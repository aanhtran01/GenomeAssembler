# -*- coding: utf-8 -*-
"""
Project goal: build a genome assembler

The input to the project are genome sequencing reads and the output is a predicted genome sequence.

Project builds a spectrum from the sequencing reads and arranges them into a de bruijn graph. The graph is then traversed using DFS. This returns fragments of the predicted reference genome. A simplified overlap-layout-consensus method is used to reconstruct a predicted reference genome. Reads are then aligned onto the reference genome.

the output is a file contains the headers of the reads sorted in the order that they appear in the genome
"""


import random
import argparse
from copy import deepcopy


#set up a command line interface with argparse, create argparse objects
parser = argparse.ArgumentParser(description='Process sequencing reads.')
parser.add_argument('read_file', type=str, help='path to reference file')


args = parser.parse_args()

read_file= args.read_file


#function to read in fasta files 
def read_fasta_file(filename):
    with open(filename, 'r') as f:
        read_id = None
        read_seq = ""
        for line in f:
            if line.startswith(">"):
                if read_id is not None:
                    yield (read_id, read_seq)
                read_id = line.strip()[1:]
                read_seq = ""
            else:
                read_seq += line.strip()
        if read_id is not None:
            yield (read_id, read_seq)

# Generate a spectrum
def generate_spectrum(reads, k):
    spectrum = {}  # Dictionary to store sequence information
    for read_id, read_seq in reads:
        # Step 1: Break the read into k-mers
        kmers = [read_seq[i:i + k] for i in range(len(read_seq) - k + 1)]
        for kmer in kmers:
            if kmer is not None:
                if kmer in spectrum:
                    spectrum[kmer]["count"] += 1
                    spectrum[kmer]["read_ids"].append(read_id)
                else:
                    spectrum[kmer] = {"count": 1, "read_ids": [read_id]}
    return spectrum

def filter_spectrum(spectrum):
    filtered_spectrum = {kmer: data for kmer, data in spectrum.items() if data["count"] >= 6}
    return filtered_spectrum


#Build a de Bruijn graph:
def de_bruijn(k, patterns):
        adjacency_db = {}
        for p in patterns:
            prefix = p[:k - 1]
            suffix = p[1:]
            adjacency_db.setdefault(prefix, []).append(suffix)
            if suffix not in adjacency_db:
                adjacency_db[suffix] = []
        return adjacency_db
        
def reconstruct_from_path(path):
        return path[0] + ''.join(seq[-1] for seq in path[1:])
        

def dfs(adjacency_db, start_node, visited):
    stack = [start_node]  # Use a stack to keep track of nodes to visit

    while stack:
        node = stack.pop()  # Pop the last node from the stack

        if node not in visited:
            visited.append(node)  # Mark the node as visited

            # Traverse the neighbors (suffixes) of the current node
            for neighbor in adjacency_db.get(node, []):
                stack.append(neighbor)  # Add neighbors to the stack

def generate_reference_kmers(reference_genome, kmer_size):
    reference_kmers = {}
    for i in range(len(reference_genome) - kmer_size + 1):
        kmer = reference_genome[i:i+kmer_size]
        position = i  # Adjust the position to be 1-indexed
        if kmer in reference_kmers:
            reference_kmers[kmer].append(position)
        else:
            reference_kmers[kmer] = [position]
    return reference_kmers

def HammingDistance(seq1, seq2):
    return len([i for i in range(len(seq1)) if seq1[i] != seq2[i]])

#function to align a read to the reference genome
def align_read_to_genome(read, reference_genome, k, reference_kmers):

    # Step 1: Break the read into k-mers
    kmers = [read[i:i+k] for i in range(len(read)-k+1)]
    #first_kmer = kmers[0]


    # Step 2: Search for matches in the BurrowsWheeler index and extend the alignment
    best_match = None
    best_score = float('inf')
    for i, kmer in enumerate(kmers):
        positions = reference_kmers.get(kmer)

        if positions is not None:  # Check if positions exist
            for pos in positions:
                # extend the alignment
                offset = i
                alignment_start = pos - offset
                alignment_end = alignment_start + len(read)
                if alignment_start < 0 or alignment_end > len(reference_genome):
                    continue  # alignment out of bounds, skip to next position
                ref_sequence = reference_genome[alignment_start:alignment_end]
                score = HammingDistance(read, ref_sequence)

                # check if this is the best match so far
                if score < best_score:
                    best_score = score
                    best_match = alignment_start

    return best_match, best_score

#function to alogn all reads to the genome
def align_all_reads_to_genome(donor_reads, reference_genome, k):
    results = []
    for read_id, read_seq in donor_reads:
        best_match, best_score = align_read_to_genome(read_seq, reference_genome, k, reference_kmers)
        results.append({'donor_read_id': read_id,'sequence' : read_seq ,'best_match': best_match, 'best_score': best_score})
    return results


k_spectrum = 20
reads = read_fasta_file(read_file)
spectrum = generate_spectrum(reads, k_spectrum)
filtered_spectrum = filter_spectrum(spectrum)

de_bruijn = de_bruijn(k_spectrum, filtered_spectrum)
#de_bruijn(k, filtered_spectrum)

# Set the seed for the random number generator
random.seed(123)

visited_dict = {}
fragments = []

# Get a random sample of 1000 keys from de_bruijn dictionary
sample_keys = int(len(de_bruijn.keys())*0.05)
#random_keys = random.sample(de_bruijn.keys(), 100)
random_keys = random.sample(list(de_bruijn.keys()), sample_keys)

sequence_lengths = []  # List to store sequence lengths

for key in random_keys:
    visited = []
    dfs(de_bruijn, key, visited)
    visited_dict[key] = visited

    reconstructed_sequence = reconstruct_from_path(visited)
    sequence_length = len(reconstructed_sequence)

    if sequence_length >= 0:
        fragments.append(reconstructed_sequence)
        sequence_lengths.append(sequence_length)


# Filter fragments with less than desired threshold length
#threshold_length = 18250
#threshold_length = 850
#filtered_fragments = [fragment for fragment in fragments if len(fragment) >= threshold_length]

# Calculate the lengths of the filtered fragments
#fragment_lengths = [len(fragment) for fragment in filtered_fragments]

# Sort the fragments based on their size and only use the top longest sequences
sorted_fragments = sorted(fragments, key=len, reverse=True)
top_longest = sorted_fragments[:1]

# Store the longest sequence as a string
longest_sequence = ''.join(top_longest)



reference_genome = longest_sequence

#calling the functions
k = 16

#filename = "/content/project3b_20000_reads_without_positions.fasta"

donor_reads = read_fasta_file(read_file)

reference_kmers = generate_reference_kmers(reference_genome, k)
results = align_all_reads_to_genome(donor_reads, reference_genome, k)
#print(results)

# Filter out results with None as the best_match value
filtered_results = [result for result in results if result['best_match'] is not None]
#print(filtered_results)

# Sort the filtered results based on the 'best_match' value
sorted_results = sorted(filtered_results, key=lambda x: x['best_match'])


# Open a text file for writing
with open('predictions.txt', 'w') as file:
    for result in sorted_results:
        file.write(f'>{result["donor_read_id"]}\n')
