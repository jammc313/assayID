#!/usr/bin/env python
# coding: utf-8

# ### This python script provides a framework to systematically find optimal window sizes and regions of mt-genomes for discriminating between a target and non-target sequences. 


# Import Libraries:
import os
import multiprocessing
import subprocess
from Bio import AlignIO
from Bio.Seq import Seq
import numpy as np
import pandas as pd
import primer3
from itertools import combinations
import matplotlib.pyplot as plt
import seaborn as sns
import random
from tqdm import tqdm
from Bio import Entrez, SeqIO
from Bio.Graphics import GenomeDiagram
from reportlab.lib.units import cm
from reportlab.lib import colors
from Bio.Graphics.GenomeDiagram import Diagram, Track, FeatureSet
from Bio.SeqFeature import SeqFeature, FeatureLocation
import configparser


# Read input from configuration file 
config = configparser.ConfigParser()
config.read('config.ini')


# Optionally run clustalo multiple sequence alignment. If "aligned_fasta" field present in config.ini file, skip this step.
unaligned_fasta = config['FILES']['unaligned_fasta']
aligned_fasta = config['FILES']['aligned_fasta']

# Retrieve number of cores available for clustalo threading
num_cores = os.cpu_count()
num_threads = max(1, num_cores - 2)

# Function to run Clustal Omega
def run_clustalo(unaligned_fasta, output_path="clustalo_aligned.fasta"):
    clustal_cmd = ["clustalo", "-i", unaligned_fasta, "-o", output_path, f"--threads={num_threads}"]
    subprocess.run(clustal_cmd)
    return output_path

# Check which option the user has provided and act accordingly
if unaligned_fasta:  # If the path to unaligned FASTA is provided
    print("Running Clustal Omega to align sequences...")
    aligned_fasta = "clustalo_aligned.fasta"  # Explicitly set the aligned_fasta path
    run_clustalo(unaligned_fasta, aligned_fasta)
    print(f"Alignment completed\nOutputted to file: {target}/clustalo_aligned.fasta")
elif aligned_fasta:  # If the path to an aligned FASTA file is provided
    print("Skipping Clustalo alignment. Using provided aligned FASTA file.")
else:
    print("No FASTA file provided. Please check your config.ini.")

working_alignment = AlignIO.read(aligned_fasta, "fasta")


# Retrieve scientific name of target taxon from config.ini
target_taxon = config['SPECIES']['target_taxon']
target = target_taxon.replace(" ", "_")


## Produce an output file to show what species are present in the multiple sequence alignment file
species_accessions = {} # Dictionary to hold species and their accession numbers

# Parse the records to extract species names and accession numbers
for record in working_alignment:
    # This assumes the record description format is 'GenBankAccessionNumber SpeciesName'
    parts = record.description.split()
    accession = parts[0]
    species_name = ' '.join(parts[1:3])  # Take the first two elements after the accession number for the species name

    # Add the accession number to the list for this species in the dictionary
    if species_name in species_accessions:
        species_accessions[species_name].add(accession)
    else:
        species_accessions[species_name] = {accession}

# Convert the species_accessions dictionary to a DataFrame
species_accessions_list = [{
    'Species': species, 
    'Accession Numbers': ';'.join(accessions),
    'Count': len(accessions)  # Count the number of accessions
} for species, accessions in species_accessions.items()]

species_accessions_df = pd.DataFrame(species_accessions_list)

# Save the DataFrame to an Excel file
species_accessions_df.to_excel(f'{target}/species_accessions.xlsx', index=False)
print(f"Outputted to file: {target}/species_accessions.xlsx\nContains list of all species and accessions present in multiple sequence alignment file")        


# Functions for extracting windows

def is_target_taxon(record, target_taxon):
    """
    Function to pass if sequence in MSA file belongs to target taxon
    """
    return target_taxon in record.description

def extend_window_to_include_nucleotides(seq, start_pos, window_size):
    """
    Extend the window from the start position until it includes window_size valid nucleotides.
    Returns the end position of the extended window.
    """
    valid_nucleotides = 0
    current_pos = start_pos
    while valid_nucleotides < window_size and current_pos < len(seq):
        if seq[current_pos] != '-':
            valid_nucleotides += 1
        current_pos += 1
    return current_pos

def check_gap_pattern_consistency(target_windows):
    """
    Ensure all target windows share the exact same gap pattern.
    """
    # Initialize a reference set of gap positions from the first sequence
    reference_gap_positions = {idx for idx, char in enumerate(target_windows[0]) if char == '-'}

    # Compare the gap positions of the remaining sequences to the reference
    for window in target_windows[1:]:
        current_gap_positions = {idx for idx, char in enumerate(window) if char == '-'}
        
        # If the gap positions differ, the pattern is inconsistent
        if current_gap_positions != reference_gap_positions:
            return False
    return True

def extract_windows(alignment, window_size, target_taxon, step_size):
    """
    Function to find all permissible windows in the MSA file.
    Windows are defined based on a random target sequence, ensuring it contains the specified number of nucleotides.
    Windows pass only if all target sequences within that window have the same pattern of gaps.
    """
    permissible_windows = []
    target_indices = [i for i, record in enumerate(alignment) if is_target_taxon(record, target_taxon)]

    total_length = len(alignment[0].seq)
    progress_bar = tqdm(total=total_length, desc="Processing windows")
    
    i = 0  # Start position in the alignment
    while i < total_length:
        # Calculate the end position for each target sequence, extending as needed to include enough valid nucleotides
        extended_end_positions = [extend_window_to_include_nucleotides(alignment[idx].seq, i, window_size) for idx in target_indices]   
        # Determine the maximum end position to ensure all required nucleotides are included
        max_end_pos = max(extended_end_positions)

        # If the maximum end position is the same as the total length, we've hit the end without finding enough valid nucleotides
        if max_end_pos >= total_length:
            print(f"No valid windows of target sequence after position {i} in alignment file.")
            break  

        # Extract sequences for the current window across all target sequences
        target_windows = [alignment[idx].seq[i:max_end_pos] for idx in target_indices]

        # Check if all target windows have the same pattern of gaps
        if check_gap_pattern_consistency(target_windows):
            # If the window is permissible, record it
            permissible_windows.append((i + 1, max_end_pos))  # Adjust start position to be 1-based
            
            # Find the first valid nucleotide position in the window for the next start
            first_valid_nucleotide_pos = next((idx for idx, char in enumerate(target_windows[0]) if char != '-'), None)
            if first_valid_nucleotide_pos is not None:
                next_start = i + first_valid_nucleotide_pos + step_size
            else:
                next_start = i + step_size  # Default case if no valid nucleotide found

            if next_start >= total_length:
                # We're at or beyond the end of the sequence; we're done
                break
        else:
            # If the window is not permissible, only increment the start position by step size
            next_start = i + step_size
        
        progress_update = next_start - i
        progress_bar.update(progress_update)  # Manually update the progress bar based on the actual progress
        i = next_start  # Update start position for next window
    
    progress_bar.close()  # Ensure to close the progress bar after completion
    return permissible_windows

# Functions for calculating diversity and distance statistics  

def calculate_shannon_entropy(sequences):
    """
    Calculate Shannon Entropy for a group of sequences, excluding gap positions
    """
    base_frequencies = {}
    total_bases = 0

    for seq in sequences:
        for base in seq:
            if base != '-': #Ignore gaps
                base_frequencies[base] = base_frequencies.get(base, 0) + 1
                total_bases += 1

    shannon_entropy = 0
    for base, freq in base_frequencies.items():
        proportion = freq / total_bases
        shannon_entropy -= proportion * np.log2(proportion)

    return shannon_entropy

def calculate_sequence_similarity(target_sequences, non_target_sequences):
    """
    Calculate sequence similarity, considering only nucleotide positions in target sequences
    and applying a partial similarity score for gaps in non-target sequences.
    """
    def similarity_score(target_seq, non_target_seq):
        nucleotide_positions = [i for i, base in enumerate(target_seq) if base != '-']
        score = 0
        gap_penalty = 0.25  # Similarity score for gap positions in non-target sequences

        for pos in nucleotide_positions:
            if non_target_seq[pos] == target_seq[pos]:
                score += 1  # Full score for a match
            elif non_target_seq[pos] == '-':
                score += gap_penalty  # Partial score for a gap in non-target sequence

        # Normalize by the number of nucleotide positions in the target sequence
        return score / len(nucleotide_positions) if nucleotide_positions else 0

    similarity_scores = []

    # Calculate similarity score for each target vs. non-target sequence pair
    for target_seq in target_sequences:
        scores = [similarity_score(target_seq, non_target_seq) for non_target_seq in non_target_sequences]
        average_score = np.mean(scores) if scores else 0
        similarity_scores.append(average_score)

    # Return the average similarity score across all target vs. non-target comparisons
    return np.mean(similarity_scores) if similarity_scores else 0

def calculate_nucleotide_diversity(sequences):
    """
    Calculate nucleotide diversity within a group of sequences, excluding positions with gaps.
    """
    pair_counts = 0
    mismatch_counts = 0

    for seq1, seq2 in combinations(sequences, 2):
        for b1, b2 in zip(seq1, seq2):
            if b1 != '-' and b2 != '-':  # Only consider positions with nucleotides in both sequences
                pair_counts += 1
                if b1 != b2:
                    mismatch_counts += 1

    nucleotide_diversity = mismatch_counts / pair_counts if pair_counts > 0 else 0
    return nucleotide_diversity

# Functions for inferring consensus sequences and adding to dataframe

def add_consensus_sequences_to_df(dataframe, alignment_file, target_taxon):
    """
    Adds a column to dataframe containing consensus sequence of target taxa sequences within that window 
    Progress indicated with a tqdm bar
    """
    alignment = AlignIO.read(alignment_file, "fasta")
    target_indices = [i for i, record in enumerate(alignment) if target_taxon in record.description]

    # Wrap dataframe.iterrows() with tqdm for progress tracking
    for index, row in tqdm(dataframe.iterrows(), total=dataframe.shape[0], desc="Generating Consensus sequences for all windows"):
        start, end = int(row['Start']), int(row['End'])
        window_sequences = [alignment[i].seq[(start-1):end] for i in target_indices]
        consensus_seq = get_consensus_sequence(window_sequences)
        dataframe.at[index, 'Consensus Sequence'] = consensus_seq

    return dataframe

def get_consensus_sequence(sequences):
    """
    Generates a consensus sequence from a list of sequences, excluding gap positions.
    Chooses the most frequent nucleotide (A, C, G, T) at each position, ignoring gaps.
    """
    consensus = ""
    length = len(sequences[0])
    for i in range(length):
        pos_counts = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
        for seq in sequences:
            base = seq[i].upper()
            if base in pos_counts:  # Only count A, C, G, T
                pos_counts[base] += 1
        
        if sum(pos_counts.values()) == 0:  # If all are gaps, consider as a gap in consensus
            consensus += '-'
        else:
            consensus += max(pos_counts, key=pos_counts.get)
    
    consensus = consensus.replace('-', '') # exclude gaps
    return consensus


# Parse window sizes and step size from config.ini
min_window_size, max_window_size = [int(size) for size in config['PARAMETERS']['window_sizes'].split(',')]
step_size = int(config['PARAMETERS']['step_size'])

# Construct the range of window sizes
window_sizes = list(range(min_window_size, max_window_size + step_size, step_size))
print(f"Window sizes to assess for assay design: {window_sizes}")

window_data = []
# Iterate across window sizes and store permissible windows in DataFrame
for size in window_sizes:
    windows = extract_windows(working_alignment, size, target_taxon, step_size)
    for start, end in windows:
        window_data.append({'Window Size': size, 'Start': start, 'End': end})
window_df = pd.DataFrame(window_data) # Convert list to dataframes


# Add consensus sequence for windows that pass. 
window_df = add_consensus_sequences_to_df(window_df, aligned_fasta, target_taxon)


# Use Primer3 to design optimal primer/probes for windows

def design_primer(seq_id, sequence):
    """
    Set criteria for Primer3 assay design.
    Retrieve parameters from config.ini.
    """
    return primer3.bindings.design_primers(
        seq_args={
            'SEQUENCE_ID': seq_id,
            'SEQUENCE_TEMPLATE': sequence,
        },
        global_args={
            'PRIMER_MIN_SIZE': int(config['PRIMER3']['PRIMER_MIN_SIZE']),
            'PRIMER_OPT_SIZE': int(config['PRIMER3']['PRIMER_OPT_SIZE']),
            'PRIMER_MAX_SIZE': int(config['PRIMER3']['PRIMER_MAX_SIZE']),
            'PRIMER_PICK_LEFT_PRIMER': int(config['PRIMER3']['PRIMER_PICK_LEFT_PRIMER']),
            'PRIMER_PICK_INTERNAL_OLIGO': int(config['PRIMER3']['PRIMER_PICK_INTERNAL_OLIGO']),
            'PRIMER_PICK_RIGHT_PRIMER': int(config['PRIMER3']['PRIMER_PICK_RIGHT_PRIMER']),
            'PRIMER_INTERNAL_MAX_SELF_END': int(config['PRIMER3']['PRIMER_INTERNAL_MAX_SELF_END']),
            'PRIMER_INTERNAL_MIN_SIZE': int(config['PRIMER3']['PRIMER_INTERNAL_MIN_SIZE']),
            'PRIMER_INTERNAL_MAX_SIZE': int(config['PRIMER3']['PRIMER_INTERNAL_MAX_SIZE']),
            'PRIMER_PRODUCT_SIZE_RANGE': eval(config['PRIMER3']['PRIMER_PRODUCT_SIZE_RANGE']),
            'PRIMER_PRODUCT_OPT_SIZE': int(config['PRIMER3']['PRIMER_PRODUCT_OPT_SIZE']),
            'PRIMER_OPT_TM': float(config['PRIMER3']['PRIMER_OPT_TM']),
            'PRIMER_MIN_TM': float(config['PRIMER3']['PRIMER_MIN_TM']),
            'PRIMER_MAX_TM': float(config['PRIMER3']['PRIMER_MAX_TM']),
            'PRIMER_INTERNAL_MIN_TM': float(config['PRIMER3']['PRIMER_INTERNAL_MIN_TM']),
            'PRIMER_INTERNAL_MAX_TM': float(config['PRIMER3']['PRIMER_INTERNAL_MAX_TM']),
            'PRIMER_INTERNAL_OPT_TM': float(config['PRIMER3']['PRIMER_INTERNAL_OPT_TM']),
            'PRIMER_MIN_GC': float(config['PRIMER3']['PRIMER_MIN_GC']),
            'PRIMER_MAX_GC': float(config['PRIMER3']['PRIMER_MAX_GC']),
            'PRIMER_OPT_GC_PERCENT': float(config['PRIMER3']['PRIMER_OPT_GC_PERCENT']),
            'PRIMER_MAX_POLY_X': int(config['PRIMER3']['PRIMER_MAX_POLY_X']),
            'PRIMER_INTERNAL_OPT_GC_PERCENT': float(config['PRIMER3']['PRIMER_INTERNAL_OPT_GC_PERCENT']),
            'PRIMER_INTERNAL_MIN_GC': float(config['PRIMER3']['PRIMER_INTERNAL_MIN_GC']),
            'PRIMER_INTERNAL_MAX_GC': float(config['PRIMER3']['PRIMER_INTERNAL_MAX_GC']),
            'PRIMER_INTERNAL_MAX_POLY_X': int(config['PRIMER3']['PRIMER_INTERNAL_MAX_POLY_X']),
            'PRIMER_MAX_END_GC': int(config['PRIMER3']['PRIMER_MAX_END_GC']),
            'PRIMER_SALT_MONOVALENT': float(config['PRIMER3']['PRIMER_SALT_MONOVALENT']),
            'PRIMER_DNA_CONC': float(config['PRIMER3']['PRIMER_DNA_CONC']),
            'PRIMER_MAX_NS_ACCEPTED': int(config['PRIMER3']['PRIMER_MAX_NS_ACCEPTED']),
            'PRIMER_MAX_SELF_ANY': int(config['PRIMER3']['PRIMER_MAX_SELF_ANY']),
            'PRIMER_MAX_SELF_END': int(config['PRIMER3']['PRIMER_MAX_SELF_END']),
            'PRIMER_PAIR_MAX_COMPL_ANY': int(config['PRIMER3']['PRIMER_PAIR_MAX_COMPL_ANY']),
            'PRIMER_PAIR_MAX_COMPL_END': int(config['PRIMER3']['PRIMER_PAIR_MAX_COMPL_END'])
        }
    )


def design_and_assess_primers(dataframe):
    """
    Design primer/probes for windows using Primer3
    Handle cases where primers/probe design fails
    Assess designed assays for tendency to form secondary structures
    """
    for index, row in tqdm(dataframe.iterrows(), total=dataframe.shape[0], desc="Designing Primers"):
        seq_id = str(index)
        sequence = row['Consensus Sequence']

        # Design primers
        primers = design_primer(seq_id, sequence)

        # Check if primers were successfully designed
        if 'PRIMER_LEFT_0_SEQUENCE' in primers and 'PRIMER_RIGHT_0_SEQUENCE' in primers:
            left_primer = primers['PRIMER_LEFT_0_SEQUENCE']
            right_primer = primers['PRIMER_RIGHT_0_SEQUENCE']
            internal_oligo = primers['PRIMER_INTERNAL_0_SEQUENCE']
            primer_L_Tm = primers['PRIMER_LEFT_0_TM']
            primer_R_Tm = primers['PRIMER_RIGHT_0_TM']
            oligo_Tm = primers['PRIMER_INTERNAL_0_TM']
            primer_L_GC = primers['PRIMER_LEFT_0_GC_PERCENT']
            primer_R_GC = primers['PRIMER_RIGHT_0_GC_PERCENT']
            oligo_GC = primers['PRIMER_INTERNAL_0_GC_PERCENT']

            # Assess secondary structures
            left_hairpin = primer3.bindings.calc_hairpin(left_primer)
            right_hairpin = primer3.bindings.calc_hairpin(right_primer)
            oligo_hairpin = primer3.bindings.calc_hairpin(internal_oligo)

            dataframe.at[index, 'L Primer'] = left_primer
            dataframe.at[index, 'R Primer'] = right_primer
            dataframe.at[index, 'Internal oligo'] = internal_oligo
            dataframe.at[index, 'L Hairpin'] = float(left_hairpin.structure_found)
            dataframe.at[index, 'R Hairpin'] = float(right_hairpin.structure_found)
            dataframe.at[index, 'Oligo Hairpin'] = float(oligo_hairpin.structure_found)
            dataframe.at[index, 'L Primer Tm'] = primer_L_Tm
            dataframe.at[index, 'R Primer Tm'] = primer_R_Tm
            dataframe.at[index, 'Oligo Tm'] = oligo_Tm
            dataframe.at[index, 'L Primer GC'] = primer_L_GC
            dataframe.at[index, 'R Primer GC'] = primer_R_GC
            dataframe.at[index, 'Oligo GC'] = oligo_GC
            
        else:
            # Handle the case where primers are not found
            dataframe.at[index, 'L Primer'] = 'No Primer Found'
            dataframe.at[index, 'R Primer'] = 'No Primer Found'
            dataframe.at[index, 'Internal oligo'] = 'No Oligo Found'
            dataframe.at[index, 'L Hairpin'] = 'N/A'
            dataframe.at[index, 'R Hairpin'] = 'N/A'
            dataframe.at[index, 'Oligo Hairpin'] = 'N/A'
            dataframe.at[index, 'L Primer Tm'] = 'N/A'
            dataframe.at[index, 'R Primer Tm'] = 'N/A'
            dataframe.at[index, 'Oligo Tm'] = 'N/A'
            dataframe.at[index, 'L Primer GC'] = 'N/A'
            dataframe.at[index, 'R Primer GC'] = 'N/A'
            dataframe.at[index, 'Oligo GC'] = 'N/A'

    return dataframe


# Primer3 design primers and probes for all windows
primerdesign_windows = design_and_assess_primers(window_df)


# Filter out windows were primer/probe design failed, and sort based on secondary structure formation.
remove_NA_df = primerdesign_windows[(primerdesign_windows['L Primer'] != 'No Primer Found') & 
                 (primerdesign_windows['R Primer'] != 'No Primer Found') & 
                 (primerdesign_windows['Internal oligo'] != 'No Oligo Found')]

# drop duplicated sets of primers/probe
filtered_df = remove_NA_df.drop_duplicates(subset=['L Primer', 'R Primer', 'Internal oligo']).copy()

# Filter out rows where 'Internal oligo' has a 'G' at 5' end
filtered_df = filtered_df[~filtered_df['Internal oligo'].str.startswith('G')]

# Further filter out rows where either 'L Hairpin' or 'R Hairpin' is "1"
filtered_df = filtered_df[(filtered_df['L Hairpin'] != 1) & (filtered_df['R Hairpin'] != 1)]

# Sort by hairpin scores
sorted_df = filtered_df.sort_values(by=['L Hairpin', 'R Hairpin', 'Oligo Hairpin'], ascending=True)


# Get complement of reverse primer, and pull out the start, stop and product length of region covered by designed assays

# Define a function to calculate the reverse complement
def reverse_complement(seq):
    return str(Seq(seq).reverse_complement())

# Define a function to find the region covered by the primers
def find_primer_region(consensus, l_primer, r_primer, window_start):
    l_start = consensus.find(l_primer)
    r_primer_revcomp = reverse_complement(r_primer)
    r_start = consensus.find(r_primer_revcomp)
    
    # Calculate the absolute positions in the alignment file
    product_start = window_start + l_start
    product_end = window_start + r_start + len(r_primer) - 1
    product_length = product_end - product_start + 1
    
    return product_start, product_end, product_length

# Now apply the updated function to the dataframe
sorted_df[['Product Start', 'Product End', 'Product Length']] = sorted_df.apply(
    lambda row: find_primer_region(
        row['Consensus Sequence'], row['L Primer'], row['R Primer'], row['Start']
    ), axis=1, result_type='expand')


def calculate_metrics_for_region(row, alignment, target_taxon):
    """
    Calculate diversity & distance metrics for assay covered regions of MSA file
    """
    start_index = row['Product Start'] - 1 # adjust for python's 0-based indexing
    end_index = row['Product End']
    
    # Extract sequences for target and non-target groups within the region
    target_sequences = [str(seq.seq[start_index:end_index]) for seq in alignment if is_target_taxon(seq, target_taxon)]
    non_target_sequences = [str(seq.seq[start_index:end_index]) for seq in alignment if not is_target_taxon(seq, target_taxon)]

    # Calculate metrics
    shannon_entropy_target = calculate_shannon_entropy(target_sequences)
    shannon_entropy_non_target = calculate_shannon_entropy(non_target_sequences)
    sequence_similarity = calculate_sequence_similarity(target_sequences, non_target_sequences)            
    nucleotide_diversity_target = calculate_nucleotide_diversity(target_sequences)
    nucleotide_diversity_non_target = calculate_nucleotide_diversity(non_target_sequences)

    return pd.Series([shannon_entropy_target, shannon_entropy_non_target, sequence_similarity, 
                      nucleotide_diversity_target, nucleotide_diversity_non_target])

# Apply the function to each row in the DataFrame
metrics_columns = ['Shannon Entropy Target', 'Shannon Entropy Non-Target', 'Sequence similarity', 
                   'Nucleotide Diversity Target', 'Nucleotide Diversity Non-Target']

sorted_df[metrics_columns] = sorted_df.apply(calculate_metrics_for_region, axis=1, args=(working_alignment, target_taxon))


# Plot distribution of Shannon Entropy for Target Sequences
plt.figure(figsize=(10, 6))
sns.histplot(sorted_df['Shannon Entropy Target'], color='green', kde=True)
plt.title('Distribution of Shannon Entropy for Target Sequences')
plt.xlabel('Shannon Entropy Target')
plt.ylabel('Frequency')
plt.savefig(f'{target}/shannon_entropy_target_distribution.png')
plt.close()


# Plot distribution of Shannon Entropy for Non-Target Sequences
plt.figure(figsize=(10, 6))
sns.histplot(sorted_df['Shannon Entropy Non-Target'], color='blue', kde=True)
plt.title('Distribution of Shannon Entropy for Non-Target Sequences')
plt.xlabel('Shannon Entropy Non-Target')
plt.ylabel('Frequency')
plt.savefig(f'{target}/shannon_entropy_non_target_distribution.png')
plt.close()


# Plot distribution of Sequence similarity
plt.figure(figsize=(10, 6))
sns.histplot(sorted_df['Sequence similarity'], color='red', kde=True)
plt.title('Distribution of Target vs NonTarget Sequence Similarity')
plt.xlabel('Sequence Similarity')
plt.ylabel('Frequency')
plt.savefig(f'{target}/sequence_similarity_distribution.png')
plt.close()


# Plot Nucleotide Diversity Target
plt.figure(figsize=(10, 6))
sns.histplot(sorted_df['Nucleotide Diversity Target'], color='green', kde=True)
plt.title('Nucleotide Diversity Target')
plt.xlabel('Nucleotide Diversity Target')
plt.ylabel('Frequency')
plt.savefig(f'{target}/nucleotide_diversity_target_distribution.png')
plt.close()


# Plot Nucleotide Diversity Non-Target
plt.figure(figsize=(10, 6))
sns.histplot(sorted_df['Nucleotide Diversity Non-Target'], color='blue', kde=True)
plt.title('Nucleotide Diversity Non-Target')
plt.xlabel('Nucleotide Diversity Non-Target')
plt.ylabel('Frequency')
plt.savefig(f'{target}/nucleotide_diversity_non_target_distribution.png')
plt.close()


# Plot and save pairwise correlations of all metrics.
metric_columns = ['Shannon Entropy Target','Shannon Entropy Non-Target','Sequence similarity', 
                  'Nucleotide Diversity Target', 'Nucleotide Diversity Non-Target']
# Create a subset dataframe with these columns
metrics_df = sorted_df[metric_columns]
# Create pairplot
pairplot = sns.pairplot(metrics_df)
plt.savefig(f'{target}/pairwise_comparisons.png')
plt.close()


# Multi-objective optimization
# Here implementing a TOPSIS-based ranking of assays. Uses the calculated statistics for regions covered by each assay to see how distant each assay is from hypothetical "ideal" and "non-ideal" assays.

# Define your ideal statistics
ideal_stats = {
    'Shannon Entropy Target': 0,
    'Shannon Entropy Non-Target': 2,  # Assuming entropy is measured in bits
    'Sequence similarity': 0,
    'Nucleotide Diversity Target': 0,
    'Nucleotide Diversity Non-Target': 1
}

worst_stats = {
    'Shannon Entropy Target': 2,
    'Shannon Entropy Non-Target': 0,  
    'Sequence similarity': 1,
    'Nucleotide Diversity Target': 1,
    'Nucleotide Diversity Non-Target': 0
}


# Preferentially define weights for each metric, with higher weight for "Sequence similarity"
weights = {
    'Shannon Entropy Target': 1.0,  # Important, but standard weight
    'Shannon Entropy Non-Target': 1.0,  # Less important
    'Sequence similarity': 2.0,  # Most important
    'Nucleotide Diversity Target': 1.0,
    'Nucleotide Diversity Non-Target': 1.0  # Less important
}

# Adjusted normalization using weights
def weighted_normalize(df, ideal, worst, weights):
    normalized_df = pd.DataFrame()
    for col in df.columns:
        if col in ideal and col in worst and col in weights:
            # Apply weights during normalization
            if ideal[col] < worst[col]:  # Maximize
                normalized_df[col] = ((df[col] - ideal[col]) / (worst[col] - ideal[col])) * weights[col]
            else:  # Minimize
                normalized_df[col] = ((df[col] - worst[col]) / (ideal[col] - worst[col])) * weights[col]
    return normalized_df

normalized_df = weighted_normalize(sorted_df, ideal_stats, worst_stats, weights)

# Calculate the Euclidean distance to the ideal and worst profiles
ideal_vector = np.array([ideal_stats[col] for col in normalized_df.columns])
worst_vector = np.array([worst_stats[col] for col in normalized_df.columns])

distance_to_ideal = np.sqrt(((normalized_df - ideal_vector) ** 2).sum(axis=1))
distance_to_worst = np.sqrt(((normalized_df - worst_vector) ** 2).sum(axis=1))

# Calculate the similarity to the ideal solution
similarity_to_ideal = distance_to_worst / (distance_to_ideal + distance_to_worst)

# Add the similarity score to the original DataFrame
sorted_df['Similarity to Ideal'] = similarity_to_ideal

# Rank the assays based on their similarity to the ideal solution
sorted_df['TOPSIS Rank'] = sorted_df['Similarity to Ideal'].rank(method='max', ascending=False)

# Sort the DataFrame based on the TOPSIS Rank
ranked_df = sorted_df.sort_values(by='TOPSIS Rank')

# Save the DataFrame to an Excel file
ranked_df.to_excel(f'{target}/ranked_assays.xlsx', index=False)
print(f"Outputted to file: {target}/ranked_assays.xlsx\nContains TOPSIS ranked assays")   


# Optionally produce diagram showing where top ranked assays bind in mt-genome
# Note: this only happens if USER has specified an accession number in conig.ini of an annotated mitochondrial genome belonging to the target species.

# Check if the accession number is provided in the config.ini
reference_genome_accession = config['SPECIES']['reference_genome_accession']
if reference_genome_accession:  # If an accession number is provided
    Entrez.email = config['DEFAULT']['email']

    def fetch_genome(accession):
        handle = Entrez.efetch(db="nucleotide", id=accession, rettype="gb", retmode="text")
        record = SeqIO.read(handle, "genbank")
        handle.close()
        return record
    
    mt_genome_record = fetch_genome(reference_genome_accession)
    
    # Extract the top 10 assays
    top_10_assays = ranked_df.head(10)
    
    # Extract their start and end positions as tuples (start, end)
    top_10_assay_positions = [(row['Product Start'], row['Product End']) for index, row in top_10_assays.iterrows()]
    
    def create_circular_genome_diagram(record, top_10_assay_positions):
        gd_diagram = GenomeDiagram.Diagram("Mitochondrial Genome")
        gd_track_for_features = gd_diagram.new_track(1, name="Annotated Features")
        gd_feature_set = gd_track_for_features.new_set()
    
        # Add genome features
        for feature in record.features:
            if feature.type == "gene":
                gd_feature_set.add_feature(feature, color=colors.skyblue, label=True, label_size=10, label_angle=0)
            elif feature.type in ["tRNA", "rRNA"]:
                gd_feature_set.add_feature(feature, color=colors.lightgreen, label=True, label_size=8, label_angle=0)
            elif feature.type == "D-loop":
                gd_feature_set.add_feature(feature, color=colors.darkblue, label=True, label_size=10, label_angle=0)

        # Highlight top 10 assays
        gd_track_for_assays = gd_diagram.new_track(3, name="Assay Regions", greytrack=True)
        gd_assays = gd_track_for_assays.new_set()  
        for i, (start, end) in enumerate(top_10_assay_positions, start=1):
            # Correctly create a feature for each assay region, including strand information
            feature_location = FeatureLocation(start, end, strand=+1)
            assay_feature = SeqFeature(feature_location)
            # Add the feature with a specific style
            gd_assays.add_feature(assay_feature, color=colors.red, label=True, label_position="middle", label_size=10, label_color=colors.black, name=f"Assay{i}",
                                  sigil="BIGARROW", arrowshaft_height=1, arrowhead_length=1)
   
        # Draw the diagram
        gd_diagram.draw(format="circular", circular=True, pagesize=(25*cm, 25*cm),
                        start=0, end=len(record), circle_core=0.7)
        gd_diagram.write(f'{target}/mtgenome_assays.png', "PNG")
        print(f"Outputted to file: {target}/mtgenome_assays.png\nContains a diagram of an annotated mitochondrial genome showing where top 10 ranked assays target")

    # Example usage (you need to define top_15_assay_positions based on your data)
    create_circular_genome_diagram(mt_genome_record, top_10_assay_positions)
else:
    print("No reference genome accession number provided. Skipping the mt-genome diagram generation.")


print("Assay design completed")

