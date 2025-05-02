#!/usr/bin/env python3

import h5py
import pandas as pd
import numpy as np
import argparse
import os
import re
import glob
import sys

def extract_gene_id(id_string):
    """Extract gene ID before the first pipe character, or ENSEMBL ID if present"""
    # First, split by pipe and take the first element
    id_part = id_string.split('|')[0].strip()
    
    # Then look for ENSEMBL transcript pattern in that part
    pattern = r'(ENS[A-Z]*\d+(?:\.\d+)?)'
    match = re.search(pattern, id_part)
    if match:
        return match.group(1)
    else:
        return id_part

def parse_gtf_for_transcript_lengths(gtf_file):
    """Parse GTF file to get transcript lengths"""
    transcript_exons = {}
    
    print(f"Parsing GTF file: {gtf_file}")
    with open(gtf_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
                
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue
                
            feature_type = fields[2]
            if feature_type != "exon":
                continue
                
            attributes = fields[8]
            transcript_match = re.search(r'transcript_id "([^"]+)"', attributes)
            if not transcript_match:
                continue
                
            transcript_id = transcript_match.group(1)
            # Remove version number if present for better matching
            transcript_id = re.sub(r'(\.[0-9]+)$', '', transcript_id)
            
            start_pos = int(fields[3])
            end_pos = int(fields[4])
            
            if transcript_id not in transcript_exons:
                transcript_exons[transcript_id] = []
            transcript_exons[transcript_id].append((start_pos, end_pos))
    
    # Calculate transcript lengths
    transcript_lengths = {}
    for transcript_id, exons in transcript_exons.items():
        # Sum of exon lengths
        length = sum(end - start + 1 for start, end in exons)
        transcript_lengths[transcript_id] = length
    
    print(f"Parsed lengths for {len(transcript_lengths)} transcripts")
    
    # Check if we have any lengths at all
    if len(transcript_lengths) == 0:
        print("ERROR: No transcript lengths were parsed from the GTF file!")
        sys.exit(1)
        
    # Print a sample of transcript lengths
    sample_transcripts = list(transcript_lengths.items())[:5]
    print("Sample transcript lengths:")
    for tid, length in sample_transcripts:
        print(f"  {tid}: {length} bp")
        
    return transcript_lengths

def process_h5_file(h5_file, transcript_lengths):
    """Process a single H5 file and return counts, TPM and library size"""
    with h5py.File(h5_file, 'r') as f:
        # Get counts
        counts = f['est_counts'][:]
        counts_int = np.array([int(count) for count in counts])
        
        # Get transcript IDs
        gene_ids_raw = f['aux/ids'][:]
        gene_ids = np.array([extract_gene_id(id.decode('utf-8')) for id in gene_ids_raw])
        
        # Calculate library size
        library_size = np.sum(counts_int)
        
        # Create a dataframe for calculations
        df = pd.DataFrame({
            'transcript_id': gene_ids,
            'count': counts_int
        })
        
        # Strip version numbers for better matching with GTF
        df['transcript_id_no_ver'] = df['transcript_id'].str.replace(r'(\.[0-9]+)$', '', regex=True)
        
        # Add transcript lengths with both full ID and ID without version
        df['length'] = df['transcript_id'].map(lambda x: transcript_lengths.get(x, 0))
        
        # If length is 0, try without version number
        df.loc[df['length'] == 0, 'length'] = df.loc[df['length'] == 0, 'transcript_id_no_ver'].map(
            lambda x: transcript_lengths.get(x, 0)
        )
        
        # Count transcripts without length information
        missing_lengths = df[df['length'] == 0].shape[0]
        if missing_lengths > 0:
            print(f"Warning: {missing_lengths} out of {len(df)} transcripts ({missing_lengths/len(df)*100:.2f}%) have no length information")
            print(f"Sample of transcripts without length: {df[df['length'] == 0]['transcript_id'].head(5).tolist()}")
        
        # Set a minimum length to avoid division by zero
        df['length'] = df['length'].apply(lambda x: max(x, 100))  # Assuming minimum transcript length of 100bp
        
        # Calculate TPM
        # Step 1: Calculate count / length for each transcript
        df['count_per_length'] = df['count'] / df['length']
        
        # Step 2: Calculate the scaling factor
        sum_cpl = df['count_per_length'].sum()
        if sum_cpl > 0:
            scaling_factor = 1e6 / sum_cpl
            df['tpm'] = df['count_per_length'] * scaling_factor
            print(f"TPM scaling factor: {scaling_factor:.2f}")
            print(f"Max TPM value: {df['tpm'].max():.2f}")
        else:
            print("ERROR: Sum of count/length is zero, cannot calculate TPM")
            df['tpm'] = 0
            
        return gene_ids, counts_int, df['tpm'].values, library_size

def main():
    parser = argparse.ArgumentParser(description='Convert H5 files to gene expression matrices (GEM)')
    parser.add_argument('-f', '--file', help='Path to single HDF5 file')
    parser.add_argument('-d', '--directory', help='Directory containing multiple HDF5 files')
    parser.add_argument('-g', '--gtf', required=True, help='Path to GTF annotation file for gene lengths')
    parser.add_argument('-o', '--output', help='Base name for output files (default: derived from input)')
    parser.add_argument('-l', '--list', action='store_true', help='List structure of HDF5 file and exit')
    
    args = parser.parse_args()
    
    # Check that at least one input option is specified
    if not args.file and not args.directory:
        print("Error: Either a single file (-f) or a directory (-d) must be specified")
        parser.print_help()
        return
    
    # List structure of a single file if requested
    if args.list and args.file:
        with h5py.File(args.file, 'r') as f:
            print(f"Structure of {args.file}:")
            for key in f.keys():
                print(f"- {key}")
                if isinstance(f[key], h5py.Group):
                    for subkey in f[key].keys():
                        print(f"  - {key}/{subkey}")
        return
    
    # Parse GTF for transcript lengths
    transcript_lengths = parse_gtf_for_transcript_lengths(args.gtf)
    
    # Determine input files
    h5_files = []
    if args.file:
        h5_files = [args.file]
    elif args.directory:
        h5_files = glob.glob(os.path.join(args.directory, "*.h5"))
        print(f"Found {len(h5_files)} H5 files in directory {args.directory}")
    
    if not h5_files:
        print("No H5 files found to process.")
        return
    
    # Set base output name
    if args.output:
        base_output = args.output
    elif args.file:
        base_output = os.path.splitext(os.path.basename(args.file))[0]
    else:
        base_output = os.path.basename(args.directory.rstrip("/"))
    
    # Process all files
    all_gene_ids = None
    counts_matrix = []
    tpm_matrix = []
    library_specs = []
    
    for h5_file in h5_files:
        filename = os.path.basename(h5_file)
        print(f"Processing {filename}...")
        
        gene_ids, counts, tpm, library_size = process_h5_file(h5_file, transcript_lengths)
        
        # Store the first gene_ids as reference
        if all_gene_ids is None:
            all_gene_ids = gene_ids
        
        # Check that gene IDs match across files
        if not np.array_equal(gene_ids, all_gene_ids):
            print(f"Warning: Gene IDs in {filename} do not match reference gene set. Skipping.")
            continue
        
        # Add to matrices
        counts_matrix.append(counts)
        tpm_matrix.append(tpm)
        
        # Add to library specs
        library_specs.append((filename, library_size))
        
        print(f"  Library size: {library_size:,} reads")
    
    # Convert lists to dataframes
    counts_df = pd.DataFrame(np.column_stack(counts_matrix), index=all_gene_ids)
    counts_df.columns = [os.path.basename(f) for f in h5_files]
    
    tpm_df = pd.DataFrame(np.column_stack(tpm_matrix), index=all_gene_ids)
    tpm_df.columns = [os.path.basename(f) for f in h5_files]
    
    # Save the matrices
    counts_output = f"{base_output}.count.gem"
    counts_df.to_csv(counts_output, sep='\t')
    print(f"Saved count matrix to {counts_output}")
    
    tpm_output = f"{base_output}.tpm.gem"
    tpm_df.to_csv(tpm_output, sep='\t')
    print(f"Saved TPM matrix to {tpm_output}")
    
    # Save library specs
    specs_output = f"{base_output}.library.specs"
    specs_df = pd.DataFrame(library_specs, columns=["filename", "library_size"])
    specs_df.to_csv(specs_output, index=False)
    print(f"Saved library specifications to {specs_output}")

if __name__ == "__main__":
    main()
