#!/usr/bin/env python3

import synapseclient
import os
import argparse
import sys
import shutil

def main():
    # Set up command line argument parsing
    parser = argparse.ArgumentParser(description='Download files from Synapse using ID list and ~/.synapseConfig credentials')
    parser.add_argument('-f', '--file', required=True, help='Path to file containing Synapse IDs (one per line)')
    parser.add_argument('-d', '--directory', default='raw_data', help='Directory to save downloaded files')
    parser.add_argument('-p', '--prefix', default=True, action='store_true', help='Prefix filenames with Synapse IDs')
    parser.add_argument('-r', '--replace', action='store_true', help='Replace filename completely with Synapse ID')
    
    args = parser.parse_args()
    
    # Create download directory
    os.makedirs(args.directory, exist_ok=True)
    
    # Load Synapse IDs from file
    try:
        with open(args.file, 'r') as file:
            synapse_ids = [line.strip() for line in file if line.strip() and not line.startswith('#')]
        
        if not synapse_ids:
            sys.exit(f"Error: No valid Synapse IDs found in {args.file}")
            
    except FileNotFoundError:
        sys.exit(f"Error: File {args.file} not found")
    
    # Login to Synapse using ~/.synapseConfig
    try:
        print("Using credentials from ~/.synapseConfig...")
        syn = synapseclient.login()
    except Exception as e:
        sys.exit(f"Authentication error: {e}")
    
    # Download each dataset
    for i, syn_id in enumerate(synapse_ids, 1):
        try:
            print(f"Downloading {i}/{len(synapse_ids)}: {syn_id}...")
            
            # Download file with its original name
            entity = syn.get(entity=syn_id, downloadLocation=args.directory)
            original_path = os.path.join(args.directory, entity.name)
            
            # Determine new filename based on options
            file_ext = os.path.splitext(entity.name)[1]
            if args.replace:
                new_name = f"{syn_id}{file_ext}"
            else:
                new_name = f"{syn_id}_{entity.name}"
            
            new_path = os.path.join(args.directory, new_name)
            
            # Rename the file
            shutil.move(original_path, new_path)
            print(f"✓ Downloaded: {new_name}")
            
        except Exception as e:
            print(f"× Failed to download {syn_id}: {e}")
    
    print(f"\nDownload complete. Files saved to {args.directory}")

if __name__ == "__main__":
    main()