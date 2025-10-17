#!/usr/bin/env python3

import argparse
import json
import csv
import os
import sys
from pathlib import Path

def parse_json_file(json_file):
    """Parse a JSON file and return the data as a dictionary"""
    try:
        with open(json_file, 'r') as f:
            data = json.load(f)
        return data
    except Exception as e:
        print(f"Error parsing {json_file}: {e}", file=sys.stderr)
        return None

def collect_json_data(input_dir):
    """Collect all JSON data from files in the input directory"""
    json_files = Path(input_dir).glob("*.json")
    all_data = []
    
    for json_file in json_files:
        data = parse_json_file(json_file)
        if data:
            all_data.append(data)
    
    if not all_data:
        print("No valid JSON files found!", file=sys.stderr)
        sys.exit(1)
    
    return all_data

def write_csv(data_list, output_file):
    """Write the collected data to a CSV file"""
    if not data_list:
        return
    
    # Get all possible fieldnames from all records
    fieldnames = set()
    for data in data_list:
        fieldnames.update(data.keys())
    
    # Sort fieldnames for consistent output (sample_id first if present)
    fieldnames = sorted(list(fieldnames))
    if 'sample_id' in fieldnames:
        fieldnames.remove('sample_id')
        fieldnames.insert(0, 'sample_id')
    
    with open(output_file, 'w', newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        
        for data in data_list:
            writer.writerow(data)
    
    print(f"Successfully compiled {len(data_list)} JSON files into {output_file}")

def main():
    parser = argparse.ArgumentParser(
        description="Compile JSON files into a single CSV file"
    )
    parser.add_argument(
        '--input_dir', 
        required=True,
        help='Directory containing JSON files'
    )
    parser.add_argument(
        '--output', 
        required=True,
        help='Output CSV file name'
    )
    
    args = parser.parse_args()
    
    # Collect data from all JSON files
    all_data = collect_json_data(args.input_dir)
    # Sort rows by sample_id if present
    if all_data and 'sample_id' in all_data[0]:
        try:
            all_data.sort(key=lambda x: int(x['sample_id']))
        except (ValueError, KeyError):
            all_data.sort(key=lambda x: str(x.get('sample_id', '')))
    # Write to CSV
    write_csv(all_data, args.output)

if __name__ == "__main__":
    main()