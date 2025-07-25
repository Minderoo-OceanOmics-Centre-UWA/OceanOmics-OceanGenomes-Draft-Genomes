#!/usr/bin/env python3

import psycopg2
import pandas as pd
import numpy as np
import configparser
import sys
import re
from pathlib import Path

# -------------------------------
# Load DB credentials from .cfg
# -------------------------------
def load_db_config(config_file):
    if not Path(config_file).exists():
        raise FileNotFoundError(f"❌ Config file '{config_file}' does not exist.")
    
    config = configparser.ConfigParser()
    config.read(config_file)

    if not config.has_section('postgres'):
        raise ValueError("❌ Missing [postgres] section in config file.")

    required_keys = ['dbname', 'user', 'password', 'host', 'port']
    for key in required_keys:
        if not config.has_option('postgres', key):
            raise ValueError(f"❌ Missing '{key}' in [postgres] section of config file.")

    return {
        'dbname': config.get('postgres', 'dbname'),
        'user': config.get('postgres', 'user'),
        'password': config.get('postgres', 'password'),
        'host': config.get('postgres', 'host'),
        'port': config.getint('postgres', 'port')    
    }

# -------------------------------
# Main logic
# -------------------------------
if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage:\n  push_mtdna_assm_results.py <config_file> <assembly_prefix> <log_file> <fasta_file>")
        sys.exit(1)

    config_file = sys.argv[1]
    assembly_prefix = sys.argv[2]
    log_path = Path(sys.argv[3])
    fasta_path = Path(sys.argv[4])

    # Read log file
    try:
        log_text = log_path.read_text()
    except Exception as e:
        print(f"❌ Failed to read log file: {e}")
        sys.exit(1)

    # Extract result status
    match_stats = re.findall(r"Result status of animal_mt:\s*(.+)", log_text)
    stats = match_stats[-1].strip() if match_stats else None

    # Extract average coverage
    match_avg_coverage = re.findall(r"Average animal_mt coverage =\s*([^\s]+)", log_text)
    avg_coverage = match_avg_coverage[-1].strip() if match_avg_coverage else None

    # Extract average base coverage
    match_avg_base_coverage = re.findall(r"Average animal_mt base-coverage =\s*([^\s]+)", log_text)
    avg_base_coverage = match_avg_base_coverage[-1].strip() if match_avg_base_coverage else None

    # Compute sequence length from FASTA
    try:
        with open(fasta_path) as f:
            length = sum(len(line.strip()) for line in f if not line.startswith(">"))
    except Exception as e:
        print(f"❌ Failed to read FASTA file: {e}")
        sys.exit(1)

    print(f"Stats: {stats}")
    print(f"Length: {length}")
    print(f"Avg Coverage: {avg_coverage}")
    print(f"Avg Base Coverage: {avg_base_coverage}")

    # Parse assembly_prefix
    try:
        og_id, tech, seq_date, code = assembly_prefix.split(".")
    except ValueError:
        print(f"❌ Failed to split assembly_prefix: {assembly_prefix}")
        sys.exit(1)

    try:
        db_params = load_db_config(config_file)
        conn = psycopg2.connect(**db_params)
        cursor = conn.cursor()

        upsert_query = """
        INSERT INTO mitogenome_data (
            og_id, tech, seq_date, code, stats, length, avg_coverage, avg_base_coverage
        )
        VALUES (
            %(og_id)s, %(tech)s, %(seq_date)s, %(code)s, %(stats)s, %(length)s, %(avg_coverage)s, %(avg_base_coverage)s
        )
        ON CONFLICT (og_id, tech, seq_date, code) 
        DO UPDATE SET
            stats = CASE WHEN mitogenome_data.stats IS NULL THEN EXCLUDED.stats ELSE mitogenome_data.stats END,
            length = CASE WHEN mitogenome_data.length IS NULL THEN EXCLUDED.length ELSE mitogenome_data.length END,
            avg_coverage = CASE WHEN mitogenome_data.avg_coverage IS NULL THEN EXCLUDED.avg_coverage ELSE mitogenome_data.avg_coverage END,
            avg_base_coverage = CASE WHEN mitogenome_data.avg_base_coverage IS NULL THEN EXCLUDED.avg_base_coverage ELSE mitogenome_data.avg_base_coverage END
        RETURNING stats, length, avg_coverage, avg_base_coverage
        """

        params = {
            "og_id": og_id,
            "tech": tech,
            "seq_date": seq_date,
            "code": code,
            "stats": stats,
            "length": int(length),
            "avg_coverage": float(avg_coverage) if avg_coverage is not None else None,
            "avg_base_coverage": float(avg_base_coverage) if avg_base_coverage is not None else None,
        }

        cursor.execute(upsert_query, params)
        returned = cursor.fetchone()
        conn.commit()

        print("📌 Final stored values:")
        print(dict(zip(["stats", "length", "avg_coverage", "avg_base_coverage"], returned)))

        # Field-by-field comparison
        field_names = ["stats", "length", "avg_coverage", "avg_base_coverage"]
        preserved_fields = []

        for idx, field in enumerate(field_names):
            if params[field] is None:
                continue  # Skipped intentionally
            elif returned[idx] != params[field]:
                preserved_fields.append(field)

        if preserved_fields:
            print(f"⚠️ Existing values preserved for: {', '.join(preserved_fields)}")
        else:
            print(f"✅ Success: Inserted/Updated mitogenome_data for {assembly_prefix}")


    except Exception as e:
        if 'conn' in locals():
            conn.rollback()
        print(f"❌ Database error: {e}")

    finally:
        if 'cursor' in locals():
            cursor.close()
        if 'conn' in locals():
            conn.close()
