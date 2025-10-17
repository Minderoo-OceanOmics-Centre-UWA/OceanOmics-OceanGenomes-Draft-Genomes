#!/usr/bin/env python3
import psycopg2
import configparser
import sys

def load_db_config(config_file):
    config = configparser.ConfigParser()
    config.read(config_file)

    return {
        'dbname': config.get('postgres', 'dbname'),
        'user': config.get('postgres', 'user'),
        'password': config.get('postgres', 'password'),
        'host': config.get('postgres', 'host'),
        'port': config.getint('postgres', 'port')
    }

def get_species_info(db_params, og_id):
    query = """
        SELECT s.nominal_species_id, sp.ncbi_taxon_id, sp.class
        FROM sample s
        JOIN species sp ON s.nominal_species_id = sp.species
        WHERE s.og_id = %s
    """
    try:
        conn = psycopg2.connect(**db_params)
        with conn.cursor() as cur:
            cur.execute(query, (og_id,))
            result = cur.fetchone()
            if result:
                nominal_species_id, taxon_id, tax_class = result
                print("og_id,nominal_species_id,taxon_id,class")
                print(f"{og_id},{nominal_species_id},{taxon_id},{tax_class}")
            else:
                print(f"[WARN] No matching species data found for OG ID: {og_id}")
    except Exception as e:
        print(f"[ERROR] Database query failed: {e}")
    finally:
        if conn:
            conn.close()

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage:\n  python lookup_species_info.py <db.cfg> <OG_ID>")
        sys.exit(1)

    config_file = sys.argv[1]
    og_id = sys.argv[2]

    db_params = load_db_config(config_file)
    get_species_info(db_params, og_id)