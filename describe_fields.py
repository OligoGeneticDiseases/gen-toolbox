import hail as hl
import sys
from pathlib import Path

def load_and_describe_vcf(file_path, destination):
    # Initialize Hail with temporary directory
    hl.init(tmp_dir=destination)

    # Load the VCF or GVCF file
    if file_path.endswith(".gvcf"):
        table = hl.import_vcf(file_path, force=True, reference_genome="GRCh37")
    else:
        print(f"Unsupported file format: {file_path}")
        sys.exit(1)

    # Describe the table to see available fields
    print("\n--- Table Description ---")
    table.describe()
    
    # Check if the CSQ field is present in rows
    if 'CSQ' in table.row:
        print("\nCSQ field is present. Sample data from CSQ field:")
        # Print the first few values of CSQ to understand its structure
        csq_data = table.select_rows(table.CSQ).take(5)
        for entry in csq_data:
            print(entry.CSQ)
    else:
        print("\nCSQ field is not present in the data.")
    
    # Stop Hail
    hl.stop()

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python describe_fields.py <input_file> <temp_dir>")
        sys.exit(1)
    
    input_file = sys.argv[1]
    temp_dir = sys.argv[2]

    load_and_describe_vcf(input_file, temp_dir)

