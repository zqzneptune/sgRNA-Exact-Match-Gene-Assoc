# -*- coding: utf-8 -*-
import os
import argparse
import pandas as pd
import matplotlib.pyplot as plt
from Bio import SeqIO
from Bio.Seq import Seq
import sys
import time # To measure execution time
import numpy as np # For NaN checking potentially
import pyranges as pr # Import pyranges

# --- Configuration ---
DEFAULT_SGRNA_INPUT_FILE = 'sgrnas.txt'

# BW25113
# DEFAULT_TARGET_GENOME_FASTA = 'E.-coli-K-12-substr.-BW25113_Genome.fasta'
# DEFAULT_GENE_TABLE_FILE = 'E.-coli-K-12-substr.-BW25113_All-genes.txt'
# DEFAULT_OUTPUT_DIR = 'E.-coli-K-12-substr.-BW25113_output'

# MG1655
DEFAULT_TARGET_GENOME_FASTA = 'E.-coli-K-12-substr.-MG1655_Genome.fasta'
DEFAULT_GENE_TABLE_FILE = 'E.-coli-K-12-substr.-MG1655_All-genes.txt'
DEFAULT_OUTPUT_DIR = 'E.-coli-K-12-substr.-MG1655_output'


# --- Helper Functions ---
# (read_sgrnas, read_genome, read_gene_table, find_exact_matches remain the same)
def read_sgrnas(filepath):
    """Reads sgRNA sequences from a file, returns a list of unique, valid sequences."""
    print(f"\nReading sgRNA sequences from: {filepath}")
    sequences = set()
    valid_count = 0
    skipped_count = 0
    try:
        with open(filepath, 'r') as f:
            for i, line in enumerate(f):
                seq = line.strip().upper()
                if seq:
                    if all(c in 'ATCGN' for c in seq):
                        sequences.add(seq)
                        valid_count += 1
                    else:
                        print(f"WARNING: Skipping invalid sequence (line {i+1}, non-ATCGN): {seq}")
                        skipped_count += 1
    except FileNotFoundError:
        print(f"ERROR: sgRNA input file not found: {filepath}")
        sys.exit(1)

    unique_sequences = sorted(list(sequences))
    print(f"Read {valid_count} valid sequences, skipped {skipped_count} invalid sequences.")
    print(f"Found {len(unique_sequences)} unique valid sgRNA sequences.")
    if not unique_sequences:
        print("ERROR: No valid sgRNA sequences to process.")
        sys.exit(1)
    return unique_sequences

def read_genome(filepath):
    """Reads a FASTA file, returns dict {id: seq}, and checks for single sequence assumption."""
    print(f"\nReading target genome from: {filepath}")
    genome = {}
    fasta_chrom_names = []
    total_len = 0
    try:
        for record in SeqIO.parse(filepath, "fasta"):
            print(f"  - Reading record: {record.id} ({len(record.seq):,} bp)")
            genome[record.id] = str(record.seq).upper()
            fasta_chrom_names.append(record.id)
            total_len += len(record.seq)
        print(f"Finished reading genome. Total length: {total_len:,} bp across {len(genome)} record(s).")
        if not genome:
             print("ERROR: No sequences found in the genome FASTA file.")
             sys.exit(1)
        # --- Assumption Check ---
        if len(fasta_chrom_names) > 1:
            print("\nWARNING: Genome FASTA contains multiple sequences/chromosomes:")
            print(f"         {fasta_chrom_names}")
            print("         This script assumes all gene coordinates in the table belong to the *first* sequence.")
            print(f"         Using '{fasta_chrom_names[0]}' as the chromosome name for gene association.")
            print("         Results might be incorrect if genes map to other sequences.")
        elif len(fasta_chrom_names) == 0:
             print("ERROR: Could not identify any sequence IDs in the FASTA file.")
             sys.exit(1)
        assumed_chromosome = fasta_chrom_names[0]
        return genome, assumed_chromosome
    except FileNotFoundError:
        print(f"ERROR: Target genome FASTA file not found: {filepath}")
        sys.exit(1)
    except Exception as e:
        print(f"ERROR: Failed to parse genome FASTA file {filepath}: {e}")
        sys.exit(1)

def read_gene_table(filepath, expected_chromosome):
    """Reads a custom gene table, prepares it, and returns a PyRanges object."""
    print(f"\nReading gene table from: {filepath}")
    try:
        # Try reading with auto-detected separator (robust for TSV/CSV)
        try:
            # Use low_memory=False for potentially mixed types or large files
            genes_df = pd.read_csv(filepath, sep=None, engine='python', comment='#')
        except Exception as e_read:
             print(f"ERROR: Failed to read gene table '{filepath}'. Check format/separator. Error: {e_read}")
             sys.exit(1)

        print(f"Read {len(genes_df)} rows initially.")
        print(f"Columns found: {list(genes_df.columns)}") # Less verbose

        # --- Define required original column names ---
        original_required_cols = ["Gene Name", "Accession-1", "Left-End-Position", "Right-End-Position", "Product"]
        coord_cols_original = ["Left-End-Position", "Right-End-Position"]

        # --- Validate required columns exist ---
        missing_cols = [orig for orig in original_required_cols if orig not in genes_df.columns]
        if missing_cols:
            print(f"ERROR: Missing required columns in gene table '{filepath}': {missing_cols}")
            print(f"       Expected columns based on example: {original_required_cols}")
            sys.exit(1)

        # --- Filter rows with invalid/blank coordinates BEFORE renaming ---
        print("Filtering rows with invalid/blank coordinate values...") # Less verbose
        # Convert coordinate columns to numeric, coercing errors (blanks, non-numeric) to NaN
        for col in coord_cols_original:
             # Make a copy to avoid SettingWithCopyWarning if genes_df is a slice
             genes_df = genes_df.copy()
             genes_df[col] = pd.to_numeric(genes_df[col], errors='coerce')

        # Drop rows where *any* coordinate conversion failed (resulting in NaN)
        original_rows = len(genes_df)
        genes_df.dropna(subset=coord_cols_original, inplace=True)
        filtered_rows = len(genes_df)

        if original_rows > filtered_rows:
             print(f"Filtered out {original_rows - filtered_rows} rows due to invalid/blank coordinate values.")

        if genes_df.empty:
             print("ERROR: No valid gene entries remaining after filtering coordinates.")
             sys.exit(1)
        else: # Less verbose
            print(f"{filtered_rows} rows remaining after coordinate filtering.")

        # --- Now safe to convert coordinates to integers and proceed ---
        try:
            for col in coord_cols_original:
                 # Ensure column type is suitable forastype(int) - handle floats if necessary
                 if pd.api.types.is_float_dtype(genes_df[col]):
                     # Check if float has decimals - if so, warn or decide policy (round? floor? error?)
                     if not genes_df[col].eq(genes_df[col].round()).all():
                         print(f"WARNING: Column '{col}' contains non-integer numbers after filtering. They will be truncated to integers.")
                 genes_df[col] = genes_df[col].astype(int)
        except ValueError as e:
             print(f"ERROR: Could not convert coordinate columns to integer after filtering: {e}")
             print("       This might indicate very large numbers exceeding standard integer limits, or a filtering issue.")
             sys.exit(1)


        # --- Proceed with renaming ---
        rename_map = {
            "Gene Name": "Gene_Symbol",
            "Accession-1": "Locus_Tag",
            "Left-End-Position": "Start", # Now guaranteed to be integer
            "Right-End-Position": "End",   # Now guaranteed to be integer
            "Product": "Product_Description"
        }
        # Select only the needed original columns and rename them
        genes_df_renamed = genes_df[original_required_cols].copy()
        genes_df_renamed.rename(columns=rename_map, inplace=True)


        # Assign expected chromosome name
        print(f"Assigning chromosome name '{expected_chromosome}' to all genes.") # Less verbose
        genes_df_renamed['Chromosome'] = expected_chromosome

        # --- Coordinate Conversion (1-based inclusive -> 0-based start, 1-based end) ---
        # Start is already integer
        genes_df_renamed['Start'] = genes_df_renamed['Start'] - 1
        # End is already integer - PyRanges uses it as exclusive relative to 0-based start

        # Add dummy Strand column (required by pyranges, but info is missing)
        genes_df_renamed['Strand'] = '+' # Or '.'
        print("WARNING: Gene table lacks Strand information. Using dummy strand '+'. Nearest gene direction cannot be determined relative to gene orientation.") # Less verbose

        # Select final columns for PyRanges object
        pr_cols = ['Chromosome', 'Start', 'End', 'Strand', 'Gene_Symbol', 'Locus_Tag', 'Product_Description']
        # Ensure Start and End are non-negative after conversion
        genes_final_df = genes_df_renamed[pr_cols].copy()
        genes_final_df = genes_final_df[genes_final_df['Start'] >= 0] # Filter out negative Start coords
        # Also ensure Start < End
        genes_final_df = genes_final_df[genes_final_df['Start'] < genes_final_df['End']]

        rows_after_coord_check = len(genes_final_df)
        if rows_after_coord_check < filtered_rows:
            print(f"Filtered out {filtered_rows - rows_after_coord_check} rows with invalid coordinates (Start < 0 or Start >= End).")

        if genes_final_df.empty:
            print("ERROR: No valid gene entries remaining after final coordinate checks.")
            sys.exit(1)

        # Convert to PyRanges object
        genes_pr = pr.PyRanges(genes_final_df)
        print(f"Prepared {len(genes_pr)} gene features for analysis from table.")
        return genes_pr

    except FileNotFoundError:
        print(f"ERROR: Gene table file not found: {filepath}")
        sys.exit(1)
    except Exception as e:
        print(f"ERROR: Failed to read or process gene table file {filepath}: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

def find_exact_matches(sgrna_seq, genome_sequences):
    """ Finds the first exact match (fwd or revcomp) and returns location. """
    sgrna_len = len(sgrna_seq)
    if sgrna_len == 0: return None
    try:
        rc_sgrna_seq = str(Seq(sgrna_seq).reverse_complement())
    except Exception: rc_sgrna_seq = None

    for chrom_id, chrom_seq in genome_sequences.items():
        index_fwd = chrom_seq.find(sgrna_seq)
        if index_fwd != -1:
            start = index_fwd + 1
            end = start + sgrna_len - 1
            return {'chromosome': chrom_id, 'start': start, 'end': end, 'strand': '+'} # Return match strand
        if rc_sgrna_seq:
            index_rev = chrom_seq.find(rc_sgrna_seq)
            if index_rev != -1:
                start = index_rev + 1
                end = start + sgrna_len - 1
                return {'chromosome': chrom_id, 'start': start, 'end': end, 'strand': '-'} # Return match strand
    return None


# --- Main Execution ---
def main():
    parser = argparse.ArgumentParser(description="Find exact matches of sgRNA sequences in a target genome and associate with overlapping/nearest genes using a custom table.")
    parser.add_argument("-s", "--sgrna_file", default=DEFAULT_SGRNA_INPUT_FILE,
                        help=f"Input file containing sgRNA sequences, one per line (default: {DEFAULT_SGRNA_INPUT_FILE})")
    parser.add_argument("-g", "--genome_file", default=DEFAULT_TARGET_GENOME_FASTA,
                        help=f"Input FASTA file for the target genome (default: {DEFAULT_TARGET_GENOME_FASTA})")
    parser.add_argument("-a", "--gene_table_file", default=DEFAULT_GENE_TABLE_FILE, # Changed argument name
                        help=f"Input gene annotation table file (CSV/TSV) for the target genome (default: {DEFAULT_GENE_TABLE_FILE})")
    parser.add_argument("-o", "--output_dir", default=DEFAULT_OUTPUT_DIR,
                        help=f"Directory to save output files (default: {DEFAULT_OUTPUT_DIR})")

    args = parser.parse_args()

    output_dir = args.output_dir
    sgrna_input_file = args.sgrna_file
    target_genome_fasta = args.genome_file
    gene_table_file = args.gene_table_file # Use new argument

    print("--- Starting sgRNA Exact Match and Gene Association Analysis (Table Input) ---")
    start_time = time.time()

    # --- Setup Output Directory ---
    if not os.path.exists(output_dir):
        print(f"Creating output directory: {output_dir}")
        os.makedirs(output_dir)
    else:
        print(f"Output directory already exists: {output_dir}")

    # --- Define Output File Paths ---
    final_csv_report = os.path.join(output_dir, 'sgrna_exact_match_gene_report.csv')
    summary_plot_file = os.path.join(output_dir, 'sgrna_exact_match_summary.png')

    # 1. Read Inputs
    sgrna_list = read_sgrnas(sgrna_input_file)
    genome, assumed_chromosome = read_genome(target_genome_fasta)
    genes_pr = read_gene_table(gene_table_file, assumed_chromosome) # Call new reader function

    # 2. Perform Sequence Search
    print(f"\nSearching for {len(sgrna_list)} unique sgRNAs in the genome...")
    match_results = []
    all_search_results = [] # Keep track of all sgRNAs for final merge
    found_count = 0
    not_found_count = 0

    search_start_time = time.time()
    for i, sgrna in enumerate(sgrna_list):
        if (i + 1) % 5000 == 0: print(f"  Searched {i+1}/{len(sgrna_list)} sgRNAs...") # Progress

        match_info = find_exact_matches(sgrna, genome)
        base_row = {'sgRNA_Sequence': sgrna, 'Match_Found': False}

        if match_info:
            found_count += 1
            base_row.update({
                'Match_Found': True,
                'Chromosome': match_info['chromosome'],
                'Start': match_info['start'], # 1-based
                'End': match_info['end'],     # 1-based
                'Strand_Match': match_info['strand']
            })
            match_results.append(base_row) # Add only if match found for gene association step
        else:
            not_found_count += 1

        all_search_results.append(base_row) # Store result for every sgRNA

    search_end_time = time.time()
    print(f"Search complete. Found matches for {found_count} sgRNAs, {not_found_count} not found.")
    print(f"Search duration: {search_end_time - search_start_time:.2f} seconds.")

    # Convert all search results to DataFrame for easy merging later
    all_results_df = pd.DataFrame(all_search_results)

    # 3. Perform Gene Association (only for sgRNAs with matches)
    print("\nAssociating matched sgRNAs with genes...")
    # Initialize columns in all_results_df to prevent KeyErrors later if no matches found
    gene_cols = ['Overlapping_Gene_Symbol', 'Overlapping_Locus_Tag', 'Overlapping_Product_Description',
                 'Nearest_Upstream_Gene_Symbol', 'Nearest_Upstream_Locus_Tag', 'Nearest_Upstream_Distance',
                 'Nearest_Downstream_Gene_Symbol', 'Nearest_Downstream_Locus_Tag', 'Nearest_Downstream_Distance']
    for col in gene_cols:
         if col not in all_results_df.columns: # Check if not already added
            all_results_df[col] = pd.NA # Use pandas NA marker

    if not match_results:
        print("No sgRNA matches found, skipping gene association.")
    else:
        matches_df = pd.DataFrame(match_results)

        # Prepare sgRNA matches for PyRanges (0-based start)
        sgrna_pr_df = matches_df[['Chromosome', 'Start', 'End', 'sgRNA_Sequence']].copy()
        sgrna_pr_df['Start'] = pd.to_numeric(sgrna_pr_df['Start'], errors='coerce')
        sgrna_pr_df['End'] = pd.to_numeric(sgrna_pr_df['End'], errors='coerce')
        sgrna_pr_df.dropna(subset=['Start', 'End'], inplace=True)
        if sgrna_pr_df.empty:
             print("WARNING: No valid coordinates found for matched sgRNAs. Skipping gene association.")
        else:
            sgrna_pr_df['Start'] = sgrna_pr_df['Start'] - 1
            sgrna_pr = pr.PyRanges(sgrna_pr_df)

            # --- Find Overlapping Genes ---
            print("  - Finding overlapping genes...")
            overlaps_pr = sgrna_pr.join(genes_pr, how="left", apply_strand_suffix=False)
            overlaps_df = overlaps_pr.df

            # Aggregate results for sgRNAs overlapping multiple genes
            overlaps_grouped = overlaps_df.groupby('sgRNA_Sequence').agg(
                Overlapping_Gene_Symbol=('Gene_Symbol', lambda x: ';'.join(x.dropna().astype(str).unique())),
                Overlapping_Locus_Tag=('Locus_Tag', lambda x: ';'.join(x.dropna().astype(str).unique())),
                Overlapping_Product_Description=('Product_Description', lambda x: ';'.join(x.dropna().astype(str).unique()))
            ).reset_index()

            # Keep only rows where an actual overlap occurred (check locus tag)
            valid_overlap_sgrnas = overlaps_df.dropna(subset=['Locus_Tag'])['sgRNA_Sequence'].unique()
            overlaps_final_df = overlaps_grouped[overlaps_grouped['sgRNA_Sequence'].isin(valid_overlap_sgrnas)].copy()
            overlaps_final_df.replace({'': pd.NA}, inplace=True) # Use pandas NA

            # --- Identify sgRNAs needing nearest search ---
            sgRNAs_with_match = set(matches_df['sgRNA_Sequence'])
            sgRNAs_with_overlaps = set(overlaps_final_df['sgRNA_Sequence'])
            non_overlapping_sgrnas = list(sgRNAs_with_match - sgRNAs_with_overlaps)

            # Prepare df for nearest search
            nearest_upstream_final_df = pd.DataFrame()
            nearest_downstream_final_df = pd.DataFrame()

            if non_overlapping_sgrnas:
                print(f"  - Finding nearest genes for {len(non_overlapping_sgrnas)} non-overlapping sgRNAs...")
                non_overlapping_df = matches_df[matches_df['sgRNA_Sequence'].isin(non_overlapping_sgrnas)].copy()

                if not non_overlapping_df.empty:
                    non_overlapping_pr_df = non_overlapping_df[['Chromosome', 'Start', 'End', 'sgRNA_Sequence']].copy()
                    non_overlapping_pr_df['Start'] = pd.to_numeric(non_overlapping_pr_df['Start'], errors='coerce')
                    non_overlapping_pr_df['End'] = pd.to_numeric(non_overlapping_pr_df['End'], errors='coerce')
                    non_overlapping_pr_df.dropna(subset=['Start', 'End'], inplace=True)

                    if not non_overlapping_pr_df.empty:
                        non_overlapping_pr_df['Start'] = non_overlapping_pr_df['Start'] - 1
                        non_overlapping_pr = pr.PyRanges(non_overlapping_pr_df)

                        # --- Find ONE Nearest (Upstream OR Downstream) ---
                        print("    - Finding single nearest gene (position will determine up/downstream)...")
                        # Run nearest once without direction constraint
                        # Nearest returns distance >= 0 or -1 if none found on chromosome
                        nearest_pr = non_overlapping_pr.nearest(genes_pr, suffix="_gene", overlap=False)
                        nearest_df = nearest_pr.df

                        # Filter out results where no nearest was found on chromosome (Distance == -1)
                        nearest_df_filtered = nearest_df[nearest_df['Distance'] >= 0].copy()

                        # Initialize empty DataFrames for results
                        nearest_upstream_final_df = pd.DataFrame()
                        nearest_downstream_final_df = pd.DataFrame()

                        if not nearest_df_filtered.empty:
                            # Determine Upstream vs Downstream based on coordinates
                            # Upstream: gene ends before sgRNA starts (End_gene < Start)
                            # Downstream: gene starts after sgRNA ends (Start_gene > End)
                            # Note: Coordinates are 0-based from PyRanges df
                            is_upstream = nearest_df_filtered['End_gene'] < nearest_df_filtered['Start']
                            is_downstream = nearest_df_filtered['Start_gene'] > nearest_df_filtered['End']

                            # Separate into upstream and downstream results
                            nearest_up_df = nearest_df_filtered[is_upstream].copy()
                            nearest_down_df = nearest_df_filtered[is_downstream].copy()

                            # Aggregate Upstream Results - *** ADD EMPTY CHECK ***
                            if not nearest_up_df.empty:
                                nearest_upstream_final_df = nearest_up_df.groupby('sgRNA_Sequence').agg(
                                    Nearest_Upstream_Gene_Symbol=('Gene_Symbol', lambda x: ';'.join(x.dropna().astype(str).unique())),
                                    Nearest_Upstream_Locus_Tag=('Locus_Tag', lambda x: ';'.join(x.dropna().astype(str).unique())),
                                    Nearest_Upstream_Distance=('Distance', 'min')
                                ).reset_index()
                                nearest_upstream_final_df.replace({'': pd.NA}, inplace=True)
                            else:
                                print("    - No valid nearest upstream genes found.")


                            # Aggregate Downstream Results - *** ADD EMPTY CHECK ***
                            if not nearest_down_df.empty:
                                nearest_downstream_final_df = nearest_down_df.groupby('sgRNA_Sequence').agg(
                                    Nearest_Downstream_Gene_Symbol=('Gene_Symbol', lambda x: ';'.join(x.dropna().astype(str).unique())),
                                    Nearest_Downstream_Locus_Tag=('Locus_Tag', lambda x: ';'.join(x.dropna().astype(str).unique())),
                                    Nearest_Downstream_Distance=('Distance', 'min')
                                ).reset_index()
                                nearest_downstream_final_df.replace({'': pd.NA}, inplace=True)
                            else:
                                print("    - No valid nearest downstream genes found.")
                        else:
                            print("    - No nearest genes found (Distance >= 0) for non-overlapping sgRNAs.")

                    else:
                         print("    - No valid coordinates for non-overlapping sgRNAs after filtering.")
                else:
                     print("    - Non-overlapping dataframe empty.") # Should not happen if list wasn't empty
            else:
                print("  - No sgRNAs required nearest gene search.")


            # --- Combine All Results ---
            print("  - Combining overlap, upstream, and downstream results...")
            # Start with the full list of sgRNAs and their match status
            final_results_df = all_results_df.copy()

            # Merge overlap info
            if not overlaps_final_df.empty:
                 final_results_df = pd.merge(final_results_df, overlaps_final_df, on='sgRNA_Sequence', how='left')

            # Merge nearest upstream info
            if not nearest_upstream_final_df.empty:
                 final_results_df = pd.merge(final_results_df, nearest_upstream_final_df, on='sgRNA_Sequence', how='left')

            # Merge nearest downstream info
            if not nearest_downstream_final_df.empty:
                 final_results_df = pd.merge(final_results_df, nearest_downstream_final_df, on='sgRNA_Sequence', how='left')

            all_results_df = final_results_df # Assign back to the main df variable


    # 4. Format Final Report
    print(f"\nGenerating final CSV report: {final_csv_report}")
    try:
        # Convert potentially numeric columns (that might have NAs) to nullable Ints first
        numeric_cols_to_convert = ['Start', 'End', 'Nearest_Upstream_Distance', 'Nearest_Downstream_Distance']
        for col in numeric_cols_to_convert:
            if col in all_results_df.columns:
                # Use pd.NA for consistency
                all_results_df[col] = pd.to_numeric(all_results_df[col], errors='coerce').astype('Int64')

        # Define dictionary for filling NAs in string columns
        string_fill_dict = {
            'Chromosome': 'N/A', 'Strand_Match': 'N/A',
            'Overlapping_Gene_Symbol': 'N/A', 'Overlapping_Locus_Tag': 'N/A', 'Overlapping_Product_Description': 'N/A',
            'Nearest_Upstream_Gene_Symbol': 'N/A', 'Nearest_Upstream_Locus_Tag': 'N/A',
            'Nearest_Downstream_Gene_Symbol': 'N/A', 'Nearest_Downstream_Locus_Tag': 'N/A',
             # Note: Int64 columns (Start, End, Distances) are handled separately below
        }
        # Apply fillna only for columns that exist in the dataframe AND are in our string fill dict
        cols_to_fill_strings = {k: v for k, v in string_fill_dict.items() if k in all_results_df.columns}
        # Use pd.NA marker for Int64 columns before converting to string
        na_marker = pd.NA
        all_results_df.fillna(cols_to_fill_strings, inplace=True)


        # Convert Int64 columns to string *now*, replacing pd.NA marker explicitly
        int_cols_to_str = ['Start', 'End', 'Nearest_Upstream_Distance', 'Nearest_Downstream_Distance']
        for col in int_cols_to_str:
             if col in all_results_df.columns:
                  # Convert to string, replacing pandas' NA representation (<NA>) with 'N/A'
                  all_results_df[col] = all_results_df[col].astype(str).replace('<NA>', 'N/A')


        # Define final column order - include the new gene info columns
        report_columns = [
            'sgRNA_Sequence', 'Match_Found', 'Chromosome', 'Start', 'End', 'Strand_Match',
            'Overlapping_Gene_Symbol', 'Overlapping_Locus_Tag', 'Overlapping_Product_Description',
            'Nearest_Upstream_Gene_Symbol', 'Nearest_Upstream_Locus_Tag', 'Nearest_Upstream_Distance',
            'Nearest_Downstream_Gene_Symbol', 'Nearest_Downstream_Locus_Tag', 'Nearest_Downstream_Distance'
        ]
        final_columns = [col for col in report_columns if col in all_results_df.columns]
        # Add any columns present in df but not explicitly listed (shouldn't happen ideally)
        final_columns.extend([col for col in all_results_df.columns if col not in final_columns])

        all_results_df = all_results_df[final_columns]

        all_results_df.to_csv(final_csv_report, index=False, na_rep='N/A') # na_rep acts as a final backup
        print("CSV report generated successfully.")
    except Exception as e:
        print(f"ERROR: Failed to generate CSV report: {e}")
        import traceback
        traceback.print_exc()

    # 5. Generate Summary Plot (same as before)
    print(f"\nGenerating summary plot: {summary_plot_file}")
    # (Plotting code remains unchanged, shows Found vs Not Found)
    try:
        plt.figure(figsize=(6, 5))
        categories = ['Found', 'Not Found']
        counts = [found_count, not_found_count]
        bars = plt.bar(categories, counts, color=['#4CAF50', '#F44336'])
        plt.xlabel("Exact Match Status")
        plt.ylabel("Number of Unique sgRNAs")
        plt.title("Exact Match Summary of sgRNAs in Target Genome")
        plt.ylim(bottom=0)
        plt.grid(axis='y', linestyle='--', alpha=0.7)
        for bar in bars:
            yval = bar.get_height()
            if yval > 0: plt.text(bar.get_x() + bar.get_width()/2.0, yval, int(yval), va='bottom', ha='center', fontsize=10)
        plt.tight_layout()
        plt.savefig(summary_plot_file, dpi=300)
        print("Summary plot generated successfully.")
    except Exception as e:
        print(f"ERROR: Failed to generate summary plot: {e}")


    end_time = time.time()
    print(f"\n--- Analysis Complete ---")
    print(f"Total execution time: {end_time - start_time:.2f} seconds.")
    print(f"Results saved in directory: {output_dir}")

if __name__ == "__main__":
    # --- Example sgRNAs --- (Same as before)
    example_sgrnas_text = """
    GGAAGGGGTGGCTTCGAGCGT
    GATCGCCGGAACGTTCACACA
    CCGGAGGTTCGCCGGGAAGG
    ACTTCCGTGCCATCAATAAA
    GCCACAGCCACATTCATTCT
    CACATTCATTCTGGGCTTTA
    GAGATTAATTAGCGACTGTT
    ATTTTACCTGAACCATAAT
    GATACGAAACCATTGTTGAC
    GCTGGGAGTTGGTGCTGGAT
    TCGCTTCCGGCGTGGCAAAG
    TGATCTGCTCGGCCTGTTCC
    GAAGTTGTAGAGACGCACAC
    TTCCATTATCGAACGACAAT
    GGAACATCCAGATGGAAATA
    CAGCCGCATCGCGCCGGCA
    TGCCAGCTTTTCGGCATTTG
    GCAAAAGCCTGCGGGAGAAA
    CAGCGTACCGAAGCGCAAGC
    INTERGENIC_SEQUENCE_1
    INTERGENIC_SEQUENCE_2
    NOT_IN_GENOME_SEQUENCE
    INVALIDSEQUENCE###INVALID
    """

    # --- Dummy Gene Table Content ---
    # Using tab separation (.tsv) as default
    dummy_table_text = """Gene Name\tAccession-1\tLeft-End-Position\tRight-End-Position\tProduct
geneA\tBW_001\t150\t250\tFirst dummy gene
geneB\tBW_002\t300\t400\tSecond dummy gene, opposite strand conceptually
#geneX\tBW_XXX\tBAD\tBAD\tCommented out gene with bad coords
geneE\tBW_005\t\t600\tGene with blank start
geneF\tBW_006\t700\t\tGene with blank end
geneG\tBW_007\t800.5\t900\tGene with float start
geneH\tBW_008\t1000\t950\tGene with start > end
geneI\tBW_009\t0\t50\tGene with 0 start (becomes -1)
geneC\tBW_003\t5070\t5170\tThird dummy gene on chr2
geneD\tBW_004\t5220\t5320\tFourth dummy gene on chr2
""" # Added more invalid/edge cases

    # Create dummy input files if they don't exist
    if not os.path.exists(DEFAULT_SGRNA_INPUT_FILE):
        print(f"Creating dummy input file: {DEFAULT_SGRNA_INPUT_FILE}")
        with open(DEFAULT_SGRNA_INPUT_FILE, "w") as f: f.write(example_sgrnas_text)

    # Create dummy genome (adjust sequence names if necessary)
    # Make sure chromosome names match what read_gene_table will use (e.g., 'Dummy_Chr1')
    if not os.path.exists(DEFAULT_TARGET_GENOME_FASTA):
         print(f"Creating dummy target genome file: {DEFAULT_TARGET_GENOME_FASTA}")
         with open(DEFAULT_TARGET_GENOME_FASTA, "w") as f:
              f.write(">Dummy_Chr1\n") # MUST match chromosome name assumption
              # ACGCTCGAAGCCACCCCTTCC -> revcomp GGAAGG... (coords 121-141 approx) -> Nearest upstream: None, Nearest downstream: geneA (dist ~ 8)
              # TGTGTGAACGTTCCGGCGATC -> revcomp GATCGC... (coords 182-202 approx) -> Nearest upstream: geneA (dist ~ 31), Nearest downstream: geneB (dist ~ 97)
              # INTERGENIC_SEQUENCE_1 (coords 162-181 approx) -> Nearest upstream: geneA (dist ~ 11), Nearest downstream: geneB (dist ~ 118)
              dummy_genome_part1 = "N"*120 + "ACGCTCGAAGCCACCCCTTCC" + "N"*20 + \
                                   "INTERGENIC_SEQUENCE_1" + "N"*20 + \
                                   "TGTGTGAACGTTCCGGCGATC" + "N"*4500 # Pad to ensure coordinates for Chr2 are distinct
              f.write(dummy_genome_part1 + "\n")
              f.write(">Dummy_Chr2\n") # Script currently assumes single chrom from FASTA, gene association won't use this unless modified
              # GCCACAGCCACATTCATTCT -> (coords 5091-5110) -> Overlaps geneC (5070-5170) YES
              # TAAAGCCCAGAATGAATGTG -> revcomp CACATT... (coords 5161-5180) -> Nearest upstream: geneC (dist ~ 0), Nearest downstream: geneD (dist ~ 39)
              # INTERGENIC_SEQUENCE_2 (coords 5231-5250) -> Nearest upstream: geneD (dist ~ 10), Nearest downstream: None
              dummy_genome_part2 = "N"*5090 + "GCCACAGCCACATTCATTCT" + "N"*50 + \
                                   "TAAAGCCCAGAATGAATGTG" + "N"*50 + \
                                   "INTERGENIC_SEQUENCE_2" + "N"*50
              f.write(dummy_genome_part2 + "\n")

    # Create dummy gene table file
    if not os.path.exists(DEFAULT_GENE_TABLE_FILE):
        print(f"Creating dummy gene table file: {DEFAULT_GENE_TABLE_FILE}")
        with open(DEFAULT_GENE_TABLE_FILE, "w") as f: f.write(dummy_table_text)

    main()