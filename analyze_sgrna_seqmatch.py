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
                    # Basic validation: Allow only A, T, C, G
                    if all(c in 'ATCG' for c in seq): # Stricter: only ATCG
                        sequences.add(seq)
                        valid_count += 1
                    elif all(c in 'ATCGN' for c in seq):
                         print(f"WARNING: Skipping sequence with 'N' (line {i+1}): {seq}")
                         skipped_count += 1
                    else:
                        print(f"WARNING: Skipping invalid sequence (line {i+1}, non-ATCGN): {seq}")
                        skipped_count += 1
    except FileNotFoundError:
        print(f"ERROR: sgRNA input file not found: {filepath}")
        sys.exit(1)

    unique_sequences = sorted(list(sequences))
    print(f"Read {valid_count} valid sequences (ATCG only), skipped {skipped_count} invalid or N-containing sequences.")
    print(f"Found {len(unique_sequences)} unique valid sgRNA sequences for searching.")
    if not unique_sequences:
        print("ERROR: No valid ATCG-only sgRNA sequences to process.")
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
            genome[record.id] = str(record.seq).upper() # Ensure uppercase
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
        for col in coord_cols_original:
             genes_df = genes_df.copy()
             genes_df[col] = pd.to_numeric(genes_df[col], errors='coerce')

        original_rows = len(genes_df)
        genes_df.dropna(subset=coord_cols_original, inplace=True)
        filtered_rows = len(genes_df)

        if original_rows > filtered_rows:
             print(f"Filtered out {original_rows - filtered_rows} rows due to invalid/blank coordinate values.")

        if genes_df.empty:
             print("ERROR: No valid gene entries remaining after filtering coordinates.")
             sys.exit(1)
        else:
            print(f"{filtered_rows} rows remaining after coordinate filtering.")

        # --- Now safe to convert coordinates to integers ---
        try:
            for col in coord_cols_original:
                 if pd.api.types.is_float_dtype(genes_df[col]):
                     if not genes_df[col].eq(genes_df[col].round()).all():
                         print(f"WARNING: Column '{col}' contains non-integer numbers after filtering. They will be truncated to integers.")
                 genes_df[col] = genes_df[col].astype(int)
        except ValueError as e:
             print(f"ERROR: Could not convert coordinate columns to integer after filtering: {e}")
             sys.exit(1)


        # --- Proceed with renaming ---
        rename_map = {
            "Gene Name": "Gene_Symbol",
            "Accession-1": "Locus_Tag",
            "Left-End-Position": "Start",
            "Right-End-Position": "End",
            "Product": "Product_Description"
        }
        genes_df_renamed = genes_df[original_required_cols].copy()
        genes_df_renamed.rename(columns=rename_map, inplace=True)


        # Assign expected chromosome name
        print(f"Assigning chromosome name '{expected_chromosome}' to all genes.")
        genes_df_renamed['Chromosome'] = expected_chromosome

        # --- Coordinate Conversion (1-based inclusive -> 0-based start, 1-based end) ---
        genes_df_renamed['Start'] = genes_df_renamed['Start'] - 1
        genes_df_renamed['Strand'] = '+' # Or '.'

        # Select final columns for PyRanges object
        pr_cols = ['Chromosome', 'Start', 'End', 'Strand', 'Gene_Symbol', 'Locus_Tag', 'Product_Description']
        genes_final_df = genes_df_renamed[pr_cols].copy()
        genes_final_df = genes_final_df[genes_final_df['Start'] >= 0] # Filter out negative Start coords
        genes_final_df = genes_final_df[genes_final_df['Start'] < genes_final_df['End']] # Ensure Start < End

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
    """
    Finds the first exact match (fwd or revcomp) that has a valid NGG PAM
    immediately downstream on the non-target strand.
    Returns dict with match info or None.
    """
    sgrna_len = len(sgrna_seq)
    if sgrna_len == 0: return None
    try:
        rc_sgrna_seq = str(Seq(sgrna_seq).reverse_complement())
    except Exception as e:
         print(f"ERROR: Could not reverse complement sgRNA '{sgrna_seq}': {e}")
         return None

    for chrom_id, chrom_seq in genome_sequences.items():
        chrom_len = len(chrom_seq)

        # --- Search Forward Strand ('+' match) ---
        current_pos = 0
        while True:
            index_fwd = chrom_seq.find(sgrna_seq, current_pos)
            if index_fwd == -1: break
            pam_start_0based = index_fwd + sgrna_len
            pam_end_0based = pam_start_0based + 3
            pam_sequence = "N/A"
            pam_found = False
            if pam_end_0based <= chrom_len:
                pam_candidate = chrom_seq[pam_start_0based:pam_end_0based]
                pam_sequence = pam_candidate
                if len(pam_candidate) == 3 and pam_candidate[1] == 'G' and pam_candidate[2] == 'G' and all(c in 'ATCG' for c in pam_candidate):
                    pam_found = True
            else: pam_sequence = "Out_of_Bounds"
            if pam_found:
                return {'chromosome': chrom_id, 'start': index_fwd + 1, 'end': index_fwd + sgrna_len,
                        'strand': '+', 'pam_sequence': pam_sequence, 'pam_found': True}
            current_pos = index_fwd + 1

        # --- Search Reverse Strand ('-' match) ---
        current_pos = 0
        while True:
            index_rev = chrom_seq.find(rc_sgrna_seq, current_pos)
            if index_rev == -1: break
            pam_start_0based = index_rev - 3
            pam_end_0based = index_rev
            pam_sequence = "N/A"
            pam_found = False
            if pam_start_0based >= 0:
                pam_candidate = chrom_seq[pam_start_0based:pam_end_0based]
                pam_sequence = pam_candidate
                if len(pam_candidate) == 3 and pam_candidate[0] == 'C' and pam_candidate[1] == 'C' and all(c in 'ATCG' for c in pam_candidate):
                     pam_found = True
            else: pam_sequence = "Out_of_Bounds"
            if pam_found:
                return {'chromosome': chrom_id, 'start': index_rev + 1, 'end': index_rev + sgrna_len,
                        'strand': '-', 'pam_sequence': pam_sequence, 'pam_found': True}
            current_pos = index_rev + 1
    return None


# --- Main Execution ---
def main():
    parser = argparse.ArgumentParser(description="Find exact matches of sgRNA sequences in a target genome, check for downstream NGG PAM, and associate with overlapping genes using a custom table.")
    parser.add_argument("-s", "--sgrna_file", required=True,
                        help=f"Input file containing sgRNA sequences (ATCG only), one per line")
    parser.add_argument("-g", "--genome_file", required=True,
                        help=f"Input FASTA file for the target genome")
    parser.add_argument("-a", "--gene_table_file", required=True,
                        help=f"Input gene annotation table file (CSV/TSV) for the target genome")
    parser.add_argument("-o", "--output_dir", required=True,
                        help=f"Directory to save output files")

    args = parser.parse_args()

    output_dir = args.output_dir
    sgrna_input_file = args.sgrna_file
    target_genome_fasta = args.genome_file
    gene_table_file = args.gene_table_file

    print("--- Starting sgRNA Exact Match, PAM Check, and Overlapping Gene Association Analysis ---")
    start_time = time.time()

    # --- Setup Output Directory ---
    if not os.path.exists(output_dir):
        print(f"Creating output directory: {output_dir}")
        os.makedirs(output_dir)
    else:
        print(f"Output directory already exists: {output_dir}")

    # --- Define Output File Paths ---
    final_csv_report = os.path.join(output_dir, 'sgrna_overlap_gene_pam_report.csv')
    summary_plot_file = os.path.join(output_dir, 'sgrna_exact_match_summary.png')

    # 1. Read Inputs
    sgrna_list = read_sgrnas(sgrna_input_file)
    genome, assumed_chromosome = read_genome(target_genome_fasta)
    genes_pr = read_gene_table(gene_table_file, assumed_chromosome)

    # 2. Perform Sequence Search and PAM Check
    print(f"\nSearching for {len(sgrna_list)} unique sgRNAs in the genome and checking for valid PAM...")
    match_results_for_genes = []
    all_search_results = []
    found_count = 0
    pam_found_count = 0
    not_found_count = 0

    search_start_time = time.time()
    for i, sgrna in enumerate(sgrna_list):
        if (i + 1) % 5000 == 0: print(f"  Processed {i+1}/{len(sgrna_list)} sgRNAs...")

        match_info = find_exact_matches(sgrna, genome)

        base_row = {
            'sgRNA_Sequence': sgrna, 'Match_Found': False, 'Chromosome': pd.NA,
            'Start': pd.NA, 'End': pd.NA, 'Strand_Match': pd.NA,
            'PAM_Sequence': pd.NA, 'PAM_Found': False
        }

        if match_info:
            pam_found_count += 1
            base_row.update({
                'Match_Found': True, 'Chromosome': match_info['chromosome'],
                'Start': match_info['start'], 'End': match_info['end'],
                'Strand_Match': match_info['strand'],
                'PAM_Sequence': match_info['pam_sequence'], 'PAM_Found': True
            })
            match_results_for_genes.append({
                'sgRNA_Sequence': sgrna, 'Chromosome': match_info['chromosome'],
                'Start': match_info['start'], 'End': match_info['end']
            })
            found_count +=1
        else:
            temp_found = False
            sgrna_rc = str(Seq(sgrna).reverse_complement())
            for chrom_seq in genome.values():
                if chrom_seq.find(sgrna) != -1 or chrom_seq.find(sgrna_rc) != -1:
                    temp_found = True
                    break
            if temp_found:
                base_row['Match_Found'] = True
                found_count += 1
            else:
                 not_found_count += 1

        all_search_results.append(base_row)

    search_end_time = time.time()
    print(f"Search complete.")
    print(f"  Found at least one match site for {found_count} sgRNAs.")
    print(f"  Reported match site has valid PAM for {pam_found_count} sgRNAs.")
    print(f"  {not_found_count} sgRNAs not found in genome.")
    print(f"Search and PAM check duration: {search_end_time - search_start_time:.2f} seconds.")

    all_results_df = pd.DataFrame(all_search_results)

    # 3. Perform Overlapping Gene Association
    print("\nAssociating matched sgRNAs (with valid PAM) with overlapping genes...")
    # Define expected gene columns JUST for later reference in formatting
    gene_cols = ['Overlapping_Gene_Symbol', 'Overlapping_Locus_Tag', 'Overlapping_Product_Description']

    if not match_results_for_genes:
        print("No sgRNA matches with valid PAM found, skipping gene association.")
        # Ensure gene columns exist in the final df even if no overlaps are found/checked
        for col in gene_cols:
             if col not in all_results_df.columns:
                  all_results_df[col] = pd.NA
    else:
        matches_df_for_genes = pd.DataFrame(match_results_for_genes)
        sgrna_pr_df = matches_df_for_genes[['Chromosome', 'Start', 'End', 'sgRNA_Sequence']].copy()
        sgrna_pr_df['Start'] = pd.to_numeric(sgrna_pr_df['Start'], errors='coerce').astype('Int64') - 1
        sgrna_pr_df['End'] = pd.to_numeric(sgrna_pr_df['End'], errors='coerce').astype('Int64')
        sgrna_pr_df.dropna(subset=['Start', 'End'], inplace=True)

        overlaps_final_df = pd.DataFrame() # Initialize empty dataframe for overlap results

        if sgrna_pr_df.empty:
             print("WARNING: No valid coordinates found for matched sgRNAs with PAM. Skipping gene association.")
        else:
            sgrna_pr_df = sgrna_pr_df[sgrna_pr_df['Start'] >= 0]
            if sgrna_pr_df.empty:
                print("WARNING: No valid non-negative coordinates after 0-based conversion. Skipping gene association.")
            else:
                sgrna_pr = pr.PyRanges(sgrna_pr_df)
                print("  - Finding overlapping genes...")
                overlaps_pr = sgrna_pr.join(genes_pr, how="left", apply_strand_suffix=False)
                overlaps_df = overlaps_pr.df

                # Check if overlaps_df is empty or has no valid overlaps
                if overlaps_df.empty or overlaps_df['Locus_Tag'].isnull().all():
                     print("  - No overlaps found.")
                else:
                    overlaps_grouped = overlaps_df.groupby('sgRNA_Sequence').agg(
                        Overlapping_Gene_Symbol=('Gene_Symbol', lambda x: ';'.join(x.dropna().astype(str).unique())),
                        Overlapping_Locus_Tag=('Locus_Tag', lambda x: ';'.join(x.dropna().astype(str).unique())),
                        Overlapping_Product_Description=('Product_Description', lambda x: ';'.join(x.dropna().astype(str).unique()))
                    ).reset_index()

                    valid_overlap_sgrnas = overlaps_df.dropna(subset=['Locus_Tag'])['sgRNA_Sequence'].unique()
                    overlaps_final_df = overlaps_grouped[overlaps_grouped['sgRNA_Sequence'].isin(valid_overlap_sgrnas)].copy()
                    overlaps_final_df.replace({'': pd.NA}, inplace=True) # Replace empty strings potentially generated by join with NA
                    print(f"  - Found overlaps for {len(overlaps_final_df)} sgRNAs.")

        # --- Combine Overlap Results ---
        print("  - Combining overlap results with main table...")
        # Perform the merge. This adds the gene columns from overlaps_final_df.
        # If overlaps_final_df is empty, the merge still works but adds no data.
        # Columns existing only in all_results_df are kept.
        # Columns existing only in overlaps_final_df are added with NaNs where no match.
        # Shared column 'sgRNA_Sequence' is used for joining.
        # Crucially, gene_cols are NOT pre-initialized in all_results_df anymore.
        all_results_df = pd.merge(all_results_df, overlaps_final_df, on='sgRNA_Sequence', how='left')

        # Ensure gene columns exist after merge, even if overlaps_final_df was empty
        for col in gene_cols:
             if col not in all_results_df.columns:
                  all_results_df[col] = pd.NA


    # 4. Format Final Report
    print(f"\nGenerating final CSV report: {final_csv_report}")
    try:
        numeric_cols_to_convert = ['Start', 'End']
        for col in numeric_cols_to_convert:
            if col in all_results_df.columns:
                all_results_df[col] = pd.to_numeric(all_results_df[col], errors='coerce').astype('Int64')

        # Use the gene_cols list defined earlier for filling NAs
        string_fill_dict = {
            'Chromosome': 'N/A', 'Strand_Match': 'N/A', 'PAM_Sequence': 'N/A',
            **{col: 'N/A' for col in gene_cols} # Dynamically add gene columns to fill dict
        }
        cols_to_fill_strings = {k: v for k, v in string_fill_dict.items() if k in all_results_df.columns}
        all_results_df.fillna(cols_to_fill_strings, inplace=True)

        if 'Match_Found' in all_results_df.columns:
             all_results_df['Match_Found'] = all_results_df['Match_Found'].fillna(False).astype(str)
        if 'PAM_Found' in all_results_df.columns:
             all_results_df['PAM_Found'] = all_results_df['PAM_Found'].fillna(False).astype(str)

        int_cols_to_str = ['Start', 'End']
        for col in int_cols_to_str:
             if col in all_results_df.columns:
                  all_results_df[col] = all_results_df[col].astype(str).replace('<NA>', 'N/A')

        # Define final column order using the gene_cols list
        report_columns = [
            'sgRNA_Sequence', 'Match_Found', 'PAM_Found',
            'Chromosome', 'Start', 'End', 'Strand_Match', 'PAM_Sequence',
            *gene_cols # Unpack the gene column names here
        ]
        # Ensure only existing columns are selected and preserve order
        final_columns = [col for col in report_columns if col in all_results_df.columns]
        # Add any other columns that might exist but weren't explicitly listed
        final_columns.extend([col for col in all_results_df.columns if col not in final_columns])

        all_results_df = all_results_df[final_columns]

        all_results_df.to_csv(final_csv_report, index=False, na_rep='N/A')
        print("CSV report generated successfully.")
    except Exception as e:
        print(f"ERROR: Failed to generate CSV report: {e}")
        import traceback
        traceback.print_exc()

    # 5. Generate Summary Plot
    print(f"\nGenerating summary plot: {summary_plot_file}")
    try:
        plt.figure(figsize=(6, 5))
        categories = ['Found', 'Not Found']
        counts = [found_count, not_found_count]
        bars = plt.bar(categories, counts, color=['#4CAF50', '#F44336'])
        plt.xlabel("Exact Match Site Status")
        plt.ylabel("Number of Unique sgRNAs")
        plt.title("Exact Match Site Summary of sgRNAs in Target Genome")
        plt.ylim(bottom=0)
        for bar in bars:
            yval = bar.get_height()
            if yval > 0: plt.text(bar.get_x() + bar.get_width()/2.0, yval, int(yval), va='bottom', ha='center', fontsize=10)
        plt.grid(axis='y', linestyle='--', alpha=0.7)
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
    main()