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
# import re # No longer needed for simple PAM check

# --- Helper Functions ---
# (read_sgrnas, read_genome, read_gene_table remain the same)
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
    """
    Finds the first exact match (fwd or revcomp) that has a valid NGG PAM
    immediately downstream on the non-target strand.

    Follows the logic similar to the example:
     - Forward match ('+'): Genome has [sgRNA_seq][NGG]
     - Reverse match ('-'): Genome has [CCN][sgRNA_rc]

    Returns:
        dict: Information about the first valid match including PAM sequence
              and validity, or None if no valid match found.
              Keys: 'chromosome', 'start', 'end', 'strand', 'pam_sequence', 'pam_found'
              'pam_sequence' is the sequence found in the genome (NGG or CCN).
    """
    sgrna_len = len(sgrna_seq)
    if sgrna_len == 0: return None
    try:
        # Ensure only ATCG - Ns were filtered earlier
        rc_sgrna_seq = str(Seq(sgrna_seq).reverse_complement())
    except Exception as e:
         # This shouldn't happen if input validation worked, but good practice
         print(f"ERROR: Could not reverse complement sgRNA '{sgrna_seq}': {e}")
         return None

    for chrom_id, chrom_seq in genome_sequences.items():
        chrom_len = len(chrom_seq)

        # --- Search Forward Strand ('+' match: sgRNA binds '-' strand, PAM is NGG on '+' strand) ---
        # Look for sgRNA_seq followed by NGG
        current_pos = 0
        while True:
            # Find the next occurrence of the sgrna sequence
            index_fwd = chrom_seq.find(sgrna_seq, current_pos)
            if index_fwd == -1:
                break # No more occurrences found in this chromosome

            # Check for NGG PAM immediately following the match
            pam_start_0based = index_fwd + sgrna_len
            pam_end_0based = pam_start_0based + 3
            pam_sequence = "N/A"
            pam_found = False

            if pam_end_0based <= chrom_len: # Check boundary
                pam_candidate = chrom_seq[pam_start_0based:pam_end_0based]
                pam_sequence = pam_candidate # Store candidate for reporting

                # Validate NGG PAM (N = A, T, C, or G)
                if len(pam_candidate) == 3 and pam_candidate[1] == 'G' and pam_candidate[2] == 'G' and all(c in 'ATCG' for c in pam_candidate):
                    pam_found = True
            else:
                pam_sequence = "Out_of_Bounds" # Indicate boundary issue

            # If PAM is valid, we found our first valid hit
            if pam_found:
                start_1based = index_fwd + 1
                end_1based = index_fwd + sgrna_len
                return {
                    'chromosome': chrom_id,
                    'start': start_1based,
                    'end': end_1based,
                    'strand': '+', # Indicates sgRNA sequence matches the '+' strand
                    'pam_sequence': pam_sequence, # The NGG sequence
                    'pam_found': True
                }

            # If PAM was not valid, continue searching from the position after the current match
            current_pos = index_fwd + 1

        # --- Search Reverse Strand ('-' match: sgRNA binds '+' strand, PAM is NGG on '-' strand / CCN on '+' strand) ---
        # Look for CCN followed by sgRNA_rc
        current_pos = 0
        while True:
            # Find the next occurrence of the reverse complement sequence
            index_rev = chrom_seq.find(rc_sgrna_seq, current_pos)
            if index_rev == -1:
                break # No more occurrences found in this chromosome

            # Check for CCN PAM immediately preceding the rc_sgrna_seq match
            pam_start_0based = index_rev - 3
            pam_end_0based = index_rev
            pam_sequence = "N/A"
            pam_found = False

            if pam_start_0based >= 0: # Check boundary
                pam_candidate = chrom_seq[pam_start_0based:pam_end_0based]
                pam_sequence = pam_candidate # Store candidate for reporting

                # Validate CCN PAM (N = A, T, C, or G) - Reverse complement of NGG
                if len(pam_candidate) == 3 and pam_candidate[0] == 'C' and pam_candidate[1] == 'C' and all(c in 'ATCG' for c in pam_candidate):
                     pam_found = True
            else:
                pam_sequence = "Out_of_Bounds" # Indicate boundary issue

            # If PAM is valid, we found our first valid hit
            if pam_found:
                start_1based = index_rev + 1
                end_1based = index_rev + sgrna_len
                return {
                    'chromosome': chrom_id,
                    'start': start_1based, # Start of the rc_sgrna_seq match
                    'end': end_1based,     # End of the rc_sgrna_seq match
                    'strand': '-', # Indicates sgRNA binds the '+' strand (rc match found)
                    'pam_sequence': pam_sequence, # The CCN sequence
                    'pam_found': True
                }

            # If PAM was not valid, continue searching from the position after the current rc match
            current_pos = index_rev + 1

    # If loops complete without returning a valid hit
    return None


# --- Main Execution ---
def main():
    parser = argparse.ArgumentParser(description="Find exact matches of sgRNA sequences in a target genome, check for downstream NGG PAM, and associate with overlapping/nearest genes using a custom table.")
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
    gene_table_file = args.gene_table_file # Use new argument

    print("--- Starting sgRNA Exact Match, PAM Check, and Gene Association Analysis (Table Input) ---")
    start_time = time.time()

    # --- Setup Output Directory ---
    if not os.path.exists(output_dir):
        print(f"Creating output directory: {output_dir}")
        os.makedirs(output_dir)
    else:
        print(f"Output directory already exists: {output_dir}")

    # --- Define Output File Paths ---
    final_csv_report = os.path.join(output_dir, 'sgrna_exact_match_gene_pam_report.csv') # Updated filename
    summary_plot_file = os.path.join(output_dir, 'sgrna_exact_match_summary.png')

    # 1. Read Inputs
    sgrna_list = read_sgrnas(sgrna_input_file)
    genome, assumed_chromosome = read_genome(target_genome_fasta)
    genes_pr = read_gene_table(gene_table_file, assumed_chromosome) # Call new reader function

    # 2. Perform Sequence Search and PAM Check
    print(f"\nSearching for {len(sgrna_list)} unique sgRNAs in the genome and checking for valid PAM...")
    match_results_for_genes = [] # Store location details ONLY for matched sgRNAs (PAM irrelevant here) for gene association
    all_search_results = [] # Keep track of all sgRNAs (match status, PAM status) for final merge
    found_count = 0 # Counts sgRNAs with at least one match site found (regardless of PAM)
    pam_found_count = 0 # Counts sgRNAs where the *first reported match* has a valid PAM
    not_found_count = 0 # Counts sgRNAs with NO match site found at all

    search_start_time = time.time()
    for i, sgrna in enumerate(sgrna_list):
        if (i + 1) % 5000 == 0: print(f"  Processed {i+1}/{len(sgrna_list)} sgRNAs...") # Progress

        # find_exact_matches now returns the first hit *with* a valid PAM, or None
        match_info = find_exact_matches(sgrna, genome)

        # Initialize base row with defaults for all sgRNAs
        base_row = {
            'sgRNA_Sequence': sgrna,
            'Match_Found': False, # Will be set to True if *any* match was potentially found (even without PAM, see below)
            'Chromosome': pd.NA,
            'Start': pd.NA,
            'End': pd.NA,
            'Strand_Match': pd.NA, # Strand of the reported match (if any)
            'PAM_Sequence': pd.NA, # PAM of the reported match (if any with PAM)
            'PAM_Found': False     # Specifically refers to the *reported* match
        }

        # Check if the specific function found a match *with* a PAM
        if match_info:
            pam_found_count += 1
            base_row.update({
                'Match_Found': True, # A valid site (match + PAM) was found
                'Chromosome': match_info['chromosome'],
                'Start': match_info['start'],
                'End': match_info['end'],
                'Strand_Match': match_info['strand'],
                'PAM_Sequence': match_info['pam_sequence'],
                'PAM_Found': match_info['pam_found'] # Should always be True if match_info is not None
            })
            # Add details to list used for gene association
            match_results_for_genes.append({
                'sgRNA_Sequence': sgrna,
                'Chromosome': match_info['chromosome'],
                'Start': match_info['start'],
                'End': match_info['end']
            })
            found_count +=1 # Count this sgRNA as having found a site

        else:
            # If find_exact_matches returned None, it means no site *with* a valid PAM was found.
            # We still need to know if the sgRNA sequence *existed* at all in the genome
            # for the "Found" vs "Not Found" summary plot.
            # Let's do a quick check *without* PAM requirement just for counting.
            temp_found = False
            sgrna_rc = str(Seq(sgrna).reverse_complement())
            for chrom_seq in genome.values():
                if chrom_seq.find(sgrna) != -1 or chrom_seq.find(sgrna_rc) != -1:
                    temp_found = True
                    break
            if temp_found:
                # Match site exists, but none had a valid PAM (or the first one didn't)
                base_row['Match_Found'] = True # Update status
                found_count += 1
                # Keep PAM_Found as False and other details as NA/default
            else:
                # No match site found at all
                 not_found_count += 1
                 # base_row defaults are already correct (Match_Found=False, etc.)

        all_search_results.append(base_row) # Store result for every sgRNA

    search_end_time = time.time()
    print(f"Search complete.")
    # Note: found_count now reflects sgRNAs with *any* match site.
    # pam_found_count reflects sgRNAs where the first reported site had a valid PAM.
    print(f"  Found at least one match site for {found_count} sgRNAs.")
    print(f"  Reported match site has valid PAM for {pam_found_count} sgRNAs.") # This count comes from find_exact_matches succeeding
    print(f"  {not_found_count} sgRNAs not found in genome.")
    print(f"Search and PAM check duration: {search_end_time - search_start_time:.2f} seconds.")


    # Convert all search results to DataFrame for easy merging later
    # This df now contains the result for the *first valid hit with PAM* found,
    # or indicates if a match was found but without a valid PAM at the first site checked,
    # or indicates no match was found at all.
    all_results_df = pd.DataFrame(all_search_results)


    # 3. Perform Gene Association (using only sgRNAs for which a match was reported by find_exact_matches)
    print("\nAssociating matched sgRNAs (with valid PAM) with genes...")
    # Initialize gene columns in all_results_df
    gene_cols = ['Overlapping_Gene_Symbol', 'Overlapping_Locus_Tag', 'Overlapping_Product_Description',
                 'Nearest_Upstream_Gene_Symbol', 'Nearest_Upstream_Locus_Tag', 'Nearest_Upstream_Distance',
                 'Nearest_Downstream_Gene_Symbol', 'Nearest_Downstream_Locus_Tag', 'Nearest_Downstream_Distance']
    for col in gene_cols:
         if col not in all_results_df.columns:
            all_results_df[col] = pd.NA # Use pandas NA marker

    # Use match_results_for_genes which ONLY contains sgRNAs where find_exact_matches returned a valid hit
    if not match_results_for_genes:
        print("No sgRNA matches with valid PAM found, skipping gene association.")
    else:
        matches_df_for_genes = pd.DataFrame(match_results_for_genes)

        # Prepare sgRNA matches for PyRanges (0-based start)
        sgrna_pr_df = matches_df_for_genes[['Chromosome', 'Start', 'End', 'sgRNA_Sequence']].copy()
        # Convert Start/End which should be integers from find_exact_matches
        sgrna_pr_df['Start'] = pd.to_numeric(sgrna_pr_df['Start'], errors='coerce').astype('Int64') - 1 # Convert to 0-based
        sgrna_pr_df['End'] = pd.to_numeric(sgrna_pr_df['End'], errors='coerce').astype('Int64')
        sgrna_pr_df.dropna(subset=['Start', 'End'], inplace=True)

        if sgrna_pr_df.empty:
             print("WARNING: No valid coordinates found for matched sgRNAs with PAM. Skipping gene association.")
        else:
            # Filter out any potential negative starts after conversion
            sgrna_pr_df = sgrna_pr_df[sgrna_pr_df['Start'] >= 0]
            if sgrna_pr_df.empty:
                print("WARNING: No valid non-negative coordinates after 0-based conversion. Skipping gene association.")
            else:
                sgrna_pr = pr.PyRanges(sgrna_pr_df)

                # --- Find Overlapping Genes ---
                print("  - Finding overlapping genes...")
                # (Rest of gene association logic remains the same as before)
                overlaps_pr = sgrna_pr.join(genes_pr, how="left", apply_strand_suffix=False)
                overlaps_df = overlaps_pr.df

                overlaps_grouped = overlaps_df.groupby('sgRNA_Sequence').agg(
                    Overlapping_Gene_Symbol=('Gene_Symbol', lambda x: ';'.join(x.dropna().astype(str).unique())),
                    Overlapping_Locus_Tag=('Locus_Tag', lambda x: ';'.join(x.dropna().astype(str).unique())),
                    Overlapping_Product_Description=('Product_Description', lambda x: ';'.join(x.dropna().astype(str).unique()))
                ).reset_index()

                valid_overlap_sgrnas = overlaps_df.dropna(subset=['Locus_Tag'])['sgRNA_Sequence'].unique()
                overlaps_final_df = overlaps_grouped[overlaps_grouped['sgRNA_Sequence'].isin(valid_overlap_sgrnas)].copy()
                overlaps_final_df.replace({'': pd.NA}, inplace=True) # Use pandas NA

                sgRNAs_with_match_and_pam = set(matches_df_for_genes['sgRNA_Sequence'])
                sgRNAs_with_overlaps = set(overlaps_final_df['sgRNA_Sequence'])
                non_overlapping_sgrnas = list(sgRNAs_with_match_and_pam - sgRNAs_with_overlaps)

                nearest_upstream_final_df = pd.DataFrame()
                nearest_downstream_final_df = pd.DataFrame()

                if non_overlapping_sgrnas:
                    print(f"  - Finding nearest genes for {len(non_overlapping_sgrnas)} non-overlapping sgRNAs (with valid PAM)...")
                    non_overlapping_df = matches_df_for_genes[matches_df_for_genes['sgRNA_Sequence'].isin(non_overlapping_sgrnas)].copy()

                    if not non_overlapping_df.empty:
                        # Prepare for PyRanges (use 0-based Start from sgrna_pr_df if available, otherwise recalculate)
                        non_overlapping_pr_df = non_overlapping_df[['Chromosome', 'Start', 'End', 'sgRNA_Sequence']].copy()
                        non_overlapping_pr_df['Start'] = pd.to_numeric(non_overlapping_pr_df['Start'], errors='coerce').astype('Int64') - 1 # Convert to 0-based
                        non_overlapping_pr_df['End'] = pd.to_numeric(non_overlapping_pr_df['End'], errors='coerce').astype('Int64')
                        non_overlapping_pr_df.dropna(subset=['Start', 'End'], inplace=True)
                        non_overlapping_pr_df = non_overlapping_pr_df[non_overlapping_pr_df['Start'] >= 0] # Filter negative starts


                        if not non_overlapping_pr_df.empty:
                            non_overlapping_pr = pr.PyRanges(non_overlapping_pr_df)

                            # --- Find ONE Nearest (Upstream OR Downstream) ---
                            print("    - Finding single nearest gene (position will determine up/downstream)...")
                            nearest_pr = non_overlapping_pr.nearest(genes_pr, suffix="_gene", overlap=False)
                            nearest_df = nearest_pr.df
                            nearest_df_filtered = nearest_df[nearest_df['Distance'] >= 0].copy()

                            if not nearest_df_filtered.empty:
                                is_upstream = nearest_df_filtered['End_gene'] < nearest_df_filtered['Start']
                                is_downstream = nearest_df_filtered['Start_gene'] > nearest_df_filtered['End']
                                nearest_up_df = nearest_df_filtered[is_upstream].copy()
                                nearest_down_df = nearest_df_filtered[is_downstream].copy()

                                if not nearest_up_df.empty:
                                    nearest_upstream_final_df = nearest_up_df.groupby('sgRNA_Sequence').agg(
                                        Nearest_Upstream_Gene_Symbol=('Gene_Symbol', lambda x: ';'.join(x.dropna().astype(str).unique())),
                                        Nearest_Upstream_Locus_Tag=('Locus_Tag', lambda x: ';'.join(x.dropna().astype(str).unique())),
                                        Nearest_Upstream_Distance=('Distance', 'min')
                                    ).reset_index()
                                    nearest_upstream_final_df.replace({'': pd.NA}, inplace=True)
                                else: print("    - No valid nearest upstream genes found.")

                                if not nearest_down_df.empty:
                                    nearest_downstream_final_df = nearest_down_df.groupby('sgRNA_Sequence').agg(
                                        Nearest_Downstream_Gene_Symbol=('Gene_Symbol', lambda x: ';'.join(x.dropna().astype(str).unique())),
                                        Nearest_Downstream_Locus_Tag=('Locus_Tag', lambda x: ';'.join(x.dropna().astype(str).unique())),
                                        Nearest_Downstream_Distance=('Distance', 'min')
                                    ).reset_index()
                                    nearest_downstream_final_df.replace({'': pd.NA}, inplace=True)
                                else: print("    - No valid nearest downstream genes found.")
                            else: print("    - No nearest genes found (Distance >= 0) for non-overlapping sgRNAs.")
                        else: print("    - No valid coordinates for non-overlapping sgRNAs after filtering.")
                    else: print("    - Non-overlapping dataframe empty.")
                else: print("  - No sgRNAs required nearest gene search.")


                # --- Combine All Results ---
                print("  - Combining overlap, upstream, and downstream results with main table...")
                final_results_df = all_results_df.copy() # Start with the full table

                # Merge gene info based on sgRNA_Sequence. Only rows where find_exact_matches succeeded will have gene info added.
                if not overlaps_final_df.empty:
                     final_results_df = pd.merge(final_results_df, overlaps_final_df, on='sgRNA_Sequence', how='left')
                if not nearest_upstream_final_df.empty:
                     final_results_df = pd.merge(final_results_df, nearest_upstream_final_df, on='sgRNA_Sequence', how='left')
                if not nearest_downstream_final_df.empty:
                     final_results_df = pd.merge(final_results_df, nearest_downstream_final_df, on='sgRNA_Sequence', how='left')

                all_results_df = final_results_df


    # 4. Format Final Report
    print(f"\nGenerating final CSV report: {final_csv_report}")
    try:
        # Convert coordinates and distances to nullable Ints first
        numeric_cols_to_convert = ['Start', 'End', 'Nearest_Upstream_Distance', 'Nearest_Downstream_Distance']
        for col in numeric_cols_to_convert:
            if col in all_results_df.columns:
                # Coerce errors, then convert to Int64
                all_results_df[col] = pd.to_numeric(all_results_df[col], errors='coerce').astype('Int64')


        # Define dictionary for filling NAs in string columns
        string_fill_dict = {
            'Chromosome': 'N/A', 'Strand_Match': 'N/A', 'PAM_Sequence': 'N/A',
            'Overlapping_Gene_Symbol': 'N/A', 'Overlapping_Locus_Tag': 'N/A', 'Overlapping_Product_Description': 'N/A',
            'Nearest_Upstream_Gene_Symbol': 'N/A', 'Nearest_Upstream_Locus_Tag': 'N/A',
            'Nearest_Downstream_Gene_Symbol': 'N/A', 'Nearest_Downstream_Locus_Tag': 'N/A',
        }
        cols_to_fill_strings = {k: v for k, v in string_fill_dict.items() if k in all_results_df.columns}
        all_results_df.fillna(cols_to_fill_strings, inplace=True)

        # Handle Boolean columns (Match_Found, PAM_Found) - fill NA, convert to String
        if 'Match_Found' in all_results_df.columns:
             all_results_df['Match_Found'] = all_results_df['Match_Found'].fillna(False).astype(str)
        if 'PAM_Found' in all_results_df.columns:
             all_results_df['PAM_Found'] = all_results_df['PAM_Found'].fillna(False).astype(str)


        # Convert Int64 columns to string *after* filling NAs in other cols
        int_cols_to_str = ['Start', 'End', 'Nearest_Upstream_Distance', 'Nearest_Downstream_Distance']
        for col in int_cols_to_str:
             if col in all_results_df.columns:
                  all_results_df[col] = all_results_df[col].astype(str).replace('<NA>', 'N/A')


        # Define final column order
        report_columns = [
            'sgRNA_Sequence', 'Match_Found', 'PAM_Found', # Moved PAM_Found earlier
            'Chromosome', 'Start', 'End', 'Strand_Match', 'PAM_Sequence',
            'Overlapping_Gene_Symbol', 'Overlapping_Locus_Tag', 'Overlapping_Product_Description',
            'Nearest_Upstream_Gene_Symbol', 'Nearest_Upstream_Locus_Tag', 'Nearest_Upstream_Distance',
            'Nearest_Downstream_Gene_Symbol', 'Nearest_Downstream_Locus_Tag', 'Nearest_Downstream_Distance'
        ]
        final_columns = [col for col in report_columns if col in all_results_df.columns]
        final_columns.extend([col for col in all_results_df.columns if col not in final_columns])

        all_results_df = all_results_df[final_columns]

        all_results_df.to_csv(final_csv_report, index=False, na_rep='N/A')
        print("CSV report generated successfully.")
    except Exception as e:
        print(f"ERROR: Failed to generate CSV report: {e}")
        import traceback
        traceback.print_exc()

    # 5. Generate Summary Plot (Shows Found vs Not Found based on *any* match site)
    print(f"\nGenerating summary plot: {summary_plot_file}")
    try:
        plt.figure(figsize=(6, 5))
        # Use found_count (any match) and not_found_count for the plot
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