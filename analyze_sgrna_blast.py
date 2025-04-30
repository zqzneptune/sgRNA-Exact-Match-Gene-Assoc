# -*- coding: utf-8 -*-
import os
import subprocess
import pandas as pd
import matplotlib.pyplot as plt
import sys
import shutil
import argparse
import glob # Needed for checking BLAST db files

# --- Configuration ---
DEFAULT_SGRNA_INPUT_FILE = 'sgrnas.txt'
DEFAULT_TARGET_GENOME_FASTA = 'target_genome.fasta'
DEFAULT_OUTPUT_DIR = 'blast_analysis_output'
DEFAULT_BLASTN_PATH = 'blastn'
DEFAULT_MAKEBLASTDB_PATH = 'makeblastdb'

# NOTE: Similarity thresholds are less relevant now as we enforce 100% identity in BLAST
# HIGH_SIMILARITY_THRESHOLD = 100.0
# MEDIUM_SIMILARITY_THRESHOLD = 85.0 # Kept for reference, but not used for BLAST filtering directly

# BLAST output format:
# qseqid: Query Seq-id (our sgRNA ID)
# qseq: Query sequence (the sgRNA itself)
# score: Raw score
# evalue: Expect value
# pident: Percentage of identical matches
# gaps: Total number of gaps
# sstrand: Subject strand (plus/minus)
# sstart: Start of alignment in subject (target genome)
# send: End of alignment in subject (target genome)
# qlen: Query length
# length: Alignment length
BLAST_OUTFMT = "6 qseqid qseq score evalue pident gaps sstrand sstart send qlen length"
BLAST_COLUMNS = ['sgRNA_ID', 'sgRNA_Sequence', 'Score', 'Expect',
                 'Identities_perc', 'Gaps', 'Strand',
                 'Subject_Start', 'Subject_End',
                 'Query_Length', 'Alignment_Length'] # Added qlen, length

# --- Helper Functions ---

def check_command(command_path, command_name):
    """Checks if a command exists and is executable."""
    resolved_path = shutil.which(command_path)
    if not resolved_path:
        print(f"ERROR: '{command_name}' command not found using path '{command_path}'.")
        print(f"Please ensure NCBI BLAST+ is installed and '{command_path}' is correct or '{command_name}' is in your system's PATH.")
        sys.exit(1)
    print(f"Found {command_name} at: {resolved_path}")
    return resolved_path

def run_subprocess(command_list, error_message):
    """Runs a subprocess command and handles errors."""
    print(f"\nRunning command: {' '.join(map(str, command_list))}") # Ensure all elements are strings
    try:
        command_list_str = [str(item) for item in command_list]
        process = subprocess.run(command_list_str, check=True, capture_output=True, text=True, errors='replace')
        print("Command executed successfully.")
        stdout_preview = process.stdout[:500] + ("..." if len(process.stdout) > 500 else "")
        stderr_preview = process.stderr[:500] + ("..." if len(process.stderr) > 500 else "")
        if stdout_preview.strip():
             print("STDOUT Preview:\n" + stdout_preview)
        if stderr_preview.strip():
             print("STDERR Preview:\n" + stderr_preview)
    except FileNotFoundError:
        print(f"ERROR: Command not found: {command_list_str[0]}. Is it installed and in PATH?")
        sys.exit(1)
    except subprocess.CalledProcessError as e:
        print(f"ERROR: {error_message} failed with exit code {e.returncode}")
        cmd_str = ' '.join(map(str, e.cmd))
        print(f"Command: {cmd_str}")
        print(f"Stderr:\n{e.stderr}")
        print(f"Stdout:\n{e.stdout}")
        sys.exit(1)
    except Exception as e:
        print(f"An unexpected error occurred while running {' '.join(map(str, command_list_str))}: {e}")
        sys.exit(1)

def prepare_sgrna_fasta(sgrna_input_file, output_fasta_path, force_overwrite=False):
    """Reads a file with one sgRNA per line and writes a FASTA file."""
    if os.path.exists(output_fasta_path) and not force_overwrite:
        print(f"\nINFO: Using existing sgRNA FASTA file: {output_fasta_path}")
        sgrnas = []
        try:
            with open(sgrna_input_file, 'r') as f:
                sgrnas = [line.strip().upper() for line in f if line.strip()]
        except FileNotFoundError:
            print(f"ERROR: sgRNA input file not found: {sgrna_input_file}")
            sys.exit(1)
        if not sgrnas:
            print(f"ERROR: No sgRNA sequences found in {sgrna_input_file}")
            sys.exit(1)
        unique_sgrnas = sorted(list(set(sgrnas)))
        unique_sgrnas_valid = [seq for seq in unique_sgrnas if all(c in 'ATCGN' for c in seq)]
        print(f"Read {len(sgrnas)} total sgRNAs, {len(unique_sgrnas_valid)} unique & valid.")
        return unique_sgrnas_valid

    print(f"\nPreparing sgRNA FASTA file from {sgrna_input_file}...")
    sgrnas = []
    try:
        with open(sgrna_input_file, 'r') as f:
            sgrnas = [line.strip().upper() for line in f if line.strip()]
    except FileNotFoundError:
        print(f"ERROR: sgRNA input file not found: {sgrna_input_file}")
        sys.exit(1)

    if not sgrnas:
        print(f"ERROR: No sgRNA sequences found in {sgrna_input_file}")
        sys.exit(1)

    print(f"Read {len(sgrnas)} total sgRNA sequences.")
    unique_sgrnas = sorted(list(set(sgrnas)))
    print(f"Processing {len(unique_sgrnas)} unique sgRNA sequences.")

    valid_sgrnas_count = 0
    unique_sgrnas_written = []
    try:
        with open(output_fasta_path, 'w') as outfile:
            for i, seq in enumerate(unique_sgrnas):
                if not all(c in 'ATCGN' for c in seq):
                    print(f"WARNING: Skipping invalid sequence in FASTA generation (non-ATCGN): {seq}")
                    continue
                outfile.write(f">sgRNA_{i+1}\n")
                outfile.write(f"{seq}\n")
                valid_sgrnas_count += 1
                unique_sgrnas_written.append(seq)
        print(f"Successfully created sgRNA FASTA file: {output_fasta_path} with {valid_sgrnas_count} sequences.")
        return unique_sgrnas_written
    except IOError as e:
        print(f"ERROR: Could not write sgRNA FASTA file to {output_fasta_path}: {e}")
        sys.exit(1)


def create_blast_db(makeblastdb_path, genome_fasta, db_name, force_overwrite=False):
    """Creates a BLAST database from the target genome, skipping if it exists."""
    print("\nChecking for existing BLAST database...")
    db_files = glob.glob(f"{db_name}.n*")

    if db_files and not force_overwrite:
        print(f"Found existing BLAST database files for prefix: {db_name}")
        print("Skipping makeblastdb step. Use --force to regenerate.")
        return
    elif db_files and force_overwrite:
        print(f"Found existing BLAST database files, but --force specified. Regenerating...")
        pass

    print("Creating BLAST database...")
    if not os.path.exists(genome_fasta):
        print(f"ERROR: Target genome FASTA file not found: {genome_fasta}")
        sys.exit(1)

    command = [
        makeblastdb_path,
        "-in", genome_fasta,
        "-dbtype", "nucl",
        "-out", db_name,
        "-parse_seqids"
    ]
    run_subprocess(command, "makeblastdb")
    print(f"BLAST database created/updated: {db_name}")


def run_blastn(blastn_path, query_fasta, db_name, output_tsv, force_overwrite=False):
    """Runs blastn alignment for PERFECT matches only, skipping if the output file exists."""
    print("\nChecking for existing BLAST results file...")
    if os.path.exists(output_tsv) and not force_overwrite:
        # Check if parameters used previously match current strict criteria? Hard to do reliably.
        # Assume user knows if existing file used the strict parameters or uses --force.
        print(f"Found existing BLAST results file: {output_tsv}")
        print("Skipping blastn step. Use --force to re-run alignment (ensure strict parameters are used).")
        if os.path.getsize(output_tsv) == 0:
             print("WARNING: Existing BLAST results file is empty.")
        return
    elif os.path.exists(output_tsv) and force_overwrite:
        print(f"Found existing BLAST results file, but --force specified. Re-running BLASTn with strict parameters...")

    print("Running BLASTn search for PERFECT matches (100% identity, 100% query coverage, 0 gaps)...")
    if not os.path.exists(query_fasta):
         print(f"ERROR: Query FASTA file not found for BLAST: {query_fasta}")
         sys.exit(1)
    if not glob.glob(f"{db_name}.n*"):
         print(f"ERROR: BLAST database files not found for prefix: {db_name}. Cannot run BLASTn.")
         sys.exit(1)

    num_threads = os.cpu_count()
    if num_threads is None:
        print("Warning: Could not detect number of CPUs. Using 1 thread for BLAST.")
        num_threads = 1

    command = [
        blastn_path,
        "-query", query_fasta,
        "-db", db_name,
        "-out", output_tsv,
        "-outfmt", BLAST_OUTFMT,      # Use the updated format string
        "-task", "blastn-short",     # Optimized for short sequences
        "-perc_identity", "100",     # Require 100% identity
        "-qcov_hsp_perc", "100",     # Require 100% query coverage
        "-evalue", "10",             # E-value threshold (less critical with strict identity/coverage)
        "-num_threads", str(num_threads),
        "-max_target_seqs", "1"      # Report only the best hit per query
        # Note: We rely on pident=100 and qcov=100 to implicitly ensure no gaps in the reported HSP.
        # Setting high gap penalties (-gapopen, -gapextend) is another way, but the above is more direct.
    ]
    run_subprocess(command, "blastn")
    print(f"BLASTn search complete (strict criteria). Results saved to: {output_tsv}")

def process_blast_results(blast_output_tsv, unique_sgrna_list):
    """Parses STRICT BLAST output, selects best hit per sgRNA, and categorizes."""
    print("\nProcessing STRICT BLAST results...")
    results_columns = BLAST_COLUMNS + ['Category']

    try:
        if not os.path.exists(blast_output_tsv):
             print(f"ERROR: BLAST output file not found: {blast_output_tsv}. Cannot process results.")
             results_df = pd.DataFrame({'sgRNA_Sequence': unique_sgrna_list})
             results_df = results_df.reindex(columns=results_columns, fill_value=pd.NA)
             results_df['Category'] = "No Hit"
             print("Created empty results structure assuming no hits.")
             return results_df

        # Specify dtypes for all columns
        dtype_spec = {
            'sgRNA_ID': str, 'sgRNA_Sequence': str, 'Score': float,
            'Expect': object, 'Identities_perc': float, 'Gaps': float,
            'Strand': str, 'Subject_Start': float, 'Subject_End': float,
            'Query_Length': float, 'Alignment_Length': float
        }
        try:
            blast_df = pd.read_csv(blast_output_tsv, sep='\t', names=BLAST_COLUMNS, dtype=dtype_spec, na_values=['N/A'])
        except pd.errors.EmptyDataError:
            print("WARNING: BLAST output file is empty. No perfect hits found.")
            blast_df = pd.DataFrame(columns=BLAST_COLUMNS)
        except Exception as e:
            print(f"ERROR: Failed to read BLAST output file {blast_output_tsv}. Error: {e}")
            print("Check the file for formatting issues or consider regenerating it with --force.")
            results_df = pd.DataFrame({'sgRNA_Sequence': unique_sgrna_list})
            results_df = results_df.reindex(columns=results_columns, fill_value=pd.NA)
            results_df['Category'] = "Processing Error"
            return results_df

        print(f"Read {len(blast_df)} total perfect BLAST hits (meeting 100% id, 100% cov).")

        # Data Cleaning & Type Conversion
        blast_df['Expect'] = pd.to_numeric(blast_df['Expect'], errors='coerce')
        # Convert relevant columns to nullable integers
        int_cols = ['Gaps', 'Subject_Start', 'Subject_End', 'Query_Length', 'Alignment_Length']
        for col in int_cols:
            # Check column exists before conversion
            if col in blast_df.columns:
                blast_df[col] = pd.to_numeric(blast_df[col], errors='coerce').astype('Int64')

        # Ensure start <= end for consistency
        if not blast_df.empty and 'Subject_Start' in blast_df.columns and 'Subject_End' in blast_df.columns:
            mask = blast_df['Subject_Start'].notna() & blast_df['Subject_End'].notna() & (blast_df['Subject_Start'] > blast_df['Subject_End'])
            start_temp = blast_df.loc[mask, 'Subject_Start'].copy()
            end_temp = blast_df.loc[mask, 'Subject_End'].copy()
            blast_df.loc[mask, 'Subject_Start'] = end_temp
            blast_df.loc[mask, 'Subject_End'] = start_temp

        # Select best hit (should be redundant with max_target_seqs=1, but safe)
        best_hits_df = pd.DataFrame(columns=BLAST_COLUMNS) # Initialize empty
        if not blast_df.empty:
            # Sorting is less critical here as all hits are 'perfect', but keep for consistency
            blast_df.sort_values(
                by=['sgRNA_ID', 'Score', 'Expect'],
                ascending=[True, False, True],
                na_position='last',
                inplace=True
            )
            best_hits_df = blast_df.drop_duplicates(subset='sgRNA_ID', keep='first').copy()
            print(f"Selected {len(best_hits_df)} best perfect hits (one per unique query sgRNA with hits).")


        sgrna_master_df = pd.DataFrame({'sgRNA_Sequence': unique_sgrna_list})
        sgrna_master_df['sgRNA_ID'] = [f'sgRNA_{i+1}' for i in range(len(unique_sgrna_list))]

        # Merge BLAST best hits with the master list
        merged_df = pd.merge(sgrna_master_df, best_hits_df, on=['sgRNA_ID', 'sgRNA_Sequence'], how='left')

        # --- Categorize Similarity (Simplified for Strict BLAST) ---
        def categorize(row):
            # If a score exists, BLAST found a perfect match based on the strict criteria
            if pd.notna(row['Score']):
                # Sanity check (optional, should always be true if BLAST worked as expected)
                is_perfect = (
                    pd.notna(row['Identities_perc']) and row['Identities_perc'] == 100.0 and
                    pd.notna(row['Gaps']) and row['Gaps'] == 0 and
                    pd.notna(row['Query_Length']) and pd.notna(row['Alignment_Length']) and
                    row['Query_Length'] == row['Alignment_Length']
                )
                if is_perfect:
                    return "High" # Perfect match found
                else:
                    # This case indicates an unexpected BLAST result despite the parameters
                    print(f"WARNING: Hit found for {row['sgRNA_ID']} but failed strict sanity check:"
                          f" Pident={row['Identities_perc']}, Gaps={row['Gaps']},"
                          f" Qlen={row['Query_Length']}, Alen={row['Alignment_Length']}. Categorizing as Error.")
                    return "Processing Error"
            else:
                # No hit reported by BLAST means no perfect match found
                return "No Hit"

        merged_df['Category'] = merged_df.apply(categorize, axis=1)
        print("Categorization complete (based on perfect matches).")
        print("\nCategory Counts:")
        category_order = ['High', 'No Hit', 'Processing Error'] # Medium/Low are not expected now
        print(merged_df['Category'].value_counts().reindex(category_order, fill_value=0))

        return merged_df

    except Exception as e:
        print(f"ERROR: Failed severely during BLAST result processing: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


def generate_report_and_plot(results_df, csv_path, plot_path):
    """Generates the final CSV report and bar plot."""
    print("\nGenerating CSV report...")

    # Select and rename columns for the final report
    # Include Query_Length and Alignment_Length for clarity
    report_columns = ['sgRNA_Sequence', 'Score', 'Expect', 'Identities (%)',
                      'Gaps', 'Strand', 'Subject_Start', 'Subject_End',
                      'Query_Length', 'Alignment_Length', 'Category']
    # Ensure all expected columns exist in results_df before selecting
    cols_to_select = [col for col in ['sgRNA_Sequence', 'Score', 'Expect', 'Identities_perc',
                                     'Gaps', 'Strand', 'Subject_Start', 'Subject_End',
                                     'Query_Length', 'Alignment_Length', 'Category'] if col in results_df.columns]
    report_df = results_df[cols_to_select].copy()

    # Rename if the column exists
    if 'Identities_perc' in report_df.columns:
        report_df.rename(columns={'Identities_perc': 'Identities (%)'}, inplace=True)

    # Format numeric columns, handling potential NA values
    if 'Score' in report_df.columns:
        report_df['Score'] = report_df['Score'].round(2)
    if 'Expect' in report_df.columns:
        report_df['Expect'] = report_df['Expect'].apply(lambda x: f"{x:.2e}" if pd.notna(x) and isinstance(x, (int, float)) else x)
    if 'Identities (%)' in report_df.columns:
        report_df['Identities (%)'] = report_df['Identities (%)'].round(2)

    # Handle NAs for final report - fill numeric NAs suitable for CSV
    fill_values = {
        'Score': 0.0, 'Expect': 'N/A', 'Identities (%)': 0.0,
        'Gaps': pd.NA, 'Strand': 'N/A', 'Subject_Start': pd.NA,
        'Subject_End': pd.NA, 'Query_Length': pd.NA, 'Alignment_Length': pd.NA
    }
    # Apply fillna only for columns that exist in the dataframe
    report_df.fillna({k: v for k, v in fill_values.items() if k in report_df.columns}, inplace=True)


    # Convert nullable Int64 columns back to string for CSV output with 'N/A' for missing values
    int_like_cols = ['Gaps', 'Subject_Start', 'Subject_End', 'Query_Length', 'Alignment_Length']
    for col in int_like_cols:
         if col in report_df.columns:
             # Ensure the column is treated as string, replace '<NA>' which is pandas' default string rep for Int64 NA
             report_df[col] = report_df[col].astype(str).replace('<NA>', 'N/A')


    try:
        # Ensure column order in output CSV, using the potentially modified list
        final_report_columns = [col for col in report_columns if col in report_df.columns]
        report_df = report_df[final_report_columns]
        report_df.to_csv(csv_path, index=False, na_rep='N/A')
        print(f"CSV report saved to: {csv_path}")
    except IOError as e:
        print(f"ERROR: Could not write CSV report to {csv_path}: {e}")
    except Exception as e:
        print(f"An unexpected error occurred generating CSV: {e}")
        import traceback
        traceback.print_exc()


    print("\nGenerating similarity category plot...")
    try:
        # Categories expected now are primarily High, No Hit, maybe Processing Error
        category_order = ['High', 'No Hit', 'Processing Error'] # Adjusted expected categories
        if 'Category' not in results_df.columns:
             print("ERROR: 'Category' column not found in results DataFrame. Cannot generate plot.")
             return

        valid_categories = results_df['Category'][results_df['Category'].isin(category_order)]
        category_counts = valid_categories.value_counts().reindex(category_order, fill_value=0)

        plt.figure(figsize=(8, 6))
        # Adjust colors if Medium/Low are definitely gone
        colors = {'High': '#d73027', 'No Hit': '#4575b4', 'Processing Error': '#bdbdbd'}
        bar_colors = [colors.get(cat, '#cccccc') for cat in category_counts.index] # Get color by name

        bars = plt.bar(category_counts.index, category_counts.values, color=bar_colors)

        plt.xlabel("Similarity Category (Perfect Match)")
        plt.ylabel("Number of Unique sgRNAs")
        plt.title("sgRNA Perfect Match Analysis vs Target Genome")
        plt.xticks(rotation=0)
        plt.grid(axis='y', linestyle='--', alpha=0.7)
        plt.ylim(bottom=0)

        for bar in bars:
            yval = bar.get_height()
            if yval > 0:
                plt.text(bar.get_x() + bar.get_width()/2.0, yval, int(yval), va='bottom', ha='center', fontsize=9)

        plt.tight_layout()
        plt.savefig(plot_path, dpi=300)
        print(f"Similarity plot saved to: {plot_path}")
    except IOError as e:
        print(f"ERROR: Could not save plot to {plot_path}: {e}")
    except Exception as e:
        print(f"An unexpected error occurred generating the plot: {e}")
        import traceback
        traceback.print_exc()
    # plt.show()

# --- Main Execution ---
def main():
    parser = argparse.ArgumentParser(description="Analyze sgRNA **perfect matches** against a target genome using BLAST.") # Updated description
    parser.add_argument("-s", "--sgrna_file", default=DEFAULT_SGRNA_INPUT_FILE,
                        help=f"Input file containing sgRNA sequences, one per line (default: {DEFAULT_SGRNA_INPUT_FILE})")
    parser.add_argument("-g", "--genome_file", default=DEFAULT_TARGET_GENOME_FASTA,
                        help=f"Input FASTA file for the target genome (default: {DEFAULT_TARGET_GENOME_FASTA})")
    parser.add_argument("-o", "--output_dir", default=DEFAULT_OUTPUT_DIR,
                        help=f"Directory to save output files (default: {DEFAULT_OUTPUT_DIR})")
    parser.add_argument("--blastn", default=DEFAULT_BLASTN_PATH,
                        help=f"Path to the blastn executable (default: '{DEFAULT_BLASTN_PATH}' or found in PATH)")
    parser.add_argument("--makeblastdb", default=DEFAULT_MAKEBLASTDB_PATH,
                        help=f"Path to the makeblastdb executable (default: '{DEFAULT_MAKEBLASTDB_PATH}' or found in PATH)")
    parser.add_argument("--force", action='store_true',
                        help="Force regeneration of BLAST database and BLAST results, even if files exist.")

    args = parser.parse_args()

    output_dir = args.output_dir
    sgrna_input_file = args.sgrna_file
    target_genome_fasta = args.genome_file
    force_overwrite = args.force

    print("--- Starting sgRNA PERFECT MATCH Analysis ---") # Updated title
    makeblastdb_path = check_command(args.makeblastdb, "makeblastdb")
    blastn_path = check_command(args.blastn, "blastn")

    if not os.path.exists(output_dir):
        print(f"Creating output directory: {output_dir}")
        os.makedirs(output_dir)
    else:
        print(f"Output directory already exists: {output_dir}")

    blast_db_name = os.path.join(output_dir, os.path.splitext(os.path.basename(target_genome_fasta))[0] + '_db')
    sgrna_fasta_file = os.path.join(output_dir, 'query_sgrnas.fasta')
    blast_output_file = os.path.join(output_dir, 'blast_results_perfect.tsv') # Changed filename
    final_csv_report = os.path.join(output_dir, 'sgrna_perfect_match_report.csv') # Changed filename
    similarity_plot_file = os.path.join(output_dir, 'sgrna_perfect_match_categories.png') # Changed filename

    unique_sgrnas = prepare_sgrna_fasta(sgrna_input_file, sgrna_fasta_file, force_overwrite)
    if not unique_sgrnas:
         print("Exiting: No valid unique sgRNAs found or processed.")
         sys.exit(1)

    create_blast_db(makeblastdb_path, target_genome_fasta, blast_db_name, force_overwrite)

    # Run BLAST with strict parameters
    run_blastn(blastn_path, sgrna_fasta_file, blast_db_name, blast_output_file, force_overwrite)

    # Process results expecting only perfect hits or no hits
    results_df = process_blast_results(blast_output_file, unique_sgrnas)

    # Generate report and plot (now showing High/No Hit)
    generate_report_and_plot(results_df, final_csv_report, similarity_plot_file)

    print("\n--- Analysis Complete ---")
    print(f"Results summary saved in: {output_dir}")

if __name__ == "__main__":
    # example_sgrnas_text = """
    # GGAAGGGGTGGCTTCGAGCGT
    # GATCGCCGGAACGTTCACACA
    # CCGGAGGTTCGCCGGGAAGG
    # ACTTCCGTGCCATCAATAAA
    # GCCACAGCCACATTCATTCT
    # CACATTCATTCTGGGCTTTA
    # GAGATTAATTAGCGACTGTT
    # ATTTTACCTGAACCATAAT
    # GATACGAAACCATTGTTGAC
    # GCTGGGAGTTGGTGCTGGAT
    # TCGCTTCCGGCGTGGCAAAG
    # TGATCTGCTCGGCCTGTTCC
    # GAAGTTGTAGAGACGCACAC
    # TTCCATTATCGAACGACAAT
    # GGAACATCCAGATGGAAATA
    # CAGCCGCATCGCGCCGGCA
    # TGCCAGCTTTTCGGCATTTG
    # GCAAAAGCCTGCGGGAGAAA
    # CAGCGTACCGAAGCGCAAGC
    # ACTTCCGTGCAATCAATAAA
    # CACATTCATTCTGGGCGTTA
    # ONLYPARTIALMATCHGGGCTTTA
    # INVALIDSEQUENCE###INVALID
    # """ # Added sequences with 1 mismatch for testing

    # if not os.path.exists(DEFAULT_SGRNA_INPUT_FILE):
    #     print(f"Creating dummy input file: {DEFAULT_SGRNA_INPUT_FILE}")
    #     with open(DEFAULT_SGRNA_INPUT_FILE, "w") as f:
    #         f.write(example_sgrnas_text)

    # if not os.path.exists(DEFAULT_TARGET_GENOME_FASTA):
    #      print(f"Creating dummy target genome file: {DEFAULT_TARGET_GENOME_FASTA}")
    #      print("WARNING: Using a dummy genome. Results will not be biologically meaningful.")
    #      with open(DEFAULT_TARGET_GENOME_FASTA, "w") as f:
    #           f.write(">Dummy_Target_Chromosome_BW25113_like\n")
    #           # Genome contains perfect matches for some, but NOT the 1-mismatch versions
    #           dummy_genome = "N"*1000 + \
    #                          "GGAAGGGGTGGCTTCGAGCGT" + "N"*500 + \
    #                          "GATCGCCGGAACGTTCACACA" + "N"*500 + \
    #                          "CCGGAGGTTCGCCGGGAAGG" + "N"*500 + \
    #                          "ACTTCCGTGCCATCAATAAA" + "N"*500 + # The original perfect match
    #                          "GCCACAGCCACATTCATTCT" + "N"*500 + \
    #                          "CACATTCATTCTGGGCTTTA" + "N"*500 + # The original perfect match
    #                          "GAGATTAATTAGCGACTGTT" + "N"*500 + \
    #                          "TCGCTTCCGGCGTGGCAAAG" + "N"*500 + \
    #                          "AAAAAAAAAAAAAAAAAAAA" + "N"*1000 + \
    #                          "GGGCTTTA" + "N"*100 # Match for the partial one, but won't meet qcov=100
    #           f.write(dummy_genome + "\n")

    main()