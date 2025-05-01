# sgRNA Exact Match and Gene Association Finder

## Summary

This Python script identifies exact sequence matches for a list of sgRNA sequences within a target genome FASTA file. Unlike alignment tools like BLAST, this script performs **only exact string matching**, meaning no mismatches or gaps are allowed relative to the provided sgRNA sequences.

For sgRNAs that have an exact match in the genome, the script further associates them with gene features based on a provided gene annotation table. It identifies:
1.  **Overlapping Genes:** Genes whose coordinates directly overlap the sgRNA match location.
2.  **Nearest Genes:** For sgRNAs landing in intergenic regions, it finds the single closest gene upstream and the single closest gene downstream along the chromosome coordinates.

The script outputs a detailed CSV report summarizing the match status and gene associations for each input sgRNA, as well as a summary plot showing the number of sgRNAs found in the genome.

## Prerequisites

1.  **Python 3:** Tested with Python 3.9+, but should work with most recent Python 3 versions.
2.  **Required Python Libraries:** Install the necessary libraries using pip:
    ```bash
    pip install pandas matplotlib biopython pyranges numpy
    ```
3.  **Input Files:** You must provide the following files:
    *   A text file containing your sgRNA sequences.
    *   A FASTA file containing the target genome sequence(s).
    *   A gene annotation table file (CSV or TSV format).

## Usage

Run the script from your command line:

```bash
python analyze_sgrna_seqmatch.py -s <sgrna_file> -g <genome_fasta> -a <gene_table_file> -o <output_directory>
```

**Required Arguments:**

*   `-s`, `--sgrna_file`: Path to the input file containing sgRNA sequences (one sequence per line).
*   `-g`, `--genome_file`: Path to the input FASTA file for the target genome sequence(s).
*   `-a`, `--gene_table_file`: Path to the input gene annotation table file (CSV or TSV format). See "Input File Formats" below for required columns.
*   `-o`, `--output_dir`: Path to the directory where output files will be saved. The directory will be created if it doesn't exist.

## Input File Formats

1.  **sgRNA File (`-s`)**:
    *   Plain text file.
    *   Each line should contain exactly one sgRNA sequence.
    *   Sequences should consist only of A, T, C, G, N characters (case-insensitive, will be converted to uppercase). Invalid sequences will be skipped with a warning.

2.  **Genome FASTA File (`-g`)**:
    *   Standard FASTA format.
    *   Can contain one or multiple sequences (e.g., chromosomes, contigs).
    *   **Important Assumption:** If the FASTA file contains multiple sequences, the script currently assumes that all coordinates in the `gene_table_file` refer to the **first sequence entry** in this FASTA file. The ID of this first sequence will be used as the 'Chromosome' name for all gene associations. Ensure this matches your data context.

3.  **Gene Table File (`-a`)**:
    *   A text file, typically tab-separated (TSV) or comma-separated (CSV). The script attempts to automatically detect the separator. Lines starting with `#` are ignored as comments.
    *   **Required Columns** (Header names must match exactly):
        *   `Gene Name`: The common name or symbol for the gene.
        *   `Accession-1`: A locus tag or other unique identifier for the gene.
        *   `Left-End-Position`: The starting coordinate of the gene.
        *   `Right-End-Position`: The ending coordinate of the gene.
        *   `Product`: A description of the gene product.
    *   **Coordinate System:** The script assumes `Left-End-Position` and `Right-End-Position` are **1-based inclusive** coordinates (standard biological representation).
    *   **Data Validity:** Rows with non-numeric, blank, or invalid coordinates (e.g., Left >= Right) will be filtered out with a warning.
    *   **Strand Information:** The script currently does not use strand information from this table (as it's not present in the example). A dummy '+' strand is assigned.

## Output Files

The script generates the following files in the specified output directory (`-o`):

1.  **CSV Report (`<sgrna_basename>_match_gene_report.csv`)**:
    *   A comprehensive report detailing the analysis for each unique input sgRNA.
    *   Columns include:
        *   `sgRNA_Sequence`: The unique input sgRNA sequence.
        *   `Match_Found`: Boolean (True/False) indicating if an exact match was found in the genome.
        *   `Chromosome`: The ID of the genome sequence where the match was found (typically the first sequence ID from the FASTA). 'N/A' if no match.
        *   `Start`: 1-based start coordinate of the exact match. 'N/A' if no match.
        *   `End`: 1-based end coordinate of the exact match. 'N/A' if no match.
        *   `Strand_Match`: Strand of the genome where the match occurred ('+' for forward, '-' for reverse complement). 'N/A' if no match.
        *   `Overlapping_Gene_Symbol`: Symbol(s) of gene(s) directly overlapping the sgRNA match. Multiple separated by ';'. 'N/A' if no overlap.
        *   `Overlapping_Locus_Tag`: Locus tag(s) of overlapping gene(s). Multiple separated by ';'. 'N/A' if no overlap.
        *   `Overlapping_Product_Description`: Product description(s) of overlapping gene(s). 'N/A' if no overlap.        

2.  **Summary Plot (`<sgrna_basename>_match_summary.png`)**:
    *   A simple bar chart visualizing the total number of unique sgRNAs that had an exact match ('Found') versus those that did not ('Not Found') in the target genome.

## Core Logic and Rationale

1.  **Input Parsing:** Uses `argparse` to handle command-line arguments, ensuring all required inputs are provided.
2.  **sgRNA Reading (`read_sgrnas`):** Reads the sgRNA file line by line, converts sequences to uppercase, validates characters (ATCGN), and stores unique valid sequences in a sorted list. This ensures clean input for matching.
3.  **Genome Reading (`read_genome`):** Uses `BioPython.SeqIO` to parse the FASTA file. Stores sequences in a dictionary for quick access. It identifies the sequence ID of the *first* record, assuming this corresponds to the chromosome relevant to the gene table.
4.  **Gene Table Reading (`read_gene_table`):** Uses `pandas.read_csv` with separator auto-detection. It validates required columns, explicitly filters rows with non-numeric or blank coordinate values, converts coordinates to integers, assigns the assumed chromosome name, renames columns, converts 1-based coordinates to 0-based for `pyranges`, adds a dummy strand, performs final coordinate validation (Start >= 0, Start < End), and finally creates a `pyranges.PyRanges` object for efficient genomic interval operations.
5.  **Exact Sequence Matching (`find_exact_matches`):** For each unique sgRNA, this function:
    *   Calculates its reverse complement using `Bio.Seq`.
    *   Iterates through the sequences stored from the genome FASTA.
    *   Uses Python's highly optimized built-in `string.find()` method to search for the forward sequence *and* the reverse complement sequence within each genome sequence string.
    *   Returns the location details (chromosome, 1-based start/end, strand) of the *first* match found. This is significantly faster than alignment-based methods for finding only exact matches.
6.  **Gene Association (`main` function logic using `pyranges`):**
    *   Filters the search results to include only sgRNAs that were found in the genome.
    *   Converts the sgRNA match locations into a `pyranges.PyRanges` object.
    *   **Overlaps:** Uses `sgrna_pr.join(genes_pr)` to efficiently find all genes that overlap with each sgRNA's genomic location. Results are aggregated if an sgRNA overlaps multiple genes.
    *   **Nearest:**
        *   Identifies sgRNAs that did *not* overlap any gene.
        *   Uses `non_overlapping_pr.nearest(genes_pr, overlap=False)` to find the single closest gene feature (regardless of direction initially).
        *   Filters out results where no gene was found on the same chromosome (`Distance = -1`).
        *   **Post-processes** the nearest results by comparing the sgRNA coordinates (`Start`, `End`) with the nearest gene's coordinates (`Start_gene`, `End_gene`) to determine if the closest gene is upstream (`End_gene <= Start`) or downstream (`Start_gene >= End`).
        *   Separately aggregates the upstream and downstream nearest gene information.
7.  **Report Generation (`main` function logic):**
    *   Merges the initial search results (all sgRNAs) with the overlap results and the separate nearest upstream/downstream results using `pandas.merge`.
    *   Initializes and formats all output columns, ensuring consistent handling of missing data (represented as 'N/A' in the final CSV). Uses `pandas` nullable integer types (`Int64`) internally to handle missing numeric data correctly before final string conversion.
    *   Saves the final, ordered DataFrame to a CSV file.
8.  **Plotting (`main` function logic):** Uses `matplotlib` to generate a simple bar plot summarizing the count of found vs. not-found sgRNAs.

## Assumptions and Limitations

*   **Exact Matches Only:** This script finds only 100% identical sequence matches. It does not perform alignments and will not find matches with mismatches or gaps.
*   **Single Chromosome Assumption for Gene Table:** The script currently assumes all coordinates in the gene table file refer to the *first* sequence entry found in the provided genome FASTA file. If your genome has multiple chromosomes/contigs *and* your gene table contains coordinates relative to different sequences, the gene association results will be incorrect for genes not on that first sequence.
*   **Gene Table Format:** The script strictly expects the specified column headers in the gene table file. The coordinate columns (`Left-End-Position`, `Right-End-Position`) are assumed to be 1-based inclusive.
*   **Gene Table Strand:** The current implementation assumes the gene table lacks strand information and assigns a dummy '+' strand. Therefore, the determination of "upstream" and "downstream" is purely based on genomic coordinates, not gene orientation.
*   **Performance:** Exact string matching is generally fast. Reading large genomes and performing `pyranges` operations (join, nearest) are the potentially more time-consuming steps, but `pyranges` is highly optimized for these tasks.