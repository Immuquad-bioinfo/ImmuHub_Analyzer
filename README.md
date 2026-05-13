ImmuHub Analyzer is a lightweight immune repertoire analysis toolkit built on top of the NCBI IgBLAST (v1.22) engine and IMGT germline reference databases for TCR/BCR repertoire sequence analysis. ￼

The project serves as a modern wrapper around IgBLAST while preserving the original alignment logic and result consistency of the official IgBLAST implementation. It provides:

* Simplified command-line interface
* Automatic input format detection
* Automatic FASTQ → FASTA conversion
* Automatic germline database management
* AIRR standard output support
* Unified multi-chain TCR/BCR analysis
* High-throughput repertoire analysis workflows

ImmuHub Analyzer supports:

BCR Chains

* IGH
* IGK
* IGL

TCR Chains

* TRA
* TRB
* TRG
* TRD

Supported input formats:

* FASTA
* FASTQ / FASTQ.GZ
* CDR3 tables (TSV/CSV)

Supported output formats:

* AIRR TSV
* IgBLAST TXT
* CSV

⸻

Features

1. Built on the Original IgBLAST Engine

ImmuHub Analyzer does not reimplement V(D)J alignment algorithms. Instead, it directly wraps and invokes:

* igblastn
* igblastp
* makeblastdb

This ensures:

* Full compatibility with official IgBLAST behavior
* Publication-grade and clinical research consistency
* AIRR-compliant outputs

⸻

2. IMGT Germline Database Support

By default, ImmuHub Analyzer supports:

* IMGT germline references
* Human (preset)
* Mouse and other species

Automatic support for:

* V
* D
* J

germline databases.

⸻

3. Automatic FASTQ Processing

Users do not need to:

* Manually convert FASTQ to FASTA
* Decompress gzip files

ImmuHub Analyzer automatically performs:

FASTQ(.gz)↓FASTA↓IgBLAST↓AIRR

⸻

4. AIRR Standard Output

Supports:

--format airr

Outputs AIRR Rearrangement TSV format.

Compatible with downstream analysis tools such as:

* MiXCR
* VDJtools
* immuneML
* scRepertoire
* immunarch

⸻

5. Full IgBLAST Argument Compatibility

Supports:

--extra

to pass through native IgBLAST arguments.

Example:

--extra -domain_system imgt

Users may continue using original IgBLAST options such as:

* domain_system
* show_translation
* num_alignments_V
* evalue
* penalty

and many others.

⸻

Installation

1. Download the Project

Download and extract the package under Linux:

cd ImmuHub_Analyzer

⸻

2. Install the Python Package

python3 -m pip install -e .

⸻

Database Preparation

1. Download IMGT Germline FASTA Files

Example:

database/mouse_IGHV.fastadatabase/mouse_IGHD.fastadatabase/mouse_IGHJ.fasta

Recommended source:

IMGT Reference Directory￼

⸻

2. Build BLAST Databases

makeblastdb \ -parse_seqids \ -dbtype nucl \ -blastdb_version 4 \ -in mouse_IGHV.clean.fasta \ -out mouse_IGHV

⸻

Usage

FASTQ Analysis

IGH Example

immuhub-analyzer run \ -i sample.fastq.gz \ -o sample.airr.tsv \ --locus IGH \ --v-db database/human_IGHV \ --d-db database/human_IGHD \ --j-db database/human_IGHJ \ --aux-file third_party/optional_file/human_gl.aux \ --threads 32 \ --format airr

⸻

TRB Example

immuhub-analyzer run \ -i sample.fastq.gz \ -o sample.airr.tsv \ --locus TRB \ --v-db database/human_TRBV \ --d-db database/human_TRBD \ --j-db database/human_TRBJ \ --aux-file third_party/optional_file/human_gl.aux \ --threads 32 \ --format airr

⸻

Project Goals

The goals of ImmuHub Analyzer are to:

* Simplify IgBLAST workflows
* Provide modern repertoire analysis pipelines
* Improve usability for TCR/BCR analysis
* Maintain full compatibility with official IgBLAST
* Support both clinical and research applications

⸻

License

This project wraps and depends on:

* NCBI IgBLAST
* IMGT germline references

Please also comply with:

* NCBI IgBLAST license
* IMGT usage policy

⸻

Citation

If you use ImmuHub Analyzer in publications, please cite:

* NCBI IgBLAST
* IMGT database
* ImmuHub Analyzer
