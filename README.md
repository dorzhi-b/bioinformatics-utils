# bioinformatics-utils

![Python Version](https://img.shields.io/badge/python-3.x-blue.svg)

This repository contains a collection of bioinformatics tools and scripts for working with DNA/RNA sequences and protein sequences. These tools can be used for various tasks such as sequence manipulation, filtering, and analysis. Also this package provides functionalities to process biological files: converting a multi-line FASTA file into a one-line FASTA file and shifting the start position of the sequence in a FASTA file by a specified number of positions.

## Table of Contents

- [Installation](#installation)
- [Usage](#usage)
- [Author](#author)

## Installation

Clone this repository to your local machine or download it manualy:

   ```bash
   git clone git@github.com:dorzhi-b/bioinformatics-utils.git
   ```

## Usage
### DNA/RNA Sequence Tools

run_dna_rna_tools.py
This script allows you to perform various operations on DNA/RNA sequences, including transcription, reversal, complementation, and reverse complementation.
Requires one or a list of nucleic acid sequences, and as the last element the command.

Example usage:

```python
run_dna_rna_tools("ATGCTA", "transcribe")
"AUGCUA"
```

```python
run_dna_rna_tools("ATGCTA", "TTAGCG", "reverse")
["ATCGTA", "GCGATT"]
```

```python
run_dna_rna_tools("ATGCTA", "complement")
"TACGAT"
```

```python
run_dna_rna_tools("ATGCTA", "reverse_complement")
"TAGCAT"
```

### FASTQ Data Filtering

filter_fastq.py
This script provides a way to filter FASTQ data based on GC-content, sequence length, and quality threshold.
Takes as input a path to FASTQ file and converts the file to dictionary of the structure key - string, sequence name; value is a tuple of two strings: sequence and quality. Return a similar dictionary consisting only of those sequences that pass all the conditions. As output save filtered dictionary as new FASTQ file.

Example usage:

```python
 #Filtering sequences with GC-content upper bound of 40%, length between 50 and 100 and quality threshold of 20   
filter_fastq(input_path, gc_bounds = 40, length_bounds = (50, 100), quality_threshold = 20)
```

### Protein Sequence Tools

protein_tool.py
This script allows you to manipulate and analyze protein sequences, including converting between one-letter and three-letter amino acid codes, converting to RNA/DNA, counting charges, and counting polar/nonpolar amino acids. 


Example usage:

```python
protein_tool('ASDRKHDE', 'charge')
{'Positive': 3, 'Negative': 3, 'Neutral': 2}
```

```python
protein_tool('FM', 'RNA')
'UUYAUG'
```

```python
protein_tool('ASDR', 'polarity')
[{'Polar': 3, 'Nonpolar': 1}]
```

```python
protein_tool('AsDr', 'DNA')
'GCN(TCN or AGY)GAYAGY'
protein_tool('YWNGAS', 'DNA')
'TAY(CGN or AGR)AAYGGNGCN(TCN or AGY)'
```

```python
protein_tool('aLa-CyS', 'one letter') #input ignore letter's size
'AC'
protein_tool('Ala-Cys', 'Ala', 'one letter')
['AC', 'A']
```

### FASTQ Data Filtering

bio_files_processor.py
Converts a multi-line FASTA file into a one-line FASTA file and shifts the start position of the sequence in a FASTA file by a specified number of positions.

Example usage:

```python
 from bio_files_processor import convert_multiline_fasta_to_oneline

input_fasta = "input_multiline.fasta"  # Specify the path to the multi-line fasta file
output_fasta = "output_oneline.fasta"  # Specify the path to save the one-line fasta file

convert_multiline_fasta_to_oneline(input_fasta, output_fasta)
```

```python
from bio_files_processor import change_fasta_start_pos

input_fasta = "input.fasta"  # Specify the existing path to the fasta file
output_fasta = "output.fasta"  # Specify the path to save the output file

change_fasta_start_pos(input_fasta, 2, output_fasta)
```

## Author
- Dorzhi Badmadashiev
