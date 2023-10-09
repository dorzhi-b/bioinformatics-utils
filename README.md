# bioinformatics-utils

![Python Version](https://img.shields.io/badge/python-3.x-blue.svg)

This repository contains a collection of bioinformatics tools and scripts for working with DNA/RNA sequences and protein sequences. These tools can be used for various tasks such as sequence manipulation, filtering, and analysis.

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
Requires one or a list of nucleic acid sequences, and as the last element the command
Example usage:

```python
run_dna_rna_tools("ATGCTA", "transcribe")
"AUGCUA"
```

```python
run_dna_rna_tools("ATGCTA", "TTAGCG", "transcribe")
["AUGCUA", "UUAGCG"]
```

### FASTQ Data Filtering

filter_fastq.py
This script provides a way to filter FASTQ data based on GC-content, sequence length, and quality threshold.
Takes as input a dictionary of the structure key - string, sequence name; value is a tuple of two strings: sequence and quality. Return a similar dictionary consisting only of those sequences that pass all the conditions.

Example usage:

```python
 #Filtering sequences with GC-content upper bound of 40%, length between 50 and 100 and quality threshold of 20   
filter_fastq(seqs, gc_bounds = 40, length_bounds = (50, 100), quality_threshold = 20)
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


## Author
- Dorzhi Badmadashiev
