from typing import Union, List, Dict, Tuple
from bioinformatics_toolkit import dna_rna_tools
from bioinformatics_toolkit import filter_fastq
from bioinformatics_toolkit import protein_tool

ALPHABET = {'a', 't', 'u', 'g', 'c', 'A', 'T', 'U', 'G', 'C'}

def run_dna_rna_tools(*args: Union[str, List[str]], ALPHABET: set = ALPHABET) -> Union[str, List[str]]:
    """
    Perform the specified operation on the given DNA/RNA sequences.
    
    This function performs the following operations: transcribtion, reversation, complementation, 
    reversal complementation.
    
    Parameters:
        *args (Union[str, List[str]]): DNA/RNA sequence(s).
        ALPHABET (set): Default set of valid nucleotide characters.
    
    Returns:
        Union[str, List[str]]: The result of the specified operation on the sequence(s).
            If multiple sequences are provided, a list of results is returned.

    Examples:
        # Transcribing a single DNA sequence
        run_dna_rna_tools("ATGCTA", "transcribe")
        # Result: "AUGCUA"

        # Transcribing multiple DNA sequences
        run_dna_rna_tools("ATGCTA", "TTAGCG", "transcribe")
        # Result: ["AUGCUA", "UUAGCG"]
    """
    operations = {
        'transcribe': transcribe,
        'reverse': reverse,
        'complement': complement,
        'reverse_complement': reverse_complement
    }
    
    seqs, operation = args[:-1], args[-1]
    
    if  not all(base in ALPHABET for seq in seqs for base in seq):
        return "Not DNA/RNA. Use only standard nucleotide sequences."
    if not operation in operations:
        return "Invalid operation. Supported operations: transcribe, reverse, complement, reverse_complement."
    results = [operations[operation](seq) for seq in seqs]
    if len(results) == 1:
        return results[0]
    else:
        return results

def filter_fastq(seqs: Dict[str, Tuple[str, str]],
                 gc_bounds: Union[float, Tuple[Union[float, int], Union[float, int]]] = (0, 100),
                 length_bounds: Union[int, Tuple[int, int]] = (0, 2**32),
                 quality_threshold: Union[float, int] = 0) -> Dict[str, Tuple[str, str]]:
    """
    Filters sequences based on GC-content, length, and quality threshold.

    Parameters:
    -----------
    ### seqs (Dict[str, Tuple[str, str]]):
        A dictionary where each value is a tuple containing two strings,
        representing DNA sequences and their quality.
    ### gc_bounds (Union[float, Tuple[Union[float, int], Union[float, int]]]):
        GC-content filtering bounds.
        - Default is 100 (no filtering).
        - If a single float is provided, it's considered as the upper bound.
        - If an integer is provided, it's considered as the upper bound.
        - If a tuple of two numbers is provided, they are treated as the lower and upper bounds.
    ### length_bounds (Union[int, Tuple[int, int]]):
        Length filtering bounds.
        - Default is (0, 2^32) (no filtering).
        - If a single integer is provided, it's considered as the upper bound.
        - If a tuple of two integers is provided, they are treated as the lower and upper bounds.
    ### quality_threshold (Union[float, int]):
        Quality threshold value.
        - Default is 0 (no filtering).
        - Sequences with an average quality below the threshold are discarded.

    Returns:
    --------
    Dict[str, Tuple[str, str]]:
        A filtered dictionary where each value is a tuple containing two strings,
        representing DNA sequences and their quality.

    Examples:
    ---------
    ### Filtering sequences with a GC-content upper bound of 40%
    filtered_seqs = filter_fastq(seqs, gc_bounds = 40)

    ### Filtering sequences with a length between 50 and 100
    filtered_seqs = filter_fastq(seqs, length_bounds = (50, 100))

    ### Filtering sequences with a quality threshold of 20
    filtered_seqs = filter_fastq(seqs, quality_threshold = 20)
    """
                   
    filtered_seqs = filter_gc_content(seqs, gc_bounds)
    filtered_seqs = filter_length(filtered_seqs, length_bounds)
    filtered_seqs = filter_quality(filtered_seqs, quality_threshold)
                   
    return filtered_seqs

def protein_tool(*args: str) -> Union[str, List[Union[Dict[str, int], str]]]:
    """
    Receives a request from the user and runs the desired function.

    Parameters
    ----------
    *args : str
        Amino acid sequences and operation type.

    Returns
    -------
    str
        If a single sequence is supplied, outputs the result as a string or
        identifies a problem with a specific sequence.
    list
        If several sequences are supplied, outputs the result as a list.
    """
    *seqs, operation = args
    operations = {
        'one letter': change_abbreviation, 'RNA': to_rna, 'DNA': to_dna,
        'charge': define_charge, 'polarity': define_polarity
    }
    output = []
    for seq in seqs:
        answer = is_correct_seq(seq.upper())
        if answer:
            function_output = operations[operation](seq.upper())
            output.append(function_output)
        else:
            print(f'Something wrong with {seq}', file=sys.stderr)
            continue
    if len(output) == 1 and (operation == 'RNA' or operation == 'DNA' or operation == 'one letter'):
        return ''.join(output)
    else:
        return output
