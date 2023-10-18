from typing import List, Union

COMPLEMENT_DICT_DNA = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'a': 't', 't': 'a', 'c': 'g', 'g': 'c'}
COMPLEMENT_DICT_RNA = {'A': 'U', 'U': 'A', 'C': 'G', 'G': 'C', 'a': 'u', 'u': 'a', 'c': 'g', 'g': 'c'}
ALPHABET = {'a', 't', 'u', 'g', 'c', 'A', 'T', 'U', 'G', 'C'}

def transcribe(seq: str) -> Union[str, List[str]]:
    """
    Transcribe a DNA sequence to RNA.
    
    Parameters:
        seq (str): The DNA sequence to transcribe.
    
    Returns:
        Union[str, List[str]]: The transcribed RNA sequence.
            If multiple sequences are provided, a list of transcribed sequences is returned.
    """
    if 'U' in seq or 'u' in seq:
        return "Invalid sequence. Contains 'U' in DNA sequence."
    
    return seq.replace('T', 'U').replace('t', 'u')

def reverse(seq: str) -> str:
    """
    Reverse a DNA/RNA sequence.
    
    Parameters:
        seq (str): The DNA/RNA sequence to reverse.
    
    Returns:
        str: The reversed DNA/RNA sequence.
    """
    seq_u = seq.upper()
    if 'T' in seq_u and 'U' in seq_u:
        return "Invalid sequence. Contains both 'T' and 'U'."
    
    return seq[::-1]

def complement(seq: str) -> str:
    """
    Find the complement of a DNA sequence.
    
    Parameters:
        seq (str): The DNA sequence to find the complement for.
    
    Returns:
        str: The complement DNA sequence.
    """
    seq_u = seq.upper()
    if 'T' in seq_u and 'U' in seq_u:
        return "Invalid sequence. Contains both 'T' and 'U'."
    
    if 'T' in seq or "t" in seq:
        return ''.join(COMPLEMENT_DICT_DNA[base] for base in seq)
    else:
        return ''.join(COMPLEMENT_DICT_RNA[base] for base in seq)

def reverse_complement(seq: str) -> str:
    """
    Find the reverse complement of a DNA sequence.
    
    Parameters:
        seq (str): The DNA sequence to find the reverse complement for.
    
    Returns:
        str: The reverse complement DNA sequence.
    """
    return reverse(complement(seq))

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
