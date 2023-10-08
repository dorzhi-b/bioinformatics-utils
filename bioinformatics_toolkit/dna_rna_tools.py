#!/usr/bin/env python
# coding: utf-8

# In[175]:


def transcribe(seq):
    """
    Transcribe a DNA sequence to RNA.
    
    Parameters
    ----------
    seq : str
        DNA sequence.

    Returns
    -------
    str
        Transcribed RNA sequence.
    """
    if 'U' and 'u' in seq:
        return "Can't transcribe RNA."
    
    return seq.replace('T', 'U').replace('t', 'u')

def reverse(seq):
    """
    Reverse a DNA/RNA sequence.
    
    Parameters
    ----------
    seq : str
        DNA/RNA sequence.

    Returns
    -------
    str
        Reversed DNA/RNA sequence.
    """
    seq_u = seq.upper()
    if 'T' in seq_u and 'U' in seq_u:
        return "Invalid sequence. Contains both 'T' and 'U'."
    
    return seq[::-1]

def complement(seq):
    """
    Find the complement of a DNAsequence.
    
    Parameters
    ----------
    seq : str
        DNA sequence.

    Returns
    -------
    str
        Complement DNA sequence.
    """
    seq_u = seq.upper()
    if 'T' in seq_u and 'U' in seq_u:
        return "Invalid sequence. Contains both 'T' and 'U'."
    
    if 'T' in seq or "t" in seq:
        complement_dict_dna = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'a': 't', 't': 'a', 'c': 'g', 'g': 'c'}
        return ''.join(complement_dict_dna[base] for base in seq)
    else:
        complement_dict_rna = {'A': 'U', 'U': 'A', 'C': 'G', 'G': 'C', 'a': 'u', 'u': 'a', 'c': 'g', 'g': 'c'}
        return ''.join(complement_dict_rna[base] for base in seq)
    

def reverse_complement(seq):
    """
    Find the reverse complement of a DNA sequence.
    
    Parameters
    ----------
    seq : str
        DNA sequence.

    Returns
    -------
    str
        Compliment reversed DNA sequence.
    """
    return reverse(complement(seq))

def run_dna_rna_tools(*args, alphabet = {'a', 't', 'u', 'g', 'c', 'A', 'T', 'U', 'G', 'C'}):
    """
    Perform the specified operation on the given DNA/RNA sequences.
    
    This function performs the following operations: transcribtion, reversation, complementation, 
    reversal complementation.
    
    Parameters
    ----------
    *args : str
        DNA/RNA sequence.
    alphabet : set
        Default set of valid nucleotide characters: {'a', 't', 'u', 'g', 'c', 'A', 'T', 'U', 'G', 'C'}

    Returns
    -------
    str
        Reversed DNA sequence.
    """
    ops = {
        'transcribe': transcribe,
        'reverse': reverse,
        'complement': complement,
        'reverse_complement': reverse_complement
    }
    
    seqs = args[:-1]
    op = args[-1]
    
    if all(base in alphabet for seq in seqs for base in seq):
        if op in ops:
            results = [ops[op](seq) for seq in seqs]
            if len(results) == 1:
                return results[0]
            else:
                return results
        else:
            return "Invalid operation. Supported operations: transcribe, reverse, complement, reverse_complement"
    else:
        return "Not DNA/RNA. Use only standard nucleotide sequences"

