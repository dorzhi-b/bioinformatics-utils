import os
import sys
import time
import datetime
import requests
from functools import wraps
from io import StringIO
from dotenv import load_dotenv
from typing import Dict, Tuple, Union
from Bio import SeqIO
from Bio.SeqUtils import GC
from abc import ABC, abstractmethod


load_dotenv() 

TG_API_TOKEN = os.getenv("TG_API_TOKEN") 

class BiologicalSequence(ABC):
    """
    Abstract base class representing a biological sequence.
    """

    @abstractmethod
    def __len__(self):
        """Return the length of the biological sequence."""
        pass

    @abstractmethod
    def __getitem__(self, index):
        """Return the item at the specified index in the sequence."""
        pass
    
    def __str__(self):
        """Return the string representation of the sequence."""
        return self.sequence
    
    @abstractmethod
    def is_valid_alphabet(self):
        """Check if the sequence contains valid alphabet characters."""
        pass
    
    @abstractmethod
    def __repr__(self):
        """Return the canonical string representation of the sequence."""
        pass

class NucleicAcidSequence(BiologicalSequence):
    """
    Abstract base class representing a nucleic acid sequence.
    """

    @abstractmethod
    def complement_dict(self):
        """Return the complement dictionary for the sequence."""
        pass

    @abstractmethod
    def valid_nucleotides(self):
        """Return the set of valid nucleotides for the sequence."""
        pass

    def __init__(self, sequence):
        """Initialize the nucleic acid sequence."""
        self.sequence = sequence.upper()
    
    def __len__(self):
        """Return the length of the nucleic acid sequence."""
        return len(self.sequence)
    
    def __getitem__(self, index):
        """Return the item at the specified index in the sequence."""
        return self.sequence[index]
    
    def is_valid_alphabet(self):
        """Check if the sequence contains valid nucleotide alphabet characters."""
        return all(base in self.valid_nucleotides() for base in self.sequence)
    
    def complement(self):
        """Return the complement of the sequence."""
        return ''.join(self.complement_dict()[base] for base in self.sequence)
    
    @abstractmethod
    def gc_content(self):
        """Calculate the GC content of the sequence."""
        pass
    
    def __repr__(self):
        """Return the canonical string representation of the sequence."""
        return f"{self.__class__.__name__}('{self.sequence}')"

class DNASequence(NucleicAcidSequence):
    """
    Class representing a DNA sequence.
    """

    def complement_dict(self):
        """Return the complement dictionary for DNA."""
        return {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    
    def valid_nucleotides(self):
        """Return the set of valid nucleotides for DNA."""
        return {'A', 'T', 'G', 'C'}
    
    def gc_content(self):
        """Calculate the GC content of the DNA sequence."""
        gc_count = sum(base in {'G', 'C'} for base in self.sequence)
        return gc_count / len(self.sequence)
    
    def transcribe(self):
        """Transcribe the DNA sequence into an RNA sequence."""
        return RNASequence(self.sequence.replace('T', 'U'))

class RNASequence(NucleicAcidSequence):
    """
    Class representing an RNA sequence.
    """

    def complement_dict(self):
        """Raise an error as RNA complement is not defined."""
        raise NotImplementedError("RNA complement is not defined.")
    
    def valid_nucleotides(self):
        """Return the set of valid nucleotides for RNA."""
        return {'A', 'U', 'G', 'C'}
    
    def gc_content(self):
        """Raise an error as GC content calculation is not defined for RNA."""
        raise NotImplementedError("GC content calculation is not defined for RNA.")

class AminoAcidSequence(BiologicalSequence):
    """
    Class representing an amino acid sequence.
    """

    def __init__(self, sequence):
        """Initialize the amino acid sequence."""
        self.sequence = sequence
    
    def __len__(self):
        """Return the length of the amino acid sequence."""
        return len(self.sequence)
    
    def __getitem__(self, index):
        """Return the item at the specified index in the sequence."""
        return self.sequence[index]
    
    def is_valid_alphabet(self):
        """Check if the sequence contains valid amino acid characters."""
        valid_amino_acids = {'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'}
        return all(aa in valid_amino_acids for aa in self.sequence)
    
    def __repr__(self):
        """Return the canonical string representation of the sequence."""
        return f"{self.__class__.__name__}('{self.sequence}')"

    def count_hydrophobic_residues(self):
        """Count the number of hydrophobic residues in the sequence."""
        hydrophobic_aa = {'A', 'V', 'L', 'I', 'P', 'F', 'W', 'M'}
        return sum(aa in hydrophobic_aa for aa in self.sequence)

def filter_fastq(file_path: str,
                 gc_bounds: Union[float, Tuple[Union[float, int], Union[float, int]]] = (0, 100),
                 length_bounds: Union[int, Tuple[int, int]] = (0, 2**32),
                 quality_threshold: Union[float, int] = 0) -> Dict[str, Tuple[str, str]]:
    """
    Filters sequences based on GC-content, length, and quality threshold.

    Parameters:
    -----------
    file_path (str):
        Path to the FastQ file.
    gc_bounds (Union[float, Tuple[Union[float, int], Union[float, int]]]):
        GC-content filtering bounds.
        - Default is 100 (no filtering).
        - If a single float is provided, it's considered as the upper bound.
        - If an integer is provided, it's considered as the upper bound.
        - If a tuple of two numbers is provided, they are treated as the lower and upper bounds.
    length_bounds (Union[int, Tuple[int, int]]):
        Length filtering bounds.
        - Default is (0, 2^32) (no filtering).
        - If a single integer is provided, it's considered as the upper bound.
        - If a tuple of two integers is provided, they are treated as the lower and upper bounds.
    quality_threshold (Union[float, int]):
        Quality threshold value.
        - Default is 0 (no filtering).
        - Sequences with an average quality below the threshold are discarded.

    Returns:
    --------
    Dict[str, Tuple[str, str]]:
        A filtered dictionary where each value is a tuple containing two strings,
        representing DNA sequences and their quality.
    """
    filtered_seqs = {}

    with open(file_path, "r") as handle:
        for record in SeqIO.parse(handle, "fastq"):
            gc_content = GC(record.seq)

            if isinstance(gc_bounds, (float, int)):
                lower_bound = 0
                upper_bound = gc_bounds
            elif isinstance(gc_bounds, tuple) and len(gc_bounds) == 2:
                lower_bound, upper_bound = gc_bounds
            else:
                raise ValueError("gc_bounds should be a single float/int or a tuple of two floats/ints.")

            if lower_bound <= gc_content <= upper_bound:
                seq_length = len(record.seq)

                if isinstance(length_bounds, int):
                    lower_bound = 0
                    upper_bound = length_bounds
                elif isinstance(length_bounds, tuple) and len(length_bounds) == 2:
                    lower_bound, upper_bound = length_bounds
                else:
                    raise ValueError("length_bounds should be a single integer or a tuple of two integers.")

                if lower_bound <= seq_length <= upper_bound:
                    quality_values = [chr(qual) for qual in record.letter_annotations["phred_quality"]]
                    quality_count = sum(record.letter_annotations["phred_quality"])
                    seq_quality = quality_count / seq_length

                    if seq_quality >= quality_threshold:
                        filtered_seqs[record.id] = (str(record.seq), "".join(quality_values))

    return filtered_seqs

def send_telegram_message(chat_id, message):
    """"
    Sends a message to Telegram.

    :param chat_id: Chat ID where the message will be sent.
    :param message: The text of the message.
    :return: Response object from the request.
    """

    token = os.getenv("TG_API_TOKEN")
    url = f"https://api.telegram.org/bot{token}/sendMessage"
    data = {
        "chat_id": chat_id,
        "text": message,
        "parse_mode": "Markdown"
    }
    response = requests.post(url, data=data)
    return response

def telegram_logger(chat_id):
    """
    Logging the execution of functions with notifications sent to Telegram.

    :param chat_id: Chat ID where the messages will be sent.
    :return: Decorator wrapper for the function, adding logging and notifications.
    """

    def decorator(func):
        @wraps(func)
        def wrapper(*args, **kwargs):
            old_stdout, old_stderr = sys.stdout, sys.stderr
            stdout = StringIO()
            stderr = StringIO()
            sys.stdout, sys.stderr = stdout, stderr

            start_time = datetime.datetime.now()
            try:
                result = func(*args, **kwargs)
                execution_time = datetime.datetime.now() - start_time
                if execution_time.days < 1:
                    execution_time_str = (datetime.datetime.min + execution_time).strftime('%H:%M:%S.%f')[:-3]
                else:
                    execution_time_str = f"{execution_time.days} days, " + (datetime.datetime.min + execution_time).strftime('%H:%M:%S')
                send_telegram_message(chat_id, f"`Function '{func.__name__}' executed successfully in {execution_time_str}.`")
            except Exception as e:
                send_telegram_message(chat_id, f"`Function '{func.__name__}' encountered an error: {type(e).__name__}: {str(e)}`")
                raise
            finally:
                sys.stdout, sys.stderr = old_stdout, old_stderr 
                stdout_value = stdout.getvalue()
                stderr_value = stderr.getvalue()
                if stdout_value:
                    send_telegram_message(chat_id, f"`==== {func.__name__} stdout log ====`\n\n{stdout_value}")
                if stderr_value:
                    send_telegram_message(chat_id, f"`==== {func.__name__} stderr log ====`\n\n{stderr_value}")
                stdout.close()
                stderr.close()
                
            return result
        return wrapper
    return decorator
