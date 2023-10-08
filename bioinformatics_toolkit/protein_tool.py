from typing import Dict, List, Union

# Dorzhi
def to_rna(seq: str, rna_dict: Dict[str, str] = {'F': 'UUY', 'L': 'YUN', 'I': 'AUH', 'M': 'AUG',
                              'V': 'GUN', 'S': 'WSN', 'P': 'CCN', 'T': 'ACN',
                              'A': 'GCN', 'Y': 'UAY', 'H': 'CAY', 'Q': 'CAR',
                              'N': 'AAY', 'K': 'AAR', 'D': 'GAY', 'E': 'GAR',
                              'C': 'UGY', 'R': 'MGN', 'G': 'GGN', 'W': 'UGG'}) -> str:
  """
  Converts an amino acid sequence into an RNA sequence.

  Parameters
  ----------
  seq : str
      Amino acid sequence.
  rna_dict : dict
      Dictionary defining the correspondence of amino acids 
      to RNA triplets (default, standard code).
  Returns
  -------
  str
      RNA sequence.

  """
  result = ''.join(rna_dict[base] for base in seq)
  return result


def define_charge(seq: str, positive_charge: List[str] = ['R', 'K', 'H'],
                  negative_charge: List[str] = ['D', 'E']) -> Dict[str, int]:
  """
  Counts the number of amino acids with positive charge, negative charge, 
  and neutral amino acids in the sequence.

  Parameters
  ----------
  seq : str
      Amino acid sequence (string).
  positive_charge : list
      List of amino acids with positive charge (default is ['R', 'K', 'H']).
  negative_charge : list
      List of amino acids with negative charge (default is ['D', 'E']).

  Returns
  -------
  dict
      A dictionary containing the counts of amino acids and their labels:
      - 'Positive' for amino acids with positive charge.
      - 'Negative' for amino acids with negative charge.
      - 'Neutral' for neutral amino acids.
  """
  positive_count = 0
  negative_count = 0
  neutral_count = 0

  for aa in seq:
    if aa in positive_charge:
      positive_count += 1
    elif aa in negative_charge:
      negative_count += 1
    else:
      neutral_count += 1
    
    result = {
        'Positive': positive_count,
        'Negative': negative_count,
        'Neutral': neutral_count
    }
  return result


# Ustin
POLAR_AA = {'D', 'E', 'R', 'K', 'H', 'N', 'Q', 'S', 'T', 'Y', 'C'}
NONPOLAR_AA = {'A', 'G', 'V', 'L', 'I', 'P', 'F', 'M', 'W'}
DNA_AA = {'F': 'TTY', 'L': '(TTR or CTN)', 'I': 'ATH', 'M': 'ATG', 'V': 'GTN', 'S': '(TCN or AGY)', 'P': 'CCN', 'T': 'ACN', 'A': 'GCN',
          'Y': 'TAY', 'H': 'CAY', 'Q': 'CAR', 'N': 'AAY', 'K': 'AAR', 'D': 'GAY', 'E': 'GAR', 'C': 'TGY', 'W': '(CGN or AGR)', 'R': 'AGY', 'G': 'GGN'}


def define_polarity(seq: str) -> Dict[str, int]:
    """
    Counts polar and nonpolar aminoacids in aminoacid sequences.

    Arguments:
     str: sequence to count polar and nonpolar aminoacids.

    Return:
     Dict[str, int]: 
          Dictionary with keys 'Polar', 'Nonpolar' and values of quantity of according groups in sequence. 
    """
    polarity_count = {'Polar': 0, 'Nonpolar': 0}
    for aminoacid in seq:
        if aminoacid in POLAR_AA:
            polarity_count['Polar'] += 1
        else:
            polarity_count['Nonpolar'] += 1
    return polarity_count


def to_dna(seq: str) -> str:
    """
    Transforms aminoacid sequence to DNA sequence

    Arguments
    ---------
     str: aminoacid sequence to transform to DNA sequence.

    Return
    ------
     str: according DNA sequence.
    """
    sequence_dna = []
    for aminoacid in seq:
        sequence_dna.append(DNA_AA[aminoacid])
    return ''.join(sequence_dna)


#Margarita
ABBREVIATION_THREE_TO_ONE = {'ALA':'A', 'CYS':'C', 'ASP':'D', 'GLU':'E', 'PHE':'F',
                             'GLY':'G', 'HIS':'H', 'ILE':'I', 'LYS':'K', 'LEU':'L',
                             'MET':'M', 'ASN':'N', 'PRO':'P', 'GLN':'Q', 'ARG':'R',
                             'SER':'S', 'TRE':'T', 'VAL':'V', 'TRP':'W', 'TYR':'Y'}
AMINO_ACIDS_ONE_LETTER = {'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K',
                          'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V',
                          'W', 'Y'}
AMINO_ACIDS_THREE_LETTER = {'ALA', 'CYS', 'ASP', 'GLU', 'PHE',
                             'GLY', 'HIS', 'ILE', 'LYS', 'LEU',
                             'MET', 'ASN', 'PRO', 'GLN', 'ARG',
                             'SER', 'TRE', 'VAL', 'TRP', 'TYR'}

import sys

def change_abbreviation(seq: str) -> str:
  """
  Changes the amino acid abbreviation from three-letter to one-letter.

  Parametrs
  ----------
  seq : str
      Amino acid sequence in three-letter form.
  Returns
  -------
  str
      Amino acid sequence in one-letter form
  
  """
  one_letter_seq = [ABBREVIATION_THREE_TO_ONE[amino_acid] for amino_acid in seq.split("-")]
  return "".join(one_letter_seq)

def is_correct_seq(seq: str) -> bool:
  """
  Check the sequence for extraneous characters.

  Parametrs
  ----------
  seq : str
      Amino acid sequence.
  Returns
  -------
  bool
      TRUE - if there is no extraneous characters, FALSE - if there is extraneous characters.
  
  """
  unique_amino_acids = set(seq)
  unique_amino_acids_three = set(seq.split("-"))
  check = unique_amino_acids <= AMINO_ACIDS_ONE_LETTER or unique_amino_acids_three <= AMINO_ACIDS_THREE_LETTER
  return check

def protein_tool(*args: str) -> Union[str, List[Union[Dict[str, int], str]]]:
  """
  Receives a request from the user and runs the desired function.

  Parametrs
  ----------
  seq : str
      Amino acid sequences.
  operation : str
      Type of user's request.
  Returns
  -------
  str
      If a single sequence is supplied, outputs the result as a string or or identify a problem with a specific sequence. 
  list
      If several sequences are supplied, outputs the result as a list.
  
  """
  *seqs, operation = args
  operations = {'one letter':change_abbreviation, 'RNA':to_rna, 'DNA':to_dna, 'charge':define_charge, 'polarity':define_polarity}
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
