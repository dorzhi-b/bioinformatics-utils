def convert_multiline_fasta_to_oneline(input_fasta: str, output_fasta: str = None):
    """
    Converts a multi-line FASTA file into a one-line FASTA file.

    Parameters:
    -----------
    input_fasta (str):
        The path to the input multi-line FASTA file.
    output_fasta (str, optional):
        The path to the output one-line FASTA file. If not provided, the default
        output filename will be created based on the input filename.

    """

    with open(input_fasta, 'r') as input_file:
        lines = input_file.readlines()

    if output_fasta is None:
        output_fasta = input_fasta.replace('.fasta', '_oneline.fasta')

    with open(output_fasta, 'w') as output_file:
        current_sequence = ''
        for line in lines:
            if line.startswith('>'):
                if current_sequence:
                    output_file.write(current_sequence + '\n')
                current_sequence = f'{line}'
            else:
                current_sequence += line.strip()
        if current_sequence:
            output_file.write(current_sequence + '\n')

def change_fasta_start_pos(input_fasta: str, shift: int, output_fasta: str = None):
    """
    Shifts the start position of the sequence in a FASTA file by a specified number of positions.

    Parameters:
    -----------
    input_fasta (str):
        The path to the input FASTA file.

    shift (int):
        The number of positions by which the start position of the sequence should be shifted. 
        A positive value shifts to the right, and a negative value shifts to the left.

    output_fasta (str, optional):
        The path to the output FASTA file. If not provided, the default output filename will be
        created based on the input filename by replacing '.fasta' with '_.fasta'.
    """
    
    if output_fasta is None:
        output_fasta = input_fasta.replace('.fasta', '_.fasta')

    with open(input_fasta, 'r') as f:
        data = f.read().strip().split('\n')
        header = data[0]
        sequence = data[1]

    shifted_sequence = sequence[shift:] + sequence[:shift]

    with open(output_fasta, 'w') as f:
        f.write(f"{header}\n{shifted_sequence}\n")
