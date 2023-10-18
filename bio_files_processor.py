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
