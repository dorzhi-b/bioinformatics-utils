from dataclasses import dataclass


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




@dataclass
class FastaRecord:
    """
    Represents a FASTA record with an ID, description, and sequence.
    """

    id: str
    description: str
    seq: str

    def __repr__(self):
        return f"id = {self.id}, description = {self.description}, seq = {self.seq}"


class OpenFasta:
    """
    Context manager for opening and reading a FASTA file.
    """

    def __init__(self, file_path):
        self.file_path = file_path

    def __enter__(self):
        self.file = open(self.file_path)
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.file.close()

    def __next__(self):
        record = self.read_record()
        if record:
            return record
        else:
            raise StopIteration

    def __iter__(self):
        record = self.read_record()
        while record:
            yield record
            record = self.read_record()

    def read_record(self) -> FastaRecord:
        """
        Read and return the next FASTA record from the file.

        Returns:
            FastaRecord: The parsed FASTA record.
        """

        record_id = ''
        description = ''
        sequence = ''

        line = self.file.readline()
        while line:
            if line.startswith('>'):
                if record_id:
                    return FastaRecord(record_id, description.strip(), sequence.strip())
                else:
                    parts = line.strip().split(maxsplit=1)
                    record_id = parts[0][1:]
                    description = parts[1] if len(parts) > 1 else ''
            else:
                sequence += line.strip()
            line = self.file.readline()

        if record_id:  
            return FastaRecord(record_id, description.strip(), sequence.strip())
        else:
            return None  

    def read_records(self) -> list:
        """
        Read and return all FASTA records from the file.

        Returns:
            list: A list of parsed FastaRecord objects.
        """
        records = []
        record_id, description, sequence = '', '', ''
        for line in self.file:
            line = line.strip()
            if line.startswith('>'):
                if record_id:
                    records.append(FastaRecord(record_id, description.strip(), sequence.strip()))
                    record_id, description, sequence = '', '', ''
                parts = line[1:].split(maxsplit=1)
                record_id = parts[0]
                description = parts[1] if len(parts) > 1 else ''
            else:
                sequence += line

        if record_id:
            records.append(FastaRecord(record_id, description.strip(), sequence.strip()))

        return records
