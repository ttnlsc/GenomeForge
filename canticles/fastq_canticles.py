import os


def calculate_gc_content(sequence: str) -> float:
    """
    Calculates the GC content (Guanine-Cytosine content) of the input sequence.

    Args:
    - sequence (str): The input sequence.

    Returns:
    - float: The calculated GC content in percent.
    """
    gc_count = 0

    for nucleotide in sequence:
        if nucleotide in 'GCgc':
            gc_count += 1

    return (gc_count / len(sequence)) * 100


def calculate_seq_length(sequence: str) -> int:
    """
    Calculate the length of a DNA, RNA, or protein sequence.

    Args:
    - sequence (str): A string representing a sequence.

    Returns:
    - int: The length of the input sequence.
    """
    return len(sequence)


def calculate_seq_quality(quality: str) -> float:
    """
    Calculate the mean quality score for a sequence quality encoded as ASCII characters (Phred+33),
    where the ASCII character is converted to a quality score by subtracting 33 from the ASCII code.

    Args:
    - quality (str): A string containing quality values.

    Returns:
    - float: The mean quality score calculated as the sum of quality values divided by the length of the input string.
    """
    score = 0

    for char in quality:
        qual = ord(char) - 33
        score += qual

    mean_quality = score / len(quality)
    return mean_quality


def parse_fastq(path_to_seqs: str) -> dict:
    """
    Parse Fastq sequences from a file and store them in a dictionary.

    Args:
    - path_to_seqs (str): The path to the Fastq file to parse.

    Returns:
    - dict: A dictionary containing Fastq sequences in the format:
        {'name': ('sequence', 'quality')}
    """
    fastq_dict = {}
    with open(path_to_seqs, mode='r') as file:
        fastq_file = file.readlines()
    record = None

    for line in fastq_file:
        line = line.strip()
        if not record:
            record = [line]
        else:
            record.append(line)

        if len(record) == 4:
            name = record[0]
            sequence = record[1]
            quality = record[3]
            fastq_dict[name] = (sequence, quality)
            record = None
            
    return fastq_dict


def write_filtered_fastq(filtered_seqs: dict, output_file_name: str) -> None:
    """
    Writes filtered Fastq sequences to a file.

    Args:
    - filtered_seqs (dict): A dictionary containing Fastq sequences in the format:
        {'name': ('sequence', 'quality')}
    - output_file_name (str): The name of the output file to write the sequences.

    Returns:
    - None
    """
    output_folder = './fastq_filtrator_results'
    
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    
    record = []
    
    for name, (sequence, quality) in filtered_seqs.items():
        third_line = ''.join(['+', name[1:]])
        record.append(name)
        record.append(sequence)
        record.append(third_line)
        record.append(quality)
    path_to_output = os.path.join(output_folder, output_file_name)
    
    with open(path_to_output, mode='w') as file:
        file.writelines('\n'.join(record))
    return None
