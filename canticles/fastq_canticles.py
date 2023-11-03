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
