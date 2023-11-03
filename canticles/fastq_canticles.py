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


def calculate_seq_quality(quality_str: str) -> float:
    """
    Calculate the mean quality score for a sequence quality encoded as ASCII characters (Phred+33),
    where the ASCII character is converted to a quality score by subtracting 33 from the ASCII code.

    Args:
    - quality_str (str): A string containing quality values.

    Returns:
    - float: The mean quality score calculated as the sum of quality values divided by the length of the input string.
    """
    score = 0
    
    for char in quality_str:
        quality = ord(char) - 33
        score += quality

    mean_quality = score / len(str)
    return mean_quality


def filter_fastq(seqs: dict, gc_bounds: tuple | int = (0, 100), length_bounds: tuple | int = (0, 0.2**32),
                 quality_threshold: int = 0) -> dict:
    """
    Filters fastq sequences based on the GC-content, length and quality parameters.

    Args:
    - seqs (dict): Dictionary of fastq sequences in the format {'name': ('sequence', 'quality')}
    - gc_bounds (tuple, int): GC content range (in percentages) for filtering
    - length_bounds (tuple, int): Length range for filtering
    - quality_threshold: Threshold value for average read quality filtering

    Returns:
    - dict: Dictionary with filtered reads
    """
    filtered_seqs = {}

    for name, (sequence, quality) in seqs.items():
        gc_content = calculate_gc_content(sequence)
        if isinstance(gc_bounds, tuple):
            if not (gc_bounds[0] <= gc_content <= gc_bounds[1]):
                continue
        elif gc_content >= gc_bounds:
            continue

        sequence_length = calculate_seq_length(sequence)
        if isinstance(length_bounds, tuple):
            if not (length_bounds[0] <= sequence_length <= length_bounds[1]):
                continue
        elif sequence_length >= length_bounds:
            continue

        mean_quality = calculate_seq_quality(quality)
        if mean_quality <= quality_threshold:
            continue
        filtered_seqs[name] = (sequence, quality)

    return filtered_seqs
