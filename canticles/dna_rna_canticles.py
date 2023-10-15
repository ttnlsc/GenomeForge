DNA_NUCLEOTIDES = set('ATGCatgc')
RNA_NUCLEOTIDES = set('AUGCaugc')
DNA_COMPLEMENT_MAP: dict[str, str] = {
    'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A',
    'a': 't', 'c': 'g', 'g': 'c', 't': 'a'
    }
RNA_COMPLEMENT_MAP: dict[str, str] = {
    'A': 'U', 'C': 'G', 'G': 'C', 'U': 'A',
    'a': 'u', 'c': 'g', 'g': 'c', 'u': 'a'
    }


def is_dna(sequence: str) -> bool:
    """
    Checks if the input sequence is DNA.

    The function compares the unique characters of the input sequence with the set DNA_NUCLEOTIDES.

    Args:
    - sequence(str): The sequence to be checked.

    Returns:
    - bool: The result of the check.
    """
    unique_chars = set(sequence)
    return unique_chars <= DNA_NUCLEOTIDES


def is_rna(sequence: str) -> bool:
    """
    Checks if the input sequence is RNA.

    The function compares the unique characters of the input sequence with the set RNA_NUCLEOTIDES.

    Args:
    - sequence(str): The sequence to be checked.

    Returns:
    - bool: The result of the check.
    """
    unique_chars = set(sequence)
    return unique_chars <= RNA_NUCLEOTIDES


def transcribe(sequence: str) -> str:
    """
    Transcribes the input nucleotide sequence.

    The function transcribes the input sequence from DNA to RNA or vice versa, depending on the type of sequence.

    Args:
    - sequence (str): The nucleotide sequence to be transcribed.

    Returns:
    - str: The transcribed nucleotide sequence.
    """
    transcript = ''

    if is_dna(sequence):
        transcript = sequence.replace('T', 'U').replace('t', 'u')
    elif is_rna(sequence):
        transcript = sequence.replace('U', 'T').replace('u', 't')

    return transcript


def get_reverse(sequence: str) -> str:
    """
    Returns the reverse of the input sequence.

    Args:
    - sequence (str): The input sequence to be reversed.

    Returns:
    - str: The reversed sequence.
    """
    return sequence[::-1]


def get_complement(sequence: str) -> str:
    """
    Returns the complement of the input sequence.

    Args:
    - sequence (str): The input sequence to be complemented.

    Returns:
    - str: The complemented sequence.
    """
    complement = ''
    
    if is_dna(sequence):
        complement = ''.join([DNA_COMPLEMENT_MAP[base] for base in sequence])
    elif is_rna(sequence):
        complement = ''.join([RNA_COMPLEMENT_MAP[base] for base in sequence])

    return complement


def get_reverse_complement(sequence: str) -> str:
    """
    Returns the reverse complement of the input sequence.

    Args:
    - sequence (str): The input sequence to be reversed and complemented.

    Returns:
    - str: The reversed and complemented sequence.
    """
    return get_complement(get_reverse(sequence))


def get_sequence_lenght(sequence: str) -> int:
    """
    Returns the length of the input sequence.

    Args:
    - sequence (str): The input sequence.

    Returns:
    - int: The length of the sequence.
    """
    return len(sequence)


def calculate_gc_content(sequence: str) -> float:
    """
    Calculates the GC content (Guanine-Cytosine content) of the input sequence.

    Args:
    - sequence (str): The input sequence.

    Returns:
    - float: The calculated GC content as a float value.
    """
    gc_count = 0

    for nucleotide in sequence:
        if nucleotide in 'GCgc':
            gc_count += 1

    return gc_count / len(sequence)


ACTION_MAP = {
    'transcribe': transcribe,
    'reverse': get_reverse,
    'complement': get_complement,
    'reverse_complement': get_reverse_complement,
    'lenght': get_sequence_lenght,
    'gc_content': calculate_gc_content,
    }


def process_single_sequence(sequence: str, action: str) -> str | int | float:
    """
    Processes a single nucleotide sequence based on the specified action.

    Args:
    - sequence (str): The input nucleotide sequence to be processed.
    - action (str): The action to perform on the sequence.

    Returns:
    - str | int | float: The result of the specified action on the input sequence.
    """
    if is_dna(sequence) or is_rna(sequence):
        if action not in ACTION_MAP:
            raise ValueError(f'Invalid operation: {action} is not a valid operation.')
        result = ACTION_MAP[action](sequence)

    else:
        raise ValueError(f'Invalid sequence: {sequence} is not a valid DNA or RNA sequence.')

    return result


def process_multiple_sequences(sequences: list[str], action: str) -> list[str] | list[int] | list[float]:
    """
    Processes a list of nucleotide sequences based on the specified action.

    Args:
    - sequences (List[str]): A list of input nucleotide sequences to be processed.
    - action (str): The action to perform on each sequence.

    Returns:
    - list[str] | list[int] | list[float]: A list of results, where each result corresponds to a processed
    input sequence.
    """
    results = []

    for sequence in sequences:
        if is_dna(sequence) or is_rna(sequence):
            if action not in ACTION_MAP:
                raise ValueError(f'Invalid operation: {action} is not a valid operation.')

            result = ACTION_MAP[action](sequence)
            results.append(result)

        else:
            raise ValueError(f'Invalid sequence: {sequence} is not a valid DNA or RNA sequence.')

    return results


def run_dna_rna_tools(*args) -> str | int | float | list[str] | list[int] | list[float]:
    """
    Runs DNA/RNA tools based on the specified arguments.

    Args:
    - args: Variable-length argument list. Should contain one or more sequences followed by an action.

    Returns:
    - str | int | float | list[str] | list[int] | list[float]: The result of the specified action(s).
    If a single sequence is provided, a string, integer or float is returned.
    If multiple sequences are provided, a list of strings, integers or floats is returned.
    """
    if len(args) < 2:
        raise ValueError('Insufficient arguments to execute operation')

    sequences = args[:-1]
    action = args[-1]

    if len(args) == 2:
        sequence = sequences[0]
        result = process_single_sequence(sequence, action=action)
    else:
        result = process_multiple_sequences(sequences, action=action)

    return result
