DNA_NUCLEOTIDES = set('ATGCatgc')
RNA_NUCLEOTIDES = set('AUGCaugc')
DNA_DICT = {
    'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A',
    'a': 't', 'c': 'g', 'g': 'c', 't': 'a'
    }
RNA_DICT = {
    'A': 'U', 'C': 'G', 'G': 'C', 'U': 'A',
    'a': 'u', 'c': 'g', 'g': 'c', 'u': 'a'
    }


def is_dna(sequense):
    unique_chars = set(sequense)
    return unique_chars <= DNA_NUCLEOTIDES


def is_rna(sequense):
    unique_chars = set(sequense)
    return unique_chars <= RNA_NUCLEOTIDES


def transcribe(sequence):
    if is_dna(sequence):
        transcript = sequence.replace('T', 'U').replace('t', 'u')
    elif is_rna(sequence):
        transcript = sequence.replace('U', 'T').replace('u', 't')

    return transcript


def get_reverse(sequence):
    reverse = sequence[::-1]

    return reverse


def get_complement(sequence):
    if is_dna(sequence):
        complement = ''.join([DNA_DICT[base] for base in sequence])
    elif is_rna(sequence):
        complement = ''.join([RNA_DICT[base] for base in sequence])

    return complement


def get_reverse_complement(sequence):
    if is_dna(sequence):
        reverse_complement = ''.join([DNA_DICT[base] for base in sequence[::-1]])
    elif is_rna(sequence):
        reverse_complement = ''.join([RNA_DICT[base] for base in sequence[::-1]])

    return reverse_complement


def get_sequence_lenght(sequence):
    return str(len(sequence))


def calculate_gc_content(sequence):
    gc_count = 0

    for nucleotide in sequence:
        if nucleotide in 'GCgc':
            gc_count += 1

    return str(gc_count / len(sequence))


actions = {
    'transcribe': transcribe,
    'reverse': get_reverse,
    'complement': get_complement,
    'reverse_complement': get_reverse_complement,
    'lenght': get_sequence_lenght,
    'gc_content': calculate_gc_content,
    }


def run_mono(sequence, action):
    if is_dna(sequence) or is_rna(sequence):
        if action not in actions:
            raise ValueError('Недопустимая операция')
        result = actions[action](sequence)

    else:
        raise ValueError('Недопустимая последовательность')

    return result


def run_multiple(sequences, action):

    results = []

    for sequence in sequences:
        if is_dna(sequence) or is_rna(sequence):
            if action not in actions:
                raise ValueError('Недопустимая операция')

            result = actions[action](sequence)
            results.append(result)

        else:
            raise ValueError('Недопустимая последовательность')

    return results


def run_dna_rna_tools(*args):
    if len(args) < 2:
        raise ValueError('Недостаточно аргументов для выполнения операции')

    sequences = args[:-1]
    action = args[-1]

    if len(args) == 2:
        sequence = sequences[0]
        result = run_mono(sequence, action=action)
    else:
        result = run_multiple(sequences, action=action)

    return result
