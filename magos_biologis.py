import os
import sys
sys.path.append('canticles')
from canticles import dna_rna_canticles
from canticles import protein_canticles
from canticles import fastq_canticles


def run_dna_rna_canticles(*args) -> str | int | float | list[str] | list[int] | list[float]:
    """
    Runs DNA/RNA canticles based on the specified arguments.

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
        result = dna_rna_canticles.process_single_sequence(sequence, action=action)
    else:
        result = dna_rna_canticles.process_multiple_sequences(sequences, action=action)

    return result


def run_protein_canticles(command,
                               inp,
                               *args,
                               **kwargs):
    """
    Accepts command and runs it on input data with params

    Args:
    - command (str): Valid command from command_dct
    - inp (str): Input in form of path, seq, seq list or seq dct

    Returns:
    - output_dct (dict): dict where keys are number or name of seq and values are results of command run
    """
    output_dct = {}
    input_dct = protein_canticles.parse_input(inp, **kwargs)

    command_dct = {
        'find_sites': protein_canticles.find_sites,
        'get_protein_rnas': protein_canticles.get_protein_rnas,
        'get_protein_rnas_number': protein_canticles.get_protein_rnas_number,
        'get_frameshift_proteins': protein_canticles.get_frameshift_proteins,
        'is_protein_valid': protein_canticles.is_protein_valid,
        'get_length_of_protein': protein_canticles.get_length_of_protein,
        'count_aa': protein_canticles.count_aa,
        'get_fracture_of_aa': protein_canticles.get_fracture_of_aa,
        'calculate_protein_mass': protein_canticles.calculate_protein_mass,
        'get_atomic_mass': protein_canticles.get_atomic_mass,
        'convert_aa_name': protein_canticles.convert_aa_name,
        }

    for name in input_dct:
        if command in command_dct and command != 'get_atomic_mass':
            if protein_canticles.is_protein_valid(input_dct[name]):
                output_dct[name] = command_dct[command](input_dct[name], *args, **kwargs)
            else:
                raise ValueError('Invalid protein sequence')
        elif command == 'get_atomic_mass':
            output_dct[name] = command_dct[command](input_dct[name], *args, **kwargs)
        else:
            raise ValueError('Invalid command')
    if len(output_dct) == 1:
        return output_dct[list(output_dct.keys())[0]]

    return output_dct


def filter_fastq(path_to_seqs: str, output_file_name: str = None, gc_bounds: tuple | int = (0, 100),
                length_bounds: tuple | int = (0, 2**32), quality_threshold: int = 0) -> None:
    """
    Filters fastq sequences based on the GC-content, length and quality parameters.

    Args:
    - path_to_seqs (str): The path to the Fastq file to be filtered.
    - output_file_name (str): The name of the file where the filtered Fastq sequences will be saved.
    - gc_bounds (tuple, int): GC content range (in percentages) for filtering. If you pass a single number 
    to the argument, it is assumed to be an upper bound.
    - length_bounds (tuple, int): Length range for filtering. If you pass a single number to the argument,
    it is assumed to be an upper bound.
    - quality_threshold (int): Threshold value for average read quality filtering

    Returns:
    - None
    """
    if output_file_name is None:
        output_file_name = os.path.basename(path_to_seqs)

    input_seqs = fastq_canticles.parse_fastq(path_to_seqs)
    filtered_seqs = {}

    for name, (sequence, quality) in input_seqs.items():
        gc_content = fastq_canticles.calculate_gc_content(sequence)
        if isinstance(gc_bounds, tuple):
            if not (gc_bounds[0] <= gc_content <= gc_bounds[1]):
                continue
        elif gc_content >= gc_bounds:
            continue
        sequence_length = fastq_canticles.calculate_seq_length(sequence)
        if isinstance(length_bounds, tuple):
            if not (length_bounds[0] <= sequence_length <= length_bounds[1]):
                continue
        elif sequence_length >= length_bounds:
            continue
        mean_quality = fastq_canticles.calculate_seq_quality(quality)
        if mean_quality <= quality_threshold:
            continue
        filtered_seqs[name] = (sequence, quality)

    fastq_canticles.write_filtered_fastq(filtered_seqs=filtered_seqs, output_file_name=output_file_name)

    return None
