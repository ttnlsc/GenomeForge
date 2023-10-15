H2O_WEIGHT: float = 18.01468

AA_MASS_DICT = {
    'G': 75.0659, 'g': 75.0659, 'L': 131.17262, 'l': 131.17262,
    'Y': 181.18894, 'y': 181.18894, 'S': 105.09158, 's': 105.09158,
    'E': 147.12826, 'e': 147.12826, 'Q': 146.1438, 'q': 146.1438,
    'D': 133.10158, 'd': 133.10158, 'N': 132.11712, 'n': 132.11712,
    'F': 165.18994, 'f': 165.18994, 'A': 89.09258, 'a': 89.09258,
    'K': 146.18716, 'k': 146.18716, 'R': 174.20056, 'r': 174.20056,
    'H': 155.15466, 'h': 155.15466, 'C': 121.15758, 'c': 121.15758,
    'V': 117.14594, 'v': 117.14594, 'P': 115.13026, 'p': 115.13026,
    'W': 204.22648, 'w': 204.22648, 'I': 131.17262, 'i': 131.17262,
    'M': 149.21094, 'm': 149.21094, 'T': 119.11826, 't': 119.11826,
}

ATOMIC_MASS = {
    'C': 12.011, 'H': 1.00784, 'O': 15.999, 'N': 14.0067, 'S': 32.065, 'Fe': 55.845, 'Se': 78.96,
}

AA_NAME_DICT = {
    'G': 'Gly', 'g': 'Gly',
    'L': 'Leu', 'l': 'Leu',
    'Y': 'Tyr', 'y': 'Tyr',
    'S': 'Ser', 's': 'Ser',
    'E': 'Glu', 'e': 'Glu',
    'Q': 'Gln', 'q': 'Gln',
    'D': 'Asp', 'd': 'Asp',
    'N': 'Asn', 'n': 'Asn',
    'F': 'Phe', 'f': 'Phe',
    'A': 'Ala', 'a': 'Ala',
    'K': 'Lys', 'k': 'Lys',
    'R': 'Arg', 'r': 'Arg',
    'H': 'His', 'h': 'His',
    'C': 'Cys', 'c': 'Cys',
    'V': 'Val', 'v': 'Val',
    'P': 'Pro', 'p': 'Pro',
    'W': 'Trp', 'w': 'Trp',
    'I': 'Ile', 'i': 'Ile',
    'M': 'Met', 'm': 'Met',
    'T': 'Thr', 't': 'Thr'
}

AA_RNA_TABLE = {
    'F': ['UUU', 'UUC'],
    'L': ['UUA', 'UUG', 'CUU', 'CUC', 'CUA', 'CUG'],
    'S': ['UCU', 'UCC', 'UCA', 'UCG', 'AGU', 'AGC'],
    'Y': ['UAU', 'UAC'],
    '*': ['UAA', 'UAG', 'UGA', 'uaa', 'uag', 'uga'],
    'C': ['UGU', 'UGC'],
    'W': ['UGG'],
    'P': ['CCU', 'CCC', 'CCA', 'CCG'],
    'H': ['CAU', 'CAC'],
    'Q': ['CAA', 'CAG'],
    'R': ['CGU', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],
    'I': ['AUU', 'AUC', 'AUA'],
    'M': ['AUG'],
    'T': ['ACU', 'ACC', 'ACA', 'ACG'],
    'N': ['AAU', 'AAC'],
    'K': ['AAA', 'AAG'],
    'V': ['GUU', 'GUC', 'GUA', 'GUG'],
    'A': ['GCU', 'GCC', 'GCA', 'GCG'],
    'D': ['GAU', 'GAC'],
    'E': ['GAA', 'GAG'],
    'G': ['GGU', 'GGC', 'GGA', 'GGG'],
    'f': ['uuu', 'uuc'],
    'l': ['uua', 'uug', 'cuu', 'cuc', 'cua', 'cug'],
    's': ['ucu', 'ucc', 'uca', 'ucg', 'agu', 'agc'],
    'y': ['uau', 'uac'],
    'c': ['ugu', 'ugc'],
    'w': ['ugg'],
    'p': ['ccu', 'ccc', 'cca', 'ccg'],
    'h': ['cau', 'cac'],
    'q': ['caa', 'cag'],
    'r': ['cgu', 'cgc', 'cga', 'cgg', 'aga', 'agg'],
    'i': ['auu', 'auc', 'aua'],
    'm': ['aug'],
    't': ['acu', 'acc', 'aca', 'acg'],
    'n': ['aau', 'aac'],
    'k': ['aaa', 'aag'],
    'v': ['guu', 'guc', 'gua', 'gug'],
    'a': ['gcu', 'gcc', 'gca', 'gcg'],
    'd': ['gau', 'gac'],
    'e': ['gaa', 'gag'],
    'g': ['ggu', 'ggc', 'gga', 'ggg']
}

RNA_CODON_TABLE = {
    'UUU': 'F', 'UUC': 'F', 'UUA': 'L', 'UUG': 'L', 'UCU': 'S', 'UCC': 'S', 'UCA': 'S', 'UCG': 'S',
    'UAU': 'Y', 'UAC': 'Y', 'UAA': '*', 'UAG': '*', 'UGU': 'C', 'UGC': 'C', 'UGA': '*', 'UGG': 'W',
    'CUU': 'L', 'CUC': 'L', 'CUA': 'L', 'CUG': 'L', 'CCU': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'CAU': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q', 'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'AUU': 'I', 'AUC': 'I', 'AUA': 'I', 'AUG': 'M', 'ACU': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'AAU': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K', 'AGU': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GUU': 'V', 'GUC': 'V', 'GUA': 'V', 'GUG': 'V', 'GCU': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'GAU': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E', 'GGU': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
    'uuu': 'f', 'uuc': 'f', 'uua': 'l', 'uug': 'l', 'ucu': 's', 'ucc': 's', 'uca': 's', 'ucg': 's',
    'uau': 'y', 'uac': 'y', 'uaa': '*', 'uag': '*', 'ugu': 'c', 'ugc': 'c', 'uga': '*', 'ugg': 'w',
    'cuu': 'l', 'cuc': 'l', 'cua': 'l', 'cug': 'l', 'ccu': 'p', 'ccc': 'p', 'cca': 'p', 'ccg': 'p',
    'cau': 'h', 'cac': 'h', 'caa': 'q', 'cag': 'q', 'cgu': 'r', 'cgc': 'r', 'cga': 'r', 'cgg': 'r',
    'auu': 'i', 'auc': 'i', 'aua': 'i', 'aug': 'm', 'acu': 't', 'acc': 't', 'aca': 't', 'acg': 't',
    'aau': 'n', 'aac': 'n', 'aaa': 'k', 'aag': 'k', 'agu': 's', 'agc': 's', 'aga': 'r', 'agg': 'r',
    'guu': 'v', 'guc': 'v', 'gua': 'v', 'gug': 'v', 'gcu': 'a', 'gcc': 'a', 'gca': 'a', 'gcg': 'a',
    'gau': 'd', 'gac': 'd', 'gaa': 'e', 'gag': 'e', 'ggu': 'g', 'ggc': 'g', 'gga': 'g', 'ggg': 'g'
}


def read_seq_from_fasta(path_to_seq: str,
                        use_full_name: bool = False,
                        **_) -> dict:
    """
    Reads sequences from a FASTA file and returns a dictionary.

    Args:
    - path_to_seq (str): The path to the FASTA file.
    - use_full_name (bool): Set to True if you want to use the full sequence name with properties.

    Returns:
    - dict: A dictionary with sequence names as keys and sequences as values, e.g. {'seq_name': 'sequence'}
    """
    sequences = {}
    current_name = None

    with open(path_to_seq) as file:
        for line in file:
            line = line.strip()
            if line.startswith('>'):
                if current_name is not None:
                    sequences[current_name] = ''.join(current_sequence)
                if use_full_name:
                    current_name = line[1:]
                else:
                    current_name = line[1:].split()[0]
                current_sequence = []
            else:
                current_sequence.append(line)

        sequences[current_name] = ''.join(current_sequence)

    return sequences


def is_protein_valid(seq: str) -> bool:
    """
    Checks if protein is valid.

    Args:
    - seq (str): seq to be checked

    Returns:
    - bool, the result of the check
    """

    if set(seq).issubset(RNA_AA_TABLE):
        return True
    return False


def find_sites(seq: str,
               *sites: str,
               is_one_based: bool = False,
               **_) -> dict:
    """
    Finds indexes of given sites.

    Args:
    - seq (str): seq to be checked
    - *args (str): sites to be found
    - is_one_based (bool): whether result should be 0- (False) or 1-indexed (True). Default False

    Returns:
    - dict: dictionary of sites as keys and lists of indexes for the site where it's been found
    """
    found_sites = {}

    for site in sites:
        site_length = len(site)
        for idx in range(len(seq) - site_length + 1):
            scatter = seq[idx:idx + site_length]
            if scatter == site:
                if is_one_based:
                    index = idx + 1
                else:
                    index = idx
                found_sites.setdefault(site, []).append(index)

    return found_sites


def get_protein_rnas(seq: str,
                     check_if_user_conscious: bool = False,
                     **_) -> list or None:
    """
    Returns list of all possible RNA's from which can serve as matrix for protein synthesis.

    WARNING: can be computationally intensive on longer sequences,
    will NOT start unless check_if_user_conscious is True!

    Args:
    - seq (str): seq to be checked
    - check_if_user_conscious (bool): checks user's consciousness. Default False

    Returns:
    - list: list of possible RNA's as str or None if not conscious (set check_if_user_conscious to True)
    """

    if check_if_user_conscious:
        kmers = ['']  # set initial kmers
        for amino_acid in seq:  # iterate AAs
            current_kmers = []
            codons = AA_RNA_TABLE[amino_acid]  # get list of codons for AA
            if amino_acid not in AA_RNA_TABLE:
                raise ValueError(f"Amino acid '{amino_acid}' is not found in the AA_RNA_TABLE.")
            for codon in codons:
                for kmer in kmers:
                    current_kmers.append(''.join([kmer + codon]))  # append every codon to existing kmers
            kmers = current_kmers  # re-write k-mers for next iteration

        return kmers
    else:
        print("You don't know what you're doing!")  # politely ask user to reconsider their actions
        return None


def get_protein_rnas_number(seq: str, **_) -> int:
    """
    Get number of all possible RNA's for a given protein.

    Args:
    - seq (str): seq to be checked

    Returns:
    - int: number of possible RNA's for seq
    """
    rnas_num = 1

    for amino_acid in seq:
        rnas_num *= len(RNA_AA_TABLE[amino_acid])
    return rnas_num


def check_all_upper(codon: str) -> bool:
    """
    Checks whether all letters in colon are upper

    Args:
    - codon (str): codon to be checked

    Returns:
    - bool: if all letters are uppercase
    """
    check_upper = True

    for letter in (set(codon)):
        check_upper = letter.isupper() and check_upper
    return check_upper


def get_frameshift_proteins(seq: str,
                            check_if_user_conscious: bool = False,
                            is_stop_codon_termination_enabled: bool = False,
                            **_) -> dict or None:
    """
    Returns list of all possible proteins from all possible frames in peptide.

    WARNING: can be computationally intensive on longer sequences,
    will NOT start unless check_if_user_conscious is True!

    Args:
    - seq (str): seq to be checked
    - check_if_user_conscious (bool): checks user's consciousness. Default False
    - is_stop_codon_termination_enabled (bool): terminate translation when reached stop-codon. Default False.

    Returns:
    - dict: dict of lists of all possible frames proteins:
    {frame_0: ['protein_seqs'], frame_1: ['protein_seqs'], frame_2: ['protein_seqs']}
    """

    if check_if_user_conscious:
        frameshift_dct = {'frame_0': [seq]}  # set current seq as frame_0 (protein from not-shifted frame)
        rnas = get_protein_rnas(seq, check_if_user_conscious=check_if_user_conscious)
        for frame_number in [1, 2]:
            frames_list = []
            for rna in rnas:
                frame = ''
                for i in range(frame_number, len(rna) - (frame_number + 1), 3):  # set frame-dependent range to iterate
                    frame_codon = rna[i:i + 3]  # extract codon
                    if not check_all_upper(frame_codon):  # check if all letters in codon uppercase
                        frame_codon = frame_codon.lower()  # if not change all to lowercase
                    frame += RNA_CODON_TABLE[frame_codon]
                    if (is_stop_codon_termination_enabled and
                            RNA_CODON_TABLE[frame_codon] == '*'):  # stop writing if meet stop-codon
                        break
                frames_list.append(frame)  # append frame to frames list
            frameshift_dct[f'frame_{frame_number}'] = list(set(frames_list))  # clean duplicates and write to dict
        return frameshift_dct
    else:
        print("You don't know what you're doing!")  # politely ask user to reconsider their actions
        return None


def get_length_of_protein(seq: str, **_) -> int:
    """
    Calculates the length of a protein.

    Argument:
    - seq (str): sequence to calculate the length

    Return:
    - int: sequence length
    """
    return len(seq)


def count_aa(seq: str, aminoacids: str = None, **_) -> dict:
    """
    Counts the number of given or all amino acids in a protein sequence.

    Arguments:
    - seq (str): sequence to count amino acids
    - aminoacids (str): which amino acids to count in sequence

    Return:
    - dict: a dictionary with amino acids and its count
    """

    aa_dict_count = {}
    if (aminoacids is None) or (aminoacids == ''):
        '''
        I added an additional condition for user-friendly experience.
        E.g., we can want to find specific aminoacid, look on result and then look on all aminoacids.
        Without this condition we have to delete keyword argument, but with it we can only make it empty.
        '''
        aminoacids = ''.join(set(seq))
    for aa in aminoacids:
        aa_dict_count[aa] = seq.count(aa)
    return aa_dict_count


def get_fracture_of_aa(seq: str, show_as_percentage: bool = False, aminoacids: str = None, **_) -> dict:
    """
    Calculates the fracture or percentage of amino acids in a protein sequence.

    Args:
    - seq (str): sequence in which you need to calculate the fracture of amino acids
    - show_as_percentage (bool): change it to True, if you want to get results with percentages
    - aminoacids (str): the fracture of which amino acids to count in the sequence

    Returns:
    - dict: a dictionary with amino acids and its fracture or percentage
    """

    if show_as_percentage:
        mult = 100
        round_var = 2
    else:
        mult = 1
        round_var = 4
    aa_dict_count = count_aa(seq, aminoacids=aminoacids)
    aa_dict_percent = {}
    len_of_protein = get_length_of_protein(seq)
    for aa, count in aa_dict_count.items():
        aa_dict_percent[aa] = round(count / len_of_protein * mult, round_var)
    return aa_dict_percent


def calculate_protein_mass(sequence: str, **_) -> float:
    """
    Calculates the molecular mass of a protein based on its amino acid sequence and a dictionary of amino acid masses.

    Arguments:
    - sequence(str or list): A string or list of characters representing the amino acid sequence.

    Return:
    - float: The molecular mass of a protein in atomic mass units, rounded to the third decimal place.
    """
    total_mass = 0.0

    for aa in sequence:
        if aa in AA_MASS_DICT:
            total_mass += AA_MASS_DICT[aa]
        else:
            raise ValueError(f'Unknown amino acid: {aa}')
    total_mass = total_mass - H2O_WEIGHT * (len(sequence) - 1)

    return round(total_mass, 3)


def get_atomic_mass(chem: str, **_) -> float:
    """
    Calculates the molecular mass of a biological molecule, primarily an amino acid, based on a simple chemical formula.
    A dictionary linking the chemical elements Carbon, Hydrogen, Oxygen, Nitrogen, Sulfur, Iron and Selenium
    with their masses in atomic mass units.

    Args:
    - chem (str): String representing a simple chemical formula, e.g. C2H5OH

    Returns:
    - float: Molecular mass of a biological molecule in atomic mass units.
    """
    total_mass = 0
    element = []
    number = []

    for character in chem:
        if character.isalpha():
            if character.isupper():
                if element and number:
                    if ''.join(element) not in ATOMIC_MASS:
                        raise ValueError(f'Unexpected character: {character}')
                    total_mass += ATOMIC_MASS[''.join(element)] * int(''.join(number))
                    number.clear(), element.clear()
                    element.append(character)
                elif element:
                    if ''.join(element) not in ATOMIC_MASS:
                        raise ValueError(f'Unexpected character: {character}')
                    total_mass += ATOMIC_MASS[''.join(element)]
                    element.clear()
                    element.append(character)
                else:
                    element.append(character)
            else:
                element.append(character)
        elif character.isdigit():
            number.append(character)
        else:
            raise ValueError(f'Unexpected character: {character}')

    if ''.join(element) not in ATOMIC_MASS:
        raise ValueError(f'Unexpected character: {character}')
    else:
        if element and number:
            total_mass += ATOMIC_MASS[''.join(element)] * int(''.join(number))
        elif element:
            total_mass += ATOMIC_MASS[''.join(element)]

    return round(total_mass, 3)


def convert_aa_name(sequence: str, name_dict: dict = None, sep: str = '',
                    use_default_register: bool = True) -> str:
    """
    Converts a sequence of one-letter amino acid codes to three-letter designations.

    Arguments:
    - sequence (str): String with one-letter amino acid codes.
    - name_dict (dict[str, str], optional): A dictionary linking one-letter codes to three-letter designations.
    If not provided, the standard AA_NAME_DICT dictionary is used.
    - sep (str, optional): Separator between three-letter amino acid designations. There is no delimiter by default.
    - use_default_register(bool, optional): Determines whether to preserve letter case in three-letter designations.
    If True, the letters will be converted to upper or lower case depending on the case of the depending
    on the case of the one-letter code. The default is False.

    Return:
    - str: A string of three-letter amino acid designations separated by the specified delimiter.
    """

    new_name = ''
    if name_dict is None:
        name_dict = AA_NAME_DICT
    for i, aa in enumerate(sequence):
        if aa in name_dict:
            if use_default_register is False:
                new_name += name_dict[aa]
            elif use_default_register is True:
                if aa.isupper():
                    new_name += name_dict[aa].upper()
                else:
                    new_name += name_dict[aa].lower()
            else:
                if aa.isupper():
                    new_name += name_dict[aa].lower()
                else:
                    new_name += name_dict[aa].upper()
            if sep and (i + 1) < len(sequence):
                new_name += sep
        else:
            raise ValueError(f'Unknown amino acid: {aa}')
    return new_name

# defined later to let all funcs be initialized before passed here
command_dct = {
    'find_sites': find_sites,
    'get_protein_rnas': get_protein_rnas,
    'get_protein_rnas_number': get_protein_rnas_number,
    'get_frameshift_proteins': get_frameshift_proteins,
    'is_protein_valid': is_protein_valid,
    'get_length_of_protein': get_length_of_protein,
    'count_aa': count_aa,
    'get_fracture_of_aa': get_fracture_of_aa,
    'calculate_protein_mass': calculate_protein_mass,
    'get_atomic_mass': get_atomic_mass,
    'convert_aa_name': convert_aa_name,
    }



def parse_input(inp: str, **kwargs) -> dict:
    """
    Parses input and returns dict of seqs.

    Arguments:
    - inp (str): Input path or seq or dict of seqs or list of seqs
    - **kwargs: Additional keyword arguments to be passed to input reader (e.g. )

    Return:
    - parsed_dct (dict): dict where keys are number or name of seq and value of seq
    """

    parsed_dct = {}
    inp_type = type(inp)  # get input type
    if inp_type == list:
        for i, seq in enumerate(inp):
            parsed_dct |= {i: seq}
    elif inp_type == dict:
        parsed_dct = inp
    elif inp_type == str and '.' in inp:  # check whether input has file extension symbols
        parsed_dct = read_seq_from_fasta(inp, **kwargs)
    elif inp_type == str:
        parsed_dct = {0: inp}

    return parsed_dct


def run_ultimate_protein_tools(command,
                               inp,
                               *args,
                               **kwargs):
    """
    Accepts command and runs it on input data with params

    Arguments:
    - command (str): Valid command from command_dct
    - inp (str): Input in form of path, seq, seq list or seq dct

    Return:
    - output_dct (dict): dict where keys are number or name of seq and values are results of command run
    """
    output_dct = {}
    input_dct = parse_input(inp, **kwargs)
    for name in input_dct:
        if command in command_dct and command != 'get_atomic_mass':
            if is_protein_valid(input_dct[name]):
                output_dct[name] = command_dct[command](input_dct[name], *args, **kwargs)
            else:
                output_dct[name] = is_protein_valid(input_dct[name])
        elif command == 'get_atomic_mass':
            output_dct[name] = command_dct[command](input_dct[name], *args, **kwargs)
        else:
            print('Command invalid')
    if len(output_dct) == 1:
        return output_dct[list(output_dct.keys())[0]]
    return output_dct
