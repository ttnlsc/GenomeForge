from __future__ import annotations
import os
import sys
import re
import datetime
import requests
from abc import ABC, abstractmethod
from Bio import SeqIO, SeqUtils
from dataclasses import dataclass
from io import BytesIO, StringIO
from bs4 import BeautifulSoup
from typing import Callable, Optional, List
from dotenv import load_dotenv

load_dotenv()


class BiologicalSequence(ABC):
    """
    Abstract base class representing a biological sequence.

    Defines methods that must be implemented by concrete subclasses:
    - __len__: Returns the length of the sequence.
    - __getitem__: Returns the item at the specified index.
    - __str__: Returns the string representation of the sequence.
    - is_valid_alphabet: Checks if the sequence contains valid symbols according to its alphabet.

    Subclasses must implement these methods to provide functionality specific to the type of biological sequence.
    """

    @abstractmethod
    def __len__(self):
        pass

    @abstractmethod
    def __getitem__(self, item):
        pass

    @abstractmethod
    def __str__(self):
        pass

    @abstractmethod
    def is_valid_alphabet(self):
        pass


class NucleicAcidSequence(BiologicalSequence):
    """
    Represents a nucleic acid sequence.
    Methods:
    - complement(): Returns the complemented sequence.
    - gc_content(percentage: bool = True) -> int | float: Calculates the GC content of the sequence.
    Raises:
    - NotImplementedError: If complement method is not implemented for this class.
    """
    alphabet = None
    complement_map = None

    def __init__(self, sequence: str):
        self.sequence = sequence

    def __len__(self) -> int:
        return len(self.sequence)

    def __getitem__(self, item: int) -> str:
        return self.sequence[item]

    def __str__(self) -> str:
        return f"{self.sequence}"

    def is_valid_alphabet(self) -> bool:
        """
        Checks if the sequence uses a valid alphabet.
        """
        if self.alphabet is None:
            raise NotImplementedError("Is valid alphabet method not implemented for this class.")
        return set(self.sequence).issubset(self.alphabet)

    def complement(self) -> NucleicAcidSequence:
        """
        Returns the complement of the sequence based on the provided complement map.
        Args:
        complement_map (dict or None): A dictionary specifying the mapping of nucleotide bases to their complements.
            If None, raises NotImplementedError.
        Returns:
        NucleicAcidSequence: A new instance of the same class representing the complemented sequence.
        Raises:
        NotImplementedError: If complement_map is None, indicating that the complement method is not implemented for
            this class.
        """
        if self.complement_map is None:
            raise NotImplementedError("Complement method not implemented for this class.")
        complemented_sequence = "".join([self.complement_map[base] for base in self.sequence])

        return self.__class__(complemented_sequence)

    def gc_content(self, percentage: bool = True) -> int | float:
        """
        Calculates the GC content of the sequence.
        Args:
        - percentage (bool): If True, returns the result as a percentage.
        Returns:
        - int | float: GC content value.
        """
        gc_symbols = set("GCgc")
        gc_count = sum(1 for nucleotide in self.sequence if nucleotide in gc_symbols)
        if percentage:
            return (gc_count / self.sequence.__len__()) * 100
        else:
            return gc_count / self.sequence.__len__()


class DNASequence(NucleicAcidSequence):
    """
    Represents a DNA sequence.
    Methods:
    - transcribe(): Returns the transcribed RNA sequence.
    """

    alphabet = set("ATGCatgc")
    complement_map = {
        "A": "T",
        "C": "G",
        "G": "C",
        "T": "A",
        "a": "t",
        "c": "g",
        "g": "c",
        "t": "a",
    }

    def __init__(self, sequence: str):
        super().__init__(sequence)

    def transcribe(self) -> RNASequence:
        """
        Returns the transcribed RNA sequence.
        Returns:
        - RNASequence: Transcribed RNA sequence.
        """
        transcribed_sequence = self.sequence.replace("T", "U").replace("t", "u")
        return RNASequence(transcribed_sequence)


class RNASequence(NucleicAcidSequence):
    alphabet = set("AUGCaugc")
    complement_map = {
        "A": "U",
        "C": "G",
        "G": "C",
        "U": "A",
        "a": "u",
        "c": "g",
        "g": "c",
        "u": "a",
    }

    def __init__(self, sequence: str):
        super().__init__(sequence)


class AminoAcidSequence(BiologicalSequence):
    """
    Represents an amino acid sequence.
    """
    alphabet = set("GgLlYySsEeQqDdNnFfAaKkRrHhCcVvPpWwIiMmTt")

    def __init__(self, sequence: str):
        self.sequence = sequence

    def __len__(self) -> int:
        return len(self.sequence)

    def __getitem__(self, item: int) -> str:
        return self.sequence[item]

    def __str__(self) -> str:
        return f"> {self.sequence}"

    def is_valid_alphabet(self) -> bool:
        """
        Checks if the sequence uses a valid alphabet.
        """
        return set(self.sequence).issubset(self.alphabet)

    def count_aa(self) -> dict:
        """
        Counts the number of given or all amino acids in a protein sequence.
        Arguments:
        - seq (str): sequence to count amino acids
        - aminoacids (str): which amino acids to count in sequence
        Return:
        - dict: a dictionary with amino acids and its count
        """

        aa_dict_count = {}
        for aa in set(self.sequence):
            aa_dict_count[aa] = self.sequence.count(aa)
        return aa_dict_count


def filter_fastq(
        path_to_seqs: str,
        output_file_name: str = None,
        gc_bounds: tuple | int = (0, 100),
        length_bounds: tuple | int = (0, 2 ** 32),
        quality_threshold: int = 0,
) -> None:
    """
    Filters FASTQ sequences based on the GC-content, length and quality parameters.

    Args:
    - path_to_seqs (str): the path to the FASTQ file to be filtered.
    - output_file_name (str): the name of the file where the filtered FASTQ sequences will be saved.
    - gc_bounds (tuple, int): GC content range (in percentages) for filtering. If you pass a single number
    to the argument, it is assumed to be an upper bound.
    - length_bounds (tuple, int): length range for filtering. If you pass a single number to the argument,
    it is assumed to be an upper bound.
    - quality_threshold (int): threshold value for average read quality filtering

    Returns:
    - None: the function doesn't return a value but writes the filtered FASTQ to the output file.
    """
    if output_file_name is None:
        output_file_name = os.path.basename(path_to_seqs)

    records = SeqIO.parse(path_to_seqs, "fastq")

    filtered_records = []

    for record in records:
        gc_content = SeqUtils.gc_fraction(record.seq) * 100
        if isinstance(gc_bounds, tuple):
            if not (gc_bounds[0] <= gc_content <= gc_bounds[1]):
                continue
        elif gc_content >= gc_bounds:
            continue
        sequence_length = len(record.seq)
        if isinstance(length_bounds, tuple):
            if not (length_bounds[0] <= sequence_length <= length_bounds[1]):
                continue
        elif sequence_length >= length_bounds:
            continue
        mean_quality = sum(record.letter_annotations["phred_quality"]) / sequence_length
        if mean_quality <= quality_threshold:
            continue
        filtered_records.append(record)

    with open(output_file_name, "w") as output_handle:
        SeqIO.write(filtered_records, output_handle, "fastq")


def send_telegram_message(chat_id: int,
                          text: str,
                          document_name: Optional[str] = None,
                          file: Optional[BytesIO] = None) -> dict:
    """
    Sends a message to Telegram.
    Args:
     - chat_id (int): ID of the chat room or user where the message will be sent
     - text (str): Message text
     - document_name (str): The name of the document being sent
     - file (Optional[BytesIO]): File to send (optional)
    Returns:
    - Response from Telegram API in JSON format
    """
    telegram_token = os.getenv("TG_API_TOKEN")
    url = f"https://api.telegram.org/bot{telegram_token}/sendMessage"
    params = {"chat_id": chat_id, "text": text, "parse_mode": "Markdown"}
    response = requests.post(url, params=params)
    if file:
        files = {"document": (document_name, file)}
        response = requests.post(
            f"https://api.telegram.org/bot{telegram_token}/sendDocument",
            params=params,
            files=files
        )
    return response.json()


def telegram_logger(chat_id: int) -> Callable:
    """
    Decorator function for logging and sending messages to Telegram
    Args:
    - chat_id (int): The chat ID to which messages will be sent
    Returns:
    - Callable: Decorator function that can be applied to other functions

    """

    def decorator(func):
        def wrapper(*args, **kwargs):
            start_time = datetime.datetime.now()
            stdout = sys.stdout
            stderr = sys.stderr
            sys.stdout = StringIO()
            sys.stderr = StringIO()
            try:
                result = func(*args, **kwargs)
                end_time = datetime.datetime.now()
                duration = end_time - start_time
                message = (
                    f"🥳 Function `{func.__name__}` completed successfully.\n"
                    f"Execution time: {duration}."
                )
                send_telegram_message(chat_id, message)
            except Exception as e:
                message = (
                    f"☠️ Function `{func.__name__}` failed with an exception:\n"
                    f"Type: {type(e).__name__}\n"
                    f"Error: {str(e)}\n"
                )
                send_telegram_message(chat_id, message)
                raise e
            finally:
                sys.stdout.seek(0)
                sys.stderr.seek(0)
                stdout_content = sys.stdout.read()
                stderr_content = sys.stderr.read()
                sys.stdout = stdout
                sys.stderr = stderr
                if stdout_content or stderr_content:
                    log_content = f"STDOUT:\n{stdout_content}\nSTDERR:\n{stderr_content}"
                    log_content_bytes = log_content.encode('utf-8')
                    send_telegram_message(chat_id, "", f"{func.__name__}.log", file=BytesIO(log_content_bytes))

            return result

        return wrapper

    return decorator


@dataclass
class Intron:
    gene: int
    number: int
    start: int
    end: int

    def __repr__(self):
        return f"Intron(gene='{self.gene}', number={self.number}, start={self.start}, end={self.end})"


@dataclass
class Exon:
    gene: int
    number: int
    start: int
    end: int

    def __repr__(self):
        return f"Exon(gene='{self.gene}', number={self.number}, start={self.start}, end={self.end})"


@dataclass
class GenscanOutput:
    status: str
    cds_list: List[str] | None
    intron_list: List['Intron'] | None
    exon_list: List['Exon'] | None

    def __repr__(self):
        return (
            f"GenscanOutput(status='{self.status}',\n"
            f"               cds_list={self.cds_list},\n"
            f"               intron_list={self.intron_list},\n"
            f"               exon_list={self.exon_list})"
        )


def parse_genscan_output(lines: List[str]) -> tuple[List[str], List[Exon], List[Intron]]:
    """
    Parse the output of Genscan and extract exons, introns, and protein sequences.
    Introns are considered to be regions lying between the initial, internal, and terminating exon.
    Promoters and poly-A regions are not considered exons.
    Args:
    - lines (List[str]): List of lines from the Genscan output
    Returns:
    Tuple containing:
    - List of protein sequences
    - List of Exon objects
    - List of Intron objects
    """

    pattern = r'^\d{1,2}\.\d\d'
    gene = ''
    exon_list = []
    intron_list = []
    exon_number, exon_start, exon_end = 0, 0, 0
    intron_number, intron_start, intron_end = 0, 0, 0
    proteins = {}
    protein_name = ''
    protein = []
    sequence_length = 0

    for line in lines:
        if line:
            line = line.strip()
            if line.startswith('Sequence'):
                sequence_length = int(line.split(':')[3].strip().split(' ')[0])
            if re.match(pattern, line):
                string = line.split()
                if string[1] not in ['Prom', 'PlyA']:
                    gene = int(string[0].split('.')[0])
                    exon_number = int(string[0].split('.')[1])
                    exon_start = int(string[3])
                    exon_end = int(string[4])
                    exon = Exon(gene=gene, number=exon_number, start=exon_start, end=exon_end)
                    exon_list.append(exon)
                if string[1] == 'Sngl':
                    continue
                elif string[1] == 'Init':
                    if string[2] == '+':
                        intron_number = exon_number
                        intron_start = exon_end + 1
                    elif string[2] == '-':
                        if not intron_number:
                            intron_number = exon_number - 1
                        intron_start = exon_end - 1
                        if not intron_end:
                            intron_end = 1
                        intron = Intron(gene=gene, number=intron_number, start=intron_start, end=intron_end)
                        intron_list.append(intron)
                        intron_number, intron_start, intron_end = 0, 0, 0
                elif string[1] == 'Intr':
                    if string[2] == '+':
                        if not intron_start:
                            intron_start = 1
                        intron_end = exon_start - 1
                        intron = Intron(gene=gene, number=intron_number, start=intron_start, end=intron_end)
                        intron_list.append(intron)
                        intron_number = exon_number
                        intron_start = exon_end + 1
                    elif string[2] == '-':
                        intron_number = exon_number
                        intron_start = exon_end - 1
                        if not intron_end:
                            intron_end = 1
                        intron = Intron(gene=gene, number=intron_number, start=intron_start, end=intron_end)
                        intron_list.append(intron)
                        intron_end = exon_start + 1
                elif string[1] == 'Term':
                    if string[2] == '+':
                        if not intron_number:
                            intron_number = exon_number - 1
                        if not intron_start:
                            intron_start = 1
                        intron_end = exon_start - 1
                        intron = Intron(gene=gene, number=intron_number, start=intron_start, end=intron_end)
                        intron_list.append(intron)
                        intron_number, intron_start, intron_end = 0, 0, 0
                    elif string[2] == '-':
                        intron_number = exon_number - 1
                        intron_end = exon_start - 1

            if line.startswith('>') and not protein_name:
                protein_name = line
            elif line.startswith('>') and protein_name:
                proteins[protein_name] = str(''.join(protein))
                protein_name = line
                protein = []
            elif protein_name:
                protein.append(line)

    proteins[protein_name] = str(''.join(protein))

    if intron_start or intron_end:
        if string[2] == '+':
            intron_end = sequence_length
        elif string[2] == '-':
            intron_start = sequence_length
        intron = Intron(gene=gene, number=intron_number, start=intron_start, end=intron_end)
        intron_list.append(intron)

    cdc_list = list(proteins.values())
    return cdc_list, exon_list, intron_list


def run_genscan(sequence: Optional[str] = None,
                sequence_file: Optional[str] = None,
                organism: str = "Vertebrate",
                exon_cutoff: Optional[float] = 1.00,
                sequence_name: Optional[str] = None) -> GenscanOutput:
    """
    Run Genscan prediction.
    Args:
    - sequence (str, optional): DNA sequence (upper or lower case, spaces/numbers ignored)
    - sequence_file (str, optional): Path to a file with DNA sequence (upper or lower case, spaces/numbers ignored)
    - organism (str): Organism type. Allowed values are 'Vertebrate', 'Arabidopsis' or 'Maize'
    - exon_cutoff (float, optional): Exon cutoff value
    - sequence_name (str, optional): Name of the DNA sequence
    Returns:
    - GenscanOutput: An object containing the status of the prediction, predicted coding sequences, introns, and exons
    Introns are considered to be regions lying between the initial, internal, and terminating exon
    Promoters and poly-A regions are not considered exons
    """

    allowed_organisms = ["Vertebrate", "Arabidopsis", "Maize"]
    if organism not in allowed_organisms:
        raise ValueError(f"Invalid value for organism. Allowed values are: {', '.join(allowed_organisms)}")

    url = "http://hollywood.mit.edu/cgi-bin/genscanw_py.cgi"
    params = {
        "-o": organism,
        "-e": exon_cutoff,
        "-n": sequence_name,
        "-p": "Predicted peptides only"
    }
    files = None
    if sequence_file:
        files = {"-u": open(sequence_file, "rb")}
    else:
        params["-s"] = sequence

    response = requests.post(url, data=params, files=files)

    if response.status_code == 200:
        soup = BeautifulSoup(response.content, 'html.parser')
        lines = soup.find('pre').text.split('\n')
        cds, exons, introns = parse_genscan_output(lines)
        return GenscanOutput(status="Success", cds_list=cds, intron_list=introns, exon_list=exons)
    else:
        return GenscanOutput(status="Failed", cds_list=None, intron_list=None, exon_list=None)
