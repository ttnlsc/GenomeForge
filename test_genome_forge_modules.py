import os
import pytest
from magos_biologis import DNASequence, AminoAcidSequence, NucleicAcidSequence
from bio_files_processor import convert_multiline_fasta_to_oneline
class TestDNASequence:
    """
    Test methods of the DNASequence class.
    """
    def test_complement(self):
        """
        Test the complement method of the DNASequence class.
        """
        dna_seq = DNASequence("ATGC")
        assert str(dna_seq.complement(DNASequence.complement_map)) == "TACG"

    def test_gc_content(self):
        """
        Test the gc_content method of the DNASequence class.
        """
        dna = DNASequence("ATGC")
        gc_content = dna.gc_content(percentage=True)
        assert gc_content == 50.0

@pytest.mark.parametrize(
    "input_sequence, expected",
    [
        ("ATGC", "AUGC"),
        ("atgc", "augc"),
        ("TTTGG", "UUUGG")
    ]
)
def test_dna_sequence_transcribe(input_sequence, expected):
    """
    Test the transcribe method of the DNASequence class.
    """
    dna_seq = DNASequence(input_sequence)
    rna_seq = dna_seq.transcribe()
    assert str(rna_seq) == expected


class TestAminoAcidSequence:
    """
    Test methods of the AminoAcidSequence class.
    """
    def test_count_aa(self):
        """
        Test the count_aa method of the AminoAcidSequence class.
        """
        amino_acid_seq = AminoAcidSequence("LISITSA")
        count = amino_acid_seq.count_aa()
        assert count == {'L': 1, 'A': 1, 'T': 1, 'S': 2, 'I': 2}

@pytest.fixture(scope="module")
def input_fasta_path():
    return "test_input.fasta"

@pytest.fixture(scope="module")
def output_fasta_path():
    return "test_output.fasta"

def test_convert_multiline_fasta_to_oneline(input_fasta_path, output_fasta_path):
    with open(input_fasta_path, "w") as f:
        f.write(">Seq1\nACGT\nTGCA\n>Seq2\nGGGG\nCCCC")

    convert_multiline_fasta_to_oneline(input_fasta_path, output_fasta_path)

    assert os.path.exists(output_fasta_path)

    with open(output_fasta_path, "r") as f:
        content = f.read()
    assert content == ">Seq1\nACGTTGCA\n>Seq2\nGGGGCCCC"

    os.remove(input_fasta_path)
    os.remove(output_fasta_path)


def test_complement_method_not_implemented():
    sequence = NucleicAcidSequence("ATGC")
    with pytest.raises(NotImplementedError):
        sequence.complement(None)