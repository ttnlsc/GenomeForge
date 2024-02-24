# GenomeForge

In this virtual temple, we have gathered sacred texts and tools developed by our Magos to explore and manipulate the code of life.
To successfully carry out various manipulations, it is necessary to perform certain canticles.

## Installation

Follow the steps below to install GenomeForge:

1. Download the GenomeForge folder to your local computer.
2. Make sure the path to the GenomeForge folder is specified in the system PATH variable. 
3. Import the GenomeForge main script magos_biologis.py as a module

```python
import magos_biologis as mb
```

## Datascrolls

### DNA RNA Canticles

The `dna_rna_canticles` offers rites for the sacred manipulations of DNA and RNA sequences. It grants the ability to perform rites on single or on multiple sequences, accepting a variable number of sequences (*str*) as input and the designation of the desired rite to be executed (always the final designation, *str*, as seen in the usage example). The designated rite is then executed on all the sequences provided.

**Catalogue of Rites:**

1. `transcribe` — transcribes a DNA sequence into an RNA sequence, and also performs reverse transcription. 

2. `reverse` — returns the reversed sequence

3. `complement` — returns the complementary sequence

4. `reverse_complement` — returns the reverse complementary sequence
5. `lenght` — returns the length of the sequence
6. `gc_content` — calculates the proportion of G and C nucleotides in the sequence and returns the result as a fraction of one

**Usage example:**

```python
mb.run_dna_rna_canticles('ATG', 'transcribe') # 'AUG'
mb.run_dna_rna_canticles('ATG', 'reverse') # 'GTA'
mb.run_dna_rna_canticles('AtG', 'complement') # 'TaC'
mb.run_dna_rna_canticles('ATg', 'reverse_complement') # 'cAT'
mb.run_dna_rna_canticles('ATG', 'aT', 'reverse') # ['GTA', 'Ta']
mb.run_dna_rna_canticles('AActG', 'ATCT', 'lenght') # [5, 4]
mb.run_dna_rna_canticles('AgCT', 'gc_content') # 0.5
```

**Note:**

- This canticle saves the character case

- This canticle **only** for nucleic acid sequences. Should the data be offered in a manner inconsistent with the holy protocol, an apt error shall manifest as divine retribution.

### Protein Canticles

This Omnissiah-approved endeavor encompasses the sacred program known as `protein_canticles`. Within its sanctified code resides the mighty `run_protein_canticles` function, a conduit for channeling protein sequences and the divine actions to be invoked upon them. Enter the name of the action (*str*), the protein sequence or path to the protein sequence file (*str*), and additional arguments if necessary.

**Catalogue of Rites:**

1. `find_sites` - finds positions of given sites.

*Args*:

- seq (*str*): sequence or path to file in fasta format with sequences to be checked
- *args (*str*): sites to be found
- is_one_based (*bool*): whether result should be 0- (False) or 1-indexed (True). Default False

*Returns*:

- *dict*: dictionary of sites as keys and lists of positions for the site where it's been found. In case you entered the path to a fasta file with multiple sequences, you will get back a dictionary containing the names of sequences as keys and the dictionary of sites as keys and lists of positions.

2. `get_protein_rnas_number` - get number of all possible RNA's for a given protein.

*Args*:

- seq (*str*): sequence or path to file in fasta format with sequences to be checked

*Returns*:

- *int* or *dict*: number of possible RNA's for sequence. In case you entered the path to a fasta file with multiple sequences, you will get back a dictionary containing the names of sequences as keys and number of possible RNA's.

3. `get_protein_rnas` - returns list of all possible RNA's from which can serve as matrix for protein synthesis.

**WARNING**: can be computationally intensive on longer sequences, will NOT start unless check_if_user_conscious is True!

*Args*:

- seq (*str*): sequence or path to file in fasta format with sequences to be checked
- check_if_user_conscious (*bool*): checks user's consciousness. Default False

*Returns*:

- *list*: list of possible RNA's as string. In case you entered the path to a fasta file with multiple sequences, you will get back a dictionary containing the names of sequences as keys and list of possible RNA's

4. `get_frameshift_proteins` - returns list of all possible proteins from all possible frames in peptide.

**WARNING**: can be computationally intensive on longer sequences, will NOT start unless check_if_user_conscious is True!

*Args*:

- seq (*str*): sequence or path to file in fasta format with sequences to be checked
- check_if_user_conscious (*bool*): checks user's consciousness. Default False
- is_stop_codon_termination_enabled (*bool*): terminate translation when reached stop-codon. Default False.

*Returns*:

- *dict*: dictionary of lists of all possible frames proteins. In case you entered the path to a fasta file with multiple sequences, you will get back a dictionary containing the names of sequences as keys and the dictionary of lists of all possible frames proteins.

5. `get_length_of_protein` - calculates the length of a protein.

*Args*:

- seq (*str*): sequence or path to file in fasta format with sequences to be checked

*Returns*:

- *int*: sequence length. In case you entered the path to a fasta file with multiple sequences, you will get back a dictionary containing the names of sequences as keys and sequences lengths.

6. `count_aa` - counts the number of given or all amino acids in a protein sequence.

*Args*:

- seq (*str*): sequence or path to file in fasta format with sequences to count amino acids
- aminoacids (*str*): which amino acids to count in sequence. If you want to count all amino acids in the whole sequence, you can provide empty string to this argument or just don't provide this keyword

*Returns*:

- *dict*: dictionary with amino acids and its count. In case you entered the path to a fasta file with multiple sequences, you will get back a dictionary containing the names of sequences as keys and the dictionary with amino acids and its count.

7. `get_fracture_of_aa` - calculates the fracture or percentage of amino acids in a protein sequence.

*Args*:

- seq (*str*): sequence or path to file in fasta format with sequences in which you need to calculate the fracture of amino acids
- show_as_percentage (*bool*): change it to True, if you want to get results with percentages
- aminoacids (*str*): the fracture of which amino acids to count in the sequence

*Returns*:

- *dict*: dictionary with amino acids and its fracture or percentage. In case you entered the path to a fasta file with multiple sequences, you will get back a dictionary containing the names of sequences as keys and the dictionary with amino acids and its  fracture.

8. `calculate_protein_mass` - calculates the molecular mass of a protein based on its amino acid sequence and a dictionary of amino acid masses.

*Args*:

- sequence (*str*): sequence or path to file in fasta format with sequences to be processed

*Returns*:

- *float* or *dict*: The molecular mass of a protein in atomic mass units, rounded to the third decimal place.  In case you entered the path to a fasta file with multiple sequences, you will get back a dictionary containing the names of sequences as keys and the molecular mass of a protein in atomic mass units.

9. `get_atomic_mass` - calculates the molecular mass of a biological molecule, primarily an amino acid, based on a simple chemical formula.

*Args*:

- chem (*str*): String representing a simple chemical formula, e.g. C2H5OH

*Returns*:

- *float*: Molecular mass of a biological molecule in atomic mass units, rounded to the third decimal place..

10. `convert_aa_name` - converts a sequence of one-letter amino acid codes to three-letter designations.

*Args*:

- sequence (*str*): sequence or path to file in fasta format with sequences with one-letter amino acid codes
- sep (*str*, optional): Separator between three-letter amino acid designations. There is no delimiter by default.
- use_default_register( *bool*, optional): Determines whether to preserve letter case in three-letter designations. If True, the letters will be converted to upper or lower case depending on the case of the depending on the case of the one-letter code. The default is False.

*Returns*:

- *str*: string of three-letter amino acid designations separated by the specified delimiter. In case you entered the path to a fasta file with multiple sequences, you will get back a dictionary containing the names of sequences as keys and the string of three-letter amino acid designations .

**Usage example**

```python
mb.run_protein_canticles('find_sites', '../ins.fasta', 'LLALL', is_one_based=True) -> {'LLALL': [10]}
mb.run_protein_canticles('get_protein_rnas_number', 'MALWMR') -> 144
mb.run_protein_canticles('get_length_of_protein', '../ins.fasta') -> 110
mb.run_protein_canticles('count_aa', '../ins.fasta', 'DRG') -> {'D': 2, 'R': 5, 'G': 12}
mb.run_protein_canticles('get_fracture_of_aa', '../ins.fasta', aminoacids='DRG', show_as_percentage=True) -> {'D': 1.82, 'R': 4.55, 'G': 10.91}
mb.run_protein_canticles('calculate_protein_mass', '../ins.fasta') -> 11980.804
mb.run_protein_canticles('get_atomic_mass', 'C5H5N5') -> 135.128
mb.run_protein_canticles('convert_aa_name', 'MALWMR', sep='~') -> 'Met~Ala~Leu~Trp~Met~Arg'
```

**Note:**

- The program devoutly scrutinizes protein sequences comprised of the **20 canonical amino acids**.

- Should the data be offered in a manner inconsistent with the holy protocol, an apt error shall manifest as divine retribution.

## FASTQ Canticles

The `fastq_canticles` is a sacred tool endowed with the ability to perform complex Fastq sequence filtering rituals. At its core, this divine utility takes four sacred offerings: **seqs**, **gc_bounds**, **length_bounds**, and **quality_threshold**.

**Catalogue of Rites:**

1. `filter_fastq` -  filters fastq sequences based on the GC-content, length and quality parameters.

*Args*:

- seqs (*dict*): dictionary containing fastq sequences. Each sequence is represented as a dictionary item with the sequence name as the key and a tuple containing the sequence and its quality as the value. 
- gc_bounds (*tuple* or *int*): GC content range (in percentages) for filtering. The default value is (0, 100), meaning all reads are retained.  If you pass a single number to the argument, it is assumed to be an upper bound.
- length_bounds (*tuple* or *int*): length range for filtering. The default value is (0, 2^32), meaning all reads are retained. If you pass a single number to the argument, it is assumed to be an upper bound.
- quality_threshold (*int*): Threshold value for average read quality filtering. The default value is 0 (Phred33 scale). Reads with an average quality score below this threshold are filtered out.

*Returns*:

- dict: dictionary with filtered reads

**Usage example**

```python
mb.filter_fastq(EXAMPLE_FASTQ, gc_bounds=(40, 50), length_bounds=85, quality_threshold=30) ->
{'@SRX079804:1:SRR292678:1:1101:149302:149302': ('TAGGGTTGTATTTGCAGATCCATGGCATGCCAAAAAGAACATCGTCCCGTCCAATATCTGCAACATACCAGTTGGTTGGTA',
  '@;CBA=:@;@DBDCDEEE/EEEEEEF@>FBEEB=EFA>EEBD=DAEEEEB9)99>B99BC)@,@<9CDD=C,5;B::?@;A'),
 '@SRX079804:1:SRR292678:1:1101:170868:170868': ('CTGCCGAGACTGTTCTCAGACATGGAAAGCTCGATTCGCATACACTCGCTGAGTAAGAGAGTCACACCAAATCACAGATT',
  'E;FFFEGFGIGGFBG;C6D<@C7CDGFEFGFHDFEHHHBBHHFDFEFBAEEEEDE@A2=DA:??C3<BCA7@DCDEG*EB')}
```

![Magos_Biologis](./img/img1.png)

In the sacred repository of GenomForge, you have discovered a set of invaluable tools bestowed upon you by the Omnissiah.  As an adept of the digital arts, you can now invoke the `dna_rna_canticles` for handling genetic scripts, chant the verses of `protein_canticles` to fathom the mysteries of proteins, and, with the `fastq_canticles` utility, undertake intricate purification rituals on Fastq sequences.

May these digital sigils and functions guide you on your quest for knowledge, unlocking the secrets of life encrypted within the digital scrolls.
