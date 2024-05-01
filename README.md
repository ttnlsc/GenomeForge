# GenomeForge

In this virtual temple, we have gathered sacred texts and tools to explore and manipulate the code of life developed by our Magos during the Python course (2023 - 2024) at the Bioinformatics Institute. To successfully carry out various manipulations, it is necessary to perform certain canticles. 

## Installation

Follow the steps below to install GenomeForge:

1. Download the GenomeForge folder to your local computer.
2. Make sure the path to the GenomeForge folder is specified in the system PATH variable. 

## Datascrolls

### magos_biologis.py

**Classes of Power**:

* *DNASequence*
* *RNASequence*
* *AminoAcidSequence*

**Functions of Invocation**

* *filter_fastq*: Filters FASTQ sequences based on GC-content, length, and quality parameters.
* *telegram_logger*: Decorator function for logging and sending messages to Telegram.

**Rituals of Revelation**

- *run_genscan*: Runs Genscan prediction, returning predicted coding sequences, introns, and exons.

Invoke these scripts and rituals with reverence and caution.

### bio_files_processor.py

Within the hallowed halls of our digital sanctum, we wield the sacred scripts bestowed upon us by the Omnissiah. These incantations allow us to unravel the mysteries of the genetic code and shape the very essence of life itself.

**Scripts of Transformation**

* *convert_multiline_fasta_to_oneline*
* *select_genes_from_gbk_to_fasta*
* *change_fasta_start_pos*
* *parse_blast_output*

**Sacred Rituals**

- *OpenFasta*: A context manager for reading FASTA files, allowing seamless iteration over records.

Invoke these scripts and rituals with reverence and caution.

## custom_random_forest.py

**Arcane Scrolls of Forest Conjuration**

Deep within the forest of algorithms, where the trees whisper secrets of classification, lies the esoteric knowledge of the RandomForestClassifierCustom. This mystical construct, a creation of human and machine collaboration, holds the power to unveil patterns hidden within data.

**Invocation Parameters**

- **n_estimators** (int, default=10): The number of trees in the forest.
- **max_depth** (int or None, default=None): The maximum depth of each tree.
- **max_features** (int or None, default=None): The maximum number of features to consider when looking for the best split.
- **random_state** (int or None, default=None): Controls the randomness of the forest.

**Manifestations**

1. **fit(X, y, n_jobs)**
   - Fit the RandomForestClassifierCustom to the training data.
2. **predict(X)**
   - Predict class labels for the input samples.
3. **predict_proba(X, n_jobs)**
   - Predict class probabilities for the input samples.

Invoke these incantations with care, for they possess the ability to unveil the hidden truths encoded within the data.

![Magos_Biologis](./img/img1.png)

In the sacred repository of GenomForge, you have discovered a set of invaluable tools bestowed upon you by the Omnissiah.  As an adept of the digital arts, you can now invoke the `dna_rna_canticles` for handling genetic scripts, chant the verses of `protein_canticles` to fathom the mysteries of proteins, and, with the `fastq_canticles` utility, undertake intricate purification rituals on Fastq sequences.

May these digital sigils and functions guide you on your quest for knowledge, unlocking the secrets of life encrypted within the digital scrolls.

Praise the Omnissiah! 
