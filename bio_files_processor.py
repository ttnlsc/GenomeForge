import os


def convert_multiline_fasta_to_oneline(input_fasta: str, output_fasta: str = None) -> None:
    """
    Converts a multiline FASTA file to a oneline FASTA file.

    Args:
    - input_fasta (str): path to the input FASTA file.
    - output_fasta (str): path to the output oneline FASTA file. If not provided, it will be generated
      using the input file name.

    Returns:
    - None: The function doesn't return a value but writes the oneline FASTA to the output file.
    """
    if output_fasta is None:
        output_fasta = 'oneline_' + os.path.basename(input_fasta)
    
    with open(input_fasta, mode='r') as infile:
        multiline = infile.readlines()
        
    output_str = []
    name = ''
    record = []
    
    for line in multiline:
        line = line.strip()
        if line.startswith('>') and not name:
            name = line
        elif not line.startswith('>'):
            record.append(line)
        elif name and record:
            sequence = ''.join(record)
            output_str.append(name)
            output_str.append(sequence)
            name = line
            record = []
            
    if name and record:
        sequence = ''.join(record)
        output_str.append(name)
        output_str.append(sequence)
    
    with open(output_fasta, mode='w') as outfile:
        outfile.write('\n'.join(output_str))
    return None


def select_genes_from_gbk_to_fasta(input_gbk: str, genes_to_find: list, n_before: int = 1, n_after: int = 1,
                                   output_fasta: str = None) -> None:
    """
    Extracts the gene name and protein sequences from the GenBank file and creates a FASTA file containing 
    the gene names and protein sequences flanking the specified genes 

    Args:
    - input_gbk (str): path to the GenBank file.
    - genes_to_find (list of str): list of gene names (str) to extract from the GenBank file.
    - n_before (int): number of genes before the target gene to include in the output. Default value is 1.
    - n_after (int): number of genes after the target gene to include in the output. Default value is 1.
    - output_fasta (str, optional): path to the output FASTA file. If not provided, a default name is generated.

    Returns:
    - None: The function writes the selected gene and protein sequences to the output FASTA file.
    """
    if output_fasta is None:
        name = os.path.basename(input_gbk).split('.')[0]
        output_fasta = name + '.fasta'
    
    with open(input_gbk, mode='r') as file:
        gbk = file.readlines()
    gene_protein_list = []
    current_gene = ''
    current_sequence = []
    for line in gbk:
        line = line.strip()
        if line.startswith('/gene'):
            current_gene = line.split('"')[1]
        elif current_gene and line.startswith('/translation='):
            protein = line.split('"')[1]
            current_sequence.append(protein)
        elif current_gene and current_sequence and not line.startswith('/') and not line.endswith('"'):
            current_sequence.append(protein)
        elif current_gene and current_sequence and line.endswith('"'):
            protein = line.replace('"', '')
            current_sequence.append(protein)
            current_protein = ''.join(current_sequence)
            gene_protein_list.append([current_gene, current_protein])
            current_gene = ''
            current_sequence = []
    if current_gene and current_sequence:
        current_protein = ''.join(current_sequence)
        gene_protein_list.append([current_gene, current_protein])
    idxs = []
    for gene_name in genes_to_find:
        for idx, (gene, protein) in enumerate(gene_protein_list):
            if gene_name in gene:
                idxs.append(idx)
    flanks = []
    for idx in idxs:
        for i in range(n_before, 0, -1):
            flanks.append(idx - i)
        for i in range(1, n_after + 1):
            flanks.append(idx + i)
    
    selected_records = [gene_protein_list[i] for i in flanks]
    with open(output_fasta, mode='w') as outfile:
        for gene, protein in selected_records:
            outfile.write(gene + '\n' + protein + '\n')
    return None
