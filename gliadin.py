import pandas as pd
import requests
import io
def sequence_breaker(sequence):
    """Breaks input sequence into a string of 15 mer peptides.
      Args:
      protein sequence
      Returns:
      dictionary of 15 mer peptides."""
    peptide_dict={}
    n=len(sequence)

    for i in range (0 ,n-15 ):
        label=f"peptide_{i+1}"
        peptide=sequence[i:i+15]
        peptide_dict[label]=peptide
    return peptide_dict
def protein_fetcher(uniprot_id):
    """Fetches a sequence from uniprot using rest API and converts it into a string format.
    Args: 
    Uniprot ID
    Returns:
    protein sequence in string format""" 
    response = requests.get(f"https://www.uniprot.org/uniprot/{uniprot_id}.fasta")
    uniprot_seq=response.text
    lines = uniprot_seq.split("\n")
    uniprot_lines = lines[1:]
    uniprot_sequencef = "".join(uniprot_lines)
    return uniprot_sequencef
gliadin_seq=protein_fetcher("Q402I5") 
gliadin_peptides=sequence_breaker(gliadin_seq)
peptide_list=list(gliadin_peptides.items())
df=pd.DataFrame.from_dict(gliadin_peptides, orient="index")
df.index.name="Peptides"
df = df.rename(columns={0: "Sequence"})
print(df)
def list_chomper (peptide_list, i, chunk_size):
    """ Takes a list and chomps it in predefined chunks
       Args:   
       list
       Returns: 
       A chunk of predefined size"""
    chunk=peptide_list[i:i+chunk_size]
    return chunk
def mhcii_predictor(fasta_sequence):
    """ Takes a single fasta formatted string and predicts the mhcii binding affinity
         Args:   
         list
         Returns: 
         A list of predicted binding affinities"""
    x= requests.post( "https://tools-cluster-interface.iedb.org/tools_api/mhcii/", data={"sequence_text":fasta_sequence,"method":'netmhciipan_el', "allele": 'HLA-DQA1*05:01/DQB1*02:01,HLA-DQA1*03:01/DQB1*03:02', "length":'asis'})
    prediction_df=pd.read_table(io.StringIO(x.text))
    return prediction_df
results=[]
for i in range(0, len(peptide_list), 8):
    fasta_list = []
    peptide_chunk = list_chomper(peptide_list, i, 8)
    print(f"Processing chunk starting at peptide {i}")
    for peptide in peptide_chunk:
        fasta_list.append(f">{peptide[0]}\n{peptide[1]}")
    fasta_sequence = "\n".join(fasta_list)
    prediction = mhcii_predictor(fasta_sequence)
    results.append(prediction)
table = pd.concat(results)
print(table)
strong_binders = table[table['rank'] < 10]
print(strong_binders)
