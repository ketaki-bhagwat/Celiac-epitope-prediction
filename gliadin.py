import pandas as pd
import requests
import io

# ── Functions ──────────────────────────────────────────────

def protein_fetcher(uniprot_id):
    """Fetches a sequence from UniProt using REST API and converts it into string format.
    Args:
        uniprot_id: UniProt accession ID
    Returns:
        Protein sequence as a string
    """
    response = requests.get(f"https://www.uniprot.org/uniprot/{uniprot_id}.fasta")
    uniprot_seq = response.text
    lines = uniprot_seq.split("\n")
    uniprot_lines = lines[1:]
    uniprot_sequencef = "".join(uniprot_lines)
    return uniprot_sequencef

def sequence_breaker(sequence):
    """Breaks input sequence into a dictionary of 15-mer peptides.
    Args:
        sequence: Protein sequence as a string
    Returns:
        Dictionary of 15-mer peptides
    """
    peptide_dict = {}
    n = len(sequence)
    for i in range(0, n - 15):
        label = f"peptide_{i+1}"
        peptide = sequence[i:i+15]
        peptide_dict[label] = peptide
    return peptide_dict

def list_chomper(peptide_list, i, chunk_size):
    """Takes a list and returns a chunk of predefined size.
    Args:
        peptide_list: List of peptides
        i: Starting index
        chunk_size: Size of chunk
    Returns:
        A chunk of predefined size
    """
    chunk = peptide_list[i:i+chunk_size]
    return chunk

def mhcii_predictor(fasta_sequence):
    """Takes a FASTA-formatted string and predicts MHC-II binding affinity.
    Args:
        fasta_sequence: FASTA-formatted string
    Returns:
        DataFrame of predicted binding affinities
    """
    x = requests.post(
        "https://tools-cluster-interface.iedb.org/tools_api/mhcii/",
        data={
            "sequence_text": fasta_sequence,
            "method": 'netmhciipan_el',
            "allele": 'HLA-DQA1*05:01/DQB1*02:01,HLA-DQA1*03:01/DQB1*03:02',
            "length": 'asis'
        }
    )
    prediction_df = pd.read_table(io.StringIO(x.text))
    return prediction_df

# ── Execution ──────────────────────────────────────────────

gliadin_seq = protein_fetcher("Q402I5")
gliadin_peptides = sequence_breaker(gliadin_seq)
peptide_list = list(gliadin_peptides.items())

df = pd.DataFrame.from_dict(gliadin_peptides, orient="index")
df.index.name = "Peptides"
df = df.rename(columns={0: "Sequence"})
print(df)

results = []
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
