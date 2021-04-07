import os
from Bio import SeqIO
from Bio import Entrez

def download(id):

    '''
        Downloads .fasta file from NCBI database

        Input:
            id: string - NCBI id of .fasta file

        Output:
            string: path to downloaded file
    '''

    Entrez.email = ""

    net_handle = Entrez.efetch(
        db="nucleotide", id=id, rettype="fasta", retmode="text"
    )

    out_handle = open(id + ".fasta", "w")
    out_handle.write(net_handle.read())

    out_handle.close()
    net_handle.close()

    return id + ".fasta"
