import sys
from Bio import SeqIO

def genbank2fasta(in_gbk, out_fasta):
    SeqIO.convert(in_gbk, "genbank", out_fasta, "fasta")

if __name__ == "__main__":
    in_gbk = sys.argv[1]
    out_fasta = sys.argv[2]
    genbank2fasta(in_gbk, out_fasta)
