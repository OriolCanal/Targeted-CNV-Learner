

from pyfaidx import Fasta

from Bio import SeqIO
fasta_file = "/home/ocanal/REF_DIR/hg19/ucsc.hg19.fasta"
# for record in SeqIO.parse(fasta_file, "fasta"):
#     print(record.id)

fasta = Fasta(fasta_file)
sequence = str(fasta[1][1234765:1324767])

if not sequence:
    ValueError(f"Sequence not foud for  in {fasta_file}")

fasta.close()
print(sequence)