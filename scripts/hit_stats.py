from Bio import SeqIO
import argparse

parser = argparse.ArgumentParser(
    description='Process some fasta file to return the stats on the matches at some e-value')
parser.add_argument('-i', '--input', type=str, help="This in the fasta file to get the stats on")
parser.add_argument('-r', '--reference', type=str, help="This is the list of files to compare the stats on")

args = parser.parse_args()
input_file = args.input
ref_file = args.reference

# Parse and count the reads from each of the genomes
dict = {}
fasta_sequences = SeqIO.parse(open(input_file), 'fasta')
for fasta in fasta_sequences:
    name, sequence = fasta.id, str(fasta.seq)
    split = name.split('__')
    if split[0] in dict.keys():
        dict[split[0]] += 1
    else:
        dict[split[0]] = 1

# Construct a list of them all, and print those that have more then 1
hit_list = []
multi_hits = []
print("+++ MULTI MATCHES +++")
for item in dict:
    hit_list.append(item)
    if dict[item] > 1:
        multi_hits.append(item)
        print(dict[item], item)
print("+++ TOTAL MULTI MATCHES ++++")
print(len(multi_hits), " ".join(multi_hits))

# Get list from the ls *.fasta > to compare to? Which arent in both lists
ref_list = []
with open(ref_file, 'r') as ref:
    for line in ref:
        line = line.strip('\n')
        line = line.strip('.fna.fasta')
        ref_list.append(line)

hit_ref = list(set(hit_list) - set(ref_list))
ref_hit = list(set(ref_list) - set(hit_list))

# print(len(hit_ref), " ".join(hit_ref))
print("+++ REFERENCE - HITS +++")
print(len(ref_hit), " ".join(ref_hit))

