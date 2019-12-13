import itertools
import sys
import argparse

parser = argparse.ArgumentParser(
    description='Calculate interset of all genes at some e-value cutoff')
parser.add_argument('-i', '--input', type=str, help="This in the fasta file to get the stats on")

args = parser.parse_args()
input_file = args.input

master_list = []
for line in open(input_file, 'r'):
    if "+++ REFERENCE - HITS +++" in line:
        list = input_file.readline().split(" ")[1:]
        master_list.append(list)
        list = []

for item in master_list:
    print(item)

result = (set(master_list[0]) & set(master_list[1]) & set(master_list[2]) & set(master_list[3]) & set(master_list[4]))
# print ("A: ", len(master_list[0]))
# print ("B: ", len(master_list[1]))
# print ("C: ", len(master_list[2]))
# print ("D: ", len(master_list[3]))
# print ("E: ", len(master_list[4]))

for i in itertools.combinations([[0,'A'],[1,'B'],[2,'C'],[3,'D'],[4,'E']], 2):
    iter = (set(master_list[i[0][0]]) & set(master_list[i[1][0]]))
    print (i[0][1] + i[1][1] + ' :', len(iter))

print("ABCDE: ", len(result), result)

result2 = (set(master_list[0]) & set(master_list[1]) & set(master_list[2]) & set(master_list[4]))
print("ABCE: ", len(result2), result2)

for i in itertools.combinations([[0,'A'],[1,'B'],[2,'C'],[3,'D'],[4,'E']], 3):
    iter = (set(master_list[i[0][0]]) & set(master_list[i[1][0]]) & set(master_list[i[2][0]]))
    print (i[0][1] + i[1][1] + i[2][1] + ' :', len(iter))