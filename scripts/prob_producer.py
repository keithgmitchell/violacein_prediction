import argparse

parser = argparse.ArgumentParser(
    description='Calculate the prob of a producer based on a list of hmm result files')
parser.add_argument('-i', '--input', type=str, help="List of hmm files from build and test hmm script")

args = parser.parse_args()
input_files = args.input
input_file_list = input_files.split(',')

dict_of_e_values = {}
for i in input_file_list:
	with open(i, 'r') as current_file:
		for line in current_file:
			gene = line.split()[8]
			microbe = gene.split('__')[0]
			e_value = line.split()[0]
			if microbe in dict_of_e_values.keys():
				dict_of_e_values[microbe].append(e_value)
			else:
				dict_of_e_values[microbe] = [e_value, ]

print (dict_of_e_values)
# calculate the p_value based on the constructed dictionary
number_of_genes = len(input_file_list)
for organism in dict_of_e_values.keys():
	offset = number_of_genes - len(dict_of_e_values[organism])
	prob = 1
	# small e_value means high chance of producing... therefore 1-e-value for all genes to calc intersection of them all occuring
	for i in dict_of_e_values[organism]:
		prob *= 1-float(i)
	# for those that did not return a hit lets assume some base value of them being a producer (0.05)
	if offset > 0:
		prob *= float(0.05)**offset
	print(organism, prob)
