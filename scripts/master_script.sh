#!/usr/bin/env bash

module load diamond
module load prodigal
# TODO be sure the snakemake conda is not active

# First lets build the inferred proteins
#cat *.fasta | prodigal > gene.coords.gbk -a protein.translations.faa

# Next we want to make a database from the proteins
#diamond makedb --in protein.translations.faa -d master

for num in 0.1 0.01 0.000001 0.00000000001 0.000000000000000000001;
do
	for gene in VioA VioB VioC VioD VioE;
	do
		matches="../matches_output/matches_${gene}_${num}.m8"
		echo "$matches"

		echo query the vio gene with certain threshold against the inferred proteins from the all fasta files
		diamond blastx -d master -q "../Vio_genes_NCBI/${gene}_Chromobacterium_violaceum.fasta" -o $matches --max-target-seqs 0 --evalue $num -f sam

		echo take the mateches format and create a fasta file from it
		tail -n +6 $matches | awk '{print ">"$3,$1"\n"$10}' > "${matches}.sfasta"

		echo check the results from the matches by looking at the hit stats
		python ../notes/hit_stats.py -i "${matches}.sfasta" -r ../notes/ref_genomes.txt > "../matches_output/${gene}_${num}.output"

		echo finally lets output the results from the hit stats based on all matches and vio genes combinations
		echo --- >> final_results.txt
		echo "${matches}.sfasta", "../matches_output/${gene}_${num}.ouput" >> final_results.txt
		echo --- >> final_results.txt
		cat "../matches_output/${gene}_${num}.output" >> final_results.txt
	done
done


