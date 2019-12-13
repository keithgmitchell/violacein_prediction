#!/usr/bin/env bash
module load hmmer
module load muscle
module load prodigal


# update non producers variable for the list.txt and same for the hmm_seqs.txt

all_seqs_dir="../new_all_seqs"
non_prod_dir="../notes/non_producers"
hmm_results="../hmm_results"


# first lets create a master fasta file for the known producers
#for i in `cat ../notes/hmm_seqs.txt`; do cat "${all_seqs_dir}/${i}.fna.fasta" >> "${all_seqs_dir}/all_producers.fasta"; done

# next lets create an amino acid inference of the known procuders master fasta file
#cat "${all_seqs_dir}/all_producers.fasta" | prodigal > "${all_seqs_dir}/gene.all_producers.gbk" -a "${all_seqs_dir}/protein.all_producers.faa"

# then lets do the same thing for the known non producers
for i in `cat "${non_prod_dir}/list.txt"`; do cat "${non_prod_dir}/${i}.faa.fasta" >> "${non_prod_dir}/non_producers.fasta"; done
#cat "${non_prod_dir}/non_producers.fasta" | prodigal > "${non_prod_dir}/gene.non_producers.gbk" -a "${non_prod_dir}/protein.non_producers.faa"

# then create a master file of the known procuders and known non producers to test against the HMM
cat "${all_seqs_dir}/protein.all_producers.faa" "${non_prod_dir}/non_producers.fasta" > "${hmm_results}/protein.hmm_test.faa"

for i in A B C D E;

do
	prefix="hmm_matches_Vio${i}_0.00000001.m8.sfasta"
	muscle -in "${prefix}" -out "${prefix}.msa"
	hmmbuild "${prefix}.hmm" "${prefix}.msa"
	hmmsearch "${prefix}.hmm" "${hmm_results}/protein.hmm_test.faa" > "${prefix}.results"
	for i in `cat "${non_prod_dir}/list.txt"`; do cat "${prefix}.results" | grep $i | head -1; done
	for i in `cat ../notes/hmm_seqs.txt`; do cat "${prefix}.results" | grep $i | head -1; done
done
