

# Part of a master script
- Create the database of the inferred proteins:

    `cat *.fasta | prodigal > gene.coords2.gbk -a protein.translations.faa`

- Make the database of proteins:
    `diamond makedb --in protein.translations.faa -d all_proteins` 

- Align the gene of interest (Fasta) to the inferred protein database created above with certain e-value cutoffs:
    ```
    diamond blastx -d nr2 -q ../Vio_genes_NCBI/VioA_Chromobacterium\ violaceum.fasta -o matches_vioA.m8 --max-target-seqs 0 --evalue 0.0000000001
    diamond blastx -d nr2 -q ../Vio_genes_NCBI/VioB_Chromobacterium\ violaceum.fasta -o matches_vioB.m8 --max-target-seqs 0 --evalue 0.0000000001
    diamond blastx -d nr2 -q ../Vio_genes_NCBI/VioC_Chromobacterium\ violaceum.fasta -o matches_vioC.m8 --max-target-seqs 0 --evalue 0.0000000001
    diamond blastx -d nr2 -q ../Vio_genes_NCBI/VioD_Chromobacterium\ violaceum.fasta -o matches_vioD.m8 --max-target-seqs 0 --evalue 0.0000000001
    diamond blastx -d nr2 -q ../Vio_genes_NCBI/VioE_Chromobacterium\ violaceum.fasta -o matches_vioE.m8 --max-target-seqs 0 --evalue 0.0000000001
    ```
    ```
    diamond blastx -d nr2 -q ../Vio_genes_NCBI/VioA_Chromobacterium\ violaceum.fasta -o matches_vioA.m8 --max-target-seqs 0 --evalue 0.00001
    diamond blastx -d nr2 -q ../Vio_genes_NCBI/VioB_Chromobacterium\ violaceum.fasta -o matches_vioB.m8 --max-target-seqs 0 --evalue 0.00001
    diamond blastx -d nr2 -q ../Vio_genes_NCBI/VioC_Chromobacterium\ violaceum.fasta -o matches_vioC.m8 --max-target-seqs 0 --evalue 0.00001
    diamond blastx -d nr2 -q ../Vio_genes_NCBI/VioD_Chromobacterium\ violaceum.fasta -o matches_vioD.m8 --max-target-seqs 0 --evalue 0.00001
    diamond blastx -d nr2 -q ../Vio_genes_NCBI/VioE_Chromobacterium\ violaceum.fasta -o matches_vioE.m8 --max-target-seqs 0 --evalue 0.00001
    ```
    ```
    diamond blastx -d nr2 -q ../Vio_genes_NCBI/VioA_Chromobacterium\ violaceum.fasta -o matches_vioA.m8 --max-target-seqs 0 --evalue 0.001
    diamond blastx -d nr2 -q ../Vio_genes_NCBI/VioB_Chromobacterium\ violaceum.fasta -o matches_vioB.m8 --max-target-seqs 0 --evalue 0.001
    diamond blastx -d nr2 -q ../Vio_genes_NCBI/VioC_Chromobacterium\ violaceum.fasta -o matches_vioC.m8 --max-target-seqs 0 --evalue 0.001
    diamond blastx -d nr2 -q ../Vio_genes_NCBI/VioD_Chromobacterium\ violaceum.fasta -o matches_vioD.m8 --max-target-seqs 0 --evalue 0.001
    diamond blastx -d nr2 -q ../Vio_genes_NCBI/VioE_Chromobacterium\ violaceum.fasta -o matches_vioE.m8 --max-target-seqs 0 --evalue 0.001
    ```

- Does this match the number of expected where every genome produces a hit?

- Lets check the results:
    `for file in matches* ; do echo $file ; cat $file | wc -l ; done`
    - 0.001
    ```
        matches_vioA.m8
        55
        matches_vioB.m8
        53
        matches_vioC.m8
        116
        matches_vioD.m8
        73
        matches_vioE.m8
        53
    ```
    - 0.00001
    ```
        matches_vioA.m8
        53
        matches_vioB.m8
        53
        matches_vioC.m8
        82
        matches_vioD.m8
        73
        matches_vioE.m8
        53
    ```
    - 0.0000000001
    ```
        matches_vioA.m8
        53
        matches_vioB.m8
        53
        matches_vioC.m8
        60
        matches_vioD.m8
        73
        matches_vioE.m8
        53
    ```
    
    
- Now in order to get the sequence alginments add `-f sam` to the end of the diamond arguments for "Align the gene of interest (Fasta) to the inferred protein database created above with certain e-value cutoffs"
    `diamond blastx -d nr2 -q ../Vio_genes_NCBI/VioA_Chromobacterium\ violaceum.fasta -o matches_vioA.m8 --max-target-seqs 0 --evalue 0.001 -f sam`
 
 
 
 
- Then create a fasta file for muscle->fastree
    `tail -n +6 matches_vioA.m8 | awk '{print ">"$3,$1"\n"$10}' > matches_vioA.m8.fasta`
    `module load muscle`
    `muscle -in matches_vioA.m8.fasta -out matches_vioA.m8.afa -maxiters 2`

- Fasttree is not installed on the cluster so I downloaded using a conda environment (local or on cluster up to you) (more details on this later if necessary):
    `FastTree matches_vioA.m8.afa > tree_vioA_faa.afa`
   
- Viewed the trees using FigTree v 1.4.4



#SCRIPT STUFF
`for f in *\ *; do mv "$f" "${f// /_}"; done` 


- Make sure all genes are shown in all organisms 
- HMM then use Hammer.
- How many hits?  (vioc have two nodes? or seperate groups in the trees)
- Being non-purple does not mean no violacein 
- What are the other genomes HMM 
- What is the range, maybe only in proteobacteria? Beta and gamma. 


ALL: e20
Collimonas
Cfungivor_Ter331_GCF_000221045.1.fna.fasta
Cfungivor_Ter6_GCF_001584145.1.fna.fasta


3/5
Collimonas
Csp_OK307_GCF_900113705.1.fna.fasta
Csp_OK412_GCF_900112135.1.fna.fasta
Csp_OK607_GCF_900111995.1.fna.fasta



Chryseobacterium
Cindologenes_assembly.fna.fasta

Collimonas
Carenae_Cal35_GCF_000786695.1.fna.fasta
Carenae_Ter10_GCF_001584165.1.fna.fasta
Carenae_Ter282_GCF_001584205.1.fna.fasta

Collimonas
Cpratensis_Ter291_GCF_001584225.1.fna.fasta
Cpratensis_Ter91_GCF_001584185.1.fna.fasta


# 
for i in `cat hmm_seqs.txt`; do cat matches_VioA_0.000000000000000000001.m8.sfasta | grep -A1 $i; done 

module load muscle

muscle -in hmm_seqs_vioA_0.00000000000000000001.fasta -out hmm_seqs_vioA_0.00000000000000000001.msa
hmmbuild hmm_seqs_vioA_0.00000000000000000001.hmm hmm_seqs_vioA_0.00000000000000000001.msa
hmmsearch hmm_seqs_vioA_0.00000000000000000001.hmm matches_VioA_0.000000000000000000001.m8.sfasta


muscle -in hmm_seqs_vioD_0.00000000000000000001.fasta -out hmm_seqs_vioD_0.00000000000000000001.msa
hmmbuild hmm_seqs_vioD_0.00000000000000000001.hmm hmm_seqs_vioD_0.00000000000000000001.msa
hmmsearch hmm_seqs_vioD_0.00000000000000000001.hmm matches_VioD_0.1.m8.sfasta

muscle -in hmm_seqs_vioC_0.00000000000000000001.fasta -out hmm_seqs_vioC_0.00000000000000000001.msa
hmmbuild hmm_seqs_vioC_0.00000000000000000001.hmm hmm_seqs_vioC_0.00000000000000000001.msa
hmmsearch hmm_seqs_vioC_0.00000000000000000001.hmm matches_VioC_0.1.m8.sfasta


# Update names of the faa negative controls
        for f in *.faa ; 
            do awk '/>/{sub(">","&"FILENAME"__");sub(/\.faa/,x)}1' "$f" > "$f".fasta ; 
        done

# Create a FAA to query test e-value cutoff in the HMM

1. Create one fasta for all the pos controls assemblies:
    `for i in `cat ../notes/hmm_seqs.txt`; do cat $i.fna.fasta >> all_producers.fasta; done`
2. Create a protein database from the positive controls:
    `cat all_producers.fasta | prodigal > gene.all_producers.gbk -a protein.all_producers.faa`
3. Create one FAA for all the neg controls
    ``
    
4. Query all the FAA neg controls + all the positive controls and then check at which value each of the positive controls return a hit?
   `cat protein.all_producers.faa protein.non_producers.faa > protein.hmm_test.faa`
   `hmmsearch hmm_seqs_vioC_0.00000000000000000001.hmm protein.hmm_test.faa > hmmsearch_results_vioC.txt`
    
    `for i in `cat non_producers/list.txt`; do cat hmmsearch_results_vioC.txt | grep $i | head -1; done`
    `for i in `cat ../notes/hmm_seqs.txt`; do cat hmmsearch_results_vioC.txt | grep $i | head -1; done`
    
    
    Do all positives return a hit before any of the neg controls do?
    

