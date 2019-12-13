# Predicting Violacein producers

## Components needed 
- directory of all the genomes (.fna ) (/share/biocore/keith/eisenlab/new_all_seqs/)
- directory with the outgroups (in the .faa.fasta form) (/share/biocore/keith/eisenlab/notes/non_producers/)
- hmm_seqs (the known positives to create the model in a list) should be in directory with all genomes (/share/biocore/keith/eisenlab/notes/hmm_seqs.txt)
- directory with all of the Vio genes to use as a gold standard for getting the other known producer's genes for the hmm



![](./workflow.png)

First we need to do some general set up.
```
    # make a directory somewhere convenient.. mine is at /share/biocore/keith/eisenlab/ then cd there
    mkdir new_all_seqs
    
    # copy over the directory 
    cp  /share/eisenlab/gjospin/misc/Marina/Marina_test/genomes_Aug19/Vio_assemblies/Vio_genomes/* new_all_seqs/ 
    
    # make a directory that will store our hmms
    mkdir hmm_results/
    
    # make a directory to build and store trees in
    mkdir trees/
    
    # make a directory that will serve to store some extra information like non_producers and names of our known producers
    mkdir notes/
    cp /share/biocore/keith/eisenlab/notes/* notes/
     
    # this genome already in that directory is incorrect
    rm new_all_seqs/Dviolaceinigra_DSM15887.*
```

## Steps
1. Lets first add the file name to the beginning of all of the seqs in `new_all_seqs`, `notes/non_producers/`, and `notes/test_seqs`
    - Add the filename to the start of all the sequence for ease when building the tree later:
        ```
            for f in *.fna ; # .faa for the known non producers 
                do awk '/>/{sub(">","&"FILENAME"__");sub(/\.fna/,x)}1' "$f" > "$f".fasta ; 
            done
        ```
 
 2. Load the necessary packages:
    - `module load diamond`
    - `module load prodigal`
    - `module load hmmer`
    - `module list` (Currently Loaded Modulefiles):
       - -1) slurm/latest   2) anaconda3/4.5.12   3) snakemake/5.6.0   4) diamond/0.9.24   5) prodigal/2.6.3 

3. Create the database of the inferred proteins:
    ```
    cat *.fasta | prodigal > gene.coords2.gbk -a protein.translations.faa
    ```

4. Make the database of proteins:
    ```
    diamond makedb --in protein.translations.faa -d all_proteins
    ``` 

5. Run master_script.sh
    - TODO 
        - add parameters for directories of the script
        - Find a cutoff that can include the Dvio and Mvio so we can grab these. 
    ```
    ./master_script.sh
    ```

6. Check to make sure that all of the genes returned the number of hits (14) so 28 should be returned for all files since it includes the sequences as well. 
    ```
    for i in hmm_matches_Vio?_0.00000001.m8.sfasta; do cat $i | wc -l; done`
    ```
7. 
    ```
    for i in `cat ../notes/hmm_seqs.txt`; do cat matches_VioE_0.00000001.m8.sfasta | grep $i -A 1 | head -2; done > hmm_matches_VioE_0.00000001.m8.sfasta
    ```
    
8. Run build_and_test_hmms.sh, this will also test the query sequences against the non_producers and known producers 
    8. Run build_and_test_hmms.sh, this will also test the query sequences against the non_producers and known producers 
    ```
    ./../matches_output/build_and_test.hmms.sh
    ```
    - the parameters here can be changed based on the locations of the components mentioned at the start. 
    - then lets generate a "Probability" that the results in the HMM are a producer of Violacein 
    ```
    python prob_producer.py -i hmm_matches_VioA_0.00000001.m8.sfasta.results.final,hmm_matches_VioB_0.00000001.m8.sfasta.results.final,hmm_matches_VioC_0.00000001.m8.sfasta.results.final,hmm_matches_VioD_0.00000001.m8.sfasta.results.final,hmm_matches_VioE_0.00000001.m8.sfasta.results.final
    ```
    Small e_value means high chance of that gene... therefore 1-e-value for all genes to calc intersection of them all occuring
    ```
    Cviolaceum_12472_GCF_000007705.1 1.0
    jliv_NFR18_GCF_900119665.1 1.0
    JLiv_1522_LRHW01.1 1.0
    JLiv_MP5059B 1.0
    JLiv_HH100 1.0
    JLiv_HH102 1.0
    JLiv_HH103 1.0
    JLiv_HH104 1.0
    JLiv_HH106 1.0
    JLiv_HH107 1.0
    DugHH01_GCF_001758795.1_ASM175879v1_genomic 1.0
    DugHH105_GCF_001758545.1_ASM175854v1_genomic 1.0
    Mviolaceinigra_GCF_002752675.1_ASM275267v1_genomic 1.0
    PsuedoDugViolaceinigra_GCF_000425385.1_ASM42538v1_genomic 1.0
    D_thermophilium_GCF_000020965.1_ASM2096v1_protein 4.750000000000001e-06
    E_minutum_GCF_000020145.1_ASM2014v1_protein 6.246187500000001e-06
    ```
    - TODO for some reasone my query sequences are not showing up in hmm queries at all. 
    
9. Get the 16s sequences for building the tree. 
    ```
    diamond blastx -d master -q ../16s/e_coli_16s.fasta -o testing.txt --max-target-seqs 0 --evalue 0.1 -f sam
    diamond blastx -d master -q ../16s/e_coli_16s.fasta -o testing.txt --max-target-seqs 0 --evalue 0.1 -f sam
    tail -n +6 testing.txt | awk '{print ">"$3,$1"\n"$10}' > testing.txt.fasta
    for i in `cat ../notes/hmm_seqs.txt`; do cat testing.txt.fasta | grep $i -A 1 | head -2; done > hmm_16s_seqs.fasta
    ```    
10. Build trees for the 16s and genes. /share/biocore/keith/eisenlab/trees/

    ```
    module load anaconda3
    source activate ~/synapse
    for i in *.sfasta.msa; do FastTree $i > "${i}.tre"; done
    ```
    


