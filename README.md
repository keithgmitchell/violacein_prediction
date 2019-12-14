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

5. Run master_script.sh. This script will grab the matches for each of the genes of interest at any e-value of interest.  
    ```
    ./master_script.sh
    ```
    The output here will be some fasta files with the hits of interest.

7.  Grab the sequences of interest that we want to build the hmm from.
    ```
    for i in `cat ../notes/hmm_seqs.txt`; do cat matches_VioE_0.00000001.m8.sfasta | grep $i -A 1 | head -2; done > hmm_matches_VioE_0.00000001.m8.sfasta
    ```

6. Check to make sure that all of the genes returned the number of hits (14) so 28 should be returned for all files since it includes the sequences as well. 
    ```
    for i in hmm_matches_Vio?_0.00000001.m8.sfasta; do cat $i | wc -l; done`
    ```
    
8. Run build_and_test_hmms.sh, this will also test the query sequences against the non_producers and known producers 
    ```
    ./../matches_output/build_and_test.hmms.sh
    ```
    The output here will be a file like the following. (hmm_matches_VioA_0.00000001.m8.sfasta.results.final)
    
    ```
          0 2432.4  17.5          0 2432.2  17.5    1.0  1  Cviolaceum_12472_GCF_000007705.1__NC_005085.1_3182                                # 3561961 # 3564957 # -1 # ID=1_
          0 1450.7   7.8          0 1450.6   7.8    1.0  1  jliv_NFR18_GCF_900119665.1__NZ_FPKH01000001.1_1012                                # 1155693 # 1158710 # -1 # ID=11
          0 1456.6   7.0          0 1456.5   7.0    1.0  1  JLiv_1522_LRHW01.1__gi|1078206172|gb|LRHW01000061.1|_19                           # 16357 # 19377 # -1 # ID=82_19;
          0 1445.2   6.5          0 1445.0   6.5    1.0  1  JLiv_MP5059B__gi|1078198150|gb|LRHX01000012.1|_113                                # 130375 # 133395 # 1 # ID=158_1
          0 1434.5   7.3          0 1434.3   7.3    1.0  1  JLiv_HH100__gi|1078222769|gb|LRHY01000006.1|_364                                  # 380888 # 383902 # 1 # ID=2271_
          0 1434.5   7.3          0 1434.3   7.3    1.0  1  JLiv_HH102__gi|1078226882|gb|LRHZ01000004.1|_390                                  # 493670 # 496684 # -1 # ID=325_
          0 1434.5   7.3          0 1434.3   7.3    1.0  1  JLiv_HH103__gi|1078236093|gb|LRIA01000027.1|_184                                  # 206314 # 209328 # -1 # ID=473_
          0 1445.7   6.8          0 1445.5   6.8    1.0  1  JLiv_HH104__gi|1078235330|gb|LRIB01000027.1|_250                                  # 265396 # 268413 # 1 # ID=3755_
          0 1434.9   8.2          0 1434.8   8.2    1.0  1  JLiv_HH106__gi|1078243536|gb|LRIC01000002.1|_386                                  # 493596 # 496610 # -1 # ID=3795
          0 1437.5   8.1          0 1437.3   8.1    1.0  1  JLiv_HH107__gi|1078255902|gb|LRID01000023.1|_67                                   # 71830 # 74844 # 1 # ID=1795_67
          0 1441.9  12.0          0 1441.7  12.0    1.0  1  DugHH01_GCF_001758795.1_ASM175879v1_genomic__NZ_LRON01000001.1_481                # 631843 # 634863 # 1 # ID=2936_
          0 1437.7  12.1          0 1437.6  12.1    1.0  1  DugHH105_GCF_001758545.1_ASM175854v1_genomic__NZ_LRHV01000001.1_231               # 309508 # 312528 # -1 # ID=916_
          0 1462.4   9.4          0 1462.2   9.4    1.0  1  Mviolaceinigra_GCF_002752675.1_ASM275267v1_genomic__NZ_CP024608.1_4234            # 5021665 # 5024682 # -1 # ID=20
          0 1435.5  10.1          0 1435.3  10.1    1.0  1  PsuedoDugViolaceinigra_GCF_000425385.1_ASM42538v1_genomic__NZ_AUDI01000004.1_430  # 450837 # 453863 # 1 # ID=2056_
       0.24   14.4   0.2       0.39   13.7   0.2    1.2  1  D_thermophilium_GCF_000020965.1_ASM2096v1_protein__WP_012547829.1                 NAD(P)/FAD-dependent oxidoreduct
    0.00061   22.9   0.0        1.8   11.5   0.0    2.3  2  E_minutum_GCF_000020145.1_ASM2014v1_protein__WP_012414458.1                       oxidoreductase [Elusimicrobium m
    ```
    
    - the parameters here can be changed based on the locations of the components mentioned at the start. 
    - then lets generate a "Probability" that the results in the HMM are a producer of Violacein 
    ```
    python prob_producer.py -i hmm_matches_VioA_0.00000001.m8.sfasta.results.final,hmm_matches_VioB_0.00000001.m8.sfasta.results.final,hmm_matches_VioC_0.00000001.m8.sfasta.results.final,hmm_matches_VioD_0.00000001.m8.sfasta.results.final,hmm_matches_VioE_0.00000001.m8.sfasta.results.final
    ```
    Small e_value means high chance of that gene... therefore 1-e-value multiplied for all genes to calc intersection of them all occuring
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
    - TODO find some query sequences that return the result  
    
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
    source activate /home/keithgmitchell/synapse
    for i in *.sfasta.msa; do FastTree $i > "${i}.tre"; done
    ```
    


