

### Table of contents

[Import raw data to iridis](#import-raw-data-to-iridis)

[Pre processing](#pre-processing)
   - [Trimmomatic](#trimmomatic)
   - [Normalise reads ](#normalise-reads)

[Trininty assemblies](#trininty-assemblies)

[Annotate each assembly using a blast search](#annotate-each-assembly-using-a-blast-search)


### Import raw data to iridis

```
# in iridis
mkdir /scratch/oww1c19/argyranthemum_transcriptomics
mkdir /scratch/oww1c19/argyranthemum_transcriptomics/raw_data

# in local terminal
rsync -avz --partial /drives/d/Trancriptomics/C101HW16120181/raw_data oww1c19@iridis4_c.soton.ac.uk:/scratch/oww1c19/argyranthemum_transcriptomics
```

The samples are named A1 to A24, with forward reads denoted by "_1.fq.gz" ending and reverse reads with "_2.fq.gz" ending.

A table with the taxonomy of each sample is below.

A1  A. broussonetii
A2  A. broussonetii
A3  A. broussonetii
A4  A. broussonetii
A5  A. broussonetii
A6  A. broussonetii
A7  A. sundingii
A8  A. sundingii
A9  A. sundingii
A10 A. sundingii
A11 A. sundingii
A12 A. sundingii
A13 A. lemsiii
A14 A. lemsiii
A15 A. lemsiii
A16 A. lemsiii
A17 A. lemsiii
A18 A. lemsiii
A19 A. frutescens
A20 A. frutescens
A21 A. frutescens
A22 A. frutescens
A23 A. frutescens
A24 A. frutescens


### Pre processing

#### Trimmomatic
```
# set up dir
mkdir /scratch/oww1c19/argyranthemum_transcriptomics/trimmomatic
cd /scratch/oww1c19/argyranthemum_transcriptomics/trimmomatic

# cp script to dir
cp /scratch/oww1c19/argyranthemum_transcriptomics/scripts/script_trimmomatic.pbs .

# run trimmomatic
qsub script_trimmomatic.pbs
```

#### Normalise reads 
```
# set up dir
mkdir /scratch/oww1c19/argyranthemum_transcriptomics/normalise_reads
cd /scratch/oww1c19/argyranthemum_transcriptomics/normalise_reads

# create normalisation input lists for each species
for i in {1..6}
do
   echo /scratch/oww1c19/argyranthemum_transcriptomics/trimmomatic/trim_paired_A${i}_1.fq.gz >> input_bro_left
   echo /scratch/oww1c19/argyranthemum_transcriptomics/trimmomatic/trim_paired_A${i}_2.fq.gz >> input_bro_right
done

for i in {7..12}
do
   echo /scratch/oww1c19/argyranthemum_transcriptomics/trimmomatic/trim_paired_A${i}_1.fq.gz >> input_sun_left
   echo /scratch/oww1c19/argyranthemum_transcriptomics/trimmomatic/trim_paired_A${i}_2.fq.gz >> input_sun_right
done

for i in {13..18}
do
   echo /scratch/oww1c19/argyranthemum_transcriptomics/trimmomatic/trim_paired_A${i}_1.fq.gz >> input_lem_left
   echo /scratch/oww1c19/argyranthemum_transcriptomics/trimmomatic/trim_paired_A${i}_2.fq.gz >> input_lem_right
done

for i in {19..24}
do
   echo /scratch/oww1c19/argyranthemum_transcriptomics/trimmomatic/trim_paired_A${i}_1.fq.gz >> input_fru_left
   echo /scratch/oww1c19/argyranthemum_transcriptomics/trimmomatic/trim_paired_A${i}_2.fq.gz >> input_fru_right
done

# run normalisation for each species
# note each script produces large temporary files so run idividually and check memory
qsub script_normalise_bro.pbs
qsub script_normalise_sun.pbs
qsub script_normalise_lem.pbs
qsub script_normalise_fru.pbs

# get input lists for normalisation of all species together
ls -1 normalise_*/trim*1.fq.gz_ext_all_reads.normalized_K25_maxC30_minC0_maxCV10000.fq | sed 's:normalise:/scratch/oww1c19/argyranthemum_transcriptomics/normalise_reads/normalise:g' > input_all_left
ls -1 normalise_*/trim*2.fq.gz_ext_all_reads.normalized_K25_maxC30_minC0_maxCV10000.fq | sed 's:normalise:/scratch/oww1c19/argyranthemum_transcriptomics/normalise_reads/normalise:g' > input_all_right

# run normalisation for all species together
qsub script_normalise_all.pbs
```

### Trininty assemblies
```
# set up dir
mkdir /scratch/oww1c19/argyranthemum_transcriptomics/trinity_all
mkdir /scratch/oww1c19/argyranthemum_transcriptomics/trinity_separate

# species specific assemblies first
cd /scratch/oww1c19/argyranthemum_transcriptomics/trinity_separate

# cp scripts to dir
for i in bro sun lem fru; do
   cp /scratch/oww1c19/argyranthemum_transcriptomics/scripts/script_trinity_${i}.pbs .
done

# run each assembly
qsub script_trinity_bro.pbs
qsub script_trinity_fru.pbs
qsub script_trinity_lem.pbs
qsub script_trinity_sun.pbs

# copy each assembly and edit such that the read headers (rh) containg a sp ref
sed -e 's/TRINITY_/TRINITY_bro_/g' trinity_bro.Trinity.fasta > trinity_bro_rh.Trinity.fasta
sed -e 's/TRINITY_/TRINITY_fru_/g' trinity_fru.Trinity.fasta > trinity_fru_rh.Trinity.fasta
sed -e 's/TRINITY_/TRINITY_lem_/g' trinity_lem.Trinity.fasta > trinity_lem_rh.Trinity.fasta
sed -e 's/TRINITY_/TRINITY_sun_/g' trinity_sun.Trinity.fasta > trinity_sun_rh.Trinity.fasta

# cat species assemblies (required for pipelie 3)
cat trinity_bro_rh.Trinity.fasta trinity_fru_rh.Trinity.fasta trinity_lem_rh.Trinity.fasta trinity_sun_rh.Trinity.fasta > trinity_cat_rh.Trinity.fasta

# assemble all species together in a single assembly
cd /scratch/oww1c19/argyranthemum_transcriptomics/trinity_all

# cp script to dir
cp /scratch/oww1c19/argyranthemum_transcriptomics/scripts/script_trinity_all.pbs .

# run assembly
qsub script_trinity_all.pbs

# edit read headers as above
sed -e 's/TRINITY_/TRINITY_all_/g' trinity_all.Trinity.fasta > trinity_all_rh.Trinity.fasta
```


### Annotate each assembly using a blast search
```
# set up dir
mkdir /scratch/oww1c19/argyranthemum_transcriptomics/blast_separate
mkdir /scratch/oww1c19/argyranthemum_transcriptomics/blast_all

# blast species specific assemblies
cd /scratch/oww1c19/argyranthemum_transcriptomics/blast_separate

# cp scripts to dir
cp /scratch/oww1c19/argyranthemum_transcriptomics/scripts/script_blast_separate.pbs .
cp /scratch/oww1c19/argyranthemum_transcriptomics/scripts/top_hit.R .

# run blast
qsub script_blast_separate.pbs


# blast single species assembly
cd /scratch/oww1c19/argyranthemum_transcriptomics/blast_all

# cp scripts to dir
cp /scratch/oww1c19/argyranthemum_transcriptomics/scripts/script_blast_all.pbs .
cp /scratch/oww1c19/argyranthemum_transcriptomics/scripts/top_hit.R .

# run blast
qsub script_blast_all.pbs
```



### Implement pipelines 1-6
```
# set up dirs for each pipeline
for i in {1..6}; do  mkdir /scratch/oww1c19/argyranthemum_transcriptomics/pipeline_${i}; done
```

#### Pipeline 1

##### Pipeline 1 - Prep reference
```
# set up dir
mkdir /scratch/oww1c19/argyranthemum_transcriptomics/pipeline_1/prep_reference
cd /scratch/oww1c19/argyranthemum_transcriptomics/pipeline_1/prep_reference

# cp script to dir
cp /scratch/oww1c19/argyranthemum_transcriptomics/scripts/script_prep_ref.pbs .

# run prep reference
qsub -v REF_PATH=/scratch/oww1c19/argyranthemum_transcriptomics/trinity_all/,REF=trinity_all_rh.Trinity.fasta script_prep_ref.pbs
````blast_separate/ blast_separate/
##### Pipeline 1 - Align and estimate abundance
```
# set up dir
mkdir /scratch/oww1c19/argyranthemum_transcriptomics/pipeline_1/align_est
cd /scratch/oww1c19/argyranthemum_transcriptomics/pipeline_1/align_est

# cp script to dir
cp /scratch/oww1c19/argyranthemum_transcriptomics/scripts/script_align_est.pbs .

# run align and estimate abundance
qsub -t 1-24 -v REF=/scratch/oww1c19/argyranthemum_transcriptomics/pipeline_1/prep_reference/trinity_all_rh.Trinity.fasta script_align_est.pbs

# delete bam files to save space
rm align_est_*/*.bam
```
##### Pipeline 1 - Build expresison matrix
```
# set up dir
mkdir /scratch/oww1c19/argyranthemum_transcriptomics/pipeline_1/build_exp_matrix
cd /scratch/oww1c19/argyranthemum_transcriptomics/pipeline_1/build_exp_matrix

# cp script to dir
cp /scratch/oww1c19/argyranthemum_transcriptomics/scripts/script_build_exp_matrix.pbs .

# run build expression matrix
qsub script_build_exp_matrix.pbs
````
##### Pipeline 1 - Differential expression
```
# set up dir
mkdir /scratch/oww1c19/argyranthemum_transcriptomics/pipeline_1/diff_exp
cd /scratch/oww1c19/argyranthemum_transcriptomics/pipeline_1/diff_exp

# cp script to dir
cp /scratch/oww1c19/argyranthemum_transcriptomics/scripts/script_diff_exp.pbs .

# cp table of input samples required for de
cp /scratch/oww1c19/argyranthemum_transcriptomics/scripts/de_samples_file.txt .

# run build expression matrix
qsub script_diff_exp.pbs
```
##### Pipeline 1 - PtR
```
# set up dir
mkdir /scratch/oww1c19/argyranthemum_transcriptomics/pipeline_1/ptr
cd /scratch/oww1c19/argyranthemum_transcriptomics/pipeline_1/ptr

# cp script to dir
cp /scratch/oww1c19/argyranthemum_transcriptomics/scripts/script_ptr.pbs .

# cp table of input samples
cp /scratch/oww1c19/argyranthemum_transcriptomics/scripts/de_samples_file.txt .

# run ptr
qsub script_ptr.pbs
```
##### Pipeline 1 - Expression phenotypes and go terms
```
# set up dir
mkdir /scratch/oww1c19/argyranthemum_transcriptomics/pipeline_1/go_terms
cd /scratch/oww1c19/argyranthemum_transcriptomics/pipeline_1/go_terms

# cp scripts to dir
cp /scratch/oww1c19/argyranthemum_transcriptomics/scripts/script_go_terms.pbs .
cp /scratch/oww1c19/argyranthemum_transcriptomics/scripts/go_terms.R .

# run go terms
qsub -v REF=/scratch/oww1c19/argyranthemum_transcriptomics/pipeline_1/prep_reference/trinity_all_rh.Trinity.fasta,BLAST=/scratch/oww1c19/argyranthemum_transcriptomics/blast_all/out_blast_all_Athaliana_top_hits.txt  script_go_terms.pbs
```
##### Pipeline 1 - Reference all-by-all blast
```
# set up dir
mkdir /scratch/oww1c19/argyranthemum_transcriptomics/pipeline_1/ref_all_by_all_blast
cd /scratch/oww1c19/argyranthemum_transcriptomics/pipeline_1/ref_all_by_all_blast

# cp script to dir
cp /scratch/oww1c19/argyranthemum_transcriptomics/scripts/script_ref_all_by_all_blast.pbs .

# run all-by-allblast 
qsub -v REF=/scratch/oww1c19/argyranthemum_transcriptomics/pipeline_1/prep_reference/trinity_all_rh.Trinity.fasta script_ref_all_by_all_blast.pbs
```
##### Pipeline 1 - Chi-squared test
```
# set up dir
mkdir /scratch/oww1c19/argyranthemum_transcriptomics/pipeline_1/chi_squared_test
cd /scratch/oww1c19/argyranthemum_transcriptomics/pipeline_1/chi_squared_test

# cp scripts to dir
cp /scratch/oww1c19/argyranthemum_transcriptomics/scripts/script_chi_squared_test.pbs .
cp /scratch/oww1c19/argyranthemum_transcriptomics/scripts/chi_squared_test.R .

# run script to get chi-squared test data
qsub script_chi_squared_test.pbs
```

#### Pipeline 2

##### Pipeline 2 - cd-hit-est
```
# set up dir
mkdir /scratch/oww1c19/argyranthemum_transcriptomics/pipeline_2/cd_hit_est
cd /scratch/oww1c19/argyranthemum_transcriptomics/pipeline_2/cd_hit_est

# cp script to dir
cp /scratch/oww1c19/argyranthemum_transcriptomics/scripts/script_cd_hit_est_0.95.pbs .

# run cd-hit-est
qsub -v REF=/scratch/oww1c19/argyranthemum_transcriptomics/trinity_all/trinity_all_rh.Trinity.fasta,OUT=trinity_all_rh_cd_hit_0.95.Trinity.fasta script_cd_hit_est_0.95.pbs
```
##### Pipeline 2 - Prep reference
```
# set up dir
mkdir /scratch/oww1c19/argyranthemum_transcriptomics/pipeline_2/prep_reference
cd /scratch/oww1c19/argyranthemum_transcriptomics/pipeline_2/prep_reference

# cp script to dir
cp /scratch/oww1c19/argyranthemum_transcriptomics/scripts/script_prep_ref.pbs .

# run prep reference
qsub -v REF_PATH=/scratch/oww1c19/argyranthemum_transcriptomics/pipeline_2/cd_hit_est/,REF=trinity_all_rh_cd_hit_0.95.Trinity.fasta script_prep_ref.pbs
```
##### Pipeline 2 - Align and estimate abundance
```
# set up dir
mkdir /scratch/oww1c19/argyranthemum_transcriptomics/pipeline_2/align_est
cd /scratch/oww1c19/argyranthemum_transcriptomics/pipeline_2/align_est

# cp script to dir
cp /scratch/oww1c19/argyranthemum_transcriptomics/scripts/script_align_est.pbs .

# run align and estimate abundance
qsub -t 1-24 -v REF=/scratch/oww1c19/argyranthemum_transcriptomics/pipeline_2/prep_reference/trinity_all_rh_cd_hit_0.95.Trinity.fasta script_align_est.pbs

# delete bam files to save space
rm align_est_*/*.bam
```
##### Pipeline 2 - Build expresison matrix
```
# set up dir
mkdir /scratch/oww1c19/argyranthemum_transcriptomics/pipeline_2/build_exp_matrix
cd /scratch/oww1c19/argyranthemum_transcriptomics/pipeline_2/build_exp_matrix

# cp script to dir
cp /scratch/oww1c19/argyranthemum_transcriptomics/scripts/script_build_exp_matrix.pbs .

# run build expression matrix
qsub script_build_exp_matrix.pbs
````
##### Pipeline 2 - Differential expression
```
# set up dir
mkdir /scratch/oww1c19/argyranthemum_transcriptomics/pipeline_2/diff_exp
cd /scratch/oww1c19/argyranthemum_transcriptomics/pipeline_2/diff_exp

# cp script to dir
cp /scratch/oww1c19/argyranthemum_transcriptomics/scripts/script_diff_exp.pbs .

# cp table of input samples required for de
cp /scratch/oww1c19/argyranthemum_transcriptomics/scripts/de_samples_file.txt .

# run build expression matrix
qsub script_diff_exp.pbs
```
##### Pipeline 2 - PtR
```
# set up dir
mkdir /scratch/oww1c19/argyranthemum_transcriptomics/pipeline_2/ptr
cd /scratch/oww1c19/argyranthemum_transcriptomics/pipeline_2/ptr

# cp script to dir
cp /scratch/oww1c19/argyranthemum_transcriptomics/scripts/script_ptr.pbs .

# cp table of input samples
cp /scratch/oww1c19/argyranthemum_transcriptomics/scripts/de_samples_file.txt .

# run ptr
qsub script_ptr.pbs
```
##### Pipeline 2 - Expression phenotypes and go terms
```
# set up dir
mkdir /scratch/oww1c19/argyranthemum_transcriptomics/pipeline_2/go_terms
cd /scratch/oww1c19/argyranthemum_transcriptomics/pipeline_2/go_terms

# cp scripts to dir
cp /scratch/oww1c19/argyranthemum_transcriptomics/scripts/script_go_terms.pbs .
cp /scratch/oww1c19/argyranthemum_transcriptomics/scripts/go_terms.R .

# run go terms
qsub -v REF=/scratch/oww1c19/argyranthemum_transcriptomics/pipeline_2/prep_reference/trinity_all_rh_cd_hit_0.95.Trinity.fasta,BLAST=/scratch/oww1c19/argyranthemum_transcriptomics/blast_all/out_blast_all_Athaliana_top_hits.txt  script_go_terms.pbs
```
##### Pipeline 2 - Reference all-by-all blast
```
# set up dir
mkdir /scratch/oww1c19/argyranthemum_transcriptomics/pipeline_2/ref_all_by_all_blast
cd /scratch/oww1c19/argyranthemum_transcriptomics/pipeline_2/ref_all_by_all_blast

# cp script to dir
cp /scratch/oww1c19/argyranthemum_transcriptomics/scripts/script_ref_all_by_all_blast.pbs .

# run all-by-allblast 
qsub -v REF=/scratch/oww1c19/argyranthemum_transcriptomics/pipeline_2/prep_reference/trinity_all_rh_cd_hit_0.95.Trinity.fasta script_ref_all_by_all_blast.pbs
```
##### Pipeline 2 - Chi-squared test
```
# set up dir
mkdir /scratch/oww1c19/argyranthemum_transcriptomics/pipeline_2/chi_squared_test
cd /scratch/oww1c19/argyranthemum_transcriptomics/pipeline_2/chi_squared_test

# cp scripts to dir
cp /scratch/oww1c19/argyranthemum_transcriptomics/scripts/script_chi_squared_test.pbs .
cp /scratch/oww1c19/argyranthemum_transcriptomics/scripts/chi_squared_test.R .

# run script to get chi-squared test data
qsub script_chi_squared_test.pbs
```


#### Pipeline 3

##### Pipeline 3 - cd_hit_est
```
# set up dir
mkdir /scratch/oww1c19/argyranthemum_transcriptomics/pipeline_3/cd_hit_est
cd /scratch/oww1c19/argyranthemum_transcriptomics/pipeline_3/cd_hit_est

# cp script to dir
cp /scratch/oww1c19/argyranthemum_transcriptomics/scripts/script_cd_hit_est_0.95.pbs .

# run cd-hit-est
qsub -v REF=/scratch/oww1c19/argyranthemum_transcriptomics/trinity_separate/trinity_cat_rh.Trinity.fasta,OUT=trinity_cat_rh_cd_hit_0.95.Trinity.fasta script_cd_hit_est_0.95.pbs
```
##### Pipeline 3 - Prep reference
```
# set up dir
mkdir /scratch/oww1c19/argyranthemum_transcriptomics/pipeline_3/prep_reference
cd /scratch/oww1c19/argyranthemum_transcriptomics/pipeline_3/prep_reference

# cp script to dir
cp /scratch/oww1c19/argyranthemum_transcriptomics/scripts/script_prep_ref.pbs .

# run prep reference
qsub -v REF_PATH=/scratch/oww1c19/argyranthemum_transcriptomics/pipeline_3/cd_hit_est/,REF=trinity_cat_rh_cd_hit_0.95.Trinity.fasta script_prep_ref.pbs
```
##### Pipeline 3 - Align and estimate abundance
```
# set up dir
mkdir /scratch/oww1c19/argyranthemum_transcriptomics/pipeline_3/align_est
cd /scratch/oww1c19/argyranthemum_transcriptomics/pipeline_3/align_est

# cp script to dir
cp /scratch/oww1c19/argyranthemum_transcriptomics/scripts/script_align_est.pbs .

# run align and estimate abundance
qsub -t 1-24 -v REF=/scratch/oww1c19/argyranthemum_transcriptomics/pipeline_3/prep_reference/trinity_cat_rh_cd_hit_0.95.Trinity.fasta script_align_est.pbs

# delete bam files to save space
rm align_est_*/*.bam
```
##### Pipeline 3 - Build expresison matrix
```
# set up dir
mkdir /scratch/oww1c19/argyranthemum_transcriptomics/pipeline_3/build_exp_matrix
cd /scratch/oww1c19/argyranthemum_transcriptomics/pipeline_3/build_exp_matrix

# cp script to dir
cp /scratch/oww1c19/argyranthemum_transcriptomics/scripts/script_build_exp_matrix.pbs .

# run build expression matrix
qsub script_build_exp_matrix.pbs
````
##### Pipeline 3 - Differential expression
```
# set up dir
mkdir /scratch/oww1c19/argyranthemum_transcriptomics/pipeline_3/diff_exp
cd /scratch/oww1c19/argyranthemum_transcriptomics/pipeline_3/diff_exp

# cp script to dir
cp /scratch/oww1c19/argyranthemum_transcriptomics/scripts/script_diff_exp.pbs .

# cp table of input samples required for de
cp /scratch/oww1c19/argyranthemum_transcriptomics/scripts/de_samples_file.txt .

# run build expression matrix
qsub script_diff_exp.pbs
```
##### Pipeline 3 - PtR
```
# set up dir
mkdir /scratch/oww1c19/argyranthemum_transcriptomics/pipeline_3/ptr
cd /scratch/oww1c19/argyranthemum_transcriptomics/pipeline_3/ptr

# cp script to dir
cp /scratch/oww1c19/argyranthemum_transcriptomics/scripts/script_ptr.pbs .

# cp table of input samples
cp /scratch/oww1c19/argyranthemum_transcriptomics/scripts/de_samples_file.txt .

# run ptr
qsub script_ptr.pbs
```
##### Pipeline 3 - Expression phenotypes and go terms
```
# set up dir
mkdir /scratch/oww1c19/argyranthemum_transcriptomics/pipeline_3/go_terms
cd /scratch/oww1c19/argyranthemum_transcriptomics/pipeline_3/go_terms

# cp scripts to dir
cp /scratch/oww1c19/argyranthemum_transcriptomics/scripts/script_go_terms.pbs .
cp /scratch/oww1c19/argyranthemum_transcriptomics/scripts/go_terms.R .

# run go terms
# note the blast input for each indivdual assembly i.e. "out_blast_*_Athaliana_top_hits.txt"
qsub -v REF=/scratch/oww1c19/argyranthemum_transcriptomics/pipeline_3/prep_reference/trinity_cat_rh_cd_hit_0.95.Trinity.fasta,BLAST=/scratch/oww1c19/argyranthemum_transcriptomics/blast_separate/out_blast_*_Athaliana_top_hits.txt  script_go_terms.pbs
```
##### Pipeline 3 - Reference all-by-all blast
```
# set up dir
mkdir /scratch/oww1c19/argyranthemum_transcriptomics/pipeline_3/ref_all_by_all_blast
cd /scratch/oww1c19/argyranthemum_transcriptomics/pipeline_3/ref_all_by_all_blast

# cp script to dir
cp /scratch/oww1c19/argyranthemum_transcriptomics/scripts/script_ref_all_by_all_blast.pbs .

# run all-by-allblast 
qsub -v REF=/scratch/oww1c19/argyranthemum_transcriptomics/pipeline_3/prep_reference/trinity_cat_rh_cd_hit_0.95.Trinity.fasta script_ref_all_by_all_blast.pbs
```
##### Pipeline 3 - Chi-squared test
```
# set up dir
mkdir /scratch/oww1c19/argyranthemum_transcriptomics/pipeline_3/chi_squared_test
cd /scratch/oww1c19/argyranthemum_transcriptomics/pipeline_3/chi_squared_test

# cp scripts to dir
cp /scratch/oww1c19/argyranthemum_transcriptomics/scripts/script_chi_squared_test.pbs .
cp /scratch/oww1c19/argyranthemum_transcriptomics/scripts/chi_squared_test.R .

# run script to get chi-squared test data
qsub script_chi_squared_test.pbs
```


#### Preparation for pipelines 4, 5 and 6

##### Transdecoder
```
# set up dir
mkdir /scratch/oww1c19/argyranthemum_transcriptomics/transdecoder
cd /scratch/oww1c19/argyranthemum_transcriptomics/transdecoder

# cp script to dir
cp /scratch/oww1c19/argyranthemum_transcriptomics/scripts/script_transdecoder.pbs .

# run transdecoder
qsub script_transdecoder.pbs
```

##### cd-hit
```
# set up dir
mkdir /scratch/oww1c19/argyranthemum_transcriptomics/cd_hit
cd /scratch/oww1c19/argyranthemum_transcriptomics/cd_hit

# cp across transdecoder output files excluding the longer header name
awk ' { print $1 } ' /scratch/oww1c19/argyranthemum_transcriptomics/transdecoder/trinity_bro_rh.Trinity.fasta.transdecoder.pep  > trinity_pep_bro.fasta
awk ' { print $1 } ' /scratch/oww1c19/argyranthemum_transcriptomics/transdecoder/trinity_fru_rh.Trinity.fasta.transdecoder.pep  > trinity_pep_fru.fasta
awk ' { print $1 } ' /scratch/oww1c19/argyranthemum_transcriptomics/transdecoder/trinity_sun_rh.Trinity.fasta.transdecoder.pep  > trinity_pep_sun.fasta
awk ' { print $1 } ' /scratch/oww1c19/argyranthemum_transcriptomics/transdecoder/trinity_lem_rh.Trinity.fasta.transdecoder.pep  > trinity_pep_lem.fasta

# cp script to dir
cp /scratch/oww1c19/argyranthemum_transcriptomics/scripts/script_cd_hit_0.995.pbs .

# run cd hit
qsub script_cd_hit_0.995.pbs
```

##### Concatenate the species specific assemblies

The reference used for pipelines 4,5 and 6 is the species-specific cd-hit files concatenated. However we need nucleotide sequences for Trinity.
 
```
# set up dir
mkdir /scratch/oww1c19/argyranthemum_transcriptomics/concat_ref
cd /scratch/oww1c19/argyranthemum_transcriptomics/concat_ref

# get Trinity transcript ids from cd-hit assemblies 
grep -e "^>" -h /scratch/oww1c19/argyranthemum_transcriptomics/cd_hit/*cd_hit.fasta  | sed 's/>//g' > transcripts_cd_hit.txt

# concate transdecoder cds files remvoving additional info in header lines added by transdecoder
cat /scratch/oww1c19/argyranthemum_transcriptomics/transdecoder/*cds | awk ' { print $1 } ' > concat_transdecoder_cds.fasta

# use Trinity transcript ids from cd-hit assemblies ids to pull out transdecoder cds sequences
perl -ne 'if(/^>(\S+)/){$c=$i{$1}}$c?print:chomp;$i{$_}=1 if @ARGV' transcripts_cd_hit.txt concat_transdecoder_cds.fasta > concat_cd_hit_cds.fasta
```

##### Blast concatenated assembly
```
# set up dir
mkdir /scratch/oww1c19/argyranthemum_transcriptomics/blast_concat
cd /scratch/oww1c19/argyranthemum_transcriptomics/blast_concat

# cp scripts to dir
cp /scratch/oww1c19/argyranthemum_transcriptomics/scripts/script_blast_concat.pbs .
cp /scratch/oww1c19/argyranthemum_transcriptomics/scripts/top_hit.R .

# run blast
qsub script_blast_concat.pbs
```







#### Pipeline 4

##### Pipeline 4 - OrthoFinder
```
# set up dir
mkdir /scratch/oww1c19/argyranthemum_transcriptomics/pipeline_4/orthofinder
mkdir /scratch/oww1c19/argyranthemum_transcriptomics/pipeline_4/orthofinder/input_fasta
cd /scratch/oww1c19/argyranthemum_transcriptomics/pipeline_4/orthofinder

# cp across argyranthemum fasta files
cp /scratch/oww1c19/argyranthemum_transcriptomics/cd_hit/*cd_hit.fasta /scratch/oww1c19/argyranthemum_transcriptomics/pipeline_4/orthofinder/input_fasta/

# cp script to dir
cp /scratch/oww1c19/argyranthemum_transcriptomics/scripts/script_orthofinder_pipeline_4.pbs .

# run orthofinder
qsub script_orthofinder_pipeline_4.pbs
```
##### Pipeline 4 - Write gene to trans map
```
# set up dir
mkdir /scratch/oww1c19/argyranthemum_transcriptomics/pipeline_4/gene_trans_map
cd /scratch/oww1c19/argyranthemum_transcriptomics/pipeline_4/gene_trans_map

# cp scripts to dir
cp /scratch/oww1c19/argyranthemum_transcriptomics/scripts/write_gene_trans_map_pipeline_4.* .

# write gene_trans_map
bash write_gene_trans_map_pipeline_4.sh
```
##### Pipeline 4 - Prep reference
```
# set up dir
mkdir /scratch/oww1c19/argyranthemum_transcriptomics/pipeline_4/prep_reference
cd /scratch/oww1c19/argyranthemum_transcriptomics/pipeline_4/prep_reference

# cp script to dir
cp /scratch/oww1c19/argyranthemum_transcriptomics/scripts/script_prep_ref_gene_trans_map.pbs .

# run prep reference
qsub -v REF_PATH=/scratch/oww1c19/argyranthemum_transcriptomics/concat_ref,REF=concat_cd_hit_cds.fasta script_prep_ref_gene_trans_map.pbs
```
##### Pipeline 4 - Align and estimate abundance
```
# set up dir
mkdir /scratch/oww1c19/argyranthemum_transcriptomics/pipeline_4/align_est
cd /scratch/oww1c19/argyranthemum_transcriptomics/pipeline_4/align_est

# cp script to dir
cp /scratch/oww1c19/argyranthemum_transcriptomics/scripts/script_align_est_gene_trans_map.pbs .

# run align and estimate abundance
# note I had to run this in smaller batches to avoid exceeding memory limit
qsub -t 1-24 -v REF=/scratch/oww1c19/argyranthemum_transcriptomics/pipeline_4/prep_reference/concat_cd_hit_cds.fasta script_align_est_gene_trans_map.pbs

# delete bam files to save space
rm align_est_*/*.bam
```
##### Pipeline 4 - Build expresison matrix
```
# set up dir
mkdir /scratch/oww1c19/argyranthemum_transcriptomics/pipeline_4/build_exp_matrix
cd /scratch/oww1c19/argyranthemum_transcriptomics/pipeline_4/build_exp_matrix

# cp script to dir
cp /scratch/oww1c19/argyranthemum_transcriptomics/scripts/script_build_exp_matrix_gene_trans_map.pbs .

# run build expression matrix
qsub script_build_exp_matrix_gene_trans_map.pbs
````
##### Pipeline 4 - Differential expression
```
# set up dir
mkdir /scratch/oww1c19/argyranthemum_transcriptomics/pipeline_4/diff_exp
cd /scratch/oww1c19/argyranthemum_transcriptomics/pipeline_4/diff_exp

# cp script to dir
cp /scratch/oww1c19/argyranthemum_transcriptomics/scripts/script_diff_exp_gene_trans_map.pbs .

# cp table of input samples required for de
cp /scratch/oww1c19/argyranthemum_transcriptomics/scripts/de_samples_file.txt .

# run build expression matrix
qsub script_diff_exp_gene_trans_map.pbs
```
##### Pipeline 4 - PtR
```
# set up dir
mkdir /scratch/oww1c19/argyranthemum_transcriptomics/pipeline_4/ptr
cd /scratch/oww1c19/argyranthemum_transcriptomics/pipeline_4/ptr

# cp script to dir
cp /scratch/oww1c19/argyranthemum_transcriptomics/scripts/script_ptr_pipelines_4_5.pbs .

# cp table of input samples
cp /scratch/oww1c19/argyranthemum_transcriptomics/scripts/de_samples_file.txt .

# run ptr
qsub script_ptr.pbs
```
##### Pipeline 4 - Expression phenotypes and go terms
```
# set up dir
mkdir /scratch/oww1c19/argyranthemum_transcriptomics/pipeline_4/go_terms
cd /scratch/oww1c19/argyranthemum_transcriptomics/pipeline_4/go_terms

# cp scripts to dir
cp /scratch/oww1c19/argyranthemum_transcriptomics/scripts/script_go_terms_gene_trans_map.pbs .
cp /scratch/oww1c19/argyranthemum_transcriptomics/scripts/go_terms_gene_trans_map.R .

# run go terms
# note the blast input for each indivdual assembly i.e. "out_blast_*_Athaliana_top_hits.txt"
qsub -v BLAST=/scratch/oww1c19/argyranthemum_transcriptomics/blast_concat/out_blast_concat_Athaliana_top_hits.txt  script_go_terms_gene_trans_map.pbs
```
##### Pipeline 4 - Reference all-by-all blast
```
# set up dir
mkdir /scratch/oww1c19/argyranthemum_transcriptomics/pipeline_4/ref_all_by_all_blast
cd /scratch/oww1c19/argyranthemum_transcriptomics/pipeline_4/ref_all_by_all_blast

# cp scripts to dir
cp /scratch/oww1c19/argyranthemum_transcriptomics/scripts/script_ref_all_by_all_blast_gene_trans_map.pbs .
cp /scratch/oww1c19/argyranthemum_transcriptomics/scripts/identify_reference_transcript.R .

# run all-by-all blast 
qsub script_ref_all_by_all_blast_gene_trans_map.pbs
```
##### Pipeline 4 - Chi-squared test
```
# set up dir
mkdir /scratch/oww1c19/argyranthemum_transcriptomics/pipeline_4/chi_squared_test
cd /scratch/oww1c19/argyranthemum_transcriptomics/pipeline_4/chi_squared_test

# cp scripts to dir
cp /scratch/oww1c19/argyranthemum_transcriptomics/scripts/script_chi_squared_test_gene_trans_map.pbs .
cp /scratch/oww1c19/argyranthemum_transcriptomics/scripts/chi_squared_test_gene_trans_map.R .

# run script to get chi-squared test data
qsub script_chi_squared_test_gene_trans_map.pbs
```



### pipeline 5

##### Pipeline 5 - OrthoFinder
```
# set up dir
mkdir /scratch/oww1c19/argyranthemum_transcriptomics/pipeline_5/orthofinder
mkdir /scratch/oww1c19/argyranthemum_transcriptomics/pipeline_5/orthofinder/input_fasta
cd /scratch/oww1c19/argyranthemum_transcriptomics/pipeline_5/orthofinder

# cp across argyranthemum fasta files
cp /scratch/oww1c19/argyranthemum_transcriptomics/cd_hit/*cd_hit.fasta /scratch/oww1c19/argyranthemum_transcriptomics/pipeline_4/orthofinder/input_fasta/

# cp across reference fasta files
cp /scratch/oww1c19/argyranthemum_transcriptomics/phytozome_protein_downloads/* /scratch/oww1c19/argyranthemum_transcriptomics/pipeline_5/orthofinder/input_fasta/


# cp scripts to dir
cp /scratch/oww1c19/argyranthemum_transcriptomics/scripts/script_orthofinder_pipeline_5_* .

# run orthofinder
# needed to split this up on iridis
qsub script_orthofinder_pipeline_5_1.pbs
qsub script_orthofinder_pipeline_5_2.pbs
```
##### Pipeline 5 - check orthogroup monophyly
```
# set up dir
mkdir /scratch/oww1c19/argyranthemum_transcriptomics/pipeline_5/check_monophyly
cd /scratch/oww1c19/argyranthemum_transcriptomics/pipeline_5/check_monophyly

# cp scripts to dir
cp /scratch/oww1c19/argyranthemum_transcriptomics/scripts/*check_monophyly_pipeline_5* .

# run script
qsub script_check_monophyly_pipeline_5.pbs
```
##### Pipeline 5 - Write gene to trans map
```
# set up dir
mkdir /scratch/oww1c19/argyranthemum_transcriptomics/pipeline_5/gene_trans_map
cd /scratch/oww1c19/argyranthemum_transcriptomics/pipeline_5/gene_trans_map

# cp scripts to dir
cp /scratch/oww1c19/argyranthemum_transcriptomics/scripts/*write_gene_trans_map_pipeline_5* .

# write gene_trans_map
qsub script_write_gene_trans_map_pipeline_5.pbs
```
##### Pipeline 5 - Prep reference
```
# set up dir
mkdir /scratch/oww1c19/argyranthemum_transcriptomics/pipeline_5/prep_reference
cd /scratch/oww1c19/argyranthemum_transcriptomics/pipeline_5/prep_reference

# cp script to dir
cp /scratch/oww1c19/argyranthemum_transcriptomics/scripts/script_prep_ref_gene_trans_map.pbs .

# run prep reference
qsub -v REF_PATH=/scratch/oww1c19/argyranthemum_transcriptomics/concat_ref,REF=concat_cd_hit_cds.fasta script_prep_ref_gene_trans_map.pbs
```
##### Pipeline 5 - Align and estimate abundance
```
# set up dir
mkdir /scratch/oww1c19/argyranthemum_transcriptomics/pipeline_5/align_est
cd /scratch/oww1c19/argyranthemum_transcriptomics/pipeline_5/align_est

# cp script to dir
cp /scratch/oww1c19/argyranthemum_transcriptomics/scripts/script_align_est_gene_trans_map.pbs .

# run align and estimate abundance
# note I had to run this in smaller batches to avoid exceeding memory limit
qsub -t 1-24 -v REF=/scratch/oww1c19/argyranthemum_transcriptomics/pipeline_5/prep_reference/concat_cd_hit_cds.fasta script_align_est_gene_trans_map.pbs

# delete bam files to save space
rm align_est_*/*.bam
```
##### Pipeline 5 - Build expresison matrix
```
# set up dir
mkdir /scratch/oww1c19/argyranthemum_transcriptomics/pipeline_5/build_exp_matrix
cd /scratch/oww1c19/argyranthemum_transcriptomics/pipeline_5/build_exp_matrix

# cp script to dir
cp /scratch/oww1c19/argyranthemum_transcriptomics/scripts/script_build_exp_matrix_gene_trans_map.pbs .

# run build expression matrix
qsub script_build_exp_matrix_gene_trans_map.pbs
````
##### Pipeline 5 - Differential expression
```
# set up dir
mkdir /scratch/oww1c19/argyranthemum_transcriptomics/pipeline_5/diff_exp
cd /scratch/oww1c19/argyranthemum_transcriptomics/pipeline_5/diff_exp

# cp script to dir
cp /scratch/oww1c19/argyranthemum_transcriptomics/scripts/script_diff_exp_gene_trans_map.pbs .

# cp table of input samples required for de
cp /scratch/oww1c19/argyranthemum_transcriptomics/scripts/de_samples_file.txt .

# run build expression matrix
qsub script_diff_exp_gene_trans_map.pbs
```
##### Pipeline 5 - PtR
```
# set up dir
mkdir /scratch/oww1c19/argyranthemum_transcriptomics/pipeline_5/ptr
cd /scratch/oww1c19/argyranthemum_transcriptomics/pipeline_5/ptr

# cp script to dir
cp /scratch/oww1c19/argyranthemum_transcriptomics/scripts/script_ptr_pipelines_4_5.pbs .

# cp table of input samples
cp /scratch/oww1c19/argyranthemum_transcriptomics/scripts/de_samples_file.txt .

# run ptr
qsub script_ptr.pbs
```
##### Pipeline 5 - Expression phenotypes and go terms
```
# set up dir
mkdir /scratch/oww1c19/argyranthemum_transcriptomics/pipeline_5/go_terms
cd /scratch/oww1c19/argyranthemum_transcriptomics/pipeline_5/go_terms

# cp scripts to dir
cp /scratch/oww1c19/argyranthemum_transcriptomics/scripts/script_go_terms_gene_trans_map.pbs .
cp /scratch/oww1c19/argyranthemum_transcriptomics/scripts/go_terms_gene_trans_map.R .

# run go terms
# note the blast input for each indivdual assembly i.e. "out_blast_*_Athaliana_top_hits.txt"
qsub -v BLAST=/scratch/oww1c19/argyranthemum_transcriptomics/blast_concat/out_blast_concat_Athaliana_top_hits.txt  script_go_terms_gene_trans_map.pbs
```
##### Pipeline 5 - Reference all-by-all blast
```
# set up dir
mkdir /scratch/oww1c19/argyranthemum_transcriptomics/pipeline_5/ref_all_by_all_blast
cd /scratch/oww1c19/argyranthemum_transcriptomics/pipeline_5/ref_all_by_all_blast

# cp scripts to dir
cp ../../scripts/script_ref_all_by_all_blast_gene_trans_map.pbs .
cp ../../scripts/identify_reference_transcript.R .

# run all-by-all blast 
qsub script_ref_all_by_all_blast_gene_trans_map.pbs
```
##### Pipeline 5 - Chi-squared test
```
# set up dir
mkdir /scratch/oww1c19/argyranthemum_transcriptomics/pipeline_5/chi_squared_test
cd /scratch/oww1c19/argyranthemum_transcriptomics/pipeline_5/chi_squared_test

# cp scripts to dir
cp /scratch/oww1c19/argyranthemum_transcriptomics/scripts/script_chi_squared_test_gene_trans_map.pbs .
cp /scratch/oww1c19/argyranthemum_transcriptomics/scripts/chi_squared_test_gene_trans_map.R .

# run script to get chi-squared test data
qsub script_chi_squared_test_gene_trans_map.pbs
```



### pipeline 6

##### Filter orthogroups with more than one hit in Helinathus from gene trans map
```
# set up dir
mkdir /scratch/oww1c19/argyranthemum_transcriptomics/pipeline_6/orthogroups_blast
cd /scratch/oww1c19/argyranthemum_transcriptomics/pipeline_6/orthogroups_blast

# cp balst results from concat assembly
# we only need columns for qseqid, sseqid and evalue
cut -f 1,2,11 /scratch/oww1c19/argyranthemum_transcriptomics/blast_concat/out_blast_concat_Hannuus_top_hits.txt > blast.txt

# gene trans map used for pipline 5
cp /scratch/oww1c19/argyranthemum_transcriptomics/pipeline_5/gene_trans_map/gene_trans_map.txt .

# run Rscript
module add R/3.4.2
Rscript orthogroups_single_blast_hit.R

# [1] "37256 orthogroups in total"
# [1] "32333 (86.79%) orthogroups with single hit"
# [1] "4923 (13.21%) orthogroups with multiple hits"

# cp gene_trans_map_single_hit.txt to gene_trans_map/gene_trans_map.txt to work with scripts
mkdir /scratch/oww1c19/argyranthemum_transcriptomics/pipeline_6/gene_trans_map
cp /scratch/oww1c19/argyranthemum_transcriptomics/pipeline_6/orthogroups_blast/gene_trans_map_single_hit.txt /scratch/oww1c19/argyranthemum_transcriptomics/pipeline_6/gene_trans_map/gene_trans_map.txt
```
##### Pipeline 6 - Prep reference
```
# set up dir
mkdir /scratch/oww1c19/argyranthemum_transcriptomics/pipeline_6/prep_reference
cd /scratch/oww1c19/argyranthemum_transcriptomics/pipeline_6/prep_reference

# cp script to dir
cp /scratch/oww1c19/argyranthemum_transcriptomics/scripts/script_prep_ref_gene_trans_map.pbs .

# run prep reference
qsub -v REF_PATH=/scratch/oww1c19/argyranthemum_transcriptomics/concat_ref,REF=concat_cd_hit_cds.fasta script_prep_ref_gene_trans_map.pbs
```
##### Pipeline 6 - Align and estimate abundance
```
# set up dir
mkdir /scratch/oww1c19/argyranthemum_transcriptomics/pipeline_6/align_est
cd /scratch/oww1c19/argyranthemum_transcriptomics/pipeline_6/align_est

# cp script to dir
cp /scratch/oww1c19/argyranthemum_transcriptomics/scripts/script_align_est_gene_trans_map.pbs .

# run align and estimate abundance
# note I had to run this in smaller batches to avoid exceeding memory limit
qsub -t 1-24 -v REF=/scratch/oww1c19/argyranthemum_transcriptomics/pipeline_6/prep_reference/concat_cd_hit_cds.fasta script_align_est_gene_trans_map.pbs

# delete bam files to save space
rm align_est_*/*.bam
```
##### Pipeline 6 - Build expresison matrix
```
# set up dir
mkdir /scratch/oww1c19/argyranthemum_transcriptomics/pipeline_6/build_exp_matrix
cd /scratch/oww1c19/argyranthemum_transcriptomics/pipeline_6/build_exp_matrix

# cp script to dir
cp /scratch/oww1c19/argyranthemum_transcriptomics/scripts/script_build_exp_matrix_gene_trans_map.pbs .

# run build expression matrix
qsub script_build_exp_matrix_gene_trans_map.pbs
````
##### Pipeline 6 - Differential expression
```
# set up dir
mkdir /scratch/oww1c19/argyranthemum_transcriptomics/pipeline_6/diff_exp
cd /scratch/oww1c19/argyranthemum_transcriptomics/pipeline_6/diff_exp

# cp script to dir
cp /scratch/oww1c19/argyranthemum_transcriptomics/scripts/script_diff_exp_gene_trans_map.pbs .

# cp table of input samples required for de
cp /scratch/oww1c19/argyranthemum_transcriptomics/scripts/de_samples_file.txt .

# run build expression matrix
qsub script_diff_exp_gene_trans_map.pbs
```
##### Pipeline 6 - PtR
```
# set up dir
mkdir /scratch/oww1c19/argyranthemum_transcriptomics/pipeline_6/ptr
cd /scratch/oww1c19/argyranthemum_transcriptomics/pipeline_6/ptr

# cp script to dir
cp /scratch/oww1c19/argyranthemum_transcriptomics/scripts/script_ptr_pipelines_4_5.pbs .

# cp table of input samples
cp /scratch/oww1c19/argyranthemum_transcriptomics/scripts/de_samples_file.txt .

# run ptr
qsub script_ptr.pbs
```
##### Pipeline 6 - Expression phenotypes and go terms
```
# set up dir
mkdir /scratch/oww1c19/argyranthemum_transcriptomics/pipeline_6/go_terms
cd /scratch/oww1c19/argyranthemum_transcriptomics/pipeline_6/go_terms

# cp scripts to dir
cp /scratch/oww1c19/argyranthemum_transcriptomics/scripts/script_go_terms_gene_trans_map.pbs .
cp /scratch/oww1c19/argyranthemum_transcriptomics/scripts/go_terms_gene_trans_map.R .

# run go terms
# note the blast input for each indivdual assembly i.e. "out_blast_*_Athaliana_top_hits.txt"
qsub -v BLAST=/scratch/oww1c19/argyranthemum_transcriptomics/blast_concat/out_blast_concat_Athaliana_top_hits.txt  script_go_terms_gene_trans_map.pbs
```
##### Pipeline 6 - Reference all-by-all blast
```
# set up dir
mkdir /scratch/oww1c19/argyranthemum_transcriptomics/pipeline_6/ref_all_by_all_blast
cd /scratch/oww1c19/argyranthemum_transcriptomics/pipeline_6/ref_all_by_all_blast

# cp scripts to dir
cp ../../scripts/script_ref_all_by_all_blast_gene_trans_map.pbs .
cp ../../scripts/identify_reference_transcript.R .

# run all-by-all blast 
qsub script_ref_all_by_all_blast_gene_trans_map.pbs
```
##### Pipeline 5 - Chi-squared test
```
# set up dir
mkdir /scratch/oww1c19/argyranthemum_transcriptomics/pipeline_6/chi_squared_test
cd /scratch/oww1c19/argyranthemum_transcriptomics/pipeline_6/chi_squared_test

# cp scripts to dir
cp /scratch/oww1c19/argyranthemum_transcriptomics/scripts/script_chi_squared_test_gene_trans_map.pbs .
cp /scratch/oww1c19/argyranthemum_transcriptomics/scripts/chi_squared_test_gene_trans_map.R .

# run script to get chi-squared test data
qsub script_chi_squared_test_gene_trans_map.pbs
```

### Chi-squared comparison

```
mkdir /scratch/oww1c19/argyranthemum_transcriptomics/pipeline_comparison

for i in {1..6}; do 
   cp /scratch/oww1c19/argyranthemum_transcriptomics/pipeline_${i}/chi_squared_test/chi_squared_data.txt /scratch/oww1c19/argyranthemum_transcriptomics/pipeline_comparison/chi_squared_data_pipeline_${i}.txt; 
done
```


