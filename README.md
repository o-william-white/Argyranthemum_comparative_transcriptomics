
# file structure

### import raw data to iridis

# in iridis
mkdir /scratch/oww1c19/argyranthemum_transcriptomics
mkdir /scratch/oww1c19/argyranthemum_transcriptomics/raw_data

# in local terminal
rsync -avz --partial /drives/d/Trancriptomics/C101HW16120181/raw_data oww1c19@iridis4_c.soton.ac.uk:/scratch/oww1c19/argyranthemum_transcriptomics


### trimmomatic

mkdir /scratch/oww1c19/argyranthemum_transcriptomics/trimmomatic
cd /scratch/oww1c19/argyranthemum_transcriptomics/trimmomatic

qsub script_trimmomatic.pbs


### normalise reads 

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

# run normalisation
# note each script produces large temporary files so run idividually and check memory
qsub script_normalise_bro.pbs
qsub script_normalise_sun.pbs
qsub script_normalise_lem.pbs
qsub script_normalise_fru.pbs

# get input for normalisation of all species together
ls -1 normalise_*/trim*1.fq.gz_ext_all_reads.normalized_K25_maxC30_minC0_maxCV10000.fq | sed 's:normalise:/scratch/oww1c19/argyranthemum_transcriptomics/normalise_reads/normalise:g' > input_all_left
ls -1 normalise_*/trim*2.fq.gz_ext_all_reads.normalized_K25_maxC30_minC0_maxCV10000.fq | sed 's:normalise:/scratch/oww1c19/argyranthemum_transcriptomics/normalise_reads/normalise:g' > input_all_right

qsub script_normalise_all.pbs


### run trininty assemblies

mkdir /scratch/oww1c19/argyranthemum_transcriptomics/trinity_all
mkdir /scratch/oww1c19/argyranthemum_transcriptomics/trinity_separate


cd /scratch/oww1c19/argyranthemum_transcriptomics/trinity_separate

qsub script_trinity_bro.pbs
qsub script_trinity_fru.pbs
qsub script_trinity_lem.pbs
qsub script_trinity_sun.pbs

cd /scratch/oww1c19/argyranthemum_transcriptomics/trinity_all

qsub script_trinity_all.pbs






### summary stats

mkdir /scratch/oww1c19/argyranthemum_transcriptomics/summary_stats
cd /scratch/oww1c19/argyranthemum_transcriptomics/summary_stats

# create sequence id list
for i in {1..24}
do 
	echo A${i}_1 >> seq_id
	echo A${i}_2 >> seq_id
done


# create sample id list
for i in {1..6}
do
   echo bro_${i}_r1 >> sample_id
   echo bro_${i}_r2 >> sample_id
done

for i in {7..12}
do
   echo sun_${i}_r1 >> sample_id
   echo sun_${i}_r2 >> sample_id
done

for i in {13..18}
do
   echo lem_${i}_r1 >> sample_id
   echo lem_${i}_r2 >> sample_id
done

for i in {19..24}
do
   echo fru_${i}_r1 >> sample_id
   echo fru_${i}_r2 >> sample_id
done


# count raw reads
qsub script_count_raw_reads.pbs

# count trimmomatic reads
qsub script_count_trimmomatic_reads.pbs

# create summary stats table
paste seq_id sample_id raw_reads trimmomatic_reads > summary_stats

