#!/bin/bash -l
# Start bash job in the current working directory AKA where I qsubbed from
#$ -cwd
# Name of job
#$ -N kneaddata_humann_rerun
# Output directory (for log)
#$ -o /u/scratch/m/mshono/humann/all_log/
# Join errors to one log
#$ -j y
# Runtime, storage and resources

#$ -l h_rt=24:00:00,h_data=2G
#$ -pe shared 32
# Array job:  one task per sample of the manifest
#$ -t 1-69
# Stop if command throws an error
set -e

# Set up

# Load  apptainer
module load apptainer/1.2.2

# Paths to kneaddata and humann containers
KD_SIF=/u/scratch/m/mshono/humann/kneaddata_0.10.0.sif
HM3_SIF=/u/scratch/m/mshono/humann/humann_3.9.sif

# Databases, and the dog index prefix
DB=/u/scratch/m/mshono/humann/databases
DOG_INDEX=$DB/dog_genome/dog_ROS_Cfam_1.0

# Where the test writes everything
ALL_OUT=/u/scratch/m/mshono/humann/all_out

# Make scratch AND the raw-data archive visible inside the containers
BIND="-B /u/scratch/m/mshono -B /u/project/pellegrini_archive"

# Sample ran is the line with $SGE_TASK ID of the manifest
MANIFEST=/u/scratch/m/mshono/humann/manifest_rerun_1.tsv

SAMPLE=$(awk -F'\t' -v n=$SGE_TASK_ID 'NR==n{print$1}' $MANIFEST)

# Stop: when the manifest lookup came back empty
if [ -z "$SAMPLE" ]; then
  echo "ERROR: no sample on this manifest line with the task id: $SGE_TASK_ID"
  exit 1
fi

# Look up this sample's read files in the manifest
R1=$(awk -F'\t' -v s="$SAMPLE" '$1==s {print $2}' $MANIFEST)
R2=$(awk -F'\t' -v s="$SAMPLE" '$1==s {print $3}' $MANIFEST)

# Output directories of kneaddata and humann3
KD_OUT=$ALL_OUT/kneaddata
HM3_OUT=$ALL_OUT/humann
mkdir -p $KD_OUT $HM3_OUT

echo "Sample: $SAMPLE  Node: $(hostname -s) Date: $(date)"

# Running KneadData

# Clean the reads: trim low-quality bases, then remove dog DNA
apptainer exec $BIND $KD_SIF kneaddata \
  --input $R1 \
  --input $R2 \
  --reference-db $DOG_INDEX \
  --output $KD_OUT \
  --output-prefix $SAMPLE \
  --threads 32 \
  --cat-final-output \
  --remove-intermediate-output

echo "KneadData complete at data!: $(date)"
ls -lh $KD_OUT || true

# Running humann

# --nucleotide-database / --protein-database: AKA  ChocoPhlAn and UniRef location
# --metaphlan-options: flags passed to MetaPhlAn when HUMAnN calls it
#   --bowtie2db =  database location, -x =  version to use
# --threads 32: matches the 32 cores requested above
apptainer exec $BIND $HM3_SIF humann \
  --input $KD_OUT/${SAMPLE}.fastq \
  --output $HM3_OUT/${SAMPLE} \
  --nucleotide-database $DB/chocophlan \
  --protein-database $DB/uniref \
  --metaphlan-options "--bowtie2db $DB/metaphlan -x mpa_vJun23_CHOCOPhlAnSGB_202403" \
  --threads 32

echo "FINISHED! on...  $(date)"
ls -lh $HM3_OUT/${SAMPLE} || true

# Cleanup step: keep the results tables, delete reads and alignment
# Rescue the species table out of the temp folder before deleting it
cp $HM3_OUT/${SAMPLE}/${SAMPLE}_humann_temp/${SAMPLE}_metaphlan_bugs_list.tsv $HM3_OUT/${SAMPLE}/

# Delete the alignment intermediates
rm -rf $HM3_OUT/${SAMPLE}/${SAMPLE}_humann_temp

# Delete this sample's FASTQ files, keeping the .log
rm -f $KD_OUT/${SAMPLE}*.fastq

echo "Cleanup done at date : $(date)"
ls -lh $HM3_OUT/${SAMPLE} || true

