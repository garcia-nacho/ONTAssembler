#!/bin/bash

set -e

# Configuration
INPUT_DIR="/data/input"       # Mount your Nanopore data here (FAST5 or FASTQ)
OUTPUT_DIR="/data/output"     # Results will be saved here
GENOME_SIZE="4.1m"            # B. pertussis genome size
THREADS=10                    # Number of CPU threads
M=500
inputfolder=${1}
homedir=$(pwd)

mkdir ${homedir}/ont_assembly
mkdir ${homedir}/ont_assembly/fasta
mkdir ${homedir}/ont_assembly/qc
mkdir ${homedir}/ont_assembly/qc/nanoplot


cd ${inputfolder}
for dir in $(ls -d */)
do
cd ${dir}

mkdir ${homedir}/ont_assembly/fasta/${dir%/}

cat *.fastq.gz > ${dir%/}_reads.fastq.gz

NanoPlot --fastq ${dir%/}_reads.fastq.gz --outdir ${homedir}/ont_assembly/qc/nanoplot
cat ${dir%/}_reads.fastq.gz | NanoFilt -q 10 -l ${M} > ${dir%/}_filtered.fastq.gz


echo "Running Flye assembly..."
flye --nano-raw ${dir%/}_filtered.fastq.gz \
  --genome-size $GENOME_SIZE \
  --out-dir ${homedir}/ont_assembly/fasta/${dir%/} \
  --threads ${THREADS}

echo "Polishing with Racon..."
minimap2 -x map-ont ${homedir}/ont_assembly/fasta/${dir%/}/assembly/assembly.fasta ${dir%/}_filtered.fastq.gz > ${dir%/}_overlaps.paf
racon -t ${THREADS} ${dir%/}_filtered.fastq.gz ${dir%/}_overlaps.paf ${homedir}/ont_assembly/fasta/${dir%/}/assembly/assembly.fasta > ${homedir}/ont_assembly/fasta/${dir%/}/${dir%/}_racon.fasta

echo "Polishing with Medaka..."
source activate medaka
medaka_consensus -i ${dir%/}_filtered.fastq.gz -d ${homedir}/ont_assembly/fasta/${dir%/}/${dir%/}_racon.fasta -o ${homedir}/ont_assembly/fasta/${dir%/} -t ${THREADS}
mv ${homedir}/ont_assembly/fasta/${dir%/}/consensus.fasta ${homedir}/ont_assembly/fasta/${dir%/}/${dir%/}_polished.fasta 
conda deactivate

echo "Running QUAST..."
quast.py ${homedir}/ont_assembly/fasta/${dir%/}/${dir%/}_polished.fasta -o ${homedir}/ont_assembly/qc


echo "Running CheckM..."
checkm lineage_wf -x fa ${homedir}/ont_assembly/fasta/${dir%/} ${homedir}/ont_assembly/qc

#Check for R1/R2
#bwa index medaka_polished.fasta
#bwa mem -t ${THREADS} medaka_polished.fasta ${r1} ${r2} | samtools sort -@ ${THREADS} -o illumina.aligned.sorted.bam -
#pilon --genome polished.fasta --frags illumina.bam --output pilon_polished

cd ..
done

echo "Pipeline completed!"
