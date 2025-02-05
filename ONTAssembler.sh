#!/bin/bash

set -e

# Configuration to be fixed on the docker container

GENOME_SIZE="4.1m"            
THREADS=10                    
M=500
homedir=/Data

echo "Running ONT assembly..."

for dir in $(ls -d */)
do
cd ${dir}

mkdir -p ${homedir}/ont_assembly/fasta/${dir%/}
mkdir -p ${homedir}/ont_assembly/qc/nanoplot_${dir%/}
mkdir -p ${homedir}/ont_assembly/qc/quast_${dir%/}
mkdir -p ${homedir}/ont_assembly/qc/checkm_${dir%/}

ls *.fastq.gz | grep -v -E "_R1|_R2" | xargs cat > ${dir%/}_reads.fastq.gz


gzip -d ${dir%/}_reads.fastq.gz
set +e
NanoPlot --fastq ${dir%/}_reads.fastq --outdir ${homedir}/ont_assembly/qc/nanoplot_${dir%/}
set -e
cat ${dir%/}_reads.fastq | NanoFilt -q 10 -l ${M} > ${dir%/}_filtered.fastq
rm ${dir%/}_reads.fastq


echo "Running Flye assembly..."
flye --nano-raw ${dir%/}_filtered.fastq --genome-size ${GENOME_SIZE} --out-dir ${homedir}/ont_assembly/fasta/${dir%/} --threads ${THREADS}
mv ${homedir}/ont_assembly/fasta/${dir%/}/assembly.fasta ${homedir}/ont_assembly/fasta/${dir%/}_flye.fasta

echo "Polishing with Racon..."
minimap2 -x map-ont ${homedir}/ont_assembly/fasta/${dir%/}_flye.fasta ${dir%/}_filtered.fastq > ${dir%/}_overlaps.paf
racon -t ${THREADS} ${dir%/}_filtered.fastq ${dir%/}_overlaps.paf ${homedir}/ont_assembly/fasta/${dir%/}_flye.fasta > ${homedir}/ont_assembly/fasta/${dir%/}_racon.fasta
rm ${dir%/}_overlaps.paf

echo "Polishing with Medaka..."
source activate medaka
medaka_consensus -i ${dir%/}_filtered.fastq -d ${homedir}/ont_assembly/fasta/${dir%/}_racon.fasta -o ${homedir}/ont_assembly/fasta/${dir%/} -t ${THREADS}
mv ${homedir}/ont_assembly/fasta/${dir%/}/consensus.fasta ${homedir}/ont_assembly/fasta/${dir%/}_polished_raw.fasta 
fold -w 80 ${homedir}/ont_assembly/fasta/${dir%/}_polished_raw.fasta > ${homedir}/ont_assembly/fasta/${dir%/}_polished.fasta
rm ${homedir}/ont_assembly/fasta/${dir%/}_polished_raw.fasta
rm ${homedir}/ont_assembly/fasta/*.fasta.fai
rm ${homedir}/ont_assembly/fasta/*ont.mmi
conda deactivate

echo "Running QUAST..."
quast.py ${homedir}/ont_assembly/fasta/${dir%/}_polished.fasta -o ${homedir}/ont_assembly/qc/quast_${dir%/}

echo "Running CheckM..."
mkdir bins
cp ${homedir}/ont_assembly/fasta/${dir%/}_polished.fasta bins/${dir%/}_polished.fasta
checkm lineage_wf -x fasta bins ${homedir}/ont_assembly/qc/checkm_${dir%/}
rm -rf bins

#Check for R1/R2
if ls *_R1* 1> /dev/null 2>&1; then
    echo "Illumina files found found!"
    r1=$(ls *_R1*)
    r2=$(ls *_R2*)
    cp ${homedir}/ont_assembly/fasta/${dir%/}_polished.fasta ./${dir%/}_polished.fasta
    bwa index ${dir%/}_polished.fasta
    bwa mem -t ${THREADS} ${dir%/}_polished.fasta ${r1} ${r2} | samtools sort -@ ${THREADS} -o illumina.aligned.sorted.bam -
    samtools index illumina.aligned.sorted.bam    
    java -Xmx32G -jar /home/docker/miniconda3/share/pilon-1.24-0/pilon.jar --genome ${dir%/}_polished.fasta --frags illumina.aligned.sorted.bam --output pilon_corrected
    mv pilon_corrected.fasta ${homedir}/ont_assembly/fasta/${dir%/}_superpolished.fasta

    
    rm ${dir%/}_polished*  
    rm illumina.aligned.sorted.bam
    rm illumina.aligned.sorted.bam.bai      
fi

rm ${dir%/}_filtered.fastq

cd ..
done

echo "Pipeline completed!"
