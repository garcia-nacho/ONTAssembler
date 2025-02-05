# ONTAssembler   
A reproducible pipeline for de novo assembly of bacterial genomes using Oxford Nanopore (and optionally Illumina) sequencing data. 
This workflow combines long-read assembly (Flye), polishing (Racon + Medaka), and optional hybrid polishing with Illumina data (Pilon) to generate high-quality bacterial genomes.

## Key Features   
Optimized for accuracy: Two-step polishing (Racon + Medaka) resolves Nanopore-specific errors.   

Containerized workflow: Docker image for dependency management.   

Hybrid polishing: Optional integration with Illumina data using Pilon.   

QC and evaluation: Includes QUAST, CheckM, and NanoPlot for assembly metrics.   



## Prerequisites   
Docker (tested on v20.0+)   

### Computational Resources:   

CPUs: 8+ cores recommended   

RAM: 16+ GB   

Storage: 50 GB free space   

## Installation

No installation is needed since it runs from a docker image.   

<code>docker run -it --rm -v $(pwd):/Data ghcr.io/garcia-nacho/ontassembler</code>

## Input Structure
<pre>
inputfolder/   
    ├── sample1/    
    |      ├── file1.fastq.gz # Fastq Nanopore files   
    |      ├── file2.fastq.gz    
    |   
    ├── sample2/    
           ├── file1.fastq.gz # Fastq Nanopore files   
           ├── file2.fastq.gz     
</pre>

## Optional Parameters    
**Illumina Hybrid Polishing:** Place paired-end R1 and R2 fastq.gz files alongside the ONT files.    
Note that files containing the strings *"_R1"* and *"_R2"* will be exluded from the first steps of the pipeline and only used on the Illumina polishing steps    

**Threads:** Modify THREADS in run_assembly.sh (default: 16).    

Genome Size: Adjust GENOME_SIZE in run_assembly.sh (default: 4.1m).    

## Workflow Overview

<pre>graph TD
  A [FASTQ] -->|Guppy| B(FASTQ)
  A --> B{QC & Filtering}
  B -->|NanoFilt| C[Filtered Reads]
  C --> D[Flye Assembly]
  D --> E[Racon Polish]
  E --> F[Medaka Polish]
  F --> G{Illumina Data?}
  G -->|Yes| H[Pilon Polish]
  G -->|No| I[Evaluation]
 </pre>


## Expected Results
<pre>
data/output/
├── assembly/          # Flye draft assembly
├── polished/          # Racon + Medaka outputs
├── quast/             # Assembly evaluation report
├── checkm/            # Completeness/contamination metrics

</pre>

## License
This pipeline is MIT-licensed. 

## Contributing
Pull requests are welcome! For major changes, please open an issue first.

## Contact
For questions or support, contact *iggl [AT] fhi [DOT] no*
