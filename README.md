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

inputfolder/   
    ├── sample1/    
    |      ├── file1.fastq.gz # Fastq Nanopore files   
    |      ├── file2.fastq.gz    
    |   
    ├── sample2/    
           ├── file1.fastq.gz # Fastq Nanopore files   
           ├── file2.fastq.gz     


Optional Parameters
Illumina Hybrid Polishing: Place paired-end reads in data/input/illumina_R1.fastq and data/input/illumina_R2.fastq.

Threads: Modify THREADS in run_assembly.sh (default: 16).

Genome Size: Adjust GENOME_SIZE in run_assembly.sh (default: 4.1m).

Workflow Overview
mermaid
Copy
graph TD
  A[FAST5] -->|Guppy| B(FASTQ)
  B --> C{QC & Filtering}
  C -->|NanoFilt| D[Filtered Reads]
  D --> E[Flye Assembly]
  E --> F[Racon Polish]
  F --> G[Medaka Polish]
  G --> H{Illumina Data?}
  H -->|Yes| I[Pilon Polish]
  H -->|No| J[Evaluation]
  I --> J
  J --> K[Circularization]
  K --> L[Annotation]
Expected Results
Copy
data/output/
├── assembly/          # Flye draft assembly
├── polished/          # Racon + Medaka outputs
├── quast/             # Assembly evaluation report
├── checkm/            # Completeness/contamination metrics
├── circularized/      # Circularized genome
└── annotation/        # Prokka gene annotations
Troubleshooting
Common Issues
GLIBCXX Errors:
Rebuild the Docker image using the updated Dockerfile in this repo.

Guppy Installation:
Mount your licensed Guppy directory to /opt/ont/guppy in the container.

CheckM Database:
Ensure the database is downloaded to /opt/checkm_data inside the container.

Citation
If you use this pipeline, please cite:

Flye: DOI:10.1038/s41587-019-0072-8

Medaka: GitHub

Prokka: DOI:10.1093/bioinformatics/btu153

See CITATIONS.md for full tool references.

License
This pipeline is MIT-licensed. Note that Guppy requires a separate license from Oxford Nanopore Technologies.

Contributing
Pull requests are welcome! For major changes, please open an issue first.

Contact
For questions or support, contact your.email@example.com.

This README provides clear instructions for users to run your pipeline and understand its components. Adjust paths, parameters, and contact details as needed!
