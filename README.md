### STAR alignment using Snakemake

#### 1. Conda environment

- References: [Conda doc](https://docs.conda.io/projects/conda/en/latest/index.html), [sra-tools](https://github.com/ncbi/sra-tools), [Snakemake Installation Guide](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html)
- Recipe: [config/conda_env.yml](https://github.com/Mira0507/snakemake_star/blob/master/config/conda_env.yml)


#### 2. Snakemake 

- Reference: [Snakemake doc](https://snakemake.readthedocs.io/en/stable), [STAR manual](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf)

- [Snakefile](https://github.com/Mira0507/snakemake_star/blob/master/Snakefile)

```
#################################### Defined by users #################################
configfile:"config/config_paired.yaml"    # Sets path to the config file

#######################################################################################

THREADS=config["THREADS"]

shell.prefix('set -euo pipefail; ')
shell.executable('/bin/bash')

rule all:
    input:
        expand("star_output/{sample}Aligned.sortedByCoord.out.bam", 
               sample=config["INPUT_PREFIX"])



rule align_star: 
    """
    This rule aligns the reads using STAR two-pass mode
    """
    input:
        gtf=config["GTF"]
    output:
        expand("star_output/{sample}Aligned.sortedByCoord.out.bam", sample=config['INPUT_PREFIX'])  
    params:
        indexing=config["INDEX_STAR"],
        indir=config['INPUT_DIR'],
        files=config["INPUT_PREFIX"], 
        read_ends=config['INPUT_END'],
        ext=config['INPUT_EXT']
    run:
        for i in range(len(params.files)):
            p=params.files[i]
            r1= params.indir + params.files[i] + "_1" + params.ext + " " 
            r2=""
            if len(params.read_ends) == 2: 
                r2= params.indir + params.files[i] + "_2" + params.ext + " " 
            shell("STAR --runThreadN {THREADS} "  
                    "--runMode alignReads "  
                    "--readFilesCommand zcat "
                    "--genomeDir {params.indexing} " 
                    "--sjdbGTFfile {input.gtf} "  
                    "--sjdbOverhang 100 "  
                    "--readFilesIn {r1}{r2}"  
                    "--outFileNamePrefix star_output/{p} "
                    "--outFilterType BySJout "  
                    "--outFilterMultimapNmax 20 "
                    "--alignSJoverhangMin 8 "
                    "--alignSJDBoverhangMin 1 "
                    "--outFilterMismatchNmax 999 "
                    "--outFilterMismatchNoverReadLmax 0.04 "
                    "--alignIntronMin 20 "
                    "--alignIntronMax 1000000 "
                    "--outSAMunmapped None "
                    "--outSAMtype BAM "
                    "SortedByCoordinate "
                    "--quantMode GeneCounts "
                    "--twopassMode Basic "
                    "--chimOutType Junctions")
```

- [config/config_single.yaml (single-end testing)](https://github.com/Mira0507/snakemake_star/blob/master/config/config_single.yaml)


```
GTF: "reference/gencode.v36.primary_assembly.annotation.gtf" # Assigns path to .gtf file

INDEX_STAR: "/home/mira/Documents/programming/Bioinformatics/snakemake_star_hisat/reference/star_index" # Assigns path to star indexing directory (ABSOLUTE PATH NEEDED due to a potential STAR error!)

INDEX_HISAT: "reference/hisat2_index" # Assigns path to hisat2 indexing directory

THREADS: 8     # Assigns the number of threads

STAR_OUT: "star_output"  # Assigns STAR output directory

HISAT2_OUT: "hisat2_output"  # Assigns HISAT2 output directory

INPUT_DIR: 'fastq_se/'

INPUT_PREFIX:
  - 'Control'
  - 'Treatment'

INPUT_END:
  - '_1'
# - '_2' # for paired end reads

INPUT_EXT: '.fastq.gz'

```


- [config/config_paired.yaml (paired-end testing)](https://github.com/Mira0507/snakemake_star/blob/master/config/config_paired.yaml)


```
GTF: "reference/gencode.v36.primary_assembly.annotation.gtf" # Assigns path to .gtf file

INDEX_STAR: "/home/mira/Documents/programming/Bioinformatics/snakemake_star_hisat/reference/star_index" # Assigns path to star indexing directory (ABSOLUTE PATH NEEDED due to a potential STAR error!)

INDEX_HISAT: "reference/hisat2_index" # Assigns path to hisat2 indexing directory

THREADS: 8     # Assigns the number of threads

STAR_OUT: "star_output"  # Assigns STAR output directory

HISAT2_OUT: "hisat2_output"  # Assigns HISAT2 output directory

INPUT_DIR: 'fastq_pe/'   # Contains input fastq.gz files

INPUT_PREFIX:
  - 'vector'
  - 'PANX'

INPUT_END:
  - '_1'
  - '_2' # for paired end reads

INPUT_EXT: '.fastq.gz'
  
```

#### 3. Running the Snakemake workflow

- Reference: [Snakemake Command Line Arguments](https://snakemake.readthedocs.io/en/stable/executing/cli.html)

- **Dry run**


```bash
#!/bin/bash

snakemake -n

```


- **DAG visualization** 


```bash
#!/bin/bash


# The dot commend requires graphviz (downloadable via conda)
snakemake --dag | dot -Tpdf > dag.pdf


```


- **Run**


```bash
#!/bin/bash

# Either -j or --cores assignes the number of cores
snakemake -j 8

```
