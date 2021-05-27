### STAR alignment using Snakemake

#### 1. Conda environment

- References: [Conda doc](https://docs.conda.io/projects/conda/en/latest/index.html), [sra-tools](https://github.com/ncbi/sra-tools), [Snakemake Installation Guide](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html)
- Config: [config/conda_env.yaml](https://github.com/Mira0507/snakemake_star/blob/master/config/conda_env.yaml)


#### 2. Snakemake 

- Reference: [Snakemake doc](https://snakemake.readthedocs.io/en/stable), [STAR manual](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf)

- [Snakefile](https://github.com/Mira0507/snakemake_star/blob/master/Snakefile)

```


#################################### Defined by users #################################
configfile:"config/config_paired1.yaml"    # Sets path to the config file

#######################################################################################

# This workflow requires fastq.gz files in fastq directory 
# e.g. paired-end: DMSO_rep1_1.fastq.gz, DMSO_rep1_2.fastq.gz, Drug_rep1_1.fastq.gz, Drug_rep1_2.fastq.gz 
#      single-end: DMSO_rep1_1.fastq.gz, Drug_rep1_1.fastq.gz




rule all: 
    input: 
        expand("reference/{ref}", ref=list(config['REFERENCE']['FILE'].values())),   # Reference genome and annotation (GTF) files
        expand("star_output/{sample}.bam", sample=list(config['SAMPLE'].keys()))  # STAR output BAM files

rule get_reference:    
    """
    This rule downloads and decompresses reference files
    """
    params:
        reflink=config['REFERENCE']['LINK']
    output:
        "reference/{ref}"  # Decompressed reference files
    run:
        link=params.reflink + wildcards.ref
        shell("wget -c {link}.gz -O {output}.gz && " 
              "gzip -d {output}.gz")



rule index_star:
    """
    This rule constructs STAR index files
    """
    input:
        fa=expand("reference/{gen}", gen=config['REFERENCE']['FILE']['GENOME']),  # Decompressed reference genome file
        gtf=expand("reference/{anno}", anno=config['REFERENCE']['FILE']['ANNOTATION'])  # Decompressed GTF file
    output:
        expand("reference/star_index/{index}", index=config['INDEX_STAR']['FILE'])
    threads: 16
    shell:
        "STAR --runThreadN {threads} "
        "--runMode genomeGenerate "
        "--genomeDir reference/star_index "
        "--genomeFastaFiles {input.fa} "
        "--sjdbGTFfile {input.gtf}"

rule align_star:   # Creates bam files in star_output directory"
    """
    This rule aligns the reads using STAR two-pass mode
    """
    input:
        gtf=expand("reference/{anno}", anno=config['REFERENCE']['FILE']['ANNOTATION']),   # Decompressed GTF file
        fastq=expand("fastq/{{sample}}_{end}.fastq.gz", end=config['END']),    # Gzipped FASTQ files
        index=expand("reference/star_index/{index}", index=config['INDEX_STAR']['FILE']) # STAR indexing files
    output:
        "star_output/{sample}.bam"     # Bam files
    params:
        indexing=config["INDEX_STAR"]['DIR'],  # STAR indexing file directory
        ext=config['FASTQ_EXT']         # extension of the FASTQ files (e.g. fastq.gz)
    threads: 16
    run:
        r1= "fastq/" + wildcards.sample + "_1" + params.ext + " " 
        r2=""
        if len(input.fastq) == 2:   # if paired-end
            r2= "fastq/" + wildcards.sample + "_2" + params.ext + " " 
        shell("STAR --runThreadN {threads} "  
                "--runMode alignReads "  
                "--readFilesCommand zcat "
                "--genomeDir {params.indexing} " 
                "--sjdbGTFfile {input.gtf} "  
                "--sjdbOverhang 100 "  
                "--readFilesIn {r1}{r2}"  
                "--outFileNamePrefix star_output/{wildcards.sample} "
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
                "--chimOutType Junctions && "
              "mv star_output/{wildcards.sample}Aligned.sortedByCoord.out.bam {output}")   

```

- [config/config_single.yaml (single-end testing)](https://github.com/Mira0507/snakemake_star/blob/master/config/config_single.yaml)


```yaml


###################### Sample info ######################

SAMPLE:
  DMSO_rep1: SRR13190144
  DMSO_rep2: SRR13190145
  DMSO_rep3: SRR13190146
  SR0813_rep1: SRR13190150
  SR0813_rep2: SRR13190151
  SR0813_rep3: SRR13190152


END: 
  - 1


###################### Reference info ######################


REFERENCE:
  LINK: 'ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_37/'
  FILE: 
    GENOME: 'GRCh38.primary_assembly.genome.fa'
    ANNOTATION: 'gencode.v37.primary_assembly.annotation.gtf'
  
# e.g. 
# Reference genome link: ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_37/GRCh38.primary_assembly.genome.fa.gz
# Reference annotation (GTF) link: ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_37/gencode.v37.primary_assembly.annotation.gtf.gz 
#
#
###################### Extra-setting info ######################

INDEX_STAR: 
  DIR: "reference/star_index" # Assigns path to star indexing directory (assign ABSOLUTE PATH if STAR generates an error)
  FILE: 
    - "Genome"
    - "SA"
    - "SAindex"


FASTQ_EXT: '.fastq.gz'
```


- [config/config_paired1.yaml (paired-end testing)](https://github.com/Mira0507/snakemake_star/blob/master/config/config_paired1.yaml)


```yaml



###################### Sample info ######################


SAMPLE:
  Treated_rep1: SRR6461133
  Treated_rep2: SRR6461134
  Treated_rep3: SRR6461135
  Control_rep1: SRR6461139 
  Control_rep2: SRR6461140
  Control_rep3: SRR6461141



END: 
  - 1
  - 2

###################### Reference info ########################

REFERENCE:
  LINK: 'ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_37/'
  FILE: 
    GENOME: 'GRCh38.primary_assembly.genome.fa'
    ANNOTATION: 'gencode.v37.primary_assembly.annotation.gtf'
  
# e.g. 
# Reference genome link: ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_37/GRCh38.primary_assembly.genome.fa.gz
# Reference annotation (GTF) link: ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_37/gencode.v37.primary_assembly.annotation.gtf.gz 



###################### Extra-setting info ######################

INDEX_STAR: 
  DIR: "reference/star_index" # Assigns path to star indexing directory (assign ABSOLUTE PATH if STAR generates an error)
  FILE: 
    - "Genome"
    - "SA"
    - "SAindex"



FASTQ_EXT: '.fastq.gz'

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
snakemake -j 10

```
