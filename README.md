### STAR alignment using Snakemake

#### 1. Conda environment

- References: [Conda doc](https://docs.conda.io/projects/conda/en/latest/index.html), [sra-tools](https://github.com/ncbi/sra-tools), [Snakemake Installation Guide](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html)
- Recipe: [config/conda_env.yml](https://github.com/Mira0507/snakemake_star/blob/master/config/conda_env.yml)


#### 2. Snakemake 

- Reference: [Snakemake doc](https://snakemake.readthedocs.io/en/stable), [STAR manual](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf)

- [Snakefile](https://github.com/Mira0507/snakemake_star/blob/master/Snakefile)

```


#################################### Defined by users #################################
configfile:"config/config_single.yaml"    # Sets path to the config file

#######################################################################################

# This workflow requires fastq.gz files in fastq directory 
# e.g. paired-end: DMSO_rep1_1.fastq.gz, DMSO_rep1_2.fastq.gz, Drug_rep1_1.fastq.gz, Drug_rep1_2.fastq.gz 
#      single-end: DMSO_rep1_1.fastq.gz, Drug_rep1_1.fastq.gz

shell.prefix('set -euo pipefail; ')
shell.executable('/bin/bash')

rule all:
    input:
        expand("star_output/{sample}Aligned.sortedByCoord.out.bam", 
               sample=config["FASTQ_PREFIX"])

rule get_reference:    
    """
    This rule downloads reference files
    """
    params:
        gen_link=config['REFERENCE_LINK']['GENOME'][0],   # Gencode reference genome file link 
        gen_name=config['REFERENCE_LINK']['GENOME'][1],   # Output reference genome location & name 
        anno_link=config['REFERENCE_LINK']['ANNOTATION'][0],  # Gencode GTF (annotation) file link
        anno_name=config['REFERENCE_LINK']['ANNOTATION'][1]   # Output GTF file location & name
    output:
        gen=expand("reference/{gen}", gen=config['REFERENCE_LINK']['GENOME'][2]),  # Decompressed reference genome file 
        anno=expand("reference/{anno}", anno=config['REFERENCE_LINK']['ANNOTATION'][2])  # Decompressed GTF file
    shell:
        "set +o pipefail; "
        "wget -c {params.gen_link} -O reference/{params.gen_name} && "
        "wget -c {params.anno_link} -O reference/{params.anno_name} && "
        "gzip -d reference/*.gz"

rule index_star:
    """
    This rule constructs STAR index files
    """
    input:
        fa=expand("reference/{gen}", gen=config['REFERENCE_LINK']['GENOME'][2]),  # Decompressed reference genome file
        gtf=expand("reference/{anno}", anno=config['REFERENCE_LINK']['ANNOTATION'][2])  # Decompressed GTF file
    output:
        "reference/star_index/Genome",   # STAR indexing files
        "reference/star_index/SA",       # STAR indexing files
        "reference/star_index/SAindex"   # STAR indexing files
    threads: 8
    shell:
        "set +o pipefail; "
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
        gtf=expand("reference/{gen}", gen=config['REFERENCE_LINK']['ANNOTATION'][2]),  # Decompressed GTF file
        fastq=expand("fastq/{sample}_{end}.fastq.gz", sample=config['FASTQ_PREFIX'], end=config['END']),                  # Gzipped FASTQ files
        index1="reference/star_index/Genome",
        index2="reference/star_index/SA",
        index3="reference/star_index/SAindex"
    output:
        expand("star_output/{sample}Aligned.sortedByCoord.out.bam", sample=config['FASTQ_PREFIX'])  # Bam files
    params:
        indexing=config["INDEX_STAR"],  # STAR indexing file directory
        files=config["FASTQ_PREFIX"],   # e.g. Ctrl, Treatment
        ext=config['FASTQ_EXT']         # extension of the FASTQ files (e.g. fastq.gz)
    threads: 8
    run:
        for i in range(len(params.files)):
            p=params.files[i]
            r1= "fastq/" + params.files[i] + "_1" + params.ext + " " 
            r2=""
            if len(input.fastq) == 2 * len(params.files): 
                r2= "fastq/" + params.files[i] + "_2" + params.ext + " " 
            shell("STAR --runThreadN {threads} "  
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


```yaml


###################### Sample info ######################

END: 
  - 1



FASTQ_PREFIX:
  - DMSO_rep1
  - DMSO_rep2
  - DMSO_rep3
  - SR0813_rep1
  - SR0813_rep2
  - SR0813_rep3


###################### Reference info ######################


REFERENCE_LINK:
  GENOME: 
    - 'ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_37/GRCh38.primary_assembly.genome.fa.gz'
    - 'GRCh38.primary_assembly.genome.fa.gz'
    - 'GRCh38.primary_assembly.genome.fa'
  ANNOTATION: 
    - 'ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_37/gencode.v37.primary_assembly.annotation.gtf.gz'
    - 'gencode.v37.primary_assembly.annotation.gtf.gz'
    - 'gencode.v37.primary_assembly.annotation.gtf'


###################### Extra-setting info ######################

INDEX_STAR: "/home/mira/Documents/programming/Bioinformatics/snakemake_star_hisat/reference/star_index" # Assigns path to star indexing directory (ABSOLUTE PATH NEEDED due to a potential STAR error!)




STAR_OUT: "star_output"  # Assigns STAR output directory



FASTQ_EXT: '.fastq.gz'


```


- [config/config_paired.yaml (paired-end testing)](https://github.com/Mira0507/snakemake_star/blob/master/config/config_paired.yaml)


```yaml

###################### Sample info ######################

END: 
  - 1
  - 2


  
FASTQ_PREFIX:
  - Untreated1
  - Untreated2
  - Treated1
  - Treated2

###################### Reference info ######################


REFERENCE_LINK:
  GENOME: 
    - 'ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_37/GRCh38.primary_assembly.genome.fa.gz'
    - 'GRCh38.primary_assembly.genome.fa.gz'
    - 'GRCh38.primary_assembly.genome.fa'
  ANNOTATION: 
    - 'ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_37/gencode.v37.primary_assembly.annotation.gtf.gz'
    - 'gencode.v37.primary_assembly.annotation.gtf.gz'
    - 'gencode.v37.primary_assembly.annotation.gtf'


###################### Extra-setting info ######################

INDEX_STAR: "/home/mira/Documents/programming/Bioinformatics/snakemake_star_hisat/reference/star_index" # Assigns path to star indexing directory (ABSOLUTE PATH NEEDED due to a potential STAR error!)




STAR_OUT: "star_output"  # Assigns STAR output directory



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
snakemake -j 8

```
