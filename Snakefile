
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


