
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
        
