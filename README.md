### STAR alignment using Snakemake

#### 1. Conda environment

- References: [Conda doc](https://docs.conda.io/projects/conda/en/latest/index.html), [sra-tools](https://github.com/ncbi/sra-tools), [Snakemake Installation Guide](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html)
- [config/conda_env.yml](https://github.com/Mira0507/snakemake_star/blob/master/config/conda_env.yml)

```
name: snakemake_mapping
channels:
  - conda-forge
  - bioconda
  - defaults
dependencies:
  - _libgcc_mutex=0.1=conda_forge
  - _openmp_mutex=4.5=1_gnu
  - aiohttp=3.7.4=py39h3811e60_0
  - amply=0.1.4=py_0
  - appdirs=1.4.4=pyh9f0ad1d_0
  - async-timeout=3.0.1=py_1000
  - atk-1.0=2.36.0=h3371d22_4
  - attrs=20.3.0=pyhd3deb0d_0
  - azure-common=1.1.27=pyhd8ed1ab_0
  - azure-core=1.12.0=pyhd8ed1ab_0
  - azure-storage-blob=12.8.0=pyhd8ed1ab_0
  - bedtools=2.29.2=hc088bd4_0
  - blinker=1.4=py_1
  - boto=2.49.0=py_0
  - boto3=1.17.49=pyhd8ed1ab_0
  - botocore=1.20.49=pyhd8ed1ab_0
  - brotlipy=0.7.0=py39h3811e60_1001
  - bzip2=1.0.8=h7f98852_4
  - c-ares=1.17.1=h7f98852_1
  - ca-certificates=2020.12.5=ha878542_0
  - cachetools=4.2.1=pyhd8ed1ab_0
  - cairo=1.16.0=hcf35c78_1003
  - certifi=2020.12.5=py39hf3d152e_1
  - cffi=1.14.5=py39he32792d_0
  - chardet=4.0.0=py39hf3d152e_1
  - coincbc=2.10.5=hcee13e7_1
  - configargparse=1.4=pyhd8ed1ab_0
  - cryptography=3.4.7=py39hbca0aa6_0
  - datrie=0.8.2=py39h07f9747_1
  - docutils=0.17=py39hf3d152e_0
  - expat=2.3.0=h9c3ff4c_0
  - filelock=3.0.12=pyh9f0ad1d_0
  - font-ttf-dejavu-sans-mono=2.37=hab24e00_0
  - font-ttf-inconsolata=2.001=hab24e00_0
  - font-ttf-source-code-pro=2.030=hab24e00_0
  - font-ttf-ubuntu=0.83=hab24e00_0
  - fontconfig=2.13.1=hba837de_1005
  - fonts-conda-ecosystem=1=0
  - fonts-conda-forge=1=0
  - freetype=2.10.4=h0708190_1
  - fribidi=1.0.10=h36c2ea0_0
  - gawk=5.1.0=h7f98852_0
  - gdk-pixbuf=2.42.6=h04a7f16_0
  - gettext=0.19.8.1=h0b5b191_1005
  - giflib=5.2.1=h36c2ea0_2
  - gitdb=4.0.7=pyhd8ed1ab_0
  - gitpython=3.1.14=pyhd8ed1ab_0
  - glib=2.68.1=h9c3ff4c_0
  - glib-tools=2.68.1=h9c3ff4c_0
  - google-api-core=1.26.2=pyhd8ed1ab_0
  - google-auth=1.28.0=pyh44b312d_0
  - google-cloud-core=1.5.0=pyhd3deb0d_0
  - google-cloud-storage=1.19.0=py_0
  - google-crc32c=1.1.2=py39hb81f231_0
  - google-resumable-media=1.2.0=pyhd3deb0d_0
  - googleapis-common-protos=1.53.0=py39hf3d152e_0
  - graphite2=1.3.13=h58526e2_1001
  - graphviz=2.42.3=h6939c30_2
  - grpcio=1.37.0=py39hff7568b_0
  - gtk2=2.24.33=hab0c2f8_0
  - gts=0.7.6=h64030ff_2
  - harfbuzz=2.4.0=h9f30f68_3
  - hisat2=2.2.1=he1b5a44_2
  - htslib=1.11=hd3b49d5_2
  - icu=64.2=he1b5a44_1
  - idna=2.10=pyh9f0ad1d_0
  - importlib-metadata=3.10.0=py39hf3d152e_0
  - ipython_genutils=0.2.0=py_1
  - isodate=0.6.0=py_1
  - jemalloc=5.2.1=h9c3ff4c_5
  - jmespath=0.10.0=pyh9f0ad1d_0
  - jpeg=9d=h36c2ea0_0
  - jsonschema=3.2.0=pyhd8ed1ab_3
  - jupyter_core=4.7.1=py39hf3d152e_0
  - krb5=1.17.2=h926e7f8_0
  - ld_impl_linux-64=2.35.1=hea4e1c9_2
  - libblas=3.9.0=8_openblas
  - libcblas=3.9.0=8_openblas
  - libcrc32c=1.1.1=h9c3ff4c_2
  - libcurl=7.75.0=hc4aaa36_0
  - libdeflate=1.7=h7f98852_5
  - libedit=3.1.20191231=he28a2e2_2
  - libev=4.33=h516909a_1
  - libffi=3.3=h58526e2_2
  - libgcc-ng=9.3.0=h2828fa1_18
  - libgfortran-ng=9.3.0=hff62375_18
  - libgfortran5=9.3.0=hff62375_18
  - libglib=2.68.1=h3e27bee_0
  - libgomp=9.3.0=h2828fa1_18
  - libiconv=1.16=h516909a_0
  - liblapack=3.9.0=8_openblas
  - libnghttp2=1.43.0=h812cca2_0
  - libopenblas=0.3.12=pthreads_h4812303_1
  - libpng=1.6.37=h21135ba_2
  - libprotobuf=3.15.7=h780b84a_0
  - libssh2=1.9.0=ha56f1ee_6
  - libstdcxx-ng=9.3.0=h6de172a_18
  - libtiff=4.2.0=hdc55705_0
  - libtool=2.4.6=h58526e2_1007
  - libuuid=2.32.1=h7f98852_1000
  - libwebp=1.2.0=h3452ae3_0
  - libwebp-base=1.2.0=h7f98852_2
  - libxcb=1.13=h7f98852_1003
  - libxml2=2.9.10=hee79883_0
  - lz4-c=1.9.3=h9c3ff4c_0
  - msrest=0.6.21=pyh44b312d_0
  - multidict=5.1.0=py39h3811e60_1
  - nbformat=5.1.3=pyhd8ed1ab_0
  - ncurses=6.2=h58526e2_4
  - oauthlib=3.0.1=py_0
  - openssl=1.1.1k=h7f98852_0
  - packaging=20.9=pyh44b312d_0
  - pango=1.42.4=h7062337_4
  - pcre=8.44=he1b5a44_0
  - perl=5.32.0=h36c2ea0_0
  - pip=21.0.1=pyhd8ed1ab_0
  - pixman=0.38.0=h516909a_1003
  - protobuf=3.15.7=py39he80948d_0
  - psutil=5.8.0=py39h3811e60_1
  - pthread-stubs=0.4=h36c2ea0_1001
  - pulp=2.4=py39hf3d152e_0
  - pyasn1=0.4.8=py_0
  - pyasn1-modules=0.2.7=py_0
  - pycparser=2.20=pyh9f0ad1d_2
  - pyjwt=2.0.1=pyhd8ed1ab_1
  - pyopenssl=20.0.1=pyhd8ed1ab_0
  - pyparsing=2.4.7=pyh9f0ad1d_0
  - pyrsistent=0.17.3=py39h3811e60_2
  - pysocks=1.7.1=py39hf3d152e_3
  - python=3.9.2=hffdb5ce_0_cpython
  - python-dateutil=2.8.1=py_0
  - python_abi=3.9=1_cp39
  - pytz=2021.1=pyhd8ed1ab_0
  - pyyaml=5.4.1=py39h3811e60_0
  - ratelimiter=1.2.0=py_1002
  - readline=8.0=he28a2e2_2
  - requests=2.25.1=pyhd3deb0d_0
  - requests-oauthlib=1.3.0=pyh9f0ad1d_0
  - rsa=4.7.2=pyh44b312d_0
  - s3transfer=0.3.6=pyhd8ed1ab_0
  - salmon=1.4.0=hf69c8f4_0
  - samtools=1.11=h6270b1f_0
  - setuptools=49.6.0=py39hf3d152e_3
  - six=1.15.0=pyh9f0ad1d_0
  - smart_open=5.0.0=pyhd8ed1ab_0
  - smmap=3.0.5=pyh44b312d_0
  - snakemake-minimal=6.1.1=pyhdfd78af_0
  - sqlite=3.34.0=h74cdb3f_0
  - sra-tools=2.8.0=0
  - star=2.7.6a=0
  - tbb=2021.1.1=h4bd325d_0
  - tk=8.6.10=h21135ba_1
  - toposort=1.6=pyhd8ed1ab_0
  - traitlets=5.0.5=py_0
  - typing-extensions=3.7.4.3=0
  - typing_extensions=3.7.4.3=py_0
  - tzdata=2021a=he74cb21_0
  - urllib3=1.26.4=pyhd8ed1ab_0
  - wheel=0.36.2=pyhd3deb0d_0
  - wrapt=1.12.1=py39h3811e60_3
  - xorg-kbproto=1.0.7=h7f98852_1002
  - xorg-libice=1.0.10=h7f98852_0
  - xorg-libsm=1.2.3=hd9c2040_1000
  - xorg-libx11=1.7.0=h7f98852_0
  - xorg-libxau=1.0.9=h7f98852_0
  - xorg-libxdmcp=1.1.3=h7f98852_0
  - xorg-libxext=1.3.4=h7f98852_1
  - xorg-libxpm=3.5.13=h7f98852_0
  - xorg-libxrender=0.9.10=h7f98852_1003
  - xorg-libxt=1.2.1=h7f98852_2
  - xorg-renderproto=0.11.1=h7f98852_1002
  - xorg-xextproto=7.3.0=h7f98852_1002
  - xorg-xproto=7.0.31=h7f98852_1007
  - xz=5.2.5=h516909a_1
  - yaml=0.2.5=h516909a_0
  - yarl=1.6.3=py39h3811e60_1
  - zipp=3.4.1=pyhd8ed1ab_0
  - zlib=1.2.11=h516909a_1010
  - zstd=1.4.9=ha95c52a_0
```

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
