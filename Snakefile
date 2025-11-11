import pandas as pd
samples = pd.read_csv("samples.tsv", sep=",")
RUN_IDS = samples.set_index("name")["run_id"].to_dict()

rule all:
    input:
        expand("aligned_reads/{sample}.bam", sample=RUN_IDS.keys()),
        expand("QC/{sample}_1_fastqc.html",sample = RUN_IDS.keys()),
        expand("QC/{sample}_1_fastqc.zip",sample = RUN_IDS.keys()),
        expand("QC/{sample}_2_fastqc.html",sample = RUN_IDS.keys()),
        expand("QC/{sample}_2_fastqc.zip",sample = RUN_IDS.keys())


rule download_reference_genome:
    output:
        "ref_genome/hg38.fa"
    shell:
        """
        wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.38_GRCh38.p12/GCF_000001405.38_GRCh38.p12_genomic.fna.gz 
        gunzip GCF_000001405.38_GRCh38.p12_genomic.fna.gz
        mv GCF_000001405.38_GRCh38.p12_genomic.fna ref_genome/hg38.fa
        """

rule build_fa_index:
    input:
        "ref_genome/hg38.fa"
    output:
        "ref_genome/hg38.fa.bwt"
    shell:
        "bwa index ref_genome/hg38.fa"

rule download_sample:
    output:
        "data/{sample}_1.fastq",
        "data/{sample}_2.fastq"

    params:
        name = lambda wildcards: wildcards.sample,
        run_id = lambda wildcards: RUN_IDS[wildcards.sample]

    shell:
        """
        cd data
        fastq-dump {params.run_id} --split-files --maxSpotId 10000
        mv {params.run_id}_1.fastq {params.name}_1.fastq
        mv {params.run_id}_2.fastq {params.name}_2.fastq
        """

rule qc:
    input:
        "data/{sample}_1.fastq",
        "data/{sample}_2.fastq"

    output:
        "QC/{sample}_1_fastqc.html",
        "QC/{sample}_2_fastqc.html",
        "QC/{sample}_1_fastqc.zip",
        "QC/{sample}_2_fastqc.zip"

    params:
        run_id = lambda wildcards: wildcards.sample
    shell:
        """
        fastqc data/{params.run_id}_1.fastq data/{params.run_id}_2.fastq -o QC
        """

rule align_sample:
    input:
        "data/{sample}_1.fastq",
        "data/{sample}_2.fastq",
        "ref_genome/hg38.fa.bwt"
    output:
        "aligned_reads/{sample}.bam"
    
    params:
        run_id = lambda wildcards: wildcards.sample
    shell:
        """
        bwa mem ref_genome/hg38.fa data/{params.run_id}_1.fastq data/{params.run_id}_2.fastq -t 8 > aligned_reads/{params.run_id}.bam
        """
