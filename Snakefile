import pandas as pd
samples = pd.read_csv("samples.tsv", sep=",")
RUN_IDS = samples.set_index("name")["run_id"].to_dict()


rule all:
    input:
        expand("aligned_reads/{sample}.sam", sample=RUN_IDS.keys()),
        expand("ref_genome/hg38.{idx}.ht2", idx = [1,2,3,4,5,6]),
        expand("QC/{sample}_1_fastqc.html",sample = RUN_IDS.keys()),
        expand("QC/{sample}_1_fastqc.zip",sample = RUN_IDS.keys()),
        expand("QC/{sample}_2_fastqc.html",sample = RUN_IDS.keys()),
        expand("QC/{sample}_2_fastqc.zip",sample = RUN_IDS.keys())



rule download_reference_genome:
    output:
        "ref_genome/hg38.fa"
    shell:
        """
        mkdir -p ref_genome
        wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.38_GRCh38.p12/GCF_000001405.38_GRCh38.p12_genomic.fna.gz 
        gunzip GCF_000001405.38_GRCh38.p12_genomic.fna.gz
        mv GCF_000001405.38_GRCh38.p12_genomic.fna ref_genome/hg38.fa
        
        """

rule build_fa_index:
    input:
        "ref_genome/hg38.fa"
    output:
        expand("ref_genome/hg38.{idx}.ht2", idx = [1,2,3,4,5,6])
    shell:
        """
        hisat2-build ref_genome/hg38.fa ref_genome/hg38
        """

rule align_reads:
    input:
        "data/{sample}_1.fastq",
        "data/{sample}_2.fastq",
        expand("ref_genome/hg38.{idx}.ht2", idx = [1,2,3,4,5,6])
    output:
        "aligned_reads/{sample}.sam"
    shell:
        """
        hisat2 -x ref_genome/hg38 -1 {input[0]} -2 {input[1]} -S {output}
        """
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

