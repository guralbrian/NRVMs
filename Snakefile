configfile: "scripts/setup/fastq_config.json"
import os
import re

READS = ["_1", "_2"]
F_SAMPLES = config["samples_fastq"]

rule all:
    input:
        "data/raw/multiqc/multiqc_report.html",
        "data/raw/anno/Rattus_norvegicus.mRatBN7.2.cdna.all.fa.gz"

rule load_index:
    output: 
        "data/raw/anno/Rattus_norvegicus.mRatBN7.2.cdna.all.fa.gz"
    shell: 
        "wget https://ftp.ensembl.org/pub/release-111/fasta/rattus_norvegicus/cdna/Rattus_norvegicus.mRatBN7.2.cdna.all.fa.gz -P data/raw/anno"
rule load_decoy:
    output: 
        "data/raw/anno/Rattus_norvegicus.mRatBN7.2.dna.toplevel.fa.gz"
    shell: 
        "wget https://ftp.ensembl.org/pub/release-111/fasta/rattus_norvegicus/dna_index/Rattus_norvegicus.mRatBN7.2.dna.toplevel.fa.gz -P data/raw/anno"
rule salmon_index:
    input:
        fa = "data/raw/anno/Rattus_norvegicus.mRatBN7.2.cdna.all.fa.gz",
        #decoys = "data/raw/anno/Rattus_norvegicus.mRatBN7.2.dna.toplevel.fa.gz"
    output:
        directory("data/raw/anno/Rattus_norvegicus.mRatBN7.2.salmon")
    shell:
        "salmon index -t {input.fa} -k 31 --keepFixedFasta -p 16 -i {output}" #-d {input.decoys} 
rule salmon_quant:
    input:
        r1 = "data/raw/fastq/{F_SAMPLES}/{F_SAMPLES}_1.fq.gz",
        r2 = "data/raw/fastq/{F_SAMPLES}/{F_SAMPLES}_2.fq.gz",
        index = "data/raw/anno/Rattus_norvegicus.mRatBN7.2.salmon"
    output:
        "data/raw/fastq/{F_SAMPLES}/quant.sf"
    params:
        dir = "data/raw/fastq/{F_SAMPLES}"
    shell:
        "salmon quant -i {input.index} -l A -p 12 --gcBias "
        "--numGibbsSamples 20 --thinningFactor 100 "
        "-o {params.dir} -1 {input.r1} -2 {input.r2}"
rule fastqc:
    input:
        "data/raw/fastq/{F_SAMPLES}/{F_SAMPLES}{READS}.fq.gz"
    output:
        html = "data/raw/fastq/{F_SAMPLES}/{F_SAMPLES}{READS}_fastqc.html",
        zip = "data/raw/fastq/{F_SAMPLES}/{F_SAMPLES}{READS}_fastqc.zip"
    shell:
        "fastqc {input} --outdir=$(dirname {input})"
rule multiqc:
    input:
        expand(["data/raw/fastq/{sample}/quant.sf",
                "data/raw/fastq/{sample}/{sample}{read}_fastqc.html"],
                sample=F_SAMPLES, read=READS)
    output:
        "data/raw/multiqc/multiqc_report.html"
    shell:
        "multiqc . -o data/raw/multiqc"
