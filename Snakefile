configfile: "scripts/setup/fastq_config.json"
configfile: "scripts/setup/chromo_config.json"
import os
import re

READS = ["_1", "_2"]
F_SAMPLES = config["samples_fastq"]
CHROMO = config["chromo"]

rule all:
    input:
        "data/raw/multiqc/multiqc_report.html",
        "data/processed/deseq/deseq.RDS",
        "results/05_upset_plot/upset_up.png",
        "results/05_upset_plot/upset_down.png"
        "results/06_plot_de/volcano_all.png",
        "results/06_plot_de/pca.png",
        "results/07_clusterProfiler/go_clusters.png"

rule load_transcript:
    output:
        "data/raw/anno/Rattus_norvegicus.mRatBN7.2.cdna.all.fa.gz"
    shell:
        "wget https://ftp.ensembl.org/pub/release-111/fasta/rattus_norvegicus/cdna/Rattus_norvegicus.mRatBN7.2.cdna.all.fa.gz -P data/raw/anno"
rule load_toplevel:
    output:
        "data/raw/anno/Rattus_norvegicus.mRatBN7.2.dna.toplevel.fa.gz"
    shell:
        "wget https://ftp.ensembl.org/pub/release-111/fasta/rattus_norvegicus/dna/Rattus_norvegicus.mRatBN7.2.dna.toplevel.fa.gz -P data/raw/anno/"
rule make_decoy_1:
    input:
        fa="data/raw/anno/Rattus_norvegicus.mRatBN7.2.dna.toplevel.fa.gz"
    output:
        decoy="data/raw/anno/decoy.txt"
    resources:
        mem_mb=10000
    shell:
        """
        grep '^>' <(gunzip -c {input.fa}) | cut -d ' ' -f 1 > {output.decoy}
        """
rule make_decoy_2:
    input:
        "data/raw/anno/decoy.txt"
    output:
        "data/raw/anno/decoy.cleaned.txt"
    shell:
        """
        sed 's/>//g' {input} > {output}
        """
rule make_gentrome:
    input:
        top = "data/raw/anno/Rattus_norvegicus.mRatBN7.2.dna.toplevel.fa.gz",
        transc ="data/raw/anno/Rattus_norvegicus.mRatBN7.2.cdna.all.fa.gz"
    output:
        "data/raw/anno/gentrome.fa.gz"
    shell:
        "cat {input.transc} {input.top} > data/raw/anno/gentrome.fa.gz"
rule salmon_index:
    input:
        fa = "data/raw/anno/gentrome.fa.gz",
        decoy = "data/raw/anno/decoy.cleaned.txt"
    resources:
        mem_mb=32000
    output:
        directory("data/raw/anno/Rattus_norvegicus.mRatBN7.2.salmon")
    shell:
        "salmon index -t {input.fa} -k 31 --keepFixedFasta -p 16 -i {output} -d {input.decoy}"
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
rule tximport:
    input:
        "data/raw/fastq/{F_SAMPLES}/{F_SAMPLES}{READS}.fq.gz"
    output:
        "data/processed/bulk/bulk_ensembl.csv"
    shell:
        "Rscript scripts/01_tximeta.R"
rule phenos:
    output:
        "data/processed/bulk/phenotypes.csv"
    shell:
        "Rscript scripts/02_phenos.R" 
rule ens_to_gene:
    input:
        "data/processed/bulk/bulk_ensembl.csv"
    output:
        "data/processed/bulk/bulk_gene.csv"
    shell:
        "Rscript scripts/03_ens_to_gene.R"
rule deseq:
    input:
        "data/processed/bulk/phenotypes",
        "data/processed/bulk/bulk_gene.csv"
    output:
        "data/processed/deseq/deseq.RDS"
    shell:
        "Rscript scripts/04_deseq.R"
rule upset:
    input:
        "data/processed/deseq/deseq.RDS"
    output:
        "results/05_upset_plot/upset_up.png",
        "results/05_upset_plot/upset_down.png"
    shell:
        "Rscript scripts/05_upset.R"
rule plot_de:
    input:
        "data/processed/deseq/deseq.RDS"
    output:
        "results/06_plot_de/volcano_all.png",
        "results/06_plot_de/pca.png"
    shell:
        "Rscript scripts/06_plot_de.R"
rule clusterProfiler:
    input:
        "data/processed/deseq/deseq.RDS"
    output:
        "results/07_clusterProfiler/go_clusters.png"
    shell:
        "Rscript scripts/07_clusterProfiler.R"
