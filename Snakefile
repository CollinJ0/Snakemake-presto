import os

configfile: "config.yaml"

SAMPLES = [f.split('_R1_001.fastq.gz')[0] for f in os.listdir('data') if '_R1_001.fastq.gz' in f]

rule all:
    input:
        expand("results/pair_pass/{sample}_R1_primers-pass_pair-pass.fastq", sample=SAMPLES),
        expand("results/pair_pass/{sample}_R2_primers-pass_pair-pass.fastq", sample=SAMPLES)

rule filter_seq:
    input:
        "data/{sample}_R1_001.fastq.gz",
        "data/{sample}_R2_001.fastq.gz"
    output:
        temp("results/quality_pass/{sample}_R1_quality-pass.fastq"),
        temp("results/quality_pass/{sample}_R2_quality-pass.fastq"),
        "results/logs/{sample}_FS1.log",
        "results/logs/{sample}_FS2.log"
    threads: 16
    run:
        shell("FilterSeq.py quality -s {input[0]} -q 20 --outdir results/quality_pass/ --outname {wildcards.sample}_R1 --log results/logs/{wildcards.sample}_FS1.log --nproc {threads}"),
        shell("FilterSeq.py quality -s {input[1]} -q 20 --outdir results/quality_pass/ --outname {wildcards.sample}_R2 --log results/logs/{wildcards.sample}_FS2.log --nproc {threads}")

rule mask_primers:
    input:
        "results/quality_pass/{sample}_R1_quality-pass.fastq",
        "results/quality_pass/{sample}_R2_quality-pass.fastq",
        config['constant_primers'],
        config['vgene_primers']
    output:
        temp("results/primers_pass/{sample}_R1_primers-pass.fastq"),
        temp("results/primers_pass/{sample}_R2_primers-pass.fastq"),
        "results/logs/{sample}_MP1.log",
        "results/logs/{sample}_MP2.log"
    threads: 16
    run:
        shell("MaskPrimers.py score -s {input[0]} -p {input[2]} --start 12 --mode cut --barcode --outdir results/primers_pass/ --outname {wildcards.sample}_R1 --log results/logs/{wildcards.sample}_MP1.log --nproc {threads}"),
        shell("MaskPrimers.py score -s {input[1]} -p {input[3]} --start 0 --mode cut --barcode --outdir results/primers_pass/ --outname {wildcards.sample}_R2 --log results/logs/{wildcards.sample}_MP2.log --nproc {threads}")

rule pair_seq:
    input:
        "results/primers_pass/{sample}_R1_primers-pass.fastq",
        "results/primers_pass/{sample}_R2_primers-pass.fastq"
    output:
        "results/pair_pass/{sample}_R1_primers-pass_pair-pass.fastq",
        "results/pair_pass/{sample}_R2_primers-pass_pair-pass.fastq"
    threads: 1
    shell:
        "PairSeq.py -1 {input[0]} -2 {input[1]} --1f BARCODE --coord illumina --outdir results/pair_pass/"
