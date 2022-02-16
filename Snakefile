import os

configfile: "config.yaml"

SAMPLES = [f.split('_R1_001.fastq')[0] for f in os.listdir('data') if '_R1_001.fastq' in f]

rule all:
    input:
        expand("results/consensus_pass/{sample}_R1_consensus-pass.fastq", sample=SAMPLES),
        expand("results/consensus_pass/{sample}_R2_consensus-pass.fastq", sample=SAMPLES)

rule filter_seq:
    input:
        "data/{sample}_R1_001.fastq",
        "data/{sample}_R2_001.fastq"
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

rule build_consensus:
    input:
        "results/pair_pass/{sample}_R1_primers-pass_pair-pass.fastq",
        "results/pair_pass/{sample}_R2_primers-pass_pair-pass.fastq"
    output:
        "results/consensus_pass/{sample}_R1_consensus-pass.fastq",
        "results/consensus_pass/{sample}_R2_consensus-pass.fastq",
        "results/logs/{sample}_BC1.log",
        "results/logs/{sample}_BC2.log"
    threads: 16
    run:
        shell("BuildConsensus.py -s {input[0]} --bf BARCODE --pf PRIMER --prcons 0.6 --maxerror 0.1 --maxgap 0.5 --outdir results/consensus_pass/ --outname {wildcards.sample}_R1 --log results/logs/{wildcards.sample}_BC1.log --nproc {threads}"),
        shell("BuildConsensus.py -s {input[1]} --bf BARCODE --pf PRIMER --prcons 0.6 --maxerror 0.1 --maxgap 0.5 --outdir results/consensus_pass/ --outname {wildcards.sample}_R2 --log results/logs/{wildcards.sample}_BC2.log --nproc {threads}")

# rule assemble_pairs:
#     input: 
#         "results/pair_pass/{sample}_R1_primers-pass_pair-pass.fastq",
#         "results/pair_pass/{sample}_R2_primers-pass_pair-pass.fastq"
#     output: 

# AssemblePairs.py align -1 MS12_R2_consensus-pass_pair-pass.fastq \
#     -2 MS12_R1_consensus-pass_pair-pass.fastq --coord presto --rc tail \
#     --1f CONSCOUNT --2f CONSCOUNT PRCONS --outname MS12 --log AP.log