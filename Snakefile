import os

configfile: "config.yaml"

SAMPLES = [f.split('_R1_001.fastq')[0] for f in os.listdir('data') if '_R1_001.fastq' in f]

rule all:
    input:
        expand("results/annotation_tables/{sample}_FS1_table.tab", sample=SAMPLES),
        expand("results/annotation_tables/{sample}_FS2_table.tab", sample=SAMPLES),
        expand("results/annotation_tables/{sample}_MP1_table.tab", sample=SAMPLES),
        expand("results/annotation_tables/{sample}_MP2_table.tab", sample=SAMPLES),
        expand("results/annotation_tables/{sample}_BC1_table.tab", sample=SAMPLES),
        expand("results/annotation_tables/{sample}_BC2_table.tab", sample=SAMPLES),
        expand("results/annotation_tables/{sample}_AP_table.tab", sample=SAMPLES)

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
        shell("MaskPrimers.py score -s {input[1]} -p {input[3]} --start 0 --mode mask --outdir results/primers_pass/ --outname {wildcards.sample}_R2 --log results/logs/{wildcards.sample}_MP2.log --nproc {threads}")

rule pair_seq:
    input:
        "results/primers_pass/{sample}_R1_primers-pass.fastq",
        "results/primers_pass/{sample}_R2_primers-pass.fastq"
    output:
        temp("results/pair_pass/{sample}_R1_primers-pass_pair-pass.fastq"),
        temp("results/pair_pass/{sample}_R2_primers-pass_pair-pass.fastq")
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
        shell("BuildConsensus.py -s {input[1]} --bf BARCODE --pf PRIMER --maxerror 0.1 --maxgap 0.5 --outdir results/consensus_pass/ --outname {wildcards.sample}_R2 --log results/logs/{wildcards.sample}_BC2.log --nproc {threads}")

rule pair_seq2:
    input:
        "results/consensus_pass/{sample}_R1_consensus-pass.fastq",
        "results/consensus_pass/{sample}_R2_consensus-pass.fastq"
    output:
        temp("results/pair_pass/{sample}_R1_consensus-pass_pair-pass.fastq"),
        temp("results/pair_pass/{sample}_R2_consensus-pass_pair-pass.fastq")
    threads: 1
    shell:
        "PairSeq.py -1 {input[0]} -2 {input[1]} --coord presto --outdir results/pair_pass/"

rule assemble_pairs:
    input:
        "results/pair_pass/{sample}_R1_consensus-pass_pair-pass.fastq",
        "results/pair_pass/{sample}_R2_consensus-pass_pair-pass.fastq"
    output:
        "results/assemble_pass/{sample}_assemble-pass.fastq",
        "results/logs/{sample}_AP.log"
    threads: 16
    shell:
        "AssemblePairs.py align -1 {input[0]} -2 {input[1]} --coord presto --rc tail --1f CONSCOUNT --2f CONSCOUNT PRCONS --outname {wildcards.sample} --outdir results/assemble_pass --log results/logs/{wildcards.sample}_AP.log --nproc {threads}"    

rule parse_headers1:
    input:
        "results/assemble_pass/{sample}_assemble-pass.fastq"
    output:
        "results/assemble_pass/{sample}_assemble-pass_reheader.fastq"
    threads: 1
    shell:
        "ParseHeaders.py collapse -s {input} -f CONSCOUNT --act min --outdir results/assemble_pass"

rule collapse_seq:
    input:
        "results/assemble_pass/{sample}_assemble-pass_reheader.fastq"
    output:
        "results/collapse_unique/{sample}_collapse-unique.fastq"
    threads: 1
    shell:
        "CollapseSeq.py -s {input} -n 20 --inner --uf PRCONS --cf CONSCOUNT --act sum --outname {wildcards.sample} --outdir results/collapse_unique"

rule split_seq:
    input:
        "results/collapse_unique/{sample}_collapse-unique.fastq"
    output:
        "results/split_seq/{sample}_atleast-2.fastq"
    threads: 1
    shell:
        "SplitSeq.py group -s {input} -f CONSCOUNT --num 2 --outname {wildcards.sample} --outdir results/split_seq"

rule parse_headers2:
    input:
        "results/split_seq/{sample}_atleast-2.fastq"
    output:
        "results/annotation_tables/{sample}_atleast-2_headers.tab"
    threads: 1
    shell:
        "ParseHeaders.py table -s {input} -f ID PRCONS CONSCOUNT DUPCOUNT --outdir results/annotation_tables"
    
rule parse_log:
    input:
        "results/logs/{sample}_FS1.log",
        "results/logs/{sample}_FS2.log",
        "results/logs/{sample}_MP1.log",
        "results/logs/{sample}_MP2.log",
        "results/logs/{sample}_BC1.log",
        "results/logs/{sample}_BC2.log",
        "results/logs/{sample}_AP.log",
        "results/annotation_tables/{sample}_atleast-2_headers.tab"
    output:
        "results/annotation_tables/{sample}_FS1_table.tab",
        "results/annotation_tables/{sample}_FS2_table.tab",
        "results/annotation_tables/{sample}_MP1_table.tab",
        "results/annotation_tables/{sample}_MP2_table.tab",
        "results/annotation_tables/{sample}_BC1_table.tab",
        "results/annotation_tables/{sample}_BC2_table.tab",
        "results/annotation_tables/{sample}_AP_table.tab",
    threads: 4
    run:
        shell('ParseLog.py -l {input[0]} {input[1]} -f ID QUALITY --outdir results/annotation_tables'),
        shell('ParseLog.py -l {input[2]} {input[3]} -f ID PRIMER BARCODE ERROR --outdir results/annotation_tables'),
        shell('ParseLog.py -l {input[4]} {input[5]} -f BARCODE SEQCOUNT CONSCOUNT PRIMER PRCONS PRCOUNT PRFREQ ERROR --outdir results/annotation_tables'),
        shell('ParseLog.py -l {input[6]} -f ID LENGTH OVERLAP ERROR PVALUE FIELDS1 FIELDS2 --outdir results/annotation_tables')