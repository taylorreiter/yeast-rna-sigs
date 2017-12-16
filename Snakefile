import pandas as pd

t = pd.read_csv("outputs/info/transcriptomic.csv.full")
t = t[t['spots'] > 1000000]
INPUTS = t.sort_values(by='size_MB')['Run']

rule all:
    input: expand("outputs/signatures/{SRA_IDS}.sig", SRA_IDS=INPUTS)


rule download_runinfo:
    output: "outputs/info/transcriptomic.csv"
    shell: """
        mkdir -p outputs/info
        wget -O {output}.full 'http://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?save=efetch&db=sra&rettype=runinfo&term=((RNA-seq[All Fields] AND Illumina[All Fields]) AND "Saccharomyces cerevisiae"[Organism] AND "Saccharomyces cerevisiae"[Organism] AND "biomol rna"[Properties]'
        head -n -1 {output}.full > {output}
        # rm {output}.full
    """

rule run_fastq_dump:
    input: 
        bin="bin/fastq-dump"
    output: "outputs/signatures/{SRA_ID}.sig"
    params: SRA_ID="{SRA_ID}"
    # set +o pipefail - explicitly tell snakemake to ignore the pipefail. 
    shell: """
        mkdir -p outputs/signatures
        set +o pipefail
        bin/fastq-dump -A {params.SRA_ID} -Z | head -n 4000000 | trim-low-abund.py --ignore-pairs --variable-coverage -M 2e9 -o - - |sourmash compute -f --track-abundance -o {output} - 
    """
    
rule download_sratoolkit:
    output: "bin/fastq-dump"
    shell: """
        mkdir -p bin
        cd bin
        wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-mac64.tar.gz
        tar xf sratoolkit.current-mac64.tar.gz
        mv sratoolkit.2.8.2-1-mac64/bin/* .
    """   
