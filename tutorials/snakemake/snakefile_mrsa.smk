from snakemake.utils import min_version
min_version("8.0.0")

rule all:
    """
    Collect the main outputs of the workflow.
    """
    input:
        "results/tables/NCTC8325.counts.tsv",
        "results/multiqc/NCTC8325.multiqc.html"



def get_sample_url(wildcards):
    samples = config["samples"]
    return samples[wildcards.sample_id]


rule get_SRA_by_accession:
    """
    Retrieve a single-read FASTQ file
    """
    output:
        "data/{sample_id}.fastq.gz"
    params:
        url = get_sample_url,
        max_reads = config["max_reads"]
    shell:
        """
        curl -L -A "Mozilla/5.0" {params.url} | seqtk sample - {params.max_reads} | gzip -c > {output[0]}
        """

rule fastqc:
    """
    Run FastQC on a FASTQ file.
    """
    shadow: "minimal"
    output:
        "results/fastqc/{sample_id}_fastqc.html",
        "results/fastqc/{sample_id}_fastqc.zip"
    input:
        "data/{sample_id}.fastq.gz"
    shell:
        """
        # Run fastQC and save the output to the current directory
        fastqc {input} -q -o .

        # Move the files which are used in the workflow
        mv {wildcards.sample_id}_fastqc.html {output[0]}
        mv {wildcards.sample_id}_fastqc.zip {output[1]}
        """

rule multiqc:
    """
    Aggregate all FastQC reports into a MultiQC report.
    """
    output:
        html = "results/multiqc/{genome_id}.multiqc.html",
        stats = "results/multiqc/{genome_id}.multiqc_general_stats.txt"
    input:
        bams = expand("results/fastqc/{sample_id}_fastqc.zip", sample_id = config["samples"].keys())
    log:
        "results/logs/multiqc/{genome_id}.log"
    shell:
        """
        # Run multiQC and keep the html report
        multiqc -n multiqc.html {input}
        mv multiqc.html {output.html}
        mv multiqc_data/multiqc_general_stats.txt {output.stats}

        # Remove the other directory that multiQC creates
        rm -rf multiqc_data
        """

def basename_noext_outfile(output):
    outfile = output[0]
    outfile = outfile.removesuffix('.gz')
    outfile = os.path.splitext(outfile)[0]
    outfile = os.path.basename(outfile)
    return outfile


def get_genome_url(wildcards):
    return config["genomes"][wildcards.genome_id]["fasta"]

rule get_genome_fasta:
    """
    Retrieve the sequence in fasta format for a genome.
    """
    output:
        "data/ref/{genome_id}.fa.gz"
    params:
        genome_url = get_genome_url
        #genome_url = config[{wildcards.genome_id}]["fasta"]
        #of_noext = basename_noext_outfile
    log:
        "results/logs/get_genome_fasta/{genome_id}.log"
    shell:
        """
        wget -o {log} {params.genome_url} -O {output}
        """

rule get_genome_gff3:
    """
    Retrieve annotation in gff3 format for a genome.
    """
    output:
        "data/ref/{genome_id}.gff3.gz"
    log:
        "results/logs/get_genome_gff3/{genome_id}.log"
    shell:
        """
        wget -o {log} ftp://ftp.ensemblgenomes.org/pub/bacteria/release-37/gff3/bacteria_18_collection/staphylococcus_aureus_subsp_aureus_nctc_8325//Staphylococcus_aureus_subsp_aureus_nctc_8325.ASM1342v1.37.gff3.gz -O {output}
        """

rule index_genome:
    """
    Index a genome using Bowtie 2.
    """
    shadow: "minimal"
    output:        
        expand("results/bowtie2/{{genome_id}}.{substr}.bt2", 
            substr = ["1", "2", "3", "4", "rev.1", "rev.2"])
    input:
        "data/ref/{genome_id}.fa.gz"
    log:
        "results/logs/index_genome/{genome_id}.log"
    shell:
        """
        # Bowtie2 cannot use .gz, so unzip to a temporary file first
        gunzip -c {input} > tempfile
        bowtie2-build tempfile results/bowtie2/{genome_id} >{log}

        # Remove the temporary file
        rm tempfile
        """

rule align_to_genome:
    """
    Align a fastq file to a genome index using Bowtie 2.
    """
    output:
        temp("results/bam/{sample_id,\\w+}.bam")
    input:
        "data/{sample_id}.fastq.gz",
        expand("results/bowtie2/{{genome_id}}.{substr}.bt2", 
            substr = ["1", "2", "3", "4", "rev.1", "rev.2"])
    shell:
        """
        bowtie2 -x results/bowtie2/{genome_id} -U {input[0]} > {output}
        """


rule sort_bam:
    """
    Sort a bam file.
    """
    output:
        "results/bam/{sample_id}.sorted.bam"
    input:
        "results/bam/{sample_id}.bam"
    shell:
        """
        samtools sort {input} > {output}
        """

rule generate_count_table:
    """
    Generate a count table using featureCounts.
    """
    output:
        "results/tables/{genome_id}.counts.tsv"
    input:
        bams = expand("results/bam/{sample_id}.sorted.bam", sample_id = config["samples"].keys()),
        annotation = "data/ref/{genome_id}.gff3.gz"
    log:
        "results/logs/generate_count_table/{genome_id}.log"
    shell:
        """
        featureCounts -t gene -g gene_id -a {input.annotation} -o {output} {input.bams} 2>{log}
        """