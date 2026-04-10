# RD-and-Rh30-Test

Pipeline used to analyze RD and Rh30 bulk RNA-sequencing data to validate the pipeline outputs coorelate with what is seen in corresponding mass spectrometry data.

---

## Software Setup

You will need to install Miniforge into our data directory using the method described in the [Conda on Biowulf](https://hpc.nih.gov/docs/diy_installation/conda.html) documentation.


### Change Working Directory 

```bash
cd /data/mckeeka
pwd
```

### Install MiniConda 

```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
#Type "yes" to accept the license. When asked for the installation location, change the directory to "/data/mckeeka/miniconda3". Type "yes" when asked to initialize MiniConda.
```

### Activate Changes and Source bashrc File

```bash
source ~/.bashrc
conda --version
```

## Upload Data and Reference Into Directory

### Make Working Directory

```bash
mkdir bulkRNA_RMS
cd /data/mckeeka/bulkRNA_RMS
```

### Download and Reformat FASTQ Files

RD and Rh30 FASTQ files were transferred to a folder called "bulkRNA" in the bulkRNA_RMS directory using Globus.
The following code changes the file names to fit the pipeline format (<sample>.fastq.<read>.gz).

UUID for Globus Transfer: 25ae9b88-ec49-4de3-803a-1c02546cee80

```bash
cd /data/mckeeka/bulkRNA_RMS/bulkRNA
for f in *_S*_R*_001.fastq.gz; do
base=$(basename "$f")

Sample=$(echo "$base" | awk -F'_' '{print $2"_"$3}')
Read=$(echo "$base" | awk -F'_' '{print $5}')

mv "$f" "${Sample}.fastq.${Read}.gz"
done
```

### Download Human Reference Transcriptome and GTF

A copy of the reference folder from bulkRNA_sarcoma was symlinked into the bulkRNA_RMS directory. This folder contains the human reference transcriptome (Ensembl release 115, GRCh38).

```bash
cd /data/mckeeka/bulkRNA_RMS/
ln -s /data/mckeeka/bulkRNA_sarcoma/reference/
```

## Create Raw QC Pipeline Working Directory

```bash
cd /data/mckeeka/bulkRNA_RMS/
mkdir run_bulkRNA
cd /data/mckeeka/bulkRNA_RMS/run_bulkRNA
mkdir rawQC
cd /data/mckeeka/bulkRNA_RMS/run_bulkRNA/rawQC
mkdir fastqc
cd /data/mckeeka/bulkRNA_RMS/
```

## Generate Raw QC Pipeline Configuration

This pipeline was generated to perform analysis of the raw FASTQ data after sequencing.

### Install QC Tools

```bash
cd /data/mckeeka/bulkRNA_RMS/
conda create -n rawQC -c bioconda snakemake fastqc multiqc -y
conda activate rawQC
```

### Create Snakemake Raw QC Configuration File

```bash
nano rawQC_pipeline.smk

# Add the following code to the configuration file:

SAMPLES = glob_wildcards("bulkRNA/{sample}.fastq.R1.gz").sample
READS = ["R1", "R2"]

rule all:
    input:
        expand("run_bulkRNA/rawQC/fastqc/{sample}.fastq.{read}_fastqc.zip", sample=SAMPLES, read=READS),
        "run_bulkRNA/rawQC/multiqc_report.html"

rule fastqc:
  input:
    "bulkRNA/{sample}.fastq.{read}.gz"
  output:
    html="run_bulkRNA/rawQC/fastqc/{sample}.fastq.{read}_fastqc.html",
    zip="run_bulkRNA/rawQC/fastqc/{sample}.fastq.{read}_fastqc.zip"
  shell:
    """
    fastqc -o run_bulkRNA/rawQC/fastqc {input}
    """

rule multiqc:
  input:
    expand("run_bulkRNA/rawQC/fastqc/{sample}.fastq.{read}_fastqc.zip",
            sample=SAMPLES,
            read=READS)
  output:
    "run_bulkRNA/rawQC/multiqc_report.html"
  shell:
    """
    multiqc run_bulkRNA/rawQC/fastqc -o run_bulkRNA/rawQC
    """
```

### Run Raw QC Configuration File

The pipeline must be run using sbatch on the Biowulf cluster.

```bash
cd /data/mckeeka/bulkRNA_RMS/
sbatch --time=00-04:00:00 --wrap="snakemake -s rawQC_pipeline.smk"
```

## Create CutAdapt Pipeline Working Directory

The CutAdapt pipeline requires a working directory where the FASTQ files can be accessed.

```bash
cd /data/mckeeka/bulkRNA_RMS/run_bulkRNA
mkdir trimmed_FASTQ
mkdir logs
cd /data/mckeeka/bulkRNA_RMS/run_bulkRNA/logs
mkdir logs_CutAdapt
cd /data/mckeeka/bulkRNA_RMS
```

## Generate CutAdapt Pipeline Configuration

This pipeline was generated to cut the adapters from the raw FASTQ files after sequencing.

### Install CutAdapt Tools

```bash
cd /data/mckeeka/bulkRNA_RMS
conda create -n CutAdapt -c bioconda snakemake cutadapt -y
conda activate CutAdapt
```

### Create Snakemake CutAdapt Configuration File

```bash
nano CutAdapt_pipeline.smk

# Add the following code to the configuration file:

adapter = "AGATCGGAAGAG"         #Illumina Universal Adapter
nextera = "CTGTCTCTTATACACATCT"  #Nextera transposase
minimum_length = 15              #Decreased from the recommended 20 since I am interested in smORFs
quality_trimming = "20,20"       #Recommended value
overlap = 5                      #Recommended value
threads = 4 

SAMPLES = glob_wildcards("bulkRNA/{sample}.fastq.R1.gz").sample

rule all:
    input:
        expand("run_bulkRNA/trimmed_FASTQ/{sample}.fastq.R1.trimmed.gz", sample=SAMPLES),
        expand("run_bulkRNA/trimmed_FASTQ/{sample}.fastq.R2.trimmed.gz", sample=SAMPLES)

rule cutadapt_pe:
  input:
    r1 = "bulkRNA/{sample}.fastq.R1.gz",
    r2 = "bulkRNA/{sample}.fastq.R2.gz"
  output:
    r1 = "run_bulkRNA/trimmed_FASTQ/{sample}.fastq.R1.trimmed.gz",
    r2 = "run_bulkRNA/trimmed_FASTQ/{sample}.fastq.R2.trimmed.gz"
  log:
    "run_bulkRNA/logs/logs_CutAdapt/{sample}.CutAdapt.log"
  shell:
    """
    cutadapt \
        -a {adapter} \
        -A {adapter} \
        -a {nextera} \
        -A {nextera} \
        -a "G{{10}}" \
        -A "G{{10}}" \
        -m {minimum_length} \
        -q {quality_trimming} \
        -O {overlap} \
        --pair-filter=any \
        --cores {threads} \
        -o {output.r1} \
        -p {output.r2} \
        {input.r1} {input.r2} > {log} 2>&1
    """

```

### Run CutAdapt Configuration File

The pipeline must be run using sbatch on the Biowulf cluster.

```bash
cd /data/mckeeka/bulkRNA_RMS
sbatch --cpus-per-task=4 --mem=16G --time=04:00:00 \--wrap "snakemake -s CutAdapt_pipeline.smk -j 4"
```

## Create Trimmed QC Pipeline Working Directory

The Trimmed QC pipeline requires a working directory where the FASTQ files can be accessed.

```bash
cd /data/mckeeka/bulkRNA_RMS/run_bulkRNA
mkdir trimmedQC
cd /data/mckeeka/bulkRNA_RMS/run_bulkRNA/trimmedQC
mkdir fastqc
```

## Generate Trimmed QC Pipeline Configuration

This pipeline was generated to perform analysis of the raw FASTQ data after sequencing.

### Install QC Tools

```bash
cd /data/mckeeka/bulkRNA_RMS/
conda create -n trimmedQC -c bioconda snakemake fastqc multiqc -y
conda activate trimmedQC
```

### Create Snakemake Trimmed QC Configuration File

```bash
nano trimmedQC_pipeline.smk

# Add the following code to the configuration file:

SAMPLES = glob_wildcards("bulkRNA/{sample}.fastq.R1.gz").sample
READS = ["R1", "R2"]

rule all:
    input:
        expand("run_bulkRNA/trimmedQC/fastqc/{sample}.fastq.{read}.trimmed_fastqc.zip", sample=SAMPLES, read=READS),
        "run_bulkRNA/trimmedQC/trimmedQC_multiqc_report.html"

rule fastqc:
  input:
    "run_bulkRNA/trimmed_FASTQ/{sample}.fastq.{read}.trimmed.gz"
  output:
    html="run_bulkRNA/trimmedQC/fastqc/{sample}.fastq.{read}.trimmed_fastqc.html",
    zip="run_bulkRNA/trimmedQC/fastqc/{sample}.fastq.{read}.trimmed_fastqc.zip"
  shell:
    """
    fastqc -o run_bulkRNA/trimmedQC/fastqc {input}
    """

rule multiqc:
  input:
    expand("run_bulkRNA/trimmedQC/fastqc/{sample}.fastq.{read}.trimmed_fastqc.zip",
            sample=SAMPLES,
            read=READS)
  output:
    "run_bulkRNA/trimmedQC/multiqc_report.html"
  shell:
    """
    multiqc run_bulkRNA/trimmedQC/fastqc -o run_bulkRNA/trimmedQC
    """

```

### Run Trimmed QC Configuration File

The pipeline must be run using sbatch on the Biowulf cluster.

```bash
cd /data/mckeeka/bulkRNA_RMS/
sbatch --time=02:00:00 --wrap "snakemake -s trimmedQC_pipeline.smk"
```


## Generate Clean FASTQ Pipeline Configuration

This pipeline was generated to clean the FASTQs after sequencing to eliminate contamination.

### Install Clean FASTQ Tools

```bash
cd /data/mckeeka/bulkRNA_RMS
conda create -n cleanFASTQ -c bioconda snakemake kraken2 bowtie2 KrakenTools BEDTools SAMtools -y
conda activate cleanFASTQ
```

### Create Standard Databases

The Clean FASTQ pipeline requires a working directory where the standard reference databases can be accessed. You can symlink these files instead of copying them into the pipeline directory to prevent the duplication of large data files in your directory.

```bash
cd /data/mckeeka/bulkRNA_RMS/run_bulkRNA
mkdir clean_FASTQ
cd /data/mckeeka/bulkRNA_RMS/run_bulkRNA/clean_FASTQ
mkdir kraken2_output
mkdir bowtie2_output
cd /data/mckeeka/bulkRNA_RMS/run_bulkRNA/logs
mkdir logs_cleanFASTQ
cd /data/mckeeka/bulkRNA_RMS/run_bulkRNA/logs/logs_cleanFASTQ
mkdir logs_kraken2
mkdir logs_bowtie2

cd /data/mckeeka/bulkRNA_sarcoma/reference

#Copy Kraken2 Standard Reference Database from Biowulf
cp -r /fdb/kraken/20220803_standard_kraken2 kraken2_database

#Create Contaminant Reference for Bowtie2
grep -E 'gene_biotype "(artifact|Mt_rRNA|Mt_tRNA|ribozyme|rRNA|rRNA_pseudogene|scaRNA|snoRNA|snRNA|vault_RNA)"' Homo_sapiens.GRCh38.115.gtf > contaminants.gtf

awk '$3=="exon"' contaminants.gtf | \
awk '{print$1"\t"$4-1"\t"$5}' > contaminants.bed

bedtools getfasta -fi Homo_sapiens.GRCh38.dna.primary_assembly.fa -bed contaminants.bed -fo contaminants.fa

mkdir -p contaminants_index
bowtie2-build contaminants.fa contaminants_index/contaminants
```

### Create Snakemake Clean FASTQ Configuration File

```bash
cd /data/mckeeka/bulkRNA_RMS
nano cleanFASTQ_pipeline.smk

# Add the following code to the configuration file:

SAMPLES = glob_wildcards("run_bulkRNA/trimmed_FASTQ/{sample}.fastq.R1.trimmed.gz").sample

rule all:
    input:
        expand("run_bulkRNA/clean_FASTQ/{sample}.fastq.R{read}.clean.gz", sample=SAMPLES, read=[1,2])

rule kraken2:
  input:
    r1 = "run_bulkRNA/trimmed_FASTQ/{sample}.fastq.R1.trimmed.gz",
    r2 = "run_bulkRNA/trimmed_FASTQ/{sample}.fastq.R2.trimmed.gz"
  output:
    report = "run_bulkRNA/clean_FASTQ/kraken2_output/{sample}_report.txt",
    output = "run_bulkRNA/clean_FASTQ/kraken2_output/{sample}_output.txt"
  threads: 4
  log:
    "run_bulkRNA/logs/logs_cleanFASTQ/logs_kraken2/{sample}.kraken2.log"
  shell:
    """
    kraken2 \
        --db reference/kraken2_database \
        --paired {input.r1} {input.r2} \
        --report {output.report} \
        --output {output.output} \
        --threads {threads} \
        &> {log}
    """

rule extract_human_unclassified:
  input:
    kraken2 = "run_bulkRNA/clean_FASTQ/kraken2_output/{sample}_output.txt",
    report = "run_bulkRNA/clean_FASTQ/kraken2_output/{sample}_report.txt",
    r1 = "run_bulkRNA/trimmed_FASTQ/{sample}.fastq.R1.trimmed.gz",
    r2 = "run_bulkRNA/trimmed_FASTQ/{sample}.fastq.R2.trimmed.gz"
  output:
    r1 = "run_bulkRNA/clean_FASTQ/kraken2_output/{sample}.fastq.R1.kraken.gz",
    r2 = "run_bulkRNA/clean_FASTQ/kraken2_output/{sample}.fastq.R2.kraken.gz"
  log:
    "run_bulkRNA/logs/logs_cleanFASTQ/logs_kraken2/{sample}.kraken2_filter.log"
  shell:
    """
    tmp_h1=$(mktemp)
    tmp_h2=$(mktemp)
    tmp_u1=$(mktemp)
    tmp_u2=$(mktemp)

    {{
        #Extract human reads
        extract_kraken_reads.py \
            -k {input.kraken2} \
            -r {input.report} \
            -s {input.r1} \
            -s2 {input.r2} \
            -t 9606 \
            --include-children \
            --fastq-output \
            -o $tmp_h1 \
            -o2 $tmp_h2

        #Extract unclassified reads
        extract_kraken_reads.py \
            -k {input.kraken2} \
            -r {input.report} \
            -s {input.r1} \
            -s2 {input.r2} \
            -t 0 \
            --fastq-output \
            -o $tmp_u1 \
            -o2 $tmp_u2

        #Combine human and unclassified reads
        cat "$tmp_h1" "$tmp_u1" | gzip -c > {output.r1}
        cat "$tmp_h2" "$tmp_u2" | gzip -c > {output.r2}

        #Remove temporary files
        rm -f "$tmp_h1" "$tmp_h2" "$tmp_u1" "$tmp_u2"
    }} &> {log}
    """

rule bowtie2_contaminant_mapping:
  input:
    r1 = "run_bulkRNA/clean_FASTQ/kraken2_output/{sample}.fastq.R1.kraken.gz",
    r2 = "run_bulkRNA/clean_FASTQ/kraken2_output/{sample}.fastq.R2.kraken.gz"
  output:
    bam = "run_bulkRNA/clean_FASTQ/bowtie2_output/{sample}_contamination.bam"
  threads: 4
  log:
    "run_bulkRNA/logs/logs_cleanFASTQ/logs_bowtie2/{sample}.bowtie2.log"
  shell:
    """
    bowtie2 \
        -x reference/contaminants_index/contaminants \
        -1 {input.r1} \
        -2 {input.r2} \
        --sensitive \
        --threads {threads} \
        2> {log} \
        | samtools view -b -o {output.bam} -
    """

rule filter_unmapped:
  input:
    bam = "run_bulkRNA/clean_FASTQ/bowtie2_output/{sample}_contamination.bam"
  output:
    r1 = "run_bulkRNA/clean_FASTQ/{sample}.fastq.R1.clean.gz",
    r2 = "run_bulkRNA/clean_FASTQ/{sample}.fastq.R2.clean.gz"
  log:
    "run_bulkRNA/logs/logs_cleanFASTQ/logs_bowtie2/{sample}.bowtie2_filter.log"
  shell:
    """
    set -euo pipefail

    tmp_bam=$(mktemp --suffix=.bam)
    tmp_namesort=$(mktemp --suffix=.bam)
    tmp_r1=$(mktemp --suffix=.fq)
    tmp_r2=$(mktemp --suffix=.fq)

    samtools view -b -f 12 -F 256 {input.bam} > "$tmp_bam"
    samtools sort -n -o "$tmp_namesort" "$tmp_bam"

    bedtools bamtofastq \
        -i "$tmp_namesort" \
        -fq "$tmp_r1" \
        -fq2 "$tmp_r2"

    gzip -c "$tmp_r1" > {output.r1}
    gzip -c "$tmp_r2" > {output.r2}

    rm -f "$tmp_bam" "$tmp_namesort" "$tmp_r1" "$tmp_r2"
    """ + "&> {log}"
```

### Run Clean FASTQ Configuration File

The pipeline must be run using sbatch on the Biowulf cluster.

```bash
cd /data/mckeeka/bulkRNA_RMS
sbatch --cpus-per-task=8 --mem=128G --time=01-00:00:00 \--wrap "snakemake -s cleanFASTQ_pipeline.smk -j 8"
```

## Create Clean QC Pipeline Working Directory

The Clean QC pipeline requires a working directory where the FASTQ files can be accessed.

```bash
cd /data/mckeeka/bulkRNA_RMS/run_bulkRNA
mkdir cleanQC
cd /data/mckeeka/bulkRNA_RMS/run_bulkRNA/cleanQC
mkdir fastqc
```

## Generate Clean QC Pipeline Configuration

This pipeline was generated to perform analysis of the raw FASTQ data after filtering.

### Install QC Tools

```bash
cd /data/mckeeka/bulkRNA_RMS/
conda create -n cleanQC -c conda-forge -c bioconda snakemake python=3.11 fastqc multiqc -y
conda activate cleanQC
```

### Create Clean QC Configuration File

```bash
nano cleanQC_pipeline.smk

# Add the following code to the configuration file:

SAMPLES = glob_wildcards("bulkRNA/{sample}.fastq.R1.gz").sample
READS = ["R1", "R2"]

rule all:
    input:
        expand("run_bulkRNA/cleanQC/fastqc/{sample}.fastq.{read}.clean_fastqc.zip", sample=SAMPLES, read=READS),
        "run_bulkRNA/cleanQC/multiqc_report.html"

rule fastqc:
  input:
    "run_bulkRNA/clean_FASTQ/{sample}.fastq.{read}.clean.gz"
  output:
    html="run_bulkRNA/cleanQC/fastqc/{sample}.fastq.{read}.clean_fastqc.html",
    zip="run_bulkRNA/cleanQC/fastqc/{sample}.fastq.{read}.clean_fastqc.zip"
  shell:
    """
    fastqc -o run_bulkRNA/cleanQC/fastqc {input}
    """

rule multiqc:
  input:
    expand("run_bulkRNA/cleanQC/fastqc/{sample}.fastq.{read}.clean_fastqc.zip",
            sample=SAMPLES,
            read=READS)
  output:
    "run_bulkRNA/cleanQC/multiqc_report.html"
  shell:
    """
    multiqc run_bulkRNA/cleanQC/fastqc -o run_bulkRNA/cleanQC
    """

```

### Run Trimmed QC Configuration File

The pipeline must be run using sbatch on the Biowulf cluster.

```bash
cd /data/mckeeka/bulkRNA_RMS/
sbatch --time=02:00:00 --wrap "snakemake -s cleanQC_pipeline.smk --cores 4"
```


## Create Read Counts QC Pipeline Working Directory

```bash
cd /data/mckeeka/bulkRNA_RMS/run_bulkRNA/clean_FASTQ/
mkdir ReadCounts_output
```

## Generate Read Counts Pipeline Configuration

This pipeline was generated to count the reads that were filtered out of each step of the Clean FASTQ pipeline.

### Install Read Counts Tools

```bash
cd /data/mckeeka/bulkRNA_RMS
conda create -n ReadCounts -c conda-forge -c bioconda snakemake python=3.10 pandas -y
conda activate ReadCounts
```

### Create Snakemake Read Counts Configuration File

```bash
nano ReadCounts_pipeline.smk

# Add the following code to the configuration file:

SAMPLES = sorted(set(glob_wildcards(
    "run_bulkRNA/trimmed_FASTQ/{sample}.fastq.R{read}.trimmed.gz"
).sample))

rule all:
    input:
        "run_bulkRNA/clean_FASTQ/ReadCounts_output/read_summary.csv",
        expand("run_bulkRNA/trimmed_FASTQ/{sample}.fastq.R{read}.count.txt", sample=SAMPLES, read=[1,2]),
        expand("run_bulkRNA/clean_FASTQ/kraken2_output/{sample}.fastq.R{read}.kraken.count.txt", sample=SAMPLES, read=[1,2]),
        expand("run_bulkRNA/clean_FASTQ/{sample}.fastq.R{read}.clean.count.txt", sample=SAMPLES, read=[1,2])

rule count_trimmed_reads:
  input:
    "run_bulkRNA/trimmed_FASTQ/{sample}.fastq.R{read}.trimmed.gz"
  output:
    "run_bulkRNA/trimmed_FASTQ/{sample}.fastq.R{read}.count.txt"
  shell:
    """
    cat {input} | wc -l | awk '{{print $1/4}}' > {output}
    """

rule count_kraken_reads:
  input:
    "run_bulkRNA/clean_FASTQ/kraken2_output/{sample}.fastq.R{read}.kraken.gz"
  output:
    "run_bulkRNA/clean_FASTQ/kraken2_output/{sample}.fastq.R{read}.kraken.count.txt"
  shell:
    """
    cat {input} | wc -l | awk '{{print $1/4}}' > {output}
    """

rule count_bowtie_reads:
  input:
    "run_bulkRNA/clean_FASTQ/{sample}.fastq.R{read}.clean.gz"
  output:
    "run_bulkRNA/clean_FASTQ/{sample}.fastq.R{read}.clean.count.txt"
  shell:
    """
    cat {input} | wc -l | awk '{{print $1/4}}' > {output}
    """

rule summarize_read_counts:
  input:
      TRIMMED = expand("run_bulkRNA/trimmed_FASTQ/{sample}.fastq.R{read}.count.txt", sample=SAMPLES, read=[1,2]),
      KRAKEN = expand("run_bulkRNA/clean_FASTQ/kraken2_output/{sample}.fastq.R{read}.kraken.count.txt", sample=SAMPLES, read=[1,2]),
      BOWTIE = expand("run_bulkRNA/clean_FASTQ/{sample}.fastq.R{read}.clean.count.txt", sample=SAMPLES, read=[1,2])
  output:
      "run_bulkRNA/clean_FASTQ/ReadCounts_output/read_summary.csv"
  run:
      import pandas as pd
      data = []
      for sample in SAMPLES:
        row = {"sample": sample}
        for step, path_template in [
            ("TRIMMED", "run_bulkRNA/trimmed_FASTQ/{sample}.fastq.R{read}.count.txt"),
            ("KRAKEN", "run_bulkRNA/clean_FASTQ/kraken2_output/{sample}.fastq.R{read}.kraken.count.txt"),
            ("BOWTIE", "run_bulkRNA/clean_FASTQ/{sample}.fastq.R{read}.clean.count.txt")
        ]:
            for read in [1,2]:
              file_path = path_template.format(sample=sample, read=read)
              row[f"{step}_{read}"] = int(float(open(file_path).read().strip()))
        data.append(row)
      df = pd.DataFrame(data)
      df.to_csv(output[0], index=False)

```

### Run Read Counts Configuration File

The pipeline must be run using sbatch on the Biowulf cluster.

```bash
cd /data/mckeeka/bulkRNA_RMS
sbatch --cpus-per-task=4 --mem=16G --time=01:00:00 \--wrap "snakemake -s ReadCounts_pipeline.smk -j 4"
```


## Create STAR Mapping Pipeline Working Directory

The STAR Mapping pipeline requires a working directory where the FASTQ files can be accessed.

```bash
cd /data/mckeeka/bulkRNA_RMS/run_bulkRNA
mkdir STAR
cd /data/mckeeka/bulkRNA_RMS/run_bulkRNA/logs
mkdir logs_STAR
cd /data/mckeeka/bulkRNA_RMS
```

## Generate STAR Mapping Pipeline Configuration

This pipeline was generated to cut the adapters from the raw FASTQ files after sequencing.

### Install STAR Mapping Tools

```bash
cd /data/mckeeka/bulkRNA_RMS
conda create -n STARmap -c bioconda snakemake star SAMtools -y
conda activate STARmap
```

### Create Snakemake STAR Mapping Configuration File

```bash
nano STARmap_pipeline.smk

# Add the following code to the configuration file:

GENOME_FASTA = "reference/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
GTF = "reference/Homo_sapiens.GRCh38.115.gtf"
STAR_INDEX_DIR = "reference/STAR_index"

SAMPLES = glob_wildcards("bulkRNA/{sample}.fastq.R1.gz").sample

rule all:
  input:
    "reference/STAR_index",
    expand("run_bulkRNA/STAR_new/{sample}.Aligned.sortedByCoord.out.bam.bai", sample=SAMPLES)

rule star_index:
  input:
    fasta=GENOME_FASTA,
    gtf=GTF
  output:
    STAR_INDEX_DIR
  params:
    outdir=STAR_INDEX_DIR
  threads: 4
  shell:
    """
    mkdir -p {params.outdir}

    STAR \
        --runThreadN {threads} \
        --runMode genomeGenerate \
        --genomeDir {output} \
        --genomeFastaFiles {input.fasta} \
        --sjdbGTFfile {input.gtf} \
        --sjdbOverhang 99
    """

rule star_two_pass:
  input:
    r1 = "bulkRNA/{sample}.fastq.R1.gz",
    r2 = "bulkRNA/{sample}.fastq.R2.gz",
    index = STAR_INDEX_DIR
  output:
    bam = "run_bulkRNA/STAR_new/{sample}.Aligned.sortedByCoord.out.bam"
  log:
    "run_bulkRNA/logs/logs_STAR_new/{sample}.log"
  threads: 4
  shell:
      """
      STAR \
          --runThreadN {threads} \
          --genomeDir {input.index} \
          --readFilesIn {input.r1} {input.r2} \
          --readFilesCommand zcat \
          --twopassMode Basic \
          --chimSegmentMin 12 \
          --outFilterMultimapNmax 20 \
          --winAnchorMultimapNmax 50 \
          --outSAMtype BAM SortedByCoordinate \
          --outFileNamePrefix run_bulkRNA/STAR_new/{wildcards.sample}. \
          &> {log}
      """

rule samtools_index:
  input:
    bam = "run_bulkRNA/STAR_new/{sample}.Aligned.sortedByCoord.out.bam"
  output:
    bai = "run_bulkRNA/STAR_new/{sample}.Aligned.sortedByCoord.out.bam.bai"
  threads: 4
  shell:
    """
    samtools index {input.bam}
    """

```

### Run STAR Mapping Configuration File

The pipeline must be run using sbatch on the Biowulf cluster.

```bash
cd /data/mckeeka/bulkRNA_RMS
sbatch --cpus-per-task=4 --mem=64G --time=01-00:00:00 --wrap "snakemake -s STARmap_pipeline.smk --cores 4"
```

## Create STAR Mapping Pipeline Working Directory

The STAR Mapping pipeline requires a working directory where the FASTQ files can be accessed.

```bash
cd /data/mckeeka/bulkRNA_RMS/run_bulkRNA
mkdir STAR_filtered
cd /data/mckeeka/bulkRNA_RMS/run_bulkRNA/logs
mkdir logs_STAR_filtered
cd /data/mckeeka/bulkRNA_RMS
```

## Generate Filtered STAR Mapping Pipeline Configuration

This pipeline was generated to cut the adapters from the raw FASTQ files after sequencing.


### Create Filtered Snakemake STAR Mapping Configuration File

```bash
conda activate STARmap

nano STARmap_pipeline_filtered.smk

# Add the following code to the configuration file:

GENOME_FASTA = "reference/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
GTF = "reference/Homo_sapiens.GRCh38.115.gtf"
STAR_INDEX_DIR = "reference/STAR_index"

SAMPLES = glob_wildcards("run_bulkRNA/clean_FASTQ/{sample}.fastq.R1.clean.gz").sample

rule all:
  input:
    "reference/STAR_index",
    expand("run_bulkRNA/STAR_filtered/{sample}.Aligned.sortedByCoord.out.bam.bai", sample=SAMPLES)

rule star_index:
  input:
    fasta=GENOME_FASTA,
    gtf=GTF
  output:
    STAR_INDEX_DIR
  params:
    outdir=STAR_INDEX_DIR
  threads: 4
  shell:
    """
    mkdir -p {params.outdir}

    STAR \
        --runThreadN {threads} \
        --runMode genomeGenerate \
        --genomeDir {output} \
        --genomeFastaFiles {input.fasta} \
        --sjdbGTFfile {input.gtf} \
        --sjdbOverhang 99
    """

rule star_two_pass:
  input:
    r1 = "run_bulkRNA/clean_FASTQ/{sample}.fastq.R1.clean.gz",
    r2 = "run_bulkRNA/clean_FASTQ/{sample}.fastq.R2.clean.gz",
    index = STAR_INDEX_DIR
  output:
    bam = "run_bulkRNA/STAR_filtered/{sample}.Aligned.sortedByCoord.out.bam"
  log:
    "run_bulkRNA/logs/logs_STAR_filtered/{sample}.log"
  threads: 4
  shell:
      """
      STAR \
          --runThreadN {threads} \
          --genomeDir {input.index} \
          --readFilesIn {input.r1} {input.r2} \
          --readFilesCommand zcat \
          --twopassMode Basic \
          --chimSegmentMin 12 \
          --outFilterMultimapNmax 20 \
          --winAnchorMultimapNmax 50 \
          --outSAMtype BAM SortedByCoordinate \
          --outFileNamePrefix run_bulkRNA/STAR_filtered/{wildcards.sample}. \
          &> {log}
      """

rule samtools_index:
  input:
    bam = "run_bulkRNA/STAR_filtered/{sample}.Aligned.sortedByCoord.out.bam"
  output:
    bai = "run_bulkRNA/STAR_filtered/{sample}.Aligned.sortedByCoord.out.bam.bai"
  threads: 4
  shell:
    """
    samtools index {input.bam}
    """

```

### Run STAR Mapping Configuration File

The pipeline must be run using sbatch on the Biowulf cluster.

```bash
cd /data/mckeeka/bulkRNA_RMS
sbatch --cpus-per-task=4 --mem=64G --time=01-00:00:00 --wrap "snakemake -s STARmap_pipeline_filtered.smk --cores 4"
```

## Generate STAR QC Pipeline Configuration

This pipeline was generated to analyze the quality of STAR mapping on our unfiltered samples.

### Install STAR QC Tools

```bash
cd /data/mckeeka/bulkRNA_RMS
conda create -n STAR_QC -c conda-forge -c bioconda snakemake python=3.10 -y
conda activate STAR_QC
```

### Create Snakemake STAR QC Configuration File

```bash
nano STAR_QC.py


# Add the following code to the configuration file:

import csv
import glob
import os

STAR_LOG_DIR = "run_bulkRNA/STAR_new"
OUT_CSV = "run_bulkRNA/STAR_new/star_qc_summary.csv"

def parse_log_final_out(path):
    # Example filename: RD_1.Log.final.out -> RD_1
    sample = os.path.basename(path).replace(".Log.final.out", "")
    metrics = {"sample": sample}

    with open(path) as f:
        for line in f:
            if "|" not in line:
                continue
            key, val = line.split("|", 1)
            key = key.strip()
            val = val.strip()
            metrics[key] = val

    return metrics

logs = glob.glob(os.path.join(STAR_LOG_DIR, "*.Log.final.out"))
rows = [parse_log_final_out(p) for p in logs]

if not rows:
    raise SystemExit(f"No STAR log files found in {STAR_LOG_DIR}")

all_keys = sorted({k for row in rows for k in row.keys()})
os.makedirs(os.path.dirname(OUT_CSV), exist_ok=True)

with open(OUT_CSV, "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=all_keys)
    writer.writeheader()
    writer.writerows(rows)

print(f"Wrote {OUT_CSV}")

```

### Run STAR QC Configuration File

```bash
cd /data/mckeeka/bulkRNA_RMS
python STAR_QC.py
```

## Generate Filtered STAR QC Pipeline Configuration

This pipeline was generated to analyze the quality of STAR mapping on our filtered samples.

### Create Snakemake STAR QC Configuration File

```bash
conda activate STAR_QC

cd /data/mckeeka/bulkRNA_RMS
nano STAR_QC_filtered.py


# Add the following code to the configuration file:

import csv
import glob
import os

STAR_LOG_DIR = "run_bulkRNA/STAR_filtered"
OUT_CSV = "run_bulkRNA/STAR_filtered/star_qc_summary.csv"

def parse_log_final_out(path):
    # Example filename: RD_1.Log.final.out -> RD_1
    sample = os.path.basename(path).replace(".Log.final.out", "")
    metrics = {"sample": sample}

    with open(path) as f:
        for line in f:
            if "|" not in line:
                continue
            key, val = line.split("|", 1)
            key = key.strip()
            val = val.strip()
            metrics[key] = val

    return metrics

logs = glob.glob(os.path.join(STAR_LOG_DIR, "*.Log.final.out"))
rows = [parse_log_final_out(p) for p in logs]

if not rows:
    raise SystemExit(f"No STAR log files found in {STAR_LOG_DIR}")

all_keys = sorted({k for row in rows for k in row.keys()})
os.makedirs(os.path.dirname(OUT_CSV), exist_ok=True)

with open(OUT_CSV, "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=all_keys)
    writer.writeheader()
    writer.writerows(rows)

print(f"Wrote {OUT_CSV}")

```

### Run STAR QC Configuration File

```bash
cd /data/mckeeka/bulkRNA_RMS
python STAR_QC_filtered.py
```


## Create FeatureCounts Pipeline Working Directory

```bash
cd /data/mckeeka/bulkRNA_RMS/run_bulkRNA/
mkdir FeatureCounts
```

## Generate Feature Counts Pipeline Configuration

This pipeline was generated to count the features after STAR two-pass mapping on the unfiltered samples.

### Install Feature Counts Tools

```bash
cd /data/mckeeka/bulkRNA_RMS
conda create -n FeatureCounts -c conda-forge -c bioconda snakemake subread python=3.10 -y
conda activate FeatureCounts
```

### Create Snakemake Feature Counts Configuration File

```bash
nano FeatureCounts_pipeline.smk

# Add the following code to the configuration file:

from pathlib import Path
from glob import glob

GTF = "reference/Homo_sapiens.GRCh38.115.gtf"
BAM_DIR = "run_bulkRNA/STAR_new"
COUNT_DIR = "run_bulkRNA/FeatureCounts"

SAMPLES = sorted(
    Path(b).name.replace(".Aligned.sortedByCoord.out.bam", "")
    for b in glob(f"{BAM_DIR}/*.Aligned.sortedByCoord.out.bam")
)

rule all:
    input:
        f"{COUNT_DIR}/gene_counts.txt"
        
rule featurecounts:
    input:
        bams=lambda wildcards: expand(
            f"{BAM_DIR}" + "/{sample}.Aligned.sortedByCoord.out.bam",
            sample=SAMPLES
        ),
        gtf=GTF
    output:
        counts=f"{COUNT_DIR}/gene_counts.txt"
    threads: 6
    shell:
        """
        featureCounts \
            -T {threads} \
            -a {input.gtf} \
            -o {output.counts} \
            -p \
            -B \
            -C \
            {input.bams}
        """
```

### Run Feature Counts Configuration File

The pipeline must be run using sbatch on the Biowulf cluster.

```bash
cd /data/mckeeka/bulkRNA_RMS
sbatch --cpus-per-task=6 --time=01:00:00 --wrap "snakemake -s FeatureCounts_pipeline.smk --cores 6"
```

## Generate Feature Counts Filtered Pipeline Configuration

This pipeline was generated to count the features after STAR two-pass mapping on the filtered samples.


### Create Snakemake Feature Counts Filtered Configuration File

```bash
conda activate FeatureCounts
nano FeatureCounts_filtered_pipeline.smk

# Add the following code to the configuration file:

from pathlib import Path
from glob import glob

GTF = "reference/Homo_sapiens.GRCh38.115.gtf"
BAM_DIR = "run_bulkRNA/STAR_filtered"
COUNT_DIR = "run_bulkRNA/FeatureCounts"

SAMPLES = sorted(
    Path(b).name.replace(".Aligned.sortedByCoord.out.bam", "")
    for b in glob(f"{BAM_DIR}/*.Aligned.sortedByCoord.out.bam")
)

rule all:
    input:
        f"{COUNT_DIR}/gene_counts.txt"
        
rule featurecounts:
    input:
        bams=lambda wildcards: expand(
            f"{BAM_DIR}" + "/{sample}.Aligned.sortedByCoord.out.bam",
            sample=SAMPLES
        ),
        gtf=GTF
    output:
        counts=f"{COUNT_DIR}/gene_counts.txt"
    threads: 6
    shell:
        """
        featureCounts \
            -T {threads} \
            -a {input.gtf} \
            -o {output.counts} \
            -p \
            -B \
            -C \
            {input.bams}
        """
```

### Run Feature Counts Filtered Configuration File

The pipeline must be run using sbatch on the Biowulf cluster.

```bash
cd /data/mckeeka/bulkRNA_RMS
sbatch --cpus-per-task=6 --time=01:00:00 --wrap "snakemake -s FeatureCounts_filtered_pipeline.smk --cores 6"
```

### Install Feature Counts vs. Proteomics Tools

```bash
cd /data/mckeeka/bulkRNA_RMS
conda create -n r-prot -c conda-forge -c bioconda \
  r-base r-tidyverse r-readxl r-pheatmap \
  bioconductor-annotationdbi bioconductor-org.hs.eg.db bioconductor-biobase -y
conda activate r-prot
```

### Create Feature Counts vs. Proteomics Configuration File

```bash
nano proteomics.r

# Add the following code to the configuration file:

suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
  library(stringr)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
})

# -----------------------------
# File paths
# -----------------------------
rna_unfiltered_file <- "gene_counts_unfiltered.txt"
rna_filtered_file   <- "gene_counts_filtered.txt"

prot_files <- c(
  gel = "11122025_SmallProtein_gel_byband.xlsx",
  sec = "11132025_SmallProtein_SEC_byfraction.xlsx"
)

# -----------------------------
# 1) Read featureCounts files
# -----------------------------
rna_unfiltered_raw <- read.delim(rna_unfiltered_file,
                                 comment.char = "#",
                                 check.names = FALSE)

rna_filtered_raw <- read.delim(rna_filtered_file,
                               comment.char = "#",
                               check.names = FALSE)

rna_unfiltered <- rna_unfiltered_raw %>%
  dplyr::select(Geneid, starts_with("run_bulkRNA/")) %>%
  mutate(across(-Geneid, as.numeric)) %>%
  tibble::column_to_rownames("Geneid")

rna_filtered <- rna_filtered_raw %>%
  dplyr::select(Geneid, starts_with("run_bulkRNA/")) %>%
  mutate(across(-Geneid, as.numeric)) %>%
  tibble::column_to_rownames("Geneid")

map_ensembl_to_symbol <- function(df) {
  df <- as.data.frame(df)
  df$Geneid <- rownames(df)

  df$Symbol <- AnnotationDbi::mapIds(
    org.Hs.eg.db,
    keys = df$Geneid,
    column = "SYMBOL",
    keytype = "ENSEMBL",
    multiVals = "first"
  )

  df <- df %>%
    filter(!is.na(Symbol)) %>%
    dplyr::select(Symbol, where(is.numeric)) %>%
    group_by(Symbol) %>%
    summarise(across(where(is.numeric), sum), .groups = "drop")

  tibble::column_to_rownames(df, "Symbol")
}

rna_unfiltered_sym <- map_ensembl_to_symbol(rna_unfiltered)
rna_filtered_sym   <- map_ensembl_to_symbol(rna_filtered)

# -----------------------------
# 2) Read one proteomics file
# -----------------------------
read_proteomics_file <- function(path) {
  prot_raw <- read_xlsx(path)
  names(prot_raw) <- make.names(names(prot_raw))

  if (!"Description" %in% names(prot_raw)) {
    stop("No Description column found in: ", path)
  }

  prot_raw <- prot_raw %>%
    mutate(
      GeneSymbol = str_extract(Description, "GN=[^ ]+") %>%
        str_remove("^GN=")
    ) %>%
    mutate(
      GeneSymbol = ifelse(is.na(GeneSymbol), NA_character_, GeneSymbol),
      GeneSymbol = str_split_fixed(GeneSymbol, ",", 2)[, 1]
    )

  candidate_numeric <- names(prot_raw)[vapply(prot_raw, is.numeric, logical(1))]

  preferred_cols <- grep(
    "intensity|lfq|abundance|area|quant|signal|psm",
    candidate_numeric,
    ignore.case = TRUE,
    value = TRUE
  )

  if (length(preferred_cols) > 0) {
    prot_value_col <- preferred_cols[1]
  } else if ("X..PSMs" %in% names(prot_raw)) {
    prot_value_col <- "X..PSMs"
  } else if ("PSMs" %in% names(prot_raw)) {
    prot_value_col <- "PSMs"
  } else {
    stop("No usable quantitative proteomics column found in: ", path)
  }

  message("Using proteomics numeric column for ", basename(path), ": ", prot_value_col)

  prot <- prot_raw %>%
    dplyr::select(GeneSymbol, all_of(prot_value_col)) %>%
    filter(!is.na(GeneSymbol), GeneSymbol != "") %>%
    mutate(across(all_of(prot_value_col), as.numeric)) %>%
    group_by(GeneSymbol) %>%
    summarise(value = mean(.data[[prot_value_col]], na.rm = TRUE), .groups = "drop") %>%
    tibble::column_to_rownames("GeneSymbol")

  log2(as.matrix(prot) + 1)
}

# -----------------------------
# 3) Correlation function
# -----------------------------
calc_cor <- function(rna_data, prot_data, method_name) {
  common_genes <- intersect(rownames(rna_data), rownames(prot_data))

  if (length(common_genes) == 0) {
    message("No overlapping genes for: ", method_name)
    return(tibble(method = character(), sample = character(),
                  pearson = numeric(), spearman = numeric()))
  }

  rna <- rna_data[common_genes, , drop = FALSE]
  prot_vec <- prot_data[common_genes, 1, drop = TRUE]

  map_dfr(colnames(rna), function(s) {
    tibble(
      method = method_name,
      sample = s,
      pearson = cor(rna[, s], prot_vec, use = "complete.obs", method = "pearson"),
      spearman = cor(rna[, s], prot_vec, use = "complete.obs", method = "spearman")
    )
  })
}

# -----------------------------
# 4) Run for each proteomics file
# -----------------------------
all_results <- map_dfr(names(prot_files), function(prot_name) {
  prot_log <- read_proteomics_file(prot_files[[prot_name]])

  bind_rows(
    calc_cor(rna_unfiltered_sym, prot_log, "Unfiltered"),
    calc_cor(rna_filtered_sym,   prot_log, "Filtered")
  ) %>%
    mutate(proteomics = prot_name)
})

print(all_results)

write.csv(all_results, "rna_vs_proteomics_all_methods.csv", row.names = FALSE)

# -----------------------------
# 5) Clean sample labels
# -----------------------------
all_results <- all_results %>%
  mutate(
    sample_short = case_when(
      str_detect(sample, "RD_1") ~ "RD1",
      str_detect(sample, "RD_2") ~ "RD2",
      str_detect(sample, "RD_3") ~ "RD3",
      str_detect(sample, "RH30_1") ~ "RH30_1",
      str_detect(sample, "RH30_2") ~ "RH30_2",
      str_detect(sample, "RH30_3") ~ "RH30_3",
      TRUE ~ sample
    ),
    sample_short = factor(sample_short,
                          levels = c("RD1", "RD2", "RD3", "RH30_1", "RH30_2", "RH30_3"))
  )

# -----------------------------
# 6) Summary: which proteomics file correlates the most?
# -----------------------------
summary_by_file <- all_results %>%
  group_by(proteomics, method) %>%
  summarise(
    mean_pearson = mean(pearson, na.rm = TRUE),
    sd_pearson   = sd(pearson, na.rm = TRUE),
    .groups = "drop"
  )

print(summary_by_file)

best_file <- summary_by_file %>%
  group_by(proteomics) %>%
  summarise(overall_mean_pearson = mean(mean_pearson, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(overall_mean_pearson))

print(best_file)

# -----------------------------
# 7) Presentation-friendly plots
# -----------------------------
suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
  library(stringr)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
})

# -----------------------------
# File paths
# -----------------------------
rna_unfiltered_file <- "gene_counts_unfiltered.txt"
rna_filtered_file   <- "gene_counts_filtered.txt"

prot_files <- c(
  gel = "11122025_SmallProtein_gel_byband.xlsx",
  sec = "11132025_SmallProtein_SEC_byfraction.xlsx"
)

# -----------------------------
# 1) Read featureCounts files
# -----------------------------
rna_unfiltered_raw <- read.delim(rna_unfiltered_file,
                                 comment.char = "#",
                                 check.names = FALSE)

rna_filtered_raw <- read.delim(rna_filtered_file,
                               comment.char = "#",
                               check.names = FALSE)

rna_unfiltered <- rna_unfiltered_raw %>%
  dplyr::select(Geneid, starts_with("run_bulkRNA/")) %>%
  mutate(across(-Geneid, as.numeric)) %>%
  tibble::column_to_rownames("Geneid")

rna_filtered <- rna_filtered_raw %>%
  dplyr::select(Geneid, starts_with("run_bulkRNA/")) %>%
  mutate(across(-Geneid, as.numeric)) %>%
  tibble::column_to_rownames("Geneid")

map_ensembl_to_symbol <- function(df) {
  df <- as.data.frame(df)
  df$Geneid <- rownames(df)

  df$Symbol <- AnnotationDbi::mapIds(
    org.Hs.eg.db,
    keys = df$Geneid,
    column = "SYMBOL",
    keytype = "ENSEMBL",
    multiVals = "first"
  )

  df <- df %>%
    filter(!is.na(Symbol)) %>%
    dplyr::select(Symbol, where(is.numeric)) %>%
    group_by(Symbol) %>%
    summarise(across(where(is.numeric), sum), .groups = "drop")

  tibble::column_to_rownames(df, "Symbol")
}

rna_unfiltered_sym <- map_ensembl_to_symbol(rna_unfiltered)
rna_filtered_sym   <- map_ensembl_to_symbol(rna_filtered)

# -----------------------------
# 2) Read one proteomics file
# -----------------------------
read_proteomics_file <- function(path) {
  prot_raw <- read_xlsx(path)
  names(prot_raw) <- make.names(names(prot_raw))

  if (!"Description" %in% names(prot_raw)) {
    stop("No Description column found in: ", path)
  }

  prot_raw <- prot_raw %>%
    mutate(
      GeneSymbol = str_extract(Description, "GN=[^ ]+") %>%
        str_remove("^GN=")
    ) %>%
    mutate(
      GeneSymbol = ifelse(is.na(GeneSymbol), NA_character_, GeneSymbol),
      GeneSymbol = str_split_fixed(GeneSymbol, ",", 2)[, 1]
    )

  candidate_numeric <- names(prot_raw)[vapply(prot_raw, is.numeric, logical(1))]

  preferred_cols <- grep(
    "intensity|lfq|abundance|area|quant|signal|psm",
    candidate_numeric,
    ignore.case = TRUE,
    value = TRUE
  )

  if (length(preferred_cols) > 0) {
    prot_value_col <- preferred_cols[1]
  } else if ("X..PSMs" %in% names(prot_raw)) {
    prot_value_col <- "X..PSMs"
  } else if ("PSMs" %in% names(prot_raw)) {
    prot_value_col <- "PSMs"
  } else {
    stop("No usable quantitative proteomics column found in: ", path)
  }

  message("Using proteomics numeric column for ", basename(path), ": ", prot_value_col)

  prot <- prot_raw %>%
    dplyr::select(GeneSymbol, all_of(prot_value_col)) %>%
    filter(!is.na(GeneSymbol), GeneSymbol != "") %>%
    mutate(across(all_of(prot_value_col), as.numeric)) %>%
    group_by(GeneSymbol) %>%
    summarise(value = mean(.data[[prot_value_col]], na.rm = TRUE), .groups = "drop") %>%
    tibble::column_to_rownames("GeneSymbol")

  log2(as.matrix(prot) + 1)
}

# -----------------------------
# 3) Correlation function
# -----------------------------
calc_cor <- function(rna_data, prot_data, method_name) {
  common_genes <- intersect(rownames(rna_data), rownames(prot_data))

  if (length(common_genes) == 0) {
    message("No overlapping genes for: ", method_name)
    return(tibble(method = character(), sample = character(),
                  pearson = numeric(), spearman = numeric()))
  }

  rna <- rna_data[common_genes, , drop = FALSE]
  prot_vec <- prot_data[common_genes, 1, drop = TRUE]

  map_dfr(colnames(rna), function(s) {
    tibble(
      method = method_name,
      sample = s,
      pearson = cor(rna[, s], prot_vec, use = "complete.obs", method = "pearson"),
      spearman = cor(rna[, s], prot_vec, use = "complete.obs", method = "spearman")
    )
  })
}

# -----------------------------
# 4) Run for each proteomics file
# -----------------------------
all_results <- map_dfr(names(prot_files), function(prot_name) {
  prot_log <- read_proteomics_file(prot_files[[prot_name]])

  bind_rows(
    calc_cor(rna_unfiltered_sym, prot_log, "Unfiltered"),
    calc_cor(rna_filtered_sym,   prot_log, "Filtered")
  ) %>%
    mutate(proteomics = prot_name)
})

print(all_results)

write.csv(all_results, "rna_vs_proteomics_all_methods.csv", row.names = FALSE)

# -----------------------------
# 5) Clean sample labels
# -----------------------------
all_results <- all_results %>%
  mutate(
    sample_short = case_when(
      str_detect(sample, "RD_1") ~ "RD1",
      str_detect(sample, "RD_2") ~ "RD2",
      str_detect(sample, "RD_3") ~ "RD3",
      str_detect(sample, "RH30_1") ~ "RH30_1",
      str_detect(sample, "RH30_2") ~ "RH30_2",
      str_detect(sample, "RH30_3") ~ "RH30_3",
      TRUE ~ sample
    ),
    sample_short = factor(sample_short,
                          levels = c("RD1", "RD2", "RD3", "RH30_1", "RH30_2", "RH30_3"))
  )

# -----------------------------
# 6) Summary: which proteomics file correlates the most?
# -----------------------------
summary_by_file <- all_results %>%
  group_by(proteomics, method) %>%
  summarise(
    mean_pearson = mean(pearson, na.rm = TRUE),
    sd_pearson   = sd(pearson, na.rm = TRUE),
    .groups = "drop"
  )

print(summary_by_file)

best_file <- summary_by_file %>%
  group_by(proteomics) %>%
  summarise(overall_mean_pearson = mean(mean_pearson, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(overall_mean_pearson))

print(best_file)

# -----------------------------
# 7) Presentation-friendly plots
# -----------------------------
p1 <- ggplot(all_results, aes(x = sample_short, y = pearson, color = proteomics, group = proteomics)) +
  geom_point(size = 3, position = position_dodge(width = 0.25)) +
  geom_line(position = position_dodge(width = 0.25)) +
  facet_wrap(~method) +
  coord_cartesian(ylim = c(0, max(all_results$pearson, na.rm = TRUE) + 0.05)) +
  theme_minimal(base_size = 14) +
  labs(
    title = "RNA vs Proteomics Correlation",
    x = "Sample",
    y = "Pearson correlation",
    color = "Proteomics file"
  )

print(p1)

p2 <- ggplot(summary_by_file, aes(x = proteomics, y = mean_pearson, fill = method)) +
  geom_col(position = "dodge", width = 0.7) +
  geom_errorbar(
    aes(ymin = mean_pearson - sd_pearson, ymax = mean_pearson + sd_pearson),
    position = position_dodge(width = 0.7),
    width = 0.15
  ) +
  coord_cartesian(ylim = c(0, max(summary_by_file$mean_pearson + summary_by_file$sd_pearson, na.rm = TRUE) + 0.05)) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Which proteomics file matches RNA best?",
    x = "Proteomics file",
    y = "Mean Pearson correlation",
    fill = "RNA method"
  )

print(p2)

# -----------------------------
# 8) Save figures
# -----------------------------
ggsave("rna_vs_proteomics_by_file.png", p1, width = 11, height = 6, dpi = 300)
ggsave("proteomics_file_comparison.png", p2, width = 8, height = 5, dpi = 300)

p2 <- ggplot(summary_by_file, aes(x = proteomics, y = mean_pearson, fill = method)) +
  geom_col(position = "dodge", width = 0.7) +
  geom_errorbar(
    aes(ymin = mean_pearson - sd_pearson, ymax = mean_pearson + sd_pearson),
    position = position_dodge(width = 0.7),
    width = 0.15
  ) +
  coord_cartesian(ylim = c(0, max(summary_by_file$mean_pearson + summary_by_file$sd_pearson, na.rm = TRUE) + 0.05)) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Which proteomics file matches RNA best?",
    x = "Proteomics file",
    y = "Mean Pearson correlation",
    fill = "RNA method"
  )

print(p2)

# -----------------------------
# 8) Save figures
# -----------------------------
ggsave("rna_vs_proteomics_by_file.png", p1, width = 11, height = 6, dpi = 300)
ggsave("proteomics_file_comparison.png", p2, width = 8, height = 5, dpi = 300)
```

### Run Feature Counts vs. Proteomics Configuration File

The pipeline must be run using sbatch on the Biowulf cluster.

```bash
cd /data/mckeeka/bulkRNA_RMS
Rscript proteomics.r
```

### Install Microproteins Tools

```bash
cd /data/mckeeka/bulkRNA_RMS
conda create -n Microproteins -c conda-forge -c bioconda snakemake gffread transdecoder r-base python=3.10 -y
conda activate Microproteins
```

### Create Microproteins Configuration File

```bash
nano Microproteins.smk

MIN_MICROPROTEIN_AA = 100
GENOME_FASTA = "reference/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
GTF = "reference/Homo_sapiens.GRCh38.115.gtf"
COUNTS = "run_bulkRNA/FeatureCounts/gene_counts_filtered.txt"
OUT_DIR = "run_bulkRNA/Microproteins"

rule all:
    input:
        f"{OUT_DIR}/orf/transcripts.fa",
        f"{OUT_DIR}/orf/transcripts.fa.transdecoder.pep",
        f"{OUT_DIR}/orf/transcripts.fa.transdecoder.gff3",
        f"{OUT_DIR}/orf/transcripts.fa.transdecoder.cds",
        f"{OUT_DIR}/orf/orf_metadata.tsv",
        f"{OUT_DIR}/orf/microproteins.tsv",
        f"{OUT_DIR}/qc/library_sizes.pdf",
        f"{OUT_DIR}/qc/count_distribution.pdf"

rule extract_transcripts:
    input:
        genome=GENOME_FASTA,
        gtf=GTF
    output:
        fa=f"{OUT_DIR}/orf/transcripts.fa"
    shell:
        """
        mkdir -p {OUT_DIR}/orf
        gffread {input.gtf} -g {input.genome} -w {output.fa}
        """

rule transdecoder_longorfs:
    input:
        fa=f"{OUT_DIR}/orf/transcripts.fa"
    output:
        flag=f"{OUT_DIR}/orf/.longorfs.done"
    shell:
        """
        cd {OUT_DIR}/orf
        TransDecoder.LongOrfs -t transcripts.fa
        touch .longorfs.done
        """

rule transdecoder_predict:
    input:
        fa=f"{OUT_DIR}/orf/transcripts.fa",
        flag=f"{OUT_DIR}/orf/.longorfs.done"
    output:
        pep=f"{OUT_DIR}/orf/transcripts.fa.transdecoder.pep",
        gff3=f"{OUT_DIR}/orf/transcripts.fa.transdecoder.gff3",
        cds=f"{OUT_DIR}/orf/transcripts.fa.transdecoder.cds"
    shell:
        """
        cd {OUT_DIR}/orf
        TransDecoder.Predict -t transcripts.fa
        """

rule build_orf_metadata:
    input:
        pep=f"{OUT_DIR}/orf/transcripts.fa.transdecoder.pep"
    output:
        tsv=f"{OUT_DIR}/orf/orf_metadata.tsv"
    run:
        import re

        def parse_fasta(path):
            header = None
            seq_parts = []
            with open(path) as fh:
                for line in fh:
                    line = line.strip()
                    if not line:
                        continue
                    if line.startswith(">"):
                        if header is not None:
                            yield header, "".join(seq_parts)
                        header = line[1:]
                        seq_parts = []
                    else:
                        seq_parts.append(line)
                if header is not None:
                    yield header, "".join(seq_parts)

        with open(output.tsv, "w") as out:
            out.write("orf_id\ttranscript_id\tlength_aa\tsequence\tcategory\n")
            for header, seq in parse_fasta(input.pep):
                orf_id = header.split()[0]

                # Try to recover transcript ID from TransDecoder-style headers
                transcript_id = orf_id
                if ".p" in orf_id:
                    transcript_id = orf_id.rsplit(".p", 1)[0]

                m = re.search(r"len=(\d+)", header)
                length_aa = int(m.group(1)) if m else len(seq)
                category = "microprotein" if length_aa < MIN_MICROPROTEIN_AA else "protein"

                out.write(f"{orf_id}\t{transcript_id}\t{length_aa}\t{seq}\t{category}\n")

rule split_microproteins:
    input:
        tsv=f"{OUT_DIR}/orf/orf_metadata.tsv"
    output:
        tsv=f"{OUT_DIR}/orf/microproteins.tsv"
    run:
        import csv

        with open(input.tsv) as inf, open(output.tsv, "w", newline="") as outf:
            reader = csv.DictReader(inf, delimiter="\t")
            writer = csv.DictWriter(outf, fieldnames=reader.fieldnames, delimiter="\t")
            writer.writeheader()
            for row in reader:
                if row["category"] == "microprotein":
                    writer.writerow(row)

rule count_qc_plots:
    input:
        counts="run_bulkRNA/FeatureCounts/gene_counts_filtered.txt"
    output:
        lib="run_bulkRNA/Microproteins/qc/library_sizes.pdf",
        dist="run_bulkRNA/Microproteins/qc/count_distribution.pdf"
    shell:
        r"""
        mkdir -p run_bulkRNA/Microproteins/qc
        Rscript -e '
          counts <- read.delim(
            "{input.counts}",
            header = TRUE,
            sep = "\t",
            comment.char = "#",
            check.names = FALSE,
            fill = TRUE,
            quote = "",
            stringsAsFactors = FALSE
          )

          annotation_cols <- c("Geneid", "Chr", "Start", "End", "Strand", "Length")
          sample_counts <- counts[, !(names(counts) %in% annotation_cols), drop = FALSE]

          if ("Geneid" %in% names(counts)) {{
            rownames(sample_counts) <- counts$Geneid
          }}

          sample_counts <- as.data.frame(lapply(sample_counts, as.numeric))
          rownames(sample_counts) <- if ("Geneid" %in% names(counts)) counts$Geneid else rownames(counts)

          pdf("{output.lib}")
          barplot(colSums(sample_counts, na.rm = TRUE),
                  las = 2,
                  main = "Library sizes",
                  ylab = "Total counts")
          dev.off()

          pdf("{output.dist}")
          boxplot(log2(sample_counts + 1),
                  las = 2,
                  main = "Log2 count distribution",
                  ylab = "log2(counts + 1)")
          dev.off()
        '
        """
```

### Run Microproteins Configuration File

The pipeline must be run using sbatch on the Biowulf cluster.

```bash
sbatch --time=00-02:00:00  --cpus-per-task=8 --mem=32G --wrap="snakemake -s Microproteins.smk --cores 8"
```



### Run Novel Microproteins Configuration File

conda create -n smorf \
  -c conda-forge -c bioconda \
  r-base=4.4 \
  r-data.table \
  r-stringr \
  bioconductor-biostrings \
  bioconductor-rtracklayer \
  diamond \
  -y

conda activate smorf

cd /data/mckeeka/bulkRNA_RMS
nano Microproteins_Novel.r

suppressPackageStartupMessages({
  library(Biostrings)
  library(data.table)
  library(rtracklayer)
  library(stringr)
})

# =========================
# User settings
# =========================

pep_fasta <- "run_bulkRNA/Microproteins/orf/transcripts.fa.transdecoder_dir/longest_orfs.pep"
gtf_file  <- "reference/Homo_sapiens.GRCh38.115.gtf" 
diamond_db <- "reference/uniprot_full.dmnd"     
out_prefix <- "smorf_results"

max_aa_len <- 100
diamond_evalue <- "1e-5"
diamond_max_hits <- 5

# =========================
# Helper functions
# =========================

parse_transdecoder_header <- function(hdr) {
  # Example:
  # ENST00000832824.p1 type:complete gc:universal ENST00000832824:764-345(-)

  first_token <- sub(" .*", "", hdr)

  # transcript/ORF ID at the beginning
  tx_orf <- first_token
  tx_id <- sub("\\.p[0-9]+$", "", tx_orf)
  orf_id <- tx_orf

  # extract coordinate field at end of header
  coord_str <- NA_character_
  if (grepl(":[0-9]+-[0-9]+\\([+-]\\)$", hdr)) {
    coord_str <- sub(".* ([^ ]+:[0-9]+-[0-9]+\\([+-]\\))$", "\\1", hdr)
  }

  transcript_id_from_coord <- NA_character_
  tx_start <- NA_integer_
  tx_end <- NA_integer_
  strand <- NA_character_

  if (!is.na(coord_str)) {
    # ENST00000832824:764-345(-)
    m <- str_match(coord_str, "^(.+):([0-9]+)-([0-9]+)\\(([+-])\\)$")
    transcript_id_from_coord <- m[, 2]
    tx_start <- suppressWarnings(as.integer(m[, 3]))
    tx_end <- suppressWarnings(as.integer(m[, 4]))
    strand <- m[, 5]
  }

  data.table(
    header = hdr,
    orf_id = orf_id,
    tx_id = tx_id,
    tx_id_from_coord = transcript_id_from_coord,
    tx_start = tx_start,
    tx_end = tx_end,
    strand = strand
  )
}

clean_interval <- function(a, b) {
  c(min(a, b), max(a, b))
}

# Map transcript coordinates to genomic coordinates using exon structure.
# Exons must be a data.frame with columns: seqnames, start, end, strand, transcript_id
map_tx_interval_to_genome <- function(tx_start, tx_end, exons_df) {
  if (nrow(exons_df) == 0 || any(is.na(c(tx_start, tx_end)))) {
    return(list(
      chrom = NA_character_,
      genomic_start = NA_integer_,
      genomic_end = NA_integer_,
      strand = NA_character_
    ))
  }

  exons_df <- as.data.table(exons_df)

  strand <- as.character(unique(exons_df$strand))[1]
  chrom  <- as.character(unique(exons_df$seqnames))[1]

  if (strand == "+") {
    setorder(exons_df, start, end)
  } else {
    setorder(exons_df, -start, -end)
  }

  exons_df[, exon_len := end - start + 1L]
  exons_df[, tx_exon_start := cumsum(c(0L, head(exon_len, -1L))) + 1L]
  exons_df[, tx_exon_end := cumsum(exon_len)]

  map_one_pos <- function(tx_pos) {
    hit <- exons_df[tx_exon_start <= tx_pos & tx_exon_end >= tx_pos][1]
    if (nrow(hit) == 0) return(NA_integer_)

    offset <- tx_pos - hit$tx_exon_start
    if (strand == "+") {
      hit$start + offset
    } else {
      hit$end - offset
    }
  }

  g1 <- map_one_pos(tx_start)
  g2 <- map_one_pos(tx_end)

  if (is.na(g1) || is.na(g2)) {
    return(list(
      chrom = chrom,
      genomic_start = NA_integer_,
      genomic_end = NA_integer_,
      strand = as.character(strand)
    ))
  }

  list(
    chrom = chrom,
    genomic_start = min(g1, g2),
    genomic_end = max(g1, g2),
    strand = as.character(strand)
  )
}

pick_best_diamond_hit <- function(hits) {
  # Expected columns from outfmt 6:
  # qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen
  if (nrow(hits) == 0) return(NULL)

  hits <- copy(hits)  # important: do not modify .SD directly
  setDT(hits)

  hits[, qcov := as.numeric(length) / as.numeric(qlen)]
  hits[, scov := as.numeric(length) / as.numeric(slen)]

  # Best hit = smallest evalue, then highest bitscore
  setorder(hits, evalue, -bitscore, -pident, -qcov)
  hits[1]
}

classify_hit <- function(best_hit) {
  if (is.null(best_hit) || nrow(best_hit) == 0) return("no_hit")

  pident <- as.numeric(best_hit$pident)
  qcov <- as.numeric(best_hit$qcov)
  evalue <- as.numeric(best_hit$evalue)

  if (is.na(pident) || is.na(qcov) || is.na(evalue)) return("uncertain")

  if (evalue <= 1e-20 && pident >= 90 && qcov >= 0.80) {
    return("known")
  }

  if (evalue <= 1e-5 && pident >= 30 && qcov >= 0.50) {
    return("homologous")
  }

  "novel_candidate"
}

# =========================
# Read and filter peptides
# =========================

pep <- readAAStringSet(pep_fasta)
hdrs <- names(pep)
aa_len <- width(pep)

meta <- rbindlist(lapply(hdrs, parse_transdecoder_header), fill = TRUE)
meta[, aa_len := aa_len]

# Keep only ORFs under threshold
keep_idx <- aa_len <= max_aa_len
pep_smorf <- pep[keep_idx]
meta_smorf <- meta[keep_idx]

# Write smORF peptide FASTA
smorf_pep_fa <- paste0(out_prefix, ".smorfs_lt", max_aa_len, "aa.pep.fa")
names(pep_smorf) <- meta_smorf$orf_id
writeXStringSet(pep_smorf, smorf_pep_fa)

# =========================
# Optional DIAMOND search
# =========================

diamond_tsv <- paste0(out_prefix, ".diamond.tsv")
best_hits_dt <- NULL

if (!is.null(diamond_db) && file.exists(diamond_db)) {
  query_fa <- smorf_pep_fa

  diamond_cmd <- c(
    "blastp",
    "-q", query_fa,
    "-d", diamond_db,
    "-o", diamond_tsv,
    "-e", diamond_evalue,
    "-k", as.character(diamond_max_hits),
    "--quiet",
    "--outfmt", "6",
    "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
    "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen", "slen"
  )

  message("Running DIAMOND...")
  system2("diamond", diamond_cmd)

  if (file.exists(diamond_tsv) && file.info(diamond_tsv)$size > 0) {
    hits <- fread(diamond_tsv, header = FALSE)
    setnames(hits, c(
      "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
      "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen", "slen"
    ))

    best_hits_dt <- hits[, {
      best <- pick_best_diamond_hit(.SD)
      if (is.null(best)) {
        data.table(
          sseqid = NA_character_,
          pident = NA_real_,
          length = NA_integer_,
          evalue = NA_real_,
          bitscore = NA_real_,
          qlen = unique(qlen),
          slen = NA_real_,
          qcov = NA_real_,
          scov = NA_real_
        )
      } else {
        best[, .(
          sseqid, pident, length, evalue, bitscore, qlen, slen,
          qcov = length / qlen,
          scov = length / slen
        )]
      }
    }, by = qseqid]
  }
}

# =========================
# Optional GTF-based genomic mapping
# =========================

tx_anno <- NULL
tx_exons <- NULL

if (!is.null(gtf_file) && file.exists(gtf_file)) {
  message("Reading GTF START")
  flush.console()

  gtf <- import(gtf_file)

  message("Reading GTF DONE")
  flush.console()

  gtf_df <- as.data.frame(gtf)
  rm(gtf)
  gc()

  if (!("transcript_id" %in% names(gtf_df))) {
    stop("GTF does not contain a transcript_id column in metadata.")
  }

  anno_cols <- intersect(
    c("transcript_id", "gene_id", "gene_name", "gene_type", "transcript_type", "biotype"),
    names(gtf_df)
  )

  if (length(anno_cols) > 0) {
    tx_anno <- unique(gtf_df[, anno_cols, drop = FALSE])
  }

  message("Filtering exon records START")
  flush.console()

  exons_only <- gtf_df[gtf_df$type == "exon", , drop = FALSE]

  message("Filtering exon records DONE")
  flush.console()

  if (nrow(exons_only) > 0) {
    exons_only <- exons_only[, intersect(
      c("seqnames", "start", "end", "strand", "transcript_id"),
      names(exons_only)
    ), drop = FALSE]

    message("Building transcript exon index START")
    flush.console()

    tx_exons <- as.data.table(exons_only)
    setkey(tx_exons, transcript_id)

    message("Building transcript exon index DONE")
    flush.console()
  }
}
# =========================
# Build final smORF table
# =========================

smorf_dt <- copy(meta_smorf)

# Use transcript ID from coord if present, otherwise tx_id parsed from first token
smorf_dt[, transcript_id := fifelse(!is.na(tx_id_from_coord), tx_id_from_coord, tx_id)]

# Genomic mapping (if GTF available)
smorf_dt[, `:=`(genomic_chrom = NA_character_,
                genomic_start = NA_integer_,
                genomic_end = NA_integer_,
                genomic_strand = NA_character_)]

if (!is.null(tx_exons)) {
  map_res <- lapply(seq_len(nrow(smorf_dt)), function(i) {
    tx <- smorf_dt$transcript_id[i]
    if (is.na(tx) || !(tx %in% names(tx_exons))) {
      return(list(
        chrom = NA_character_,
        genomic_start = NA_integer_,
        genomic_end = NA_integer_,
        strand = NA_character_
      ))
    }

    tx_df <- tx_exons[J(tx), nomatch = 0L]
    tx_df <- as.data.frame(tx_df)
    # normalize to a minimal set of columns
    tx_df <- tx_df[, intersect(c("seqnames", "start", "end", "strand"), names(tx_df)), drop = FALSE]

    # ORF coordinates in the header may be reversed on minus strand
    interval <- clean_interval(smorf_dt$tx_start[i], smorf_dt$tx_end[i])
    map_tx_interval_to_genome(interval[1], interval[2], tx_df)
  })

  smorf_dt[, genomic_chrom := vapply(map_res, `[[`, character(1), "chrom")]
  smorf_dt[, genomic_start := vapply(map_res, `[[`, integer(1), "genomic_start")]
  smorf_dt[, genomic_end := vapply(map_res, `[[`, integer(1), "genomic_end")]
  smorf_dt[, genomic_strand := vapply(map_res, `[[`, character(1), "strand")]
}

# =========================
# Add similarity classification
# =========================

smorf_dt[, diamond_hit := NA_character_]
smorf_dt[, diamond_pident := NA_real_]
smorf_dt[, diamond_qcov := NA_real_]
smorf_dt[, diamond_evalue := NA_real_]
smorf_dt[, diamond_bitscore := NA_real_]
smorf_dt[, diamond_class := "no_hit"]

if (!is.null(best_hits_dt) && nrow(best_hits_dt) > 0) {
  setkey(best_hits_dt, qseqid)
  setkey(smorf_dt, orf_id)

  smorf_dt[best_hits_dt, on = c(orf_id = "qseqid"), `:=`(
    diamond_hit = i.sseqid,
    diamond_pident = i.pident,
    diamond_qcov = i.qcov,
    diamond_evalue = i.evalue,
    diamond_bitscore = i.bitscore
  )]

  
  smorf_dt[, diamond_class := mapply(function(h, pid, qc, ev) {
    if (is.na(h)) return("no_hit")
    if (!is.na(ev) && ev <= 1e-20 && !is.na(pid) && pid >= 90 && !is.na(qc) && qc >= 0.80) {
      return("known")
    }
    if (!is.na(ev) && ev <= 1e-5 && !is.na(pid) && pid >= 30 && !is.na(qc) && qc >= 0.50) {
      return("homologous")
    }
    "novel_candidate"
  }, diamond_hit, diamond_pident, diamond_qcov, diamond_evalue)]
}

# =========================
# Transcript annotation join if available
# =========================

if (!is.null(tx_anno)) {
  tx_anno_dt <- as.data.table(tx_anno)
  if ("transcript_id" %in% names(tx_anno_dt)) {
    setkey(tx_anno_dt, transcript_id)
    setkey(smorf_dt, transcript_id)
    smorf_dt <- tx_anno_dt[smorf_dt]
  }
}

# =========================
# Write outputs
# =========================

fwrite(smorf_dt, paste0(out_prefix, ".smorf_annotation.tsv"), sep = "\t", na = "NA")

# Also write a fasta of only the "novel_candidate" smORFs
novel_ids <- smorf_dt[diamond_class == "novel_candidate" | diamond_class == "no_hit", header]
novel_pep <- pep_smorf[names(pep_smorf) %in% novel_ids]
writeXStringSet(novel_pep, paste0(out_prefix, ".novel_candidates.pep.fa"))

message("Done.")
message("Wrote:")
message("  - ", smorf_pep_fa)
message("  - ", paste0(out_prefix, ".smorf_annotation.tsv"))
message("  - ", paste0(out_prefix, ".novel_candidates.pep.fa"))
if (!is.null(best_hits_dt) && nrow(best_hits_dt) > 0) {
  message("  - ", diamond_tsv)
}

sbatch --job-name=smorf --cpus-per-task=4 --mem=32G --time=12:00:00 --output=smorf_%j.out --error=smorf_%j.err --wrap="Rscript Microproteins_Novel.r"

sbatch --job-name=smorf --cpus-per-task=4 --mem=32G --time=12:00:00 --output=smorf_%j.out --error=smorf_%j.err --wrap="set -x; echo START $(date); source ~/.bashrc; conda activate smorf; echo ENV_OK; Rscript Microproteins_Novel.r; echo DONE $(date)"

###Second filtering

conda activate smorf

cd /data/mckeeka/bulkRNA_RMS
nano Microproteins_Novel_new.r

suppressPackageStartupMessages({
  library(Biostrings)
  library(data.table)
  library(stringr)
  library(rtracklayer)
})

# =========================
# User inputs
# =========================

pep_fasta    <- "run_bulkRNA/Microproteins/orf/transcripts.fa.transdecoder_dir/longest_orfs.pep"
gtf_file     <- "reference/Homo_sapiens.GRCh38.115.gtf"

out_prefix   <- "smorf_results"
max_aa_len   <- 100

diamond_evalue    <- "1e-5"
diamond_max_hits   <- 5

# DIAMOND databases
uniprot_full_db <- "reference/uniprot_full.dmnd"
sprot_db        <- "reference/uniprot_sprot.dmnd"
openprot_db     <- "reference/openprot_human.dmnd"
smprot_db       <- "reference/smprot2_human_ribo.dmnd"

# =========================
# Helper functions
# =========================

parse_transdecoder_header <- function(hdr) {
  # Example:
  # ENST00000832824.p1 type:complete gc:universal ENST00000832824:764-345(-)

  first_token <- sub(" .*", "", hdr)

  tx_orf <- first_token
  tx_id <- sub("\\.p[0-9]+$", "", tx_orf)
  orf_id <- tx_orf

  coord_str <- NA_character_
  if (grepl(":[0-9]+-[0-9]+\\([+-]\\)$", hdr)) {
    coord_str <- sub(".* ([^ ]+:[0-9]+-[0-9]+\\([+-]\\))$", "\\1", hdr)
  }

  transcript_id_from_coord <- NA_character_
  tx_start <- NA_integer_
  tx_end <- NA_integer_
  strand <- NA_character_

  if (!is.na(coord_str)) {
    m <- str_match(coord_str, "^(.+):([0-9]+)-([0-9]+)\\(([+-])\\)$")
    transcript_id_from_coord <- m[, 2]
    tx_start <- suppressWarnings(as.integer(m[, 3]))
    tx_end <- suppressWarnings(as.integer(m[, 4]))
    strand <- m[, 5]
  }

  data.table(
    header = hdr,
    orf_id = orf_id,
    tx_id = tx_id,
    tx_id_from_coord = transcript_id_from_coord,
    tx_start = tx_start,
    tx_end = tx_end,
    strand = strand
  )
}

clean_interval <- function(a, b) {
  c(min(a, b), max(a, b))
}

map_tx_interval_to_genome <- function(tx_start, tx_end, exons_df) {
  if (nrow(exons_df) == 0 || any(is.na(c(tx_start, tx_end)))) {
    return(list(
      chrom = NA_character_,
      genomic_start = NA_integer_,
      genomic_end = NA_integer_,
      strand = NA_character_
    ))
  }

  exons_df <- as.data.table(exons_df)

  strand <- unique(exons_df$strand)
  chrom <- unique(as.character(exons_df$seqnames))
  strand <- strand[1]
  chrom <- chrom[1]

  if (strand == "+") {
    setorder(exons_df, start, end)
  } else {
    setorder(exons_df, -start, -end)
  }

  exons_df[, exon_len := end - start + 1L]
  exons_df[, tx_exon_start := cumsum(c(0L, head(exon_len, -1L))) + 1L]
  exons_df[, tx_exon_end := cumsum(exon_len)]

  map_one_pos <- function(tx_pos) {
    hit <- exons_df[tx_exon_start <= tx_pos & tx_exon_end >= tx_pos][1]
    if (nrow(hit) == 0) return(NA_integer_)

    offset <- tx_pos - hit$tx_exon_start
    if (strand == "+") {
      hit$start + offset
    } else {
      hit$end - offset
    }
  }

  g1 <- map_one_pos(tx_start)
  g2 <- map_one_pos(tx_end)

  if (is.na(g1) || is.na(g2)) {
    return(list(
      chrom = chrom,
      genomic_start = NA_integer_,
      genomic_end = NA_integer_,
      strand = strand
    ))
  }

  list(
    chrom = chrom,
    genomic_start = min(g1, g2),
    genomic_end = max(g1, g2),
    strand = strand
  )
}

run_diamond <- function(query_fa, db, out_tsv, evalue = "1e-5", max_hits = 5) {
  diamond_cmd <- c(
    "blastp",
    "-q", query_fa,
    "-d", db,
    "-o", out_tsv,
    "-e", evalue,
    "-k", as.character(max_hits),
    "--quiet",
    "--outfmt", "6",
    "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
    "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen", "slen"
  )

  message("Running DIAMOND against: ", db)
  status <- system2("diamond", diamond_cmd)
  if (!identical(status, 0L)) {
    warning("DIAMOND returned non-zero exit status for: ", db)
  }
}

read_best_hits <- function(tsv_file, db_name) {
  if (!file.exists(tsv_file) || file.info(tsv_file)$size == 0) return(NULL)

  hits <- fread(tsv_file, header = FALSE)
  setnames(hits, c(
    "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
    "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen", "slen"
  ))

  hits[, `:=`(
    pident = as.numeric(pident),
    length = as.numeric(length),
    evalue = as.numeric(evalue),
    bitscore = as.numeric(bitscore),
    qlen = as.numeric(qlen),
    slen = as.numeric(slen)
  )]

  hits[, `:=`(
    qcov = length / qlen,
    scov = length / slen,
    db = db_name
  )]

  # Best hit = lowest evalue, then highest bitscore, then highest identity, then highest qcov
  setorder(hits, qseqid, evalue, -bitscore, -pident, -qcov)
  hits[, .SD[1], by = qseqid]
}

# =========================
# Read and filter peptides
# =========================

pep <- readAAStringSet(pep_fasta)
hdrs <- names(pep)
aa_len <- width(pep)

meta <- rbindlist(lapply(hdrs, parse_transdecoder_header), fill = TRUE)
meta[, aa_len := aa_len]

# Keep only ORFs below the threshold
keep_idx <- aa_len < max_aa_len
pep_smorf <- pep[keep_idx]
meta_smorf <- meta[keep_idx]

# Write filtered peptide FASTA
smorf_pep_fa <- paste0(out_prefix, ".smorfs_lt", max_aa_len, "aa.pep.fa")
names(pep_smorf) <- meta_smorf$orf_id
writeXStringSet(pep_smorf, smorf_pep_fa)

# =========================
# DIAMOND searches
# =========================

openprot_tsv <- paste0(out_prefix, ".openprot.tsv")
smprot_tsv   <- paste0(out_prefix, ".smprot.tsv")
sprot_tsv    <- paste0(out_prefix, ".uniprot_sprot.tsv")
uniprot_tsv  <- paste0(out_prefix, ".uniprot_full.tsv")

run_diamond(smorf_pep_fa, openprot_db, openprot_tsv, evalue = diamond_evalue, max_hits = diamond_max_hits)
run_diamond(smorf_pep_fa, smprot_db,   smprot_tsv,   evalue = diamond_evalue, max_hits = diamond_max_hits)
run_diamond(smorf_pep_fa, sprot_db,    sprot_tsv,    evalue = diamond_evalue, max_hits = diamond_max_hits)
run_diamond(smorf_pep_fa, uniprot_full_db, uniprot_tsv, evalue = diamond_evalue, max_hits = diamond_max_hits)

openprot_best <- read_best_hits(openprot_tsv, "openprot")
smprot_best   <- read_best_hits(smprot_tsv,   "smprot")
sprot_best    <- read_best_hits(sprot_tsv,    "sprot")
uniprot_best  <- read_best_hits(uniprot_tsv,  "uniprot_full")

# =========================
# GTF processing
# =========================

tx_anno <- NULL
tx_exons <- NULL

if (!is.null(gtf_file) && file.exists(gtf_file)) {
  message("Reading GTF...")
  gtf <- import(gtf_file)
  gtf_df <- as.data.frame(gtf)
  rm(gtf)
  gc()

  if (!("transcript_id" %in% names(gtf_df))) {
    stop("GTF does not contain a transcript_id column in metadata.")
  }

  anno_cols <- intersect(
    c("transcript_id", "gene_id", "gene_name", "gene_type", "transcript_type", "biotype"),
    names(gtf_df)
  )

  if (length(anno_cols) > 0) {
    tx_anno <- unique(gtf_df[, anno_cols, drop = FALSE])
  }

  exons_only <- gtf_df[gtf_df$type == "exon", , drop = FALSE]

  if (nrow(exons_only) > 0) {
    exons_only <- exons_only[, intersect(
      c("seqnames", "start", "end", "strand", "transcript_id"),
      names(exons_only)
    ), drop = FALSE]

    tx_exons <- as.data.table(exons_only)
    setkey(tx_exons, transcript_id)
  }
}

# =========================
# Build annotation table
# =========================

smorf_dt <- copy(meta_smorf)

# Use transcript ID from coordinate field if available
smorf_dt[, transcript_id := fifelse(!is.na(tx_id_from_coord), tx_id_from_coord, tx_id)]

# Initialize columns
smorf_dt[, `:=`(
  openprot_hit = NA_character_,
  openprot_pident = NA_real_,
  openprot_qcov = NA_real_,
  openprot_evalue = NA_real_,
  openprot_bitscore = NA_real_,

  smprot_hit = NA_character_,
  smprot_pident = NA_real_,
  smprot_qcov = NA_real_,
  smprot_evalue = NA_real_,
  smprot_bitscore = NA_real_,

  sprot_hit = NA_character_,
  sprot_pident = NA_real_,
  sprot_qcov = NA_real_,
  sprot_evalue = NA_real_,
  sprot_bitscore = NA_real_,

  uniprot_hit = NA_character_,
  uniprot_pident = NA_real_,
  uniprot_qcov = NA_real_,
  uniprot_evalue = NA_real_,
  uniprot_bitscore = NA_real_,

  genomic_chrom = NA_character_,
  genomic_start = NA_integer_,
  genomic_end = NA_integer_,
  genomic_strand = NA_character_
)]

# Attach DIAMOND hits by qseqid -> orf_id
if (!is.null(openprot_best) && nrow(openprot_best) > 0) {
  smorf_dt[openprot_best, on = c(orf_id = "qseqid"), `:=`(
    openprot_hit = i.sseqid,
    openprot_pident = i.pident,
    openprot_qcov = i.qcov,
    openprot_evalue = i.evalue,
    openprot_bitscore = i.bitscore
  )]
}

if (!is.null(smprot_best) && nrow(smprot_best) > 0) {
  smorf_dt[smprot_best, on = c(orf_id = "qseqid"), `:=`(
    smprot_hit = i.sseqid,
    smprot_pident = i.pident,
    smprot_qcov = i.qcov,
    smprot_evalue = i.evalue,
    smprot_bitscore = i.bitscore
  )]
}

if (!is.null(sprot_best) && nrow(sprot_best) > 0) {
  smorf_dt[sprot_best, on = c(orf_id = "qseqid"), `:=`(
    sprot_hit = i.sseqid,
    sprot_pident = i.pident,
    sprot_qcov = i.qcov,
    sprot_evalue = i.evalue,
    sprot_bitscore = i.bitscore
  )]
}

if (!is.null(uniprot_best) && nrow(uniprot_best) > 0) {
  smorf_dt[uniprot_best, on = c(orf_id = "qseqid"), `:=`(
    uniprot_hit = i.sseqid,
    uniprot_pident = i.pident,
    uniprot_qcov = i.qcov,
    uniprot_evalue = i.evalue,
    uniprot_bitscore = i.bitscore
  )]
}

# =========================
# Optional genomic mapping
# =========================

if (!is.null(tx_exons)) {
  map_res <- lapply(seq_len(nrow(smorf_dt)), function(i) {
    tx <- smorf_dt$transcript_id[i]
    if (is.na(tx) || !(tx %in% tx_exons$transcript_id)) {
      return(list(
        chrom = NA_character_,
        genomic_start = NA_integer_,
        genomic_end = NA_integer_,
        strand = NA_character_
      ))
    }

    tx_df <- tx_exons[J(tx), nomatch = 0L]
    tx_df <- as.data.frame(tx_df)
    tx_df <- tx_df[, intersect(c("seqnames", "start", "end", "strand"), names(tx_df)), drop = FALSE]

    interval <- clean_interval(smorf_dt$tx_start[i], smorf_dt$tx_end[i])
    map_tx_interval_to_genome(interval[1], interval[2], tx_df)
  })

smorf_dt[, genomic_chrom := vapply(map_res, function(x) as.character(x$chrom), character(1))]
smorf_dt[, genomic_start := vapply(map_res, function(x) as.integer(x$genomic_start), integer(1))]
smorf_dt[, genomic_end := vapply(map_res, function(x) as.integer(x$genomic_end), integer(1))]
smorf_dt[, genomic_strand := vapply(map_res, function(x) as.character(x$strand), character(1))]
}

# =========================
# Classification logic
# =========================

smorf_dt[, `:=`(
  hit_in_openprot = !is.na(openprot_hit),
  hit_in_smprot = !is.na(smprot_hit),

  strong_hit_in_sprot = !is.na(sprot_hit) &
    !is.na(sprot_evalue) & sprot_evalue <= 1e-20 &
    !is.na(sprot_pident) & sprot_pident >= 90 &
    !is.na(sprot_qcov) & sprot_qcov >= 0.80,

  weak_hit_in_uniprot = !is.na(uniprot_hit) &
    !is.na(uniprot_evalue) & uniprot_evalue <= 1e-5 &
    !is.na(uniprot_pident) & uniprot_pident >= 30 &
    !is.na(uniprot_qcov) & uniprot_qcov >= 0.50
)]

smorf_dt[, class := fifelse(
  hit_in_openprot | hit_in_smprot,
  "known_microprotein",
  fifelse(
    strong_hit_in_sprot,
    "known_canonical",
    fifelse(
      weak_hit_in_uniprot,
      "homologous",
      "novel_candidate"
    )
  )
)]

# =========================
# Final cleanup / output
# =========================

# Optional: order columns a little more nicely
setcolorder(smorf_dt, c(
  "header", "orf_id", "tx_id", "tx_id_from_coord", "transcript_id",
  "tx_start", "tx_end", "strand", "aa_len",
  "genomic_chrom", "genomic_start", "genomic_end", "genomic_strand",
  "openprot_hit", "openprot_pident", "openprot_qcov", "openprot_evalue", "openprot_bitscore",
  "smprot_hit", "smprot_pident", "smprot_qcov", "smprot_evalue", "smprot_bitscore",
  "sprot_hit", "sprot_pident", "sprot_qcov", "sprot_evalue", "sprot_bitscore",
  "uniprot_hit", "uniprot_pident", "uniprot_qcov", "uniprot_evalue", "uniprot_bitscore",
  "hit_in_openprot", "hit_in_smprot", "strong_hit_in_sprot", "weak_hit_in_uniprot",
  "class"
))

fwrite(
  smorf_dt,
  paste0(out_prefix, ".smorf_classification.tsv"),
  sep = "\t",
  na = "NA"
)

message("Done. Wrote: ", paste0(out_prefix, ".smorf_classification.tsv"))
