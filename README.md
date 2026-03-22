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
for f in *_RD_*_S*_R*_001.fastq.gz; do
  RD=$(echo "$f" | grep -o 'RD_[123]')
  Read=$(echo "$f" | grep -o 'R[12]')
  mv "$f" "${RD}.fastq.${Read}.gz"
done

for f in *_RH30_*_S*_R*_001.fastq.gz; do
  RH30=$(echo "$f" | grep -o 'RH30_[123]')
  Read=$(echo "$f" | grep -o 'R[12]')
  mv "$f" "${RH30}.fastq.${Read}.gz"
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
        cat $tmp_h1 $tmp_u1 > gzip -c {output.r1}
        cat $tmp_h2 $tmp_u2 > gzip -c {output.r2}

        #Remove temporary files
        rm -f $tmp_h1 $tmp_h2 $tmp_u1 $tmp_u2
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
    temp_bam = temp("run_bulkRNA/clean_FASTQ/bowtie2_output/{sample}_unmapped.bam"),
    r1 = "run_bulkRNA/clean_FASTQ/{sample}.fastq.R1.clean.gz",
    r2 = "run_bulkRNA/clean_FASTQ/{sample}.fastq.R2.clean.gz",
  log:
    "run_bulkRNA/logs/logs_cleanFASTQ/logs_bowtie2/{sample}.bowtie2_filter.log"
  shell:
    """
    {{
        samtools view -b -f 12 -F 256 {input.bam} > {output.temp_bam}

        tmp_r1=$(mktemp)
        tmp_r2=$(mktemp)

        bedtools bamtofastq \
            -i {output.temp_bam} \
            -fq tmp_r1 \
            -fq2 tmp_r2

        gzip -c $tmp_r1 > {output.r1}
        gzip -c $tmp_r2 > {output.r2}

        rm -f $tmp_r1 $tmp_r2

    }} &> {log}
    """

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
conda create -n cleanQC -c bioconda snakemake fastqc multiqc -y
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
sbatch --time=02:00:00 --wrap "snakemake -s cleanQC_pipeline.smk"
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

### Run CutAdapt Configuration File

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

## Create Indexing Pipeline Working Directory

The pipeline requires a working directory where the FASTQ files and reference transcriptome can be accessed. 

```bash
mkdir run_bulkRNA
cd /data/mckeeka/bulkRNA_sarcoma/run_bulkRNA

ln -s /data/mckeeka/bulkRNA_sarcoma/MCI_fastq_117_STS_FASTQ/
ln -s /data/mckeeka/bulkRNA_sarcoma/reference/
```

## Generate Pipeline Configuration

The pipeline requires a working directory where the FASTQ files and reference transcriptome can be accessed. You can symlink these files instead 
of copying them into the pipeline directory to prevent the duplication of large data files in your directory.

```bash
mkdir run_bulkRNA
cd /data/mckeeka/bulkRNA_sarcoma/run_bulkRNA

ln -s /data/mckeeka/bulkRNA_sarcoma/MCI_fastq_117_STS_FASTQ/
ln -s /data/mckeeka/bulkRNA_sarcoma/reference/
```





