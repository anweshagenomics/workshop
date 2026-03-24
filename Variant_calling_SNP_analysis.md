# 🧬 SNP & Variant Calling Pipeline

A complete bioinformatics pipeline for detecting **Single Nucleotide Polymorphisms (SNPs)** from raw Illumina paired-end reads using BWA-MEM, SAMtools, and BCFtools.

> Demonstrated on SARS-CoV-2 (GCF_009858895.2) with SRA accession `SRR37561623`

---

## 📋 Table of Contents

- [Overview](#overview)
- [Installation](#installation)
- [Input Data](#input-data)
- [Pipeline Steps](#pipeline-steps)
- [Output Files](#output-files)
- [Quick Reference](#quick-reference)
- [What's Next](#whats-next)

---

## Overview

This pipeline takes raw sequencing reads (FASTQ files) and a reference genome, aligns the reads, cleans up the alignment, and then identifies positions in the genome where your sample differs from the reference — those differences are the SNPs.

```
FASTQ reads + Reference FASTA
        │
        ▼
  [1] Index reference         bwa index, samtools faidx
        │
        ▼
  [2] Align reads             bwa mem
        │
        ▼
  [3] SAM → BAM + query sort  samtools view, samtools sort -n
        │
        ▼
  [4] Fix mate information    samtools fixmate -m
        │
        ▼
  [5] Sort by position        samtools sort
        │
        ▼
  [6] Mark duplicates         samtools markdup
        │
        ▼
  [7] Variant calling         bcftools mpileup | bcftools call
        │
        ▼
  [8] Extract SNPs            bcftools view -v snps
        │
        ▼
  [9] Filter SNPs             bcftools filter -e 'QUAL<30 || DP<10'
        │
        ▼
  covid.filtered.vcf.gz  ✅
```

---

## Installation

### Option 1 — Install Miniconda (if you don't have conda yet)

Miniconda is a lightweight installer for conda — a package and environment manager. If you already have Anaconda or Miniconda installed, skip to Option 2.

```bash
# Download the Miniconda installer (Linux x86_64)
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh

# Run the installer
bash Miniconda3-latest-Linux-x86_64.sh

# Follow the on-screen prompts, then reload your shell
source ~/.bashrc

# Verify conda is installed
conda --version
```

> For **macOS**, replace the URL with `https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh`

---

### Option 2 — Create the Environment and Install All Tools

All tools are available on the `bioconda` channel. The command below creates a fresh isolated environment called `bioset1` and installs everything in one go.

```bash
# Add the required channels (only needed once — sets the search order for packages)
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

# Create a new environment and install all tools at once
conda create -n bioset1 -c bioconda -c conda-forge \
    sra-tools \
    bwa \
    samtools \
    bcftools \
    tabix

# Activate the environment before running any pipeline commands
conda activate bioset1
```

> 💡 **What is a conda environment?** Think of it as a clean, isolated workspace where all your tools live together without interfering with anything else on your system. You activate it before working and deactivate it when done.

---

### Option 3 — Install Tools One by One (if the joint install fails)

If the combined install runs into conflicts, install each tool separately:

```bash
conda activate bioset1

conda install -c bioconda sra-tools
conda install -c bioconda bwa
conda install -c bioconda samtools
conda install -c bioconda bcftools
conda install -c bioconda tabix
```

---

### Verify Your Installation

After installing, confirm every tool is available and working:

```bash
conda activate bioset1

bwa 2>&1 | head -3
samtools --version | head -1
bcftools --version | head -1
prefetch --version
tabix --version | head -1
```

You should see version numbers printed for each tool. If any command fails with "command not found", re-run the install for that specific tool.

---

### Tool Summary

| Tool | Minimum Version | Purpose |
|------|----------------|---------|
| `sra-tools` | ≥3.0 | Download reads from NCBI SRA |
| `bwa` | ≥0.7.17 | Align reads to the reference genome |
| `samtools` | ≥1.17 | Work with alignment files (BAM/SAM) |
| `bcftools` | ≥1.17 | Call and filter variants |
| `tabix` | ≥1.17 | Index compressed VCF files |

---

## Input Data

### Reference Genome

The reference genome is the "gold standard" sequence that your reads will be compared against. Here we use the official SARS-CoV-2 reference from NCBI.

```bash
# Download the reference genome (compressed)
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/858/895/GCF_009858895.2_ASM985889v3/GCF_009858895.2_ASM985889v3_genomic.fna.gz

# Decompress it
gunzip GCF_009858895.2_ASM985889v3/GCF_009858895.2_ASM985889v3_genomic.fna.gz

# Rename it to something simpler to type
mv GCF_009858895.2_ASM985889v3/GCF_009858895.2_ASM985889v3_genomic.fna ref.fa
```

### Paired-End Reads

Paired-end sequencing reads each DNA fragment from both ends, producing two files per sample (R1 and R2). We download them directly from NCBI's Sequence Read Archive (SRA).

```bash
# Download the SRA archive
prefetch SRR37561623

# Split into two paired FASTQ files
fasterq-dump --split-files SRR37561623/SRR37561623.sra

# This produces:
#   SRR37561623_1.fastq  →  forward reads (R1)
#   SRR37561623_2.fastq  →  reverse reads (R2)
```

---

## Pipeline Steps

### Step 1 — Index the Reference Genome

**In simple terms:** Before we can search through the reference genome quickly, we need to build an "index" — similar to how a book's index lets you find a word without reading every page. We build two indexes here: one for BWA (speeds up read alignment) and one for SAMtools (enables jumping to any genome position instantly).

```bash
bwa index ref.fa
samtools faidx ref.fa
```

**Files created:**

| File | Created by | Used by |
|------|-----------|---------|
| `ref.fa.amb`, `.ann`, `.bwt`, `.pac`, `.sa` | `bwa index` | `bwa mem` in Step 2 |
| `ref.fa.fai` | `samtools faidx` | `bcftools mpileup` in Step 7 |

> ⏱️ For SARS-CoV-2 (~30 kb genome) this finishes in seconds. For large genomes like human (~3 Gb), `bwa index` can take 1–2 hours and is only done once per reference.

---

### Step 2 — Align Reads with BWA-MEM

**In simple terms:** This is the core mapping step. BWA-MEM takes each short sequencing read and finds the best matching location for it in the reference genome — like placing thousands of puzzle pieces onto a picture. We also attach a metadata "label" (called a read group tag) to every read, which carries the sample name and platform info into the final VCF output.

```bash
bwa mem \
  -t 8 \
  -R '@RG\tID:covid\tSM:covid\tPL:illumina' \
  ref.fa \
  SRR37561623_1.fastq \
  SRR37561623_2.fastq \
  > covid.sam
```

| Flag | Meaning |
|------|---------|
| `-t 8` | Use 8 CPU threads to speed up alignment |
| `-R '@RG\t...'` | Attach a read group tag to every aligned read |
| `> covid.sam` | Save the output as a SAM file (a human-readable text alignment format) |

> **Read group fields explained:** `ID` = a unique run identifier, `SM` = sample name (this appears as the column header in your VCF), `PL` = sequencing platform (Illumina).

---

### Step 3 — Convert SAM → BAM and Sort by Read Name

**In simple terms:** The SAM file is plain text and very large — it's inefficient for both storage and computation. We compress it into BAM format (about 4× smaller), then sort the reads so that each read and its paired partner appear next to each other in the file. This side-by-side pairing is needed for the next step to work correctly.

```bash
# Convert SAM (plain text) to BAM (compressed binary)
samtools view -@ 8 -bS covid.sam > covid.bam

# Sort reads by name so each read pair is grouped together
samtools sort -@ 8 -n -o covid.querysort.bam covid.bam
```

| Flag | Meaning |
|------|---------|
| `-@ 8` | Use 8 CPU threads |
| `-bS` | Read SAM format as input (`S`), write BAM format as output (`b`) |
| `-n` | Sort by read name (not genomic position) |

> ⚠️ **Why sort by name here?** The next step (`fixmate`) needs to see both reads of a pair simultaneously. Sorting by name ensures Read 1 and Read 2 of each pair are adjacent in the file.

---

### Step 4 — Fix Mate Information

**In simple terms:** In paired-end sequencing, each read should "know" the position and orientation of its partner read (its mate). After alignment, this information is sometimes missing or inconsistent. This step fills in all the mate details correctly, which is essential for the duplicate detection step that follows.

```bash
samtools fixmate -@ 8 -m covid.querysort.bam covid.fixmate.bam
```

| Flag | Meaning |
|------|---------|
| `-@ 8` | Use 8 CPU threads |
| `-m` | Add a mate score tag — this is required by `samtools markdup` in Step 6 |

> Think of it like making sure each pair of reads has each other's correct contact information before the next step checks them.

---

### Step 5 — Re-sort by Genomic Position

**In simple terms:** Now that the mate information is corrected, we re-sort the reads by their actual physical position in the genome (chromosome name and coordinate). All downstream tools — variant callers, genome browsers like IGV, GATK — expect reads to be in this order, like chapters in a book ordered by page number.

```bash
samtools sort -@ 8 -o covid.positionsort.bam covid.fixmate.bam
```

> This is the final sort. All files from this point forward stay in coordinate order.

---

### Step 6 — Mark PCR Duplicates

**In simple terms:** During library preparation in the lab, DNA is amplified by PCR, which can create multiple identical copies of the same original DNA fragment. If we don't account for this, the variant caller will count those copies as independent evidence and over-inflate confidence in a variant. This step identifies those duplicate reads and flags them so the variant caller knows to ignore them — without actually deleting the data.

```bash
# Identify and flag duplicate reads
samtools markdup -@ 8 covid.positionsort.bam covid.markdup.bam

# Create an index for fast random access (required by bcftools)
samtools index covid.markdup.bam
```

> 🗑️ Duplicates are **flagged, not deleted**. The raw data is preserved — downstream tools are simply instructed which reads to skip when making variant calls.

> ✅ `covid.markdup.bam` is your **final cleaned alignment file**. You can load it into [IGV](https://igv.org/) to visually inspect how reads pile up across the genome.

---

### Step 7 — Variant Calling

**In simple terms:** This is where SNPs are actually discovered. BCFtools scans every position across the genome, stacks up all the reads that map there, and counts how many show the reference base versus a different (alternate) base. Using a statistical model, it calculates the probability that a difference is a real biological variant — not just a sequencing error. All positions with strong evidence of a variant are written to a VCF file.

```bash
bcftools mpileup \
  --threads 8 \
  -Ou \
  -f ref.fa \
  covid.markdup.bam \
| bcftools call \
  --threads 8 \
  -mv \
  -Oz \
  -o covid.vcf.gz

bcftools index covid.vcf.gz
```

**Two commands joined by a pipe (`|`) — data flows directly from one to the other without writing to disk:**

| Command | Role |
|---------|------|
| `bcftools mpileup` | Stacks all reads at each genome position and reports the bases seen |
| `bcftools call` | Applies a Bayesian model to decide which positions are genuine variants |

| Flag | Tool | Meaning |
|------|------|---------|
| `-Ou` | mpileup | Stream output directly to `bcftools call` — no intermediate file |
| `-f ref.fa` | mpileup | Compare reads against this reference genome |
| `-m` | call | Use the multiallelic caller (handles sites with more than two alleles) |
| `-v` | call | Only output positions where a variant was detected |
| `-Oz` | call | Write a bgzip-compressed `.vcf.gz` output file |

---

### Step 8 — Extract SNPs Only

**In simple terms:** The VCF file from Step 7 contains all types of variants — SNPs (single letter changes like `A→T`), indels (insertions or deletions of bases), and others. For this pipeline we focus on SNPs only, so we filter out everything else and keep just the single nucleotide changes.

```bash
# Keep only SNP-type variants; discard indels, MNPs, and structural variants
bcftools view -v snps covid.vcf.gz -Oz -o covid.snps.vcf.gz

# Index the SNP file
bcftools index covid.snps.vcf.gz
```

> **SNP vs indel in plain language:**
> - SNP: one letter changes to another — e.g., position 23403 in the genome reads `A` in the reference but `T` in your sample
> - Indel: one or more letters are inserted or deleted — e.g., `AGTC` becomes `AGC` (one base removed)

---

### Step 9 — Filter High-Quality SNPs

**In simple terms:** Not every SNP in the previous file is reliable. Some are based on very few reads (low depth), making the call uncertain. Others have a low quality score, meaning the statistical model isn't confident the variant is real. We apply two straightforward cutoffs to keep only the SNPs we can trust, and discard the rest.

```bash
# Exclude SNPs that fail either quality check
bcftools filter \
  -e 'QUAL<30 || DP<10' \
  covid.snps.vcf.gz \
  -Oz -o covid.filtered.vcf.gz

# Index the final file
bcftools index covid.filtered.vcf.gz
```

| Filter | Plain-language meaning |
|--------|----------------------|
| `QUAL<30` | Less than 99.9% confident this is a real variant — too risky, discard |
| `DP<10` | Fewer than 10 reads supporting this position — not enough evidence, discard |

> The `-e` flag means **exclude** variants where the condition is true (it's the opposite of a "keep" filter).

> ✅ `covid.filtered.vcf.gz` is your **final result** — a high-confidence set of SNPs ready for annotation, comparison, or lineage typing.

---

## Output Files

| File | Description |
|------|-------------|
| `covid.markdup.bam` | Final processed alignment — load into IGV to visualize read pileups |
| `covid.vcf.gz` | All raw called variants — SNPs, indels, and other types |
| `covid.snps.vcf.gz` | SNP-only variants, before quality filtering |
| `covid.filtered.vcf.gz` | ✅ Final high-confidence SNPs (QUAL≥30 and DP≥10) |

---

## Quick Reference

```bash
# View your final SNPs (skip comment/header lines starting with #)
bcftools view covid.filtered.vcf.gz | grep -v "^#"

# Count how many SNPs passed the filter
bcftools view -H covid.filtered.vcf.gz | wc -l

# Get a summary report (number of SNPs, ts/tv ratio, etc.)
bcftools stats covid.filtered.vcf.gz | grep "^SN"

# Check alignment stats (% mapped, duplicates, etc.)
samtools flagstat covid.markdup.bam

# Deactivate the conda environment when you're done
conda deactivate
```

---

## What's Next?

- **Annotation:** Use [snpEff](https://pcingola.github.io/SnpEff/) to label each SNP with its gene name and predicted effect (synonymous, missense, stop-gain, etc.)
- **Multi-sample calling:** Pass multiple BAM files to `bcftools mpileup` at once for joint genotyping across samples — improves accuracy
- **Lineage assignment (SARS-CoV-2):** Submit your filtered VCF to [Pangolin](https://cov-lineages.org/resources/pangolin.html) or [Nextclade](https://clades.nextstrain.org/) to identify the viral clade or Pango lineage
- **Read quality control:** Run [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) on your raw FASTQ files before alignment to catch low-quality samples early
- **Aggregate QC report:** Use [MultiQC](https://multiqc.info/) to combine FastQC, samtools, and bcftools reports into a single HTML summary page

---

## License

MIT
