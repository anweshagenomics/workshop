# 🧬 Hands-on Session on Antimicrobial Resistance: Genomic Perspective

> **Exploring AMR Databases and Genome-based Resistance Tools**  
> **Instructor:** Anwesha &nbsp;|&nbsp; **Date:** March 19, 2026 &nbsp;|&nbsp; **Duration:** 180 minutes

---

## 📖 Table of Contents

- [📋 Session Overview](#-session-overview)
- [🎯 Learning Objectives](#-learning-objectives)
- [📊 Sample Datasets](#-sample-datasets)
- [Prerequisites](#prerequisites)
- [🔧 Part 1: Environment Setup](#-part-1-environment-setup)
  - [1.1 Check System Information](#11-check-system-information)
  - [1.2 Install & Configure Conda](#12-install--configure-conda-if-not-already-done)
  - [1.3 Create Working Directories](#13-create-working-directories)
- [📥 Part 2: Download Sample Genomes](#-part-2-download-sample-genomes)
  - [2.1 Inspect Genome Files](#21-inspect-genome-files)
- [🧪 Exercise 1: AMRFinderPlus Analysis](#-exercise-1-amrfinderplus-analysis)
  - [1.1 Install AMRFinderPlus](#11-install-amrfinderplus)
  - [1.2 Update the Database](#12-update-the-amrfinderplus-database)
  - [1.3 Run AMRFinderPlus](#13-run-amrfinderplus)
  - [1.4 Inspect Results](#14-inspect-amrfinderplus-results)
- [🧪 Exercise 2: CARD RGI Analysis](#-exercise-2-card-rgi-analysis)
  - [2.1 Install RGI](#21-install-rgi)
  - [2.2 Run RGI](#22-run-rgi)
  - [2.3 Inspect Results](#23-inspect-rgi-results)
- [🧪 Exercise 3: ABRicate Multi-Database Screening](#-exercise-3-abricate-multi-database-screening)
  - [3.1 Install ABRicate](#31-install-abricate)
  - [3.2 Set Up Databases](#32-set-up-databases)
  - [3.3 Run ABRicate](#33-run-abricate-against-multiple-databases)
  - [3.4 Inspect Results](#34-inspect-abricate-results)
- [📊 Exercise 4: Comparative Analysis](#-exercise-4-comparative-analysis)
  - [4.1 Install Python Dependencies](#41-install-python-dependencies)
  - [4.2 Run the Comparison Script](#42-run-the-comparison-script)
- [💾 Export Results](#-export-results)
- [💭 Discussion Questions](#-discussion-questions)
- [🎓 Additional Exercises](#-additional-exercises-optional)
- [✅ Best Practices Summary](#-best-practices-summary)
- [📚 Additional Resources](#-additional-resources)
- [🎯 Session Summary](#-session-summary)

---

## 📋 Session Overview

This guide provides hands-on experience with four AMR analysis tools run entirely from your **local terminal** (Linux/macOS). No Colab or cloud environment is required.

| Tool | Description |
|------|-------------|
| **AMRFinderPlus** | NCBI's comprehensive AMR detection tool |
| **RGI (CARD)** | Resistance Gene Identifier from the CARD database |
| **ABRicate** | Multi-database screening tool |
| **Comparative Analysis** | Cross-validation across tools |

---

## 🎯 Learning Objectives

- Install and configure AMR detection tools via `conda`
- Download and prepare *Klebsiella pneumoniae* genomes from NCBI
- Execute AMR prediction using multiple tools
- Interpret resistance gene detection results
- Compare outputs across different databases
- Understand quality metrics (coverage, identity)

---

## 📊 Sample Datasets

Two *Klebsiella pneumoniae* isolates will be analyzed:

| Isolate | Accession | Strain |
|---------|-----------|--------|
| Isolate 1 | `GCF_051414815.1` | K. pneumoniae MVS2 |
| Isolate 2 | `GCF_051549635.1` | K. pneumoniae ASM5154963v1 |

---

## Prerequisites

- **conda** installed (Miniconda or Anaconda). If not installed, follow the [Miniconda installation guide](https://docs.conda.io/en/latest/miniconda.html).
- `wget`, `unzip`, `awk`, `tree` available on your system.
- At least **4 GB free disk space** for databases.

---

## 🔧 Part 1: Environment Setup

### 1.1 Check System Information

```bash
echo "=== System Information ==="
uname -a

echo "=== CPU Information ==="
lscpu | grep -E 'Model name|CPU\(s\)'

echo "=== Memory Information ==="
free -h

echo "=== Disk Space ==="
df -h ~
```

### 1.2 Install & Configure Conda (if not already done)

```bash
# Download Miniconda installer
wget -q https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh

# Install silently
bash miniconda.sh -b -p ~/miniconda
rm miniconda.sh

# Add conda to PATH
export PATH="$HOME/miniconda/bin:$PATH"
echo 'export PATH="$HOME/miniconda/bin:$PATH"' >> ~/.bashrc

# Restart your shell or source the profile
source ~/.bashrc

# Configure bioinformatics channels
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict

echo "✅ Conda ready:"
conda --version
```

### 1.3 Create Working Directories

```bash
mkdir -p ~/amr_analysis/genomes
mkdir -p ~/amr_analysis/results/amrfinder
mkdir -p ~/amr_analysis/results/rgi
mkdir -p ~/amr_analysis/results/abricate
mkdir -p ~/amr_analysis/comparison

cd ~/amr_analysis

echo "✅ Directory structure:"
tree -L 2 ~/amr_analysis 2>/dev/null || find ~/amr_analysis -type d
```

---

## 📥 Part 2: Download Sample Genomes

```bash
DOWNLOAD_DIR="$HOME/amr_analysis/genomes"

echo "📥 Downloading K. pneumoniae genomes from NCBI..."

# Isolate 1: K. pneumoniae MVS2
echo "Downloading Isolate 1: K. pneumoniae MVS2 (GCF_051414815.1)..."
wget -q --show-progress -P "$DOWNLOAD_DIR" \
  https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/051/414/815/GCF_051414815.1_KpMVS2/GCF_051414815.1_KpMVS2_genomic.fna.gz

# Isolate 2: K. pneumoniae ASM5154963v1
echo "Downloading Isolate 2: K. pneumoniae ASM5154963v1 (GCF_051549635.1)..."
wget -q --show-progress -P "$DOWNLOAD_DIR" \
  https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/051/549/635/GCF_051549635.1_ASM5154963v1/GCF_051549635.1_ASM5154963v1_genomic.fna.gz

# Decompress
echo "📦 Decompressing..."
gunzip "$DOWNLOAD_DIR"/*.gz

echo "✅ Download complete! Files:"
ls -lh "$DOWNLOAD_DIR"/*.fna
```

### 2.1 Inspect Genome Files

```bash
echo "=== Genome Statistics ==="

echo ""
echo "📊 K. pneumoniae MVS2:"
grep -c ">" ~/amr_analysis/genomes/GCF_051414815.1_KpMVS2_genomic.fna \
  | xargs -I {} echo "  Contigs: {}"
grep -v ">" ~/amr_analysis/genomes/GCF_051414815.1_KpMVS2_genomic.fna \
  | wc -c | awk '{print "  Genome size: ~" int($1/1000000) " Mb"}'

echo ""
echo "📊 K. pneumoniae ASM5154963v1:"
grep -c ">" ~/amr_analysis/genomes/GCF_051549635.1_ASM5154963v1_genomic.fna \
  | xargs -I {} echo "  Contigs: {}"
grep -v ">" ~/amr_analysis/genomes/GCF_051549635.1_ASM5154963v1_genomic.fna \
  | wc -c | awk '{print "  Genome size: ~" int($1/1000000) " Mb"}'

echo ""
echo "🔍 First 3 lines of MVS2 FASTA:"
head -n 3 ~/amr_analysis/genomes/GCF_051414815.1_KpMVS2_genomic.fna
```

---

## 🧪 Exercise 1: AMRFinderPlus Analysis

**AMRFinderPlus** is NCBI's flagship tool for AMR detection. It identifies:
- Acquired resistance genes
- Point mutations in chromosomal targets
- Virulence factors
- Stress response genes

### 1.1 Install AMRFinderPlus

```bash
# Create a dedicated conda environment
conda create -n amr -y

# Install AMRFinderPlus from bioconda
conda install -n amr -y -c bioconda ncbi-amrfinderplus

echo "✅ AMRFinderPlus installed:"
conda run -n amr amrfinder --version
```

### 1.2 Update the AMRFinderPlus Database

> ⚠️ **Critical step.** Always update the database before running analyses to capture the latest resistance mechanisms.

```bash
echo "📦 Updating AMRFinderPlus database (may take 2–3 minutes)..."
conda run -n amr amrfinder --update

echo "✅ Database update complete!"
conda run -n amr amrfinder --database_version 2>/dev/null || echo "Database ready"
```

### 1.3 Run AMRFinderPlus

```bash
# --- Isolate 1: K. pneumoniae MVS2 ---
echo "🔬 Running AMRFinderPlus on K. pneumoniae MVS2..."

conda run -n amr amrfinder \
  -n ~/amr_analysis/genomes/GCF_051414815.1_KpMVS2_genomic.fna \
  --organism Klebsiella_pneumoniae \
  --plus \
  --output ~/amr_analysis/results/amrfinder/kp_mvs2_amr.tsv

echo "✅ Analysis complete!"
wc -l < ~/amr_analysis/results/amrfinder/kp_mvs2_amr.tsv \
  | awk '{print $1-1 " resistance determinants found"}'


# --- Isolate 2: K. pneumoniae ASM5154963v1 ---
echo "🔬 Running AMRFinderPlus on K. pneumoniae ASM5154963v1..."

conda run -n amr amrfinder \
  -n ~/amr_analysis/genomes/GCF_051549635.1_ASM5154963v1_genomic.fna \
  --organism Klebsiella_pneumoniae \
  --plus \
  --output ~/amr_analysis/results/amrfinder/kp_asm_amr.tsv

echo "✅ Analysis complete!"
wc -l < ~/amr_analysis/results/amrfinder/kp_asm_amr.tsv \
  | awk '{print $1-1 " resistance determinants found"}'
```

### 1.4 Inspect AMRFinderPlus Results

```bash
echo "📋 Preview of AMRFinderPlus results (MVS2):"
head -n 5 ~/amr_analysis/results/amrfinder/kp_mvs2_amr.tsv

echo ""
echo "📊 Summary — Drug Classes detected:"
# Print column header then count by class (column 12 = Class in standard output)
awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) if($i=="Class") col=i}
            NR>1 && col{print $col}' \
  ~/amr_analysis/results/amrfinder/kp_mvs2_amr.tsv \
  | sort | uniq -c | sort -rn
```

---

## 🧪 Exercise 2: CARD RGI Analysis

**RGI (Resistance Gene Identifier)** uses the CARD database with three detection modes:

| Mode | Description |
|------|-------------|
| **Perfect** | 100% identity to reference |
| **Strict** | High similarity (>95% identity) |
| **Loose** | Divergent homologs — requires validation |

### 2.1 Install RGI

```bash
# Create a dedicated conda environment for RGI
conda create -n tool -y
conda install -n tool -y -c bioconda rgi

echo "✅ RGI installed:"
MPLBACKEND=Agg conda run -n tool rgi main --version
```

### 2.2 Run RGI

```bash
# --- Isolate 1: K. pneumoniae MVS2 ---
echo "🔬 Running CARD RGI on K. pneumoniae MVS2..."

MPLBACKEND=Agg conda run -n tool rgi main \
  --input_sequence ~/amr_analysis/genomes/GCF_051414815.1_KpMVS2_genomic.fna \
  --output_file ~/amr_analysis/results/rgi/kp_mvs2_rgi \
  --input_type contig \
  --clean

echo "✅ RGI complete!"
wc -l < ~/amr_analysis/results/rgi/kp_mvs2_rgi.txt \
  | awk '{print $1-1 " genes found"}'


# --- Isolate 2: K. pneumoniae ASM5154963v1 ---
echo "🔬 Running CARD RGI on K. pneumoniae ASM5154963v1..."

MPLBACKEND=Agg conda run -n tool rgi main \
  --input_sequence ~/amr_analysis/genomes/GCF_051549635.1_ASM5154963v1_genomic.fna \
  --output_file ~/amr_analysis/results/rgi/kp_asm_rgi \
  --input_type contig \
  --clean

echo "✅ RGI complete!"
wc -l < ~/amr_analysis/results/rgi/kp_asm_rgi.txt \
  | awk '{print $1-1 " genes found"}'
```

### 2.3 Inspect RGI Results

```bash
echo "📋 Preview of RGI results (MVS2):"
head -n 5 ~/amr_analysis/results/rgi/kp_mvs2_rgi.txt

echo ""
echo "🎯 Detection stringency breakdown (Cut_Off column):"
awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) if($i=="Cut_Off") col=i}
            NR>1 && col{print $col}' \
  ~/amr_analysis/results/rgi/kp_mvs2_rgi.txt \
  | sort | uniq -c | sort -rn

echo ""
echo "💡 Perfect = 100% identity | Strict = >95% | Loose = Divergent"
```

---

## 🧪 Exercise 3: ABRicate Multi-Database Screening

**ABRicate** screens genomes against multiple databases in a single command.

| Database | Contents |
|----------|----------|
| `card` | CARD resistance genes |
| `resfinder` | ResFinder resistance genes |
| `argannot` | ARG-ANNOT annotations |
| `ncbi` | NCBI AMR reference genes |
| `megares` | MEGARes resistance database |
| `vfdb` | Virulence Factor Database |

### 3.1 Install ABRicate

```bash
# Install into base conda environment
conda install -y -c bioconda abricate

echo "✅ ABRicate installed:"
abricate --version
```

### 3.2 Set Up Databases

```bash
echo "📦 Setting up ABRicate databases..."
abricate --setupdb

echo "📚 Available databases:"
abricate --list
```

### 3.3 Run ABRicate Against Multiple Databases

```bash
GENOME="$HOME/amr_analysis/genomes/GCF_051414815.1_KpMVS2_genomic.fna"
OUTDIR="$HOME/amr_analysis/results/abricate"

# CARD database
echo "🔬 ABRicate (CARD)..."
abricate --db card "$GENOME" > "$OUTDIR/mvs2_card.tab"
echo "  Genes found: $(( $(wc -l < "$OUTDIR/mvs2_card.tab") - 1 ))"

# ResFinder database
echo "🔬 ABRicate (ResFinder)..."
abricate --db resfinder "$GENOME" > "$OUTDIR/mvs2_resfinder.tab"
echo "  Genes found: $(( $(wc -l < "$OUTDIR/mvs2_resfinder.tab") - 1 ))"

# ARG-ANNOT database
echo "🔬 ABRicate (ARG-ANNOT)..."
abricate --db argannot "$GENOME" > "$OUTDIR/mvs2_argannot.tab"
echo "  Genes found: $(( $(wc -l < "$OUTDIR/mvs2_argannot.tab") - 1 ))"

# VFDB — virulence factors
echo "🔬 ABRicate (VFDB - Virulence)..."
abricate --db vfdb "$GENOME" > "$OUTDIR/mvs2_vfdb.tab"
echo "  Virulence factors found: $(( $(wc -l < "$OUTDIR/mvs2_vfdb.tab") - 1 ))"

echo "✅ All ABRicate analyses complete!"
```

### 3.4 Inspect ABRicate Results

```bash
echo "📋 CARD results (columns: GENE, PRODUCT, %COVERAGE, %IDENTITY):"
awk -F'\t' 'NR==1 || NR>1' ~/amr_analysis/results/abricate/mvs2_card.tab \
  | cut -f6,8,9,10 | column -t | head -20

echo ""
echo "📊 Database comparison:"
echo "  CARD:      $(( $(wc -l < ~/amr_analysis/results/abricate/mvs2_card.tab) - 1 )) genes"
echo "  ResFinder: $(( $(wc -l < ~/amr_analysis/results/abricate/mvs2_resfinder.tab) - 1 )) genes"
echo "  ARG-ANNOT: $(( $(wc -l < ~/amr_analysis/results/abricate/mvs2_argannot.tab) - 1 )) genes"
```

---

## 📊 Exercise 4: Comparative Analysis

Use the following Python script to compare results across all tools, compute consensus genes, and generate visualizations.

### 4.1 Install Python Dependencies

```bash
pip install pandas matplotlib seaborn openpyxl
```

### 4.2 Run the Comparison Script

Save the script below as `~/amr_analysis/compare_results.py` and run it with:

```bash
cd ~/amr_analysis
python compare_results.py
```

<details>
<summary><b>📄 compare_results.py — click to expand</b></summary>

```python
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from collections import Counter

# ── Load Results ───────────────────────────────────────────────────────────
df_mvs2      = pd.read_csv('results/amrfinder/kp_mvs2_amr.tsv', sep='\t')
df_rgi       = pd.read_csv('results/rgi/kp_mvs2_rgi.txt', sep='\t')
df_card_abr  = pd.read_csv('results/abricate/mvs2_card.tab', sep='\t')
df_resfinder = pd.read_csv('results/abricate/mvs2_resfinder.tab', sep='\t')

# ── Gene Sets ──────────────────────────────────────────────────────────────
amrfinder_genes  = set(df_mvs2['Element symbol'].str.lower().str.strip())
rgi_genes        = set(df_rgi['Best_Hit_ARO'].str.lower().str.strip())
card_abr_genes   = set(df_card_abr['GENE'].str.lower().str.strip())
resfinder_genes  = set(df_resfinder['GENE'].str.lower().str.strip())
all_genes        = amrfinder_genes | rgi_genes | card_abr_genes | resfinder_genes

print("=" * 60)
print("🔄 Cross-Tool Comparison — K. pneumoniae MVS2")
print("=" * 60)
print(f"  AMRFinderPlus : {len(amrfinder_genes)} genes")
print(f"  CARD RGI      : {len(rgi_genes)} genes")
print(f"  ABRicate-CARD : {len(card_abr_genes)} genes")
print(f"  ResFinder     : {len(resfinder_genes)} genes")
print(f"  Total unique  : {len(all_genes)} genes")

# Consensus genes (detected by 3+ tools)
gene_counts = Counter({
    g: sum([g in amrfinder_genes, g in rgi_genes,
            g in card_abr_genes, g in resfinder_genes])
    for g in all_genes
})
consensus_genes = [g for g, c in gene_counts.items() if c >= 3]

print(f"\n✅ High-confidence genes (3+ tools): {len(consensus_genes)}")
for gene in sorted(consensus_genes):
    print(f"  • {gene}")

# ── Quality Metrics ────────────────────────────────────────────────────────
print("\n📊 Quality Metrics (AMRFinderPlus):")
id_col  = '% Identity to reference sequence'
cov_col = '% Coverage of reference sequence'
if id_col in df_mvs2.columns:
    print(f"  Avg Identity : {df_mvs2[id_col].mean():.2f}%")
    print(f"  Avg Coverage : {df_mvs2[cov_col].mean():.2f}%")

# ── Resistance Classes ─────────────────────────────────────────────────────
print("\n🧬 Resistance Determinants by Class (AMRFinderPlus):")
if 'Class' in df_mvs2.columns:
    for cls, cnt in df_mvs2['Class'].value_counts().items():
        print(f"  {cls}: {cnt}")

# ── Plots ──────────────────────────────────────────────────────────────────
sns.set_style("whitegrid")

# Plot 1: Drug class bar chart
if 'Class' in df_mvs2.columns:
    fig, ax = plt.subplots(figsize=(12, 6))
    df_mvs2['Class'].value_counts().plot(kind='barh', ax=ax, color='steelblue')
    ax.set_xlabel('Number of Genes', fontsize=12)
    ax.set_title('AMR Genes by Drug Class (K. pneumoniae MVS2)', fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig('comparison/drug_classes.png', dpi=300)
    print("\n✅ Saved: comparison/drug_classes.png")

# Plot 2: Quality metrics scatter
if id_col in df_mvs2.columns:
    fig, ax = plt.subplots(figsize=(10, 8))
    sc = ax.scatter(df_mvs2[cov_col], df_mvs2[id_col],
                    c=df_mvs2[id_col], cmap='RdYlGn',
                    s=100, alpha=0.6, edgecolors='black', linewidth=0.5)
    ax.axhline(90, color='red', linestyle='--', alpha=0.5, label='90% Identity threshold')
    ax.axvline(80, color='blue', linestyle='--', alpha=0.5, label='80% Coverage threshold')
    ax.set_xlabel('Coverage (%)', fontsize=12)
    ax.set_ylabel('Identity (%)', fontsize=12)
    ax.set_title('AMR Gene Detection Quality Metrics', fontsize=14, fontweight='bold')
    ax.legend()
    plt.colorbar(sc, ax=ax, label='Identity (%)')
    plt.tight_layout()
    plt.savefig('comparison/quality_metrics.png', dpi=300)
    print("✅ Saved: comparison/quality_metrics.png")

# Plot 3: Tool comparison bar chart
tool_data = {
    'Tool':           ['AMRFinderPlus', 'CARD RGI', 'ABRicate (CARD)', 'ABRicate (ResFinder)'],
    'Genes Detected': [len(amrfinder_genes), len(rgi_genes),
                       len(card_abr_genes), len(resfinder_genes)]
}
fig, ax = plt.subplots(figsize=(10, 6))
bars = ax.bar(tool_data['Tool'], tool_data['Genes Detected'],
              color=['#3498db', '#e74c3c', '#2ecc71', '#f39c12'])
ax.set_ylabel('Number of Genes', fontsize=12)
ax.set_title('AMR Gene Detection: Tool Comparison', fontsize=14, fontweight='bold')
ax.grid(axis='y', alpha=0.3)
for bar in bars:
    h = bar.get_height()
    ax.text(bar.get_x() + bar.get_width()/2., h, f'{int(h)}',
            ha='center', va='bottom', fontsize=11, fontweight='bold')
plt.tight_layout()
plt.savefig('comparison/tool_comparison.png', dpi=300)
print("✅ Saved: comparison/tool_comparison.png")

# ── Excel Export ───────────────────────────────────────────────────────────
with pd.ExcelWriter('comparison/AMR_Results_Comparison.xlsx', engine='openpyxl') as writer:
    df_mvs2.to_excel(writer, sheet_name='AMRFinderPlus', index=False)
    df_rgi.to_excel(writer, sheet_name='CARD_RGI', index=False)
    df_card_abr.to_excel(writer, sheet_name='ABRicate_CARD', index=False)
    df_resfinder.to_excel(writer, sheet_name='ABRicate_ResFinder', index=False)
    pd.DataFrame({
        'Tool':           ['AMRFinderPlus', 'CARD RGI', 'ABRicate (CARD)', 'ABRicate (ResFinder)'],
        'Genes Detected': [len(amrfinder_genes), len(rgi_genes),
                           len(card_abr_genes), len(resfinder_genes)]
    }).to_excel(writer, sheet_name='Summary', index=False)

print("✅ Saved: comparison/AMR_Results_Comparison.xlsx")
```

</details>

---

## 💾 Export Results

### Generate a Text Summary Report

```bash
cat > ~/amr_analysis/comparison/SUMMARY_REPORT.txt << 'EOF'
================================================================================
AMR ANALYSIS SUMMARY REPORT
K. pneumoniae MVS2 (GCF_051414815.1)
================================================================================

Analysis performed using: AMRFinderPlus, CARD RGI, ABRicate (CARD + ResFinder)

--------------------------------------------------------------------------------
FILES
--------------------------------------------------------------------------------
  AMRFinderPlus : results/amrfinder/kp_mvs2_amr.tsv
  CARD RGI      : results/rgi/kp_mvs2_rgi.txt
  ABRicate-CARD : results/abricate/mvs2_card.tab
  ABRicate-ResF : results/abricate/mvs2_resfinder.tab
  Virulence     : results/abricate/mvs2_vfdb.tab

================================================================================
EOF

echo "✅ Summary report scaffold created. Run compare_results.py for full stats."
```

### Archive All Results

```bash
cd ~
zip -r AMR_Analysis_Results.zip amr_analysis/ -q

echo "✅ Archive created: AMR_Analysis_Results.zip"
echo "📁 Contents:"
unzip -l AMR_Analysis_Results.zip | head -30
```

---

## 💭 Discussion Questions

1. **Which resistance genes were consistently detected across all tools?**
   - What does this tell you about database concordance?

2. **Were there any database-specific calls?**
   - How would you investigate these discrepancies?
   - Which database would you trust more for clinical reporting?

3. **Quality Metrics:**
   - How many genes had ≥90% identity and ≥80% coverage?
   - What would you do with low-coverage or low-identity hits?

4. **What major antibiotic classes show resistance?**
   - β-lactams (penicillins, cephalosporins, carbapenems)?
   - Aminoglycosides? Fluoroquinolones?
   - What are the clinical implications?

5. **If this were a clinical sample:**
   - Which antibiotics would you predict as ineffective?
   - What treatment options might remain?

6. **Tool Selection:**
   - For routine surveillance, which tool would you choose and why?
   - For research, would your choice differ?

---

## 🎓 Additional Exercises (Optional)

### Challenge 1: Analyze the Second Sample

Repeat the full analysis for `GCF_051549635.1_ASM5154963v1` and compare findings with MVS2.

### Challenge 2: Virulence Factor Analysis

```bash
# View VFDB results
cat ~/amr_analysis/results/abricate/mvs2_vfdb.tab | column -t | head -20
```

### Challenge 3: Custom Visualizations

Ideas to try in Python:
- Venn diagram showing gene overlap between tools
- Heatmap of resistance genes × samples
- Network analysis of co-occurring resistance genes

### Challenge 4: Export for Phylogenetic Analysis

Prepare your data for integration with MLST, SNP calling, or pangenome analysis pipelines.

---

## ✅ Best Practices Summary

### DO ✓
- Always update databases before analysis
- Use multiple tools for cross-validation
- Set appropriate quality thresholds (≥90% ID, ≥80% coverage)
- Correlate with phenotypic AST when available
- Document software versions and database dates
- Consider local epidemiology when interpreting results

### DON'T ✗
- Trust low-confidence predictions without validation
- Rely on a single tool/database
- Ignore partial gene matches without investigation
- Make clinical decisions based solely on genomic predictions
- Forget to check for novel or divergent variants

### Key Principles

| Principle | Detail |
|-----------|--------|
| **Database Currency** | Regular updates capture new resistance mechanisms |
| **Cross-Validation** | Multiple tools reduce false positives/negatives |
| **Quality Control** | Coverage and identity metrics indicate confidence |
| **Clinical Context** | Interpret findings within epidemiological framework |
| **Reproducibility** | Document everything for traceable analyses |

---

## 📚 Additional Resources

### Databases
- [CARD](https://card.mcmaster.ca)
- [NCBI AMRFinderPlus](https://www.ncbi.nlm.nih.gov/pathogens/antimicrobial-resistance/)
- [ResFinder](https://cge.food.dtu.dk/services/ResFinder/)
- [VFDB](http://www.mgc.ac.cn/VFs/)

### Web Tools
- [PathogenWatch](https://pathogen.watch)
- [BV-BRC](https://www.bv-brc.org)
- [Virulent2](https://bioinfo.icgeb.res.in/virulent2/)

### Key Publications
- CARD 2023: Alcock et al., *Nucleic Acids Res.* 2023
- AMRFinderPlus: Feldgarden et al., *Antimicrob. Agents Chemother.* 2021
- ResFinder 4.0: Bortolaia et al., *J. Antimicrob. Chemother.* 2020

### Tool Documentation
- [AMRFinderPlus Wiki](https://github.com/ncbi/amr/wiki)
- [CARD RGI](https://github.com/arpcard/rgi)
- [ABRicate](https://github.com/tseemann/abricate)

---

## 🎯 Session Summary

### What We Accomplished

- [x] Set up a complete AMR analysis environment using conda
- [x] Downloaded and inspected *Klebsiella pneumoniae* genomes from NCBI
- [x] Ran AMRFinderPlus for comprehensive resistance detection
- [x] Used CARD RGI with multiple detection stringencies
- [x] Performed multi-database screening with ABRicate
- [x] Compared results across tools and generated visualizations
- [x] Identified virulence factors alongside resistance genes

### Key Takeaways

1. **No single tool is perfect** — use multiple databases for validation
2. **Quality metrics matter** — coverage and identity guide confidence
3. **Context is critical** — genomic predictions need clinical correlation
4. **Regular updates essential** — AMR databases evolve constantly
5. **Reproducibility** — document everything for traceable science

---

*Workshop material: Hands-on Session on Antimicrobial Resistance: Genomic Perspective*  
*March 19, 2026 | Instructor: Anwesha*
