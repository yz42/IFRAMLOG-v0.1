# IFRAMLOG
 IFRAMLOG: Inference of ML Tree with a outgroup from SNPs Pipeline

# IFRAMLOG: Inference of ML Tree with a outgroup from SNPs Pipeline

### **Author:**

Yangzi Wang ([wangyz.benniao@gmail.com](mailto\:wangyz.benniao@gmail.com))

---

### **Introduction**

The IFRAMLOG (InFeRence A Maximum Likelihood tree with OutGroup from SNPs) pipeline is designed for generating a Maximum Likelihood (ML) phylogenetic tree with a outgroup from SNPs stored in VCF files. The pipeline extracts flanking sequences around SNPs, identifies orthologous alleles in an outgroup reference genome, and uses ModelTest-NG and RAxML-NG to generate an ML tree. This pipeline was originally used in the great duckweed population genomics project: https://doi.org/10.1038/s42003-024-06266-7

---

### **Key Features**

1. Extracts flanking 150 bp (totaling 301 bp, or user customized) regions around SNPs.
2. Removes problematic SNPs (e.g., multi-allelic or sequences with too many ‘N’s).
3. Uses BLAST to map flanking regions to an outgroup genome.
4. Identifies orthologous alleles in the outgroup species.
5. Performs model selection using ModelTest-NG.
6. Constructs ML phylogenetic trees with RAxML-NG.

---

### **Workflow Overview**

1. **Input Data:** VCF file with SNPs, a reference genome FASTA, and an outgroup genome FASTA.
2. **Flanking Sequence Extraction:** Extracts 150 bp around each SNP from the reference genome.
3. **BLAST Search:** Maps flanking sequences to the outgroup genome to identify orthologous alleles.
4. **Model Selection:** Uses ModelTest-NG to identify the best substitution model for phylogenetic analysis.
5. **ML Tree Construction:** Constructs an ML phylogenetic tree using RAxML-NG.

---

### **Usage**

```bash
perl IFRAMLOG_v0.1.pl --vcf <FILE> --ref <FILE> --out_dir <PATH> --out <FILE> [options]
```

#### Required Arguments:

- `--vcf <FILE>`: Input VCF file containing SNP data.
- `--ref <FILE>`: Reference genome FASTA file.
- `--OG_ref <FILE>`: Outgroup reference genome FASTA file.
- `--out_dir <PATH>`: Output directory to store results.
- `--out <FILE>`: Basename for output files.

#### Optional Arguments:

- `--blast_C <INT>`: Number of CPUs to use for BLAST and RAxML (default: 4).
- `--FLK_len <INT>`: Length of flanking sequences to extract on each side of the SNP (default: 150).
- `--N_perc <FLOAT>`: Maximum fraction of ‘N’ allowed in the flanking region (default: 0.5).
- `--help`: Print help message.

#### Example Command:

```bash
perl IFRAMLOG_v0.1.pl \
    --vcf input.vcf \
    --ref genome.fa \
    --OG_ref OG_genome.fa \
    --out_dir /path/to/output/ \
    --out snp_regions.fa \
    --blast_C 4 \
    --FLK_len 150 \
    --N_perc 0.5
```

---

### **Prerequisites**

Ensure the following dependencies are installed and accessible from the command line:

#### **Perl Modules:**

- `Bio::DB::Fasta`
  - Install via CPAN:
    ```bash
    cpanm Bio::DB::Fasta
    ```

#### **BLAST Tools:**

- `makeblastdb`
- `blastn`
  - Install via `conda`:
    ```bash
    conda install -c bioconda blast
    ```

#### **ModelTest and RAxML Tools:**

- `modeltest-ng`
  - Install via `conda`:
    ```bash
    conda install -c bioconda modeltest-ng
    ```
- `raxml-ng`
  - Install via `conda`:
    ```bash
    conda install -c bioconda raxml-ng
    ```

---

### **Output Files**

#### Flanking Sequence Extraction:

- `<output_dir>/<out>.fa`: FASTA file of extracted flanking sequences.
- `<output_dir>/<out>Excluded_SNP.log`: Log file of excluded SNPs with reasons.

#### BLAST and Parsing:

- `<output_dir>/<out>.blast.fmt0.res`: Raw BLAST results.
- `<output_dir>/<out>.blast.fmt0.res.final.4outgroup.tsv`: Table of query vs outgroup genotypes.
- `<output_dir>/<out>.blast.fmt0.res.outgroup.fa`: FASTA file of outgroup sequences.

#### Model Selection and ML Tree:

- `<output_dir>/<out>.MLtree.MDtest`: ModelTest-NG output.
- `<output_dir>/<out>.MLtree.raxml.supportFBP`: RAxML tree with Fast Bootstrap Proportions.
- `<output_dir>/<out>.MLtree.raxml.supportTBE`: RAxML tree with Transfer Bootstrap Expectations.

---

### **Pipeline Steps in Detail**

#### **Step 1: Flanking Sequence Extraction**

- Extracts 150 bp upstream and downstream of each SNP, resulting in 301 bp sequences.
- Filters SNPs based on:
  - Multi-allelic SNPs.
  - Chromosome regions near edges.
  - Too many ‘N’s in the flanking region.

#### **Step 2: BLAST Search**

- Uses `blastn` to map flanking sequences to the outgroup reference genome.
- Generates a table of orthologous alleles.

#### **Step 3: Parse BLAST Results**

- Identifies orthologous alleles for each SNP in the outgroup.
- Outputs FASTA and TSV files for downstream analysis.

#### **Step 4: Model Selection**

- Runs ModelTest-NG to identify the best substitution model using the BIC criterion.

#### **Step 5: ML Tree Construction**

- Constructs the ML tree using RAxML-NG.
- Outputs tree files with bootstrap support metrics (`FBP` and `TBE`).

---

### **Error Handling**

- Ensures required files are present before execution.
- Provides detailed error messages for missing files, BLAST failures, or pipeline issues.

---

### **Contact**

For any issues or questions, please contact **Yangzi Wang** at [wangyz.benniao@gmail.com](mailto\:wangyz.benniao@gmail.com).

