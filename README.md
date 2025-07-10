# De-novo Genome Assembly Pipeline (Nextflow)

This Nextflow pipeline performs **de-novo genome assembly** from **paired-end FASTQ files**, using a sequence of filtering, trimming, assembly, and quality control steps. The pipeline leverages the SPAdes assembler and supports optional consensus sequence generation when a reference genome is provided.

---

## 📦 Overview

This pipeline processes paired-end sequencing reads and outputs high-quality de-novo assemblies along with extensive quality metrics and reports. Key steps include:

- Adapter trimming and quality filtering
- Duplicate read removal
- Repairing of FASTQ pairing
- Genome assembly using SPAdes
- Quality assessment of assemblies with QUAST
- Quality checks using FastQC at each stage
- Aggregated reporting via MultiQC
- Consensus sequence generation (optional, requires reference genome)

---

## 🧬 Input Format

The input should be provided as a CSV file (e.g., `input.csv`) with the following structure:

```csv
sample,fastq_1,fastq_2
SRR1785701,/data/projects/snsb/de_nov_assembly_reads/e_coli/SRR1785701_1.fastq.gz,/data/projects/snsb/de_nov_assembly_reads/e_coli/SRR1785701_2.fastq.gz
SRR3584989,/data/projects/snsb/de_nov_assembly_reads/e_coli/SRR3584989_1.fastq.gz,/data/projects/snsb/de_nov_assembly_reads/e_coli/SRR3584989_2.fastq.gz


Each line represents a sample with its corresponding paired-end FASTQ files.

---

## ⚙️ Step-by-Step Pipeline Breakdown

1. **Adapter Trimming & Quality Filtering**
   - Tools: `bbmap` or `fastp`
   - Filters based on:
     - Adapter content
     - Minimum average base quality
     - Minimum read length

2. **Duplicate Read Removal**
   - Tool: `pardre`

3. **FASTQ Repair**
   - Tool: `bbmap`
   - Ensures paired-end reads are synchronized

4. **De-novo Assembly**
   - Tool: `SPAdes`
   - Assembles filtered reads into contigs

5. **Assembly Quality Metrics**
   - Tool: `QUAST`
   - Evaluates assembly statistics and contiguity

6. **Quality Control**
   - Tool: `FastQC`
   - Run **after every filtering step**

7. **Summary Reporting**
   - Tool: `MultiQC`
   - Aggregates all QC and assembly reports into one

8. **Consensus Sequence Generation** (Optional)
   - Tools: `bwamem2`, `samtools`
   - Requires a reference genome
   - Aligns reads to the reference and generates consensus sequences

---

## 🚀 How to Run

Run the pipeline using the following command:

```bash
nextflow assemble_shortreads/ \
  --input input.csv \
  --outdir ./test_run \
  -profile mamba \
  -resume

