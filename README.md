# gcf/assemble_shortreads

[![GitHub Actions CI Status](https://github.com/gcf/assemble_shortreads/actions/workflows/nf-test.yml/badge.svg)](https://github.com/gcf/assemble_shortreads/actions/workflows/nf-test.yml)
[![GitHub Actions Linting Status](https://github.com/gcf/assemble_shortreads/actions/workflows/linting.yml/badge.svg)](https://github.com/gcf/assemble_shortreads/actions/workflows/linting.yml)[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.XXXXXXX)
[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)

[![Nextflow](https://img.shields.io/badge/version-%E2%89%A524.10.5-green?style=flat&logo=nextflow&logoColor=white&color=%230DC09D&link=https%3A%2F%2Fnextflow.io)](https://www.nextflow.io/)
[![nf-core template version](https://img.shields.io/badge/nf--core_template-3.3.2-green?style=flat&logo=nfcore&logoColor=white&color=%2324B064&link=https%3A%2F%2Fnf-co.re)](https://github.com/nf-core/tools/releases/tag/3.3.2)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Seqera Platform](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Seqera%20Platform-%234256e7)](https://cloud.seqera.io/launch?pipeline=https://github.com/gcf/assemble_shortreads)

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

