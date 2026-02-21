# Ribosome Profiling Analysis Pipeline

Code for processing ribosome profiling (Ribo-seq) data and computing CELP-corrected ribosome footprint counts, as described in Wu et al.

## Overview

This pipeline processes raw Ribo-seq FASTQ files through adapter trimming, barcode demultiplexing, UMI handling, rRNA filtering, alignment, and deduplication. Downstream analysis uses [Ribolog](https://github.com/goodarzilab/Ribolog) to compute P-site offsets, codon-level read counts, and CELP (Consistent Excess of Loess Predictions) bias correction. CELP correction is performed separately for each experimental condition.

## Sequencing files

Sequencing files used for this analysis can be downloaded directly from the NCBI Gene Expression Omnibus under accession [#GSE213472](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE213472)

## Dependencies

### Command-line tools
- Trimmomatic (0.39)
- Cutadapt
- fastx_toolkit (fastx_barcode_splitter)
- UMI-tools
- STAR
- samtools

### R packages
- [Ribolog](https://github.com/goodarzilab/Ribolog) (see its README for installation via conda)
- ggplot2
- dplyr

## Reference files

The following reference files are required (not included):

| File | Description |
|------|-------------|
| Adapter FASTA | Linker sequence `AGATCGGAAGAGCAC` in TruSeq3-compatible format |
| Barcode file | Tab-delimited file mapping sample IDs to barcode sequences |
| rRNA FASTA | Complete human ribosomal repeating unit ([NCBI U13369](https://www.ncbi.nlm.nih.gov/nuccore/U13369)) |
| cDNA FASTA | Human cDNA sequences (Ensembl GRCh38), one transcript per gene with longest CDS |
| Ribolog annotation | `annotation.txt` generated per Ribolog documentation (Ensembl GRCh38) |
| Longest CDS file | `cDNA_longest_CDS.txt` generated per Ribolog documentation |

The annotation, longest CDS, and cDNA FASTA files are generated using the Python script provided in the Ribolog repository (`Biomart_cDNA_fasta_to_rW_annotation_and_reheadered_longest_CDS_cDNA_fasta_v2.py`) following the instructions in the Ribolog vignette.

**Note:** Alignments are performed against cDNA (not the genome), as required by Ribolog/riboWaltz for downstream P-site analysis.

## Pipeline

### 1. Adapter trimming

Remove linker sequence from RPF reads using Trimmomatic:

```bash
java -jar trimmomatic-0.39.jar SE \
  input.fastq.gz trimmed.fastq.gz \
  ILLUMINACLIP:adapters/TruSeq3-ribo-SE.fa:2:30:10 \
  LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:6
```

### 2. Barcode demultiplexing

Split reads by sample barcode:

```bash
gunzip -c trimmed.fastq.gz | fastx_barcode_splitter.pl \
  --bcfile barcode_ids.txt --prefix ./split/ --eol --mismatches 0
```

### 3. Barcode trimming

Remove barcode sequences from demultiplexed reads:

```bash
cutadapt -a BARCODE_SEQ -m 13 --cores=6 -o sample_trim.fastq.gz sample.fastq.gz
```

Replace `BARCODE_SEQ` with the appropriate barcode for each sample. The `-m 13` flag discards reads shorter than 13 nt after trimming.

### 4. UMI extraction

Extract UMIs (2 nt at 5' end, 5 nt at 3' end) and append to read headers:

```bash
umi_tools extract \
  --stdin sample_trim.fastq.gz \
  --extract-method=regex \
  --bc-pattern='^(?P<umi_1>.{2}).+(?P<umi_2>.{5})$' \
  --log=sample.log \
  --stdout sample_umi.fastq.gz
```

### 5. rRNA filtering

Remove ribosomal RNA reads using STAR aligned to the complete ribosomal repeating unit. Retain unmapped reads (non-rRNA):

```bash
STAR --genomeDir /path/to/rRNA_star_index \
  --readFilesIn sample_umi.fastq.gz \
  --runThreadN 6 \
  --readFilesCommand "gunzip -c" \
  --outFileNamePrefix sample_ribo_filter \
  --outReadsUnmapped Fastx
```

The `Unmapped.out.mate1` output contains non-rRNA reads for subsequent alignment.

### 6. Align to cDNA

Align filtered reads to the human cDNA reference using STAR:

```bash
STAR --genomeDir /path/to/cDNA_star_index \
  --readFilesIn sample_no_ribo.fastq \
  --runThreadN 6 \
  --outSAMtype BAM SortedByCoordinate \
  --outFileNamePrefix sample_cDNA_align
```

### 7. Index and deduplicate

```bash
samtools index sample.bam

umi_tools dedup \
  -I sample.bam \
  --output-stats=sample_stats \
  -S sample_dedup.bam
```

### 8. Downstream analysis (R)

See `ribolog_analysis.R` for the complete R workflow. In brief:

1. Read annotation and load deduplicated RPF BAMs into Ribolog via `bamtolist_rW`
2. Compute P-site offsets (`psite_rW`) and assign P-site positions (`psite_info_rW`)
3. Generate QC plots (read length distribution, periodicity)
4. Generate codon-level read counts via `psite_to_codon_count` (read lengths 20–32 nt)
5. Compute CELP bias correction via `CELP_bias`, **separately for each condition**
6. Save per-gene codon-level output files (one file per gene per sample)

**Note:** When using Ensembl cDNA references, `bamtolist_rW` is called with `refseq_sep = "."` to strip version suffixes from transcript IDs (e.g., `ENST00000000233.10` → `ENST00000000233`).

## Output

For each sample, a directory of per-gene text files is generated. Each file contains codon-level data with the following columns:

| Column | Description |
|--------|-------------|
| `codon_number` | Position along the CDS |
| `codon_type` | Trinucleotide codon |
| `aa_type` | Amino acid |
| `observed_count` | Raw P-site read count at this codon |
| `bias_coefficient` | CELP stalling bias coefficient |
| `corrected_count` | Bias-corrected read count |
