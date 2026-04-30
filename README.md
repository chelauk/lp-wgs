## Introduction

**lp-wgs** is a Nextflow pipeline for shallow whole genome sequencing analysis. It can start from paired-end FASTQ files or aligned BAM files, then generate low-pass copy number and QC outputs from tools including FastQC, fastp, BWA, mosdepth, Picard, HMMcopy, ichorCNA, QDNAseq, ACE, ASCAT-style low-pass fitting, and optionally MEDICC2.

The pipeline uses Nextflow DSL2 and supports containerised execution with Docker, Singularity/Apptainer, Podman, Shifter, Charliecloud, or Conda.

## Pipeline Summary

When starting from FASTQ files (`--step mapping`), the pipeline runs:

1. Read QC with [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
2. Read trimming with [fastp](https://github.com/OpenGene/fastp)
3. Alignment with [BWA](https://github.com/lh3/bwa)
4. Lane merging and BAM indexing
5. Optional fragment-size BAM filtering with [samtools](https://www.htslib.org/)
6. Alignment and insert-size QC with [Picard](https://broadinstitute.github.io/picard/)
7. Coverage QC with [mosdepth](https://github.com/brentp/mosdepth)

When starting from BAM files (`--step calling`), or after mapping has completed, the pipeline runs copy number preparation and callers:

1. HMMcopy read counting, using a supplied GC wig or generating one with `--call_gc`
2. ichorCNA tumour fraction and copy number calling
3. QDNAseq segmentation and bin-level log2 ratios
4. ACE absolute copy number estimation
5. ASCAT-style low-pass purity, ploidy, and copy number fitting
6. Optional MEDICC2 preparation and execution
7. MultiQC reporting

The caller set is controlled with `--tools`, for example `--tools ace,ichor,ascat`.

## Quick Start

Run from paired-end FASTQ files:

```bash
nextflow run /path/to/lp_wgs \
    --input samplesheet_fastq.csv \
    --outdir results \
    --igenomes_base /path/to/igenomes \
    --step mapping \
    -profile singularity \
    -resume
```

Run from existing BAM files:

```bash
nextflow run /path/to/lp_wgs \
    --input samplesheet_bam.csv \
    --outdir results \
    --igenomes_base /path/to/igenomes \
    --step calling \
    -profile singularity \
    -resume
```

## Samplesheets

FASTQ input:

```csv
patient,sample,lane,fastq_1,fastq_2,predicted_ploidy
patient1,sample1,1,/data/patient1_sample1_L001_R1.fastq.gz,/data/patient1_sample1_L001_R2.fastq.gz,2
patient1,sample1,2,/data/patient1_sample1_L002_R1.fastq.gz,/data/patient1_sample1_L002_R2.fastq.gz,2
patient1,sample2,1,/data/patient1_sample2_L001_R1.fastq.gz,/data/patient1_sample2_L001_R2.fastq.gz,3
```

BAM input:

```csv
patient,sample,bam,bai,predicted_ploidy
patient1,sample1,/data/patient1_sample1.bam,/data/patient1_sample1.bam.bai,2
patient1,sample2,/data/patient1_sample2.bam,/data/patient1_sample2.bam.bai,3
```

`patient` and `sample` are required. `lane` is required for FASTQ input. `predicted_ploidy` is optional and defaults to `2`; it is used by ACE/ASCAT-related copy number steps.

## Common Parameters

| Parameter | Default | Description |
| --- | ---: | --- |
| `--step` | `mapping` | Use `mapping` for FASTQ input or `calling` for BAM input. |
| `--tools` | `ace,ichor,ascat` | Comma-separated copy number tools to run. |
| `--bin` | `1000` | Bin size in kb for HMMcopy/ichorCNA-style read counting. Supported values in the bundled config are `10`, `50`, `500`, and `1000`. |
| `--ploidy` | `2` | Default ploidy used by ASCAT-style fitting when sample-level ploidy is not supplied. |
| `--ascat_pcf_gamma` | `10` | Penalty passed to `copynumber::pcf()` for ASCAT low-pass segmentation. Higher values produce fewer segments. |
| `--filter_bam` | `false` | Filter BAMs by insert size before calling/QC. |
| `--filter_bam_min` | `90` | Minimum insert size when `--filter_bam` is enabled. |
| `--filter_bam_max` | `150` | Maximum insert size when `--filter_bam` is enabled. |
| `--ichor_purity` | unset | Set to `cf_dna` to use low-fraction ichorCNA defaults. |
| `--ichor_purity_manual` | unset | Manually pass ichorCNA normal fractions, for example `"c(0.95,0.99,0.995,0.999)"`. |
| `--qdnaseq_genome` | genome config | Genome label used by QDNAseq/ASCAT support code: `hg19`, `hg38`, or `mm10`. |

See [docs/usage.md](docs/usage.md) for more complete run examples and reference configuration notes.

## Reference Data

The default human configuration expects an iGenomes-style directory layout under `--igenomes_base`, including BWA indexes, FASTA, FASTA index, sequence dictionary, GC/mappability wig files, centromeres, and chromosome arm boundary files.

You can override individual reference paths directly, for example:

```bash
--fasta /refs/Homo_sapiens_assembly38.fasta \
--fasta_fai /refs/Homo_sapiens_assembly38.fasta.fai \
--dict /refs/Homo_sapiens_assembly38.dict \
--bwa /refs/BWAIndex \
--chr_arm_boundaries /refs/chrArmBoundaries_hg38.txt
```

## Output

Results are written under:

```text
<outdir>/<patient>/<patient>_<sample>/low_pass_wgs/
```

The main result areas are:

- `reports/`: FastQC, Picard, mosdepth, and other QC outputs.
- `ichorcna_<bin>/`: ichorCNA outputs for the chosen bin size.
- `ace/`: ACE outputs.
- `ascat/`: ASCAT low-pass copy number calls and plot PDFs.
- `bwa/`: mapped BAM/BAM index outputs when mapping is run.
- `<outdir>/reports/low_pass_wgs/`: MultiQC report and pipeline-level reporting files.

See [docs/output.md](docs/output.md) for details.

## Documentation

- [Usage](docs/usage.md)
- [Output](docs/output.md)

## Credits

lp-wgs was originally written by Chela James George Cresswell.

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

## Citations

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/master/LICENSE).
