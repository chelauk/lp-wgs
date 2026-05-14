# lp-wgs: Output

This page describes the main outputs produced by the pipeline. Paths are relative to `--outdir`.

Most sample-level outputs are published under:

```text
<patient>/<patient>_<sample>/low_pass_wgs/
```

Pipeline-level reports are published under:

```text
reports/low_pass_wgs/
```

## Mapping and QC Outputs

These outputs are produced when running `--step mapping`. Some QC outputs are also produced from BAM input where applicable.

### FastQC

```text
<patient>/<patient>_<sample>/low_pass_wgs/reports/fastqc/
```

Files:

- `*_fastqc.html`: per-read quality report.
- `*_fastqc.zip`: FastQC data archive used by MultiQC.

### BWA BAMs

```text
<patient>/<patient>_<sample>/low_pass_wgs/bwa/
```

Files:

- `*.bam`: mapped, merged BAM for the sample.
- `*.bam.bai`: BAM index.

When multiple lanes are supplied for the same `patient` and `sample`, the pipeline aligns each lane and then merges them to one sample-level BAM.

### Picard Metrics

```text
<patient>/<patient>_<sample>/low_pass_wgs/reports/picard/
```

Files:

- `*metrics`: alignment summary metrics.
- `*.txt`: insert size metrics.
- `*.pdf`: insert size histogram.

### mosdepth

```text
<patient>/<patient>_<sample>/low_pass_wgs/reports/mosdepth/
```

Files can include:

- `*.summary.txt`: coverage summary.
- `*.global.dist.txt`: global coverage distribution.
- `*.regions.bed.gz`: region-level coverage, when regions are configured.
- `*.regions.bed.gz.csi`: index for region-level coverage.

## Copy Number Outputs

The copy number callers are not redundant wrappers around the same algorithm. They provide complementary views of the same low-pass WGS read-depth signal:

- ichorCNA uses hidden Markov modelling.
- ACE uses circular binary segmentation.
- ASCATlp uses piecewise constant fitting.
- BayesCNA uses Bayesian inference.

Running several callers can be useful when the goal is to compare method agreement, inspect caller-specific assumptions, or prioritise robust copy number events across approaches.

### HMMcopy Read Counts

HMMcopy read-counter wig files are primarily intermediate files consumed by ichorCNA. They are collected into the workflow and may be visible in work directories, but the main published copy number outputs come from the caller directories below.

### ichorCNA

```text
<patient>/<patient>_<sample>/low_pass_wgs/ichorcna_<bin>/
```

or, when `--ichor_purity_manual` is set:

```text
<patient>/<patient>_<sample>/low_pass_wgs/ichorcna_<bin>/ichor_purity_manual/
```

Files:

- `filter*`: ichorCNA result files produced by the module, including copy number calls, parameters, and plots depending on ichorCNA output mode.

The `<bin>` suffix reflects `--bin`, for example `ichorcna_1000`. ichorCNA is the hidden Markov model caller in the pipeline.

### QDNAseq Preparation

QDNAseq produces bin-level and segment-level files used by ASCAT low-pass
fitting, plus a reusable QDNAseq RDS object used by ACE.

Key intermediate files:

- `*bins.txt`: bin-level log2 ratio data.
- `*cna_segments.txt`: QDNAseq segment file.
- `*kbp.rds`: normalized segmented QDNAseq object.
- `*.pdf`: QDNAseq plot outputs.

These are emitted by `RUN_QDNASEQ` and passed directly into `RUN_ASCAT` and
`ACE`; the main published ASCAT-facing outputs are listed in the ASCAT section
below.

### ASCAT Low-Pass

```text
<patient>/<patient>_<sample>/low_pass_wgs/ascat/ascat_ploidy_<ploidy>/
```

Files:

- `*_ascat_lp_plot.pdf`: low-pass copy number plot with segment means and integer calls.
- `*_cna_ploidy_search_calls.txt`: integer copy number calls from ASCAT-style low-pass fitting.

ASCATlp is the piecewise constant fitting caller in the pipeline. It creates one `ascat_ploidy_<ploidy>` folder per value in `--ploidy`. Segmentation uses `copynumber::pcf()` with `--ascat_pcf_gamma` controlling the breakpoint penalty.

### ACE

```text
<patient>/<patient>_<sample>/low_pass_wgs/ace/
```

Files:

- `*filter_*`: ACE output directory/files for the filtered or unfiltered input state.
- `*_sky_on_fire.pdf`: ACE square-model matrix plot.
- `*_sqmodel_minmadf.txt`: ACE square-model minima table.

ACE is the circular binary segmentation caller in the pipeline.

### BayesCNA

BayesCNA is the Bayesian inference caller in the pipeline. Its outputs are published under the BayesCNA-specific result directory when `bayes_cna` is included in `--tools`.

### MEDICC2

When `medicc` is included in `--tools`, the pipeline prepares MEDICC2 input from copy number calls and runs MEDICC2 at patient level.

Main outputs:

- MEDICC2 input TSV files from `PREP_MEDICC2`.
- `medicc2_output/`: MEDICC2 output directory.

MEDICC is currently treated as human-specific by the workflow and is blocked for `mm10` runs.

## MultiQC

```text
reports/low_pass_wgs/
```

Files:

- `multiqc_report.html`: interactive QC report.
- `multiqc_data/`: parsed MultiQC data tables.
- `multiqc_plots/`: exported plots, when generated.

MultiQC summarises supported QC outputs such as FastQC, fastp, Picard, mosdepth, and software versions.

## Pipeline Information

Depending on the run configuration, Nextflow report files are written under:

```text
reports/pipeline_info/
```

Files can include:

- `execution_report_*.html`: resource and task execution summary.
- `execution_timeline_*.html`: timeline of tasks.
- `execution_trace_*.txt`: tabular task trace.
- `pipeline_dag_*.html`: workflow graph.
- `software_versions.yml`: collected tool versions where emitted by modules.

These files are useful for audit trails, troubleshooting, and reproducing runs.
