# lp-wgs: Usage

## Running the Pipeline

The pipeline can start from paired-end FASTQ files or from existing BAM files.

FASTQ input:

```bash
nextflow run /path/to/lp_wgs \
    --input samplesheet_fastq.csv \
    --outdir results \
    --igenomes_base /path/to/igenomes \
    --step mapping \
    -profile singularity \
    -resume
```

BAM input:

```bash
nextflow run /path/to/lp_wgs \
    --input samplesheet_bam.csv \
    --outdir results \
    --igenomes_base /path/to/igenomes \
    --step calling \
    -profile singularity \
    -resume
```

Use Docker locally by replacing the profile:

```bash
-profile docker
```

Use a cluster profile or local config file with `-c` when your institution requires scheduler-specific settings:

```bash
nextflow run /path/to/lp_wgs \
    --input samplesheet_bam.csv \
    --outdir results \
    --step calling \
    -profile singularity \
    -c local.config \
    -resume
```

## Input Samplesheets

The samplesheet is a CSV file. `patient` and `sample` are always required. File paths must not contain spaces.

### FASTQ Input

Use FASTQ input with `--step mapping`.

```csv
patient,sample,lane,fastq_1,fastq_2,predicted_ploidy
patient1,sample1,1,/data/patient1_sample1_L001_R1.fastq.gz,/data/patient1_sample1_L001_R2.fastq.gz,2
patient1,sample1,2,/data/patient1_sample1_L002_R1.fastq.gz,/data/patient1_sample1_L002_R2.fastq.gz,2
patient1,sample2,1,/data/patient1_sample2_L001_R1.fastq.gz,/data/patient1_sample2_L001_R2.fastq.gz,3
```

Required columns:

| Column | Description |
| --- | --- |
| `patient` | Patient or case identifier. No spaces. |
| `sample` | Sample identifier within patient. No spaces. |
| `lane` | Sequencing lane or run identifier. Multiple rows with the same `patient` and `sample` are merged after alignment. |
| `fastq_1` | Read 1 FASTQ path. Must end in `.fastq.gz` or `.fq.gz`. |
| `fastq_2` | Read 2 FASTQ path. Must end in `.fastq.gz` or `.fq.gz`. |

Optional columns:

| Column | Description |
| --- | --- |
| `predicted_ploidy` | Integer ploidy prior for copy number tools. Defaults to `2` if absent. |

### BAM Input

Use BAM input with `--step calling`.

```csv
patient,sample,bam,bai,predicted_ploidy
patient1,sample1,/data/patient1_sample1.bam,/data/patient1_sample1.bam.bai,2
patient1,sample2,/data/patient1_sample2.bam,/data/patient1_sample2.bam.bai,3
```

Required columns:

| Column | Description |
| --- | --- |
| `patient` | Patient or case identifier. No spaces. |
| `sample` | Sample identifier within patient. No spaces. |
| `bam` | Coordinate-sorted BAM file. Must end in `.bam`. |
| `bai` | BAM index. Must end in `.bai`. |

Optional columns:

| Column | Description |
| --- | --- |
| `predicted_ploidy` | Integer ploidy prior for copy number tools. Defaults to `2` if absent. |

## Choosing Tools

The `--tools` parameter controls downstream copy number tools. It is a comma-separated list.

Default:

```bash
--tools ace,ichor,ascat
```

Examples:

```bash
--tools ichor
--tools ace,ascat
--tools ichor,ascat,medicc
```

For mouse genomes (`qdnaseq_genome = mm10`), MEDICC is currently blocked by the workflow because it remains human-specific in this pipeline.

## Bin Size

Set the bin size in kilobases with `--bin`.

```bash
--bin 1000
```

The bundled HMMcopy settings currently map these values to read-counter windows:

| `--bin` | Window |
| ---: | ---: |
| `10` | 10,000 bp |
| `50` | 50,000 bp |
| `500` | 500,000 bp |
| `1000` | 1,000,000 bp |

## BAM Filtering

Set `--filter_bam` to filter alignments by insert size before analysis.

```bash
--filter_bam \
--filter_bam_min 90 \
--filter_bam_max 150
```

The defaults are `90` and `150`. For Illumina data, the filter is applied with the standard samtools view module. For Nanopore mode, the pipeline uses the Nanopore-specific filtering module.

## ichorCNA Purity Settings

For low tumour fraction or cfDNA-style data, use:

```bash
--ichor_purity cf_dna
```

For human genomes this applies ichorCNA settings equivalent to:

```text
--normal "c(0.95, 0.99, 0.995, 0.999)"
--ploidy "c(2)"
--maxCN 3
--estimateScPrevalence FALSE
--scStates "c()"
--chrs 'paste0("chr",c(1:22))'
--chrNormalize 'paste0("chr", c(1:22))'
--chrTrain 'paste0("chr", c(1:22))'
```

To provide normal fractions manually:

```bash
--ichor_purity_manual "c(0.95,0.99,0.995,0.999)"
```

## ASCAT Low-Pass Segmentation

ASCAT low-pass fitting uses QDNAseq-style bin log2 ratios and `copynumber::pcf()` segmentation before purity/ploidy fitting.

The PCF penalty is controlled by:

```bash
--ascat_pcf_gamma 10
```

The default is `10`. Higher values penalise new breakpoints more strongly and generally produce fewer segments. Lower values allow more segments.

Example:

```bash
nextflow run /path/to/lp_wgs \
    --input samplesheet_bam.csv \
    --outdir results \
    --step calling \
    --tools ascat \
    --ascat_pcf_gamma 20 \
    -profile singularity \
    -resume
```

## Reference Configuration

The pipeline can resolve references from `conf/igenomes.config` using `--genome` and `--igenomes_base`.

Default:

```bash
--genome GATK.GRCh38
```

Common reference paths include:

| Parameter | Purpose |
| --- | --- |
| `--igenomes_base` | Base directory for bundled genome path templates. |
| `--fasta` | Reference FASTA. |
| `--fasta_fai` | FASTA index. |
| `--dict` | Picard/GATK sequence dictionary. |
| `--bwa` | BWA index directory. |
| `--gc_wig` | GC wig file for ichorCNA/HMMcopy. |
| `--map_wig_file` | Mappability wig file. |
| `--chr_bed` | Chromosome BED used by mosdepth on supported genomes. |
| `--centromere` | Centromere file for ichorCNA. |
| `--chr_arm_boundaries` | Chromosome arm boundaries for ASCAT segmentation. |
| `--qdnaseq_genome` | Genome label for QDNAseq/ASCAT support code: `hg19`, `hg38`, or `mm10`. |
| `--qdnaseq_package` | QDNAseq annotation package name. |
| `--hmmcopy_chromosomes` | Explicit chromosome list for HMMcopy modules. |

You can either provide `--igenomes_base` with a matching directory layout, or override individual paths.

## Generating GC Wig Files

If `--call_gc` is set, the pipeline runs HMMcopy `gccounter` using the supplied FASTA and bin size instead of using `--gc_wig`.

```bash
--call_gc true
```

## Profiles and Containers

Use one of the bundled software profiles:

```bash
-profile docker
-profile singularity
-profile podman
-profile shifter
-profile charliecloud
-profile conda
```

Docker or Singularity/Apptainer are preferred for reproducibility. Conda is supported but should usually be a fallback.

Multiple profiles can be combined:

```bash
-profile test,docker
```

## Resuming Runs

Use `-resume` when restarting a run:

```bash
nextflow run /path/to/lp_wgs ... -resume
```

Nextflow reuses cached work where inputs, commands, and relevant parameters are unchanged.

## Running on SLURM

Example submission wrapper:

```bash
#!/bin/bash -l
#SBATCH --job-name=lp_wgs
#SBATCH --output=nextflow_%j.out
#SBATCH --ntasks=1
#SBATCH --time=96:00:00

module load java

nextflow run /path/to/lp_wgs \
    --input samplesheet_bam.csv \
    --outdir results \
    --igenomes_base /path/to/igenomes \
    --step calling \
    -profile singularity \
    -c local.config \
    -resume
```

The Nextflow process should remain alive for the duration of the workflow. On clusters, run Nextflow inside an interactive session, `screen`/`tmux`, or a scheduler job like the example above.

## Nextflow Runtime Files

The launch directory will contain:

```text
work/              Nextflow working directory
.nextflow.log      Current Nextflow log
.nextflow/         Nextflow run metadata
<outdir>/          Published results
```

Keep `work/` if you want `-resume` to work.
