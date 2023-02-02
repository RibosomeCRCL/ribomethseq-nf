# ribomethseq-nf

`ribomethseq-nf` is a [nextflow](https://www.nextflow.io/) pipeline dedicated to RiboMethSeq data processing. It generates quality control data and counts, which can be directly used for further analyses using the rRMSAnalyzer package.

### Versions

 - 1.0 : first release

## Software requirements

[Nextflow](https://www.nextflow.io) (**21.04 or later**) is required to run this pipeline.

> **Warning**
> Older Nextflow versions before 20.10.0 will likely fail to run RiboMethSeq-nf.

The following software program/packages are also required to run this pipeline.

  - [bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
  - [samtools](https://github.com/samtools/samtools) [>=1.15]
  - [bedtools](https://github.com/arq5x/bedtools2)
  - [Trimmomatic](https://github.com/usadellab/Trimmomatic)
  - [FastQC](https://github.com/s-andrews/FastQC)
  - [MultiQC](https://github.com/ewels/MultiQC)
  - [pandoc](https://pandoc.org/)
  - [R](https://cran.r-project.org) and the following libraries :
    + ade4
    + dplyr
    + tidyr
    + tibble
    + rmarkdown
    + pheatmap


If you do not want to install them manually, you can use either docker/singularity
or conda to build the required environment (see Installation section)

Note: Earlier version of samtools may work as well but will retain a few more
SAM records due to an issue in expression handling in `samtools view` prior to
version 1.15 (See [this issue](https://github.com/samtools/samtools/issues/1476)).
In our case this is related to this filtering expression `-d 'NM' -e '![XS]'`given
to samtools. This is expected to be really minor though.

## General workflow description

The workflow is currently designed to process RiboMethSeq data generated in
single-end mode on Illumina sequencers. For each detected fastq file (one per
sample), quality control (FastQC) and reads trimming (Trimmomatic) steps are
launched. Trimmed reads are then aligned on rRNA reference sequences (human or mouse,
included in the pipeline) with Bowtie2 in end-to-end mode with `--sensitive -L 17` parameters.
Aligned files are then coordinate sorted and filtered with samtools to obtain uniquely mapped reads.
From these, we generate count (coverage) files using bedtools. A MultiQC report
gathering metrics from FastQC, Trimmomatic and Bowtie2 is provided as well as an
HTML report with custom QC metrics computed on the full dataset (all detected
fastq files).

## Installation

To have the workflow installed you simply need to clone this repository.

```sh
git clone https://github.com/EpiRnaTools/ribomethseq-nf
```

### Conda environment

Perhaps the easiest way to have a proper environment to run the pipeline is to
use conda. Two YAML files are provided : `docker/rms-processing.yml` and
`docker/rms-report.yml` corresponding to two distinct environments. The former
contains requirements for the processing part of the pipeline (samtools, bowtie2, ...)
and the latter is dedicated to the report generation. We chose to split the two
environment because solving the full conda environment could take a really long time.

You can either directly use the conda profile provided in `nextflow.config`
and the workflow will automatically build the environment at workflow
initiation or you can build it in advance (**recommended**) with for instance the
following command:

```sh
cd docker
conda env create -f rms-processing.yml
conda env create -f rms-report.yml
```

If you create the environment prior to running the workflow, you will then need
to adapt the `process.conda` directives from the conda profile in `nextflow.config`
as shown below :

```
conda {
	includeConfig "nf-config/exec.config"
	conda.enabled = true

	process {
		withName: 'fastqc|trim|bowtie2|filter|multiqc|counts' {
			conda = "/path/to/your/conda/envs/rms-processing"
		}
		withName: 'split|report' {
			conda = "/path/to/your/conda/envs/rms-report"
		}
	}
}
```

### Docker image

A docker image can also be built using the provided `docker/Dockerfile.prod`. The
image will be built upon a debian 11 (bullseye) base image.

To build the docker image :

```sh
cd docker
docker build -t ribomethseq-nf:1.0 -f Dockerfile.prod .
```

If you need to mount specific path(s) of your infrastructure, adapt the docker
profile in the configuration file (`containerOptions`).

More information on docker image content in the [docker README](docker/README.md).

### Singularity image

There is no proper singularity recipe provided at the moment, but in the meantime
you can convert the Docker image to a singularity image.

```sh
# first, create an archive of the docker image
docker save ribomethseq-nf:1.0 | gzip > ribomethseq-nf_1.0.tar.gz

# (... somewhere else ...)
# second, build a sif image from the archive
singularity build [--sandbox] ribomethseq-nf_1.0.sif docker-archive://ribomethseq-nf_1.0.tar.gz
```
If problems occur with the loop device, the option --sandbox can be used to build the singularity image.

### HPC environment

By default in `nf-config/exec.config`, the executor is set to `slurm`. You may
need to adapt it to your own computing infrastructure (e.g. `pbs`, `LSF`, or other)

You can also set `params.scheduler = 'local'` if you plan to run locally on your
computer but then pay attention to control the queue size with `--qsize` (set to
20 by default) to avoid launching to many jobs depending on your hardware
capabilities.

## Tests

Some tests are provided in the `tests` directory. Once you have set up your
software environment :

```sh
cd tests

# Test with tools installed in PATH
make test

# Test with conda profile
make test-conda

# Test with docker profile on local machine
make test-docker SCHEDULER=local

# Test with singularity profile on HPC with pbs
make test-singularity SCHEDULER=pbs

# Test with your own created profile (e.g. `custom_profile`)
make test PROFILE=custom_profile

# Clean the tests directory
make clean
```

## Running the workflow

Let's suppose that you have cloned the repository in the following directory:
`/path/to/ribomethseq-nf`

To run the workflow, you will need to specify both your execution environment
and your species (human or mouse) of interest. This is done here by combining
nextflow profiles.

`human` and `mouse` profiles are already available. They gather for each species
all the required reference data and the bowtie indexes to save you some time.

### Run the workflow for human data with docker
```sh
nextflow run /path/to/ribomethseq-nf -profile human,docker \
    --fqdir $FastqDir \
    --outdir $OutDir
```

### Run the workflow for mouse data with conda
```sh
nextflow run /path/to/ribomethseq-nf -profile mouse,conda \
    --fqdir $FastqDir \
    --outdir $OutDir
```

See `Input parameters` section for a description of all available parameters.

### Quickstart

Easiest way to get you started for the non-bioinformatician

  1. You just need nextflow and conda available on your system.
  2. Find out what is your job scheduler. e.g. `pbs`

Then:

```
nextflow EpiRnaTools/ribomethseq-nf -profile conda,human --scheduler 'pbs' --qsize 10 --fqdir '/path/to/fastq/files'
```

This will automatically retrieve the nextflow pipeline from GitHub, build the
required conda environment and finally process your data. Output files will be
located in current directory (default).

## Reference data used by the pipeline

The following reference rRNA are used :

### Human
- https://www.ncbi.nlm.nih.gov/nuccore/NR_046235 (18S, 28S and 5.8S)
- https://www.ncbi.nlm.nih.gov/nuccore/NR_023363.1 (5S)

### Mouse
- https://www.ncbi.nlm.nih.gov/nuccore/NR_030686.1 (5S)
- https://www.ncbi.nlm.nih.gov/nuccore/NR_003280.2 (5.8S)
- https://www.ncbi.nlm.nih.gov/nuccore/NR_003278.3 (18S)
- https://www.ncbi.nlm.nih.gov/nuccore/NR_003279.1 (28S)

The associated fasta sequences for human and mouse organisms are stored in the
`data/fasta` directory. It is possible to add other species, but it will be
necessary to add a new profile in `nextflow.config` file for them.

Precomputed `bowtie2` indexes are also already provided for both human and mouse
in the folder `data/bowtie`. If you need to use your own index, you can specify
it through the `--bowtie_index` parameter.

## Input parameters

```
                            +++++++++++++++++++++++++
                            +  ribomethseq-nf help  +
                            +++++++++++++++++++++++++

--fqdir              DIR    Fastq files location                       Required
--fastq_pattern      STR    Pattern for fastq file selection           Optional (.fastq.gz)

--adapters           FILE   (Trimmomatic) Path to illumina adapters    Optional ($baseDir/data/adapters/TruSeq3-SE.fa)
--leading            INT    (Trimmomatic) LEADING parameter            Optional (30)
--trailing           INT    (Trimmomatic) TRAILING parameter           Optional (30)
--slidingwindow      STR    (Trimmomatic) SLIDINGWINDOW parameter      Optional (4:15)
--avgqual            INT    (Trimmomatic) AVGQUAL parameter            Optional (30)
--minlen             INT    (Trimmomatic) MINLEN parameter             Optional (8)

--bowtie_index       FILE   (Bowtie) Path to index                     Optional ($baseDir/data/bowtie/human/human_index)
--bowtie_opts        STR    (Bowtie) additional options                Optional (--sensitive -L 17)

--samtools_opts      STR    (samtools) options to view                 Optional (--no-PG -h -u -d 'NM' -e '![XS]')

--bowtie_threads     INT    Threads for bowtie                         Optional (7)
--fastqc_threads     INT    Threads for fastqc                         Optional (2)
--trimmo_threads     INT    Threads for trimmomatic                    Optional (3)
--samtools_threads   INT    Threads for samtools                       Optional (4)

--split              FLAG   Split count files by RNA                   Optional (false)
--scheduler          STR    Job scheduler                              Optional (slurm)
--qsize              INT    Max number of parallel jobs                Optional (20)
--outdir             DIR    Output directory                           Optional (.)
--logdir             DIR    Log directory                              Optional ($outdir)
--help               FLAG   Displays this help
```

## Output files
### Overview 
Output directory tree, after running the pipeline on the test dataset : 
```
.
├── bowtie2
│   ├── logs
│   │   ├── sample1_R1_001.bowtie2.stats.log
│   │   ├── sample2_R1_001.bowtie2.stats.log
│   │   └── sample3_R1_001.bowtie2.stats.log
│   ├── sample1_R1_001.uniq.bam
│   ├── sample2_R1_001.uniq.bam
│   └── sample3_R1_001.uniq.bam
├── counts
│   ├── sample1_R1_001.3_counts.csv
│   ├── sample1_R1_001.5_counts.csv
│   ├── sample2_R1_001.3_counts.csv
│   ├── sample2_R1_001.5_counts.csv
│   ├── sample3_R1_001.3_counts.csv
│   └── sample3_R1_001.5_counts.csv
├── fastqc
│   ├── sample1_R1_001_fastqc.html
│   ├── sample2_R1_001_fastqc.html
│   ├── sample3_R1_001_fastqc.html
│   └── zip
│       ├── sample1_R1_001_fastqc.zip
│       ├── sample2_R1_001_fastqc.zip
│       └── sample3_R1_001_fastqc.zip
├── multiqc_report.html
├── rms_report.html
└── trimmomatic
    ├── logs
    │   ├── sample1_R1_001.trimmomatic.stats.log
    │   ├── sample2_R1_001.trimmomatic.stats.log
    │   └── sample3_R1_001.trimmomatic.stats.log
    ├── sample1_R1_001.trim.fastq.gz
    ├── sample2_R1_001.trim.fastq.gz
    └── sample3_R1_001.trim.fastq.gz

7 directories, 27 files
```

Five directories are generated : 
* Bowtie 2 : BAM alignment files and logs.
* Counts : 5' and 3' read-end count for each genomic position. (See [Read-end count files](#read-end-count-files)
 )
* Trimmomatic : Trimmed fastq files and logs.
* fastqc : fastqc report for each sample in both zip and html formats.

Two HTML reports are also generated : 

* multiqc_report.html : A RNA-seq report generated by MultiQC.
* rms_report.html : A RiboMethSeq report generated (See [RiboMethSeq (RMS) Report](#ribomethseq-rms-report))


### Read-end count files

The read-end count files represent the main output from this pipeline and are stored in the counts directory. **By default, only the 5'end-read counts are exported** (one file per sample). 3'end-read counts can also be exported alongside, using the `--threeandcount` parameter.

The files are in CSV format and have the following structure : 

| Ref RNA          | position on ref RNA | 5/3' read en count |
|------------------|---------------------|--------------------|
| NR_046235.3_5.8S | 1                   | 735                |
| NR_046235.3_5.8S | 2                   | 173                |
| NR_046235.3_5.8S | 3                   | 59                 |
| NR_046235.3_5.8S | 4                   | 32                 |
| NR_046235.3_5.8S | 5                   | 21                 |

The column headers have been added in the above example for clarity, but the real output files do not have them.

### RiboMethSeq (RMS) quality control Report

The RiboMethSeq QC report is stored in rms_report.html, at the root of the output directory.

It currently contains the following analyses :
* A end-read count boxplot and RLE for each samples.
* A distance heatmap to compare coverage profiles between samples.
* A correspondence analysis of the coverage profiles.

