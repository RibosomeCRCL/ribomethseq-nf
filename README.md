# ribomethseq-nf workflow

`ribomethseq-nf` is a [nextflow](https://www.nextflow.io/) pipeline dedicated to
RiboMethSeq data processing. It generates counts and quality control data.

### Versions

 - 1.0 : first release of the workflow

## Software requirements

The following software program/packages are required to run this pipeline.

  - [nextflow](https://www.nextflow.io)
  - [bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
  - [samtools](https://github.com/samtools/samtools)
  - [bedtools](https://github.com/arq5x/bedtools2)
  - [Trimmomatic](https://github.com/usadellab/Trimmomatic)
  - [FastQC](https://github.com/s-andrews/FastQC)
  - [MultiQC](https://github.com/ewels/MultiQC)
  - [R](https://cran.r-project.org) and the following libraries :
    + reshape2
    + dplyr
    + tidyr
    + rmdformats
    + rmarkdown
    + ngsReports
    + ComplexHeatmap


If you do not want to install them manually, you can use either docker/singularity
or conda to build the required environment (see next section)

## Installation

To have the workflow installed you simply need to clone this repository.

```sh
git clone https://github.com/EpiRnaTools/ribomethseq-nf
```

### Conda environment

Perhaps the easiest way to have a proper environment to run the pipeline is to
use conda. A yaml file is provided here : `docker/conda.yml`

You can either directly use the conda   profile provided in `nextflow.config`
and the workflow will automatically build the environment at workflow
initiation or you can build it in advance with for instance the following
command:

```sh
cd docker
conda create -n ribomethseq-1.0 -f conda.yml
```

If you create the environment prior to running the workflow, you will then need
to adapt the `process.conda` directive from the conda profile in `nextflow.config`.

```
conda {
  conda.enabled = true
  process.conda = '/path/to/your/conda/envs/ribomethseq-1.0'
  includeConfig "nf-config/exec.config"
}
```

### Docker image

A docker image can also be built using the provided `docker/Dockerfile`. The
image will be built upon a `continuumio/miniconda3` base image and the previous
conda environment.

To build the docker image :

```sh
cd docker
docker build -t ribomethseq-nf:1.0 .
```

If you need to mount specific path(s) of your infrastructure, adapt the docker
profile in the configuration file (`containerOptions`).

### Singularity image

There is no proper singularity recipe provided at the moment, but in the meantime
you can convert the Docker image to a singularity image.

### HPC environment

By default in `nf-config/exec.config`, the executor is set to `slurm`. You may
need to adapt it to your own computing infrastructure (e.g. `pbs`, `LSF`, or other)

You can also set `params.scheduler = 'local'` if you plan to run locally on your
computer but then pay attention to control the queue size with `--qsize` (set to
20 by default) to avoid launching to many jobs depending on your hardware
capabilities.

## Running the workflow

Let's suppose that you have cloned the repository in the following directory:
`/path/to/ribomethseq-nf`

To run the workflow, you will need to specify both your execution environment
and the species (human or mouse) of interest. This is done here by combining
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

--bowtie_threads     INT    Threads for bowtie                         Optional (7)
--fastqc_threads     INT    Threads for fastqc                         Optional (2)
--trimmo_threads     INT    Threads for trimmomatic                    Optional (3)
--samtools_threads   INT    Threads for samtools                       Optional (4)

--scheduler          STR    Job scheduler                              Optional (slurm)
--qsize              INT    Max number of parallel jobs                Optional (20)
--outdir             DIR    Output directory                           Optional (.)
--logdir             DIR    Log directory                              Optional ($outdir)
--help               FLAG   Displays this help
```

## Output files

**TODO Theo**
