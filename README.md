# RiboMethSeq-pipe

A Nextflow pipeline dedicated to RiboMethSeq counting and Quality Control. Its outputs are the 3/5 end count for each position of your samples (one CSV file per sample).

## Quickstart

There is no need to install the pipeline before running it. Nextflow will automatically download this repository.

`nextflow EPIRNATools/ribomethseq-nf --results_path "/path/to/result/" --dir "/path/to/data/*.fastq.gz" -profile conda,human`

`conda` can be replaced with `docker` (soon!) or `singularity` (soon!) depending on your configuration.

The `human` profile, for using human reference data, can be replaced by `mouse` if your samples are from this species.
## Installation

This is not necessary for simply running the pipeline.

Download this repository with : <br>
`git clone https://github.com/EPIRNAtools/ribomethseq-nf`

Then `cd ribomethseq-nf` and you are ready to run the pipeline !

__Please note that you have to replace `EPIRNATools/ribomethseq-nf` in the quickstart examples with `main.nf` in case of a manual installation.__ Otherwise, Nextflow will download the Github version.

## Parameters 

The following parameters MUST be specified before running this pipeline:

**--dir** : the path to your input fastq.gz files <br> 
Example : `--dir "/home/user/InputDir/*.fastq.gz"`

**--results_path** : the path to the directory where generated files will be stored. <br>
Example: `--results_path "/home/user/OutputDir/"`

**-profile** : Specify how the pipeline will run (Docker, Singularity or Conda) and what reference data will be used (human or mouse). Examples : `-profile docker,human` or  `-profile conda,mouse`.

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

It is possible to add other species, but it will be necessary to add a new profile in nextflow.config for them.