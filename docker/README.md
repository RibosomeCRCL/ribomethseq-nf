# Execution environments

### Dockerfile

This is the file used to test/debug the docker image. It is built upon debian 11
(bullseye) and it contains **everything** (libraries required for compilation,
all sources, ...) and it is then quite huge (~2.7Gb).

### Dockerfile.prod

All tools are built like it is done in Dockerfile but the final image (still
based on debian 11) is reduced to only (hopefully) the runtime requirements.
This image is much "lighter" than the dev one and is about 1.3Gb in size. The
compressed docker archive is ~490Mb.

### Dockerfile.conda - DEPRECATED

Original docker image built with conda.

### trimmomatic

This is an helper wrapper script that allows to call trimmomatic in the form
`trimmomatic [args]` as it is done when built with conda. This way calling
trimmomatic in conda, docker or singularity environments is done exactly the
same way.

You may need to do something similar if you are using the trimmomatic jar file
directly (copying it into your `PATH` + adapting the path to trimmotatic jar
	file should do the trick).

### rms-processing.yml and rms-report.yml

Initially we tried to create a single conda environment but due to conda solving
time issue we chose to split the environment dedicated to sequencing data processing
and the one for creating a report (R packages).

The recommended way to use them is to create both these conda envs in your system
and adapt the conda profile in `nextflow.config`.
