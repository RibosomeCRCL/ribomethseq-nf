
# Tests directory

We generated a tiny dataset (`fastq` folder) which allows to quickly test if
your workflow environment is properly working.

Running this test should take more than 1 or 2 minutes.

A Makefile is also provided to help you test the workflow :

```
$ make help
--------------------------------------------------
test : test workflow using user's local PATH

test-conda : test workflow using conda profile

test-docker : test workflow using docker profile

test-singularity : test workflow using singularity profile

clean : clean local nextflow stuff and results
--------------------------------------------------
```

## Local installation

If all the required dependencies are installed in your local environment and all
the software are available on your PATH, you can simply run `make test`.

## Conda installation

Once you have set up local conda envs (preferred over building the environments
each time the workflow is launched) and modified the `conda` profile accordingly
in the config, you can test the workflow with `make test-conda`.

## Docker installation

Similarly, test your docker configuration with `make test-docker`. If necessary,
be sure to provide necessary bind mounts _via_ the `process.containerOptions`
directive.

## Singularity installation

Test your singularity configuration with `make test-singularity`. If necessary,
be sure to provide necessary bind mounts _via_ the `process.containerOptions`
directive.
