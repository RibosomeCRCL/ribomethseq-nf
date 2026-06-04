# Container image

The pipeline ships a single container image holding **all** tools (processing +
R reporting). The recommended way to get it is to pull the prebuilt image from
GHCR (see the main README) — you normally do **not** need to build it yourself.

### Dockerfile.prod

The production image. It is built on the official **`rocker/r-ver`** base
(versioned R on Ubuntu LTS) using a multi-stage build:

- the `builder` stage holds the compilers / `-dev` headers and builds the only
  tool that is compiled (samtools, htslib bundled); every other tool is a
  prebuilt binary;
- the final stage copies the binaries in and keeps only runtime shared
  libraries, for a slimmer image.

R packages are installed from the **Posit Package Manager** binary repo pinned to
a dated snapshot, so the package set is reproducible. All tool/R versions are set
via `ARG`s at the top of the file. See the header comments in `Dockerfile.prod`
for the rationale and the noble-specific notes.

To build it locally (optional):

```sh
cd docker
docker build -t ribomethseq-nf:1.1 -f Dockerfile.prod .
```

### Dockerfile

Legacy variant (everything in a single stage). Kept for reference but
the maintained image is `Dockerfile.prod`.

### trimmomatic

A small wrapper script so Trimmomatic can be called as `trimmomatic [args]`
(matching the way it is invoked in the pipeline). It is copied into the image's
`PATH` and points to the bundled `trimmomatic-<version>.jar`.
