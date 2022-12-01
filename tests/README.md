
# Tests directory

In this directory, store a **tiny** dataset which allows to run the pipeline in
a reasonable time. Total size should not be more than a few megabytes.

If it is not possible to have (or to generate) a very small sized test dataset,
it is better to keep it external from the git repository and give instructions
in this _README.md_ on how access these data. Heavy datasets are a pain when
manipulating git repositories.

Beyond a successful pipeline completion, if possible, explain how to check that
results given by the execution are indeed what we expect from this test.
