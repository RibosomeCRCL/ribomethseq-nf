
RESDIR = ribomethseq
QSIZE = 6
FQDIR = fastq
FQ_PATTERN = '.fastq.gz'
PROFILE = human
SCHEDULER = slurm

all: help

help:
	@echo "--------------------------------------------------"
	@echo "  test : test workflow using user's local PATH"
	@echo
	@echo "  test-conda : test workflow using conda profile"
	@echo
	@echo "  test-docker : test workflow using docker profile"
	@echo
	@echo "  test-singularity : test workflow using singularity profile"
	@echo
	@echo " clean : clean local nextflow stuff and results"
	@echo "--------------------------------------------------"

test: reset
	nextflow -bg run .. -profile  $(PROFILE) --scheduler $(SCHEDULER) \
	--qsize $(QSIZE) --outdir $(RESDIR) --logdir . \
	--fqdir $(FQDIR) --fastq_pattern $(FQ_PATTERN)

test-conda: reset
	$(MAKE) test PROFILE='human,conda'

test-docker: reset
	$(MAKE) test PROFILE='human,docker'

test-singularity: reset
	$(MAKE) test PROFILE='human,singularity'

clean:
	rm -rf .nextflow* work ribomethseq-nf.* $(RESDIR)

resdir:
	mkdir -p $(RESDIR)

reset: clean resdir
