
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
	@echo " clean : clean local nextflow stuff and results"
	@echo "--------------------------------------------------"

test: reset
	nextflow -bg run .. -profile  $(PROFILE) --scheduler $(SCHEDULER) \
	--qsize $(QSIZE) --outdir $(RESDIR) --logdir . \
	--fqdir $(FQDIR) --fastq_pattern $(FQ_PATTERN)

test-docker: reset
	$(MAKE) test PROFILE='human,docker'

test-conda: reset
	$(MAKE) test PROFILE='human,conda'

clean:
	rm -rf .nextflow* work ribomethseq-nf.* $(RESDIR)

resdir:
	mkdir -p $(RESDIR)

reset: clean resdir
