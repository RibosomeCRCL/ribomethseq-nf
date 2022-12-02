
RESDIR = ribomethseq

QSIZE = 6
FQDIR = fastq
FQ_PATTERN = '_R1_001.fastq.gz'
PROFILE = human

all: help

help:
	@echo "----------------------------------------"
	@echo "  test : launch test workflows"
	@echo
	@echo " clean : clean local nextflow stuff and results"
	@echo "----------------------------------------"

test: reset
	nextflow -bg run .. -profile  $(PROFILE) \
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
