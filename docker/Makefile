version ?= v1.0

all: build

build:
	cp ../conda.yml .
	cp -R ../Rscripts .
	cp -R ../ref_human .
	docker build --output=type=docker --progress=plain .
