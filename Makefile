.PHONY: install build queries check benchmark test clean all

PYTHON ?= python3

all: build queries benchmark

install:
	$(PYTHON) -m pip install -e ".[dev]"

## Create the schema and load the Ensembl exports into a local SQLite file
build:
	$(PYTHON) -m genomedb.cli build

## Run Q1-Q7 and write results/
queries:
	$(PYTHON) -m genomedb.cli query

## Integrity and consistency checks
check:
	$(PYTHON) -m genomedb.cli check

## Time every query with and without the secondary indexes
benchmark:
	$(PYTHON) -m genomedb.cli benchmark

test:
	$(PYTHON) -m pytest -q

clean:
	rm -f data/genome.sqlite
	rm -rf results/*.tsv results/*.json results/*.md
	find . -name __pycache__ -type d -exec rm -rf {} +
