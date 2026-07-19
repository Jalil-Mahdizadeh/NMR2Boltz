.PHONY: install test example benchmark-corpus distance-check docker

install:
	python -m pip install -e '.[test]'

test:
	pytest

example:
	nmr2boltz convert examples/example.nef -o artifacts/example --hypotheses 4

benchmark-corpus:
	python validation/benchmark_corpus.py benchmark/input --output-directory benchmark/output

distance-check:
	python validation/distance_check.py benchmark/input --conversion-output benchmark/output --output-directory benchmark/distance_check

docker:
	docker build -t nmr2boltz:0.1.0 .
