.PHONY: install test example docker

install:
	python -m pip install -e '.[test]'

test:
	pytest

example:
	nmr2boltz convert examples/example.nef -o artifacts/example --hypotheses 4

docker:
	docker build -t nmr2boltz:0.1.0 .
