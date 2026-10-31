.PHONY: help test fixture run validate clean

SERVICE_NAME ?= sc-markers
PYTHON ?= python3
INPUT ?= ./data/fixture/linked

.DEFAULT: help

help:
	@echo "Make Help for $(SERVICE_NAME)"
	@echo ""
	@echo "make test     - pytest: the synthetic fixture through the processor, both input layouts (no Docker)"
	@echo "make fixture  - write the synthetic fixture to data/fixture/{linked,flat}/"
	@echo "make run      - build the image and run it in Docker over INPUT=$(INPUT), writing to data/output/"
	@echo "make clean    - remove data/output/"
	@echo ""
	@echo "PYTHON=$(PYTHON) (set PYTHON=.venv/bin/python after pip install -r requirements-dev.txt)"

test:
	$(PYTHON) -m pytest -q tests

fixture:
	rm -rf data/fixture
	$(PYTHON) -m tests.make_fixture data/fixture --layout both

run: clean
	@if [ "$(INPUT)" = "./data/fixture/linked" ]; then $(MAKE) fixture PYTHON=$(PYTHON); fi
	mkdir -p data/output
	INPUT=$(INPUT) docker compose -f docker-compose.yml down --remove-orphans
	INPUT=$(INPUT) docker compose -f docker-compose.yml build
	INPUT=$(INPUT) docker compose -f docker-compose.yml up --exit-code-from processor
	INPUT=$(INPUT) docker compose -f docker-compose.yml down --remove-orphans
	@echo "--- outputs:"
	@find data/output -type f | sort

clean:
	rm -rf data/output

validate:
	$(PYTHON) scripts/validate.py
