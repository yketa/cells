SHELL:=/bin/bash
CELLS:=src/cells

all: container.sif

container.sif: container.def requirements.txt pyproject.toml $(wildcard *.md) _custom_build.py $(CELLS)/Makefile $(wildcard $(CELLS)/*.*pp) $(wildcard $(CELLS)/*.py) $(wildcard $(CELLS)/*.sh) $(wildcard $(CELLS)/docs/*.pdf)
	sudo singularity build $@ $<

