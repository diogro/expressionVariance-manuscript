# Assume pandoc-scholar is in the parent directory.
# Should usually be overwritten or configured via an environment variable.
PANDOC_SCHOLAR_PATH   ?= $(shell dirname $(PWD))
ARTICLE_FILE=main.md
BIBLIOGRAPHY_FILE=references.bib
include $(PANDOC_SCHOLAR_PATH)/Makefile

all:
	xelatex outfile.latex
	biber outfile
	xelatex outfile.latex
	xelatex outfile.latex
	mv outfile.pdf out/
	ls outfile* | grep -v main.md | xargs rm
clean:
	ls outfile* | grep -v main.md | xargs rm