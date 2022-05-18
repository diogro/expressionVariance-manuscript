# Assume pandoc-scholar is in the parent directory.
# Should usually be overwritten or configured via an environment variable.
PANDOC_SCHOLAR_PATH ?= $(PWD)/pandoc-scholar
ARTICLE_FILE=main.md
BIBLIOGRAPHY_FILE=references.bib
DEFAULT_EXTENSIONS = latex
include $(PANDOC_SCHOLAR_PATH)/Makefile

pdf:
	xelatex outfile.latex
	biber outfile
	xelatex outfile.latex
	xelatex outfile.latex
	mv outfile.pdf out/
	ls outfile* | grep -v main.md | xargs rm
cleanall:
	ls outfile* | grep -v main.md | xargs rm
