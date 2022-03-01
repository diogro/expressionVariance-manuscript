all:
	pandoc main.md --bibliography references.bib --citeproc -o out/main.pdf --pdf-engine=xelatex