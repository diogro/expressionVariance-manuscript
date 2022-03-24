all:
	pandoc main.md -s --filter pandoc-crossref --bibliography references.bib --biblatex -o main.tex --pdf-engine=xelatex
	xelatex main
	biber main
	xelatex main
	xelatex main
	mv main.pdf out/
	ls main* | grep -v main.md | xargs rm
clean:
	ls main* | grep -v main.md | xargs rm