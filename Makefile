logfile := expressionVariance_$(shell date +%F)
all:
	pandoc main.md --lua-filter=pandoc/scholarly-metadata.lua --lua-filter=pandoc/author-info-blocks.lua  --filter pandoc-crossref --citeproc --pdf-engine=xelatex -o expressionVariance.pdf
	mv expressionVariance.pdf out/
doc:
	pandoc main.md --reference-doc=custom-reference.docx --lua-filter=pandoc/scholarly-metadata.lua --lua-filter=pandoc/author-info-blocks.lua  --filter pandoc-crossref --citeproc -o expressionVariance.docx 
	mv expressionVariance.docx out/
log:
	cp out/expressionVariance.pdf out/archive/$(logfile).pdf
	cp out/expressionVariance.docx out/archive/$(logfile).docx
cleanall:
	ls outfile* | grep -v main.md | xargs rm
