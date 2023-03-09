main := expressionVariance
logfile := $(main)_$(shell date +%F)
all:
	pandoc main.md --template=pandoc/default.latex --lua-filter=pandoc/scholarly-metadata.lua --lua-filter=pandoc/author-info-blocks.lua  --filter pandoc-crossref --citeproc --pdf-engine=xelatex -o $(main).pdf
	mv $(main).pdf out/
	pandoc main.md --reference-doc=pandoc/custom-reference.docx --template=pandoc/default.latex --lua-filter=pandoc/scholarly-metadata.lua --lua-filter=pandoc/author-info-blocks.lua  --filter pandoc-crossref --citeproc -o $(main).docx
	mv $(main).docx out/
	cp out/$(main).pdf out/archive/$(logfile).pdf
	cp out/$(main).docx out/archive/$(logfile).docx
pdf:
	pandoc main.md --lua-filter=pandoc/scholarly-metadata.lua --lua-filter=pandoc/author-info-blocks.lua  --filter pandoc-crossref --citeproc --pdf-engine=xelatex -o $(main).pdf
	mv $(main).pdf out/
doc:
	pandoc main.md --reference-doc=pandoc/custom-reference.docx --lua-filter=pandoc/scholarly-metadata.lua --lua-filter=pandoc/author-info-blocks.lua  --filter pandoc-crossref --citeproc -o $(main).docx
	mv $(main).docx out/
log:
	cp out/$(main).pdf out/archive/$(logfile).pdf
	cp out/$(main).docx out/archive/$(logfile).docx
cleanall:
	ls outfile* | grep -v main.md | xargs rm
