
.PHONY: all pdf clean cleanall watch

MAIN := report/main.tex
BUILD := build

all: pdf

pdf:
	latexmk -pdf -interaction=nonstopmode \
		-halt-on-error -file-line-error \
		-outdir=$(BUILD) $(MAIN)

watch:
	latexmk -pdf -pvc -interaction=nonstopmode \
		-halt-on-error -outdir=$(BUILD) $(MAIN)

clean:
	latexmk -c -outdir=$(BUILD) $(MAIN)

cleanall:
	latexmk -C -outdir=$(BUILD) $(MAIN)
