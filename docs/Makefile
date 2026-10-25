.PHONY: all pdf clean cleanall

MAIN := main.tex
BUILD := build

all: pdf

pdf:
	latexmk -pdf -interaction=nonstopmode -halt-on-error -file-line-error -outdir=$(BUILD) $(MAIN)

clean:
	latexmk -c -outdir=$(BUILD) $(MAIN)

cleanall:
	latexmk -C -outdir=$(BUILD) $(MAIN)
