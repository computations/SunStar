PDF=pdflatex -halt-on-error -quiet

FIGS := $(wildcard *.tex)
PDF_FIGS := $(FIGS:.tex=.pdf)

all: $(PDF_FIGS)
	rm -rf *.aux *.log

%.pdf: %.tex
	$(PDF) $^ > /dev/null 2>&1

clean:
	rm -rf $(PDF_FIGS)
