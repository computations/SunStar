DOCNAME=writeup.pdf
TITLEPAGE=titlepage.pdf
TITLEPAGE_S=titlepage.tex
MDFILE =writeup.md
BIBFILE=bib.yaml
FIGDIR=figs

all: $(DOCNAME)

$(DOCNAME): $(MDFILE) $(BIBFILE) $(TITLEPAGE)
	make -C $(FIGDIR)
	pandoc --filter=pandoc-citeproc -f markdown-superscript $(MDFILE) -o $@

$(TITLEPAGE): $(TITLEPAGE_S)
	pdflatex --halt-on-error $^ > /dev/null 2>&1

see: $(DOCNAME)
	evince $(DOCNAME)

clean: 
	make -C $(FIGDIR) clean
	rm -rf $(DOCNAME) $(TITLEPAGE)
