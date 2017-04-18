DOCNAME=writeup.pdf
MDFILE =writeup.md
BIBFILE=bib.yaml
FIGDIR=figs

all: $(DOCNAME)

$(DOCNAME): $(MDFILE) $(BIBFILE) 
	make -C $(FIGDIR)
	pandoc --filter=pandoc-citeproc -f markdown-superscript $(MDFILE) -o $@

see: $(DOCNAME)
	evince $(DOCNAME)

clean: 
	make -C $(FIGDIR) clean
	rm -rf $(DOCNAME)
