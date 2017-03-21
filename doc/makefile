DOCNAME=writeup.pdf
MDFILE =writeup.md
BIBFILE=bib.yaml

all: $(DOCNAME)

$(DOCNAME): $(MDFILE) $(BIBFILE)
	pandoc --filter=pandoc-citeproc -f markdown-superscript $(MDFILE) -o $@

see: $(DOCNAME)
	evince $(DOCNAME)

clean: 
	rm -rf $(DOCNAME)
