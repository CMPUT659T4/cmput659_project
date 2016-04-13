FILENAME = T4_Report
RM = rm -f

$(FILENAME).pdf: $(FILENAME).tex $(FILENAME).aux
	pdflatex $<
	pdflatex $<

$(FILENAME).aux: $(FILENAME).tex $(FILENAME).bib
	pdflatex $<
	pdflatex $<
	bibtex $(FILENAME).aux

.PHONY: clean
clean: 
	$(RM) $(FILENAME).pdf
	$(RM) $(FILENAME).aux
	$(RM) $(FILENAME).dvi
	$(RM) $(FILENAME).bbl
	$(RM) $(FILENAME).log
	$(RM) $(FILENAME).out
	$(RM) $(FILENAME).synctex*

