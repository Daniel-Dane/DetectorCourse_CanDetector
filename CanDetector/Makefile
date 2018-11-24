# Makefile for creating document
#------------------------------------------------------------------------------
# By default makes mydocument.pdf using target run_pdflatex.
# Replace mydocument with your main filename or add another target set.
# Adjust TEXLIVE if it is not correct, or pass it to "make new".
# Replace BIBTEX = bibtex with BIBTEX = biber if you use biber instead of bibtex.
# Adjust FIGSDIR for your figures directory tree.
# Adjust the %.pdf dependencies according to your directory structure.
# Use "make clean" to cleanup.
# Use "make cleanpdf" to delete $(BASENAME).pdf.
# "make cleanall" also deletes the PDF file $(BASENAME).pdf.
# Use "make cleanepstopdf" to rmeove PDF files created automatically from EPS files.
#   Note that FIGSDIR has to be set properly for this to work.

# If you have to run latex rather than pdflatex adjust the dependencies of %.dvi target
#   and use the command "make run_latex" to compile.
# Specify dvipdf or dvips as the run_latex dependency,
#   depending on which you want to use.

#-------------------------------------------------------------------------------
# Check which TeX Live installation you have with the command pdflatex --version
TEXLIVE  = 2016
LATEX    = latex
PDFLATEX = pdflatex
#BIBTEX   = bibtex
BIBTEX = biber
DVIPS    = dvips
DVIPDF   = dvipdf

#-------------------------------------------------------------------------------
# The main document filename
BASENAME = master
#-------------------------------------------------------------------------------
# Adjust this according to your top-level figures directory
# This directory tree is used by the "make cleanepstopdf" command
FIGSDIR  = graphics
#-------------------------------------------------------------------------------

# EPSTOPDFFILES = `find . -name \*eps-converted-to.pdf`
rwildcard=$(foreach d,$(wildcard $1*),$(call rwildcard,$d/,$2) $(filter $(subst *,%,$2),$d))
EPSTOPDFFILES = $(call rwildcard, $(FIGSDIR), *eps-converted-to.pdf)

# Default target - make mydocument.pdf with pdflatex
default: all

.PHONY: all clean cleanpdf help

all: pre-build main-build

pre-build:
	@echo PRE

post-build:
	@echo POST

# Standard pdflatex target
main-build: $(BASENAME).pdf
	@echo "Made $<"


#-------------------------------------------------------------------------------
# Specify the tex and bib file dependencies for running pdflatex
# If your bib files are not in the main directory adjust this target accordingly
#%.pdf:	%.tex *.tex bibtex/bib/*.bib
%.pdf:	%.tex
	$(PDFLATEX) $<
	-$(BIBTEX)  $(basename $<)
	$(PDFLATEX) $<
	$(PDFLATEX) $<
#-------------------------------------------------------------------------------

help:
	@echo ""
	@echo "make               to produce the output PDF"
	@echo "make clean         to clean auxiliary files (not output PDF)"
	@echo "make cleanpdf      to clean output PDF files"
	@echo "make cleanps       to clean output PS files"
	@echo "make cleanall      to clean all files"
	@echo "make cleanepstopdf to clean PDF files automatically made from EPS"
	@echo ""

clean:
	-rm *.dvi *.toc *.aux *.log *.out \
		*.bbl *.blg *.brf *.bcf *-blx.bib *.run.xml \
		*.cb *.ind *.idx *.ilg *.inx \
		*.synctex.gz *~ *.fls *.fdb_latexmk .*.lb spellTmp 

cleanpdf:
	-rm $(BASENAME).pdf 
	-rm $(BASENAME)-draft-cover.pdf $(BASENAME)-preprint-cover.pdf
	-rm $(BASENAME)-auxmat.pdf

cleanps:
	-rm $(BASENAME).ps 
	-rm $(BASENAME)-draft-cover.ps $(BASENAME)-preprint-cover.ps
	-rm $(BASENAME)-auxmat.ps

cleanall: clean cleanpdf cleanps

# Clean the PDF files created automatically from EPS files
cleanepstopdf: $(EPSTOPDFFILES)
	@echo "Removing PDF files made automatically from EPS files"
	-rm $^
