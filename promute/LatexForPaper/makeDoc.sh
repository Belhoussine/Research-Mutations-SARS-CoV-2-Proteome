#!/bin/bash

# Remove old files
rm -f *.log *.aux *.pdf *.out *.gz *.blg *.bbl *.dvi *.ps

# create new doc
pdflatex proMute2
bibtex proMute2
pdflatex proMute2
pdflatex proMute2

# Remove just created intermediate files
rm -f *.log *.aux *.out *.gz *.blg *.bbl *.dvi *.ps



