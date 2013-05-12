# use pdflatex for converting the tex file to a pdf file
TEX     = pdflatex 
# use image magick's convert in order to get a png file
CONVERT = convert 

# input are all tex files of this directory
SRC	:= $(shell egrep -l '^[^%]*\\begin\{document\}' *.tex)

# target are all png files generated from the tex files
TRG	= $(SRC:%.tex=%.png)
all: $(TRG)

# convert a tex to a pdf file
%.pdf: %.tex
		$(TEX) $<

# convert a pdf to a png file
%.png: %.pdf
	$(CONVERT) $< $@

