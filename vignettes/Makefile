# Makefile to save the result from the bumphunter example
# Modified from https://hedgehog.fhcrc.org/bioconductor/trunk/madman/Rpacks/ggbio/vignettes/Makefile

HTMLS= bumphunterExample.html regionReport.html

all: $(HTMLS)

clean:
	cp fakeVignettes/*.Rmd .

%.html: %.Rmd
	"$(R_HOME)/bin$(R_ARCH_BIN)/Rscript" -e "library(rmarkdown); render('$*.Rmd')"
