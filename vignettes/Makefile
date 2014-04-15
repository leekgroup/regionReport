all:
	"$(R_HOME)/bin$(R_ARCH_BIN)/Rscript" -e "library('knitrBootstrap'); library('rmarkdown'); render('derfinderReport.Rmd', bootstrap_document(theme = 'default', highlight = 'Brown Paper', theme.chooser = TRUE, highlight.chooser = TRUE))"
clean:
	rm -f references.bib
	rm -fr derfinderReport_files
	rm -fr report/chr21
	rm -f report/*Rdata