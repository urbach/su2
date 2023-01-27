# Makefile

# file paths and variables
main = ./main.Rmd

define generate_html
	Rscript -e 'library(rmarkdown); rmarkdown::render("$(1)", "html_document")'
	Rscript -e 'library(rmarkdown); rmarkdown::render("$(1)", "pdf_document")'
endef

all: main
	make clean

main: $(main)
	$(call generate_html,$(main))


clean:
	@rm *~ *.gz *.fls *.fdb_latexmk *.log *.aux *.xml \
	*.bbl *.blg *.out diagram.* *.bcf *.mp  *.dvi \
	2> /dev/null || true
	

