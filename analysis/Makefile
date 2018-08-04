R_OPT=--vanilla

SCRIPTS=scripts
FIG_SRC=scripts/figures
TAB_SRC=scripts/tables

MANUSCRIPT=manuscript
FIG_DIR=manuscript/figures
TAB_DIR=manuscript/tables

DATA=data
LOG=log

all: dir_tree \
	clean

# install/update the functions
install:
	Rscript -e "devtools::install('..')"

# make directory tree
dir_tree:
	test ! -d $(LOG) && mkdir $(LOG) || echo '$(LOG)/ is already there.'
	test ! -d $(DATA) && mkdir $(DATA) || echo '$(DATA)/ is already there.'
	test ! -d $(MANUSCRIPT) && mkdir $(MANUSCRIPT) || echo '$(MANUSCRIPT)/ is already there.'
	test ! -d $(TAB_DIR) && mkdir $(TAB_DIR) || echo '$(TAB_DIR)/ is already there.'
	test ! -d $(FIG_DIR) && mkdir $(FIG_DIR) || echo '$(FIG_DIR)/ is already there.'

# generate tables
#$(TAB_DIR)/%.tex: $(TAB_SRC)/%.R
#	R CMD BATCH $(R_OPT) $< $(LOG)/$(<F).Rout

# generate figures
#$(FIG_DIR)/%.png: $(FIG_SRC)/%.R
#	R CMD BATCH $(R_OPT) $< $(LOG)/$(<F).Rout

# generate workflow
#$(FIG_DIR)/workflow.png: $(FIG_SRC)/workflow.dot
#	dot -Tpng $(FIG_SRC)/workflow.dot -o $(FIG_DIR)/workflow.png

# clean up
.PHONY: clean
clean:
	rm -f *.pdf
	rm -f *.RData
