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
	get_fastq \
	align_reads \
	macs2 \
	annotate_peaks \
	clean_data \
	$(FIG_DIR)/available_data.png \
	$(FIG_DIR)/tf_tf_arcs.png \
	$(TAB_DIR)/chip_source.tex \
	$(TAB_DIR)/tf_fold_enrichment.tex \
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

# bigwig to bedgraph
#bwtobdg: $(SCRIPTS)/bw_to_bdg.sh $(DATA)/bigwig/*
#	sh $(SCRIPTS)/bw_to_bdg.sh

# run macs2
#macs2: $(SCRIPTS)/macs2.sh $(DATA)/bedgraph/*
#	sh $(SCRIPTS)/macs2.sh

# get fastq files
get_fastq: $(SCRIPTS)/download_fastq.sh $(DATA)/tf_full.urls
	sh $(SCRIPTS)/download_fastq.sh

# align reads
align_reads: $(SCRIPTS)/bowtie2.sh $(DATA)/fastq/* $(DATA)/mm10/*
	sh $(SCRIPTS)/bowtie2.sh

# sam to bam
sam_to_bam: $(SCRIPTS)/sam_to_bam.sh $(DATA)/sam/*
	sh $(SCRIPTS)/sam_to_bam.sh

# run macs2
macs2: $(SCRIPTS)/macs2.sh $(DATA)/sam/*
	sh $(SCRIPTS)/macs2.sh

# annotate peaks
annotate_peaks: $(SCRIPTS)/annotate_peaks.R $(DATA)/bed/*
	R CMD BATCH $(R_OPT) $< $(LOG)/$(<F).Rout

clean_data: $(SCRIPTS)/clean_data.R $(DATA)/annotated_peaks/*
	R CMD BATCH $(R_OPT) $< $(LOG)/$(<F).Rout

# generate tables
$(TAB_DIR)/%.tex: $(TAB_SRC)/%.R $(DATA)/chip.csv $(DATA)/tf_tf_targets.csv
	R CMD BATCH $(R_OPT) $< $(LOG)/$(<F).Rout

# generate figures
$(FIG_DIR)/%.png: $(FIG_SRC)/%.R $(DATA)/chip.csv
	R CMD BATCH $(R_OPT) $< $(LOG)/$(<F).Rout

# generate workflow
#$(FIG_DIR)/workflow.png: $(FIG_SRC)/workflow.dot
#	dot -Tpng $(FIG_SRC)/workflow.dot -o $(FIG_DIR)/workflow.png

# clean up
.PHONY: clean
clean:
	rm -f *.pdf
	rm -f *.RData
