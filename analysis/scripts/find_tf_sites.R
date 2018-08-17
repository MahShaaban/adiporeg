# load required libraries
library(tidyverse)
library(reshape2)
library(IRanges)
library(GenomicRanges)
library(rtracklayer)
library(Gviz)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(BSgenome.Hsapiens.UCSC.hg38)
library(trackViewer)

get_coords <- function(gene_id, txdb) {
  gene <- genes(txdb, filter = list(gene_id = gene_id))
  chr <- as.character(seqnames(gene))
  seqlevels(gene) <- chr

  transcript <- transcripts(txdb, filter = list(gene_id = gene_id))
  seqlevels(transcript) <- chr

  exon <- exons(txdb, filter = list(gene_id = gene_id))
  seqlevels(exon) <- chr

  promoter <- promoters(txdb, filter = list(gene_id = gene_id))
  seqlevels(promoter) <- chr

  list(gene = gene,
       transcript = transcript,
       exon = exon,
       promoter = promoter)
}

get_coverage <- function(fls, gene.gr) {
  bw.dfs <- map(fls, function(x) {
    import.bw(con = BigWigFile(x),
              selection = BigWigSelection(gene.gr)) %>%
      as.data.frame()
  })
  names(bw.dfs) <- 1:length(bw.dfs)

  bw.gr <- bw.dfs %>%
    bind_rows(.id = 'sample') %>%
    filter(score < 1) %>%
    makeGRangesFromDataFrame(keep.extra.columns = TRUE)
  bw.gr
}

plot_tracks <- function(gene, coverage) {
  chr = as.character(seqnames(gene$gene))
  start = start(gene$gene)
  end = end(gene$gene)
  genome = as.character(genome(gene$gene))

  itrack <- IdeogramTrack(chromosome = chr,
                          genome = genome)

  xtrack <- GenomeAxisTrack()

  gtrack <- GeneRegionTrack(txdb,
                            chromosome = chr,
                            start = start,
                            end = end,
                            transcriptAnnotation = 'symbol')

  atrack <- AnnotationTrack(gene$gene, name = 'Gene')
  ttrack <- AnnotationTrack(gene$transcript, name = 'Transcript')
  etrack <- AnnotationTrack(gene$exon, name = 'Exon')
  ptrack <- AnnotationTrack(gene$promoter, name = 'Promotor')

  dtrack <- DataTrack(coverage, data = coverage$score, type = 'h')

  plotTracks(list(
    itrack,
    xtrack,
    gtrack,
    atrack,
    ttrack,
    etrack,
    ptrack,
    dtrack
  ))
}

get_peaks <- function(coverage, q, return = 'gr') {
  quant <- quantile(coverage$score, q)
  sub.cov <- coverage[coverage$score > quant]
  p <- reduce(sub.cov)
  t <- paste0(as.character(seqnames(p)), ':',
              as.integer(start(p)), '-',
              as.integer(end(p)))

  switch (return,
    'gr' = p,
    'text' = t
  )
}

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
hs <- BSgenome.Hsapiens.UCSC.hg38

ecrr6_files <- list.files('data/bigwig/', pattern = 'ercc6_*', full.names = TRUE)
vezf1.files <- list.files('data/bigwig/', pattern = 'vezf1_*', full.names = TRUE)

pebp1 <- 5037
pebp1.coords <- get_coords(pebp1, txdb)
pebp1.coverage <- get_coverage(ecrr6_files, pebp1.coords$gene)
plot_tracks(pebp1.coords, pebp1.coverage)


p.pebp1 <- get_peaks(pebp1.coverage, .95)
get_peaks(pebp1.coverage, .95, return = 'text')
fa <- getSeq(hs, p.pebp1)
fa
writeXStringSet(fa, '~/Desktop/ercc6.pebp1.fa')

##

tbc1d5 <- 9779
tbc1d5.coords <- get_coords(tbc1d5, txdb)
tbc1d5.coverage <- get_coverage(ecrr6_files, tbc1d5.coords$gene)
plot_tracks(tbc1d5.coords, tbc1d5.coverage)

p.tbc1d5 <- get_peaks(tbc1d5.coverage, .99)
p.tbc1d5  <- p.tbc1d5[unlist(lapply(p.tbc1d5, width)) > 50]

resized.tbc1d5 <- resize(p.tbc1d5, 150, fix = 'center')
fa <- getSeq(hs, resized.tbc1d5)

writeXStringSet(fa, '~/Desktop/ercc6.tbc1d5.fa')


FLYWCH2 = 114984
FLYWCH2.coords <- get_coords(FLYWCH2, txdb)
FLYWCH2.coverage <- get_coverage(ecrr6_files, FLYWCH2.coords$gene)
plot_tracks(FLYWCH2.coords, FLYWCH2.coverage)

p.FLYWCH2 <- get_peaks(FLYWCH2.coverage, .95)
fa <- getSeq(hs, p.FLYWCH2)
fa
writeXStringSet(fa, '~/Desktop/ercc6.FLYWCH2.fa')


pebp1.coverage2 <- get_coverage(vezf1.files, pebp1.coords$gene)
plot_tracks(pebp1.coords, pebp1.coverage2)

p.pebp1 <- get_peaks(pebp1.coverage, .95)
fa <- getSeq(hs, p.pebp1)
fa
writeXStringSet(fa, '~/Desktop/vezf1.pebp1.fa')

pik3c3 <- 5289
pik3c3.coords <- get_coords(pik3c3, txdb)
pik3c3.coverage <- get_coverage(vezf1.files, pik3c3.coords$gene)
plot_tracks(pik3c3.coords, pik3c3.coverage)

p.pik3c3 <- get_peaks(pik3c3.coverage, .99)
fa <- getSeq(hs, p.pik3c3)
fa
writeXStringSet(fa, '~/Desktop/vezf1.pik3c3.fa')


wdr45 <- 11152
wdr45.coords <- get_coords(wdr45, txdb)
wdr45.coverage <- get_coverage(vezf1.files, wdr45.coords$gene)
plot_tracks(wdr45.coords, wdr45.coverage)

p.wdr45<- get_peaks(wdr45.coverage, .99)


resized.wdr45 <- resize(p.wdr45, 150, fix = 'center')

fa <- getSeq(hs, resized.wdr45)
fa
writeXStringSet(fa, '~/Desktop/vezf1.wdr45.fa')

tmem230 <- 29058
tmem230.coords <- get_coords(tmem230, txdb)
tmem230.coverage <- get_coverage(vezf1.files, tmem230.coords$gene)
plot_tracks(tmem230.coords, tmem230.coverage)

p.tmem230<- get_peaks(tmem230.coverage, .99)
fa <- getSeq(hs, p.tmem230)
fa
writeXStringSet(fa, '~/Desktop/vezf1.tmem230.fa')
