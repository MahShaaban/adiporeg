#!/bin/bash
# define variables
SAM=$(ls data/sam | cut -d '.' -f1)

# make output director
test ! -d data/bed && mkdir data/bed || echo 'Already exists.'

# run bowtie2 on all files
for i in $SAM; do
  PEAKS=$(printf "data/bed/%s_peaks.xls" "$i")
  if [ ! -f $PEAKS ]; then
    macs2 callpeak -f SAM --outdir data/bed -t data/sam/$i.sam -n $i
  fi
done
