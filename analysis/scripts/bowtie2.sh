#!/bin/bash
# define variables
INDEX='data/mm10/mm10'
FASTQ=$(ls data/fastq | cut -d '.' -f1)

# make directory of the alignment output
test ! -d data/sam && mkdir data/sam || echo 'Already exists'

# run bowtie2 on all files
for i in $FASTQ; do
  if [ ! -f "data/sam/$i.sam" ]; then
    bowtie2 --no-unal -P 8 -x $INDEX -U data/fastq/$i.fastq.gz -S data/sam/$i.sam
    echo "data/sam/$i.sam was created."
  fi
done
