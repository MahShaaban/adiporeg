#!/bin/bash
# define variables
SAM=$(ls data/sam | cut -d '.' -f1)

# make directory of bam output
test ! -d data/bam && mkdir data/bam || echo 'Already exists'

# run bowtie2 on all files
for i in $SAM; do
  if [ ! -f "data/bam/$i.bam" ]; then
    samtools view -Sb data/sam/$i.sam > data/bam/$i.bam
    echo "data/bam/$i.bam was created." >> log/sam_to_bam.out
  fi
done
