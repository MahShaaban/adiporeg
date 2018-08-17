#!/bin/bash

BW=$(ls data/bigwig | cut -d '.' -f1)

# make output director
test ! -d data/bedgraph && mkdir data/bedgraph || echo 'Already exists.'

# run bigwigtobedgraph
for i in $BW; do
  if [ ! -e "data/bedgraph/$i.bdg" ]; then
    bigWigToBedGraph data/bigwig/$i.bw data/bedgraph/$i.bgd
  fi
done
