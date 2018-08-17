#!/bin/bash

test ! -d data/fastq && mkdir data/fastq || echo 'data/fastq/ is already there.'

cat data/tf_full.urls | xargs -P76 -n 1 wget -q -nc -P data/fastq/
