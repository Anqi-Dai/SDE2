#!/bin/bash

#$ -P flynngrp
#$ -cwd

python rMATS_pkg/rMATS.4.0.2/rMATS-turbo-Linux-UCS4/rmats.py --b1 data/RF89.bam.txt  --b2 data/RFControl.bam.txt \
    --gtf /usr4/bs831/adai/bubhub-home/SDE2/reference/v27/gencode.v27.annotation.gtf --od output/RF89_control/  -t paired \
    --readLength  76  --cstat 0.1 --libType fr-secondstrand 
