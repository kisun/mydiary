#!/bin/bash

for i in $(find ./ -type f -name "*.fastq.gz" | while read F; do basename $F | cut -d _ -f 1,2,3 ; done | sort | uniq)

    do echo "Merging R1"

zcat "$i"_S*_L00*_R1_001.fastq.gz > "$i"_R1.fq.gz

       echo "Merging R2"

zcat "$i"_S*_L00*_R2_001.fastq.gz > "$i"_R2.fq.gz

done;


