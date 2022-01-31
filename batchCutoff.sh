#!/bin/bash
for ((i=0;i<1000;i+=10));
do
    echo $i
    for ((j=0;j<1000;j++));
    do
	seed=$((i*1000+j))
        ./build/hintonBatch configBatch.json logs/cutoff$i/sim$j $i $seed &
    done

done
