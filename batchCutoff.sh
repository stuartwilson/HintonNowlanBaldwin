#!/bin/bash
for ((i=0;i<1000;i+=10));
do
    echo $i
    for ((j=0;j<1000;j++));
    do
        ./build/hintonBatch config.json logs/cutoff$i/sim$j $i &
    done

done
