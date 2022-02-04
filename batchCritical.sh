#!/bin/bash
seed=0
for ((k=0;k<20;k++));
do
    echo $k
    for ((i=0;i<1000;i+=10));
    do
        for ((j=i+10;j<1000;j+=10));
        do
            ./build/hintonBatchCritical config.json logs/Lower${i}_Upper${j}/sim${k} $i $j $seed &
            ((seed++))
        done
    done

done
