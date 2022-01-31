#!/bin/bash
for ((i=0;i<1000;i+=10));
do
    for ((j=0;j<1000;j+=10));
    do
        echo lower=$i, upper=$j
        for ((k=0;j<100;j++));
        do
	        seed=$(((i*1000+j*1000)*100+k))
            ./build/hintonCritical config.json logs/Lower$i_Upper$j/sim$k $i $j $seed &
        done
    done

done
