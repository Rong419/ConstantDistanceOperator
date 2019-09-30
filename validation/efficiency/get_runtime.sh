#!/bin/sh


 for param in {1..2}
 do 
    for sim in {1..100}
    do
    sed -n 37p ~/Desktop/efficiency/anolis/times/error_anolis_${param}_${sim}.txt >>./anolis/time_anolis.txt
    done
  done

 for param in {1..2}
 do 
    for sim in {1..100}
    do
    sed -n 36p ~/Desktop/efficiency/RSV2/times/error_RSV2_${param}_${sim}.txt >>./RSV2/time_RSV2.txt
    done
  done