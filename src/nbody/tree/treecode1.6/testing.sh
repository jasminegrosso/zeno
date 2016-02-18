#!/bin/bash
stars=1
timesteps=1
function treeforce {
  run=$(./treecode nbody=$stars timesteps=$timesteps )
  var=0
  for x in $run; do 
    #echo $x
    var=$((var+1))
    # startrun
    if [ $var -eq 25 ]
      then
        echo startrun: $x
        printf $x >> test.csv
        printf , >> test.csv
    fi
    # Initial treeforce
    if [ $var -eq 28 ]
      then
        echo Inital treeforce $x
        printf $x >> test.csv
        printf , >> test.csv
    fi
    # treeforce loop
    if [ $var -eq 33 ]
      then
        echo treeforce loop $x
        printf $x >> test.csv
        printf , >> test.csv
    fi
    # Total
    if [ $var -eq 35 ]
      then
        echo Total $x
        printf $x >> test.csv
    fi
  done
}
# Comment out if just want to add stuff to the end.
# printf stars,timesteps,startrun,initial_treeforce,treeforce,totalp > test.csv
# echo >> test.csv

while [ $stars -lt 10000 ]; do
  while [ $timesteps -lt 10000 ]; do
    echo ---$stars stars----
    echo ---$timesteps timesteps---
    printf $stars >> test.csv
    printf , >> test.csv
    printf $timesteps >> test.csv
    printf , >> test.csv
    treeforce
    let timesteps=timesteps*10
    echo >> test.csv
  done
  timesteps=1
  let stars=stars*10
done 