#!/bin/bash
stars=1
timesteps=1
function run_cs257 {
  run=$(./cs257 $stars $timesteps 0)
  var=0
  for x in $run; do 
    #echo $x
    var=$((var+1))
    # Loop 0
    echo $x
    # if [ $var -eq 16 ]
    #   then
    #     echo Loop 0 $x
    #     printf $x >> test.csv
    #     printf , >> test.csv
    # fi
    # # Loop 1
    # if [ $var -eq 21 ]
    #   then
    #     echo Loop 1 $x
    #     printf $x >> test.csv
    #     printf , >> test.csv
    # fi
    # # Loop 2
    # if [ $var -eq 26 ]
    #   then
    #     echo Loop 2 $x
    #     printf $x >> test.csv
    #     printf , >> test.csv
    # fi
    # # Loop 3
    # if [ $var -eq 31 ]
    #   then
    #     echo Loop 3 $x
    #     printf $x >> test.csv
    #     printf , >> test.csv
    # fi
  done
}
# Comment out if just want to add stuff to the end.
printf stars,timesteps,loop0,loop1,loop2,loop3,total,gflop,gb,answer > test.csv
echo >> test.csv

while [ $stars -lt 1000 ]; do
  while [ $timesteps -lt 1000 ]; do
    echo ---$stars stars----
    echo ---$timesteps timesteps---
    printf $stars >> test.csv
    printf , >> test.csv
    printf $timesteps >> test.csv
    printf , >> test.csv
    run_cs257
    let timesteps=timesteps*10
    echo >> test.csv
  done
  timesteps=1
  let stars=stars*10
done 