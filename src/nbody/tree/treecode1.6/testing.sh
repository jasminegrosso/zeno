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
        echo startrun $x
        printf $x >> test.csv
        printf , >> test.csv
    fi
    # Initial Treeforce
    if [ $var -eq 28 ]
      then
        echo Initial Treeforce $x
        printf $x >> test.csv
        printf , >> test.csv
    fi
    # Treeforce Loop
    if [ $var -eq 33 ]
      then
        echo Treeforce loop $x
        printf $x >> test.csv
        printf , >> test.csv
    fi
    # Total
    if [ $var -eq 35 ]
      then
        echo Total $x
        printf $x >> test.csv
        printf , >> test.csv
    fi
  done
}
# Comment out if just want to add stuff to the end.
# printf stars,timesteps,startrun,initial_treeforce,treeforce,totalp > test.csv
# echo >> test.csv

# while [ $stars -lt 10000 ]; do
#   while [ $timesteps -lt 10000 ]; do
#     echo ---$stars stars----
#     echo ---$timesteps timesteps---
#     printf $stars >> test.csv
#     printf , >> test.csv
#     printf $timesteps >> test.csv
#     printf , >> test.csv
#     treeforce
#     let timesteps=timesteps*10
#     echo >> test.csv
#   done
#   timesteps=1
#   let stars=stars*10
# done 
# printf stars,timesteps,loop0,loop1,loop2,loop3,total,gflop,gb,answer > test.csv
# echo >> test.csv


# count=1

# while [ $count -lt 11 ]; do
#   #run_cs257
#   count=$((count+1))
#   stars=10
#   timesteps=10

#   while [ $stars -lt 10000 ]; do
#     # while [ $timesteps -lt 10000 ]; do
#       echo ---$stars stars----
#       echo ---$timesteps timesteps---
#       # printf $stars >> test.csv
#       #printf , >> test.csv
#       #printf $timesteps >> test.csv
#       #printf , >> test.csv
#       run_cs257
#       let timesteps=timesteps*10
#       #echo >> test.csv
#     # done
#     # timesteps=1
#     let stars=stars*10
#   done 

#   echo >> test.csv

# done

count=1

while [ $count -lt 11 ]; do
  run_treecode
  echo >> test.csv
  count=$((count+1))
done 

# while [ $stars -lt 10000 ]; do
#   while [ $timesteps -lt 10000 ]; do
#     echo ---$stars stars----
#     echo ---$timesteps timesteps---
#     printf $stars >> test.csv
#     printf , >> test.csv
#     printf $timesteps >> test.csv
#     printf , >> test.csv
#     run_cs257
#     let timesteps=timesteps*10
#     echo >> test.csv
#   done
#   timesteps=1
#   let stars=stars*10
# done 
>>>>>>> 655aa5e837c3d42d04302652f02614abb2d6497c
