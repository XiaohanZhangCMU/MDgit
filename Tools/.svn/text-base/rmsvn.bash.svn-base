#!/bin/bash
# This file removes a hidden directory .svn from every subdirectories of the current directory.
# Aug/03/2007 Keonwook Kang

for i in `ls -AR ./*` 
do
########## Find Directory ##########
  numchar=${#i}                  # number of charcther in i
  lastch=${i:$[$numchar - 1]:1}  # the last character of i
  if [ "$lastch" = ":" ]         # if the last character is ":",
  then                           # i is a directory name. 
    dir=${i%:}    # remove ":" from i and set the output as "dir" 
########## Find .svn ##########
    last4ch=${dir:$[$numchar - 5]:4} # the last 4 charcters of "dir"
#    echo $last4ch
    if [ "$last4ch" = ".svn" ]   # if i is ".svn",
    then
      if [ -d $dir ]             # if dir exists,
      then
        rm -rf $dir
        echo .svn removed from $dir directory.
      fi
    fi
  fi
done
