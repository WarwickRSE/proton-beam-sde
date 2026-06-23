#!/bin/bash

files=($(ls *.txt))
for file in "${files[@]}" ; do
  ct=$(head -n 1 $file | wc -w)
  echo $ct > $file".bak"
  cat $file >> $file".bak"
  mv $file".bak" $file
done