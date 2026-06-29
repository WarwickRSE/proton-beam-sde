#!/bin/bash

files=($(ls *_ne_rate.txt))
for file in "${files[@]}" ; do
  ct=$(head -n 1 $file | wc -w)
  echo $ct > $file".bak"
  cat $file >> $file".bak"
  mv $file".bak" $file
done

files=($(ls *_el_ruth_cross_sec.txt))
for file in "${files[@]}" ; do
  ct=$(head -n 1 $file | wc -w)
  ct2=$(head -n 3 $file | tail -n 1 | wc -w)
  echo $ct $ct2 > $file".bak"
  cat $file >> $file".bak"
  mv $file".bak" $file
done
