#!/bin/bash
# local part
l_loc=$1
echo 'l_loc = '$l_loc
echo ' Call psgen'
./psgen ncpp.ini
echo ' Call pswatch'
if [ ! -n "$1" ]; then
 ./pswatch ncpp.ini
else
 ./pswatch ncpp.ini -l $l_loc
fi 
echo ' Call linear'
../tools/linear.x
../tools/qlqgrid.x
../tools/mk_pseudoFile.x

