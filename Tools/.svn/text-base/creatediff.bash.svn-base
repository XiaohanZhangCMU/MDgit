#!/bin/bash
# This file does diff between MD++ and MD++.svn.
# $ cd MD++/src
# $ ../Tools/creatediff.bash

echo "Begin..." > diff.svn
echo "" >> diff.svn

#for i in `ls `
for i in *
do
  echo "*****************" >> diff.svn
  echo $i >> diff.svn
  echo " " >> diff.svn
  diff $i ../../MD++.svn/trunk/src/$i >> diff.svn
done

unset i

