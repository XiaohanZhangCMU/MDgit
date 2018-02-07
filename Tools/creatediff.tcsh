#!/bin/tcsh

echo "Begin...\n" > diff.svn

foreach i (*)
  echo "*****************" >> diff.svn
  echo $i >> diff.svn
  echo " " >> diff.svn
  diff $i ../../MD++.svn/trunk/src/$i >> diff.svn
end
  

