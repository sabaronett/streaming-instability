#!/bin/bash

echo "***START ORGANIZATION***"
mkdir athdf
mkdir core
mkdir dat
mkdir output
mkdir rst
mkdir video
mkdir xdmf
mv -v *.athdf athdf/
mv -v core.* core/
mv -v *.dat dat/
mv -v *.xdmf xdmf/
mv -v SI.0*.rst rst/
cp -v SI.final.rst rst/
cp -v SI.hst output/
mv -v *.o* output/
echo "***END ORGANIZATION***"
