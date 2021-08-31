#!/bin/bash

echo "***START ORGANIZATION***"
mkdir athdf
mkdir output
mkdir rst
mkdir video
mkdir xdmf
mv -v *.athdf athdf/
mv -v *.o* output/
mv -v SI.hst output/
mv -v *.rst rst/
mv -v *.xdmf xdmf/
echo "***END ORGANIZATION***"
