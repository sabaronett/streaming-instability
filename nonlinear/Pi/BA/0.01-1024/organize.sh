#!/bin/bash

echo "***START ORGANIZATION***"
mkdir athdf
mkdir output
mkdir rst
mkdir video
mkdir xdmf
mv -v *.athdf athdf/
mv -v *.xdmf xdmf/
mv -v *.rst rst/
mv -v SI.hst output/
mv -v *.o* output/
echo "***END ORGANIZATION***"
