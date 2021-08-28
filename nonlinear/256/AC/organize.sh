#!/bin/bash

echo "***START ORGANIZATION***"
mkdir rst
mkdir xdmf
mkdir athdf
mv -v *.rst rst/
mv -v *.xdmf xdmf/
mv -v *.athdf athdf/
echo "***END ORGANIZATION***"
