#!/bin/bash

echo "Organizing..."

if [[ -f *.athdf ]]
then
    if [[ ! -d athdf ]]
    then
        mkdir athdf
    fi
    mv -v *.athdf athdf/
fi

if [[ -f core.* ]]
then
    if [[ ! -d core ]]
    then
        mkdir core
    fi
    mv -v core.* core/ 
fi

if [[ -f *.dat ]]
then
    if [[ ! -d dat ]]
    then
        mkdir dat
    fi
    mv -v *.dat dat/
fi

if [[ -f *.o* ]]
then
    if [[ ! -d output ]]
    then
        mkdir output
    fi
    cp -v *.hst output/
    mv -v *.o* output/
fi

if [[ -f *.rst ]]
then
    if [[ ! -d rst ]]
    then
        mkdir rst
    fi
    mv -v 0*.rst rst/
    if [[ -f *.final.rst ]]
    then
        cp -v *.final.rst rst/
    fi
fi

if [[ ! -d video ]]
then
    mkdir video
fi

if [[ -f *.xdmf]]
then
    if [[ ! -d xdmf ]]
    then
        mkdir xdmf
    fi
    mv -v *.xdmf xdmf/
fi

echo "... Done."
