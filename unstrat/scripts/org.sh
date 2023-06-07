#!/bin/bash
#==============================================================================
# org.sh
#
# Organizes typical Athena++ and post-processing outputs into subdirectories by
# file format (extension).
# Reference: https://stackoverflow.com/a/34195247
#
# Author: Stanley A. Baronett
# Updated: 2023-06-07
#==============================================================================

echo "Organizing..."

# Athena++ HDF5 outputs
if compgen -G "*.athdf" > /dev/null; then
    if [[ ! -d athdf ]]; then
        mkdir athdf
    fi
    mv -v *.athdf athdf/
fi

# Athena++ Dust particle data
if compgen -G "*.dat" > /dev/null; then
    if [[ ! -d dat ]]; then
        mkdir dat
    fi
    mv -v *.dat dat/
fi

# Compressed post-processing outputs
if compgen -G "*.npz" > /dev/null; then
    if [[ ! -d npz ]]; then
        mkdir npz
    fi
    mv -v *.npz npz/
fi

# Athena++ restart files
if compgen -G "*.rst" > /dev/null; then
    if [[ ! -d rst ]]; then
        mkdir rst
    fi
    if compgen -G "*.[0-9][0-9][0-9][0-9][0-9].rst" > /dev/null; then
        mv -v *.[0-9][0-9][0-9][0-9][0-9].rst rst/
    fi
    if compgen -G "*.final.rst" > /dev/null; then
        cp -v *.final.rst rst/
    fi
fi

# Athena++ XDMF files (VisIt/ParaView)
if compgen -G "*.xdmf" > /dev/null; then
    if [[ ! -d xdmf ]]; then
        mkdir xdmf
    fi
    mv -v *.xdmf xdmf/
fi

# C++ crash memory dumps
if compgen -G "core.*" > /dev/null; then
    if [[ ! -d core ]]; then
        mkdir core
    fi
    mv -v core.* core/ 
fi

# PBS stdout, Athena++ history files
if compgen -G "*.o*" > /dev/null; then
    if [[ ! -d output ]]; then
        mkdir output
    fi
    mv -v *.o* output/
    if compgen -G "*.hst" > /dev/null; then
        cp -v *.hst output/
    fi
fi

echo "... Done."
