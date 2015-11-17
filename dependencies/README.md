dependencies
============

This folder contains a Makefile for installing all of the external programs and
libraries on which this project depends. They will all be installed in this
directory, so you should modify your environment as follows.

    export PATH=/path/to/thesis/dependencies/bin:$PATH
    export LIBRARY_PATH=/path/to/thesis/dependencies/bin:$LIBRARY_PATH
    export LD_LIBRARY_PATH=/path/to/thesis/dependencies/bin:$LD_LIBRARY_PATH
    export C_INCLUDE_PATH=/path/to/thesis/dependencies/include:$C_INCLUDE_PATH

Then just type "make".
