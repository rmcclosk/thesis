dependencies
============

This folder contains a Makefile for installing all of the external programs and
libraries on which this project depends. They will all be installed in this
directory, so you should modify your environment as follows.

    export PATH=/path/to/thesis/dependencies/bin:$PATH
    export LIBRARY_PATH=/path/to/thesis/dependencies/lib64:/path/to/thesis/dependencies/lib:$LIBRARY_PATH
    export LD_LIBRARY_PATH=/path/to/thesis/dependencies/lib64:/path/to/thesis/dependencies/lib:$LD_LIBRARY_PATH
    export PKG_CONFIG_PATH=/path/to/thesis/dependencies/lib/pkg-config
    export C_INCLUDE_PATH=/path/to/thesis/dependencies/include:$C_INCLUDE_PATH
    export R_LIBS_USER=/path/to/thesis/dependencies/R/library

Then just type "make". It is assumed that you have a relatively modern version
of gcc installed. 
