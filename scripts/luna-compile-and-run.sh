#!/bin/bash

LIBRARY_DIR=../cmake-build-debug/projects/Luna

../cmake-build-debug/projects/Core/luna/LUNA_RTS --fa ../projects/Luna/Scheme.ja --library $LIBRARY_DIR/libLuna.so
