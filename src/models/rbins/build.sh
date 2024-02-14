#!/bin/bash

# mkdir build
# cd build

cmake /home/flo/WORK/voet_model/fabm -DFABM_HOST=python

# TO BE RUN IN THE BUILD DIRECTORY 
cmake --build . --target install
