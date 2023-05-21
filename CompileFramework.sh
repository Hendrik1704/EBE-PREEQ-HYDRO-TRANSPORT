#!/bin/bash

echo "Build KoMPoST code for pre-equilibrium stage"
cd KoMPoST && make && cd ..

echo "Build MUSIC code for hydrodynamic stage"
cd MUSIC && mkdir build && cd build && cmake .. && make -Bj && make install && cd .. && rm -r build && cd ..

echo "Build iSS sampler to perform freeze out"
cd iSS && mkdir -p build && cd build && rm -fr * && cmake .. && make -Bj && make install && cd .. && rm -fr build && cd ..

echo "Build Pythia which is used in SMASH for the hadronic afterburner phase"
tar xf pythia8309.tgz
cd pythia8309
./configure --cxx-common='-std=c++11 -march=native -O3 -fPIC'
make -j6 && cd ..

echo "Prepare the Eigen library for SMASH"
tar -xf eigen-3.4.0.tar.gz

echo "Build SMASH as an hadronic afterburner"
cd smash && mkdir build && cd build && cmake .. -DTRY_USE_ROOT=OFF -DPythia_CONFIG_EXECUTABLE=../../pythia8309/bin/pythia8-config -DCMAKE_PREFIX_PATH=../../eigen-3.4.0/ && make -j6
echo "Run SMASH tests"
ctest -j6 && cd ../..

echo "Finished building all the modules successfully"
