#!/bin/bash

echo "Download KoMPoST from GitHub:"
git clone --depth 1 https://github.com/Hendrik1704/KoMPoST.git
cd KoMPoST
git checkout 2f2d65692f5a6f4c0fe99072f92bb6587b7aa2c4
cd ..

echo "Download MUSIC from GitHub:"
git clone --branch chun_dev --single-branch https://github.com/Hendrik1704/MUSIC
cd MUSIC
git checkout ef77326527929d8db472fa9a7b0c4f55613c6d6c
cd EOS
bash download_hotQCD.sh SMASH_binary
cd ../..

echo "Download iSS from GitHub:"
git clone --depth 1 https://github.com/chunshen1987/iSS.git
cd iSS
git checkout 8b3e43b26dade8252e2671907518a0d111625fd1
cd ..

echo "Download SMASH from GitHub:"
git clone --depth 1 https://github.com/smash-transport/smash.git --branch SMASH-3.2.2

echo "Download Pythia (used for SMASH):"
wget https://pythia.org/download/pythia83/pythia8315.tgz
tar xf pythia8315.tgz

echo "Download stable Eigen library version 3.4.0:"
wget https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz
tar xf eigen-3.4.0.tar.gz

echo "Downloaded all the modules successfully"
