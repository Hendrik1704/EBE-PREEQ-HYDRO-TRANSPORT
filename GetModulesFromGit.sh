#!/bin/bash

echo "Download KoMPoST from GitHub:"
git clone --depth 1 https://github.com/KMPST/KoMPoST.git

echo "Download MUSIC from GitHub:"
git clone --depth 1 https://github.com/MUSIC-fluid/MUSIC.git

echo "Download iSS from GitHub:"
git clone --depth 1 https://github.com/chunshen1987/iSS.git

echo "Download SMASH from GitHub:"
git clone --depth 1 https://github.com/smash-transport/smash.git --branch SMASH-3.0

echo "Download Pythia (used for SMASH):"
wget https://pythia.org/download/pythia83/pythia8309.tgz

echo "Download stable Eigen library version 3.4.0:"
wget https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz

echo "Downloaded all the modules successfully"
