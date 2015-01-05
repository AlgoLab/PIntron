#!/bin/bash


VERSION=${pintron_version:master}

git clone https://github.com/AlgoLab/PIntron.git
cd pintron
git checkout $VERSION
make dist
