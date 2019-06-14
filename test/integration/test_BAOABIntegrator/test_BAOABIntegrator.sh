#!/bin/bash

# mkdir -p data
[ -d data ] || mkdir data

echo "================================================"
echo " Running a simulation. It would take some time. "
echo "================================================"

./test_BAOABIntegrator test_BAOABIntegrator.toml

/usr/bin/env python3 test_BAOABIntegrator.py
