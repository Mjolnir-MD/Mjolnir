#!/bin/bash

# mkdir -p data
[ -d data ] || mkdir data

echo "================================================"
echo " Running a simulation. It would take some time. "
echo "================================================"

./test_G-JFIntegrator test_G-JFIntegrator.toml
/usr/bin/env python3 test_G-JFIntegrator.py
