#!/bin/bash

# mkdir -p data
[ -d data ] || mkdir data

echo "================================================"
echo " Running a simulation. It would take some time. "
echo "================================================"

./test_g-BAOABIntegrator test_g-BAOABIntegrator.toml

/usr/bin/env python3 check_constraint.py
/usr/bin/env python3 check_constraint.py
echo "passed check constraint.py"
/usr/bin/env python3 test_g-BAOABIntegrator.py
