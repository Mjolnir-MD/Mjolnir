## Performance

__NOTE__
Now, specification of __Mjolnir__ is unstable.
Currently, I have not done tuning.


### Single core performance test

Condition
* 4044 particle(1JZ8), 4 chains
* 100,000 steps with dt = 0.4
* simple Clementi-Go(and EXV) starts from native state
* try 5 times with different random seed

#### CafeMol

compiled with `ifort -cpp -O3 -xHost`

1. 350.50 user 3.28 system 5:56.63 elapsed
2. 351.30 user 3.02 system 5:54.33 elapsed
3. 350.25 user 3.07 system 5:53.32 elapsed
4. 352.35 user 3.13 system 5:55.54 elapsed
5. 350.36 user 3.39 system 5:54.82 elapsed

#### Mjolnir

compiled with `g++-5 -std=c++11 -O3 -mmmx -msse`

1. 263.79 user 2.76 system 4:26.66 elapsed
2. 261.16 user 2.92 system 4:24.13 elapsed
3. 261.78 user 3.02 system 4:24.83 elapsed
4. 269.00 user 2.86 system 4.31.88 elapsed
5. 259.67 user 2.82 system 4:22.52 elapsed

#### result

1.33x faster than CafeMol.

| value    | CafeMol  | Mjolnir |
|:---------|:---------|:--------|
| mean     | 350.952  | 263.08  |
| variance | 3.12308  | 52.555  |
| stddev   | 1.76722  | 7.24948 |
| stderr   | 0.790327 | 3.24207 |


