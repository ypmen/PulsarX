# Python Scripts

## Requirements
- multinest
- mpi4py
- pymultinest

## Installation
- python3 setup.py install

## Tools

- ### orbit_solver (solve the binary orbit using Multinest)
  - *usage*: mpiexec -n 6 orbit_solver -y [yaml.conf](examples/orbit_solve/yaml.conf)
  Here is an example input data [test.txt](examples/orbit_solve/test.txt)
  - *tips*: if the elliptic orbit can't fit the data, you can try to fit with circular orbit first to narrow the search range of PB, and then run elliptic orbit search again.