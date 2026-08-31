#!/bin/bash

FLAGS=-O3

#g++ createParameters.cpp -o parameters
#g++ createGeometry.cpp -o geometry
#g++ createModel.cpp -o velocityModel
g++ onda_acustica_2D_4_ordem.cpp -O3 -fopenmp -o runAcustica2D

time ./runAcustica2D

g++ RTM.cpp -O3 -fopenmp -o rtm

time ./rtm

python3 -i onda_acustica_2D.py