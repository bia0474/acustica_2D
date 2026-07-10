#!/bin/bash

FLAGS=-O3


g++ onda_acustica_2D_4_ordem.cpp -O3 -fopenmp -o runAcustica2D

time ./runAcustica2D

python3 -i onda_acustica_2D.py