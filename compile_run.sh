#!/bin/bash

FLAGS=-O3


g++ onda_acustica_2D_4_ordem.cpp $FLAGS -o runAcustica2D 

./runAcustica2D

python3 -i onda_acustica_2D.py