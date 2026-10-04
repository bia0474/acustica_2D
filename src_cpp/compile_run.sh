#!/bin/bash

FLAGS="-O3"

#g++ createParameters.cpp -o parameters
#g++ createGeometry.cpp -o geometry
#g++ createModel.cpp -o velocityModel

cd "$(dirname "$0")" || exit 1

export ACC_NUM_CORES=$(nproc)

nvc++ onda_acustica_2D_4_ordem.cpp $FLAGS -acc=multicore -Minfo=accel -o runAcustica2D || exit 1

time ./runAcustica2D

nvc++ RTM.cpp $FLAGS -acc=multicore -Minfo=accel -o rtm || exit 1

time ./rtm

python3 -i onda_acustica_2D.py