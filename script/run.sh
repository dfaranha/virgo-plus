#!/bin/bash

set -x

#touch ../build && rm ../build -rf
mkdir -p ../build
cd ../build
cmake -DCMAKE_BUILD_TYPE=Debug ..
make

./src/virgo_plus_run ../data/test_128_fail.pws

cd ../script
