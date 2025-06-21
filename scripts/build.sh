#!/bin/bash

RHOMBUS_ENABLE_TESTS=ON

mkdir -p build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release -DRHOMBUS_BUILD_TESTS=${RHOMBUS_ENABLE_TESTS}
make -j4