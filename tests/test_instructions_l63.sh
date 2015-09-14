#!/bin/bash
set -o nounset

cd ../
make v1
cd examples/lorenz63

./instructions.sh
