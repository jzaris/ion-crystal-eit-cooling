#!/bin/bash

set -e

./bin/main_EIT_cooling.exe \
    examples/input/input_EIT_100i_3D_092926_0.txt \
    examples/output/output_EIT_100i_3D_092926_0.txt \
    examples/output/output_EIT_100i_3D_092926_0_1.txt \
    examples/output/int_state_output_EIT_100i_3D_092926_0.txt \
    examples/output/int_state_output_EIT_100i_3D_092926_0_1.txt

