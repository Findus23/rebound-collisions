#!/bin/bash

set -x
set -e

gcc -c -fPIC -O3 -Wall heartbeat.c -o heartbeat.o
gcc -L. -shared heartbeat.o -o heartbeat.so -lrebound -Wl,-rpath='./heartbeat' -Wall
