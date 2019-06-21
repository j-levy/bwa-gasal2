#!/bin/bash
sed -i "s/NB_STREAMS ([0-9])/NB_STREAMS ($1)/" src/fastmap.c
cat src/fastmap.c | grep "NB_STREAMS ("
