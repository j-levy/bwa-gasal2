#!/bin/bash
if [ $# -eq 0 ]; then
	CUR_GPU=$(cat src/bwamem.h | grep "GPU_SELEC" | sed "s/[^01]*//" | sed "s/).*//")
	if [ $CUR_GPU -eq 0 ]; then
		echo "Currently using GPU 0 (Tesla K40c)"
	else
		echo "Currently using GPU 1 (GeForce GTX 750 Ti)"
	fi
else
if [ $# -eq 1 ] && [ $1 -eq 0 -o $1 -eq 1 ] && [ $(cat src/bwamem.h | grep "GPU_SELEC" | sed "s/[^01]*//" | sed "s/).*//") -ne $1 ]; then
	if [ $1 -eq 1 ]; then
		sed -i "s/GPU_SELECT (0)/GPU_SELECT (1)/" src/bwamem.h
		sed -i "s/sm_35/sm_50/" GASAL2/run_all.sh
		echo "GPU set to 1 (GTX 750 Ti)"
	fi
	if [ $1 -eq 0 ]; then
		sed -i "s/GPU_SELECT (1)/GPU_SELECT (0)/" src/bwamem.h
		sed -i "s/sm_50/sm_35/" GASAL2/run_all.sh
		echo "GPU set to 0 (Quadro K40c)"
	fi
	if [ "$(ls obj/ | grep "fastmap.o")" = "fastmap.o" ]; then
		make clean_light
	fi
fi
fi

