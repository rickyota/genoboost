#!/bin/bash

fwrap="./wrapper.h"
rm -rf $fwrap


#echo "#pragma once" >> $fwrap


for f in $(find . -name '*.h') ; do
	if [[ $f == "./wrapper.h" ]]; then 
		echo "skip $f"
	#elif [[ $f == "*/avx512/*.h" ]]; then 
		#echo "skip avx512 $f"
	else
		h="#include \"${f:2}\""
		echo $h >> $fwrap

	fi
done
