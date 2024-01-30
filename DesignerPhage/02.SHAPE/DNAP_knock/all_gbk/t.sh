#!/usr/bin/sh
ls CPB* | while read i;
do
	echo -n "${i}" >> result.txt 
	grep "portal" ${i} > test.txt
	if [ -s test.txt ]; then
	echo -n " 1" >> result.txt
	else
	echo -n " 0" >> result.txt
	fi

	grep "terminase" ${i} > test.txt
	if [ -s test.txt ]; then
	echo -n " 1" >> result.txt
	else
	echo -n " 0" >> result.txt
	fi

	grep "tail" ${i} > test.txt
	if [ -s test.txt ]; then
	echo -n " 1" >> result.txt
	else
	echo -n " 0" >> result.txt
	fi

	grep "fiber" ${i} > test.txt
	if [ -s test.txt ]; then
	echo -n " 1" >> result.txt
	else
	echo -n " 0" >> result.txt
	fi


	grep "spike" ${i} > test.txt
	if [ -s test.txt ]; then
	echo -n " 1" >> result.txt
	else
	echo -n " 0" >> result.txt
	fi

	grep "head" ${i} > test.txt
	if [ -s test.txt ]; then
	echo -n " 1" >> result.txt
	else
	echo -n " 0" >> result.txt
	fi

	grep "coat" ${i} > test.txt
	if [ -s test.txt ]; then
	echo -n " 1" >> result.txt
	else
	echo -n " 0" >> result.txt
	fi

	grep "capsid" ${i} > test.txt
	if [ -s test.txt ]; then
	echo " 1" >> result.txt
	else
	echo " 0" >> result.txt
	fi
done



