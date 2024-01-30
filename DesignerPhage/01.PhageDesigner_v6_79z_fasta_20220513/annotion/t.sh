#!/usr/bin/sh
ls -d * | head -79| while read i;
do
	echo -n "${i}" >> result.txt
	grep "DNA polymerase" ${i}/02.blastp/*.txt > test.txt
	if [ -s test.txt ]; then
	echo -n " 1" >> result.txt
	else
	echo -n " 0" >> result.txt
	fi
	
	grep "tRNA" ${i}/02.blastp/*.txt > test.txt
	if [ -s test.txt ]; then
	echo -n " 1" >> result.txt
	else
	echo -n " 0" >> result.txt
	fi

	grep "portal" ${i}/02.blastp/*.txt > test.txt
	if [ -s test.txt ]; then
	echo -n " 1" >> result.txt
	else
	echo -n " 0" >> result.txt
	fi

	grep "terminase" ${i}/02.blastp/*.txt > test.txt
	if [ -s test.txt ]; then
	echo -n " 1" >> result.txt
	else
	echo -n " 0" >> result.txt
	fi

	grep "tail" ${i}/02.blastp/*.txt > test.txt
	if [ -s test.txt ]; then
	echo -n " 1" >> result.txt
	else
	echo -n " 0" >> result.txt
	fi

	grep "fiber" ${i}/02.blastp/*.txt > test.txt
	if [ -s test.txt ]; then
	echo -n " 1" >> result.txt
	else
	echo -n " 0" >> result.txt
	fi


	grep "spike" ${i}/02.blastp/*.txt > test.txt
	if [ -s test.txt ]; then
	echo -n " 1" >> result.txt
	else
	echo -n " 0" >> result.txt
	fi

	grep "head" ${i}/02.blastp/*.txt > test.txt
	if [ -s test.txt ]; then
	echo -n " 1" >> result.txt
	else
	echo -n " 0" >> result.txt
	fi

	grep "coat" ${i}/02.blastp/*.txt > test.txt
	if [ -s test.txt ]; then
	echo -n " 1" >> result.txt
	else
	echo -n " 0" >> result.txt
	fi

	grep "capsid" ${i}/02.blastp/*.txt > test.txt
	if [ -s test.txt ]; then
	echo " 1" >> result.txt
	else
	echo " 0" >> result.txt
	fi
done



