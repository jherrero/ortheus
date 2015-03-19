#!/bin/sh

ortheusPath=$1

for i in chunk.tar*
do
	python ${ortheusPath}/src/python/ortheus/scripts/utils/shellScripts/splitInThree.py ${i}
done	
