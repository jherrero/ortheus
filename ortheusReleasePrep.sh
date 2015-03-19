#!/bin/sh

programName='ortheus'
baseDir="/Users/benedictpaten/sync/eclipse/"
benLibDir="BenLib"
programDir="ortheus"
buildDir="/tmp/externalBuild"



####this should fail
rm -r ${buildDir}
rm ${programName}Ext.tar
####
mkdir ${buildDir}
dir=`pwd`
cd ${baseDir}
tar -cvf benLib.tar ${benLibDir}
tar -cvf ${programName}.tar ${programDir}
mv benLib.tar ${buildDir}
mv ${programName}.tar ${buildDir}
cd ${buildDir}
tar xf ${programName}.tar 
tar xf benLib.tar 
mv benLib/src/python/benLib ${programName}/src/python/
mv benLib/src/c/benLib/*.[ch] ${programName}/src/c/${programName}/
mv ${programName}/makeFileExternal ${programName}/makefile
rm ${programName}.tar
tar -cvf ${programName}Ext.tar ${programName}
cd ${dir}
mv ${buildDir}/${programName}Ext.tar .
rm -r ${buildDir}