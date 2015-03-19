makefile:

all : clean bin/OrtheusC

#Ortheus make file
 
cxx = gcc

#Release compiler flags
#cflags = -O3 -Wall -funroll-loops -lm

#Debug flags
cflags = -Wall -W -fno-inline -lm -g

#Profile flags
#cflags = -Wall -W -pg -lm


lP = src/c/ortheus
libSources = ${lP}/bioioC.c ${lP}/commonC.c ${lP}/fastCMaths.c ${lP}/hashTableC.c ${lP}/heapC.c ${lP}/substitutionC.c
libHeaders = ${lP}/bioioC.h ${lP}/commonC.h ${lP}/fastCMaths.h ${lP}/hashTableC.h ${lP}/hashTablePrivateC.h ${lP}/heapC.h ${lP}/substitutionC.h

sources = src/c/ortheus/constraintsC.c src/c/ortheus/OrtheusC.c src/c/ortheus/sequenceGraphC.c src/c/ortheus/xyzModelC.c ${libSources}
headers = src/c/ortheus/constraintsC.h src/c/ortheus/sequenceGraphC.h src/c/ortheus/xyzModelC.h ${libHeaders}

clean :
	rm -f bin/OrtheusC
	
bin/OrtheusC : ${sources} ${headers}
	${cxx} ${cflags} -o bin/OrtheusC src/c/ortheus/*.c