#!/bin/sh

anc=${1}
clo=${2}
seq=${3}
string=${4}
tree=${5}
orthArgs1=${6}
ortheusPath=${7}
passToOrth=1
outFile=${8}


PYTHONPATH="${ortheusPath}/src/python"

echo $PYTHONPATH >&2

echo $1 >&2
echo $2 >&2
echo $3 >&2
echo $4 >&2
echo $5 >&2
echo $6 >&2
echo $7 >&2
echo $8 >&2
	
mv ${seq} scratchy
if [ $passToOrth -eq 1 ]
then
	python ${ortheusPath}/Ortheus.py -b -l '#-w ' -e *.fa -d "${tree}" -f out1.mfa -l "$orthArgs1" -h out1.score -m "java -server -Xmx1700m " > one
fi

if [ $passToOrth -eq 0 ]
then
	python ${ortheusPath}/Ortheus.py -b -l '#-w ' -e *.fa -d "${tree}" -f out1.mfa $orthArgs1 -h out1.score -m "java -server -Xmx1700m " > one
fi

seqNo=`grep '>' out1.mfa | wc -l`
python ${ortheusPath}/src/python/ortheus/utils/CutUpAlignment.py -d out1.mfa -e *.ranges -f $seqNo -g tempFile
mv tempFile out1.mfa

mv scratchy ${seq}

echo '>temp' > ancFile
sed -n "${anc},${anc}p" out1.mfa | sed 's/-//g' >> ancFile
cat ancFile >&2

echo '>clo' > cloFile
sed -n "${clo},${clo}p" out1.mfa | sed 's/-//g' >> cloFile
cat cloFile >&2

python ${ortheusPath}/Ortheus.py -d '(h:0.1, a:0.1);' -e 000.*.fa ${seq} -f tempFile
python ${ortheusPath}/src/python/ortheus/utils/CutUpAlignment.py -d tempFile -e *.ranges -f 3 -g tempFile2
echo '>held' > heldFile
sed -n "6,6p" tempFile2 | sed 's/-//g' >> heldFile
cat heldFile >&2

python ${ortheusPath}/src/python/ortheus/utils/CrossValidate.py -d heldFile -e ancFile -f cloFile -g ${string} >> ${outFile}
rm ancFile
rm cloFile
rm heldFile
rm tempFile
rm tempFile2
grep "total_time" one >> ${outFile}