#!/bin/bash

brownie=~/Documents/Tools/brownie/Release/src/brownie
bwa=~/Documents/Tools/bwa/bwa
refalign=~/Documents/Projects/sv_detection/src/align.py
samtools=samtools


ref=${1}
graph=${2}
k=${3}
dir=tempdir
perfect=0
if [ "$#" -ge 4 ]; then
        if [ "${4}" = 'perfect_graph' ]; then
                perfect=1
                if [ "$#" -ge 5 ]; then
                        dir=${5}
                fi
        else
                dir=${4}
        fi
fi

mkdir -p ${dir}/bwa
mkdir -p ${dir}/graph
cp ${ref} ${dir}/bwa

if [ ${perfect} -eq 1 ]; then
        ${brownie} graphConstruction -k ${k} -p ${dir}/graph ${graph} ${graph}
else
        ${brownie} graphCorrection -k ${k} -p ${dir}/graph ${graph}
fi
cd ${dir}
cd graph
	sed -i -e 's#NODE\t##g' DBGraph.fasta
cd ..
cd bwa
	${bwa} index *.fasta
	${bwa} mem -a *.fasta ../graph/DBGraph.fasta > graph.sam
	${samtools} view -Sb graph.sam > graph.bam
	${samtools} sort graph.bam sorted
	${samtools} view -h sorted.bam > sorted.sam
cd ..

mv graph/DBGraph.fasta .
rm -r graph
mv bwa/sorted.sam .
rm -r bwa
/usr/bin/time -v python3 ${refalign} -s sorted.sam -g DBGraph.fasta -k ${k}

cd ..
