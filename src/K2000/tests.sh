#!/bin/bash

echo "*** Running Command: ./run_K2000.sh tests/numbers.txt tests/dbg4.fa 240 tests/afac.gfa tests/afac.fasta***"

./run_K2000.sh tests/dbg_paths.txt tests/unitigs.fa 240 tests/afac.gfa tests/afac.fasta
if [ $? -ne 0 ] ; then
echo "*** Test: FAILURE on command ./run_K2000.sh tests/dbg_paths.txt tests/unitigs.fa 240 tests/afac.gfa tests/afac.fasta"
rm -f compacted_unitigs.txt tests/afac.gfa tests/afac.fasta
exit 1
fi

diff tests/dbg_paths.txt_compacted  tests/dbg_paths.txt_compacted_ref
if [ $? -ne 0 ] ; then
echo "*** Test: FAILURE on diff compacted_numbers"
rm -f tests/dbg_paths.txt_compacted tests/afac.gfa tests/afac.fasta
exit 1
fi

diff tests/afac.gfa tests/out.gfa
if [ $? -ne 0 ] ; then
echo "*** Test: FAILURE on diff gfa"
rm -f tests/dbg_paths.txt_compacted tests/afac.gfa tests/afac.fasta
exit 1
fi


diff tests/afac.fasta tests/out.fasta
if [ $? -ne 0 ] ; then
echo "*** Test: FAILURE on diff fasta"
rm -f tests/dbg_paths.txt_compacted tests/afac.gfa tests/afac.fasta
exit 1
fi

echo "*** Test: OK"
rm -f tests/dbg_paths.txt_compacted tests/afac.gfa tests/afac.fasta
exit 0