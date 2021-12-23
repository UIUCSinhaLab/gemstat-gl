#!/bin/bash
LD_LIBRARY_PATH=/home/grad/samee1/packages/gsl-1.14/lib
export LD_LIBRARY_PATH
BASE=.
DATA=$BASE/data
PAR=$BASE/par
OUT=$BASE/out
SRC=$BASE/src
LOG=$BASE/log
BS_VAL=$2
$SRC/seq2expr -s $DATA/seqs.fa -e $DATA/expr.tab -m $DATA/factors.wtmx -f $DATA/factor_expr.tab -o Direct -oo PGP -fo $OUT/$1.$BS_VAL.out -i $DATA/factor_info.txt -ct 25 -c $DATA/coop.txt -rt 250 -ff $DATA/free_fix.txt -ft $DATA/factor_thr.txt -p $PAR/$1.par -bs $BS_VAL > $LOG/$1.$BS_VAL.log
