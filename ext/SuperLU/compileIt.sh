#!/bin/bash
export SUPERLUDIR=SuperLU_4.3
cp make.inc.SuperLU_4.3.hermes $SUPERLUDIR/make.inc
cd $SUPERLUDIR
make
