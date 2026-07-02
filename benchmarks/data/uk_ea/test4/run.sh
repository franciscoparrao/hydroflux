#!/bin/bash
CASE=ea4-10m

mkdir -p $CASE

lisflood -DEMfile $CASE.dem -dirroot $CASE/acc -acceleration ea4.par
lisflood -DEMfile $CASE.dem -dirroot $CASE/fv1 -fv1 ea4.par
lisflood -DEMfile $CASE.dem -dirroot $CASE/fv1-gpu -fv1 -cuda ea4.par
lisflood -DEMfile $CASE.dem -dirroot $CASE/dg2 -dg2 ea4.par
lisflood -DEMfile $CASE.dem -dirroot $CASE/dg2-gpu -dg2 -cuda ea4.par
