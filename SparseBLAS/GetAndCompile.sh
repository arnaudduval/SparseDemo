#!/bin/sh

wget http://www.netlib.org/toms/818.gz
gunzip 818.gz
tail -n+4 818 | /bin/sh
rm 818

sed -i -e "s/FFLAGS = -g/FFLAGS = -O3/g" SOFTWARE/Makefile
sed -i -e "s/FC = f90/FC = gfortran/g" SOFTWARE/Makefile
sed -i -e "s/CC = cc/CC = gcc/g" SOFTWARE/Makefile
