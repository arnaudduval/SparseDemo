#!/bin/sh

wget http://www.netlib.org/toms/818.gz
gunzip 818.gz
tail -n+4 818 | /bin/sh
