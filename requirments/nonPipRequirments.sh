#!/usr/bin/env sh
wget -O bwa-0.7.5.tar.bz2 http://sourceforge.net/projects/bio-bwa/files/bwa-0.7.5a.tar.bz2/download && tar -xvf bwa-0.7.5.tar.bz2 && rm bwa-0.7.5.tar.bz2 && cd bwa-0.7.5a && make && cp bwa /usr/local/bin && cd .. && rm -r bwa-0.7.5a
apt-get install r-base && Rscript Rrequirements.R