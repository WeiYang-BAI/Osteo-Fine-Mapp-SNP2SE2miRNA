#!/bin/sh

python \
	./FineMapp_SNP2SE2miRNA.py \
	-SNP ./example/example_SNP.txt \
	-PCol p-value \
	-outPrefix ./example/Res
