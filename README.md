# Fine-Mapp-SNP2SE2miRNA

This is an easily-to-use tool for users to detect potential long-range effects of SNPs on miRNA via super-enhancers in human osteoblast. 


Use the follow command to download:

	git clone https://github.com/WeiYang-BAI/Fine-Mapp-SNP2SE2miRNA.git 

Please run the follow command to check the helpdoc first:

	python ./FineMapp_SNP2SE2miRNA.py -H
	
And this will output a flag list:
	
	-H/-h	Show this help-doc.

	-SNP	The INDEPENDENT SNPs file, which consists of at least 
		three tab-separated columns: chromosome, position and P-value.
		Extra columns are allowed, but the first two columns must be
		chromosome (1-22, X, Y) and position (in hg19). 
	-PCol	Columns name of P-value in the SNP input file.
	
	-outPrefix	Prefix for output results.

An example is given:

	python ./FineMapp_SNP2SE2miRNA.py \
	-SNP ./example/example_SNP.txt \
	-PCol p-value \
	-outPrefix ./example/Res
	
The script will output the following five files, run sh example.sh in command line for details:

	-mapping.tsv, results of SNPs mapping to SEs.
	-mappedSE-Bonferroni_sig.tsv, mapped SEs with P < 0.05/1224.
	-Interaction_SE2miRNA.tsv, interactions between miRNA and SE.
	-SE-interacted_pri-miRNA.txt, SE-interacted pri-miRNAs.
	-SE-interacted_mature-miRNA.txt, mature miRNAs.



 Note that the step of long-range interaction identification between Super-Enhancer and miRNA may take a while. 


