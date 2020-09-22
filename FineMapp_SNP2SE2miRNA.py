import re
import sys
import os
import time
import pandas as pd
import numpy as np
import scipy
from scipy import stats


def getArgvDict(argv):
		paraDict = {}
		c = 0
		for i in argv:
				if re.match('-',i):
						paraDict[i] = argv[c+1]
				c += 1
		return paraDict

def GetSNP(gFile):
		with open(gFile, 'r') as f:
				allInfo = f.readlines()
		tmp2DSNP = []
		for i in allInfo:
				tmp = re.match(r'(\S+)\t(\d+)', i)
				if tmp:
						tmpSNP = [i.strip(), str(tmp.group(1)), int(tmp.group(2))]
						tmp2DSNP.append(tmpSNP)
		return allInfo[0], tmp2DSNP

def Map2SE(seFile, SNPArray, outPut, gHeader):
		t = open(outPut+'-mapping.tsv', 'w')
		with open(seFile, 'r') as f:
				allInfo = f.readlines()
		t.write(gHeader.strip()+'\t'+allInfo[0])
		for i in allInfo:
				tmp = re.search(r'chr(\S+)\t(\d+)\t(\d+)', i)
				if tmp:
						seChr = str(tmp.group(1))
						seStart = int(tmp.group(2))
						seEnd = int(tmp.group(3))
						for snp in SNPArray:
								if snp[1] == seChr and snp[2] >= seStart and snp[2] <= seEnd:
										t.write(snp[0] + '\t' + i.strip() + '\n')
		t.close()

def Calculate_P(seFile, outPut, PvCol, SEHeader):
		entryDict = {}
		for entry in open(seFile, 'r'):
				match = re.match(r'(.+?K27ac.+?)\tchr(\d+\t\d+\t\d+)\t', entry)
				if match:
						entryDict[match.group(1)] = match.group(2)
		oriRes = pd.read_table(outPut+'-mapping.tsv', header = 0, sep = '\t')
		oriRes.to_csv(outPut+'-TMP', columns = [PvCol,"seID"], index = False, header = False, sep = '\t')
		sePvDict = {}
		SE = list(set(oriRes['seID'].values))
		for s in SE:
				sePvDict[s] = []
		for x in open(outPut+'-TMP', 'r'):
				match = re.match(r'(\S+)\t(\S+)', x)
				if match:
						sePvDict[match.group(2)].append(float(match.group(1)))
		t = open(outPut+'-mappedSE-Bonferroni_sig.tsv', 'w')
		t.write(SEHeader)
		for s in SE:
				X2 = sum(np.log(np.array(sePvDict[s])))*(-2)
				N = len(sePvDict[s])
				FD = len(sePvDict[s])*2
				Pv = stats.distributions.chi2.sf(X2, FD)
				if Pv <= 0.05/1224 and N >= 10:
						t.write(s + '\t' + entryDict.get(s, '-\t-\t-') + '\t' + str(N) + '\t' + str(X2) + '\t' + str(FD) + '\t' + str(Pv)+ '\n')
		t.close()
		if os.path.exists(outPut+'-TMP'):
				os.remove(outPut+'-TMP')

def FindInteraction(miRNAFile, outPut, HiCFile):
		with open(miRNAFile, 'r') as f:
				allMiRNA = f.readlines()
		with open(outPut+'-mappedSE-Bonferroni_sig.tsv', 'r') as f:
				allSE = f.readlines()
		with open(HiCFile, 'r') as f:
				allInterR = f.readlines()
		t = open(outPut+'-Interaction_SE2miRNA.tsv', 'w')
		t.write(allMiRNA[0].strip()+'\t'+allSE[0].strip()+'\t'+allInterR[0].strip()+'\t'+'distance_SE2miRNA'+'\n')
		for miRNA in allMiRNA:
				tmpMiRNA = re.match(r'(chr\S+)\t(\d+)\t(\d+)\t', miRNA)
				if tmpMiRNA:
						mirC = tmpMiRNA.group(1)
						mirS = int(tmpMiRNA.group(2))
						mirE = int(tmpMiRNA.group(3))
						for SE in allSE:
								tmpSE = re.match(r'.+?\t(\S+)\t(\d+)\t(\d+)\t', SE)
								if tmpSE:
										seC = 'chr'+tmpSE.group(1)
										seS = int(tmpSE.group(2))
										seE = int(tmpSE.group(3))
										if mirC != seC:
												continue
										else:
												for InterR in allInterR:
														tmpInterR = re.match(r'(chr\S+),(\d+),(\d+)\t(chr\S+),(\d+),(\d+)\t', InterR)
														if tmpInterR:
																r1C = tmpInterR.group(1)
																r1S = int(tmpInterR.group(2))
																r1E = int(tmpInterR.group(3))
																r2C = tmpInterR.group(4)
																r2S = int(tmpInterR.group(5))
																r2E = int(tmpInterR.group(6))
																if r1C == seC:
																		if seS > r1E or seE < r1S or mirS > r2E or mirE < r2S:
																				uFool = 'try below'
																		else:
																				if mirE <= seS:
																						distance = seS - mirE
																				elif mirS >= seE:
																						distance = mirS - seE
																				else:
																						distance = 'overlapped'
																				t.write(miRNA[0:-1]+'\t'+SE[0:-1]+'\t'+InterR[0:-1]+'\t'+str(distance)+'\n')
																if r1C == mirC:
																		if seS > r2E or seE < r2S or mirS > r1E or mirE < r1S:
																				uFool = 'try next pair'
																		else:
																				if mirE <= seS:
																						distance = seS - mirE
																				elif mirS >= seE:
																						distance = mirS - seE
																				else:
																						distance = 'overlapped'
																				t.write(miRNA[0:-1]+'\t'+SE[0:-1]+'\t'+InterR[0:-1]+'\t'+str(distance)+'\n')
		t.close()

def ExtracMiRNA(outPut, miRNAFile):
		with open(outPut+'-Interaction_SE2miRNA.tsv', 'r') as f:
				allPri = f.read()
		with open(miRNAFile, 'r') as f:
				allMirna = f.read()
		match = re.findall(r'(hsa-\S+-\S+)\t', allPri)
		t1 = open(outPut+'-SE-interacted_pri-miRNA.txt', 'w')
		t2 = open(outPut+'-SE-interacted_mature-miRNA.txt', 'w')
		if match:
				for i in list(set(match)):
						t1.write(i + '\n')
						tmpEx = re.search(r'ID=(\S+?);.*?Name='+i, allMirna)
						if tmpEx:
								miRNA = re.findall(r'Name=(\S+?);Derives_from='+tmpEx.group(1)+'\n', allMirna)
								for im in miRNA:
										t2.write(im+'\n')
		t1.close()
		t2.close()



if __name__ == '__main__':
		helpDoc = '''
-H/-h	Show this help message and exit.

-SNP	The INDEPENDENT SNPs file, which consists of at least 
	three tab-separated columns: chromosome, position and P-value.
	Extra columns are allowed, but the first two columns must be
	chromosome (1-22, X, Y) and position (in hg19). 
-PCol	Columns name of P-value in the SNP input file.
-outPrefix	Prefix for output results.

Note that the step of long-range interaction identification
between Super-Enhancer and miRNA may take a while.

The script will output the following five files:

	-mapping.tsv, results of SNPs mapping to SEs.
	-mappedSE-Bonferroni_sig.tsv, mapped SEs with P < 0.05/1224.
	-Interaction_SE2miRNA.tsv, interactions between miRNA and SE.
	-SE-interacted_pri-miRNA.txt, SE-interacted pri-miRNAs.
	-SE-interacted_mature-miRNA.txt, mature miRNAs.
		'''
		argv = sys.argv
		if re.search('-H|-h', str(argv)):
				print(helpDoc)
				sys.exit()
		try:
				paraDict = getArgvDict(argv)
				GWASFile = paraDict['-SNP']
				PvalueCol = paraDict['-PCol']
				OUTPrefix = paraDict['-outPrefix']
		except KeyError as e:
				print('ERROR: Incomplete or invalid arg '+str(e)+' ! Please check your input, or use -H/-h flag to get help.')
		else:
				SEFile = './data/GSM733697_rep1_K27ac_SE_1224.tsv'
				miRNAFile = './data/miRNA_miRBase_v22.txt'
				miRBase = './data/miRBaseV22.hsa37.liftOver.res'
				HiCFile = './data/BMP2_both_frag_washU_text.txt'
		print('Collect SNPs...')
		gHeader, SNP2DArray = GetSNP(GWASFile)
		print('Mapping to SEs...')
		Map2SE(SEFile, SNP2DArray, OUTPrefix, gHeader)
		SEHeader = 'GSM733697-SE'+'\t'+'CHR_SE'+'\t'+'startPOS_SE'+'\t'+'endPOS_SE'+'\t'+'N_SNP'+'\t'+'X2'+'\t'+'freedom_degree'+'\t'+'Pvalue'+'\n'
		print('Get significantly mapped SEs...')
		Calculate_P(SEFile, OUTPrefix, PvalueCol, SEHeader)
		print('Find interactions between miRNA and SE...')
		FindInteraction(miRNAFile, OUTPrefix, HiCFile)
		ExtracMiRNA(OUTPrefix, miRBase)
		print('Done.')
		localtime = time.asctime(time.localtime(time.time()))
		print('Done. Finish time: '+ localtime + '\n')

