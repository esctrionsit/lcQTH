#python3
import os
import sys
import math
import numpy
from scipy import stats
from hmmlearn import hmm
from scipy.stats import norm
from multiprocessing import Pool, Process, Manager


# Paras
diffpath = './diff/'
cnvpath = './CNV/'
outpath = './'
SM_accidlst = sys.argv[1].split(",")
CHRlis = sys.argv[2].split(",")
MAXTHR = int(sys.argv[3])

# Pre-defines
SM_vcfidlst = SM_accidlst
gIBDCNVCode = {
	"nor":{
		"del": "y",
		"dup": "b"
	}, "del":{
		"del": "z",
		"dup": "p",
		"nor": "x"
	}, "dup":{
		"del": "1",
		"dup": "c",
		"nor": "a"
	}
}

# HMM Model Paras
AA = 0.8127
AB = 0.1399
AC = 0.0008
AD = 0.0466
BA = 0.0465
BB = 0.9345
BC = 0.0038
BD = 0.0152
CA = 0.0056
CB = 0.0884
CC = 0.8856
CD = 0.0204
DA = 0.2023
DB = 0.1915
DC = 0.0098
DD = 0.5964
AA = math.log10(AA)
AB = math.log10(AB)
AC = math.log10(AC)
AD = math.log10(AD)
BA = math.log10(BA)
BB = math.log10(BB)
BC = math.log10(BC)
BD = math.log10(BD)
CA = math.log10(CA)
CB = math.log10(CB)
CC = math.log10(CC)
CD = math.log10(CD)
DA = math.log10(DA)
DB = math.log10(DB)
DC = math.log10(DC)
DD = math.log10(DD)
# ["A", "B", "C", "D"]
states = ["h", "m", "l", "CNV"]
n_states = 4
n_obs = 4
start_arr_sing = math.log10(0.25)
start_arr = [start_arr_sing, start_arr_sing, start_arr_sing, start_arr_sing]
trans_mat = [[AA,AB,AC,AD],[BA,BB,BC,BD],[CA,CB,CC,CD],[DA,DB,DC,DD]]
means = [3.5, 2.1, 1.1, -3]
covars = [0.21149496, 0.31481414, 0.09464676, 0.04473012]
weight = [math.log10(0.1167319), math.log10(0.6282251), math.log10(0.2550430), math.log10(1)]

# Load CNV status
CNVstatusdic = {}
for SM in SM_accidlst:
	CNVstatusdic[SM] = {}
	with open(cnvpath + "/" + SM + ".CNV") as f:
		lines = f.readlines()
	for line in lines:
		ele = line.strip().split("\t")
		if not ele[0] in CNVstatusdic[SM] :
			CNVstatusdic[SM][ele[0]] = []
		CNVstatusdic[SM][ele[0]].append(ele[1])

# Main function
def call(chrom):
	# Load Diff Matrix and Mask CNV Bins
	with open(diffpath + "/" + chrom + ".combineDiff") as f:
		header = f.readline().strip().split("\t")
		lines = f.readlines()
	# Check header
	SMPAIRidx = 0
	for i in range(len(SM_accidlst)):
		for j in range(i+1, len(SM_accidlst)):
			if (SM_accidlst[i] + "|" + SM_accidlst[j]) != header[SMPAIRidx]:
				print("Error when checking diff matrix header: " + str(i) + "\t" + str(j) + "\t" + str(SMPAIRidx))
				exit(-1)
			SMPAIRidx += 1
	# Phase Diff Matrix
	DiffData = [[] for i in range(len(header))]
	for i in range(len(lines)):
		ele = lines[i].strip().split("\t")
		for j in range(len(ele)):
			SMlst = header[j].split("|")
			if CNVstatusdic[SMlst[0]][chrom][i] != "nor" or CNVstatusdic[SMlst[1]][chrom][i] != "nor":
				DiffData[j].append(math.log10(0.001))
			else:
				DiffData[j].append(math.log10(10+float(ele[j])))
	# Use HMM Model to Phase gIBD
	gIBDStatusSeqs = []
	for SMPAIRidx in range(len(header)):
		raw_seq = DiffData[SMPAIRidx]
		# viterbi
		mat = [[0 for m in range(len(raw_seq))] for n in range(4)]
		matTB = [[0 for m in range(len(raw_seq))] for n in range(4)]
		# Fill in first column
		for m in range(0, len(states)):
			pdf_val = norm.pdf(raw_seq[0], loc=means[m], scale=covars[m])
			if pdf_val == 0.0:
				eq = -64
			else:
				eq = math.log10(pdf_val)
			mat[m][0] = start_arr[m] + eq + weight[m]
		# Fill in the rest of mat, and choose elements of matTB
		for n in range(1, len(raw_seq)):
			for m in range(0, len(states)):
				pdf_val = norm.pdf(raw_seq[n], loc=means[m], scale=covars[m])
				if pdf_val == 0.0:
					eq = -64
				else:
					eq = math.log10(pdf_val)
				mx, mxi = mat[0][n-1] + trans_mat[0][m] + eq + weight[m], 0
				for m_former in range(1, len(states)):
					pr = mat[m_former][n-1] + trans_mat[m_former][m] + eq + weight[m]
					if pr > mx:
						mx, mxi = pr, m_former
				mat[m][n], matTB[m][n] = mx, mxi
		# Find the final state which has the maximal probability
		omx, omxi = mat[0][len(raw_seq)-1], 0
		for m in range(1, len(states)):
			if mat[m][len(raw_seq)-1] > omx:
				omx, omxi = mat[m][len(raw_seq)-1], m
		# Trace for the state sequence
		m, p = omxi, [omxi]
		for n in range(len(raw_seq)-1, 0, -1):
			m = matTB[m][n]
			p.append(m)
		
		logprob = omx
		# log_s += (path_c + "/" + "chr" + str(i+1) + j) + "\n"
		# log_s += str(logprob) + "\n"
		count_lev_to_CNV = 0
		statusseq = ""
		for k in range(len(raw_seq)):
			state = states[p[len(raw_seq)-k-1]]
			if raw_seq[k] > 0:
				if state == "CNV":
					count_lev_to_CNV += 1
					if raw_seq[k] < (means[1]+means[2])/2:
						statusseq += states[2]
					elif raw_seq[k] < (means[0]+means[1])/2:
						statusseq += states[1]
					else:
						statusseq += states[0]
				else:
					statusseq += state
			else:
				SMlst = header[SMPAIRidx].split("|")
				statusseq += gIBDCNVCode[CNVstatusdic[SMlst[0]][chrom][k]][CNVstatusdic[SMlst[1]][chrom][k]]
		# log_s += str(count_lev_to_CNV) + "\n"
		gIBDStatusSeqs.append(statusseq)
	return [chrom, gIBDStatusSeqs]
		
	
proce_pool = Pool(processes = MAXTHR)
mulres = []
for chrom in CHRlis:
	mulres.append(proce_pool.apply_async(call, args=(chrom, )))

proce_pool.close()
proce_pool.join()

# Output
resdic = {}
for res in mulres:
	g = res.get()
	resdic[g[0]] = g[1]
oheader = ["Acc1_vs_Acc2"]
for chrom in CHRlis:
	for i in range(len(resdic[chrom][0])):
		oheader.append(chrom + ":" + str(i))
olines = []
for i in range(len(SM_accidlst)):
	for j in range(i+1, len(SM_accidlst)):
		olines.append(SM_accidlst[i] + "_vs_" + SM_accidlst[j] + "\t")
for chrom in CHRlis:
	SMPAIRidx = 0
	for i in range(len(SM_accidlst)):
		for j in range(i+1, len(SM_accidlst)):
			olines[SMPAIRidx] += resdic[chrom][SMPAIRidx]
			SMPAIRidx += 1
with open(outpath + "/gIBD.tsv", "w") as f:
	f.write("\t".join(oheader) + "\n")
	for line in olines:
		f.write(line + "\n")

