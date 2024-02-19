#python3
import os
import sys
import math
import numpy
import traceback
from hmmlearn import hmm
from multiprocessing import Pool, Process, Manager

MAXTHR = int(sys.argv[1])
BINwidth = int(sys.argv[2])
CNVBinWidth = int(sys.argv[3])
Parents = sys.argv[4].split(",")

BinmergeCount = 1
if CNVBinWidth == BINwidth:
	pass
elif CNVBinWidth > BINwidth:
	BinmergeCount = CNVBinWidth/BINwidth
	if math.ceil(BinmergeCount) - BinmergeCount < 0.0001:
		BinmergeCount = math.ceil(BinmergeCount)
	elif BinmergeCount - math.floor(BinmergeCount) < 0.0001:
		BinmergeCount = math.floor(BinmergeCount)
	else:
		print("Error: invaild bin width!", file=sys.stderr)
		exit(1)
else:
	print("Error: invaild bin width! Bin size for CNV should NOT smaller than status bin!", file=sys.stderr)
	exit(1)

# print("MAXTHR:", MAXTHR)
# print("BINwidth:", BINwidth)
# print("CNVBinWidth:", CNVBinWidth)
# print("BinmergeCount:", BinmergeCount)
# print("==============\n")


observes = ["+", "*", "±", "/", "-"]
observesMap = {
	"+": 0,
	"*": 1,
	"±": 2,
	"/": 3,
	"-": 4
}
states = ["+", "±", "-"]
n_states = len(states)
n_obs = len(observes)

model = hmm.MultinomialHMM(n_components=n_states, n_iter=20, tol=0.001)
model.startprob_ = numpy.array([0.33, 0.33, 0.33])
model.transmat_ = numpy.array([
	[0.9, 0.08, 0.02],
	[0.15, 0.7, 0.15],
	[0.02, 0.08, 0.9],
])
model.emissionprob_ = numpy.array([
	[0.9, 0.04, 0.04, 4e-05, 0.02],
	[0.05, 0.1, 0.7, 0.1, 0.05],
	[0.02, 4e-05, 0.04, 0.04, 0.9]
])

model.startprob_ = model.startprob_/ sum(model.startprob_)
model.transmat_ = model.transmat_ / model.transmat_.sum(axis=1)[:, numpy.newaxis]
model.emissionprob_ = model.emissionprob_ / model.emissionprob_.sum(axis=1)[:, numpy.newaxis]

# Load Sample List
SMlst = []
with open('Childs.txt') as f:
	lines = f.readlines()
for line in lines:
	SMlst.append(line.replace("\n", ""))

# Load Parents CNV
CNVReMap = {
	"del": "0",
	"nor": "1",
	"dup" : "2"
}
P1CNVDic = {}
with open("../02.CNVidentify/CNV/" + Parents[0] + ".CNV") as f:
	lines = f.readlines()
for line in lines:
	ele = line.replace("\n", "").split("\t")
	if not ele[0] in P1CNVDic:
		P1CNVDic[ele[0]] = []
	P1CNVDic[ele[0]].append(CNVReMap[ele[1]])
P2CNVDic = {}
with open("../02.CNVidentify/CNV/" + Parents[1] + ".CNV") as f:
	lines = f.readlines()
for line in lines:
	ele = line.replace("\n", "").split("\t")
	if not ele[0] in P2CNVDic:
		P2CNVDic[ele[0]] = []
	P2CNVDic[ele[0]].append(CNVReMap[ele[1]])
del lines


# ParentsCNVDiffRegions = {}
# for chrom in P1CNVDic:
# 	for i in range(len(P1CNVDic[chrom])-1):
# 		if P1CNVDic[chrom][i] != P2CNVDic[chrom][i]:
# 			if not chrom in ParentsCNVDiffRegions:
# 				ParentsCNVDiffRegions[chrom] = []
# 			ParentsCNVDiffRegions[chrom].append(i)


def smooth(SM):
	try:
		ChildCNVDic = {}
		with open("../02.CNVidentify/CNV/" + SM + ".CNV") as f:
			lines = f.readlines()
		for line in lines:
			ele = line.replace("\n", "").split("\t")
			if not ele[0] in ChildCNVDic:
				ChildCNVDic[ele[0]] = []
			ChildCNVDic[ele[0]].append(CNVReMap[ele[1]])

		ChildPhaseDic = {}
		with open("SMphase/" + SM + ".txt") as f:
			lines = f.readlines()
		for line in lines:
			ele = line.replace("\n", "").split("\t")
			ChildPhaseDic[ele[0]] = list(ele[1])

		for chrom in ChildPhaseDic:
			# CNVbin >= Statusbin
			CNVidx = 0
			Sbinwidth = 1
			for idx in range(len(ChildPhaseDic[chrom])):
				if P2CNVDic[chrom][CNVidx] != P1CNVDic[chrom][CNVidx] and ChildPhaseDic[chrom][idx] in ("."):
					if ChildCNVDic[chrom][CNVidx] == P1CNVDic[chrom][CNVidx]:
						ChildPhaseDic[chrom][idx] = "-" 
					elif ChildCNVDic[chrom][CNVidx] == P2CNVDic[chrom][CNVidx]:
						ChildPhaseDic[chrom][idx] = "+" 
				if Sbinwidth == BinmergeCount:
					CNVidx += 1
					if CNVidx >= len(ChildCNVDic[chrom]):
						break
					Sbinwidth = 0
				Sbinwidth += 1
		
		# for chrom in ParentsCNVDiffRegions:
		# 	for idx in ParentsCNVDiffRegions[chrom]:
		# 		if ChildCNVDic[chrom][idx] == P1CNVDic[chrom][idx]:
		# 			ChildPhaseDic[chrom][idx] = "+" 
		# 		elif ChildCNVDic[chrom][idx] == P2CNVDic[chrom][idx]:
		# 			ChildPhaseDic[chrom][idx] = "-" 
		# 		else:
		# 			ChildPhaseDic[chrom][idx] = "." 

		SmoothedChildPhase = ""
		for chrom in ChildPhaseDic:
			raw_seq = []
			ValidBinIdxlst = []
			for idx in range(len(ChildPhaseDic[chrom])):
				if ChildPhaseDic[chrom][idx] != ".":
					ValidBinIdxlst.append(idx)
					raw_seq.append(observesMap[ChildPhaseDic[chrom][idx]])
			
			decseq = []
			if len(ValidBinIdxlst) != 0:
				raw_seq = numpy.array([raw_seq]).T
				# print(raw_seq)

				logprob, decseq = model.decode(raw_seq, algorithm="viterbi")
			for decidx in range(len(decseq)):
				ChildPhaseDic[chrom][ValidBinIdxlst[decidx]] = states[decseq[decidx]]

			unknownStart = -1
			for idx in range(len(ChildPhaseDic[chrom])):
				if ChildPhaseDic[chrom][idx] == ".":
					if unknownStart == -1:
						unknownStart = idx
				else:
					if unknownStart != -1:
						if unknownStart == 0:
							for sidx in range(unknownStart, idx):
								ChildPhaseDic[chrom][sidx] = ChildPhaseDic[chrom][idx]
						else:
							if ChildPhaseDic[chrom][unknownStart-1] == ChildPhaseDic[chrom][idx]:
								for sidx in range(unknownStart, idx):
									ChildPhaseDic[chrom][sidx] = ChildPhaseDic[chrom][idx]
							else:
								for sidx in range(unknownStart, idx):
									rsidx = idx-(sidx-unknownStart)
									if sidx > rsidx:
										break
									elif sidx > rsidx:
										ChildPhaseDic[chrom][sidx] = "±"
										break
									else:
										ChildPhaseDic[chrom][sidx] = ChildPhaseDic[chrom][unknownStart-1]
										ChildPhaseDic[chrom][rsidx] = ChildPhaseDic[chrom][idx]
						unknownStart = -1
			if unknownStart != -1:
				if unknownStart != 0:
					for sidx in range(unknownStart, len(ChildPhaseDic[chrom])):
						ChildPhaseDic[chrom][sidx] = ChildPhaseDic[chrom][unknownStart-1]
			SmoothedChildPhase += chrom + "\t" + "".join(ChildPhaseDic[chrom]) + "\n"

		with open("SMphase.smooth/" + SM + ".txt", "w") as f:
			f.write(SmoothedChildPhase)
	except Exception as e:
		traceback.print_exc()

proce_pool = Pool(processes = int(MAXTHR))
mulres = []
for SM in SMlst:
	mulres.append(proce_pool.apply_async(smooth, args=(SM,)))
	#mklevelfile(folder)

proce_pool.close()
proce_pool.join()
