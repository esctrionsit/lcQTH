import math
import sys
from multiprocessing import Pool, Process, Manager
import traceback

MAFthreahold = 0.3 # lower loci will be deleted, 0 for ignore this filter
gIBDWindowSize = 1000  # Kb
CNVWindowSize = int(sys.argv[1])  # Kb
PhaseWindowSize = int(sys.argv[2])  # Kb
MergeWindowSize = int(sys.argv[3]) # Kb
doBinning = True
doResmoothing = True
doStrictSmooth = False
UnknownRateThreshold = 0.2 # higher loci will be deleted, 0 for ignore this filter

SMlstFile = "Childs.txt"
BinInfoPath = "../05.HMMParsing/05.Bininfo.txt"
SMphaseFolder = "../05.HMMParsing/SMphase.smooth"
OutParentalRatioPath = "ParentalRatio"
OutBinInfoPath = "02.Bininfo.binned.txt"
OutHapPath = "binned"

MAXTHR = int(sys.argv[4])

# chrlen (Mb)
with open("chrlen.txt") as f:
    lines = f.readlines()
chrlen = {}
for line in lines:
    ele = line.strip().split("\t")
    chrlen[ele[0]] = int(ele[1])

# Calculating gIBD <-> Phase scaling factor
BinmergeCount4gIBD = 1
if gIBDWindowSize == PhaseWindowSize:
	pass
elif gIBDWindowSize > PhaseWindowSize:
	BinmergeCount4gIBD = gIBDWindowSize/PhaseWindowSize
	if math.ceil(BinmergeCount4gIBD) - BinmergeCount4gIBD < 0.0001:
		BinmergeCount4gIBD = math.ceil(BinmergeCount4gIBD)
	elif BinmergeCount4gIBD - math.floor(BinmergeCount4gIBD) < 0.0001:
		BinmergeCount4gIBD = math.floor(BinmergeCount4gIBD)
	else:
		print("Error: invaild bin width! Bin size for gIBD should be the integer times of status bin.", file=sys.stderr)
		exit(1)
else:
	print("Error: invaild bin width! Bin size for gIBD should be the integer times of status bin.", file=sys.stderr)
	exit(1)

# Merge Ratio
BinmergeRatio = MergeWindowSize/PhaseWindowSize


# Load Parents gIBD
gIBDdic = {}
with open("../01.ParentsHaplotype/gIBD.txt") as f:
    lines = f.readlines()
header = lines[0].replace("\n", "").split("\t")[1:]
gIBDstatus = lines[1].replace("\n", "").split("\t")[1]
for i in range(len(header)):
    chrom = header[i].split(":")[0]
    if not chrom in gIBDdic:
        gIBDdic[chrom] = []
    gIBDdic[chrom].append(gIBDstatus[i])
gIBDSameMap = {
    "l": None,
    "c": None,
    "z": None
}

# Load Sample List
SMlst = []
with open(SMlstFile) as f:
    lines = f.readlines()
for line in lines:
    SMlst.append(line.replace("\n", ""))

# print("[+] Loading raw bin info", file=sys.stderr)
RawBininfo = {}
RawBinInfoCount = {}
gIBDCursorDic = {}
with open(BinInfoPath) as f:
    lines = f.readlines()
for line in lines:
    ele = line.replace("\n", "").split("\t")
    if not ele[0] in RawBininfo:
        RawBininfo[ele[0]] = []
        gIBDCursorDic[ele[0]] = 0
        RawBinInfoCount[ele[0]] = 1

    if gIBDCursorDic[ele[0]] < len(gIBDdic[ele[0]]) and gIBDdic[ele[0]][gIBDCursorDic[ele[0]]] in gIBDSameMap:
        # Parents are gIBD
        pass
    else:
        # Parents have diversity
        RawBininfo[ele[0]].append([int(ele[1]), int(ele[2])])

    if RawBinInfoCount[ele[0]] == BinmergeCount4gIBD:
        gIBDCursorDic[ele[0]] += 1
        RawBinInfoCount[ele[0]] = 0
    RawBinInfoCount[ele[0]] += 1
del ele
del gIBDCursorDic
del RawBinInfoCount
# print("[+] Loading raw bin info...Done", file=sys.stderr)

# print("[+] Loading raw phase data", file=sys.stderr)
SMData = {}
for SM in SMlst:
    # print("\33[2K\r" + SM, file=sys.stderr, end="")
    SMData[SM] = {}
    with open(SMphaseFolder + "/" + SM + ".txt") as f:
        lines = f.readlines()
    for line in lines:
        ele = line.replace("\n", "").split("\t")
        chrom = ele[0]
        SMData[SM][ele[0]] = ""
        PhaseBinCount = 1
        gIBDCursor = 0
        for i in range(len(ele[1])):
            if gIBDCursor < len(gIBDdic[chrom]) and gIBDdic[chrom][gIBDCursor] in gIBDSameMap:
                # Parents are gIBD
                pass
            else:
                # Parents have diversity
                SMData[SM][ele[0]] += ele[1][i]

            if PhaseBinCount == BinmergeCount4gIBD:
                gIBDCursor += 1
                PhaseBinCount = 0
            PhaseBinCount += 1
del ele
# print("\33[2K\r[+] Loading raw phase data...Done", file=sys.stderr)


# Merge
# print("[+] Start mergeing", file=sys.stderr)
def mergefunc(chrom):
    try:
        newSMData = {}
        for SM in SMData:
            newSMData[SM] = {}
        Bininfo = {}
        Sideout = {}

        # print("\r" + chrom, file=sys.stderr, end="")
        Bininfo[chrom] = []
        for SM in SMData:
            newSMData[SM][chrom] = ""
            Sideout[SM] = ["", ""]
        BinInfoCursor = 0
        for MWS in range(0, chrlen[chrom]*1000, MergeWindowSize):
            BinInfoCursor_tmp = None
            passflag = False
            for SM in SMData:
                BinInfoCursor_tmp = BinInfoCursor
                count = [0, 0, 0] # + ± -
                countmap = {"+":0, "±":1, "-":2}
                while BinInfoCursor_tmp < len(RawBininfo[chrom]):
                    if RawBininfo[chrom][BinInfoCursor_tmp][0] >= MWS and RawBininfo[chrom][BinInfoCursor_tmp][1] <= MWS+MergeWindowSize:
                        count[countmap[SMData[SM][chrom][BinInfoCursor_tmp]]] += 1
                    if RawBininfo[chrom][BinInfoCursor_tmp][0] >= MWS+MergeWindowSize:
                        break
                    BinInfoCursor_tmp += 1
                if sum(count) < BinmergeRatio/2:
                    # No enough valid markers
                    passflag = True
                    break
                else:
                    if count[0] / sum(count) > 0.8:
                        newSMData[SM][chrom] += "+" 
                    elif count[2] / sum(count) > 0.8:
                        newSMData[SM][chrom] += "-"
                    else:
                        newSMData[SM][chrom] += "±"
                Sideout[SM][0] += str(round(count[0]/sum(count)*100,1)) + ","
                Sideout[SM][1] += str(round(count[1]/sum(count)*100,1)) + ","
            if not passflag:
                Bininfo[chrom].append([MWS, (MWS+MergeWindowSize if MWS+MergeWindowSize < chrlen[chrom]*1000 else chrlen[chrom]*1000)])
            BinInfoCursor = BinInfoCursor_tmp
                

        SideoutText = ""
        for SM in Sideout:
            SideoutText += SM + "\t" + Sideout[SM][0][:-1] + "\n"
            SideoutText += SM + "\t" + Sideout[SM][1][:-1] + "\n"
        with open(OutParentalRatioPath + "/" + chrom + ".txt", "w") as f:
            f.write(SideoutText)

        return Bininfo, newSMData
    except Exception as e:
        traceback.print_exc()
        return Bininfo, newSMData
    
Bininfo = {}
newSMData = {}
with Pool(processes = MAXTHR) as proce_pool:
    mulres = []
    for i in range(7):
        for j in ["A", "B", "D"]:
            chrom = "chr" + str(i+1) + j
            mulres.append(proce_pool.apply_async(mergefunc, args=(chrom,)))

    proce_pool.close()
    proce_pool.join()

    for res in mulres:
        g,k = res.get()
        for chrom in g:
            Bininfo[chrom] = g[chrom]
        for SM in k:
            if not SM in newSMData:
                newSMData[SM] = {}
            for chrom in k[SM]:
                newSMData[SM][chrom] = k[SM][chrom]
    del g
    del k
    del mulres




# print("\33[2K\r[+] Start mergeing...Done", file=sys.stderr)
SMData = newSMData
RawBininfo = Bininfo

# Removing rag
if not doStrictSmooth:
    for SM in SMData:
        for chrom in SMData[SM]:
            tmppheseq = list(SMData[SM][chrom])
            for idx in range(len(tmppheseq)):
                if idx == 0 and tmppheseq[idx+1] != tmppheseq[idx]:
                    tmppheseq[idx] = "±"
                elif idx == len(tmppheseq)-1 and tmppheseq[idx-1] != tmppheseq[idx]:
                    tmppheseq[idx] = "±"
                elif idx > 0 and idx < len(tmppheseq)-1 and tmppheseq[idx+1] != tmppheseq[idx] and tmppheseq[idx-1] != tmppheseq[idx]:
                    tmppheseq[idx] = "±"
            SMData[SM][chrom] = "".join(tmppheseq)

# Re-smoothing heter
if doResmoothing:
    # print("[+] Start re-smooth", file=sys.stderr)
    for SM in SMData:
        for chrom in SMData[SM]:
            tmppheseq = list(SMData[SM][chrom])
            heterStart = -1
            for idx in range(len(tmppheseq)):
                if tmppheseq[idx] == "±":
                    if heterStart == -1:
                        heterStart = idx
                    elif RawBininfo[chrom][idx][0] != RawBininfo[chrom][idx-1][1]:
                        if RawBininfo[chrom][heterStart][0] != RawBininfo[chrom][heterStart-1][1]:
                            # ...±±±...
                            laternonheter = -1
                            for sidx in range(idx+1, len(tmppheseq)):
                                if tmppheseq[sidx] != "±":
                                    laternonheter = sidx
                                    break
                            formernonheter = -1
                            for sidx in range(heterStart-1, -1, -1):
                                if tmppheseq[sidx] != "±":
                                    formernonheter = sidx
                                    break
                            if laternonheter == -1 and formernonheter != -1:
                                for sidx in range(heterStart, idx):
                                    tmppheseq[sidx] = tmppheseq[formernonheter]
                            elif laternonheter != -1 and formernonheter == -1:
                                for sidx in range(heterStart, idx):
                                    tmppheseq[sidx] = tmppheseq[laternonheter]
                            elif laternonheter != -1 and formernonheter != -1:
                                if heterStart-formernonheter < laternonheter-idx:
                                    for sidx in range(heterStart, idx):
                                        tmppheseq[sidx] = tmppheseq[formernonheter]
                                elif heterStart-formernonheter < laternonheter-idx:
                                    for sidx in range(heterStart, idx):
                                        tmppheseq[sidx] = tmppheseq[laternonheter]
                            
                        else: 
                            # ???±±±...
                            if tmppheseq[heterStart-1] != "±":
                                for sidx in range(heterStart, idx):
                                    tmppheseq[sidx] = tmppheseq[heterStart-1]

                        heterStart = idx
                else:
                    if heterStart != -1:
                        if heterStart == 0:
                            for sidx in range(heterStart, idx):
                                tmppheseq[sidx] = tmppheseq[idx]
                        else:
                            if tmppheseq[heterStart-1] == tmppheseq[idx]:
                                for sidx in range(heterStart, idx):
                                    tmppheseq[sidx] = tmppheseq[idx]
                            elif RawBininfo[chrom][heterStart][0] != RawBininfo[chrom][heterStart-1][1] and RawBininfo[chrom][idx][0] != RawBininfo[chrom][idx-1][1]:
                                replacehap = "±"
                                if RawBininfo[chrom][heterStart][0]-RawBininfo[chrom][heterStart-1][1] > RawBininfo[chrom][idx][0]-RawBininfo[chrom][idx-1][1]:
                                    replacehap = tmppheseq[idx]
                                elif RawBininfo[chrom][heterStart][0]-RawBininfo[chrom][heterStart-1][1] < RawBininfo[chrom][idx][0]-RawBininfo[chrom][idx-1][1]:
                                    replacehap = tmppheseq[heterStart-1]
                                for sidx in range(heterStart, idx):
                                    tmppheseq[sidx] = replacehap
                            elif RawBininfo[chrom][heterStart][0] != RawBininfo[chrom][heterStart-1][1]:
                                for sidx in range(heterStart, idx):
                                    tmppheseq[sidx] = tmppheseq[idx]
                            elif RawBininfo[chrom][idx][0] != RawBininfo[chrom][idx-1][1]:
                                for sidx in range(heterStart, idx):
                                    tmppheseq[sidx] = tmppheseq[heterStart-1]
                        heterStart = -1
            if heterStart != -1:
                if heterStart != 0:
                    for sidx in range(heterStart, len(tmppheseq)):
                        tmppheseq[sidx] = tmppheseq[heterStart-1]
            SMData[SM][chrom] = "".join(tmppheseq)
            # if SM == "R1" and chrom == "chr5B":
            #     print(SMData[SM][chrom] )
    # print("\33[2K\r[+] Start re-smooth...Done", file=sys.stderr)

# Binning
# print("[+] Start Binning", file=sys.stderr)
def binningfunc(chrom):
    try:
        newSMData = {}
        for SM in SMData:
            newSMData[SM] = {}
        Bininfo = {}
        SMDataKeys = list(SMData.keys())
        # print("\r" + chrom, file=sys.stderr, end="")
        Bininfo[chrom] = ""
        for SM in SMData:
            newSMData[SM][chrom] = ""
        POSidx = 0
        while POSidx < len(SMData[SMDataKeys[0]][chrom]):
            binsize = 0
            while (POSidx + binsize + 1) < len(SMData[SMDataKeys[0]][chrom]):
                InBinSMCount = 0
                if doBinning:
                    for SM in SMData:
                        if SMData[SM][chrom][POSidx] == SMData[SM][chrom][POSidx + binsize + 1]:
                            InBinSMCount += 1
                if InBinSMCount == len(SMData):
                    binsize += 1
                else:
                    break
            # Summary GT
            PlusCount = 0
            MinusCount = 0
            for SM in SMData:
                if SMData[SM][chrom][POSidx] == "+":
                    PlusCount += 1
                elif SMData[SM][chrom][POSidx] == "-":
                    MinusCount += 1
            # print("\r" + str(PlusCount), end="")
            if (MAFthreahold == 0) or (PlusCount >= (PlusCount + MinusCount) * MAFthreahold and PlusCount <= (PlusCount + MinusCount) * (1-MAFthreahold)):
                if (UnknownRateThreshold == 0) or ((len(SMData)-PlusCount-MinusCount)/len(SMData) < UnknownRateThreshold):
                    for SM in SMData:
                        newSMData[SM][chrom] += SMData[SM][chrom][POSidx]
                    Bininfo[chrom] += chrom + "\t" + str(RawBininfo[chrom][POSidx][0]) + "\t" + str(RawBininfo[chrom][POSidx + binsize][1]) + "\n"
            # Update POSidx
            POSidx += binsize + 1
            # print("\r" + chrom + ":" + str(POSidx), end="")
        return Bininfo, newSMData
    except Exception as e:
        traceback.print_exc()
        return Bininfo, newSMData

Bininfo = {}
newSMData = {}
with Pool(processes = MAXTHR) as proce_pool:
    mulres = []
    for i in range(7):
        for j in ["A", "B", "D"]:
            chrom = "chr" + str(i+1) + j
            mulres.append(proce_pool.apply_async(binningfunc, args=(chrom,)))

    proce_pool.close()
    proce_pool.join()

    for res in mulres:
        g,k = res.get()
        for chrom in g:
            Bininfo[chrom] = g[chrom]
        for SM in k:
            if not SM in newSMData:
                newSMData[SM] = {}
            for chrom in k[SM]:
                newSMData[SM][chrom] = k[SM][chrom]
    del g
    del k
    del mulres

# print("\33[2K\r[+] Start Binning...Done", file=sys.stderr)


# Re-smoothing heter
if doResmoothing:
    # print("[+] Start re-smooth", file=sys.stderr)
    for SM in newSMData:
        for chrom in newSMData[SM]:
            tmppheseq = list(newSMData[SM][chrom])
            heterStart = -1
            for idx in range(len(tmppheseq)):
                if tmppheseq[idx] == "±":
                    if heterStart == -1:
                        heterStart = idx
                else:
                    if heterStart != -1:
                        if heterStart == 0:
                            for sidx in range(heterStart, idx):
                                tmppheseq[sidx] = tmppheseq[idx]
                        else:
                            if tmppheseq[heterStart-1] == tmppheseq[idx]:
                                for sidx in range(heterStart, idx):
                                    tmppheseq[sidx] = tmppheseq[idx]
                        heterStart = -1
            if heterStart != -1:
                if heterStart != 0:
                    for sidx in range(heterStart, len(tmppheseq)):
                        tmppheseq[sidx] = tmppheseq[heterStart-1]
            newSMData[SM][chrom] = "".join(tmppheseq)
            # if SM == "R1" and chrom == "chr5B":
            #     print(newSMData[SM][chrom] )
    # print("\33[2K\r[+] Start re-smooth...Done", file=sys.stderr)

with open(OutBinInfoPath, "w") as f:
    for chrom in Bininfo:
        f.write(Bininfo[chrom])
for SM in SMData:
    with open(OutHapPath + "/" + SM + ".txt", "w") as f:
        for chrom in Bininfo:
            f.write(chrom +  "\t" + newSMData[SM][chrom] + "\n")
            