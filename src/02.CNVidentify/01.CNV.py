import os
import sys
from multiprocessing import Pool, Process, Manager
from concurrent.futures import ThreadPoolExecutor, wait, FIRST_COMPLETED, ALL_COMPLETED

# Paras
bampath = sys.argv[1]
bedpath = "./bed/"
outpath = './CNV/'
tmppath = sys.argv[2]
CHRlis = sys.argv[3].split(",")
SM_accidlst = sys.argv[4].split(",")
MAXTHR = int(sys.argv[5])
MAXTHR = len(CHRlis) if len(CHRlis) < MAXTHR else MAXTHR

# Check Bam file
checkFlag = False
for SM in SM_accidlst:
    for CHR in CHRlis:
        if (not os.path.exists(bampath + "/" + SM + "." + CHR + ".bam")) or (not (os.path.exists(bampath + "/" + SM + "." + CHR + ".bam.bai") or os.path.exists(bampath + "/" + SM + "." + CHR + ".bam.csi"))):
            print("[E] Bam file of " + SM + " for chromosom " + CHR + " does not exists, or does not indexed. Please note that the bam file should be named as \"<Samplename>.<Chromosome>.bam\"", file=sys.stderr)
            checkFlag = True
if checkFlag:
    exit(1)

def DPextract(chrom, SM, Bampath, Bedpath, Tmppath):
    os.system("rm -f " + Tmppath + "/" + chrom + ".1M.tmp")
    os.system("bedtools coverage -a " + Bedpath + "/" + chrom + ".1M.bed -b " + Bampath + " -counts -sorted > " + Tmppath + "/" + chrom + ".1M.DP")
    os.system("gawk '{if($4 != 0){print int($4/100)*100}else{print 0}}' " + Tmppath + "/" + chrom + ".1M.DP > " + Tmppath + "/" + chrom + ".1M.tmp")

for k in range(len(SM_accidlst)):
    with ThreadPoolExecutor(max_workers=MAXTHR) as t: 
        all_task = []
        for chrom in CHRlis:
            all_task.append(t.submit(DPextract, chrom, SM_accidlst[k], bampath + "/" + SM_accidlst[k] + "." + CHR + ".bam", bedpath, tmppath))
        wait(all_task, return_when=ALL_COMPLETED)
    # Mode of DP
    os.system("cat " + tmppath + "/chr*.1M.tmp | sort | uniq -c | sort -rn -k1 > " + tmppath + "/DPvaluecount")
    with open(tmppath + "/DPvaluecount") as f:
        lines = f.readlines()
    maxcount = 0
    sumDP = 0
    countsum = 0
    for line in lines:
        ele = line.strip().split(" ")
        if int(ele[0]) > maxcount and ele[1] != "0":
            maxcount = int(ele[0])
            countsum = 1
            sumDP = float(ele[1])
        elif int(ele[0]) == maxcount and ele[1] != "0":
            countsum += 1
            sumDP += float(ele[1])
    modeDP = sumDP/countsum
    # Load DP & normalize & phase
    ob = ""
    for chrom in CHRlis:
        with open(tmppath + "/" + chrom + ".1M.tmp") as f:
            lines = f.readlines()
        for line in lines:
            ele = line.strip()
            normedDP = float(ele)/modeDP
            if normedDP < 0.5:
                ele = "del"
            elif normedDP > 1.5:
                ele = "dup"
            else:
                ele = "nor"
            ob += "\t".join([chrom, ele]) + "\n"
    with open(outpath + "/" + SM_accidlst[k] + ".CNV", "w") as f:
        f.write(ob)
    os.system("rm " + tmppath + "/chr*.1M.tmp")
    os.system("rm " + tmppath + "/DPvaluecount")
