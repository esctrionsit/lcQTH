import os
import sys
from multiprocessing import Pool, Process, Manager

# Paras
vcfpath = sys.argv[1]
vcfend = sys.argv[2]
bedpath = "./bed/"
outpath = './diff/'
tmppath = sys.argv[3]
SM_accidlst = sys.argv[4].split(",")
CHRlis = sys.argv[5].split(",")
MAXTHR = int(sys.argv[6])

# Pre-defines
SM_vcfidlst = SM_accidlst

# Check vcf file
checkFlag = False
for CHR in CHRlis:
    if (not os.path.exists(vcfpath + "/" + CHR + vcfend)) or ((not os.path.exists(vcfpath + "/" + CHR + vcfend + ".csi") and (not os.path.exists(vcfpath + "/" + CHR + vcfend + ".tbi")))):
        print("[E] Vcf file of " + CHR + " does not exists, or does not indexed. Please note that the vcf file should be named as \"<Chromosome>.<VcfFileSuffixes>\"", file=sys.stderr)
        checkFlag = True
if checkFlag:
    exit(1)

# Cal DSR
def calDSR(chrom):
    combineDSRheader = []
    combineDSRcontent = ""
    for i in range(len(SM_accidlst)):
        for j in range(i+1, len(SM_accidlst)):
            combineDSRheader.append(SM_accidlst[i] + "|" + SM_accidlst[j])

    with open(bedpath + "/" + chrom + ".1M.bed") as f:
        bedlines = f.readlines()
    for line in bedlines:
        bedele = line.strip().split("\t")
        # Generate SNP Matrix
        region = bedele[0] + ":" + bedele[1] + "-" + bedele[2]
        vcffilepath = vcfpath + "/" + chrom + vcfend
        shell = "bcftools view -v snps -r " + region + " " + vcffilepath + " -s " + ",".join(SM_vcfidlst) + " | bcftools query -f '[%GT:%DP:%GQ\\t]\\n'"
        SNPs = os.popen(shell).readlines()
        SNPMAT = "\t".join(SM_accidlst) + "\n"
        for SNP in SNPs:
            ele = SNP.strip().split("\t")
            tmpSNPline = []
            for SMSNP in ele:
                SNPstatus = SMSNP.split(":")
                GT = SNPstatus[0]
                DP = SNPstatus[1]
                GQ = SNPstatus[2]
                if GT == "./." or DP == "." or GQ == "." or int(DP)<3 or int(GQ)<8 or int(DP)>99:
                    tmpSNPline.append("-1")
                else:
                    GT = GT.split("/")
                    if GT[0] == "0" and GT[1] == "0":
                        tmpSNPline.append("0")
                    elif GT[0] == GT[1]:
                        tmpSNPline.append("2")
                    else:
                        tmpSNPline.append("1")
            SNPMAT += "\t".join(tmpSNPline) + "\n"
        with open(tmppath + "/" + chrom + ".SNPMAT", "w") as f:
            f.write(SNPMAT)
        # Calculate DSR
        os.system(
            "src/dsrdist " + 
            tmppath + "/" + chrom + ".SNPMAT " + 
            tmppath + "/" + chrom + ".dsr "
        )
        # Load DSR Matrix
        with open(tmppath + "/" + chrom + ".dsr") as f:
            DSRlines = f.readlines()
        # Combine DSR Matrix
        combineDSRcontentline = []
        for i in range(len(SM_accidlst)):
            ele = DSRlines[i].replace("\n", "").split("\t")
            for j in range(i+1, len(SM_accidlst)):
                combineDSRcontentline.append(ele[j])
        combineDSRcontent += "\t".join(combineDSRcontentline) + "\n"
    with open(outpath + "/" + chrom + ".combineDiff", "w") as f:
        f.write("\t".join(combineDSRheader) + "\n")
        f.write(combineDSRcontent)




proce_pool = Pool(processes = MAXTHR) # 限制同时运行的最大进程数
mulres = [] # 可保存返回值，顺序不定
for chrom in CHRlis:
    mulres.append(proce_pool.apply_async(calDSR, args=(chrom,)))

proce_pool.close()
proce_pool.join()
