import os
import random
from multiprocessing import Pool, Process, Manager
import traceback

MergeWindowSize = 1000 # Kb
SMlstFile = "../../04.HaplotypeRCR/Bamlst.Childs.txt"
SMPHFile = "../../05.HMMParsing.1M/SMphase.smooth/" # "../06.MarkerScaling/binned/"
MarkerInfoFile = "../../05.HMMParsing.1M/05.Bininfo.txt"  #'../06.MarkerScaling/02.Bininfo.binned.txt'
# MergeWindowSize = 1000 # Kb
# SMlstFile = "../../04.HaplotypeRCR/Bamlst.Childs.txt"
# SMPHFile = "../../06.MarkerScaling/binned/"
# MarkerInfoFile = '../../06.MarkerScaling/02.Bininfo.binned.txt'
Chrom = "chr4B"

def random_str(num):
    H = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789'
    salt = ''
    for i in range(num):
        salt += random.choice(H)
    return salt
    
RawBininfo = {}
with open(MarkerInfoFile) as f:
    lines = f.readlines()
for line in lines:
    ele = line.replace("\n", "").split("\t")
    if not ele[0] in RawBininfo:
        RawBininfo[ele[0]] = []
    RawBininfo[ele[0]].append([int(ele[1]), int(ele[2])])
del ele

##################
SMIDs = []
with open(SMlstFile) as f:
    lines = f.readlines()
for line in lines:
    SMIDs.append(line.replace("\n", ""))

def process_data(SMID):
    # with open("../06.MarkerScaling/binned/" + SMID + ".txt") as f:
    with open(SMPHFile + "/" + SMID + ".txt") as f:
        lines = f.readlines()

    pseudogIBD = SMID + "\t"
    for line in lines:
        ele = line.replace("\n", "").split("\t")
        if ele[0] == Chrom:
            lastend = 0
            count = 0
            for i in range(len(ele[1])):
                binrange = RawBininfo[ele[0]][i]
                while binrange[0] > lastend:
                    pseudogIBD += "."
                    count += 1
                    lastend += MergeWindowSize
                while binrange[1] > lastend:
                    pseudogIBD += ele[1][i]
                    count += 1
                    lastend += MergeWindowSize

    return pseudogIBD + "\n"


with Pool(processes = 30) as proce_pool:  # 限制同时运行的最大进程数
    mulres = []
    for SMID in SMIDs:
        mulres.append(proce_pool.apply_async(process_data, args=(SMID,)))

    proce_pool.close()
    proce_pool.join() # 执行，等待全部进程执行完毕后主进程才执行下一步

    real_res = []
    for res in mulres:
        real_res.append(res.get())
    
    with open("plotdata.txt", "w") as f:
        f.writelines(real_res)