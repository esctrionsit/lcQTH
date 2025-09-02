import os
import sys
import random
from concurrent.futures import ThreadPoolExecutor, wait, FIRST_COMPLETED, ALL_COMPLETED
import concurrent.futures
from time import sleep
import traceback

MergeWindowSize = int(sys.argv[1]) # Kb
SMlstFile = "../Childs.txt"
SMPHFile = "../../06.MarkerScaling/binned/"
MarkerInfoFile = '../../06.MarkerScaling/02.Bininfo.binned.txt'
CHRinfo = '../chrlen.txt'
tmppath = sys.argv[2]
MAXTHR = int(sys.argv[3])

os.system("bash init.sh")

chrlen = {}
with open(CHRinfo) as f:
    lines = f.readlines()
for line in lines:
    ele = line.replace("\n", "").split("\t")
    chrlen[ele[0]] = int(ele[1])

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

DataTypeColMap = {
    "+": "0",
    "±": "1",
    "-": "2",
    "N": "3"
}



##################
SMIDs = []
with open(SMlstFile) as f:
    lines = f.readlines()
for line in lines:
    SMIDs.append(line.replace("\n", ""))

def plot(SMID):
    # with open("../06.MarkerScaling/binned/" + SMID + ".txt") as f:
    with open(SMPHFile + "/" + SMID + ".txt") as f:
        lines = f.readlines()

    header = "Acc1_vs_Acc2"
    pseudogIBD = SMID + "_vs_" + "Parents" + "\t"
    for line in lines:
        ele = line.replace("\n", "").split("\t")
        lastend = 0
        count = 0
        for i in range(len(ele[1])):
            binrange = RawBininfo[ele[0]][i]
            while binrange[0] > lastend:
                header += "\t" + ele[0] + ":" + str(count)
                pseudogIBD += DataTypeColMap["N"]
                count += 1
                lastend += MergeWindowSize
            while binrange[1] > lastend:
                header += "\t" + ele[0] + ":" + str(count)
                pseudogIBD += DataTypeColMap[ele[1][i]]
                count += 1
                lastend += MergeWindowSize
        while chrlen[ele[0]] > lastend:
            header += "\t" + ele[0] + ":" + str(count)
            pseudogIBD += DataTypeColMap["N"]
            count += 1
            lastend += MergeWindowSize

    t = random_str(5)
    with open(tmppath + "/" + t + "a", "w") as f:
        f.write(header)
    with open(tmppath + "/" + t + "b", "w") as f:
        f.write(pseudogIBD)

    shell = "node "
    shell += "app.pairwise.readColorIdentify.js "
    shell += tmppath + "/" + t + "a "
    shell += tmppath + "/" + t + "b "
    shell += SMID + " "
    shell += "Inheritplot/" + SMID + " "
    shell += str(MergeWindowSize) + " " 


    os.system(shell)

    os.system("rm " + tmppath + "/" + t + "*")

with ThreadPoolExecutor(max_workers=MAXTHR) as t: 
    all_task = []
    for SMID in SMIDs:
        all_task.append(t.submit(plot, SMID))
    wait(all_task, return_when=ALL_COMPLETED)
    
