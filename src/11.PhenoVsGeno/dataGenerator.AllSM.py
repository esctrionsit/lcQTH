import math
import sys

phenoFile = sys.argv[1]
targetchrom = sys.argv[2]
targetregion = [int(sys.argv[3]), int(sys.argv[4])] # kb
highlightregion = [int(sys.argv[5])] # kb

with open(phenoFile) as f:
    lines = f.readlines()
PhenoDic = {}
for line in lines[3:]:
    ele = line.replace("\n", "").split("\t")
    PhenoDic[ele[0]] = ele[1]

# Load Sample List
SMlst = []
with open('Childs.txt') as f:
    lines = f.readlines()
for line in lines:
    SMlst.append(line.replace("\n", ""))


SMData = {}
for SM in SMlst:
    SMData[SM] = {}
    with open("../06.MarkerScaling/binned/" + SM + ".txt") as f:
        lines = f.readlines()
    for line in lines:
        ele = line.replace("\n", "").split("\t")
        SMData[SM][ele[0]] = ele[1]

Bininfo = {}
with open("../06.MarkerScaling/02.Bininfo.binned.txt") as f:
    lines = f.readlines()
for line in lines:
    ele = line.replace("\n", "").split("\t")
    if not ele[0] in Bininfo:
        Bininfo[ele[0]] = []
    Bininfo[ele[0]].append([float(ele[1]), float(ele[2])])

# init SM list
SMlst = []
Phenolst = []
for SM in PhenoDic:
    if SM in SMData:
        SMlst.append(SM)
        Phenolst.append(float(PhenoDic[SM]))
Flag = True
while Flag:
    Flag = False
    for i in range(1, len(SMlst)):
        if Phenolst[i] < Phenolst[i-1]:
            Phenolst[i-1], Phenolst[i] = Phenolst[i], Phenolst[i-1]
            SMlst[i-1], SMlst[i] = SMlst[i], SMlst[i-1]
            Flag = True



targetbinidxlst = []
for i in range(len(Bininfo[targetchrom])):
    if Bininfo[targetchrom][i][0] > targetregion[1]:
        break
    elif Bininfo[targetchrom][i][0] >= targetregion[0] and Bininfo[targetchrom][i][1] <= targetregion[1]:
        targetbinidxlst.append(i)

s = "SM"
colplotdata = ""
for idx in targetbinidxlst:
    s += "\t" + str(round((Bininfo[targetchrom][idx][0]+Bininfo[targetchrom][idx][1])/2/1000, 2))
    for hr in highlightregion:
        if Bininfo[targetchrom][idx][0] < hr and Bininfo[targetchrom][idx][1] > hr:
            colplotdata += "1\n"
        else:
            colplotdata += "0\n"

s += "\tPheno\n"
GTmap = {
    "+":"A",
    "-":"B",
    "±":"C",
    "*":"C",
    "/":"C",
    ".":"D"
}

SMlinelst = []
SMcrosslst = []
for SM in SMlst:
    ts = SM
    for idx in targetbinidxlst:
        ts += "\t" + GTmap[SMData[SM][targetchrom][idx]]
    ts += "\t" + PhenoDic[SM] + "\n"
    SMlinelst.append(ts)

# sort
Flag = True
while Flag:
    Flag = False
    for i in range(1, len(SMcrosslst)):
        if SMcrosslst[i] < SMcrosslst[i-1]:
            SMcrosslst[i-1], SMcrosslst[i] = SMcrosslst[i], SMcrosslst[i-1]
            SMlinelst[i-1], SMlinelst[i] = SMlinelst[i], SMlinelst[i-1]
            Flag = True

if len(SMlinelst) == 0:
    print("[I] No Sample left")
    exit(0)

for line in SMlinelst:
    s += line

with open("colAnno.txt", "w") as f:
    f.write(colplotdata)

with open("data.txt", "w") as f:
    f.write(s)
    