
# Load Sample List
SMlst = []
with open('../01.MergeGvcf/Bamlst.Childs.txt') as f:
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

# Generate csv
SMDataKeys = list(SMData.keys())
SMGTs = {}
POSID = []
Chrinfo = []
POSinfo = []
for SM in SMData:
    SMGTs[SM] = []
for i in range(7):
    for j in ["A", "B", "D"]:
        chrom = "chr" + str(i+1) + j
        chromid = i*3 +  ["A", "B", "D"].index(j) + 1
        for POSidx in range(len(SMData[SMDataKeys[0]][chrom])):
            for SM in SMData:
                if SMData[SM][chrom][POSidx] == "+":
                    SMGTs[SM].append("GG")
                elif SMData[SM][chrom][POSidx] == "-":
                    SMGTs[SM].append("TT")
                else:
                    SMGTs[SM].append("GT")
            Chrinfo.append(str(chromid))
            # POSinfo.append(str(POSidx*100))
            # POSID.append(chrom + str(POSidx))
            POSinfo.append(str((Bininfo[chrom][POSidx][0]+Bininfo[chrom][POSidx][1])/2/1000))
            POSID.append(chrom + str(POSidx))


s = ",".join(["wheat_ID"] + POSID) + "\n"
s += ",".join([""] + Chrinfo) + "\n"
s += ",".join([""] + POSinfo) + "\n"
for SM in SMGTs:
    s += ",".join([SM] + SMGTs[SM]) + "\n"
with open("GM.1M/data.merged.csv", "w") as f:
    f.write(s)
