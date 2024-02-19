import sys
CHRlis = sys.argv[1].split(",")

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

# Generate csv
SMDataKeys = list(SMData.keys())
SMGTs = {}
POSID = []
Chrinfo = []
POSinfo = []
ChrName = []
for SM in SMData:
    SMGTs[SM] = []
for chrom in CHRlis:
        chromid = CHRlis.index(chrom) + 1
        for POSidx in range(len(SMData[SMDataKeys[0]][chrom])):
            for SM in SMData:
                if SMData[SM][chrom][POSidx] == "+":
                    SMGTs[SM].append("G")
                elif SMData[SM][chrom][POSidx] == "-":
                    SMGTs[SM].append("T")
                else:
                    SMGTs[SM].append("K")
            Chrinfo.append(str(chromid))
            ChrName.append(chrom)
            # POSinfo.append(str(POSidx*100))
            # POSID.append(chrom + str(POSidx))
            POSinfo.append(str((Bininfo[chrom][POSidx][0]+Bininfo[chrom][POSidx][1])/2/1000))
            POSID.append(chrom + ":" + str(round(Bininfo[chrom][POSidx][0]/1000, 2)) + "-" + str(round(Bininfo[chrom][POSidx][1]/1000, 2)))


s = "\t".join(["rs#", "alleles", "chrom", "pos", "strand", "assembly#", "center", "protLSID", "assayLSID", "panelLSID", "QCcode"] + list(SMGTs.keys())) + "\n"
for i in range(len(POSID)):
    s += "\t".join([
        POSID[i],
        "T/G",
        ChrName[i],
        POSinfo[i],
        "+",
        "NA",
        "NA",
        "NA",
        "NA",
        "NA",
        "NA"
    ])
    for SM in SMGTs:
        s += "\t" + SMGTs[SM][i]
    s += "\n"
with open("02.data.merged.hmp.txt", "w") as f:
    f.write(s)
