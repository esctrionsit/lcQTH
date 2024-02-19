import sys

Resolution = int(sys.argv[1]) # k
CHRlst = sys.argv[2].split(",")

# Load Sample List
SMlst = []
with open('Childs.txt') as f:
    lines = f.readlines()
for line in lines:
    SMlst.append(line.replace("\n", ""))
SMData = {}
for SM in SMlst:
    SMData[SM] = {}
    with open("SMphase.smooth/" + SM + ".txt") as f:
        lines = f.readlines()
    for line in lines:
        ele = line.replace("\n", "").split("\t")
        SMData[SM][ele[0]] = ele[1]


# Main
newSMData = {}
for SM in SMData:
    newSMData[SM] = {}

Bininfo = {}
SMDataKeys = list(SMData.keys())
for chrom in CHRlst:
    # print("\r" + chrom, end="")
    Bininfo[chrom] = ""
    for SM in SMData:
        newSMData[SM][chrom] = ""
    for POSidx in range(len(SMData[SMDataKeys[0]][chrom])):
        Bininfo[chrom] += chrom + "\t" + str(POSidx*Resolution) + "\t" + str((POSidx + 1)*Resolution) + "\n"
# print("\n")

with open("05.Bininfo.txt", "w") as f:
    for chrom in Bininfo:
        f.write(Bininfo[chrom])
