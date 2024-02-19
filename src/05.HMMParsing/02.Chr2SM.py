import sys

CHRlst = sys.argv[1].split(",")

SMData = {}
for chrom in CHRlst:
	with open("data/" + chrom + ".txt") as f:
		lines = f.readlines()
	for line in lines:
		ele = line.replace("\n", "").split("\t")
		if not ele[0] in SMData:
			SMData[ele[0]] = ""
		SMData[ele[0]] += chrom + "\t" + ele[1] + "\n"

for SM in SMData:
	with open("SMphase/" + SM + ".txt", "w") as f:
		f.write(SMData[SM])