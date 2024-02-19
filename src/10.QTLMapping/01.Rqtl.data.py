import sys
phenoFile = sys.argv[1]

with open(phenoFile) as f:
    lines = f.readlines()
PhenoDic = {}
for line in lines[3:]:
    ele = line.replace("\n", "").split("\t")
    PhenoDic[ele[0]] = ele[1]
PhenoName = lines[2].replace("\n", "").split("\t")[1]


# Load Genetic Map
with open("../09.BuildGeneticMap/GM/GM.full.csv") as f:
    lines = f.readlines()

# Generate csv
s = PhenoName + "," + lines[0] + "," + lines[1] + "," + lines[2]
for i in range(3, len(lines)):
    ele = lines[i].split(",")
    if ele[0] in PhenoDic:
        s += PhenoDic[ele[0]] + "," + lines[i]

with open("data.merged.csv", "w") as f:
    f.write(s)
