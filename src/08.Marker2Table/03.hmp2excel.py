import sys
from openpyxl import Workbook, load_workbook
from multiprocessing import Pool, Process, Manager

Parents = sys.argv[1].split(",")

with open("02.data.merged.hmp.txt") as f:
    lines = f.readlines()
data = [
    ["T for " + Parents[0] + ", G for " + Parents[1] + ", K for heterozygous and unknown source"],
    [""],
    lines[0].replace("\n", "").split("\t")
]
for line in lines[1:]:
    ele = line.replace("\n", "").split("\t")
    ele[3] = float(ele[3])
    data.append(ele)

wb = Workbook()
ws = wb.active
for line in data:
    ws.append(line)
wb.save("03.hmp2excel.xlsx")