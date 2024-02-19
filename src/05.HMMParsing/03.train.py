import os
import sys
import traceback
from multiprocessing import Pool, Process, Manager
from concurrent.futures import ThreadPoolExecutor, wait, FIRST_COMPLETED, ALL_COMPLETED, as_completed

chrlen = {'chr1A': 595, 'chr1B': 690, 'chr1D': 496, 'chr2A': 781, 'chr2B': 802, 'chr2D': 652, 'chr3A': 751, 'chr3B': 831, 'chr3D': 616, 'chr4A': 745, 'chr4B': 674, 'chr4D': 510, 'chr5A': 710, 'chr5B': 714, 'chr5D': 567, 'chr6A': 619, 'chr6B': 721, 'chr6D': 474, 'chr7A': 737, 'chr7B': 751, 'chr7D': 639}


def cal(chrom):
	try:
		with open("data/" + chrom + ".txt") as f:
			PPs = f.readlines()
		count = -1
		for i in range(len(PPs)):
			PPs[i] = PPs[i].replace("\n", "").split("\t")
			if count == -1:
				count = len(PPs[i][1])
			if count != len(PPs[i][1]):
				print(count, len(PPs[i][1]), file=sys.stderr)
				return ""
		
		counts = [0, 0, 0, 0]
		countsmap = {"++":0, "+-":1, "-+":2, "--":3}
		
		for PP in PPs:
			for i in range(len(PP[1])):
				if PP[1][i] == ".":
					continue
				pattern = PP[1][i]
				for j in range(i+1, len(PP[1])):
					if PP[1][j] != ".":
						pattern += PP[1][j]
						break
				if len(pattern) == 2:
					counts[countsmap[pattern]] += 1
		
		SUM = sum(counts)
		s = []
		for c in counts:
			s.append(str(round(c/SUM, 5)))

		return chrom + "\t" + str(s) + "\n"
	except Exception as e:
		traceback.print_exc()
		return ""

# cal("chr1A")

with Pool(processes = 21) as proce_pool:  # 限制同时运行的最大进程数
	mulres = []
	for i in range(7):
		for j in ["A", "B", "D"]:
			chrom = "chr" + str(i+1) + j
			mulres.append(proce_pool.apply_async(cal, args=(chrom,)))

	proce_pool.close()
	proce_pool.join() # 执行，等待全部进程执行完毕后主进程才执行下一步

	s = ""
	for res in mulres:
		s += res.get()
	
	with open("03.o", "w") as f:
		f.write(s)