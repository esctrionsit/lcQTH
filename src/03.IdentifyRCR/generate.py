import os
import sys
import traceback
import subprocess
from multiprocessing import Pool, Process, Manager

SMID = sys.argv[1]
AlignPath = sys.argv[2]
tmppath = sys.argv[3]
MAXTHR = int(sys.argv[4])
CHRlis = sys.argv[5].split(",")

# Check Bam files
checkFlag = False
for CHR in CHRlis:
	if (not os.path.exists(AlignPath + "/" + SMID + "." + CHR + ".bam")) or (not (os.path.exists(AlignPath + "/" + SMID + "." + CHR + ".bam.bai") or os.path.exists(AlignPath + "/" + SMID + "." + CHR + ".bam.csi"))):
		print("[E] Bam file of " + SMID + " for chromosom " + CHR + " does not exists, or does not indexed. Please note that the bam file should be named as \"<Samplename>.<Chromosome>.bam\"", file=sys.stderr)
		checkFlag = True
if checkFlag:
	exit(1)

def process_data(SM, chrom):
	try:
		if os.path.exists("ReadCovRegionList.tmp/" + SM + "." + chrom + ".txt"):
			return
		
		ReturnData = [chrom, [], SM]
		shell = "bedtools genomecov -ibam " + AlignPath + "/" + SMID + "." + chrom + ".bam -bga"
		sp = subprocess.Popen(shell, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		shellstdout, err = sp.communicate()
		if err:
			sp = subprocess.Popen(shell, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
			shellstdout, err = sp.communicate()
		if err:
			print("[E] " + SM + "\n" + err.decode("utf-8"))
			exit(1)
		shellstdout = shellstdout.decode("utf-8").split("\n")
		for line in shellstdout:
			if len(line) > 4:
				ele = line.split("\t")
				if len(ele) < 4:
					print(ele)
					exit(1)
				if ele[3] != "0":
					if len(ReturnData[1]) > 0 and ReturnData[1][len(ReturnData[1])-1][1] == ele[1]:
						ReturnData[1][len(ReturnData[1])-1][1] = ele[2]
					else:
						ReturnData[1].append([ele[1], ele[2]])
		del shellstdout
		s = ""
		for val in ReturnData[1]:
			s += "\t".join([ReturnData[0], val[0], val[1], str(int(val[1])-int(val[0]))]) + "\n"
		with open("ReadCovRegionList.tmp/" + SM + "." + chrom + ".txt", "w") as f:
			f.write(s)
		del s
		
		print(SM, chrom)
	except Exception as e:
		traceback.print_exc()

MAXTHR = len(CHRlis) if len(CHRlis) < MAXTHR else MAXTHR
with Pool(processes = MAXTHR) as proce_pool:
	mulres = []
	if os.path.exists("ReadCovRegionList/" + SMID + ".txt") and os.path.getsize("ReadCovRegionList/" + SMID + ".txt") > 50:
		print("[I] Skiping " + SMID + ".")
		exit()

	for chrom in CHRlis:
		mulres.append(proce_pool.apply_async(process_data, args=(SMID, chrom)))

	proce_pool.close()
	proce_pool.join()



s = "Chrom\tStart\tEnd\tLength\n"
for chrom in CHRlis:
	with open("ReadCovRegionList.tmp/" + SMID + "." + chrom + ".txt") as f:
		s += f.read().replace("\n\n", "\n")
with open("ReadCovRegionList/" + SMID + ".txt", "w") as f:
	f.write(s)
os.system("rm -f ReadCovRegionList.tmp/" + SMID + "*")
	

