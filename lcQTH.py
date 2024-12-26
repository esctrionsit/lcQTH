from concurrent.futures import ThreadPoolExecutor, wait, FIRST_COMPLETED, ALL_COMPLETED
import concurrent.futures
import traceback
import hashlib
import random
import math
import json
import sys
import os
scriptPath = os.path.abspath(__file__)
projectPath = os.path.split(scriptPath)[0]

def printHelp():
	print("")
	print("Usage: python3 lcQTH.py [options] COMMAND")
	print("")
	print("Toolkit for tracing parental haplotypes and perform simple QTL mapping with aligned ultra-low-coverage sequencing data.")
	print("")
	print("Options:")
	print("  --help            Show this message and exit.")
	print("")
	print("Commands:")
	print("  init              Initialize this lcQTH project, cleaning previous results.")
	# print("  status            Print status and analysis progress of this lcQTH project.")
	# print("  loadconf          Load a configuration file into this project.")
	# print("  updateconf        Update configurations with options directly.")
	print("  run               Run the analysis program continually (with no options) or from a specific step.")
	print("  exporthmp         Export parental haplotypes of RILs in hapmap format.")
	print("  exportxlsx        Export parental haplotypes of RILs in hapmap format and saved in an Excel spreadsheet")
	print("")
	print("Run 'python3 lcQTH.py COMMAND --help' for more information on a command.")
	print("lcQTH online help: <https://github.com/esctrionsit/lcQTH/wiki/>")

def printCommandHelp(command):
	if command == "init":
		print("")
		print("Usage: python3 lcQTH.py init [OPTIONS]")
		print("")
		print("Initialize this lcQTH project, cleaning previous results.")
		print("")
		print("Options:")
		print("  --config=[path]   Initialize the project with a configuration file.")
		print("  --ParentalResolutionInKb=[int]   (Optional) Manually setting the window size in parental haplotype defining in step 1. Its default value will be 2*ResolutionInKb.")
		print("")
		print("lcQTH online help: <https://github.com/esctrionsit/lcQTH/wiki/>")
	elif command == "status" and False:
		print("")
		print("Usage: python3 lcQTH.py status [OPTIONS]")
		print("")
		print("Print status and analysis progress of this lcQTH project.")
		print("")
		print("Options:")
		print("  --html=[path]     Save status and progress to an HTML file.")
		print("")
		print("lcQTH online help: <https://github.com/esctrionsit/lcQTH/wiki/>")
	elif command == "loadconf" and False:
		print("")
		print("Usage: python3 lcQTH.py loadconf [OPTIONS]")
		print("")
		print("Load a configuration file into this project.")
		print("")
		print("Options:")
		print("  --config=[path]   Provide the path of config file.")
		print("  --html=[path]     Save status and progress to an HTML file.")
		print("")
		print("For details of config file, please visit <https://>")
		print("lcQTH online help: <https://github.com/esctrionsit/lcQTH/wiki/>")
	elif command == "updateconf" and False:
		print("")
		print("Usage: python3 lcQTH.py updateconf [OPTIONS]")
		print("")
		print("Update configurations with options directly.")
		print("")
		print("Options:")
		print("  --resolution=[int]  Change the window size in haplotype tracing. Window size should be in kilo-basepairs.")
		print("")
		print("lcQTH online help: <https://github.com/esctrionsit/lcQTH/wiki/>")
	elif command == "run":
		print("")
		print("Usage: python3 lcQTH.py run [OPTIONS]")
		print("")
		print("Run the analysis program continually (with no options) or from a specific")
		print("step.")
		print("")
		print("Options:")
		print("  --from=stepID     Start analysis from a specific step.")
		print("  --step=stepID     Run analysis of a specific step only.")
		# print("  --overwrite       Overwrite previous results in running.")
		print("")
		print("Available step IDs:")
		print("  1   Identifying parental haplotypes.")
		print("  2   Identificate CNV blocks in parents and offspring lines.")
		print("  3   Identificate reads-covered regions in offspring lines.")
		print("  4   Assign raw parental haplotypes to RCRs.")
		print("  5   Trace parental haplotypes in offspring lines.")
		print("  6   Remove markers in regions that parents have the identical backgrounts and do filtering.")
		print("  7   (optional) Plot parental haplotypes in offspring lines.")
		print("  8   Build genetic map.")
		print("  9   QTL mapping.")
		print("")
		print("lcQTH online help: <https://github.com/esctrionsit/lcQTH/wiki/>")
	elif command == "exporthmp":
		print("")
		print("Usage: python3 lcQTH.py exporthmp [OPTIONS]")
		print("")
		print("Export parental haplotypes of RILs in hapmap format.")
		print("")
		print("Options:")
		print("  --out=[path]      Provide a path of the file.")
		print("")
		print("lcQTH online help: <https://github.com/esctrionsit/lcQTH/wiki/>")
	elif command == "exportxlsx":
		print("")
		print("Usage: python3 lcQTH.py exportxlsx [OPTIONS]")
		print("")
		print("Export parental haplotypes of RILs in hapmap format and saved in an Excel spreadsheet.")
		print("")
		print("Options:")
		print("  --out=[path]      Provide a path of the file.")
		print("")
		print("lcQTH online help: <https://github.com/esctrionsit/lcQTH/wiki/>")
	else:
		print("")
		print("Unknown command: " + command)
		print("For helps, please run: python3 lcQTH.py --help ")
		print("")
		print("lcQTH online help: <https://github.com/esctrionsit/lcQTH/wiki/>")

def printHelpInfo(scene):
	if scene == "runAnalysis":
		print("[I] LcQTH will continue to use the parameter settings analyzed last time. To update parameters, please run \"python3 lcQTH.py loadconf\".")

def random_str(num):
    H = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789'
    salt = ''
    for i in range(num):
        salt += random.choice(H)
    return salt

def init(Para):
	# clean files
	cleanFiles(projectPath)
	# clean files
	newRuntime(projectPath)
	# Update config
	if "config" in Para:
		conf = loadConfigFile(Para["config"], "" if (not "ParentalResolutionInKb" in Para) else Para["ParentalResolutionInKb"])
		saveConfig(conf)

def loadConfigFile(filePath, ParentalResolutionInKb):
	# Read file
	if filePath == None:
		filePath = "config.json"
	with open(filePath) as f:
		s = f.read()
	# Parse JSON
	try:
		configs = json.loads(s)
	except Exception as e:
		print("[E] Error occured during parsing config file. The details are as follow:", file=sys.stderr)
		traceback.print_exc()
		exit(1)
	# Validate
	missingList = []
	for configTerm in ["ChromosomeLengthFile", "ChromosomeCentromereInMBFile", "BamFolder", "VcfFolder", "VcfFileSuffixes", "ParentsName", "OffspringListFile", "ResolutionInKb", "MaxThreads", "TmpFolderPath", "PhenotypeFile"]:
		if not configTerm in configs:
			missingList.append(configTerm)
	if len(missingList) != 0:
		print("[E] Key terms are missing in config file. The missing terms are: " + ", ".join(missingList), file=sys.stderr)
		exit(1)
	if not "transmat_" in configs:
		configs["transmat_"] = [
			[0.9, 0.08, 0.02],
			[0.15, 0.7, 0.15],
			[0.02, 0.08, 0.9],
	    ]
	if not "emissionprob_" in configs:
		configs["emissionprob_"] = [
	        [0.9, 0.04, 0.04, 4e-05, 0.02],
	        [0.05, 0.1, 0.7, 0.1, 0.05],
	        [0.02, 4e-05, 0.04, 0.04, 0.9]
	    ]
	# Parse ParentalResolutionInKb
	if ParentalResolutionInKb != "":
		configs["ParentalResolutionInKb"] = int(ParentalResolutionInKb)
	else:
		configs["ParentalResolutionInKb"] = 2*configs["ResolutionInKb"]
	# Convert relative path to absolute path
	if len(configs["ChromosomeLengthFile"]) != 0:
		if configs["ChromosomeLengthFile"][0] != "/":
			configs["ChromosomeLengthFile"] = os.path.normpath(projectPath + "/" + configs["ChromosomeLengthFile"])
			print("[W] Value of term \"ChromosomeLengthFile\" in config file is a relative path. It has been converted to \"" + configs["ChromosomeLengthFile"] + "\"", file=sys.stderr)
	else:
		print("[E] Value of term \"ChromosomeLengthFile\" in config file is empty.", file=sys.stderr)
		exit(1)
	if len(configs["ChromosomeCentromereInMBFile"]) != 0:
		if configs["ChromosomeCentromereInMBFile"][0] != "/":
			configs["ChromosomeCentromereInMBFile"] = os.path.normpath(projectPath + "/" + configs["ChromosomeCentromereInMBFile"])
			print("[W] Value of term \"ChromosomeCentromereInMBFile\" in config file is a relative path. It has been converted to \"" + configs["ChromosomeCentromereInMBFile"] + "\"", file=sys.stderr)
	else:
		print("[E] Value of term \"ChromosomeCentromereInMBFile\" in config file is empty.", file=sys.stderr)
		exit(1)
	if len(configs["BamFolder"]) != 0:
		if configs["BamFolder"][0] != "/":
			configs["BamFolder"] = os.path.normpath(projectPath + "/" + configs["BamFolder"])
			print("[W] Value of term \"BamFolder\" in config file is a relative path. It has been converted to \"" + os.path.normpath(projectPath + "/" + configs["BamFolder"]) + "\"", file=sys.stderr)
	else:
		print("[E] Value of term \"BamFolder\" in config file is empty.", file=sys.stderr)
		exit(1)
	if len(configs["VcfFolder"]) != 0:
		if configs["VcfFolder"][0] != "/":
			configs["VcfFolder"] = os.path.normpath(projectPath + "/" + configs["VcfFolder"])
			print("[W] Value of term \"VcfFolder\" in config file is a relative path. It has been converted to \"" + os.path.normpath(projectPath + "/" + configs["VcfFolder"]) + "\"", file=sys.stderr)
	else:
		print("[E] Value of term \"VcfFolder\" in config file is empty.", file=sys.stderr)
		exit(1)
	if len(configs["OffspringListFile"]) != 0:
		if configs["OffspringListFile"][0] != "/":
			configs["OffspringListFile"] = os.path.normpath(projectPath + "/" + configs["OffspringListFile"])
			print("[W] Value of term \"OffspringListFile\" in config file is a relative path. It has been converted to \"" + os.path.normpath(projectPath + "/" + configs["OffspringListFile"]) + "\"", file=sys.stderr)
	else:
		print("[E] Value of term \"OffspringListFile\" in config file is empty.", file=sys.stderr)
		exit(1)
	if len(configs["TmpFolderPath"]) != 0:
		if configs["TmpFolderPath"][0] != "/":
			configs["TmpFolderPath"] = os.path.normpath(projectPath + "/" + configs["TmpFolderPath"])
			print("[W] Value of term \"TmpFolderPath\" in config file is a relative path. It has been converted to \"" + os.path.normpath(projectPath + "/" + configs["TmpFolderPath"]) + "\"", file=sys.stderr)
	else:
		print("[E] Value of term \"TmpFolderPath\" in config file is empty.", file=sys.stderr)
		exit(1)
	if len(configs["PhenotypeFile"]) != 0:
		if configs["PhenotypeFile"][0] != "/":
			configs["PhenotypeFile"] = os.path.normpath(projectPath + "/" + configs["PhenotypeFile"])
			print("[W] Value of term \"PhenotypeFile\" in config file is a relative path. It has been converted to \"" + os.path.normpath(projectPath + "/" + configs["PhenotypeFile"]) + "\"", file=sys.stderr)
	else:
		print("[E] Value of term \"PhenotypeFile\" in config file is empty.", file=sys.stderr)
		exit(1)
	# Return
	return configs

def saveConfig(conf):
	# Load chrom length
	if not os.path.exists(conf["ChromosomeLengthFile"]):
		print("[E] The path for chromosome length file is invalid.", file=sys.stderr)
		exit(1)
	with open(conf["ChromosomeLengthFile"]) as f:
		lines = f.readlines()
	ChrLen = {}
	for line in lines:
		if line.strip() == "":
			continue
		else:
			ele = line.strip().split()
			ChrLen[ele[0]] = int(ele[1])
	conf["ChrLen"] = ChrLen
	# Load chrom cent pos
	if not os.path.exists(conf["ChromosomeCentromereInMBFile"]):
		print("[E] The path for chromosome centromere position file is invalid.", file=sys.stderr)
		exit(1)
	with open(conf["ChromosomeCentromereInMBFile"]) as f:
		lines = f.readlines()
	ChrCent = {}
	for line in lines:
		if line.strip() == "":
			continue
		else:
			ele = line.strip().split()
			ChrCent[ele[0]] = int(ele[1])
	conf["ChrCent"] = ChrCent
	# Save
	with open(projectPath + "/.paras", "w") as f:
		f.write(json.dumps(conf))

def loadConfig():
	if (not os.path.exists(projectPath + "/.paras")) or (not os.path.exists(projectPath + "/runtime")):
		print("[E] The project does not appear to have been initialized. Please run \"python3 lcQTH.py init --help\" for helps of initialization.", file=sys.stderr)
		exit(1)
	else:
		confText = ""
		with open(projectPath + "/.paras") as f:
			confText = f.read()
		try:
			configs = json.loads(confText)
		except Exception as e:
			print("[E] Error occured during loading project settings.", file=sys.stderr)
			traceback.print_exc()
			exit(1)
		return configs	

def cleanFiles(projectPath):
	# Execuate cleaning
	if os.path.exists(projectPath + "/runtime"):
		os.system("rm -r " + projectPath + "/runtime")

def newRuntime(projectPath):
	if not os.path.exists(projectPath + "/src"):
		print("[E] Failed in initializing the project due to incomplete code. Please re-clone the code from GitHub at https://", file=sys.stderr)
		exit(1)
	else:
		os.system("cp -r " + projectPath + "/src " + projectPath + "/runtime")
		os.system("g++ " + projectPath + "/runtime/01.ParentsHaplotype/src/dsrdist.cpp -o " + projectPath + "/runtime/01.ParentsHaplotype/src/dsrdist")
		# os.system("g++ " + projectPath + "/runtime/01.ParentsHaplotype/src/dsrdist.cpp -o " + projectPath + "/runtime/01.ParentsHaplotype/src/dsrdist")

def checkEnv():
	ErrFlag = False
	try:
		import hmmlearn
		if not hmmlearn.__version__ in ("0.2.4", "0.2.5", "0.2.6"):
			print("[W] Current version (" + hmmlearn.__version__ + ") of hmmlearn is not tested. Version 0.2.4, 0.2.5, and 0.2.6 are known to be supported.", file=sys.stderr)
	except:
		print("[E] hmmlearn is not installed in current Python environment.", file=sys.stderr)
	
def loadStatus():
	pass

def saveStatusToHtml(status):
	pass

def runStep1(conf):
	# Identifying parental haplotypes
	## Generate bed files
	os.system("mkdir -p " + projectPath + "/runtime/01.ParentsHaplotype/bed/")
	os.system("mkdir -p " + projectPath + "/runtime/01.ParentsHaplotype/CNV/")
	os.system("mkdir -p " + projectPath + "/runtime/01.ParentsHaplotype/diff/")
	for chrom in conf["ChrLen"]:
		s = ""
		for i in range(0, conf["ChrLen"][chrom], conf["ParentalResolutionInKb"]*1000):
			s += chrom + "\t" + str(i+1) + "\t" + str(i+conf["ParentalResolutionInKb"]*1000) + "\n"
		with open(projectPath + "/runtime/01.ParentsHaplotype/bed/" + chrom + ".bed", "w") as f:
			f.write(s)
	## Generate running script
	token = random_str(5)
	s = "cd " + projectPath + "/runtime/01.ParentsHaplotype/\n"
	s += "python3 01.CNV.py \"" + conf["BamFolder"] + "\" \"" + conf["TmpFolderPath"] + "\" \"" + ",".join(list(conf["ChrLen"].keys())) + "\" \"" + ",".join(conf["ParentsName"]) + "\" \"" + str(conf["MaxThreads"]) + "\"\n"
	s += "python3 02.DSR.py \"" + conf["VcfFolder"] + "\" \"" + conf["VcfFileSuffixes"] + "\" \"" + conf["TmpFolderPath"] + "\" \"" + ",".join(conf["ParentsName"]) + "\" \"" + ",".join(list(conf["ChrLen"].keys())) + "\" \"" + str(conf["MaxThreads"]) + "\"\n"
	s += "python3 03.HMM.py \"" + ",".join(conf["ParentsName"]) + "\" \"" + ",".join(list(conf["ChrLen"].keys())) + "\" \"" + str(conf["MaxThreads"]) + "\" \"" + str(conf["ParentalResolutionInKb"]) + "\"\n"
	with open(conf["TmpFolderPath"] + "/" + token, "w") as f:
		f.write(s)
	os.system("bash " + conf["TmpFolderPath"] + "/" + token)

def runStep2(conf):
	# Identificate CNV blocks in parents and offspring lines
	## Generate bed files
	os.system("mkdir -p " + projectPath + "/runtime/02.CNVidentify/bed/")
	os.system("mkdir -p " + projectPath + "/runtime/02.CNVidentify/CNV/")
	for chrom in conf["ChrLen"]:
		s = ""
		for i in range(0, conf["ChrLen"][chrom], conf["ParentalResolutionInKb"]*1000):
			s += chrom + "\t" + str(i+1) + "\t" + str(i+conf["ParentalResolutionInKb"]*1000) + "\n"
		with open(projectPath + "/runtime/02.CNVidentify/bed/" + chrom + ".bed", "w") as f:
			f.write(s)
	## Load offspring list
	if not os.path.exists(conf["OffspringListFile"]):
		print("[E] Offspring list file does not exist. It was supposed to be at: " + conf["OffspringListFile"], file=sys.stderr)
		exit(1)
	with open(conf["OffspringListFile"]) as f:
		lines = f.readlines()
	offspring = []
	for line in lines:
		if line.strip() != "":
			offspring.append(line.strip())
	## Generate running script
	token = random_str(5)
	s = "cd " + projectPath + "/runtime/02.CNVidentify/\n"
	s += "python3 01.CNV.py \"" + conf["BamFolder"] + "\" \"" + conf["TmpFolderPath"] + "\" \"" + ",".join(list(conf["ChrLen"].keys())) + "\" \"" + ",".join(conf["ParentsName"]+offspring) + "\" \"" + str(conf["MaxThreads"]) + "\"\n"
	with open(conf["TmpFolderPath"] + "/" + token, "w") as f:
		f.write(s)
	os.system("bash " + conf["TmpFolderPath"] + "/" + token)

def runStep3(conf):
	# Identificate reads-covered regions in offspring lines
	## Load offspring list
	os.system("mkdir -p " + projectPath + "/runtime/03.IdentifyRCR/ReadCovRegionList.tmp/")
	os.system("mkdir -p " + projectPath + "/runtime/03.IdentifyRCR/ReadCovRegionList/")
	if not os.path.exists(conf["OffspringListFile"]):
		print("[E] Offspring list file does not exist. It was supposed to be at: " + conf["OffspringListFile"], file=sys.stderr)
		exit(1)
	with open(conf["OffspringListFile"]) as f:
		lines = f.readlines()
	offspring = []
	for line in lines:
		if line.strip() != "":
			offspring.append(line.strip())
	# Run
	MAXTHR4Call = int(conf["MaxThreads"]/len(conf["ChrLen"].keys()))
	MAXTHR4Call = 1 if MAXTHR4Call < 1 else MAXTHR4Call
	with ThreadPoolExecutor(max_workers=16) as t: 
		all_task = []
		for SM in offspring:
			token = random_str(5)
			s = "cd " + projectPath + "/runtime/03.IdentifyRCR/\n"
			s += "python3 generate.py \"" + SM + "\" \"" + conf["BamFolder"] + "\" \"" + conf["TmpFolderPath"] + "\" \""  + str(conf["MaxThreads"]) + "\" \"" + ",".join(list(conf["ChrLen"].keys())) + "\"\n"
			s += "rm " + conf["TmpFolderPath"] + "/" + token
			with open(conf["TmpFolderPath"] + "/" + token, "w") as f:
				f.write(s)
			all_task.append(t.submit(runShell, "bash " + conf["TmpFolderPath"] + "/" + token))
		wait(all_task, return_when=ALL_COMPLETED)

def runStep4(conf):
	# Assign raw parental haplotypes to RCRs.
	## Load offspring list
	os.system("mkdir -p " + projectPath + "/runtime/04.HaplotypeRCR/Rawphase/")
	if not os.path.exists(conf["OffspringListFile"]):
		print("[E] Offspring list file does not exist. It was supposed to be at: " + conf["OffspringListFile"], file=sys.stderr)
		exit(1)
	with open(conf["OffspringListFile"]) as f:
		lines = f.readlines()
	offspring = []
	for line in lines:
		if line.strip() != "":
			offspring.append(line.strip())
	## Run
	token = random_str(5)
	retries = 0
	maxRetries = 2
	while retries <= maxRetries:
		filenames = os.listdir(projectPath + "/runtime/04.HaplotypeRCR/Rawphase/")
		SMlst = {}
		for filename in filenames:
			SM = filename.split(".")[0]
			SMlst[SM] = None
		s = ""
		count = 0
		for SM in offspring:
			if not SM in SMlst:
				s += SM + "\n"
				count += 1
		if s != "":
			if retries > 0:
				print("[W] " + str(count) + " samples are failed in step 4. Retrying..." + str(retries) + "/" + str(maxRetries), file=sys.stderr)
			# Write offspring list
			with open(projectPath + "/runtime/04.HaplotypeRCR/Childs.txt", "w") as f:
				f.write(s)
			s = "cd " + projectPath + "/runtime/04.HaplotypeRCR/\n"
			s += "go run phaseRaw.bychr.go \"" + projectPath + "/runtime/04.HaplotypeRCR/Childs.txt" + "\" \"" + conf["ParentsName"][0] + "\" \"" + conf["ParentsName"][1] + "\" \""  + str(conf["MaxThreads"]) + "\" \"" + conf["VcfFolder"] + "\" \"" + conf["VcfFileSuffixes"] + "\" \"" + conf["TmpFolderPath"] + "\" \""  + str(len(conf["ChrLen"])) + "\"\n"
			s += "rm " + conf["TmpFolderPath"] + "/" + token + "\n"
			with open(conf["TmpFolderPath"] + "/" + token, "w") as f:
				f.write(s)
			os.system("bash " + conf["TmpFolderPath"] + "/" + token)
		else:
			break
		retries += 1
	# Check result integrality
	if retries >= maxRetries:
		filenames = os.listdir(projectPath + "/runtime/04.HaplotypeRCR/Rawphase/")
		SMlst = {}
		for filename in filenames:
			SM = filename.split(".")[0]
			SMlst[SM] = None
		count = 0
		for SM in offspring:
			if not SM in SMlst:
				count += 1
		if count != 0:
			print("[E] " + str(count) + " samples still failed in step 4. Maximum retry times have been reached.", file=sys.stderr)
			exit(1)

def runStep5(conf):	
	# Trace parental haplotypes in offspring lines
	## Load offspring list
	os.system("mkdir -p " + projectPath + "/runtime/05.HMMParsing/data/")
	os.system("mkdir -p " + projectPath + "/runtime/05.HMMParsing/RawPhaseParentalRatio/")
	os.system("mkdir -p " + projectPath + "/runtime/05.HMMParsing/SMphase/")
	os.system("mkdir -p " + projectPath + "/runtime/05.HMMParsing/SMphase.smooth/")
	if not os.path.exists(conf["OffspringListFile"]):
		print("[E] Offspring list file does not exist. It was supposed to be at: " + conf["OffspringListFile"], file=sys.stderr)
		exit(1)
	with open(conf["OffspringListFile"]) as f:
		lines = f.readlines()
	offspring = []
	for line in lines:
		if line.strip() != "":
			offspring.append(line.strip())
	## Write offspring list
	with open(projectPath + "/runtime/05.HMMParsing/Childs.txt", "w") as f:
		f.write("\n".join(offspring) + "\n")
	## Write chrom length
	with open(projectPath + "/runtime/05.HMMParsing/chrlen.txt", "w") as f:
		for chrom in conf["ChrLen"]:
			f.write(chrom + "\t" + str(conf["ChrLen"][chrom]) + "\n")
	# Run
	token = random_str(5)
	tp = {
		"transmat_": conf["transmat_"],
		"emissionprob_": conf["emissionprob_"],
	}
	with open(conf["TmpFolderPath"] + "/" + token + ".para", "w") as f:
		f.write(json.dumps(tp))
	s = "cd " + projectPath + "/runtime/05.HMMParsing/\n"
	s += "rm -f data/*\n"
	s += "rm -f RawPhaseParentalRatio/*\n"
	s += "go run 01.rawphase.go \"Childs.txt\" \"chrlen.txt\" \"" + str(conf["MaxThreads"]) + "\" \"" + str(conf["ResolutionInKb"]*1000) + "\" \"" + str(conf["ResolutionInKb"]*1000*2) + "\"\n"
	s += "python3 02.Chr2SM.py  \"" + ",".join(list(conf["ChrLen"].keys())) + "\"\n"
	s += "python3 04.Smooth.py  \"" + str(conf["MaxThreads"]) + "\" \"" + str(conf["ResolutionInKb"]*1000) + "\" \"" + str(conf["ParentalResolutionInKb"]*1000) + "\" \"" + ",".join(conf["ParentsName"]) + "\" " + conf["TmpFolderPath"] + "/" + token + ".para" + "\n"
	s += "python3 05.finalRawMarkers.py  \"" + str(conf["ResolutionInKb"]) + "\" \"" + ",".join(list(conf["ChrLen"].keys())) + "\"\n"
	s += "rm " + conf["TmpFolderPath"] + "/" + token + "*\n"
	with open(conf["TmpFolderPath"] + "/" + token, "w") as f:
		f.write(s)
	os.system("bash " + conf["TmpFolderPath"] + "/" + token)

def runStep6(conf):
	# Remove markers in regions that parents have the identical backgrounts and do filtering
	## Load offspring list
	os.system("mkdir -p " + projectPath + "/runtime/06.MarkerScaling/ParentalRatio/")
	os.system("mkdir -p " + projectPath + "/runtime/06.MarkerScaling/binned/")
	if not os.path.exists(conf["OffspringListFile"]):
		print("[E] Offspring list file does not exist. It was supposed to be at: " + conf["OffspringListFile"], file=sys.stderr)
		exit(1)
	with open(conf["OffspringListFile"]) as f:
		lines = f.readlines()
	offspring = []
	for line in lines:
		if line.strip() != "":
			offspring.append(line.strip())
	## Write offspring list
	with open(projectPath + "/runtime/06.MarkerScaling/Childs.txt", "w") as f:
		f.write("\n".join(offspring) + "\n")
	## Write chrom length
	with open(projectPath + "/runtime/06.MarkerScaling/chrlen.txt", "w") as f:
		for chrom in conf["ChrLen"]:
			f.write(chrom + "\t" + str(conf["ChrLen"][chrom]) + "\n")
	# Run
	token = random_str(5)
	s = "cd " + projectPath + "/runtime/06.MarkerScaling/\n"
	s += "python3 downScaling.parallel.py \"" + str(conf["ParentalResolutionInKb"]) + "\" \"" + str(conf["ParentalResolutionInKb"]) + "\" \"" + str(conf["ResolutionInKb"]) + "\" \"" + str(conf["ResolutionInKb"]) + "\" \"" + str(conf["MaxThreads"]) + "\"\n"
	s += "rm " + conf["TmpFolderPath"] + "/" + token + "*\n"
	with open(conf["TmpFolderPath"] + "/" + token, "w") as f:
		f.write(s)
	os.system("bash " + conf["TmpFolderPath"] + "/" + token)

def runStep7(conf):
	# Plot parental haplotypes in offspring lines
	## Load offspring list
	os.system("mkdir -p " + projectPath + "/runtime/07.ChildInheritionPlot/InheritplotBySample/Inheritplot/")
	if not os.path.exists(conf["OffspringListFile"]):
		print("[E] Offspring list file does not exist. It was supposed to be at: " + conf["OffspringListFile"], file=sys.stderr)
		exit(1)
	with open(conf["OffspringListFile"]) as f:
		lines = f.readlines()
	offspring = []
	for line in lines:
		if line.strip() != "":
			offspring.append(line.strip())
	## Write offspring list
	with open(projectPath + "/runtime/07.ChildInheritionPlot/Childs.txt", "w") as f:
		f.write("\n".join(offspring) + "\n")
	## Write chrom length
	with open(projectPath + "/runtime/07.ChildInheritionPlot/chrlen.txt", "w") as f:
		for chrom in conf["ChrLen"]:
			f.write(chrom + "\t" + str(math.ceil(conf["ChrLen"][chrom]/1000)) + "\n")
	## Write chrom cent pos
	with open(projectPath + "/runtime/07.ChildInheritionPlot/chrcent.txt", "w") as f:
		for chrom in conf["ChrCent"]:
			f.write(chrom + "\t" + str(conf["ChrCent"][chrom]) + "\n")
	# Run
	token = random_str(5)
	s = "cd " + projectPath + "/runtime/07.ChildInheritionPlot/InheritplotBySample\n"
	s += "python3 01.CallPlot.py \"" + str(conf["ResolutionInKb"]) + "\" \"" + conf["TmpFolderPath"] + "\" \"" + str(conf["MaxThreads"]) + "\"\n"
	with open(conf["TmpFolderPath"] + "/" + token, "w") as f:
		f.write(s)
	os.system("bash " + conf["TmpFolderPath"] + "/" + token)

def runStep8(conf):
	# Build genetic map
	## Load offspring list
	os.system("mkdir -p " + projectPath + "/runtime/09.BuildGeneticMap/GM/")
	if not os.path.exists(conf["OffspringListFile"]):
		print("[E] Offspring list file does not exist. It was supposed to be at: " + conf["OffspringListFile"], file=sys.stderr)
		exit(1)
	with open(conf["OffspringListFile"]) as f:
		lines = f.readlines()
	offspring = []
	for line in lines:
		if line.strip() != "":
			offspring.append(line.strip())
	## Write offspring list
	with open(projectPath + "/runtime/09.BuildGeneticMap/Childs.txt", "w") as f:
		f.write("\n".join(offspring) + "\n")
	## Write chrom length
	with open(projectPath + "/runtime/09.BuildGeneticMap/chrlen.txt", "w") as f:
		for chrom in conf["ChrLen"]:
			f.write(chrom + "\t" + str(conf["ChrLen"][chrom]) + "\n")
	# Run
	token = random_str(5)
	s = "cd " + projectPath + "/runtime/09.BuildGeneticMap/\n"
	s += "python3 01.generate.Rqtl.csv.py \"" + ",".join(list(conf["ChrLen"].keys())) + "\"\n"
	s += "python3 02.parallel.py \"" + str(conf["MaxThreads"]) + "\"\n"
	s += "Rscript 03.merge.R \"" + str(len(conf["ChrLen"])) + "\"\n"
	s += "Rscript 04.plotGeneticMap.R\n"
	with open(conf["TmpFolderPath"] + "/" + token, "w") as f:
		f.write(s)
	os.system("bash " + conf["TmpFolderPath"] + "/" + token)

def runStep9(conf):
	# QTL mapping
	token = random_str(5)
	s = "cd " + projectPath + "/runtime/10.QTLMapping/\n"
	s += "python3 01.Rqtl.data.py \"" + conf["PhenotypeFile"] + "\"\n"
	s += "Rscript 02.plot.R\n"
	with open(conf["TmpFolderPath"] + "/" + token, "w") as f:
		f.write(s)
	os.system("bash " + conf["TmpFolderPath"] + "/" + token)

def runFrom(conf, stepID):
	if stepID <= 1:
		runStep1(conf)
	if stepID <= 2:
		runStep2(conf)
	if stepID <= 3:
		runStep3(conf)
	if stepID <= 4:
		runStep4(conf)
	if stepID <= 5:
		runStep5(conf)
	if stepID <= 6:
		runStep6(conf)
	if stepID <= 7:
		runStep7(conf)
	if stepID <= 8:
		runStep8(conf)
	if stepID <= 9:
		runStep9(conf)

def saveHmp(conf, savepath):
	# relative to absolute
	if savepath != "/":
		savepath = os.path.normpath(projectPath + "/" + savepath)
	print("[I] The hmp file would be saved at: " + savepath, file=sys.stderr)
	if os.path.exists(savepath):
		if not os.path.isfile(savepath):
			print("[E] Save path (\"" + savepath + "\") for hmp file is a folder.", file=sys.stderr)
			exit(1)
		else:
			print("[W] Save path (\"" + savepath + "\") for hmp file is exists. The previous file would be re-writed.", file=sys.stderr)
	# Load offspring list
	if not os.path.exists(conf["OffspringListFile"]):
		print("[E] Offspring list file does not exist. It was supposed to be at: " + conf["OffspringListFile"], file=sys.stderr)
		exit(1)
	with open(conf["OffspringListFile"]) as f:
		lines = f.readlines()
	offspring = []
	for line in lines:
		if line.strip() != "":
			offspring.append(line.strip())
	# Write offspring list
	with open(projectPath + "/runtime/08.Marker2Table/Childs.txt", "w") as f:
		f.write("\n".join(offspring) + "\n")
	# run
	token = random_str(5)
	s = "cd " + projectPath + "/runtime/08.Marker2Table/\n"
	s += "python3 02.generate.hmp.py \"" + ",".join(list(conf["ChrLen"].keys())) + "\"\n"
	s += "cp 02.data.merged.hmp.txt " + savepath + "\n"
	with open(conf["TmpFolderPath"] + "/" + token, "w") as f:
		f.write(s)
	os.system("bash " + conf["TmpFolderPath"] + "/" + token)

def saveHmpXlsx(conf, savepath):
	# relative to absolute
	if savepath != "/":
		savepath = os.path.normpath(projectPath + "/" + savepath)
	print("[I] The haplotype xlsx file would be saved at: " + savepath, file=sys.stderr)
	if os.path.exists(savepath):
		if not os.path.isfile(savepath):
			print("[E] Save path (\"" + savepath + "\") for hmp file is a folder.", file=sys.stderr)
			exit(1)
		else:
			print("[W] Save path (\"" + savepath + "\") for hmp file is exists. The previous file would be re-writed.", file=sys.stderr)
	# Load offspring list
	if not os.path.exists(conf["OffspringListFile"]):
		print("[E] Offspring list file does not exist. It was supposed to be at: " + conf["OffspringListFile"], file=sys.stderr)
		exit(1)
	with open(conf["OffspringListFile"]) as f:
		lines = f.readlines()
	offspring = []
	for line in lines:
		if line.strip() != "":
			offspring.append(line.strip())
	# Write offspring list
	with open(projectPath + "/runtime/08.Marker2Table/Childs.txt", "w") as f:
		f.write("\n".join(offspring) + "\n")
	# run
	token = random_str(5)
	s = "cd " + projectPath + "/runtime/08.Marker2Table/\n"
	s += "python3 02.generate.hmp.py \"" + ",".join(list(conf["ChrLen"].keys())) + "\"\n"
	s += "python3 03.hmp2excel.py \"" + ",".join(conf["ParentsName"]) + "\"\n"
	s += "cp 03.hmp2excel.xlsx " + savepath + "\n"
	with open(conf["TmpFolderPath"] + "/" + token, "w") as f:
		f.write(s)
	os.system("bash " + conf["TmpFolderPath"] + "/" + token)

def readCommand():
	l = sys.argv
	command = {
		"Command": l[1],
		"Parameters": {}
	}
	for i in range(2, len(l)):
		if len(l[i]) < 3 or l[i][:2] != "--":
			print("[W] Unreconized parameter: " + l[i], file=sys.stderr)
		else:
			if "=" in l[i]:
				ele = l[i].split("=", 1)
				command["Parameters"][ele[0]] = ele[1]
			else:
				command["Parameters"][l[i]] = True
	return command

def runShell(command):
	os.system(command)

if __name__ == '__main__':
	if len(sys.argv) == 1:
		printHelp()
	else:
		command = readCommand()
		if command["Command"] == "init":
			# Check parameters
			for key in command["Parameters"]:
				if not key in ("--help", "--config", "--ParentalResolutionInKb"):
					print("[W] Unused parameter: " + key, file=sys.stderr)
			# Branchs
			if "--help" in command["Parameters"]:
				printCommandHelp(command["Command"])
			else:
				if "--config" in command["Parameters"]:
					initPara = {
						"config": command["Parameters"]["--config"]
					}
					if "--ParentalResolutionInKb" in command["Parameters"]:
						initPara["ParentalResolutionInKb"] = command["Parameters"]["--ParentalResolutionInKb"]
					init(initPara)
				else:
					printCommandHelp(command["Command"])
		elif command["Command"] == "run":
			# Check parameters
			for key in command["Parameters"]:
				if not key in ("--help", "--from", "--step"):
					print("[W] Unused parameter: " + key, file=sys.stderr)
			# Branchs
			if "--help" in command["Parameters"]:
				printCommandHelp(command["Command"])
			else:
				conf = loadConfig()
				if "--from" in command["Parameters"]:
					try:
						stepID = int(command["Parameters"]["--from"])
					except:
						print("[E] Unknown step ID: " + command["Parameters"]["--from"], file=sys.stderr)
						exit(1)
					if stepID < 1 or stepID > 9:
						print("[E] Unknown step ID: " + command["Parameters"]["--from"], file=sys.stderr)
						exit(1)
					runFrom(conf, stepID)
				elif "--step" in command["Parameters"]:
					if command["Parameters"]["--step"] == "1":
						runStep1(conf)
					elif command["Parameters"]["--step"] == "2":
						runStep2(conf)
					elif command["Parameters"]["--step"] == "3":
						runStep3(conf)
					elif command["Parameters"]["--step"] == "4":
						runStep4(conf)
					elif command["Parameters"]["--step"] == "5":
						runStep5(conf)
					elif command["Parameters"]["--step"] == "6":
						runStep6(conf)
					elif command["Parameters"]["--step"] == "7":
						runStep7(conf)
					elif command["Parameters"]["--step"] == "8":
						runStep8(conf)
					elif command["Parameters"]["--step"] == "9":
						runStep9(conf)
					else:
						printCommandHelp(command["Command"])
				else:
					# from start
					runFrom(conf, 1)
		elif command["Command"] == "exporthmp":
			# Check parameters
			for key in command["Parameters"]:
				if not key in ("--help", "--out"):
					print("[W] Unused parameter: " + key, file=sys.stderr)
			if "--help" in command["Parameters"]:
				printCommandHelp(command["Command"])
			elif "--out" in command["Parameters"]:
				conf = loadConfig()
				saveHmp(conf, command["Parameters"]["--out"])
			else:
				print("[E] Please provide a path to save the hmp file with parameter \"--out=[path]\"", file=sys.stderr)
				exit(1)
		elif command["Command"] == "exportxlsx":
			# Check parameters
			for key in command["Parameters"]:
				if not key in ("--help", "--out"):
					print("[W] Unused parameter: " + key, file=sys.stderr)
			if "--help" in command["Parameters"]:
				printCommandHelp(command["Command"])
			elif "--out" in command["Parameters"]:
				conf = loadConfig()
				saveHmpXlsx(conf, command["Parameters"]["--out"])
			else:
				print("[E] Please provide a path to save the xlsx file with parameter \"--out=[path]\"", file=sys.stderr)
				exit(1)
		elif command["Command"] == "--help":
			printHelp()
		else:
			printCommandHelp(command["Command"])
