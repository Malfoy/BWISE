


# ***************************************************************************
#
#							   Bwise:
#				 High-order De Bruijn graph assembler
#
#
#
# ***************************************************************************

# ############################################################################
#									Utils functions
# ############################################################################


def printCommand(cmd,pc=True):
	if pc:
		print(cmd, flush=True)

# get the platform
def getPlatform():
	if sys.platform == "linux" or sys.platform == "linux2":
		return "linux"
	elif sys.platform == "darwin":
		return "OSX"
	else:
		print("[ERROR] BWISE is not compatible with Windows.")
		sys.exit(1);


# get the timestamp as string
def getTimestamp():
	return "[" + time.strftime("%H:%M:%S") + " " + time.strftime("%d/%m/%Y") + "] "



# check if reads files are present
def checkReadFiles(readfiles):
	if readfiles is None:
		return True
	allFilesAreOK = True
	#~ for file in readfiles:
	if not os.path.isfile(readfiles):
		print("[ERROR] File \""+file+"\" does not exist.")
		allFilesAreOK = False
	if not allFilesAreOK:
		dieToFatalError("One or several read files do not exist.")


# check if files written by BWISE are present
def checkWrittenFiles(files):
	allFilesAreOK = True
	if not os.path.isfile(files):
		print("[ERROR] There was a problem generating \"" + files + "\".")
		allFilesAreOK = False
	if not allFilesAreOK:
		dieToFatalError("One or several files could not be generated.")



# to return if an error makes the run impossible
def dieToFatalError (msg):
  print("[FATAL ERROR] " + msg)
  print("To find out why, try `Bwise --help` and/or check the logs files of the various steps of the pipeline (logs/logBloocoo, logs/logBcalm, logs/logTips, logs/logBgreat, logs/logK2000).")
  sys.exit(1);


# launch subprocess
def subprocessLauncher(cmd, argstdout=None, argstderr=None,	 argstdin=None):
	args = shlex.split(cmd)
	p = subprocess.Popen(args, stdin = argstdin, stdout = argstdout, stderr = argstderr).communicate()
	return p

def printTime(msg, seconds):
	m, s = divmod(seconds, 60)
	h, m = divmod(m, 60)
	return msg + " %d:%02d:%02d" % (h, m, s)


def printWarningMsg(msg):
	print("[Warning] " + msg)







# ############################################################################
#			   graph generation with BCALM + BTRIM + BGREAT
# ############################################################################

def graphConstruction(BWISE_MAIN, BWISE_INSTDIR, OUT_DIR, fileBcalm,k_min, k_max, kmer_solidity, Kmer_Coverage, SR_solidity, SR_Coverage, toolsArgs, fileCase, nb_cores, mappingEffort, missmatchAllowed,anchorSize, OUT_LOG_FILES,greedy_K2000,fastq_Bool):
	try:
		endLoop=False
		inputBcalm=fileBcalm
		print("\n" + getTimestamp() + "--> Starting Graph construction and Super Reads generation...")
		kmerList = ["0",str(k_min),"101","201","251","301","351","401","451","501","551","601","651","701","751","801"]
		os.chdir(OUT_DIR)
		logBcalm = "logs/logBcalm"
		logBcalmToWrite = open(logBcalm, 'w')
		logTips = "logs/logTips"
		logTipsToWrite = open(logTips, 'w')
		logBgreat = "logs/logBgreat"
		logBgreatToWrite = open(logBgreat, 'w')
		logK2000 = "logs/logK2000"
		logK2000ToWrite = open(logK2000, 'w')
		#~ os.chdir(BWISE_MAIN)
		#~ os.chdir(OUT_DIR)
		fastq_option=""
		if(fastq_Bool):
			fastq_option=" -q "

		indiceGraph = 1
		kmerSize = kmerList[indiceGraph]
		coreUsed = "20" if nb_cores == 0 else str(nb_cores)
		for indiceGraph in range(1, len(kmerList)):
			kmerSize = kmerList[indiceGraph]
			if(int(kmerList[indiceGraph]) >= k_max):
				kmerSize=str(k_max)
				greedy_K2000=1
				endLoop=True;
			if(os.path.isfile(OUT_DIR +"/dbg" + str(kmerList[indiceGraph])+".fa")):
				print("\t#Graph " + str(indiceGraph) + ": Already here ! Let us use it ", flush=True)
			else:
				print("#Graph " + str(indiceGraph) + ": Construction... ", flush=True)
				# BCALM
				cmd=BWISE_INSTDIR + "/bcalm -max-memory 10000 -in " + OUT_DIR + "/" + inputBcalm + " -kmer-size " + kmerSize + " -abundance-min " + str(kmer_solidity) + " -out " + OUT_DIR + "/out " + " -nb-cores " + coreUsed
				printCommand( "\t\t"+cmd)
				p = subprocessLauncher(cmd, logBcalmToWrite, logBcalmToWrite)
				checkWrittenFiles(OUT_DIR + "/out.unitigs.fa")

				#  Graph Cleaning
				print("\t #Graph cleaning... ", flush=True)
				# BTRIM
				if(kmer_solidity == 1):
					cmd=BWISE_INSTDIR + "/btrim -u out.unitigs.fa -k "+kmerSize+"  -t "+str((2*(int(kmerSize)-1)))+" -T 3  -c "+coreUsed+" -h 8 -o dbg"+str(kmerSize)+".fa -f 1"
				else:
					cmd=BWISE_INSTDIR + "/btrim -u out.unitigs.fa -k "+kmerSize+" -t "+str((2*(int(kmerSize)-1)))+" -T 3 -o dbg"+str(kmerSize)+".fa -c "+coreUsed+" -h 8 -f "+str(Kmer_Coverage)
				printCommand("\t\t"+cmd)
				p = subprocessLauncher(cmd, logTipsToWrite, logTipsToWrite)
				checkWrittenFiles(OUT_DIR + "/dbg"+str(kmerSize)+".fa")
				os.remove(OUT_DIR + "/out.unitigs.fa")
				os.remove(OUT_DIR + "/dbg"+str(kmerSize)+".fa1")
				os.remove(OUT_DIR + "/dbg"+str(kmerSize)+".fa2")
				for filename in glob.glob(OUT_DIR + "/out.*"):
					os.remove(filename)
				for filename in glob.glob(OUT_DIR + "/trashme*"):
					os.remove(filename)

			# Read Mapping
			if(not os.path.isfile(OUT_DIR +"/dbg" + str(kmerList[indiceGraph+1])+".fa")):

				print("\t#Read mapping with BGREAT... ", flush=True)
				# BGREAT
				cmd=BWISE_INSTDIR + "/bgreat   -k " + kmerSize + "  " + toolsArgs['bgreat'][fileCase] + " -g dbg" + str(kmerSize) + ".fa "+fastq_option+" -t " + coreUsed + "  -a "+str(anchorSize)+"   -m "+str(missmatchAllowed)+" -e "+str(mappingEffort)
				printCommand("\t\t"+cmd)
				p = subprocessLauncher(cmd, logBgreatToWrite, logBgreatToWrite)
				checkWrittenFiles(OUT_DIR + "/paths")

				print("\t#Super reads filtering... ", flush=True)
				#NUMBERFILTER
				#PREFILTER
				cmd=BWISE_INSTDIR + "/numbersFilter paths "+str(SR_Coverage)+" cleanedPaths_"+str(kmerSize)+" "+ coreUsed +" "+str(SR_solidity)
				printCommand("\t\t"+cmd)
				p = subprocessLauncher(cmd, logBgreatToWrite, logBgreatToWrite)
				#PHASE ONE
				cmd=BWISE_INSTDIR +"/compact.sh -i cleanedPaths_"+str(kmerSize)+" -u dbg" +	 str(kmerSize) + ".fa  -k "+kmerSize
				printCommand("\t\t"+cmd)
				p = subprocessLauncher(cmd, logBgreatToWrite, logBgreatToWrite)
				cmd=BWISE_INSTDIR + "/numbersFilter paths "+str(SR_Coverage)+" cleanedPaths_"+str(kmerSize)+" "+ coreUsed +" "+str(SR_solidity+1)
				printCommand("\t\t"+cmd)
				p = subprocessLauncher(cmd, logBgreatToWrite, logBgreatToWrite)
				cmd="cat compact > cleanedPaths_"+str(kmerSize)
				printCommand("\t\t"+cmd)
				#PHASE TWO
				p = subprocessLauncher(cmd, logBgreatToWrite, logBgreatToWrite)
				cmd=BWISE_INSTDIR +"/compact.sh -i cleanedPaths_"+str(kmerSize)+" -u dbg" +	 str(kmerSize) + ".fa  -k "+kmerSize
				printCommand("\t\t"+cmd)
				p = subprocessLauncher(cmd, logBgreatToWrite, logBgreatToWrite)
				p = subprocessLauncher(cmd, logBgreatToWrite, logBgreatToWrite)
				cmd=BWISE_INSTDIR + "/numbersFilter paths "+str(SR_Coverage)+" cleanedPaths_"+str(kmerSize)+" "+ coreUsed +" "+str(SR_solidity+2)
				printCommand("\t\t"+cmd)
				cmd="cat compact > cleanedPaths_"+str(kmerSize)
				printCommand("\t\t"+cmd)
				p = subprocessLauncher(cmd, logBgreatToWrite, logBgreatToWrite)
				#PHASE THREE
				p = subprocessLauncher(cmd, logBgreatToWrite, logBgreatToWrite)
				cmd=BWISE_INSTDIR +"/compact.sh -i cleanedPaths_"+str(kmerSize)+" -u dbg" +	 str(kmerSize) + ".fa  -k "+kmerSize
				printCommand("\t\t"+cmd)
				p = subprocessLauncher(cmd, logBgreatToWrite, logBgreatToWrite)
				p = subprocessLauncher(cmd, logBgreatToWrite, logBgreatToWrite)
				cmd=BWISE_INSTDIR + "/numbersFilter paths "+str(SR_Coverage)+" cleanedPaths_"+str(kmerSize)+" "+ coreUsed +" "+str(SR_solidity+3)
				printCommand("\t\t"+cmd)
				cmd="cat compact > cleanedPaths_"+str(kmerSize)
				printCommand("\t\t"+cmd)
				p = subprocessLauncher(cmd, logBgreatToWrite, logBgreatToWrite)

				print("\t#Contig generation... ", flush=True)
				#K2000
				if(greedy_K2000==0):
					cmd=BWISE_INSTDIR +"/run_K2000.sh -i cleanedPaths_"+str(kmerSize)+" -u dbg" +	 str(kmerSize) + ".fa  -k "+kmerSize+" -f  contigs_k"+kmerSize+".fa  -g  assemblyGraph_k"+kmerSize+".gfa"
				else:
					cmd=BWISE_INSTDIR +"/run_K2000.sh -i cleanedPaths_"+str(kmerSize)+" -u dbg" + str(kmerSize) + ".fa  -k "+kmerSize+" -f  contigs_k"+kmerSize+".fa  -g  assemblyGraph_k"+kmerSize+".gfa -t 1000 -c 100"
				#~ cmd=BWISE_INSTDIR +"/numbersToSequences dbg" + str(kmerList[indiceGraph]) + ".fa cleanedPaths_"+str(kmerList[indiceGraph])+" "+kmerSize+"  compacted_unitigs_k"+kmerSize+".fa"
				printCommand("\t\t"+cmd)
				p = subprocessLauncher(cmd, logK2000ToWrite, logK2000ToWrite)
				for filename in glob.glob(OUT_DIR + "/dbg_path_file_compacted*"):
					os.remove(filename)
				#~ os.remove(OUT_DIR + "/paths")
				if(endLoop):
					break;
			inputBcalm = "contigs_k"+kmerSize+".fa";
			kmer_solidity = 1

		#~ os.chdir(BWISE_MAIN)

		print(getTimestamp() + "--> Done!")
		os.remove(OUT_DIR + "/bankBcalm.txt")
		return {'indiceGraph': indiceGraph, 'kmerSize': kmerSize}
	except SystemExit:	# happens when checkWrittenFiles() returns an error
		sys.exit(1);
	except KeyboardInterrupt:
		sys.exit(1);
	except:
		print("Unexpected error during graph construction:", sys.exc_info()[0])
		dieToFatalError('')







# ############################################################################
#									Main
# ############################################################################
def main():

	wholeT = time.time()
	print("\n*** This is BWISE - High order De Bruijn graph assembler ***\n")
	#~ BWISE_MAIN = os.path.dirname(os.path.realpath(__file__))
	BWISE_INSTDIR =	 BWISE_MAIN + "/bin"  # todo : change using getPlatform()
	print("Binaries are in: " + BWISE_INSTDIR)

	# ========================================================================
	#						 Manage command line arguments
	# ========================================================================
	parser = argparse.ArgumentParser(description='Bwise - High order De Bruijn graph assembler ')

	# ------------------------------------------------------------------------
	#							 Define allowed options
	# ------------------------------------------------------------------------
	parser.add_argument("-x", action="store", dest="paired_readfiles",		type=str,					help="input fasta or (compressed .gz if -c option is != 0) paired-end read files. Several read files must be concatenated.")
	parser.add_argument("-u", action="store", dest="single_readfiles",		type=str,					help="input fasta or (compressed .gz if -c option is != 0) single-end read files. Several read files must be concatenated.")

	parser.add_argument('-s', action="store", dest="kmer_solidity",				type=int,	default = 2,	help="an integer, k-mers present strictly less than this number of times in the dataset will be discarded (default 2)")
	parser.add_argument('-S', action="store", dest="Kmer_Coverage",		type=int,	default = 5,	help="an integer, minimal unitig coverage for first cleaning (default 5)")

	parser.add_argument('-p', action="store", dest="SR_solidity",			type=int,	default = 2,	help="an integer,  super-reads present strictly less than this number of times will be discarded (default 2)")
	parser.add_argument('-P', action="store", dest="SR_Coverage",			type=int,	default = 5,	help="an integer  unitigs with less than S reads mapped is filtred (default 5)")

	parser.add_argument('-k', action="store", dest="k_min",					type=int,	default = 63,	help="an integer, smallest k-mer size (default 63)")
	parser.add_argument('-K', action="store", dest="k_max",					type=int,	default = 201,	help="an integer, largest k-mer size (default 201)")

	parser.add_argument('-e', action="store", dest="mapping_Effort",				type=int,	default = 1000,	help="Anchors to test for mapping (default 1000)")
	parser.add_argument('-a', action="store", dest="anchor_Size",				type=int,	default = 41,	help="Anchors size (default 41)")
	parser.add_argument('-m', action="store", dest="missmatch_allowed",				type=int,	default = 5,	help="missmatch allowed in mapping (default 5)")

	parser.add_argument('-g', action="store", dest="greedy_K2000",				type=int,	default = 0,	help="Greedy contig extension")

	parser.add_argument('-t', action="store", dest="nb_cores",				type=int,	default = 0,	help="number of cores used (default max)")
	parser.add_argument('-o', action="store", dest="out_dir",				type=str,	default=os.getcwd(),	help="path to store the results (default = current directory)")

	parser.add_argument('--version', action='version', version='%(prog)s 0.0.1')


	# ------------------------------------------------------------------------
	#				Parse and interpret command line arguments
	# ------------------------------------------------------------------------
	options = parser.parse_args()

	# ------------------------------------------------------------------------
	#				  Print command line
	# ------------------------------------------------------------------------
	print("The command line was: " + ' '.join(sys.argv))


	# ------------------------------------------------------------------------
	#				  Misc parameters
	# ------------------------------------------------------------------------
	k_min				= options.k_min
	k_max				= options.k_max
	kmer_solidity		= options.kmer_solidity
	Kmer_Coverage		= options.Kmer_Coverage
	SR_solidity			= options.SR_solidity
	nb_cores			= options.nb_cores
	mappingEffort		= options.mapping_Effort
	anchorSize		= options.anchor_Size
	SR_Coverage		= options.SR_Coverage
	missmatchAllowed	= options.missmatch_allowed
	greedy_K2000	= options.greedy_K2000


	# ------------------------------------------------------------------------
	#				Create output dir and log files
	# ------------------------------------------------------------------------
	OUT_DIR = options.out_dir
	try:
		if not os.path.exists(OUT_DIR):
			os.mkdir(OUT_DIR)
		else:
			printWarningMsg(OUT_DIR + " directory already exists, BWISE will use it.")

		OUT_LOG_FILES = OUT_DIR + "/logs"
		if not os.path.exists(OUT_LOG_FILES):
			os.mkdir(OUT_LOG_FILES)
		outName = OUT_DIR.split("/")[-1]
		OUT_DIR = os.path.dirname(os.path.realpath(OUT_DIR)) + "/" + outName


		print("Results will be stored in: ", OUT_DIR)
	except:
		print("Could not write in out directory :", sys.exc_info()[0])
		dieToFatalError('')

	# ------------------------------------------------------------------------
	#				  Parse input read options
	# ------------------------------------------------------------------------
	fastqFile=False
	try:
		bankBcalm = open(OUT_DIR + "/bankBcalm.txt", 'w');
	except:
		print("Could not write in out directory :", sys.exc_info()[0])

	# check if the given paired-end read files indeed exist
	paired_readfiles = None
	single_readfiles = None
	errorReadFile = 0
	if options.paired_readfiles:
		paired_readfiles =	''.join(options.paired_readfiles)
		try:
			paired_readfiles = os.path.abspath(paired_readfiles)
			checkReadFiles(options.paired_readfiles)
			if(paired_readfiles[-2:]=="fq"):
				fastqFile=True
		except:
			paired_readfiles = None
			errorReadFile = 1
	else:
		paired_readfiles = None
		errorReadFile = 1

	# check if the given single-end read files indeed exist
	if options.single_readfiles:
		single_readfiles = ''.join(options.single_readfiles)
		try:
			single_readfiles = os.path.abspath(single_readfiles)
			checkReadFiles(options.single_readfiles)
			if(single_readfiles[-2:]=="fq"):
				fastqFile=True
			errorReadFile *= 0
		except:
			single_readfiles = None
			errorReadFile *= 1
	else:
		single_readfiles = None
		errorReadFile *= 1

	if errorReadFile:
		parser.print_usage()
		dieToFatalError("BWISE requires at least a read file")


	parametersLog = open(OUT_DIR + "/ParametersUsed.txt", 'w');
	parametersLog.write("reads: "+str(paired_readfiles)+" "+ str(single_readfiles)+"	k_min: %s	k_max:%s	k-mer_solidity:%s	kmer_coverage:%s	SR_solidity:%s	SR_coverage:%s	mapping_effort:%s 	missmatch_allowed:%s	greedy_parameter:%s\n " %(k_min,k_max, kmer_solidity, Kmer_Coverage, SR_solidity, SR_Coverage, mappingEffort,missmatchAllowed,greedy_K2000 ))
	parametersLog.close()
	bloocooArg = ""
	bgreatArg = ""
	paired = '' if paired_readfiles is None else str(paired_readfiles)
	single = '' if single_readfiles is None else str(single_readfiles)
	both = paired + "," + single
	toolsArgs = {'bloocoo':{1: paired + " " , 2:  single + " " , 3: both + " "}, 'bgreat':{1:" -x "+str(paired_readfiles)+" ", 2: " -u "+str(single_readfiles)+" ", 3: " -x "+str(paired_readfiles)+"  -u "+str(single_readfiles)+" "}}




	if single_readfiles is not None and paired_readfiles is not None:  # paired end + single end
		fileCase = 3
		bankBcalm.write(str(paired_readfiles) +"\n"+str(single_readfiles))
	elif single_readfiles is None:	# paired end only
		fileCase = 1
		bankBcalm.write(str(paired_readfiles))
	else:  # single end only
		fileCase = 2
		bankBcalm.write(str(single_readfiles) )
	# bankBcalm.write(OUT_DIR + "lost_unitig.fa")
	bankBcalm.close()

	# ========================================================================
	#									RUN
	# ========================================================================


	# ------------------------------------------------------------------------
	#						   Graph construction and cleaning
	# ------------------------------------------------------------------------
	t = time.time()
	valuesGraph = graphConstruction(BWISE_MAIN, BWISE_INSTDIR, OUT_DIR, "bankBcalm.txt",k_min, k_max, kmer_solidity, Kmer_Coverage, SR_solidity, SR_Coverage,toolsArgs, fileCase, nb_cores, mappingEffort ,missmatchAllowed,anchorSize, OUT_LOG_FILES,greedy_K2000,fastqFile)
	print(printTime("Graph Construction took: ", time.time() - t))



	print(printTime("\nThe end !\nBWISE assembly took: ", time.time() - wholeT))



if __name__ == '__main__':
	main()
