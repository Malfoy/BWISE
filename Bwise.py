#!/usr/bin/env python3 
# -*- coding: utf-8 -*-

import os
import re
import sys
import time
import shlex, subprocess
import struct
import shutil
import os.path
import tempfile
import argparse
import threading
import multiprocessing
# from random import randint
from operator import itemgetter
from subprocess import Popen, PIPE, STDOUT

# ***************************************************************************
#
#							   Bwise:
#				 High order De Bruijn graph assembler
#	 
#	   
#
# ***************************************************************************

# ############################################################################
#									Utils functions
# ############################################################################

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
		dieToFatalError("One or more read files do not exist.")

		
# check if files written by BWISE are present
def checkWrittenFiles(files):
	allFilesAreOK = True
	if not os.path.isfile(files):
		print("[ERROR] There was a problem writing \"" + files + "\".")
		allFilesAreOK = False
	if not allFilesAreOK:
		dieToFatalError("One or more files could not be written.")



# to return if an error makes the run impossible
def dieToFatalError (msg):
  print("[FATAL ERROR] " + msg)
  print("Try `Bwise --help` for more information")
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
#									Correction fonction using Bloocoo
# ############################################################################

def correctionReads(BWISE_MAIN, BWISE_INSTDIR, paired_readfiles, single_readfiles, toolsArgs, fileCase, nb_correction_steps, OUT_DIR, nb_cores, OUT_LOG_FILES):
	try:
		print("\n" + getTimestamp() + "--> Starting Read Correction with Bloocoo...")
		slowParameter = " -slow "
		kmerSizeCorrection = ["31", "63", "95", "127"]
		bloocooversion = ["32", "64", "128", "128"]
		os.chdir(OUT_LOG_FILES)
		logBloocoo = "logBloocoo"
		logBloocooToWrite = open(logBloocoo, 'w')
		os.chdir(BWISE_MAIN)
		os.chdir(OUT_DIR)
		indiceCorrection = 0
		for indiceCorrection in range(min(nb_correction_steps, len(kmerSizeCorrection))):
			# print("		Correction step " + str(indiceCorrection + 1) + "... ", end='', flush=True)
			logHistoCorr = "histocorr" + str(kmerSizeCorrection[indiceCorrection])
			logHistoCorrToWrite = open(logHistoCorr, 'w')
			# Bloocoo
			cmd=BWISE_INSTDIR + "/Bloocoo" + bloocooversion[indiceCorrection] + " -file " + toolsArgs['bloocoo'][fileCase] + slowParameter + "-kmer-size " + kmerSizeCorrection[indiceCorrection] + " -nbits-bloom 24  -out reads_corrected" + str(indiceCorrection + 1) + ".fa -nb-cores " + nb_cores
			print("\tCorrection step " + str(indiceCorrection + 1), flush=True)
			print( "\t\t"+cmd)
			p = subprocessLauncher(cmd, logBloocooToWrite, logBloocooToWrite)
			# Deal with files after Bloocoo
			
			#TODO=put back the histogram creation
			# cmd=BWISE_INSTDIR + "/h5dump -y -d histogram_"+kmerSizeCorrection[indiceCorrection]+" reads_corrected" + str(indiceCorrection + 1) + ".fa.h5"
			# print("\t\t"+cmd)
			# p = subprocessLauncher(cmd, logHistoCorrToWrite, logHistoCorrToWrite)
			checkWrittenFiles(OUT_DIR + "/reads_corrected" + str(indiceCorrection + 1) + ".fa.h5")
			# if (indiceCorrection > 0):
			#	  cmd="rm -f " + OUT_DIR + "/reads_corrected" + str(indiceCorrection) + "* "
			#	  print("\t\t\t"+cmd)
			#	  p = subprocessLauncher(cmd, None, logHistoCorrToWrite)
			if fileCase == 3:
				cmd="mv reads_corrected" + str(indiceCorrection + 1) + "_0_.fasta reads_corrected" + str(indiceCorrection + 1) + "1.fa "
				print("\t\t\t"+cmd)
				p = subprocessLauncher(cmd)
				p = subprocessLauncher("mv reads_corrected" + str(indiceCorrection + 1) + "_1_.fasta reads_corrected" + str(indiceCorrection + 1) + "2.fa ")
				checkWrittenFiles(OUT_DIR + "/reads_corrected" + str(indiceCorrection + 1) + "1.fa")
				checkWrittenFiles(OUT_DIR + "/reads_corrected" + str(indiceCorrection + 1) + "2.fa")
				toolsArgs['bloocoo'][fileCase] = OUT_DIR + "/reads_corrected" + str(indiceCorrection + 1) + "1.fa," + OUT_DIR + "/reads_corrected" + str(indiceCorrection + 1) + "2.fa"
			else:
				checkWrittenFiles(OUT_DIR + "/reads_corrected" + str(indiceCorrection + 1) + ".fa")
				toolsArgs['bloocoo'][fileCase] = OUT_DIR + "/reads_corrected" + str(indiceCorrection + 1) + ".fa "
			logHistoCorrToWrite.close()
			checkWrittenFiles(OUT_DIR + "/histocorr" + str(kmerSizeCorrection[indiceCorrection]))
			
		os.chdir(BWISE_MAIN)
		# links and file check
		if nb_correction_steps == 0:
			if fileCase == 3:
				# no correction : linking raw read files to reads_corrected1.fa and reads_corrected2.fa
				cmd="ln -s " + paired_readfiles + " " + OUT_DIR + "/reads_corrected1.fa"
				print("\t\t\t"+cmd)
				p = subprocessLauncher(cmd, None, subprocess.DEVNULL)
				cmd="ln -s " + single_readfiles + " " + OUT_DIR + "/reads_corrected2.fa"
				print("\t\t\t"+cmd)
				p = subprocessLauncher(cmd, None, subprocess.DEVNULL)
				checkWrittenFiles(OUT_DIR + "/reads_corrected1.fa")
				checkWrittenFiles(OUT_DIR + "/reads_corrected2.fa")
			else:
				cmd="ln -s " + toolsArgs['bloocoo'][fileCase] + " " + OUT_DIR + "/reads_corrected.fa"
				print("\t\t\t"+cmd)
				p = subprocessLauncher(cmd, None, subprocess.DEVNULL)
				checkWrittenFiles(OUT_DIR + "/reads_corrected.fa")
		else:
			if fileCase == 3:
				# linking last corrected reads files to reads_corrected1.fa and reads_corrected2.fa
				cmd="ln -s " + OUT_DIR + "/reads_corrected" + str(indiceCorrection + 1) + "1.fa " + OUT_DIR + "/reads_corrected1.fa"
				print("\t\t\t"+cmd)
				p = subprocessLauncher(cmd, None, subprocess.DEVNULL) 
				cmd="ln -s " + OUT_DIR + "/reads_corrected" + str(indiceCorrection + 1) + "2.fa " + OUT_DIR + "/reads_corrected2.fa"
				print("\t\t\t"+cmd)
				p = subprocessLauncher(cmd, None, subprocess.DEVNULL)
				checkWrittenFiles(OUT_DIR + "/reads_corrected1.fa")
				checkWrittenFiles(OUT_DIR + "/reads_corrected2.fa")
			else:
				cmd="ln -s " + toolsArgs['bloocoo'][fileCase] + " " + OUT_DIR + "/reads_corrected.fa"
				print("\t\t\t"+cmd)
				p = subprocessLauncher(cmd)
				checkWrittenFiles(OUT_DIR + "/reads_corrected.fa")
		
		print("\n" + getTimestamp() + "--> Done!")
	except SystemExit:	# happens when checkWrittenFiles() returns an error
		sys.exit(1);
	except KeyboardInterrupt:
		sys.exit(1);
	except:
		print("Unexpected error during read correction:", sys.exc_info()[0])
		dieToFatalError('')




# ############################################################################
#			   graph generation with BCALM + kMILL + BGREAT
# ############################################################################

def graphConstruction(BWISE_MAIN, BWISE_INSTDIR, OUT_DIR, fileBcalm, k_max, solidity, unitigFilter, superReadsCleaning, toolsArgs, fileCase, nb_cores, OUT_LOG_FILES, cleanedPathsFile, newPathsFile):
	try:
		print("\n" + getTimestamp() + "--> Starting Graph construction and Super Reads generation...")
		kmerList = ["0","51","101","151","201","251","301","351","401","451","501"]	 # todo: better kmer list
		os.chdir(OUT_LOG_FILES)
		logBcalm = "logBcalm"
		logBcalmToWrite = open(logBcalm, 'w')
		logTips = "logTips"
		logTipsToWrite = open(logTips, 'w')
		logBgreat = "logBgreat"
		logBgreatToWrite = open(logBgreat, 'w')
		logK2000 = "logK2000"
		logK2000ToWrite = open(logK2000, 'w')
		os.chdir(BWISE_MAIN)
		os.chdir(OUT_DIR)
		indiceGraph = 1
		kmerSize = kmerList[indiceGraph]
		coreUsed = "20" if nb_cores == "0" else nb_cores
		for indiceGraph in range(1, len(kmerList)):
			if int(kmerList[indiceGraph]) <= k_max:
				kmerSize = kmerList[indiceGraph]
				kmerSizeTip = 2 * int(kmerSize)
				print("\tGraph " + str(indiceGraph) + ": Construction... ", flush=True)
				#  Graph construction, Bcalm
				cmd=BWISE_INSTDIR + "/bcalm -in " + OUT_DIR + "/" + fileBcalm + " -kmer-size " + kmerSize + " -abundance-min " + str(solidity) + " -out " + OUT_DIR + "/out " + " -nb-cores " + nb_cores
				print( "\t\t"+cmd)
				p = subprocessLauncher(cmd, logBcalmToWrite, logBcalmToWrite)
				checkWrittenFiles(OUT_DIR + "/out.unitigs.fa")
				#  Graph Cleaning
				print("\t\t Cleaning... ", flush=True)
				
				# kMILL + tip cleaning
				cmd=BWISE_INSTDIR + "/kMILL out.unitigs.fa " + str(int(kmerSize) - 1) + " " + str(int(kmerSize) - 2)
				print("\t\t\t"+cmd)
				p = subprocessLauncher(cmd, logBcalmToWrite, logBcalmToWrite)
				checkWrittenFiles(OUT_DIR + "/out_out.unitigs.fa.fa")
				cmd=BWISE_INSTDIR + "/tipCleaner out_out.unitigs.fa.fa "	 + str(int(kmerSize) - 1) + " " + str(kmerSizeTip)
				print("\t\t\t"+cmd)
				p = subprocessLauncher(cmd, logTipsToWrite, logTipsToWrite)
				checkWrittenFiles(OUT_DIR + "/tiped.fa")
				cmd=BWISE_INSTDIR + "/kMILL tiped.fa " + str(int(kmerSize) - 1) + " " + str(int(kmerSize) - 2)
				print("\t\t\t"+cmd)
				p = subprocessLauncher(cmd, logBcalmToWrite, logBcalmToWrite)
				checkWrittenFiles(OUT_DIR + "/out_tiped.fa.fa")
				cmd="mv out_tiped.fa.fa dbg" + str(kmerList[indiceGraph]) + ".fa"
				print("\t\t\t"+cmd)
				p = subprocessLauncher(cmd)
				checkWrittenFiles(OUT_DIR + "/dbg" + str(kmerList[indiceGraph]) + ".fa")
				# Read Mapping
				print("\tRead mapping with BGREAT... ")
				# BGREAT
				cmd=BWISE_INSTDIR + "/bgreat -k " + kmerSize + " -M -i 10 " + toolsArgs['bgreat'][fileCase] + " -g dbg" + str(kmerList[indiceGraph]) + ".fa -t " + coreUsed + " -a 63 -m 0 -e 100"
				print("\t\t"+cmd)
				p = subprocessLauncher(cmd, logBgreatToWrite, logBgreatToWrite)
				checkWrittenFiles(OUT_DIR + "/paths")
				if (indiceGraph == 1):
					cmd=BWISE_INSTDIR + "/numbersFilter paths " + str(unitigFilter) + " " + str(superReadsCleaning) + " dbg" + str(kmerList[indiceGraph]) + ".fa "	+ kmerSize
					print("\t\t"+cmd)
					p = subprocessLauncher(cmd, cleanedPathsFile, logBgreatToWrite)
					cmd=BWISE_INSTDIR + "/numbersToSequences dbg" + str(kmerList[indiceGraph]) + ".fa cleanedPaths " + str(int(kmerSize) - 1)
					print("\t\t"+cmd)
					p = subprocessLauncher(cmd, newPathsFile, logBgreatToWrite)
				else:
					if int(kmerList[indiceGraph]) < k_max:
						compactedUnitigs = open("compacted_unitigs" + str(indiceGraph) + ".txt", 'w')
						cmd=BWISE_INSTDIR + "/numbersFilter paths " + str(unitigFilter) + " " + str(superReadsCleaning) + " dbg" + str(kmerList[indiceGraph]) + ".fa " + kmerSize
						print("\t\t"+cmd)
						p = subprocessLauncher(cmd, cleanedPathsFile, logBgreatToWrite)
						# cmd="python3 " + BWISE_INSTDIR + "/K2000.py cleanedPaths dbg" + str(kmerList[indiceGraph]) + ".fa " + kmerSize + " -e"
						
						cmd=BWISE_INSTDIR +"/run_K2000.sh cleanedPaths dbg" + str(kmerList[indiceGraph]) + ".fa "+kmerSize+" compacted_unitigs_k"+kmerSize+".gfa compacted_unitigs_k"+kmerSize+".fa"
						print("\t\t"+cmd)
						p = subprocessLauncher(cmd, logK2000ToWrite, logK2000ToWrite)
						system.exit(1)
				
						# checkWrittenFiles(OUT_DIR + "/" + "compacted_unitigs" + str(indiceGraph) + ".txt")
						# cmd=BWISE_INSTDIR + "/numbersToSequences dbg" + str(kmerList[indiceGraph]) + ".fa compacted_unitigs" + str(indiceGraph) + ".txt " + str(int(kmerSize) - 1)
						# print("\t\t"+cmd)
						# p = subprocessLauncher(cmd, newPathsFile, logBgreatToWrite)
						cmd="ln -s " + toolsArgs['bloocoo'][fileCase] + " " + OUT_DIR + "/reads_corrected.fa"
						compactedUnitigs.close()
					else:
						cmd=BWISE_INSTDIR + "/numbersFilter paths " + str(unitigFilter) + " " + str(superReadsCleaning) + " dbg" + str(kmerList[indiceGraph]) + ".fa " + kmerSize
						print("\t\t"+cmd)
						p = subprocessLauncher(cmd, cleanedPathsFile, logBgreatToWrite)
			else:
				break
			fileBcalm = "newPaths";
			solidity = 1
		os.chdir(BWISE_MAIN)
		 
		print(getTimestamp() + "--> Done!")
		return {'indiceGraph': indiceGraph, 'kmerSize': kmerSize}
	except SystemExit:	# happens when checkWrittenFiles() returns an error
		sys.exit(1);
	except KeyboardInterrupt:
		sys.exit(1);
	except:
		print("Unexpected error during graph construction:", sys.exc_info()[0])
		dieToFatalError('')

	

# ############################################################################
#			   Super Reads compaction with k2000
# ############################################################################

def srCompaction(BWISE_MAIN, BWISE_INSTDIR, OUT_DIR, valuesGraph, OUT_LOG_FILES):
	try:
		print("\n" + getTimestamp() + "--> Starting Super Reads Compaction...")
		# K2000
		os.chdir(OUT_LOG_FILES)
		logkmill = "logkmill"
		logSRC = "logSRC"
		logK2000ToWrite = open("logK2000", 'a')
		logSRCToWrite = open(logSRC, 'a')
		logkmillToWrite = open(logkmill, 'a')
		os.chdir(OUT_DIR)
		#~ print(BWISE_INSTDIR + "/run_K2000.sh cleanedPaths dbg" + str(valuesGraph['indiceGraph'] - 1) + ".fa " + valuesGraph['kmerSize'] + " contigs.gfa contigs.fa")
		cmd=BWISE_INSTDIR + "/run_K2000.sh cleanedPaths dbg" + str(valuesGraph['indiceGraph'] - 1) + ".fa " + valuesGraph['kmerSize'] + " contigs.gfa contigs.fa"
		print("\t\t"+cmd)
		p = subprocessLauncher(cmd, None, logK2000ToWrite)
		checkWrittenFiles(OUT_DIR + "/contigs.gfa")
		checkWrittenFiles(OUT_DIR + "/contigs.fa")
		checkWrittenFiles(OUT_DIR + "/cleanedPaths_compacted")
		cmd="mv cleanedPaths_compacted compacted_unitigs.txt"
		print("\t\t\t"+cmd)
		p = subprocessLauncher(cmd)
		checkWrittenFiles(OUT_DIR + "/compacted_unitigs.txt")
		cmd="rm -rf trashme* *.h5 out.unitigs.fa notAligned.fa bankBready bankBcalm.txt maximalSuperReads.fa newPaths out_out.unitigs.fa.fa tiped.fa paths"
		print("\t\t\t"+cmd)
		p = subprocessLauncher(cmd, logkmillToWrite, logkmillToWrite)
		print(getTimestamp() + "--> Done!")
	except SystemExit:	# happens when checkWrittenFiles() returns an error
		sys.exit(1);
	except KeyboardInterrupt:
		sys.exit(1);
	except:
		print("Unexpected error during super reads compaction:", sys.exc_info()[0])
		dieToFatalError('')


		

# ############################################################################
#									Main
# ############################################################################
def main():

	wholeT = time.time()
	print("\n*** This is BWISE - High order De Bruijn graph assembler ***\n")
	BWISE_MAIN = os.path.dirname(os.path.realpath(__file__))
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
	parser.add_argument('-s', action="store", dest="min_cov",				type=int,	default = 2,	help="an integer, k-mers present strictly less than this number of times in the dataset will be discarded (default 2)")
	parser.add_argument('-S', action="store", dest="min_cov_uni",						default = 2,	help="an integer, unitigs present strictly less than this number of times in the dataset will be discarded (default 2)")
	parser.add_argument('-o', action="store", dest="out_dir",				type=str,	default=".",	help="path to store the results (default = .)")
	parser.add_argument('-k', action="store", dest="k_max",					type=int,	default = 301,	help="an integer, largest k-mer size (default=301)")
	parser.add_argument('-p', action="store", dest="min_cov_SR",			type=int,	default = 2,	help="an integer,  super-reads present strictly less than this number of times will be discarded(default 2)")
	parser.add_argument('-c', action="store", dest="nb_correction_steps",	type=int,	default = 4,	help="an integer, number of steps of read correction (default max=4)")
	parser.add_argument('-t', action="store", dest="nb_cores",				type=str,	default = "0",	help="number of cores used (default max)")
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
	k_max				= options.k_max
	min_cov				= options.min_cov
	min_cov_uni			= options.min_cov_uni
	min_cov_SR			= options.min_cov_SR
	nb_correction_steps = options.nb_correction_steps
	nb_cores			= options.nb_cores
	
	if nb_correction_steps > 4:
		dieToFatalError("Please use value <= 4 for correction steps.") 

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
		parametersLog = open(OUT_DIR + "/ParametersUsed.txt", 'w');
		parametersLog.write("k_max:%s k-mer_solidity:%s unitig_solidity:%s SR_cleaning:%s correction_steps:%s" %(k_max, min_cov, min_cov_uni, min_cov_SR, nb_correction_steps))
		print("Results will be stored in: ", OUT_DIR)
	except:
		print("Could not write in out directory :", sys.exc_info()[0])
		dieToFatalError('')
		
	# ------------------------------------------------------------------------
	#				  Parse input read options
	# ------------------------------------------------------------------------
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

	bloocooArg = ""
	bgreatArg = ""
	paired = '' if paired_readfiles is None else str(paired_readfiles) 
	single = '' if single_readfiles is None else str(single_readfiles) 
	both = paired + "," + single
	toolsArgs = {'bloocoo':{1: paired + " " , 2:  single + " " , 3: both + " "}, 'bgreat':{1:" -x reads_corrected.fa ", 2: " -u reads_corrected.fa ", 3: " -x reads_corrected1.fa  -u reads_corrected2.fa "}}

	try:
		cleanedPathsFile =	open(OUT_DIR + "/cleanedPaths", 'w')
		newPathsFile	 =	open(OUT_DIR + "/newPaths", 'w')
		checkWrittenFiles(OUT_DIR + "/cleanedPaths")
		checkWrittenFiles(OUT_DIR + "/newPaths")
	except:
		print("Could not write in out directory :", sys.exc_info()[0])
		dieToFatalError('')
	
	
	if single_readfiles is not None and paired_readfiles is not None:  # paired end + single end
		fileCase = 3
		bankBcalm.write(OUT_DIR + "/reads_corrected1.fa\n" + OUT_DIR + "/reads_corrected2.fa\n")
	elif single_readfiles is None:	# paired end only
		fileCase = 1
		bankBcalm.write(OUT_DIR + "/reads_corrected.fa\n")
	else:  # single end only
		fileCase = 2
		bankBcalm.write(OUT_DIR + "/reads_corrected.fa\n")
	bankBcalm.close()

	# ========================================================================
	#									RUN
	# ========================================================================

	
	# ------------------------------------------------------------------------
	#						   Correction
	# ------------------------------------------------------------------------
	t = time.time()
	correctionReads(BWISE_MAIN, BWISE_INSTDIR, paired_readfiles, single_readfiles, toolsArgs, fileCase, nb_correction_steps, OUT_DIR, nb_cores, OUT_LOG_FILES)
	print(printTime("Correction took: ", time.time() - t))


	# ------------------------------------------------------------------------
	#						   Graph construction and cleaning
	# ------------------------------------------------------------------------
	t = time.time()
	valuesGraph = graphConstruction(BWISE_MAIN, BWISE_INSTDIR, OUT_DIR, "bankBcalm.txt", k_max, min_cov, min_cov_uni, min_cov_SR, toolsArgs, fileCase, nb_cores, OUT_LOG_FILES, cleanedPathsFile, newPathsFile)
	print(printTime("Graph Construction took: ", time.time() - t))
	cleanedPathsFile.close()
	newPathsFile.close()

	# ------------------------------------------------------------------------
	#						   Super Reads Compaction
	# ------------------------------------------------------------------------
	t = time.time()
	srCompaction(BWISE_MAIN, BWISE_INSTDIR, OUT_DIR, valuesGraph, OUT_LOG_FILES)
	print(printTime("Super Reads Compaction took:", time.time() - t))

	
	print(printTime("\nThe end !\nBWISE assembly took: ", time.time() - wholeT))
	


if __name__ == '__main__':
	main()
