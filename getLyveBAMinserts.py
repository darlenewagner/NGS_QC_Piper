#!/usr/bin/python

import sys
import os.path
import argparse
import re
import logging
import warnings
import subprocess
from os import listdir
from os.path import isfile, join

## Convert all .sorted.bam files in folder to .sam files
## Calculate median insert lengths for each .sam file
## RECOMMENDED for SMALT mapping of reads to their own draft assembly

## Requires samtools-0.1.x or higher
## also requires perl script: revis_gather_PE_read_stats.pl

def readable_dir(prospective_dir):
	if not os.path.isdir(prospective_dir):
    		raise argparse.ArgumentTypeError("readable_dir:{0} is not a valid path".format(prospective_dir))
	if os.access(prospective_dir, os.R_OK):
		if( not prospective_dir.endswith("/") ):
			prospective_dir = prospective_dir + "/"
		return prospective_dir
	else:
		raise argparse.ArgumentTypeError("readable_dir:{0} is not a readable dir".format(prospective_dir))


def getParentDir(filePathString):
	splitStr = re.split(pattern='/', string=filePathString)
	folderCount = len(splitStr) - 2
	isolatePath=splitStr[0]
	ii = 1
	while ii < folderCount:
		isolatePath = isolatePath + "/" + splitStr[ii]
		ii = ii + 1
	return isolatePath


logger = logging.getLogger("getLyveBAMinserts.py")
logger.setLevel(logging.INFO)

parser = argparse.ArgumentParser(description='compute insert sizes from .sorted.bam', usage="getLyveBAMinserts.py --dir Project/bam/ --verbose (Y/N)")

parser.add_argument('dir', type=readable_dir, action='store')

parser.add_argument("--verbose", default='N', type = lambda s : s.lower(), choices=['Y','N'])

args = parser.parse_args()

inFolder = args.dir

parentFolder = getParentDir(inFolder)

#print(inFolder + "\n" + parentFolder)

onlyFiles = [f for f in listdir(inFolder) if isfile(join(inFolder, f))]

samMedian = []
samAverage = []

samFolder = parentFolder + "/smSAM_BAM"

referencePath = parentFolder + "/reference/"

referenceFiles = [f for f in listdir(referencePath) if isfile(join(referencePath, f))]

referenceFasta = None

rr = 0

while rr < len(referenceFiles):
	if(re.search(r'\.fasta$', referenceFiles[rr], flags=re.IGNORECASE)):
		referenceFasta = referencePath + referenceFiles[rr]
		break
	rr = rr + 1

print(referenceFasta)

os.mkdir(samFolder)

ii = 0

while ii < len(onlyFiles):
	if((re.search(r'\.sorted\.bam$', onlyFiles[ii], flags=re.IGNORECASE)) and (re.match(r'C6472-M3235-14-007', onlyFiles[ii]) is None) and (re.match(r'CDC03-98-DC-M3235-15-014', onlyFiles[ii]) is None)):
		inputPath = inFolder + onlyFiles[ii]
		outName = re.sub(r'sorted\.bam', 'sam', onlyFiles[ii])
		outputSam = samFolder + "/" + outName
		os.system("samtools view {} > {}".format(inputPath, outputSam))
		outDepth = re.sub(r'\.sorted\.bam', '_depth.txt', onlyFiles[ii])
		outputDepth = samFolder + "/" + outDepth
		os.system("samtools depth {} > {}".format(inputPath, outputDepth))
		os.system("perl revis_gather_PE_read_stats.pl {} {} {}".format(referenceFasta, outputSam, outputDepth))
	ii = ii + 1

#os.rmdir(samFolder)

