#!/usr/bin/python

import sys
import os

#print("script name is", sys.argv[0])

fas = sys.argv[1]

#print(fas)

#fasDNA = open(fas, "r")

#fastaDNA = fasDNA.read()

#print(fastaDNA)

os.system("head {}".format(fas))

os.system("perl /scicomp/home/ydn3/exorcise_endline.pl < {} | perl /scicomp/home/ydn3/Sort_and_screen_contig_lengths.pl | head".format(fas))
