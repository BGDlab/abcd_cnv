#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  4 16:38:04 2018

@author: maude
"""

import time, argparse, subprocess



parser = argparse.ArgumentParser(description="Pipeline for CNV calling and their annotation using genotyping data (final reports)")



parser.add_argument("outputDir", help="the directory for the results.")
parser.add_argument("cohortID", help="A cohort identifier. It will be inserted in the name of the outputs.")

args = parser.parse_args()

 
def formatForCNVision():
    start_time3 = time.time()
    print("format and filter files to CNVision requirements")
    
    formatScript = "Scripts/callCNVisionFormat_ukbb.sh"
    pennOption = "--PNformat"
    QSoption = "--QTformat"
    
    allPenn = baseDir + "/CNVisionManip/PC_allCNV.txt"
    allQS = baseDir + "/CNVisionManip/QS_allCNV.txt"
    
    
    proc1=subprocess.Popen([formatScript, allPenn, baseDir, tag, pennOption])
    proc2=subprocess.Popen([formatScript, allQS, baseDir, tag, QSoption])
    
    proc1.wait()
    proc2.wait()

    print("**** time format file for CNVision: %f" % (time.time() - start_time3)) 




def merge2algo():
    start_time4 = time.time()
    print("in CNVision merge 2 algo")
    mergeScript = "Scripts/callCNVisionMerge_ukbb.sh"
    
    PennFile = baseDir + "/CNVisionManip/PC_" + tag + "_CNVisionFormated.txt"
    QSfile = baseDir + "/CNVisionManip/QS_" + tag + "_CNVisionFormated.txt"
    
    
    proc1=subprocess.Popen([mergeScript, PennFile, QSfile, baseDir, tag])
    proc1.wait()

    print("**** time for CNVision merge 2 algo: %f" % (time.time() - start_time4)) 
    

    
             
start = time.time()

#global variables 

#PennDir = args.PC
#PennParentDir, resultDir = os.path.split(PennDir)
#QuantiDir = args.QS
#QuantiParentDir, resultDir = os.path.split(QuantiDir) 
baseDir = args.outputDir
tag = args.cohortID



#mergeRawResults()
formatForCNVision()
merge2algo()


print("the whole thing takes: %f" % (time.time() - start)) 
