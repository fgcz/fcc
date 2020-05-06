# Script for converting Thermo raw files to scpectrum files with Proteome Discoverer Daemon
# Prerequisite is a running Proteome Dicoverer process
# Simon Barkow-Oesterreicher, Functional Genomic Center Zurich
# 20140226 

from random import randint
import sys
import shutil
import os
import stat
import fileinput
import re

fileExtension=sys.argv[1]
pdParamFile=sys.argv[2]
inputFile=sys.argv[3]
outputFile=sys.argv[4]

tempPath = "E:\pd_temp\\"
randomName=str(randint(0,100000))
intermediateFile = tempPath + randomName + ".raw"
try:
    shutil.copy(inputFile,intermediateFile) 
except OSError as e:
    print("{} could not be copied into {}! Error {} occured".format(inputFile,tempPath,e))
    sys.exit()
            
myCommand = r'"D:\Program Files\Thermo\Discoverer Daemon 1.4\System\Release\DiscovererDaemon.exe" -p C:\FGCZ\fcc\proteomeDiscoverer\{} '.format(pdParamFile) + intermediateFile
try:
    os.system(myCommand)
    print("{} succeeded. {} has been converted.".format(myCommand, inputFile))
except OSError as e:
    print("{} failed!. {} has NOT been converted. Error {} occured".format(myCommand, inputFile, e))
    sys.exit()


def replaceTitle(spectrumFile):
    for line in fileinput.input(spectrumFile, inplace=True):
        line=line.strip()
        if (line.startswith("TITLE=")):
            print(line.replace(line, line + " [" + inputFile + "]"))
        else:
            print(line)

outFiles = []
#find OutputFile in tempDirectory
for root, dirs, files in os.walk(tempPath):
    #count number of mgf files exeeding 10KB
    for name in files:
        fullPath = os.path.join(root, name)
        if name.startswith(randomName) and name.endswith(fileExtension) and os.path.getsize(fullPath) > 10000:
            outFiles.append(name)

for fileName in outFiles:
    fullPath = os.path.join(root, fileName)
    replaceTitle(fullPath)
    outFileName=outputFile
    if (len(outFiles) > 1):
        if not os.path.isfile(outputFile):
            with open(outputFile, 'w') as outFile: print('dummy mgf file', file=outFile,  end='\r')
        #change output filename according to node
        if re.search("\[Node_10\]",fileName):
            outFileName = os.path.splitext(outputFile)[0] + '__ETDIT' + os.path.splitext(outputFile)[1]
        elif re.search("\[Node_27\]",fileName):
            outFileName = os.path.splitext(outputFile)[0] + '__HCDIT' + os.path.splitext(outputFile)[1]
        elif re.search("\[Node_31\]",fileName):
            outFileName = os.path.splitext(outputFile)[0] + '__ETDFT' + os.path.splitext(outputFile)[1]
        elif re.search("\[Node_14\]",fileName):
            outFileName = os.path.splitext(outputFile)[0] + '__HCDFT' + os.path.splitext(outputFile)[1]
        elif re.search("\[Node_34\]",fileName):
            outFileName = os.path.splitext(outputFile)[0] + '__CIDFT' + os.path.splitext(outputFile)[1]
        elif re.search("\[Node_04\]",fileName):
            outFileName = os.path.splitext(outputFile)[0] + '__CIDIT' + os.path.splitext(outputFile)[1]
    try:
        shutil.copy(fullPath,outFileName)
        print(fileName + " copied to " + outFileName)
    except OSError as e:
        print("{} could not be copied into {}! Error {} occured".format(fullPath,outFileName, e))
        pass

print("clean up:")
for root, dirs, files in os.walk(tempPath):
    for name in files:
        if (name.startswith(randomName)):
            fullPath = os.path.join(root, name)
            try:
                os.chmod(fullPath, stat.S_IWRITE)
            except OSError as e:
                print("chmod did not work, {}".format(e))
                pass
            try:
                os.remove(fullPath)
                print("{} removed".format(fullPath))
            except OSError as e:
                print("{} was not removed! Error {} occured".format(fullPath,e))
                pass
            
                



    






 
