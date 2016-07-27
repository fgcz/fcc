# python script for converting Thermo raw files to scpectrum files with Proteome Discoverer Daemon
# Prerequisite is a running Proteome Dicoverer process
# Simon Barkow-Oesterreicher, Functional Genomic Center Zurich, 20140226 
#   
# $HeadURL: http://fgcz-svn.uzh.ch/repos/fgcz/stable/proteomics/fcc/fgcz_pd2distiller_wrapper.py $
# $Id: fgcz_pd2distiller_wrapper.py 7385 2015-04-16 11:46:29Z cpanse $
# $Date: 2015-04-16 13:46:29 +0200 (Thu, 16 Apr 2015) $
# $Author: cpanse $


from random import randint
import sys
import shutil
import os
import stat
import fileinput
import re
import signal
import tempfile


"""
Goes through all the lines and adapt the TITLE field working for FGCZ 
BFabric mapping.
Input and Output files are the same.
"""
def replace_title(spectrumFile, para):


    pattern="(^TITLE=File:).(\"[-_:\.\\a-zA-Z]{3,}\")(;\ Spe.+)$"
    p = re.compile(pattern)
    patternFile = ".*_(ETD|EThcD|HCD)mgf.*"
    pFile = re.compile(patternFile)
    
    try:
        fr = open(os.path.normpath(spectrumFile), 'r')
    except:
        raise

    try:
        fw0 = open(os.path.normpath(para['outputFile']), 'w')
        fw0.close()
    except:
        raise

    try:
        matchFile = pFile.search(spectrumFile)
        if matchFile is None:
            return
        
        fw = open(os.path.normpath(para['outputFile'].replace(".mgf", "_{0}.mgf".format(matchFile.group(1)))), 'w')
        
        for line in fr:
            if line.startswith("TITLE="):
                line = line.strip()
                match = p.search(line)
                fw.write("""{0} "{1}" {2}\n""".format(match.group(1), para['inputFile'], match.group(3)))
            else:
                fw.write(line)
               
    except:
        raise


def post_process(para):
    outFiles = []

    # find OutputFile in tempDirectory
    for root, dirs, files in os.walk(os.path.normpath(para['tempPath'])):
        #count number of mgf files exeeding 10KB
        
        for name in files:
            fullPath = os.path.join(root, name)
           
            if name.endswith(para['fileExtension']) and os.path.getsize(fullPath) > 10000:
                outFiles.append(name)

    for fileName in outFiles:
        fullPath = os.path.join(root, fileName)
        print fullPath
        replace_title(fullPath, para)
                                     
def compose_batch_file(para):
    myCommand = """
echo HELLO
copy {0} {1}

"{2}" -c Rawfiles 
"{2}" -a Rawfiles {1}
"{2}" -e Rawfiles ANY "{3};{4}"

""".format(os.path.normpath(para['inputFile']),
           os.path.normpath(para['tempFileName']),
           os.path.normpath(para['pdCmd']),
           os.path.normpath("c:/FGCZ/fcc/pd/p1352_AllFT_mgfs.pdProcessingWF"),
           os.path.normpath("c:/FGCZ/fcc/pd/CWF_minimal.pdConsensusWF")) 


    batchFilename = os.path.normpath("{0}/{1}".format(para['tempPath'], 'pd21_temp.bat'))

    try:
        f = open(batchFilename, 'w')
        f.write(myCommand)
        f.close()
    except:
        raise
        
  
    try:
        os.system(batchFilename)
        print "\tDONE"
        #print("{} succeeded. {} has been converted.".format(myCommand, inputFile))
    except OSError as e:
        print("{} failed!. {} has NOT been converted. Error {} occured".format(myCommand, inputFile, e))
        raise

if __name__ == "__main__":
    tempDir = None
    try:
        tempDir = tempfile.mkdtemp("-pd21", dir=os.path.normpath('d:/tmp/'))
    except:
        raise
    
    para = {
        'pdParamFile': sys.argv[1],
        'inputFile': sys.argv[2],
        'outputFile': sys.argv[3],
        'pdCmd': "c:/Program Files/Thermo\Proteome Discoverer Daemon 2.1/System/Release/DiscovererDaemon.exe",
        'tempFileName': os.path.normpath("{0}/temp.raw".format(tempDir)),
        'tempPath': tempDir,
        'fileExtension': 'mgf'
        }

    #para['tempFileName']='D:\\tmp\\tmpownkuk-pd21\\temp.raw'
    #para['tempPath']='D:\\tmp\\tmpownkuk-pd21'
    compose_batch_file(para)
    post_process(para)
    
