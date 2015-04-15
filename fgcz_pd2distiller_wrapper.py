# $HeadURL: http://fgcz-svn.uzh.ch/repos/fgcz/stable/proteomics/fcc/fgcz_pd2distiller_wrapper.py $
# $Id: fgcz_pd2distiller_wrapper.py 7374 2015-04-14 06:53:01Z cpanse $
# $Date: 2015-04-14 08:53:01 +0200 (Tue, 14 Apr 2015) $
# $Author: cpanse $


from random import randint
import sys
import shutil
import os
import stat
import fileinput
import re

class FgczProteinDiscovererWrapper(object):
    """
    python script for converting Thermo raw files to scpectrum files with Proteome Discoverer Daemon
    Prerequisite is a running Proteome Dicoverer process

    Author:
        Simon Barkow-Oesterreicher, Functional Genomic Center Zurich, 20140226 

    Maintainer:
        Christian Panse <cp@fgcz.ethz.ch>
    """

    para = dict()

    def __init__(self):
        self.para['tempPath'] = "E:/pd_temp/"
        # todo(cp); use tf = tempfile.NamedTemporaryFile()
        self.para['randomName'] = str(randint(0, 100000))
        try:
            pass
        except:
            raise

    def set_para(self, key, value):
        """ class parameter setting """
        self.para[key] = value

    def print_para(self):
        """ print class parameter setting """
        for k, v in self.para.items():
            sys.stdout.write("{0}\t=\t{1}\n".format(k, v))

    def run(self):
        self.para['intermediateFile'] = "{0}/{1}.raw".format(self.para['tempPath'], self.para['randomName']) 
        self.para['cmd'] = "D:/Program Files/Thermo/Discoverer Daemon 1.4/System/Release/DiscovererDaemon.exe" 
        self.para['cmdPara'] = "-p C:\\FGCZ\\fcc\\proteomeDiscoverer/{0} {1}".format(self.para['pdParamFile'], self.para['intermediateFile'])

#    try:
#        shutil.copy(inputFile, intermediateFile) 
#    except OSError as e:
#        print("{} could not be copied into {}! Error {} occured".format(inputFile, tempPath, e))
#        sys.exit()
#            
#    myCommand = r'"D:\\Program Files\\Thermo\Discoverer Daemon 1.4\\System\\Release\\DiscovererDaemon.exe" -p C:\\FGCZ\\fcc\\proteomeDiscoverer\\{1} {2}'.format(pdParamFile, intermediateFile)
#
#    try:
#        os.system(myCommand)
#        print("{} succeeded. {} has been converted.".format(myCommand, inputFile))
#    except OSError as e:
#        print("{} failed!. {} has NOT been converted. Error {} occured".format(myCommand, inputFile, e))
#        sys.exit(1)

    def replaceTitle(self, spectrumFile):
        """
        Goes through all the lines and adapt the TITLE field working for FGCZ 
        BFabric mapping.
        Input and Output files are the same.
        """
        for line in fileinput.input(spectrumFile, inplace=True):
            line = line.strip()
            #if (line.startswith("TITLE=")):j
            #    print(line.replace(line, "{0} [{1}]".format(line, inputFile))
            #else:
            #    print(line)

    def reshape(self):
        outFiles = []
        #find OutputFile in tempDirectory

        for root, dirs, files in os.walk(tempPath):
            # count number of mgf files exeeding 10KB
            # outFiles = map(lambda f: os.path.join(root, f), files)
            # outFiles = filter(lambda f: f.startswith(self.para['randomName']), outFiles)
            # outFiles = filter(lambda f: f.endswith(self.para['fileExtension']), outFiles)
            # outFiles = filter(lambda f: os.path.getsize(self.para['fullPath']) > 10000, outFiles)
            for name in files:
                fullPath = os.path.join(root, name)
                if name.startswith(self.para['randomName']) \
                    and name.endswith(self.para['fileExtension']) \
                    and os.path.getsize(self.para['fullPath']) > 10000:
                    outFiles.append(name)

            for fileName in outFiles:
                fullPath = os.path.join(root, fileName)
#        replaceTitle(fullPath)
#        outFileName=outputFile
#        if (len(outFiles) > 1):
#            if not os.path.isfile(outputFile):
#            with open(outputFile, 'w') as outFile: print('dummy mgf file', file=outFile,  end='\r')
#            #change output filename according to node
#            if re.search("\[Node_10\]",fileName):
#            outFileName = os.path.splitext(outputFile)[0] + '__ETDIT' + os.path.splitext(outputFile)[1]
#            elif re.search("\[Node_27\]",fileName):
#                outFileName = os.path.splitext(outputFile)[0] + '__HCDIT' + os.path.splitext(outputFile)[1]
#            elif re.search("\[Node_31\]",fileName):
#                outFileName = os.path.splitext(outputFile)[0] + '__ETDFT' + os.path.splitext(outputFile)[1]
#            elif re.search("\[Node_14\]",fileName):
#                outFileName = os.path.splitext(outputFile)[0] + '__HCDFT' + os.path.splitext(outputFile)[1]
#            elif re.search("\[Node_34\]",fileName):
#                outFileName = os.path.splitext(outputFile)[0] + '__CIDFT' + os.path.splitext(outputFile)[1]
#            elif re.search("\[Node_04\]",fileName):
#                outFileName = os.path.splitext(outputFile)[0] + '__CIDIT' + os.path.splitext(outputFile)[1]
#        try:
#            shutil.copy(fullPath,outFileName)
#            print(fileName + " copied to " + outFileName)
#        except OSError as e:
#            print("{} could not be copied into {}! Error {} occured".format(fullPath,outFileName, e))
#            pass
#

    def clean(self):
        for root, dirs, files in os.walk(self.para['tempPath']):
            for name in files:
                if (name.startswith(randomName)):
                    fullPath = os.path.join(root, name)
                    try:
                        os.chmod(fullPath, stat.S_IWRITE)
                    except OSError as e:
                        print("chmod did not work, {}".format(e))
                    try:
                        os.remove(fullPath)
                        print("{} removed".format(fullPath))
                    except OSError as e:
                        print("{} was not removed! Error {} occured".format(fullPath,e))


if __name__ == "__main__":
    #fileExtension=sys.argv[1]
    #pdParamFile=sys.argv[2]
    #inputFile=sys.argv[3]
    #outputFile=sys.argv[4]
    PD = FgczProteinDiscovererWrapper()
    PD.set_para('fileExtension', sys.argv[1])
    PD.set_para('pdParamFile', sys.argv[2])
    PD.set_para('inputFile', sys.argv[3])
    PD.set_para('outputFile', sys.argv[4])
    PD.run()
    PD.print_para()

    print "HELLO"
    sys.exit(0)

#
#"""
#Example:
#    python fgcz_pd2distiller_wrapper.py
#"""
#
#
#
