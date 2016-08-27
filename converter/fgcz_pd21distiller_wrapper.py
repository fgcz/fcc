#!/usr/bin/python
# -*- coding: latin1 -*-

# python script for converting Thermo raw files to scpectrum
# files with Proteome Discoverer Daemon 2.1
# Prerequisite is a running Proteome Dicoverer process

# AUTHORS:
# Simon Barkow-Oesterreicher, Functional Genomic Center Zurich, 20140226
# Christian Panse <cp@fgcz.ethz.ch>
# 2016-07-27
# 2016-08-24,26
#
# SOURCE
# https://github.com/fgcz/fcc/blob/fgcz/converter/fgcz_pd21distiller_wrapper.py


# This file is part of the fcc package.
# https://github.com/fgcz/fcc

# fcc is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# fcc is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with fcc.  If not, see <http://www.gnu.org/licenses/>.

import sys
import shutil
import os
import stat
import fileinput
import re
import signal
import tempfile


def replace_title(spectrumFile, para):
    """
    Goes through all the lines and adapt the TITLE field working for FGCZ
    BFabric mapping.
    Input and Output files are the same.
    """
    pattern = "(^TITLE=File:).(\"[-_:\.\\a-zA-Z]{3,}\")(;\ Spe.+)$"
    p = re.compile(pattern)
    patternFile = ".*_(FTETD|FTEThcD|FTHCD|ITETD|ITCID|ITHCD|ITEThcD)_.*mgf"
    pFile = re.compile(patternFile)

    try:
        fr = open(os.path.normpath(spectrumFile), 'r')
    except:
        raise

    try:
        # create an empty file
        fw0 = open(os.path.normpath(para['outputFile']), 'w')
        fw0.close()
    except:
        raise

    try:
        match_file = pFile.search(spectrumFile)
        if match_file is None:
            return

        mgf_filename = os.path.normpath(para['outputFile'] \
            .replace(".mgf", "_{0}.mgf".format(match_file.group(1))))

        print "mgf_filename =", mgf_filename
        fw = open(mgf_filename, 'w')
        for line in fr:
            if line.startswith("TITLE="):
                line = line.strip()
                match = p.search(line)
                fw.write("""{0} "{1}" {2}\n""".format(match.group(1), \
                    para['inputFile'], match.group(3)))
            else:
                fw.write(line)
    except:
        raise


def post_process(para):
    """
    reshape PD output;
    1. replaces all TITLE fields by copying the original input
    2. copy the mgf to its final destination
    """
    result_mgf_file_list = []

    # find OutputFile in tempDirectory
    for root, dirs, files in os.walk(os.path.normpath(para['tempPath'])):
        #count number of mgf files exeeding 10KB
        for name in files:
            full_path = os.path.join(root, name)
            if name.endswith(para['fileExtension']) \
                and os.path.getsize(full_path) > 10000:
                print
                print "adding file", name, "for post processing ..."
                result_mgf_file_list.append(name)

    for file_name in result_mgf_file_list:
        full_path = os.path.join(root, file_name)
        print full_path
        replace_title(full_path, para)

def compose_batch_file(para):
    """
    composes a Microsoft windows batch file.
    """

    pd_batch_command = """
echo HELLO
copy {0} {1}

"{2}" -c Rawfiles 
"{2}" -a Rawfiles {1}
"{2}" -e Rawfiles ANY "{3};{4}"

""".format(os.path.normpath(para['inputFile']),
           os.path.normpath(para['tempFileName']),
           os.path.normpath(para['pdCmd']),
           os.path.normpath("c:/FGCZ/fcc/pd/p1352_AllIT_FT_mgfs.pdProcessingWF"),
           os.path.normpath("c:/FGCZ/fcc/pd/CWF_minimal.pdConsensusWF"))

    batch_filename = os.path.normpath("{0}/{1}".format(para['tempPath'], \
        'pd21_temp.bat'))

    try:
        f = open(batch_filename, 'w')
        f.write(pd_batch_command)
        f.close()
    except:
        raise
    try:
        os.system(batch_filename)
        print "\tDONE"
    except OSError as e:
        print("{} failed!. {} has NOT been converted. Error {} occured" \
            .format(pd_batch_command, inputFile, e))
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
        'pdCmd': "c:/Program Files/Thermo/Proteome Discoverer Daemon 2.1/System/Release/DiscovererDaemon.exe",
        'tempFileName': os.path.normpath("{0}/temp.raw".format(tempDir)),
        'tempPath': tempDir,
        'fileExtension': 'mgf'
        }

    #para['tempFileName']='D:\\tmp\\tmpownkuk-pd21\\temp.raw'
    #para['tempPath']='D:\\tmp\\tmpownkuk-pd21'
    compose_batch_file(para)
    post_process(para)

    print "python script done."
# END
