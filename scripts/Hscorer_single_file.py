import sys
import getopt
import numpy as np
import os
import re



class IonFound:
    def __init__(self,mz,intensity):
        """
        @brief simple object to hold the measured ions data
        
        """
        self.intensity = intensity 
        self.mz = mz

class SeriesObject:
    def __init__(self):
        """
        @brief simple object to different series (and generator for theoretical 
        fragments) 
        """
        self.name = 'n' #name used to index series later on
        self.version = 0

class IonProcessing:
    """
    @brief class holding all processing methods to match theoretical 
    fragments to measured ions
    """
    def __init__(self, masses):
        self.masses = masses
        bs = SeriesObject()
        bs.series = 'b'
        bs.name = 'b'
        ys = SeriesObject()
        ys.series = 'y'
        ys.name = 'y'
        self.standardseries = [bs,ys]

    def mergeCleavageList(self, cleavagelist, seq):
        """
        @brief merges the cleavagepatterns generated from the b & y series (in cleavagelist)
        @param cleavagelist <list>: list of cleavage patterns from the 2 series
        @param seq <string>: peptide sequence
        @return fullcleavagepattern <string>: complete cleavage pattern after merging 
        """
        cleavageline = ''
        fullcleavagepattern = ''
        runningindex = 1
        for idx in range(0,len(seq)):
            cleave = 0 
            for cl in cleavagelist:
                if cl[idx+runningindex] == '|':
                    cleave = 1
                    break
            runningindex+=1
            
            if cleave == 1:
                fullcleavagepattern+=seq[idx]+'|'
            else:
                fullcleavagepattern+=seq[idx]
        return fullcleavagepattern
    def doubleB(self,mass):
        
        x = mass+self.masses['Proton']
        
        x = x/2.0
        x[x<=(self.masses['Proton']/2.0)] = 0
        return x
    def calculateFragmentIons(self,seq, modstring):
        """
        @brief calculates theoretical b & y ion series for sequence 'seq' 
        taking into account modifications encoded within modsstring
        @param seq <string>: peptide 
        @param modstring <string>: series of numbers indicating type & position 
        of modification on peptide sequence
        @return allbseries+allyseries <list>: list of python generators to make
        theoretical fragments for sequence given
        """

        modstring_noterm = modstring[1:]
        modstring_noterm = modstring_noterm[:-1]
        startmass = 0 
        #b-series should start with mass of one proton:
        startmass+=float(self.masses['Proton']) 
        if int(modstring[0]): #if there's an N-terminal modificaiton
            startmass+=self.masses[str(modstring[0])]

        currentIdx = 0
        theseseries = self.standardseries #standard series is just 'y' & 'b'

        allbseries = [x for x in theseseries if x.series =='b']
        for bs in allbseries:
            mass = startmass
            bs.generator = self.getSeqMass(mass,currentIdx,modstring_noterm, seq)
            bs.seriesindex =currentIdx/len(seq)
            currentIdx+=len(seq)            

         ######## And now for the Y-series ##################                
        startmass = 0
        #y-series always starts with 1 water
        startmass+=(self.masses['Oxygen']+(self.masses['Hydrogen']*2))        
        
        startmass+=float(self.masses['Proton']) 
        seq = list(seq)
        seq.reverse()
        modstring_noterm = list(modstring)
        modstring_noterm.reverse()
        modstring_noterm = ''.join(modstring_noterm[1:])

        modstring_noterm = modstring_noterm[:-1]        

        seq = ''.join([x for x in seq])
        allyseries = [x for x in theseseries if x.series =='y']
        for ys in allyseries:
            mass = startmass
            ys.generator = self.getSeqMass(mass,currentIdx,modstring_noterm, seq)
            ys.seriesindex =currentIdx/len(seq)
            
            currentIdx+=len(seq)
        return allbseries+allyseries






    def getSeqMass(self,startmass, idxstart,modstring_noterm ,sequence):
        """
        @brief creates python generator to calculate cumulative mass of a 
        peptide sequence, adding masses of modifictions where appropriate
        @param startmass <float>: mass before adding residue specific masses
        @param idxstart <int> : current index
        @param modstring_noterm <string>: series of numbers indicating type & position 
        of modification on peptide sequence
        @param sequence <string>: peptide sequence for which mass list is being generated
        @return mass, index for use in the generator
        """
        seqlen = len(sequence) -1
        idx = 0 
        modcounter = 0 
        mass = startmass
        while idx < len(sequence)-1:
            aa = sequence[idx]
            if int(modstring_noterm[idx]):
                modmass = self.masses[str(modstring_noterm[idx])]

                mass+=modmass           
            mass+=self.masses[aa]

            yield mass, idx+idxstart
                    
            idx+=1

        yield 0, idx+idxstart
    
    
    def getCleavagePattern(self, foundresidueset, seq, series):
        """
        @brief creates the cleaveage pattern based on list of residue 
        positions calculated using peptide sequence.
        @foudnresidueset <list> : list of matched residue positions
        @seq <string>: peptide sequence
        @series <string>: whether foundresidueset from 'b' or 'y' series
        @return <string>: peptide sequence with mapped cleavage pattern
        """
        if series.series == 'y':
            seq = list(seq)
            seq.reverse()            
        cleavageline = ''
        for idx, r in enumerate(seq):
            #for each position in foundresidueset add a pipe if found
            if idx in foundresidueset:        
                cleavageline+=r+'|'
            else:
                cleavageline+=r+'-'
        if series.series == 'y':
            cleavageline = list(cleavageline[:-1])
            cleavageline.reverse()
            cleavageline.append('-')
            return ''.join(cleavageline)
        else:
            return cleavageline

    def matchFragmentIons(self,ionslist,  fragindex, fragments,tolerance,units):
        """        
        @brief Theoretical fragment arrays created above are converted into a matrix.
        For each, a second matrix is created representing the  same masses, but doubly charged
        All theoretical ions matching to the real ion within given tolerance are recorded
        in the temp matrix.  The intersection of matches for all ions is then mapped back to the 
        corresponding residues and returned
        @param ionslist : list of IonFound objects 
        @param fragindex: index mapping fragment to corresponding residue in peptide
        @param fragments: array of theoretical masses for the series
        @param tolerance : user-supplied tolerance
        @param units: user-supplied tolerance units (ppm or da ONLY)
        @return matchedresidues: all theoretical residue positions explained by ions 
        """
        
        basic_matrix = np.array(fragments)
        double = self.doubleB(basic_matrix)
        double = basic_matrix*0

        allposs = np.vstack((basic_matrix,double)) #stack doubly-charged matrix on top of basic (singly-charged masses) matrix
        
        multifragidx = []
        
        for r in range(0,2):  #add fragmentidex for both levels (double & single charged)
            multifragidx.append(fragindex)
        allpossidx = np.vstack(multifragidx)
        true = basic_matrix == -123 #true-match matrix: same size as basic matrix.  Filled all false as no ions can be minus
        for i in ionslist:
            
            if units == 'ppm':
                #ie need to convert to ppm first then filter
                #temp is a boolean matrix where the tested ion matches theoretical masses within tolerance
                temp = abs((i.mz-allposs) / i.mz)*1000000 < tolerance
            elif units =='da':
                #no need to convert to ppm first, just test
                #temp is a boolean matrix  where the tested ion matches theoretical masses within tolerance
                temp = abs(i.mz-allposs)  < tolerance
            else:
                sys.exit('dont know what to do with units %s exiting' %units)
            true = true|temp            
        matchedresidues = list(allpossidx[true])

        return matchedresidues
            
            
def getPeptideData(file):
    """
    @brief extracts data from a tab-delimited file of format ID\tPeptide sequence\tmodlist
    @param <open file handle>: file handle from which to extract the data
    @retun <dict>: dictionary indexed by spectrum ID, containing peptide string, mod locations 
    """
    dataDict = dict()
    for line in file:
        items = line.strip().split('\t')
        
        id = items[0]
        peptide = items[1]
        if len(items) > 2:
            modstring = items[2]
        else:
            modstring = '0'*(len(peptide))
            #print modstring
        dataDict[id] = (peptide,modstring)
    sys.stdout.write('There are %i records for which to create H-Score\n' %len(dataDict))
    
    return dataDict

        
def doDeisotopeAll(totlist,maxcharge):
    """
    
    @brief set intensity to zero for all but the monoisotopic peak of an 
    isotopic cluster for all chargestates up to that of the parent ion
    @param totlist <list>: list of tuples : mz/intensity
    @param maxcharge <int>: parent chargestate
    @return totlist <list>: list where intensity zero pairs filtered away
    """
    maxstart = maxcharge
    neutron = 1.0033548
    newlist = [] #add extra column to track what starting chargestate was
    for x in totlist:
        x.extend([0])
        newlist.append(x)
    totlist = newlist

    totlist.reverse()

    for index, x in enumerate(totlist):# deisotoping

        for y in totlist[index:]:
            charge = maxstart

            while charge:
                if abs(x[0] - y[0] - neutron  / charge)<0.02 and y[1]>(x[1]/2.0) and x[0] > 132:
                    x[1] = 0
                    x[2] = charge
                    y[2] = charge   
                charge-=1
            if x[0]-y[0]>2:
                break
    totlist = [x for x in totlist if x[1] > 0]        #remove all datapoints where intensity zero or less               
    return totlist

def removeReporterIons(totlist,repions):
    """
    @brief removes spectra from totlist if they are from reporter ions (within
    0.02 Da tolerance)
    @param totlist <list>: list of all mz/intensity tuples
    @param repions <list>: list of reporter ion masses
    @return totlist <list>: new list of mz/intensity tuples without reporter ions 
    """
    for repion in repions:
        totlist = [x for x in totlist if not (abs(float(x[0]) - repion ) <= 0.02 )]
    return totlist
def doDeconvolute(totlist,proton, parentmass):
    """
    @brief recalcultes mz to mz of chargestate 1 when chargestate > 1
    @param totlist <list>: list of m/z intensity and parent chargestate
    @param proton <float>: value for proton
    @param parentmass <float>: mass of parent ion
    @retrun 
    """
    totlistnew = []
    for mz, intensity, chargestate in totlist:
        if int(chargestate) > 1:
            monoisotope = mz * float(chargestate) - ((float(chargestate) -1)  * proton) 
            #if newly calculated monoistopic mass is greater than parent mass, 
            #something went wrong, leave as original mz
            if monoisotope - parentmass > -1:
                monoisotope = mz
        else:
            monoisotope = mz
        totlistnew.append((monoisotope,intensity,chargestate))
    return totlistnew    
def doHCDdeconvolution(myMGF,newMGF,reporterions,masses,dataIN,dataOUT,tolerance,units):
    """
    @brief main function of the script.  Take existing MGF file, extract the m/z 
    intensity pairs and filter / deconvolute these.  A new MGF file is created with 
    these amended data.  For ions where the spectrum ID matches a spectrum ID
    given in the dataIN file, the Hscore & cleaveage pattern will be generated
    and written to a file.
    @param myMGF <string>: name of file being processed
    @param newMGF <string>: name of new MGF file
    @param reporterions <list>: list of reporterion masses to filter away
    @param masses <dict>: dictionary of residues / modification and masses
    @param dataIN <string>: name of file containing peptide & modstring with \
    spectrum id as first field
    @param dataOUT <string>: name of file to which to write cleaveage pattern / 
    hsocre for data given in dataIN file
    @param tolerance <float>: tolerance for matching theretical fragment to ions 
    @param units <string>: units which tolerance is measured in (da OR ppm)
    @return titlecounter <int>: total number of spectra in MGF file processed 
    (title parsed)
    """
    fh = open(myMGF,'r')
    fho = open(newMGF,'w')
    try:
        datain = open(dataIN,'r')
        peptideData = getPeptideData(datain)
        datain.close()
    except IOError:
        sys.stderr.write('cannot find file %s, will NOT calculate Hscore' %dataIN)
        peptideData = dict() 
    
    
    dato = open(dataOUT, 'w')
    dato.write('spectrum id\tpeptide\tmodstring\tcleavagepattern\thscore\n')
    proton = masses['Proton']
    ch = 0
    datastr = []
    totlist = []
    ignore = 0
    titlecounter =0
    lastid = None
    lineID = ""
    specid = ""
    scanNumber = ""
    rt = ""
    mz = 0
    for line in fh:
        if line != '\n':
            lineID = line.strip().split('=')[0]
            values = line.replace(lineID + "=", "");
        if lineID == 'TITLE': 
            titlecounter+=1
            myline = line.strip()
            # The following lines do only work with mascot distiller mgfs
            #name = values.split('[')[0]
            #specid = values.split('[')[1][:-1]
            #specid = specid[:-1]
            # This line should work with all mgfs
            specid = values
            scanNumber = line.split('.')[1]
            p = re.compile( '\\\\mgf__proteowizard_.+\\\\')
            myMGF = re.sub(p, '\\\\', myMGF)
            # write title line after getting the rtinseconds value
        if lineID == 'RTINSECONDS':
            rt = float(values)
            datastr.append(myline + ' Scan ' + scanNumber +  ' (rt='+str(rt/60)+') [' + myMGF.replace(".mgf.tmp",".raw") + ']')
            datastr.append(line.strip())   

        if lineID == 'CHARGE':
            chp = values.split('+')
            ch = float(chp[0])
            datastr.append(line.strip())
            # include scanNumber here:
            datastr.append("SCANS=" + scanNumber)
        if lineID == 'PEPMASS':
            mz = float((values.split(' '))[0])
            datastr.append(line.strip())
        if lineID == 'END IONS':
            if len(totlist) > 0 :
                sys.stdout.flush()
                if (mz*ch < 16000) :
                    fho.write('BEGIN IONS\n')
                    fho.write('\n'.join(datastr))
                    fho.write('\n')
                    totlist = removeReporterIons(totlist,reporterions)
                    totlist = doDeisotopeAll(totlist,ch)
                    parentmass = mz * ch - (proton * (ch -1)) 
                    totlist = doDeconvolute(totlist,proton,parentmass)
                    totlist.sort()
                    for v,y,z in totlist:
                        fho.write('%f %f %i\n' %(v,y,z))

                    fho.write('END IONS\n\n')
                else:
                    print("PEPTIDE REMOVED" + " MZ:" + str(mz) + " CHARGE:" + str(ch))
                sys.stdout.flush()
        
                if specid in peptideData:
                    sequence = peptideData[specid][0]
                    ptmstring = peptideData[specid][1]
                    cleavagepattern, hscore = calcHScore(sequence, ptmstring,totlist,masses,tolerance,units)
                    dato.write(specid+'\t'+sequence+'\t'+ptmstring+'\t'+cleavagepattern+'\t'+str(hscore)+'\n')
        
                #only add mz / intensity data, other lines not added

            totlist = []             
            datastr = []
            lastid = specid
        expmint = line.split()
        try:
            totlist.append([float(expmint[0]), float(expmint[1])])
        except ValueError:
            #don't add non float data to the totlist
            pass
    #lest we forget the final spectrum!
    if len(totlist) > 0 :

        sys.stdout.flush()
        fho.write('BEGIN IONS\n')
        fho.write('\n'.join(datastr))
        fho.write('\n')
        totlist = removeReporterIons(totlist,reporterions)
        totlist = doDeisotopeAll(totlist,ch)
        parentmass = mz * ch - (proton * (ch -1)) 
        totlist = doDeconvolute(totlist,proton,parentmass) 
       
        for v,y,z in totlist:
            fho.write('%f %f %i\n' %(v,y,z))
        fho.write('END IONS\n\n')            
        if specid in peptideData:
            sequence = peptideData[specid][0]
            ptmstring = peptideData[specid][1]
            cleavagepattern, hscore = calcHScore(sequence, ptmstring,totlist,masses,tolerance,units)
            dato.write(specid+'\t'+sequence+'\t'+ptmstring+'\t'+cleavagepattern+'\t'+str(hscore)+'\n')
        
    fho.close()
    dato.close()
    return titlecounter

def calcHScore(sequence, ptmstring, ionlist, masses,tolerance,units):
    """
    @brief combines data from sequence and ptmstring to create
    theoetical peptide fragments and then mathces filtered ions
    to these within the given tolerance to generate cleavage pattern 
    and H-score for peptides in list
    @param sequence <string>: peptide sequence
    @param ptmstring <string>: string denoting modification type and position
    on sequence
    @param ionlist <list>: list containing tuples of m/z / intensity after filtering
    @param masses <dict>: dictionary with residue/modification and corresponding mass
    @param tolerance <float>: tolerance for matching theretical fragment to ions 
    @param units <string>: units which tolerance is measured in (da OR ppm)
    @return cleavagepattern <string>: modified version of the peptide showing which
    cleavage sites are explained by found ions
    @return hscore <int>: H-score based on number of explained cleavage sites
    """
    ionObjlist = []
    ionp = IonProcessing(masses) #this object has all methods to create theoretical fragments
    
    #convert ionlist tuples into IonFound object (makes later handling easer & can be extended)
    for i in ionlist:
        
        j,k = i[0],i[1]
        if float(k) > 0: #only add if intensity > 0!
            f = IonFound(float(j),float(k))
            ionObjlist.append(f)
    cleavageList = []
    #create python generators for all series
    fragmentgenerators = ionp.calculateFragmentIons(sequence,ptmstring)
    for series in fragmentgenerators:
        fragments = []
        fragindex = []
        foundresidueset = set()
        for mass, idx in series.generator:
            #generator gives fragment mass & index on the peptide sequence
            fragments.append(mass)
            fragindex.append(idx)
            
        #create cleavage pattern for each series
        foundresidues = ionp.matchFragmentIons(ionObjlist,fragindex,fragments,tolerance,units)
        foundresidueset.update([x%len(sequence) for x in foundresidues])
        cleavageList.append(ionp.getCleavagePattern(foundresidueset,sequence,series))
        

    cleavagepattern = ionp.mergeCleavageList(cleavageList, sequence)
    hscore = cleavagepattern.count('|')
    
    if hscore +1 == len(sequence): # bonus for all explained cleavage site
        hscore+=3
    elif hscore+2 ==len(sequence):# bonus for 1 non explained cleavage sites
        hscore+=1      
    return cleavagepattern,  hscore
   
def getMasses(file=None):
    """
    @brief extracts data from user-supplied tab-delimited file containing 
    modifications / AA masses.
    @parma file <string>: filename from which to extract the masses
    @return <dict>: key value pairs of residue / mod : mass
    
    """
    itraq_masses=[114.11120,115.10830,
                  116.11160,117.11500]
    tmt_masses=[126.12830,127.13160,128.13500,
                129.13830,130.14170,131.13870]
    proton = 1.0072760000000001

    if not file:


        masses = {
            'Oxygen': 15.994915000000001,'N_term': 1.007825, 
            'Electron': 0.00054900000000000001,'Carbon': 12.0, 
            'A': 71.037114000000003,'C': 160.03064900000001,
            'B': 114.53494000000001, 'E': 129.04259300000001,
            'D': 115.026943, 'G': 57.021464000000002,
            'F': 147.06841399999999,'I': 113.084064,
            'H': 137.05891199999999, 'K':357.257895, 
            'J': 0.0, 'M': 131.04048499999999, 
            'L': 113.084064, 'O': 0.0, 
            'N': 114.04292700000001, 'Q': 128.05857800000001,
            'P': 97.052763999999996, 'S': 87.032027999999997,
            'R': 156.101111, 'U': 150.95364000000001,
            'T': 101.047679, 'W': 186.07931300000001,
            'V': 99.068414000000004, 'Y': 163.06332900000001,
            'X': 111.0, 'Z': 128.55059, 
            'C_term': 17.002739999999999,
            'Nitrogen': 14.003074, 'Hydrogen': 1.007825}
    else:
        masses = dict()
        fh = open(file, 'r')
        for line in fh: 
            items = line.strip().split('\t')
            try: 
                float(items[1])
                residue, mass = items[0], items[1]
                masses[residue] = float(mass)
            except ValueError:
                #line not contining valid data, ignore
                pass
    masses['itraq'] = itraq_masses
    masses['tmt'] = tmt_masses
    masses['Proton'] = proton
    return masses

def usage():
    print ("""
USAGE:
  Hscorer.py [OPTION]

takes options from command line otherwise uses in-built defaults (given inside '<>')
OPTION:
    --inp : input mgf file 
    --outp : output hmgf file 
    --tolerance <20.0> :tolerance for ion to theoretical fragment matching
    --units <ppm>: unit of tolerance given above (ppm or da only)
    --quantmeth <>: quantification method used (tmt/itraq only); these reporter ions are 
    removed from the MGF file
    --massfile <>: file containing masses of amino acids & any modifications (loosely based on 
    Mascot's dat file 'mods' sections.
""")

def hasOption( option, options):
    for o, a in options:
        if o == option:
            return 1
    else:
        return 0

def getOption( option, options):
    for o, a in options:
        
        if o == option:
            
            return a   


def getconfig(opts=None):
    """Build config container
    switch between test/production mode
    """
    configops = dict(inp=None,outp=None,tolerance=None,
                              units = None,quantmeth=None,
                              massfile=None,reporterions=None)
    for configItem in configops.keys():
        
        if hasOption('--' + configItem, opts):
            configops[configItem] = getOption('--' + configItem,
                    opts)
    
    return configops
        
    
if __name__=="__main__":
    

    try:
        opts, args = getopt.getopt(sys.argv[1:], 'help',['inp=','outp=','tolerance=',
                                                         'units=','quantmeth=',
                                                         'massfile=','reporterions='])
    except getopt.GetoptError:
        usage()
        sys.exit(2)
    
    cfg = getconfig(opts)
    
    if cfg['inp']:
        inp = cfg['inp']
    
    if cfg['outp']:
        outp = cfg['outp']
        
    if cfg['tolerance']:
        tolerance = float(cfg['tolerance'])
    else:
        tolerance = 20.0 #mass accuracy in PPM
    if cfg['units']:
        units = cfg['units']
    else:
        units = 'ppm'
    if cfg['quantmeth']:
        if cfg['quantmeth'] not in ('tmt', 'itraq'):
            sys.exit('quantmethod %s not supported ( I dont know which ions to remove ), terminating' % cfg['quantmeth'])
        quantmeth = cfg['quantmeth']
    else:
        quantmeth = None
    if cfg['massfile']:
        massfile = cfg['massfile']
    else:
        print('No massfile provided using default settings (AA values ONLY - non modification masses)')
        massfile =None
    reporterions = []
        
    masses = getMasses(massfile)
    if quantmeth in masses:
        reporterions = masses[quantmeth]   
    print('running with following parameters:\non file %s\noutput file %s\nwith tolerance %f\nunits %s\nquantmeth %s\nmassfile %s\nreporterions to remove %s\n'   %(inp,outp,tolerance,units,quantmeth,massfile,','.join([str(x) for x in reporterions])))
    print(os.path.abspath(inp))
    
    processingList = []
    processingList.append(os.path.abspath(inp)) 
    outpDir = os.path.dirname(outp)

    
    for file in processingList:
        
        #newfilename = os.path.join(outpDir, os.path.abspath(file.name.replace('.mgf.tmp','.hmgf')))
        print(os.path.basename(file))
        datain = os.path.join(outpDir,os.path.basename(file).replace('.mgf.tmp','.data'))
        dataout = os.path.join(outpDir, os.path.basename(file).replace('.mgf.tmp','_hscore.data'))
        #def doHCDdeconvolution(myMGF,newMGF,reporterions,masses,dataIN,dataOUT,tolerance,units):
        doHCDdeconvolution(file,outp,reporterions,masses,datain,dataout,tolerance,units)
        

