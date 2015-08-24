# fcc - [FGCZ](http://www.fgcz.ethz.ch) Converter Control

# description
fcc is a minimalistic workflow engine.
The specification/properties of the program are as follow:
- converting instrument files (e.g. RAW-files) to all kinds of formats
- being generic
- follows fgcz granularity: project, user, instrument, time range
- multi platform, host, task
- configurable through xml file which means new converter by new tag in xml file NO CODE CHANGE!
- stdout and errout logging

# synopsis
## on Microsoft OS
```
python fcc.py --dir=s:\ --output=runme.bat
python fcc.py --dir=s:\ --pattern=".*VELOS_2.*" --exec --loop
C:\FGCZ\fcc>c:\Python32\python.exe fcc.py --dir s:\ --output test.bat --pattern "^[Ss]:\\p[0-9]+\\.+"
```

## on UNIX

```
# use python 2.7.x
python fcc_run_linux.py --dir=/srv/www/htdocs/Data2San/p720 --output=runme.bat
```
## arguments:
-output writes a batch file to run later manually
-dir specifies the directory to crawl
-exec automatically triggers the execution of the generated converter commands
-loop the FCC automatically restarts after it has finished one crawling round

#author
Simon Barkow-Oesterreicher <simon@uberchord.com> and Christian Panse <cp@fgcz.ethz.ch>
    
#maintainer
<cp@fgcz.ethz.ch>

# see also / references
- [doi:10.1186/1751-0473-8-3](http://www.scfbm.org/content/8/1/3/abstract), [PMC3614436](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3614436/),
PMID: 23311610

# configuration

```xml
 <controllerRuleSet>
    <rule converterID='000' 
    project='p103' 
    omics='Proteomics' 
    user='' instrument='ORBI_2' 
    beginDate='20080901' 
    endDate='20990101' 
    keyword="iTRAQ">
    </rule>
</controllerRuleSet>
    
<converterList>
    <converter converterID='000' 
    converterDir= 'mgf__low_res_MS2_iTRAQ' 
    converterCmd='cscript "C:\FGCZ\fgcz-proteomics\stage\mascot_distiller\fgczRaw2Mgf.vbs"'         converterOptions='"C:\FGCZ\fgcz-proteomics\stage\generalRawFileConverterRobot\MascotDistillerOPTs\Orbitrap_low_res_MS2_iTRAQ.opt"' 
    toFileExt='.mgf' 
    hostname='fgcz-s-034'> 
    </converter>
</converterList>
```
    
# example

run fcc on a linux box and fake the hostname
```bash
./fcc.py --output runme.bat --hostname fgcz-s-021
```

# TODO

## use lxml 
## replace XML fcc_config by something easier to config, e.g. see rc files in obsd.
## provide cmd argv for log file name
## unit test


# branches
## fgcz

should run out of the box on a linux machine.
