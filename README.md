# FCC - [FGCZ](http://www.fgcz.ethz.ch) Converter Control

# Description
fcc is a minimalistic workflow engine.
The specification/properties of the program are as follow:
- converting instrument files (e.g. RAW-files) to all kinds of formats
- being generic
- follows fgcz granularity: project, user, instrument, time range
- multi platform, host, task
- configurable through xml file which means new converter by new tag in xml file NO CODE CHANGE!
- stdout and errout logging

# Synopsis
## on Microsoft OS
```
c:\Python27\python.exe fcc_run_windows.py --output __runme.bat --ncpu 2 --exec
```

## on UNIX

```
# use python 2.7.x
$ python fcc_run_linux.py --hostname fgcz-s-021 --output out-20150824.bash --ncpu 1 --exec --loop
```
## Arguments:

- `--output <file>` writes a batch file to run later manually
- `--exec` automatically triggers the execution of the generated converter commands
- `--loop` the FCC automatically restarts after it has finished one crawling round
- `--ncpu` number of jobs should be used
- `--hostname <host configuration>` tells fcc which host configuration should be used.

# Authors

Simon Barkow-Oesterreicher <simon@uberchord.com> and Christian Panse <cp@fgcz.ethz.ch>
    
#maintainer

<cp@fgcz.ethz.ch>


# Configuration

[the current FCGZ configuration can be found here.](http://fgcz-data.uzh.ch/config/fcc_config.xml)

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

# TODO

- use [lxml](http://lxml.de/) or even better use [YAML](http://yaml.org/) to ease implementation/configuration
- provide cmd argv for log file name
- unit test - no glue how to do that 
- ensure that no process is killed by CTRL-C. add a `pool wait` at the right position on the code when CTRL-C is pressed or give at least a warning.
- write a  pid file to ensure that only one fcc is running

# branches
- [fgcz](https://github.com/fgcz/fcc) (default)

# see also 

- FCC - An automated rule-based processing tool for life science data.
Barkow-Oesterreicher S1, Türker C, Panse C. Source Code Biol Med. 2013 Jan 11;8(1):3. doi: 10.1186/1751-0473-8-3. [doi:10.1186/1751-0473-8-3](http://www.scfbm.org/content/8/1/3/abstract), [PMC3614436](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3614436/),
PMID: 23311610

- [FCC @ fgcz configuration](http://fgcz-data.uzh.ch/config/fcc_config.xml)
- 
# Have Fun
![Mt Pico](https://lh4.googleusercontent.com/-Y71pZMTpDHo/VVBtoGacseI/AAAAAAAACJA/kw6Lli7wTlA/w1052-h884-no/15%2B-%2B1)
