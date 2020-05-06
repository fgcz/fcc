#!/usr/bin/python
import os
import getopt
import sys

import tempfile
import fcc


def create_pidfile():
    try:
        pidfile = "{0}/fcc.pid".format(tempfile.gettempdir())
        if os.path.isfile(pidfile):
            print ("{0} already exists.  exit.".format(pidfile))
            sys.exit(1)
        else:
            with open(pidfile, 'w') as f:
                f.write("fcc is running")
    except: 
        print ("creating {0} failed.".format(pidfile))
        sys.exit(1)

def unlink_pidfile():
    pidfile = "{0}/fcc.pid".format(tempfile.gettempdir())
    try:
        os.unlink(pidfile)
    except:
        print ("removing '{0}' failed".format(pidfile))

if __name__ == "__main__":
    create_pidfile()
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hoepl", [
                                   "help", "output=", "exec", "pattern=", "loop", "hostname=", "ncpu="])
    except getopt.GetoptError as err:
        print (str(err))
        sys.exit(2)

    fcc = fcc.Fcc()

    for o, value in opts:
        if o == "--output":
            fcc.set_para('myOutputFile', value)
        elif o == "--exec":
            fcc.set_para('exec', True)
        elif o == "--loop":
            fcc.set_para('loop', True)
        elif o == "--pattern":
            fcc.set_para('myPattern', value)
        elif o == "--hostname":
            fcc.set_para('hostname', value)
        elif o == "--ncpu":
            fcc.set_para('ncpu',  int(value))
        elif o in ("--help"):
            usage()
            sys.exit(0)
        else:
            usage()
            sys.exit(1)

    crawl_pattern = ['/srv/www/htdocs/Data2San/',
        'p[0-9]{2,4}', '(Proteomics|Metabolomics)',
        '(TIMSTOF|EXPLORIS|QEXACTIVEHF|FUSION|QEXACTIVEHFX|LUMOS|QDA|GCT|G2HD|G2SI)_[0-9]',
        '[a-z]{2,18}_[0-9]{8}(_[-a-zA-Z0-9_]{0,100}){0,1}',
        '[-a-zA-Z0-9_]+.(d|raw|RAW|wiff|wiff\.scan)']

    fcc.set_para('crawl_pattern', crawl_pattern)
    fcc.set_para('min_time_diff', 3600)
    fcc.set_para('max_time_diff', 3600 * 24 * 30 * 6)

    fcc.run()
    unlink_pidfile()
