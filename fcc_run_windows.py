
import getopt
import sys
from fgcz import fcc

if __name__ == "__main__":
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hoepl", [
                                   "help", "output=", "exec", "pattern=", "loop", "hostname=", "ncpu="])
    except getopt.GetoptError as err:
        print (str(err))
        usage()
        raise

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
            myHostname = value
        elif o == "--ncpu":
            fcc.set_para('ncpu',  int(value))
        elif o in ("--help"):
            usage()
            sys.exit(0)
        else:
            usage()
            sys.exit(1)

    crawl_pattern = ['S:', 'p[0-9]+',
                     'Proteomics',
                     '(EXTERNAL|LUMOS|FUSION|QEXACTIVE|QEXACTIVEHF|QEXACTIVEHFX)_[1-9]',
                     '[a-z]{3,18}_[0-9]{8}(_[-a-zA-Z0-9_]+){0,1}',
                     '[-a-zA-Z0-9_]+.(RAW|raw)$']

    fcc.set_para('crawl_pattern', crawl_pattern)
    fcc.set_para('max_time_diff', 60 * 60 * 24 * 7 * 2) 
    fcc.run()
