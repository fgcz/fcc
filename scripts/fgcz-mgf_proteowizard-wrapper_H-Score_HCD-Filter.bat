echo %1
echo %2


"C:\Program Files\ProteoWizard\ProteoWizard 3.0.3860\msconvert.exe"  -o - --mgf --filter "peakPicking true [2,3]" --filter "activation HCD" %1 > %2.tmp
python C:\FGCZ\fcc\scripts\Hscorer_single_file.py --massfile C:\FGCZ\fcc\mascotmasses.txt --tolerance 50 --inp %2.tmp --outp %2
del %2.tmp

exit 1

