#! /bin/sh
../bin/msroll -f ptci.sel -m ../pdb/2ptc.pdb -q ptci.pqms -n ptci
../bin/msroll -f ptce.sel -m ../pdb/2ptc.pdb -q ptce.pqms -n ptce
../bin/msdraw -i ptc.script -b 18.0 -r ptc.sun sun -z -0.05 -s 600
