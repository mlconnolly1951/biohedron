#! /bin/sh
../bin/msroll -m ../pdb/1crb.pdb -f crb.sel -q crb.pqms
../bin/msroll -m ../pdb/1crb.pdb -p 0.0 -f retinol.sel -q retinol.pqms
../bin/msdraw -i crb.script -z 0.0 -r crb.sgi

