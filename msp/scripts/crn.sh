#! /bin/sh
../bin/msdraw -i bas.script -r bas.bmp bmp
../bin/msroll -m ../pdb/1crn.pdb -q crn.pqms
../bin/msdraw -i ms.script -r ms.sun sun
../bin/msdraw -i lucent.script -r lucent.sgi

