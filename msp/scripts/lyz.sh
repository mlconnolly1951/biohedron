#! /bin/sh
../bin/msroll -m ../pdb/2lyz.pdb -q lyz.pqms -f lyz.sel
../bin/msdraw -i lyz.script -r lyz.rgb -z 0.5
