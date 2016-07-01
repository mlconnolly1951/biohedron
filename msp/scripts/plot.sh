#! /bin/sh
../bin/msroll -m ../pdb/4pti.pdb -f plot.sel -t pti.vet
../bin/msdraw -i plot.script -p plot.ps ps
