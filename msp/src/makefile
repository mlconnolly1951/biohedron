# Molecular Surface Package
# Copyright 1986 by Michael L. Connolly
# All Rights Reserved
#
# February 17, 2006
# Version 3.9.3

CC = cc
#uncomment the gcc and comment the cc for using gcc for compiling
#CC = gcc

LD = cc
#uncomment the gcc and comment the cc for using gcc for linking
#LD = gcc

CFLAGS = -c
INCLUDES = ms.h msenum.h mslimit.h msdefault.h msstruct.h msglobal.h msfunction.h
BIN = ../bin/

all: msroll msdraw msform mstran

# General routines

ms.o: ms.c ${INCLUDES}
	${CC} ${CFLAGS} ms.c
msmath.o: msmath.c ${INCLUDES}
	${CC} ${CFLAGS} msmath.c
msmem.o: msmem.c ${INCLUDES}
	${CC} ${CFLAGS} msmem.c
msdata.o: msdata.c ${INCLUDES}
	${CC} ${CFLAGS} msdata.c
msarray.o: msarray.c ${INCLUDES}
	${CC} ${CFLAGS} msarray.c
msrecord.o: msrecord.c ${INCLUDES}
	${CC} ${CFLAGS} msrecord.c
mscept.o: mscept.c ${INCLUDES}
	${CC} ${CFLAGS} mscept.c
	
# SGI Open Inventor 3D file format

msinventor.o: msinventor.c ${INCLUDES}
	${CC} ${CFLAGS} msinventor.c

# Vector Output
msvector.o: msvector.c ${INCLUDES}
	${CC} ${CFLAGS} msvector.c

# Image Output
msimage.o: msimage.c ${INCLUDES}
	${CC} ${CFLAGS} msimage.c

# Molecule

msmolecule.o: msmolecule.c ${INCLUDES}
	${CC} ${CFLAGS} msmolecule.c
msatom.o: msatom.c ${INCLUDES}
	${CC} ${CFLAGS} msatom.c
mspdb.o: mspdb.c ${INCLUDES}
	${CC} ${CFLAGS} mspdb.c

# Coloring and Opacity

mscolor.o: mscolor.c ${INCLUDES}
	${CC} ${CFLAGS} mscolor.c

# Piecewise Quartic surface

msquartic.o: msquartic.c ${INCLUDES}
	${CC} ${CFLAGS} msquartic.c
mscircle.o: mscircle.c ${INCLUDES}
	${CC} ${CFLAGS} mscircle.c
mscycle.o: mscycle.c ${INCLUDES}
	${CC} ${CFLAGS} mscycle.c
msinside.o: msinside.c ${INCLUDES}
	${CC} ${CFLAGS} msinside.c


# Areas and Volumes

mschunk.o: mschunk.c ${INCLUDES}
	${CC} ${CFLAGS} mschunk.c
msvolume.o: msvolume.c ${INCLUDES}
	${CC} ${CFLAGS} msvolume.c

# Solid Angle (omega)
	
mssolid.o: mssolid.c ${INCLUDES}
	${CC} ${CFLAGS} mssolid.c
mscritical.o: mscritical.c ${INCLUDES}
	${CC} ${CFLAGS} mscritical.c


# Polyhedra

msphn.o: msphn.c ${INCLUDES}
	${CC} ${CFLAGS} msphn.c
msphnio.o: msphnio.c ${INCLUDES}
	${CC} ${CFLAGS} msphnio.c
mshedron.o: mshedron.c ${INCLUDES}
	${CC} ${CFLAGS} mshedron.c
mstess.o: mstess.c ${INCLUDES}
	${CC} ${CFLAGS} mstess.c
mscontour.o: mscontour.c ${INCLUDES}
	${CC} ${CFLAGS} mscontour.c
msalpha.o: msalpha.c ${INCLUDES}
	${CC} ${CFLAGS} msalpha.c
msbisect.o: msbisect.c ${INCLUDES}
	${CC} ${CFLAGS} msbisect.c
msdegen.o: msdegen.c ${INCLUDES}
	${CC} ${CFLAGS} msdegen.c
mspurge.o: mspurge.c ${INCLUDES}
	${CC} ${CFLAGS} mspurge.c
mssimple.o: mssimple.c ${INCLUDES}
	${CC} ${CFLAGS} mssimple.c
mstent.o: mstent.c ${INCLUDES}
	${CC} ${CFLAGS} mstent.c
msedger.o: msedger.c ${INCLUDES}
	${CC} ${CFLAGS} msedger.c

# Dot Surface

msdot1.o: msdot1.c ${INCLUDES}
	${CC} ${CFLAGS} msdot1.c
msdot2.o: msdot2.c ${INCLUDES}
	${CC} ${CFLAGS} msdot2.c
msdot3.o: msdot3.c ${INCLUDES}
	${CC} ${CFLAGS} msdot3.c
msbox.o: msbox.c ${INCLUDES}
	${CC} ${CFLAGS} msbox.c

# Density

msdenio.o: msdenio.c ${INCLUDES}
	${CC} ${CFLAGS} msdenio.c
msdensity.o: msdensity.c ${INCLUDES}
	${CC} ${CFLAGS} msdensity.c
msphnden.o: msphnden.c ${INCLUDES}
	${CC} ${CFLAGS} msphnden.c
mscube.o: mscube.c ${INCLUDES}
	${CC} ${CFLAGS} mscube.c
msjoin.o: msjoin.c ${INCLUDES}
	${CC} ${CFLAGS} msjoin.c
msborder.o: msborder.c ${INCLUDES}
	${CC} ${CFLAGS} msborder.c
msshape.o: msshape.c ${INCLUDES}
	${CC} ${CFLAGS} msshape.c
msback.o: msback.c ${INCLUDES}
	${CC} ${CFLAGS} msback.c
msoringe.o: msoringe.c ${INCLUDES}
	${CC} ${CFLAGS} msoringe.c


# Rendering

msrender.o: msrender.c ${INCLUDES}
	${CC} ${CFLAGS} msrender.c
msbuffer.o: msbuffer.c ${INCLUDES}
	${CC} ${CFLAGS} msbuffer.c
msmodel.o: msmodel.c ${INCLUDES}
	${CC} ${CFLAGS} msmodel.c

# Plotting

msplot.o: msplot.c ${INCLUDES}
	${CC} ${CFLAGS} msplot.c

msplotio.o: msplotio.c ${INCLUDES}
	${CC} ${CFLAGS} msplotio.c

# Protein Shape Library
	
msset.o: msset.c ${INCLUDES}
	${CC} ${CFLAGS} msset.c


# MSRoll

msroll0.o: msroll0.c ${INCLUDES}
	${CC} ${CFLAGS} msroll0.c
msroll.o: msroll.c ${INCLUDES}
	${CC} ${CFLAGS} msroll.c
msstream.o: msstream.c ${INCLUDES}
	${CC} ${CFLAGS} msstream.c
msgrid.o: msgrid.c ${INCLUDES}
	${CC} ${CFLAGS} msgrid.c
mstorus.o: mstorus.c ${INCLUDES}
	${CC} ${CFLAGS} mstorus.c
msquartet.o: msquartet.c ${INCLUDES}
	${CC} ${CFLAGS} msquartet.c
msface.o: msface.c ${INCLUDES}
	${CC} ${CFLAGS} msface.c
mscusp.o: mscusp.c ${INCLUDES}
	${CC} ${CFLAGS} mscusp.c
mscomponent.o: mscomponent.c ${INCLUDES}
	${CC} ${CFLAGS} mscomponent.c
msconvert.o: msconvert.c ${INCLUDES}
	${CC} ${CFLAGS} msconvert.c
msrollio.o: msrollio.c ${INCLUDES}
	${CC} ${CFLAGS} msrollio.c

msroll: msroll0.o msmath.o ms.o mscept.o \
	msroll.o msquartic.o msrollio.o mscircle.o \
	mscycle.o msstream.o msmolecule.o mscusp.o msquartet.o \
	mstorus.o msface.o msvolume.o msconvert.o msmem.o mstess.o \
	msalpha.o msbisect.o msdegen.o mspurge.o mssimple.o mstent.o \
	msinside.o mshedron.o msphn.o msphnio.o msgrid.o \
	msdot1.o msdot2.o msdot3.o msbox.o mscomponent.o \
	mschunk.o mscolor.o msset.o msatom.o mspdb.o msarray.o
	${LD} -o ${BIN}msroll \
	msroll0.o msmath.o ms.o mscept.o \
	msroll.o msquartic.o msrollio.o mscircle.o \
	mscycle.o msstream.o msmolecule.o mscusp.o msquartet.o \
	mstorus.o msface.o msvolume.o msconvert.o msmem.o mstess.o \
	msalpha.o msbisect.o msdegen.o mspurge.o mssimple.o mstent.o \
	msinside.o mshedron.o msphn.o msphnio.o msgrid.o \
	msdot1.o msdot2.o msdot3.o msbox.o mscomponent.o \
	mschunk.o mscolor.o msset.o msatom.o mspdb.o msarray.o -lm


# MSDraw

msdraw0.o: msdraw0.c ${INCLUDES}
	${CC} ${CFLAGS} msdraw0.c
msdraw.o: msdraw.c ${INCLUDES}
	${CC} ${CFLAGS} msdraw.c
msdrawio.o: msdrawio.c ${INCLUDES}
	${CC} ${CFLAGS} msdrawio.c
msscript.o: msscript.c ${INCLUDES}
	${CC} ${CFLAGS} msscript.c

msdraw: msdraw0.o msdraw.o msdrawio.o ms.o mscept.o msmath.o \
	msmolecule.o msphn.o msedger.o msphnio.o mscontour.o \
	msmem.o msset.o msatom.o mspdb.o msarray.o \
	msrender.o msbuffer.o msplot.o msplotio.o msinside.o msinventor.o \
	msimage.o msbox.o msquartic.o mschunk.o \
	mscolor.o msmodel.o mscircle.o msscript.o
	${LD} -o ${BIN}msdraw msdraw0.o msdraw.o msdrawio.o \
	ms.o mscept.o msmath.o msphn.o msedger.o msphnio.o mscontour.o \
	msmem.o msset.o msatom.o mspdb.o msarray.o \
	msrender.o msbuffer.o msplot.o msplotio.o msimage.o msquartic.o mschunk.o \
	msinside.o msinventor.o msbox.o msmolecule.o \
	mscolor.o msmodel.o mscircle.o msscript.o -lm


# MSForm

msform0.o: msform0.c ${INCLUDES}
	${CC} ${CFLAGS} msform0.c

msform.o: msform.c ${INCLUDES}
	${CC} ${CFLAGS} msform.c

msform: msform0.o msform.o ms.o mscept.o msmath.o msdata.o \
	mscircle.o msphn.o msedger.o msphnio.o \
	msmem.o msinside.o msinventor.o mscycle.o mssolid.o mscritical.o \
	msoringe.o msrecord.o mscolor.o msquartic.o mschunk.o \
	msdensity.o msdenio.o msphnden.o mscube.o msjoin.o \
	msborder.o msshape.o msback.o
	${LD} -o ${BIN}msform msform0.o msform.o ms.o mscept.o msmath.o \
	msdata.o msphn.o msedger.o msphnio.o msmem.o mscircle.o \
	msinside.o msinventor.o mscycle.o mssolid.o mscritical.o msoringe.o \
	mscolor.o msquartic.o mschunk.o \
	msrecord.o msdensity.o msdenio.o msphnden.o mscube.o \
	msjoin.o msborder.o msshape.o msback.o -lm

# MSTran

mstran0.o: mstran0.c ${INCLUDES}
	${CC} ${CFLAGS} mstran0.c
mstran.o: mstran.c ${INCLUDES}
	${CC} ${CFLAGS} mstran.c

mstran: mstran0.o mstran.o ms.o mscept.o \
	msmath.o mscolor.o msphn.o msedger.o \
	msphnio.o mscontour.o msquartic.o msatom.o mspdb.o msarray.o \
	msmolecule.o msbox.o msmem.o msset.o msrecord.o msinventor.o \
	msoringe.o msdenio.o msdensity.o msvector.o mschunk.o mscircle.o
	${LD} -o ${BIN}mstran mstran0.o mstran.o ms.o mscept.o msmath.o \
	mscolor.o msquartic.o msphn.o msedger.o msphnio.o mscontour.o \
	msatom.o mspdb.o msarray.o msmolecule.o msbox.o msmem.o msset.o \
	msrecord.o msinventor.o msoringe.o msdenio.o msdensity.o \
	msvector.o mschunk.o mscircle.o -lm


nomsdraw:
	rm msdraw0.? msdraw.? msdrawio.? msscript.? msplot.? msplotio.? msmodel.? msrender.? msbuffer.?
	
nomsform:
	rm msform0.? msform.? msdensity.? msphnden?.? mssolid.? msoringe.?
	

clean:
	/bin/rm -f *.o
