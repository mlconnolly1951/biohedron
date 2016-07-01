/*
DS

DOT SURFACE

DOT SURFACE COMPUTER PROGRAM WRITTEN IN THE C LANGUAGE

Copyright 1986 by Michael L. Connolly
All Rights Reserved

WRITTEN BY MICHAEL L. CONNOLLY

OWNER OF BIOHEDRON

APRIL 23, 1986

Revised:  February 8, 2000

*/

/* DS SUBROUTINE PACKAGE */

#include "ms.h"
#include "msenum.h"
#include "mslimit.h"
#include "msdefault.h"
#include "msstruct.h"
#include "msglobal.h"
#include "msfunction.h"

/* global variables for ds */


/* leave these two strings here so that
   there is a copyright statement in the object code */

char *cpyrit = "Copyright 1986 by Michael L. Connolly";
char *arr = "All Rights Reserved";


/* highest level function in DS subroutine package:
                   calculates dot surface of molecule */

int ds (struct dsdesc  *dsd)
{

    /* determine number of entries in table */
    determine_ntab ();
    if (dsd == NULL) {
        set_error1 ("null pointer passed to ds");
        return (550);
    }

 /* initialization and transfering records to global variables */

    dsdini (dsd);
    dsdget (dsd);
    dsgetir (dsd -> intern);

    dsstage ("input checking");
    dsinch ();
    if (errflg)
        goto error;

 /* allocate memory for problem array */
    problem = allocate_shorts (natom);
    if (problem == NULL) {
        dserr (500, "ds");
        goto error;
    }

    dsstage ("neighbor list creation");
    dslocal ();
    if (errflg)
        goto error;
    if (connected && natom < 3)
        connected = 0;            /* automatic reset */

    dsstage ("torus creation");
    dsgtor ();
    if (errflg)
        goto error;
    if (connected && ntori <= 0)
        connected = 0;            /* automatic reset */

    dsstage ("probe placement");
    dsplace ();
    if (errflg)
        goto error;

    /* check for trouble in starting connected rolling */

    if (connected && nprobe <= 0) {
        connected = 0;            /* automatic reset */
        if (ntori > 0) {
        /* get rid of torus */
            if (ntori != 1) {
                dserr (855, "ds");
                goto error;
            }
            dsfcir (hedtor -> central);
            if (errflg)
                goto error;
            free_object (TORUS, (short *) hedtor);
            ntori = 0;
            hedtor = NULL;
            taltor = NULL;
        }
    /* try again */
        dsstage ("torus creation (2nd try)");
        dsgtor ();
        if (errflg)
            goto error;
        dsstage ("probe placement (2nd try)");
        dsplace ();
        if (errflg)
            goto error;
    }

    if (connected) {
        dsstage ("connected rolling");
        dsplin ();
        if (errflg)
            goto error;
        dsstep (hedprb);
    }
    if (errflg)
        goto error;

    dsstage ("concave surface generation");
    dsprbs ();
    if (errflg)
        goto error;

    dsstage ("convex and toroidal surface generation");
    dsall ();
    if (errflg)
        goto error;

    dsstage ("cleanup");
    dsclean ();
    if (errflg)
        goto error;

    dsstage ("surface point count");
    dscount ();
    if (errflg)
        goto error;

error:

/* store records from global variables */

    dsputir (dsd -> intern);
    dsdput (dsd);
    return (errflg);
}


/* FUNCTIONS CALLED BY DS, IN ALPHABETICAL ORDER */

/* see: dsa.c, dsb.c, dsc.c, dsd.c, dse.c, dsf.c, dsg.c, dsh.c, dsi.c,
        dsl.c, dsm.c, dsn.c, dso.c, dsp.c, dsr.c, dss.c, dst.c, dsu.c */

/* End of DS subroutine package.

Copyright 1986 by Michael L. Connolly
All Rights Reserved

*/
