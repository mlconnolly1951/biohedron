/*
DS

DOT SURFACE

DOT SURFACE COMPUTER PROGRAM WRITTEN IN THE C LANGUAGE

Copyright 1986 by Michael L. Connolly
All Rights Reserved

WRITTEN BY MICHAEL L. CONNOLLY

Revised:  March 6, 2000

*/

/* DS SUBROUTINE PACKAGE */

#include "ms.h"
#include "msenum.h"
#include "mslimit.h"
#include "msdefault.h"
#include "msstruct.h"
#include "msglobal.h"
#include "msfunction.h"



void dsaac (struct circle *cir, int ia, double ang1, double ang2)     /* And Arc with Circle */
{
    struct endpnt  *ept1, *ept2, *ept3, *ept4;

    atom = ia;
    if (cir -> head == NULL) {
    /* start new endpoint list */
        ept1 = dsnwept (1, 0.0, ia);
        if (errflg)
            return;
        ept2 = dsnwept (0, ang2, ia);
        if (errflg)
            return;
        cir -> head = ept1;
        ept1 -> next = ept2;
        return;
    }

 /* we already have an endpoint list */

 /* AND with current accessibility list */
    if (ang2 >= ang1) {
        for (ept1 = cir -> head; ept1 != NULL; ept1 = ept2 -> next) {
            ept2 = ept1 -> next;
            if (!ept1 -> begin) {
                dserr (610, "dsaac");
                return;
            }

            if (ept1 -> angle >= ept2 -> angle) {
                dserr (611, "dsaac");
                return;
            }

            if (ang2 < ept1 -> angle || ang1 > ept2 -> angle) {
            /* no overlap, mark to delete */
                ept1 -> angle = -1.0;
                ept2 -> angle = -1.0;
            }
            else {
            /* just like the real line */
                if (ang1 > ept1 -> angle) {
                    ept1 -> angle = ang1;
                    ept1 -> atom = (unsigned short) ia;
                }
                if (ang2 < ept2 -> angle) {
                    ept2 -> angle = ang2;
                    ept2 -> atom = (unsigned short) ia;
                }
            }
        }                       /* end of for loop */
    }                           /* end of ang2 >= ang1 */
    else {
    /* ang2 < ang1 */
        for (ept1 = cir -> head; ept1 != NULL; ept1 = ept2 -> next) {
            if (!ept1 -> begin) {
                dserr (610, "dsaac");
                return;
            }
            ept2 = ept1 -> next;
            if (ept1 -> angle >= ept2 -> angle) {
                dserr (611, "dsaac");
                return;
            }
            if (ang2 < ept1 -> angle && ang1 > ept2 -> angle) {
            /* no overlap, mark to delete */
                ept1 -> angle = -1.0;
                ept2 -> angle = -1.0;
            }
            else
                if (ang2 >= ept1 -> angle && ang1 > ept2 -> angle) {
                /* connected overlap at beginning */
                    if (ang2 < ept2 -> angle) {
                        ept2 -> angle = ang2;
                        ept2 -> atom = (unsigned short) ia;
                    }
                }
                else
                    if (ang2 < ept1 -> angle && ang1 <= ept2 -> angle) {
                    /* connected overlap at end */
                        if (ang1 > ept1 -> angle) {
                            ept1 -> angle = ang1;
                            ept1 -> atom = (unsigned short) ia;
                        }
                    }
                    else {      /* disconnected overlap */
                    /* need new endpoints */
                        ept3 = dsnwept (0, ang2, ia);
                        if (errflg)
                            return;
                        ept4 = dsnwept (1, ang1, ia);
                        if (errflg)
                            return;
                        ept1 -> next = ept3;
                        ept3 -> next = ept4;
                        ept4 -> next = ept2;
                    }
        }                       /* end of for loop */
    }                           /* end of ang2 < ang1 */

    dsccir (cir);               /* clean circle of deleted endpoints */

    return;
}

void dsall ()                        /* gather or create all shapes of surface */
{
    int    k;
    int     a, which, ntorsp, ncavsp, ncvxsp;
    struct cluster *clu, *headclu, *tailclu;
    struct cluster *reenclu, *cavclu, *torclu;
    struct torus   *tor;
    struct tlist   *tl, *tlnext;
    struct tlist  **tlhead, **tltail;

    atom = EOL;                 /* no relevant atom number if error now */

 /* form atom -> torus lists */

    tlhead = (struct tlist **) allocate_pointers (TLIST, natom);
    if (tlhead == NULL) {
        dserr (500, "dsall");
        return;
    }

    tltail = (struct tlist **) allocate_pointers (TLIST, natom);
    if (tltail == NULL) {
        dserr (500, "dsall");
        return;
    }

    for (tor = hedtor; tor != NULL; tor = tor -> next)
        for (k = 0; k < 2; k++) {
            a = tor -> atom[k];
            if (*(atmatt + a) < AREA)
                continue;
            tl = (struct tlist *) allocate_object (TLIST);
            if (tl == NULL) {
                dserr (500, "dsall");
                return;
            }
        /* link into appropriate list */
            if (*(tlhead + a) == NULL)
                *(tlhead + a) = tl;
            else
                (*(tltail + a)) -> next = tl;
            *(tltail + a) = tl;
            tl -> tor = tor;
        }

 /* for each atom, gather convex, toroidal, concave surface points */

    for (a = 0; a < natom; a++) {
        atom = a;               /* store atom number for possible error */
        if (*(atmatt + a) < AREA)
            continue;           /* ignored or blocker */
        cavclu = *(hedcav + a); /* concave surface cluster */
        if (cavclu == NULL)
            ncavsp = 0;
        else
            ncavsp = cavclu -> nmem;
        if (*(atmatt + a) < POINTS && ncavsp > 0) {
            dserr (780, "dsall");
            return;
        }
        tl = *(tlhead + a);
        /* cannot have concave w/o toroidal */
        if (cavclu && tl == NULL && ! *(problem + a)) {
            dserr (782, "dsall");
            return;
        }
        /* atom with neighbors but no accessible tori must be buried */
        if ((*(hedrun + a) != NULL) && (tl == NULL) && ! *(problem + a)) {
            *(cvxsp + a) = NULL;
            *(rensp + a) = NULL;
            continue;
        }
        *(cvxsp + a) = dscvx (a, tl);   /* calculate convex surface */
        if (errflg)
            return;
        if (*(cvxsp + a) == NULL)
            ncvxsp = 0;
        else
            ncvxsp = (*(cvxsp + a)) -> nmem;
        if (*(atmatt + a) < POINTS && ncvxsp > 0) {
            dserr (780, "dsall");
            return;
        }
    /* reentrant */
        headclu = NULL;         /* initialize linked list of clusters */
        tailclu = NULL;
        for (tl = *(tlhead + a); tl != NULL; tl = tl -> next) {
            tor = tl -> tor;
            which = (a == tor -> atom[1]);      /* flag for which side */
            clu = dstor (tor, which);           /* compute toroidal srf */
            if (errflg)
                return;
            if (clu == NULL)
                continue;
            if (headclu == NULL)        /* link into list */
                headclu = clu;
            else
                tailclu -> next = clu;
            tailclu = clu;
        }

        if (headclu == NULL)
            torclu = NULL;
        else
            torclu = dsmerge (headclu); /* merge list into one cluster */
        if (errflg)
            return;

        if (torclu == NULL)
            ntorsp = 0;
        else
            ntorsp = torclu -> nmem;

        if (*(atmatt + a) < POINTS && ntorsp > 0) {
            dserr (780, "dsall");
            return;
        }

        /* combine toroidal and concave surface point clusters */
        if (torclu == NULL && cavclu == NULL)
            reenclu = NULL;
        else
            if (torclu != NULL && cavclu == NULL)
                reenclu = torclu;
            else
                if (torclu == NULL && cavclu != NULL)
                    reenclu = cavclu;
                else {
                /* merge toroidal and concave */
                    torclu -> next = cavclu;
                    reenclu = dsmerge (torclu);
                    if (errflg)
                        return;
                }
        *(rensp + a) = reenclu; /* store pointer to reentrant cluster */

    /* reentrant areas are unreliable (possibly 0 or low) if the atom
       has an att number < POINTS */

        if (*(atmatt + a) < POINTS && reenclu != NULL)
            reenclu -> area = 0.0;

    /* free tl memory */
        for (tl = *(tlhead + a), tlnext = NULL; tl != NULL; tl = tlnext) {
            tlnext = tl -> next;
            free_object (TLIST, (short *) tl);
        }
    }                           /* end of atom loop */
    free_pointers (TLIST, (void *) tlhead);
    free_pointers (TLIST, (void *) tltail);
}

struct run *dsand (struct run *f1, struct run *f2) /* intersection (and) of two atom lists */
{
    int     b1, b2, e1, e2, b3, e3;     /* begin and end of ranges */
    struct run *r1, *r2, *p, *r, *f;

    if (f1 == NULL || f2 == NULL)
        return (NULL);

    /* initialization */
    r1 = f1;
    r2 = f2;
    f = NULL;
    r = NULL;

    while (r1 != NULL && r2 != NULL) {
        b1 = r1 -> atom[0];
        b2 = r2 -> atom[0];
        e1 = r1 -> atom[1];
        e2 = r2 -> atom[1];
        if (e1 < b2)
            r1 = r1 -> next;
        else
            if (e2 < b1)
                r2 = r2 -> next;
            else {              /* overlap between ranges of two runs */
                b3 = b1 > b2 ? b1 : b2;
                e3 = e1 < e2 ? e1 : e2;
                p = r;
                /* new run */
                r = (struct run *) allocate_object (RUN);
                if (r == NULL) {
                    dserr (500, "dsand");
                    return (NULL);
                }
                /* link into function output list */
                if (p == NULL)
                    f = r;
                else
                    p -> next = r;

                /* store begin and end of new range */
                r -> atom[0] = (unsigned short) b3;
                r -> atom[1] = (unsigned short) e3;
                if (e1 <= e2)
                    r1 = r1 -> next;
                if (e2 <= e1)
                    r2 = r2 -> next;
            }
    }
    return (f);                 /* return pointer to first run of new list */
}

int dsans (atomnum *list)                    /* check whether atom numbers are sorted */
{
    int    n, i;

    if (list == NULL)
        return (1);
    /* determine length of list */
    for (n = 0; *(list + n) != EOL; n++);
    if (n <= 1)
        return (1);

    for (i = 0; i < n - 1; i++)
        if (*(list + i) > *(list + i + 1))
            return (0);         /* inverted order */

    return (1);
}

/* compute arbitrary perpendicular unit vector to given unit vector */
void dsarbp (double   given[3], double perp[3])
{
    int    k;
    double  dott;
    double   try[3];

    /* my own weird algorithm */
    try[0] = given[1] * given[1] + given[2] * given[2];
    try[1] = given[0] * given[0] + given[2] * given[2];
    try[2] = given[0] * given[0] + given[1] * given[1];

    if (!dsnize (try)) {
        dserr (800, "dsarbp");
        return;
    }

    /* subtract parallel component of first try */
    dott = dsdot (given, try);
    for (k = 0; k < 3; k++)
        try[k] = try[k] - dott * given[k];

    if (!dsnize (try)) {
        dserr (800, "dsarbp");
        return;
    }

 /* check correctness */
    if (fabs (dsdot (given, try)) > EPSILON) {
        dserr (810, "dsarbp");
        return;
    }

    for (k = 0; k < 3; k++)
        perp[k] = try[k];
}



/* clean circle of endpoints marked for deletion */
void dsccir (struct circle  *cir)
{
    struct endpnt  *ept1, *prev, *next;

    if (cir == NULL) {
        dserr (760, "dsccir");
        return;
    }
    if (cir->head == NULL) return;
    prev = NULL;
    for (ept1 = cir -> head, next = NULL; ept1 != NULL; ept1 = next) {
        next = ept1 -> next;
        if (ept1 -> angle >= 0.0) {     /* keep this one */
            prev = ept1;
            continue;
        }
        /* fix up pointers for deleted endpoint */
        if (prev == NULL)
            cir -> head = next;
        else
            prev -> next = next;
    /* free memory for endpoint */
        free_object (ENDPNT, (short *) ept1);
    }
    /* mark circle inaccessible if no endpoints left */
    if (cir -> head == NULL)
        cir -> gone = 1;
}

/* count the total number of surface points and sum the areas */
void dscount ()
{
    int     a;
    char message[MAXLINE];
    struct cluster *clu;

    /* the total variables have already been initialized elsewhere */

    for (a = 0; a < natom; a++) {
        clu = *(cvxsp + a);
        if (clu != NULL) {
            nsrfpnt += clu -> nmem;
            cvxarea += clu -> area;
        }
        clu = *(rensp + a);
        if (clu != NULL) {
            nsrfpnt += clu -> nmem;
            renarea += clu -> area;
        }
    }

    if (verbose) {
        sprintf(message, "%8d tori",ntori);
		inform (message);
        sprintf(message, "%8d probes",nprobe);
		inform (message);
        sprintf(message, "%8d low probes",nlow);
		inform (message);
        if (connected) {
            sprintf(message, "%8d dead ends",ndead);
            inform (message);
        }
    }
}

/* circle - plane intersection  (used for concave surface) */
void dscpi (struct circle *cir, double pcen[3], double paxis[3])
{
    int    k;
    double   crad, l, h, ang1, ang2, d2, r, dot1, dot2, dot3;
    double   caxis[3], ccen[3], lvect[3], cpvect[3], svect[3];
    double   mco[3], end1[3], end2[3], ve1[3], ve2[3], base[3];

    if (cir == NULL) {
        dserr (760, "dscpi");
        return;
    }
    if (cir -> gone)
        return;

    atom = cir -> atom;

 /* transfer circle data to local arrays */
    for (k = 0; k < 3; k++) {
        ccen[k] = cir -> center[k];
        caxis[k] = cir -> axis[k];
    }
    crad = cir -> radius;
    if (crad <= 0.0) {
        dserr (859, "dscpi");
        return;
    }

    /* some solid geometry */

    for (k = 0; k < 3; k++)
        cpvect[k] = pcen[k] - ccen[k];
    dot1 = dsdot (cpvect, paxis);
    h = dsdot (cpvect, caxis);
    dscross (paxis, caxis, lvect);
    l = dsnorm (lvect);
    if (l <= 0.0) {
        if (dot1 < 0.0)
            cir -> gone = 1;
        return;
    }

    if (!dsnize (lvect)) {
        dserr (800, "dscpi");
        return;
    }

    dscross (paxis, lvect, svect);

    if (!dsnize (svect)) {
        dserr (800, "dscpi");
        return;
    }

    for (k = 0; k < 3; k++)
        mco[k] = pcen[k] + (h / l) * svect[k];
    d2 = dsdis2 (mco, ccen);
    r = crad * crad - d2;
    if (r <= 0.0) {                     /* no intersection */
        if (dot1 < 0.0)
            cir -> gone = 1;            /* entire circle inaccessible */
        return;
    }

    r = sqrt (r);
    /* set up intersection point coordinates */
    for (k = 0; k < 3; k++) {
        end1[k] = mco[k] - r * lvect[k];
        end2[k] = mco[k] + r * lvect[k];
    }

    /* unit vectors from circle center */
    for (k = 0; k < 3; k++) {
        ve1[k] = (end1[k] - ccen[k]) / crad;
        ve2[k] = (end2[k] - ccen[k]) / crad;
    }

    if (cir -> head == NULL) {          /* first endpoints of circle */
        for (k = 0; k < 3; k++)
            cir -> base[k] = ve1[k];
        ang1 = 0.0;
    }
    else {                              /* not first: base vector exists */
        for (k = 0; k < 3; k++)
            base[k] = cir -> base[k];
        dot2 = dsdot (base, ve1);
        if (dot2 < -1.0 - EPSILON || dot2 > 1.0 + EPSILON) {
            dserr (862, "dscpi");
            return;
        }
        /* don't give math library a bad argument */
        if (dot2 > 1.0)
            dot2 = 1.0;
        else
            if (dot2 < -1.0)
                dot2 = -1.0;
        ang1 = acos (dot2);
        /* inverse cosine function multiple valued: determine quadrant */
        if (dstrip (base, ve1, caxis) <= 0.0)
            ang1 = 2 * PI - ang1;
    }

    /* determine angle between endpoints */
    for (k = 0; k < 3; k++)
        base[k] = cir -> base[k];
    dot3 = dsdot (base, ve2);
    if (dot3 < -1.0 - EPSILON || dot3 > 1.0 + EPSILON) {
        dserr (862, "dscpi");
        return;
    }
    if (dot3 > 1.0)
        dot3 = 1.0;
    else
        if (dot3 < -1.0)
            dot3 = -1.0;
    ang2 = acos (dot3);
    if (dstrip (base, ve2, caxis) <= 0.0)
        ang2 = 2 * PI - ang2;

    if (ang2 <= 0.0)
        return;

 /* intersection exists */

    /* call function to handle linked list of endpoints */
    dsaac (cir, EOL, ang1, ang2);
    if (errflg)
        return;
}

/* cross product of two vectors */
void dscross (double a[3], double b[3], double c[3])
{
    c[0] = a[1] * b[2] - a[2] * b[1];
    c[1] = a[2] * b[0] - a[0] * b[2];
    c[2] = a[0] * b[1] - a[1] * b[0];
}

void dscsi (struct circle *cir, int ia)                 /* circle-sphere intersection */
{
    int    k;
    double   rad1, rad2, rad3, rad4, plus13, minus13;
    double   ang1, ang2, d13, root1, root2;
    double   dot1, dot2, dot3;
    double   cen1[3], cen2[3], cen3[3], cen4[3], axis[3];
    double   u2[3], v1[3], v2[3], v3[3], v4[3], v5[3];

    if (cir -> gone)
        return;
    atom = ia;

 /* transfer circle data to local arrays */
    for (k = 0; k < 3; k++) {
        cen1[k] = cir -> center[k];
        axis[k] = cir -> axis[k];
    }
    rad1 = cir -> radius;
    if (rad1 <= 0.0) {
        dserr (859, "dscsi");
        return;
    }

 /* transfer atom data to local array */
    rad2 = *(atmrad + ia) + pradius;
    for (k = 0; k < 3; k++)
        cen2[k] = *(atmco + 3 * ia + k);

 /* calculate plane-sphere intersection */
    for (k = 0; k < 3; k++)
        v1[k] = cen2[k] - cen1[k];

    dot1 = dsdot (axis, v1);

    rad3 = rad2 * rad2 - dot1 * dot1;
    if (rad3 <= 0.0)
        return;                 /* expanded atom does not intersect plane */

    rad3 = sqrt (rad3);
    plus13 = rad1 + rad3;
    minus13 = rad1 - rad3;
    for (k = 0; k < 3; k++)
        cen3[k] = cen2[k] - dot1 * axis[k];
    for (k = 0; k < 3; k++)
        v2[k] = cen3[k] - cen1[k];
    d13 = dsnorm (v2);
    if (d13 >= plus13)
        return;                 /* expanded atom does not intersect circle */

    root1 = plus13 * plus13 - d13 * d13;
    if (root1 <= 0.0)
        return;

    root1 = sqrt (root1);
    root2 = d13 * d13 - minus13 * minus13;

    if (root2 <= 0.0) {
        if (rad1 >= rad3)
            return;             /* expanded atom inside circle */
        cir -> gone = 1;        /* expanded atom contains circle */
        dsnoend (cir);          /* delete endpoints of circle */
        if (errflg)
            return;
        return;
    }

 /* we have an intersection */
    root2 = sqrt (root2);
    rad4 = 0.5 * root1 * root2 / d13;
    for (k = 0; k < 3; k++)
        u2[k] = v2[k] / d13;
 /* center of s0 (two intersection points) */
    for (k = 0; k < 3; k++)
        cen4[k] = 0.5 * (cen1[k] + cen3[k]) +
            0.5 * u2[k] * (rad1 * rad1 - rad3 * rad3) / d13;
    dscross (u2, axis, v3);

 /* calculate endpoints */
 /* the sign ! the sign ! */
    for (k = 0; k < 3; k++) {
        v4[k] = (cen4[k] - rad4 * v3[k] - cen1[k]) / rad1;
        v5[k] = (cen4[k] + rad4 * v3[k] - cen1[k]) / rad1;
    }
    if (cir -> head == NULL) {          /* first endpoints of circle */
        for (k = 0; k < 3; k++)
            cir -> base[k] = v4[k];
        ang1 = 0.0;
    }
    else {                              /* circle only partially accessible */
        dot2 = dsdot (cir -> base, v4);
        if (dot2 < -1.0 - EPSILON || dot2 > 1.0 + EPSILON) {
            dserr (862, "dscsi");
            return;
        }
        if (dot2 > 1.0)
            dot2 = 1.0;
        else
            if (dot2 < -1.0)
                dot2 = -1.0;
        ang1 = acos (dot2);
        if (dstrip (cir -> base, v4, axis) <= 0.0)
            ang1 = 2 * PI - ang1;
    }
    dot3 = dsdot (cir -> base, v5);
    if (dot3 < -1.0 - EPSILON || dot3 > 1.0 + EPSILON) {
        dserr (862, "dscsi");
        return;
    }

    /* be nice to math library */
    if (dot3 > 1.0)
        dot3 = 1.0;
    else
        if (dot3 < -1.0)
            dot3 = -1.0;
    ang2 = acos (dot3);
    /* get the right quadrant */
    if (dstrip (cir -> base, v5, axis) <= 0.0)
        ang2 = 2 * PI - ang2;

    if (ang2 <= 0.0)
        return;

 /* intersection exists */

    /* do the algebra now that the geometry is done */
    dsaac (cir, ia, ang1, ang2);
    if (errflg)
        return;
}

struct cluster *dscvx (int ia, struct tlist *ftl)         /* convex surface points */
{
    int    k;
    int     nlat, ilat, nnt, nnbr, onecyc;
    int     pass, all, n, i, ja, jnb, kn;
    int     loopb[2], loopi[2];
    unsigned    size;
    double   width, zone, irad, erad, ratio, z, zl, rsmall, rlarge;
    double  *fn, *f, *f1;
    double   vectl[3], ico[3], zsmall[3], zlarge[3], north[3];
    atomnum * nbrs, *jnbr, *nbrtor, *nbrnt;
    struct endpnt   ept1, ept2;
    struct circle   mercir;
    struct run *irun;
    struct circle  *cir;
    struct cluster *merclu, *fclu, *clu, *pclu;
    struct torus   *tor;
    struct tlist   *tl;

    atom = ia;                          /* store atom number of err msg */

    /* move atom data to local variables */

    irad = *(atmrad + ia);              /* radius of atom */
    erad = irad + pradius;              /* radius of expanded sphere */
 /* ratio is used to convert between contact and accessible surfaces */
    ratio = erad / irad;
    density = *(atmden + ia);
    for (k = 0; k < 3; k++)
        ico[k] = *(atmco + 3 * ia + k);

 /* set up lists for neighboring atoms */

    irun = *(hedrun + ia);              /* pointer to first nbr run */
    nbrs = dsupk (irun);                /* unpack neighbors to list */
    if (errflg)
        return (NULL);

    nnbr = dslen (irun);                /* count number of neighbors */
    if (errflg)
        return (NULL);

    nnt = nnbr;                         /* initialize number of non-tori nbrs */
         if (nnbr > 0) {
        nbrtor = allocate_atomnums ((unsigned) nnbr);
        if (nbrtor == NULL) {
            dserr (500, "dscvx");
            return (NULL);
        }
    /* mark neighbors of ia that share tori with ia */
        for (tl = ftl; tl != NULL; tl = tl -> next) {
            tor = tl -> tor;
            /* find atom sharing torus with ia */
            if (ia == tor -> atom[0])
                ja = tor -> atom[1];
            else
                ja = tor -> atom[0];
            /* look for it in the neighbor list */
            for (jnb = 0; jnb < nnbr; jnb++)
                if (*(nbrs + jnb) == ja && *(nbrtor + jnb) == 0) {
                    *(nbrtor + jnb) = 1;        /* mark as torus neighbor */
                    nnt--;                      /* decrement non-tori counter */
                                 }
        }
    }
    else
        nbrtor = NULL;

    if (nnt < 0) {
        dserr (712, "dscvx");
        return (NULL);
    }

 /* set up array of non-torus neighbors */
    if (nnt > 0) {
        nbrnt = allocate_atomnums ((unsigned) (nnt + 1));
        if (nbrnt == NULL) {
            dserr (500, "dscvx");
            return (NULL);
        }
        kn = 0;
        for (jnb = 0; jnb < nnbr; jnb++)
            if (!*(nbrtor + jnb))                       /* non-torus check */
                *(nbrnt + kn++) = *(nbrs + jnb);        /* store atom num */
        if (kn != nnt) {
            dserr (712, "dscvx");
            return (NULL);
        }
        *(nbrnt + nnt) = EOL;                           /* mark end of list */
    }
    else
        nbrnt = NULL;


 /* speed-up code for atoms with one cycle of edges */

    /* this flag is 1 if dscyc finds exactly one cycle and no problems */
    onecyc = dscyc (ia, nbrs, ftl, north, &zone);
    if (errflg)
        return (NULL);

    if (!onecyc) {
        /* north vector points towards positive z */
        for (k = 0; k < 3; k++)
            north[k] = (k == 2) ? 1.0 : 0.0;
        zone = 0.0;     /* ignored */
    }

 /* set up meridian semi-circle from north to south pole */
    for (k = 0; k < 3; k++) {
        mercir.center[k] = ico[k];
        mercir.base[k] = north[k];
    }
    dsarbp (north, mercir.axis);        /* arbitrary perpendicular */
    if (errflg)
        return (NULL);
    mercir.radius = irad;               /* circle radius = atom radius */
    mercir.width = 0.0;
    mercir.gone = 0;
    mercir.atom = (unsigned short) ia;
    mercir.head = &ept1;
    ept1.angle = 0.0;
    ept1.begin = 1;
    ept1.next = &ept2;
    ept2.angle = PI;
    ept2.begin = 0;
    ept2.next = NULL;
    merclu = dsdivid (&mercir);         /* points for latitudes */
    if (errflg)
        return (NULL);
    nlat = merclu -> nmem;
    if (nlat <= 0) {
        dserr (865, "dscvx");
        return (NULL);
    }
    /* initialize latitude cluster linked list */
    pclu = NULL;                        /* previous cluster pointer */
    fclu = NULL;                        /* first cluster pointer */

    /* adjust width for more accurate spherical zone area formula */
    width = 2 * irad * sin (PI / (2 * nlat));
    f1 = merclu -> points;

    loopb[0] = 0;
    loopi[0] = 1;
    loopb[1] = nlat - 1;
    loopi[1] = -1;
    all = 0;

    /* two passes in opposite directions:

                    0: from partially accessible to inaccessible
                    1: from partially accessible to completely accessible */

    for (pass = 0; pass < 2; pass++)
        for (ilat = loopb[pass]; loopi[pass] * (ilat - loopb[pass]) < nlat;
                ilat += loopi[pass]) {
         /* project point on meridian onto axis */
            for (k = 0; k < 3; k++)
                vectl[k] = *(f1 + 3 * ilat + k) - ico[k];
            z = dsdot (north, vectl);
         /* expand for expanded sphere for accessibility check */
            zl = z * ratio;
            if (pass == 0 && zl > zone)         /* ignore northern hemisphere */
                                continue;
            if (pass == 1 && zl <= zone)        /* ignore southern hemisphere */
                                continue;
            rsmall = irad * irad - z * z;
            if (rsmall < 0.0) {
                dserr (861, "dscvx");
                return (NULL);
            }
            rsmall = sqrt (rsmall);             /* contact latitude radius */
            if (rsmall < 0.0) {
                dserr (859, "dscvx");
                return (NULL);
            }
            if (rsmall == 0.0) continue;
            rlarge = rsmall * ratio;            /* accessible latitude radius */
                         for (k = 0; k < 3; k++)
                zsmall[k] = ico[k] + z * north[k];      /* contact center */
            for (k = 0; k < 3; k++)
                zlarge[k] = ico[k] + z * ratio * north[k]; /* access. center */
        /* use large for accessibility */
            cir = dsnwcir (zlarge, rlarge, north);
            if (errflg)
                return (NULL);

        /* cut away inaccessible parts of circle: */
            if (nbrs != NULL && !all) {
            /* loop thru neighbors that share a torus with ia: */
                for (jnbr = nbrs; *jnbr != EOL; jnbr++) {
                    jnb = jnbr - nbrs;
                    ja = *jnbr;
                    if (!*(nbrtor + jnb))
                        continue;
                    dscsi (cir, ja);
                    if (errflg)
                        return (NULL);
                    if (cir -> gone)
                        break;
                }
            /* handle all nbrs that don't share a torus with ia: */
                if (nnt > 0 && !cir -> gone)
                    dsnti (cir, nbrnt);
                if (errflg)
                    return (NULL);
            }
            /* check for circle rendered inaccessible by neighbors */
            if (!cir -> gone) {
            /* at least partially accessible */
            /* convert to surface points */
            /* use contact radius, center */
                cir -> radius = rsmall;
                cir -> width = width;
                cir -> atom = (unsigned short) ia;
                for (k = 0; k < 3; k++)
                    cir -> center[k] = zsmall[k];
                clu = dsdivid (cir);
                if (errflg)
                    return (NULL);
                /* link into latitude cluster list */
                if (clu != NULL) {
                    if (pclu == NULL)
                        fclu = clu;
                    else
                        pclu -> next = clu;
                    pclu = clu;
                }
            }
            /* speed up code for accessible area bounded by 1 cycle */
            if (onecyc) {
                /* check for remainder of latitudes inaccessible */
                if (pass == 0 && zl < zone && cir -> gone) {
                    dsfcir (cir);       /* free circle */
                    if (errflg)
                        return (NULL);
                    break;              /* finished with first pass */
                }
                /* check for remainder of latitudes all accessible */
                if (pass == 1 && zl > zone && !all &&
                        !cir -> gone && cir -> head == NULL)
                    all = 1;            /* turn on flag to skip cuts */
            }
            dsfcir (cir);               /* free circle */
            if (errflg)
                return (NULL);
        }                       /* end of latitude loop */

    /* free temporary memory */
    free_atomnums (nbrs);
    free_atomnums (nbrtor);
    free_atomnums (nbrnt);
    dsfclu (merclu);
    if (errflg)
        return (NULL);

 /* merge latitude clusters */
    if (fclu == NULL)
        return (NULL);
    clu = dsmerge (fclu);
    if (errflg)
        return (NULL);
    n = clu -> nmem;
    clu -> normals = NULL;

    /* create normal vectors if requested */
    if (*(atmatt + ia) >= NORMALS && n > 0) {
        size = n * 3 * sizeof (double);
        clu -> normals = allocate_doubles (n * 3, 0, NORMALS);
        if (clu -> normals == NULL) {
            dserr (500, "dscvx");
            return (NULL);
        }
        fn = clu -> normals;
        f = clu -> points;
        for (i = 0; i < n; i++)
            for (k = 0; k < 3; k++)
                *fn++ = (*f++ - ico[k]) / irad;
    }
    return (clu);
}

/* check for one cycle on atom ia */
int dscyc (int ia, atomnum *nbrs, struct tlist *ftl, double north[3], double *z)
{
    int    k, ie, je;
    int     nedg,  nfree,  ja, nlink, ntrav, a1, a2;
    double   ang, sinang, cosang, dt;
    atomnum * jnbr;
    double   vect1[3], tp[3], ept1co[3], ept2co[3], zone[3];
    double   jvec[3], ico[3], zenith[3];
    struct tlist   *tl;
    struct torus   *tor;
    struct circle  *cir;
    struct endpnt  *ept1, *ept2;
    struct edg *edges, *eptr;

    atom = ia;
    /* transfer coordinates to local array */
    for (k = 0; k < 3; k++)
        ico[k] = *(atmco + 3 * ia + k);


    /* cycle not possible without tori */
    if (ftl == NULL)
        return (0);
    /* initialization */
    nedg = 0;
    nfree = 0;

    /* count number of edges and number of complete (free) circles */
    /* we need these numbers for memory allocation */
    for (tl = ftl; tl != NULL; tl = tl -> next) {
        tor = tl -> tor;
        cir = tor -> central;
        if (cir -> head == NULL) {
            nedg++;
            nfree++;
        }
        else {
            for (ept1 = cir -> head, ept2 = NULL; ept1 != NULL; ept1 = ept2 -> next) {
                if (!ept1 -> begin) {
                    dserr (610, "dscyc");
                    return (0);
                }
                ept2 = ept1 -> next;
                nedg++;
            }
        }
    }

    if (nedg <= 0) {
        dserr (713, "dscyc");
        return (0);
    }

    /* allocate memory for edges */
    edges = (struct edg *) allocate_objects (EDG, nedg);
    if (edges == NULL) {
        dserr (500, "dscyc");
        return (0);
    }

    /* initialization of index into newly allocated memory */
    ie = 0;
 /* store info: */
    for (tl = ftl; tl != NULL; tl = tl -> next) {
        tor = tl -> tor;
        cir = tor -> central;

        if (tor -> atom[0] == ia)
            ja = tor -> atom[1];
        else
            ja = tor -> atom[0];

        dscross (cir -> axis, cir -> base, zenith);

        if (cir -> head == NULL) {
            (edges + ie) -> ept1 = NULL;
            (edges + ie) -> ept2 = NULL;
            (edges + ie) -> atom = ja;
            (edges + ie) -> next = NULL;
            for (k = 0; k < 3; k++)
                (edges + ie) -> ept1co[k] = 0.0;
            for (k = 0; k < 3; k++)
                (edges + ie) -> ept2co[k] = 0.0;
            ie++;
        }
        else {
            for (ept1 = cir -> head; ept1 != NULL; ept1 = ept2 -> next) {
                ept2 = ept1 -> next;
                if ((ia < ja) == !tor -> reverse) {
                    (edges + ie) -> ept1 = ept1;
                    (edges + ie) -> ept2 = ept2;
                }
                else {
                    (edges + ie) -> ept1 = ept2;
                    (edges + ie) -> ept2 = ept1;
                }
                (edges + ie) -> atom = ja;
                (edges + ie) -> next = NULL;
            /* calculate endpoint coordinates */
                ang = (edges + ie) -> ept1 -> angle;
                cosang = cos (ang);
                sinang = sin (ang);
                for (k = 0; k < 3; k++)
                    tp[k] = cir -> radius *
                        (cosang * cir -> base[k] + sinang * zenith[k]);
                for (k = 0; k < 3; k++)
                    (edges + ie) -> ept1co[k] = cir -> center[k] + tp[k];
                ang = (edges + ie) -> ept2 -> angle;
                cosang = cos (ang);
                sinang = sin (ang);
                for (k = 0; k < 3; k++)
                    tp[k] = cir -> radius *
                        (cosang * cir -> base[k] + sinang * zenith[k]);
                for (k = 0; k < 3; k++)
                    (edges + ie) -> ept2co[k] = cir -> center[k] + tp[k];
                ie++;
            }
        }
    }

 /* calculate north vector from all neighbors */
 /* we want this vector to point towards the solvent */
 /* towards the solvent means away from neighboring atoms */
    for (k = 0; k < 3; k++)
        north[k] = 0.0;
    for (jnbr = nbrs; *jnbr != EOL; jnbr++) {
        ja = *jnbr;
        for (k = 0; k < 3; k++)
            jvec[k] = (ico[k] - *(atmco + 3 * ja + k));
        for (k = 0; k < 3; k++)
            north[k] += jvec[k];
    }

    /* normalize to unit length */
    if (!dsnize (north)) {
        free_objects (EDG, (short *) edges);
        return (0);
    }

    /* one edge that is a complete circle is a special case */
    if (nedg == 1) {
        if (nfree != 1) {
            dserr (714, "dscyc");
            free_objects (EDG, (short *) edges);
            return (0);
        }
        for (k = 0; k < 3; k++)
            zone[k] = cir -> center[k];
    /* project zone onto axis */
        for (k = 0; k < 3; k++)
            vect1[k] = zone[k] - *(atmco + 3 * ia + k);
        dt = dsdot (vect1, north);
        *z = dt;
        /* free temporary memory */
        free_objects (EDG, (short *) edges);
        return (1);
    }

    /* more than one free circle means more than one cycle */
    if (nfree > 0) {
        free_objects (EDG, (short *) edges);
        return (0);
    }

 /* set up pointers from edge to edge */
    nlink = 0;
    for (ie = 0; ie < nedg; ie++)
        for (je = 0; je < nedg; je++) {
            if (ie == je)
                continue;
            ept2 = (edges + ie) -> ept2;        /* 2nd endpnt of ie */
            ept1 = (edges + je) -> ept1;        /* 1st endpnt of je */
            a1 = (edges + ie) -> atom;          /* ie lies along a1 */
            a2 = (edges + je) -> atom;          /* je lies along a2 */
            /* check that atoms match */
            if (a1 != ept1 -> atom)
                continue;
            if (a2 != ept2 -> atom)
                continue;

            /* check that endpoints are close in space */
            for (k = 0; k < 3; k++)
                ept2co[k] = (edges + ie) -> ept2co[k];
            for (k = 0; k < 3; k++)
                ept1co[k] = (edges + je) -> ept1co[k];
            if (dsdis (ept1co, ept2co) > EPSILON)
                continue;
            /* edge je follows edge ie */
            nlink++;
            (edges + ie) -> next = (edges + je);
            break;
        }

    /* check for missing link */
    if (nlink != nedg) {
        free_objects (EDG, (short *) edges);
        return (0);
    }

 /* traverse links starting at first to form 1 cycle */
    ntrav = 0;
    for (eptr = edges; eptr != NULL; eptr = eptr -> next) {
        if (ntrav > 0 && eptr == edges)
            break;              /* back to start */
        ntrav++;
    }
    /* finished with first cycle */
    if (ntrav > nedg) {
        dserr (716, "dscyc");
        return (0);
    }
    /* if edges remain, there is more than one cycle */
    if (ntrav < nedg) {
        free_objects (EDG, (short *) edges);
        return (0);
    }

 /* exactly one, non-free cycle */
 /* calculate zone by averaging endpoints of edges */

    for (k = 0; k < 3; k++)
        zone[k] = 0.0;
    for (ie = 0; ie < nedg; ie++)
        for (k = 0; k < 3; k++)
            zone[k] += (edges + ie) -> ept1co[k];
    for (k = 0; k < 3; k++)
        zone[k] /= nedg;
 /* and projecting zone onto axis */
    for (k = 0; k < 3; k++)
        vect1[k] = zone[k] - *(atmco + 3 * ia + k);
    dt = dsdot (vect1, north);
    *z = dt;                    /* store value in calling function */
    /* free temporary memory */
    free_objects (EDG, (short *) edges);
    return (1);                 /* success */
}



/* transfer dot surface descriptor to global variables */
void dsdget (struct dsdesc *dsd)
{
    maxatom = dsd -> maxatom;
    natom = dsd -> natom;
    atmco = dsd -> atmco;
    atmrad = dsd -> atmrad;
    atmden = dsd -> atmden;
    atmatt = dsd -> atmatt;
    pradius = dsd -> pradius;
    connected = dsd -> connected;

    errflg = dsd -> errflg;
    errstr = dsd -> errstr;
    strcpy (stage, dsd -> stage);
    atom = dsd -> atom;
    verbose = dsd -> verbose;
    nsrfpnt = dsd -> nsrfpnt;
    cvxarea = dsd -> cvxarea;
    renarea = dsd -> renarea;
    cvxsp = dsd -> cvxsp;
    rensp = dsd -> rensp;
    intern = dsd -> intern;
}

void dsdini (struct dsdesc *dsd)                    /* initialize dot surface descriptor */
{
    dsd -> errflg = 0;
    dsd -> errstr = NULL;
    strcpy (dsd -> stage, " ");
    dsd -> atom = EOL;
    dsd -> nsrfpnt = 0;
    dsd -> cvxarea = 0.0;
    dsd -> renarea = 0.0;
    dsd -> cvxsp = (struct cluster **) allocate_pointers (CLUSTER, dsd -> maxatom);
    if (dsd -> cvxsp == NULL) {
        dserr (500, "dsdini");
        return;
    }
    dsd -> rensp = (struct cluster **) allocate_pointers (CLUSTER, dsd -> maxatom);
    if (dsd -> rensp == NULL) {
        dserr (500, "dsdini");
        return;
    }
    dsd -> intern = (struct hidden *) allocate_object(HIDDEN);
    if (dsd -> intern == NULL) {
        dserr (500, "dsdini");
        return;
    }
    dsirini (dsd -> intern);
}

double dsdis (double a[3], double b[3])            /* distance from a to b */
{
    double  d2;

    d2 = dsdis2 (a, b);
    if (d2 < 0.0)
        d2 = 0.0;
    d2 = sqrt (d2);
    return (d2);
}

double dsdis2 (double a[3], double b[3])           /* distance from a to b, squared */
{
    int    k;
    double  d2, r;

    d2 = 0.0;
    for (k = 0; k < 3; k++) {
        r = *a++ - *b++;
        d2 += (r * r);
    }
    return (d2);
}

struct cluster *dsdivid (struct circle *cir)   /* subdivide arcs of circle into points */
{
    int    k, ndiv, ndiv2, i;
    int     genpnts;
    double   ainc, crad, a, ca, sa, circum, ang1, ang2, ang12, vk;
    double   base[3], zenith[3];
    double  *f, *fptr, *cptr, *bptr, *zptr;
    struct cluster *fclu, *clu, *tclu, *pclu;
    struct endpnt  *ept1, *ept2;

 /* check for accessibility */
    if (cir -> gone)
        return (NULL);

    atom = cir -> atom;

    /* determining whether to generate points is not so simple */

    if (cir -> width <= 0.0)
        genpnts = 1;            /* not area call */
    else
        if (cir -> atom == EOL)
            genpnts = 1;        /* atom unknown */
        else
            genpnts = (*(atmatt + cir -> atom) >= POINTS);

    if (fabs (dsnorm (cir -> axis) - 1.0) > EPSILON) {
        dserr (850, "dsdivid");
        return (NULL);
    }
    if (density <= 0.0) {
        dserr (864, "dsdivid");
        return(NULL);
    }
    /* transfer circle radius to local variable */
    crad = cir -> radius;

 /* check for complete circle */

    if (cir -> head == NULL) {
        dsarbp (cir -> axis, base);
        if (errflg)
            return (NULL);
        dscross (cir -> axis, base, zenith);
        circum = 2 * PI * cir -> radius;
        /* generate an even number of points for the circle */
        ndiv2 = 0.5 * circum * sqrt (density) + 0.25;
        if (ndiv2 <= 0)
            ndiv2 = 1;
        ndiv = 2 * ndiv2;
        clu = (struct cluster  *) allocate_object (CLUSTER);
        if (clu == NULL) {
            dserr (500, "dsdivid");
            return (NULL);
        }

        /* compute area associated with circle */
        clu -> area = 2 * PI * cir -> radius * cir -> width;

        if (genpnts) {
            clu -> nmem = ndiv;
            /* allocate memory for subdivision point coordinates */
            f = allocate_doubles (ndiv * 3, 0, F);
            if (f == NULL) {
                dserr (500, "dsdivid");
                return (NULL);
            }
            /* angular increment between subdivision points */
            ainc = 2 * PI / ndiv;
            for (i = 0, a = 0.0, fptr = f; i < ndiv2;
                    i++, a += ainc) {
                ca = cos (a) * crad;
                sa = sin (a) * crad;
                /* pointer arithmetic is faster than array arithmetic */
                cptr = cir -> center;
                bptr = base;
                zptr = zenith;
                /* calculate two points */
                for (k = 0; k < 3; k++) {
                    vk = ca * *bptr++ + sa * *zptr++;
                    *(fptr + ndiv2 * 3) = *cptr - vk;   /* first point */
                    *fptr++ = *cptr++ + vk;             /* opposite point */
                }
            }
        }
        else {                  /* area call -- no points */
            f = NULL;
            clu -> nmem = 0;
        }

        clu -> points = f;
        clu -> normals = NULL;
        return (clu);
    }

 /* remaining case: one or more arcs */

    /* initialization of cluster linked list */
    fclu = NULL;
    pclu = NULL;

    for (ept1 = cir -> head, ept2 = NULL; ept1 != NULL; ept1 = ept2 -> next) {
        ept2 = ept1 -> next;
        ang1 = ept1 -> angle;
        ang2 = ept2 -> angle;
        ang12 = ang2 - ang1;
        dscross (cir -> axis, cir -> base, zenith);
        circum = ang12 * cir -> radius;         /* arc length */
        ndiv = circum * sqrt (density) + 0.5;   /* number of subdivisions */
        if (ndiv <= 0)
            ndiv = 1;
        clu = (struct cluster  *) allocate_object (CLUSTER);
        if (clu == NULL) {
            dserr (500, "dsdivid");
            return (NULL);
        }

        clu -> area = circum * cir -> width;    /* area of band */

        if (genpnts) {
            clu -> nmem = ndiv;
            f = allocate_doubles ((unsigned) ndiv * 3, 0, F);
            if (f == NULL) {
                dserr (500, "dsdivid");
                return (NULL);
            }
            ainc = ang12 / ndiv;                /* angular increment */
            for (i = 0, a = ang1 + 0.5 * ainc, fptr = f; i < ndiv;
                    i++, a += ainc) {
                ca = cos (a) * crad;
                sa = sin (a) * crad;
                /* pointer arithmetic is faster than array arithmetic */
                cptr = cir -> center;
                bptr = cir -> base;
                zptr = zenith;
                /* calculate subdivision point coordinates */
                for (k = 0; k < 3; k++)
                    *fptr++ = *cptr++ +
                        (ca * *bptr++ + sa * *zptr++);
            }
        }
        else {
            f = NULL;
            clu -> nmem = 0;
        }

        clu -> points = f;
        clu -> normals = NULL;

        /* link into list of clusters for arcs on this circle */
        if (fclu == NULL)
            fclu = clu;
        else
            pclu -> next = clu;
        pclu = clu;
    }

    /* merge the clusters into one cluster */
    tclu = dsmerge (fclu);
    if (errflg)
        return (NULL);
    return (tclu);
}

double dsdot (double *a, double  *b)            /* dot product of two vectors */
{
    int k;
    double  dp;

    dp = 0.0;
    for (k = 0; k < 3; k++)
        dp += (*a++ * *b++);
    return (dp);
}

void dsdput (struct dsdesc *dsd)    /* transfer global variables to dot surface descriptor */
{
    dsd -> maxatom = maxatom;
    dsd -> natom = natom;
    dsd -> atmco = atmco;
    dsd -> atmrad = atmrad;
    dsd -> atmden = atmden;
    dsd -> atmatt = atmatt;
    dsd -> pradius = pradius;
    dsd -> connected = connected;

    dsd -> errflg = errflg;
    dsd -> errstr = errstr;
    strcpy (dsd -> stage, stage);
    dsd -> atom = atom;
    dsd -> verbose = verbose;
    dsd -> nsrfpnt = nsrfpnt;
    dsd -> cvxarea = cvxarea;
    dsd -> renarea = renarea;
    dsd -> cvxsp = cvxsp;
    dsd -> rensp = rensp;
    dsd -> intern = intern;
}



void dsfclu (struct cluster *clu)                    /* free cluster */
{
    if (clu == NULL) {
        dserr (760, "dsfclu");
        return;
    }
    /* check for consistent opinion about presence of points */
    if ((clu -> nmem == 0) != (clu -> points == NULL)) {
        dserr (625, "dsfclu");
        return;
    }
    if (clu -> nmem > 0)
        free_doubles (clu -> points, 0, POINTS);
    if (clu -> normals != NULL)
        free_doubles (clu -> normals, 0, NORMALS);
    free_object (CLUSTER, (short *) clu);
}

int dsfree (struct dsdesc *dsd)                    /* free surface calculated by ds */
{
    int     i;
    struct cluster *clu;

    if (dsd == NULL) {
        set_error1 ("null pointer passed to dsfree");
        return (550);
    }

    atom = EOL;
    dsdget (dsd);                       /* move record to global */

    /* be careful */
    if (cvxsp == NULL) {
        dserr (555, "dsfree");
        goto error;
    }
    if (rensp == NULL) {
        dserr (555, "dsfree");
        goto error;
    }
    if (intern == NULL) {
        dserr (555, "dsfree");
        goto error;
    }

    /* for each atom, free convex and reentrant clusters */
    for (i = 0; i < natom; i++) {
        atom = i;
        clu = *(cvxsp + i);
        if (clu != NULL)
            dsfclu (clu);
        if (errflg)
            break;
        clu = *(rensp + i);
        if (clu != NULL)
            dsfclu (clu);
        if (errflg)
            break;
    }
    if (errflg)
        goto error;

    /* free arrays of pointers to clusters */
    free_pointers (CLUSTER, cvxsp);
    free_pointers (CLUSTER, rensp);

    /* free internal record */
    free_object (HIDDEN, (short *) intern);

    /* make sure pointers are null */
    cvxsp = NULL;
    rensp = NULL;
    intern = NULL;

error:

    dsdput (dsd);               /* store global into record */
    return (errflg);
}


void dsgetir (struct hidden *ir)                    /* move internal record to global */
{
    ntori = ir -> ntori;
    nprobe = ir -> nprobe;
    nlow = ir -> nlow;
    ndead = ir -> ndead;
    hedtor = ir -> hedtor;
    taltor = ir -> taltor;
    hedprb = ir -> hedprb;
    talprb = ir -> talprb;
    lowp = ir -> lowp;
    hedrun = ir -> hedrun;
    hedcav = ir -> hedcav;
    hedpl = ir -> hedpl;
    problem = ir -> problem;
    nplist = ir -> nplist;
    minwid = ir -> minwid;
    rootbox = ir -> rootbox;
}

void dsgtor ()                       /* generate tori */
{
    int    ia, k, ja, ma;
    double   irad, jrad, dij2, dij, root1, root2, radius;
    double   ico[3], jco[3], vij[3], axis[3], center[3];
    atomnum * nbrs, *jnbr, *mut, *m;
    struct run *irun;
    struct circle  *cir;
    struct torus   *tor;

    if (natom < 2)
        return;

    for (ia = 0; ia < natom; ia++) {
        atom = ia;
        irun = *(hedrun + ia);          /* get packed neighbors */
        if (irun == NULL)
            continue;
        nbrs = dsupk (irun);            /* unpack neighbors */
        if (errflg)
            return;
        if (nbrs == NULL)
            continue;
        /* transfer atomic info to local variables */
        irad = *(atmrad + ia) + pradius;
        for (k = 0; k < 3; k++)
            ico[k] = *(atmco + 3 * ia + k);

        for (jnbr = nbrs; *jnbr != EOL; jnbr++) {
            ja = *jnbr;                 /* retrieve atom number from list */
            if (ja <= ia)               /* for complete surface: ia < ja */
                continue;
            /* if neither torus has even area attention, skip this neighbor */
            if (*(atmatt + ia) < AREA && *(atmatt + ja) < AREA)
                continue;
            /* transfer atomic info to local variables */
            jrad = *(atmrad + ja) + pradius;
            for (k = 0; k < 3; k++)
                jco[k] = *(atmco + 3 * ja + k);

            /* geometric calculations for torus central circle */
            /* see M.L. Connolly, J. Appl. Crystallogr. 16, 548-558 (1983) */
            for (k = 0; k < 3; k++)
                vij[k] = jco[k] - ico[k];
            dij2 = vij[0] * vij[0] + vij[1] * vij[1] + vij[2] * vij[2];
            if (dij2 <= 0.0) {
                dserr (852, "dsgtor");
                return;
            }
            dij = sqrt (dij2);
            if (dij <= 0.0) {
                dserr (852, "dsgtor");
                return;
            }
            for (k = 0; k < 3; k++)
                axis[k] = vij[k] / dij;
            for (k = 0; k < 3; k++)
                center[k] = 0.5 * (ico[k] + jco[k])
                    + 0.5 * vij[k] * (irad * irad - jrad * jrad) / dij2;
            root1 = (irad + jrad) * (irad + jrad) - dij2;
            /* the two atoms ought to be neighbors */
            if (root1 < 0.0)
                continue;               /* for some reason, they aren't */
            root1 = sqrt (root1);
            root2 = dij2 - (irad - jrad) * (irad - jrad);
            if (root2 < 0.0)
                continue;               /* one atom inside the other */
            root2 = sqrt (root2);
            radius = 0.5 * root1 * root2 / dij;
        /* allocate torus central circle */
            cir = dsnwcir (center, radius, axis);
            if (errflg)
                return;
        /* get list of mutual neighbors of both atoms */
            mut = dsmut (ia, ja);
            if (errflg)
                return;
            if (mut != NULL) {
                /* for each mutual neighbor:
                    cut torus central circle with its expanded sphere */
                for (m = mut; *m != EOL; m++) {
                    ma = *m;
                    dscsi (cir, ma);
                    if (errflg)
                        return;
                    if (cir -> gone)    /* circle is inaccessible */
                        break;
                }
                free_atomnums (mut);
            }
            if (cir -> gone) {
				/* no endpoints to free */
                free_object (CIRCLE, (short *) cir);
                continue;
            }

        /* we have a torus that is entirely or partially accessible */
            tor = dsnwtor (ia, ja, cir);
            if (errflg)
                return;
            if (tor == NULL) {
                dserr (750, "dsgtor");
                return;
            }
            /* for connected rolling,
               we need only one torus to initially position the probe */
            if (connected) {
                free_atomnums (nbrs);
                return;
            }
        }                       /* end of jnbr loop */
        free_atomnums (nbrs);
    }                           /* end of ia loop */
    return;
}


/* height of probe center above plane passing through three atom centers */
double dshei (struct probe *prb)
{
    int    k;
    int     i1, i2, i3;
    double   h;
    double   v1[3], v2[3], v3[3], v4[3];
    double  *a1, *a2, *a3;

    i1 = prb -> atom[0];
    i2 = prb -> atom[1];
    i3 = prb -> atom[2];
    atom = i1;
    a1 = atmco + 3 * i1;
    a2 = atmco + 3 * i2;
    a3 = atmco + 3 * i3;
    /* vectors from atom i1 to i2 and from i1 to i3 */
    for (k = 0; k < 3; k++) {
        v1[k] = *(a2 + k) - *(a1 + k);
        v2[k] = *(a3 + k) - *(a1 + k);
    }
    /* vector perpendicular to plane */
    dscross (v1, v2, v3);
    if (!dsnize (v3)) {
        dserr (800, "dshei");
        return (0.0);
    }

    /* vector from plane to probe */
    for (k = 0; k < 3; k++)
        v4[k] = prb -> center[k] - *(a1 + k);
    h = dsdot (v3, v4);         /* component along normal vector */
    return (h);
}


void dsinch ()                        /* input checking */
{
    int     a, att, k;
    double   rad, den;

    atom = EOL;
    if (connected != 0 && connected != 1) {
        dserr (574, "dsinch");
        return;
    }
    if (pradius < MINPRADIUS || pradius > MAXPRADIUS) {
        dserr (572, "dsinch");
        return;
    }
    if (maxatom < MINATOM || maxatom > MAXATOM) {
        dserr (575, "dsinch");
        return;
    }
    if (natom < MINATOM || natom > maxatom) {
        dserr (576, "dsinch");
        return;
    }
    if (atmco == NULL) {
        dserr (775, "dsinch");
        return;
    }
    if (atmrad == NULL) {
        dserr (775, "dsinch");
        return;
    }
    if (atmden == NULL) {
        dserr (775, "dsinch");
        return;
    }
    if (atmatt == NULL) {
        dserr (775, "dsinch");
        return;
    }

    for (a = 0; a < natom; a++) {
        atom = a;
        for (k = 0; k < 3; k++)
            if (fabs (*(atmco + 3 * a + k)) > MAXCOORDINATE) {
                dserr (577, "dsinch");
                return;
            }
        rad = *(atmrad + a);
        den = *(atmden + a);
        att = *(atmatt + a);
        if (den < MINDENSITY || den > MAXDENSITY) {
            dserr (573, "dsinch");
            return;
        }
        if (rad < MINARADIUS || rad > MAXARADIUS) {
            dserr (571, "dsinch");
            return;
        }
        if (att < MINATTENTION || att > MAXATTENTION) {
            dserr (570, "dsinch");
            return;
        }
    }
}

void dsirini (struct hidden *ir)                    /* internal record initialization */
{
    ir -> ntori = 0;
    ir -> nprobe = 0;
    ir -> nlow = 0;
    ir -> ndead = 0;
    ir -> hedtor = NULL;
    ir -> hedprb = NULL;
    ir -> taltor = NULL;
    ir -> talprb = NULL;
    ir -> lowp = NULL;
    ir -> hedrun = NULL;
    ir -> hedcav = NULL;
    ir -> hedpl = NULL;
    ir -> problem = NULL;
    ir -> nplist = 0;
    ir -> minwid = 0.0;
    ir -> rootbox = NULL;
}


int dslen (struct run *first)                   /* number of atoms in linked list of runs */
{
    int     n;
    struct run *r;

    n = 0;
    if (first == NULL)
        return (0);
    for (r = first; r != NULL; r = r -> next)
        n += (r -> atom[1] - r -> atom[0] + 1);
    return (n);
}

void dslfa (struct circle *cir)             /* leave first arc of circle, delete the rest */
{
    struct endpnt  *ept1, *ept2, *next;

    atom = cir -> atom;
    ept1 = cir -> head;
    if (ept1 == NULL) {
        dserr (612, "dslfa");
        return;
    }
    ept2 = ept1 -> next;
    if (ept2 == NULL) {
        dserr (612, "dslfa");
        return;
    }
 /* free the rest */
    next = NULL;
    for (ept1 = ept2 -> next; ept1 != NULL; ept1 = next) {
        next = ept1 -> next;
        free_object (ENDPNT, (short *) ept1);
    }
    ept2 -> next = NULL;        /* 2nd endpoint is now end of list */
}

void dslocal ()                       /* form neighbor lists of all atoms */
{
    int    ia, ja;
    double  d2, sum2;
    int     k;
    double   maxrad, extra, minsep2, erad;
    double   rect[3][2];
    double  *atomi, *atomj;
    atomnum * nbrs, *this, *jnbr, *inrect;

 /* calculate maximum atomic radius */
    maxrad = 0.0;
    for (ia = 0; ia < natom; ia++)
        if (*(atmrad + ia) > maxrad)
            maxrad = *(atmrad + ia);
    extra = 2 * pradius + maxrad;

    dsbtree ();                 /* create box octree */
    if (errflg)
        return;

 /* allocate memory for head pointers of run-length encoded nbr lists */
    hedrun = (struct run  **) allocate_pointers (RUN, natom);
    if (hedrun == NULL) {
        dserr (500, "dslocal");
        return;
    }

 /* allocate temporary storage to be used for each list of atom nbrs */
    nbrs = allocate_atomnums ((unsigned) natom);
    if (nbrs == NULL) {
        dserr (500, "dslocal");
        return;
    }

    /* minimum atom center sep^2 */
    minsep2 = MINSEPARATION * MINSEPARATION;

    for (ia = 0; ia < natom; ia++) {            /* loop through all atoms */
        atom = ia;
        if (*(atmatt + ia) <= IGNORE) {         /* ignored atoms get no list */
            *(hedrun + ia) = NULL;
            continue;
        }
        this = nbrs;                    /* initialize pointer for storing nbr */
                 atomi = atmco + 3 * ia;         /* pointer to atomic coordinates */
        /* set up rectangular range for searching box tree */
        /* make it big enough to include any possible nbr of atom ia */
        for (k = 0; k < 3; k++) {
            rect[k][0] = *(atomi + k) - *(atmrad + ia) - extra;
            rect[k][1] = *(atomi + k) + *(atmrad + ia) + extra;
        }
        /* get atoms in rectangle */
        inrect = dsgair (rect, rootbox);
        if (errflg)
            return;
        if (inrect == NULL) {
            dserr (750, "dslocal");
            return;
        }
        atom = ia;              /* restore error variable after dsgair call */
        erad = *(atmrad + ia) + 2 * pradius;    /* doubly expanded radius */

        /* loop through atoms in rectangular range */
        for (jnbr = inrect; *jnbr != EOL; jnbr++) {
            ja = *jnbr;         /* retrieve atom number from list */
            if (ia == ja)
                continue;
            if (*(atmatt + ja) <= IGNORE)       /* ignore atom */
                continue;
            atomj = atmco + 3 * ja;     /* pointer to ja center */
            d2 = dsdis2 (atomi, atomj); /* distance between centers, squared */
            sum2 = erad + *(atmrad + ja);
            sum2 = sum2 * sum2;         /* sum of expanded radii, squared */
            if (d2 >= sum2) continue;   /* skip atom if too far away */

            /* the two atoms are close enough for probe sphere
               to touch both of them simultaneously, that is,
               they are neighbors */

            *this++ = (unsigned short) ja;               /* store atom number in new list */

         /* some error checking */
            if (d2 <= 0.0) {
                dserr (852, "dslocal");
                return;
            }
            if (d2 <= minsep2) {
                dserr (585, "dslocal");
                return;
            }
        }
        free_atomnums (inrect);         /* free list from box tree search */
        *this = EOL;                    /* mark End Of List */
        *(hedrun + ia) = dspak (nbrs);  /* pack neighbors and store head ptr */
        if (errflg)
            return;
    }

    free_atomnums (nbrs);               /* free temporary storage */

 /* recursively free box tree */
    dsfbox (rootbox);
    if (errflg)
        return;
	free_cache(BOX);
}

/* 

Copyright 1986 by Michael L. Connolly
All Rights Reserved

*/
