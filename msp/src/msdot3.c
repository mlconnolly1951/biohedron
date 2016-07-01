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


struct cluster *dsmerge (struct cluster *fclu)          /* merge clusters */
{
    int    k, i, nmem;
    double *f, *g;
    int     nrmls, n;
    unsigned    size;
    double   totarea;
    double  *fn, *gn;
    struct cluster *clu, *tclu, *next;

    if (fclu == NULL) {
        dserr (760, "dsmerge");
        return (NULL);
    }

    if (fclu -> next == NULL)
        return (fclu);
    nrmls = (fclu -> normals != NULL);

    /* count total number of points and accumulate total area */
    n = 0;
    totarea = 0.0;
    for (clu = fclu; clu != NULL; clu = clu -> next) {
        n += clu -> nmem;
        totarea += clu -> area;
    }

    /* create a new cluster */
    tclu = (struct cluster *) allocate_object (CLUSTER);
    if (tclu == NULL) {
        dserr (500, "dsmerge");
        return (NULL);
    }
    tclu -> nmem = n;
    tclu -> area = totarea;

    if (n > 0) {
        /* allocate space for points */
        size = n * 3 * sizeof (double);
        tclu -> points = allocate_doubles (n * 3, 0, POINTS);
        if (tclu -> points == NULL) {
            dserr (500, "dsmerge");
            return (NULL);
        }
        f = tclu -> points;
        if (nrmls) {
            /* allocate space for normals */
            tclu -> normals = allocate_doubles (n * 3, 0, NORMALS);
            if (tclu -> normals == NULL) {
                dserr (500, "dsmerge");
                return (NULL);
            }
            fn = tclu -> normals;
        }

        /* transfer points and normals from old clusters to new cluster */
        for (clu = fclu, next = NULL; clu != NULL; clu = next) {
            nmem = clu -> nmem;
            g = clu -> points;
            if (nmem != 0 && g != NULL) {
                for (i = 0; i < nmem; i++)
                    for (k = 0; k < 3; k++)
                        *f++ = *g++;
                free_doubles (clu -> points, 0, POINTS);          /* free old memory */
                if (nrmls) {
                    gn = clu -> normals;
                    if (gn == NULL) {
                        dserr (627, "dsmerge");
                        return (NULL);
                    }
                    for (i = 0; i < nmem; i++)
                        for (k = 0; k < 3; k++)
                            *fn++ = *gn++;
                    free_doubles (clu -> normals, 0, NORMALS);     /* free old memory */
                }
            }
            next = clu -> next;
            free_object (CLUSTER, (short *) clu);                        /* free old cluster */
        }
    }
    else {
        tclu -> points = NULL;
        tclu -> normals = NULL;
    }
    return (tclu);
}


atomnum *dsmut (int ia, int ja)        /* return mutual neighbor list of two atoms */
{
    atomnum * nbrs;
    struct run *f1, *f2, *and12;

    atom = ia;

    /* set pointers to first run of each atom,
                        return null if either empty */

    f1 = *(hedrun + ia);
    if (f1 == NULL)
        return (NULL);
    f2 = *(hedrun + ja);
    if (f2 == NULL)
        return (NULL);

    and12 = dsand (f1, f2);     /* call function to compute intersection */
    if (errflg)
        return (NULL);
    if (and12 == NULL)
        return (NULL);

    nbrs = dsupk (and12);       /* unpack linked list of runs into array */
    if (errflg)
        return (NULL);
    dsfrun (and12);             /* free temporary memory */
    if (errflg)
        return (NULL);

    return (nbrs);
}


int dsnize (double a[3])                      /* normalize vector */
{
    int    k;
    double  na;

    na = dsnorm (a);
    if (na <= 0.0)
        return (0);             /* failure: zero norm */
    for (k = 0; k < 3; k++)
        a[k] = a[k] / na;
    return (1);                 /* success */
}

void dsnoend (struct circle *cir)                   /* eliminate endpoints of circle */
{
    struct endpnt *first, *r, *next;

    atom = cir -> atom;
    first = cir -> head;
    if (first == NULL)
        return;
    for (r = first, next = NULL; r != NULL; r = next) {
        next = r -> next;
        free_object (ENDPNT, (short *) r);
    }
    cir -> head = NULL;
}

double dsnorm (double a[3])              /* return norm of vector */
{
    double  dp;

    dp = a[0] * a[0] + a[1] * a[1] + a[2] * a[2];
    if (dp < 0.0)
        return (0.0);
    dp = sqrt (dp);
    return (dp);
}

void dsnti (struct circle *cir, atomnum *nbrnt) /* non-tori neighbor intersection with cir */
{
    int    k;
    int     ja;
    double   a, ca, sa, ang1, ang2, ang12, jrad, jrad2, dj2;
    double   base[3], zenith[3], midpnt[3];
    double  *jco;
    atomnum * jnbr, *knbr;
    struct endpnt  *ept1, *ept2;

 /* check for accessibility */
    if (cir -> gone)
        return;

 /* check for complete circle */

    if (cir -> head == NULL) {
    /* calculate an arbitrary point on circle */
        dsarbp (cir -> axis, base);
        if (errflg)
            return;
        for (k = 0; k < 3; k++)
            midpnt[k] = cir -> center[k] + cir -> radius * base[k];
    /* check whether midpoint inside expanded sphere of non-torus nbr */
        for (jnbr = nbrnt; *jnbr != EOL; jnbr++) {
            ja = *jnbr;
            jco = atmco + 3 * ja;
            jrad = *(atmrad + ja) + pradius;
            jrad2 = jrad * jrad;
            dj2 = dsdis2 (midpnt, jco);
            if (dj2 < jrad2) {
                cir -> gone = 1;
            /* shift to first if not first (for next latitude) */
                if (jnbr != nbrnt) {
                    for (knbr = jnbr - 1; knbr >= nbrnt; knbr--)
                        *(knbr + 1) = *knbr;
                    *nbrnt = (unsigned short) ja;
                }
                return;
            }
        }
        return;
    }

 /* one or more arcs */
    for (ept1 = cir -> head, ept2 = NULL; ept1 != NULL; ept1 = ept2 -> next) {
        ept2 = ept1 -> next;
        ang1 = ept1 -> angle;
        ang2 = ept2 -> angle;
        ang12 = ang2 - ang1;
        if (ang12 < 0.0) {
            dserr (611, "dsnti");
            return;
        }
        /* calculate midpoint of arc */
        dscross (cir -> axis, cir -> base, zenith);
        a = ang1 + 0.5 * ang12;
        ca = cos (a);
        sa = sin (a);
        for (k = 0; k < 3; k++)
            midpnt[k] = cir -> center[k] + cir -> radius
                * (ca * cir -> base[k] + sa * zenith[k]);
    /* check whether midpoint inside expanded sphere of non-torus nbr */
        for (jnbr = nbrnt; *jnbr != EOL; jnbr++) {
            ja = *jnbr;
            atom = ja;
            jco = atmco + 3 * ja;
            jrad = *(atmrad + ja) + pradius;
            dj2 = dsdis2 (midpnt, jco);
            jrad2 = jrad * jrad;
            if (dj2 < jrad2) {
            /* midpoint eaten, therefore arc eaten */
            /* mark arc endpoints to be removed */
                ept1 -> angle = -1.0;
                ept2 -> angle = -1.0;
            /* swap with first if not first (for next latitude) */
                if (jnbr != nbrnt) {
                    for (knbr = jnbr - 1; knbr >= nbrnt; knbr--)
                        *(knbr + 1) = *knbr;
                    *nbrnt = (unsigned short) ja;
                }
                break;
            }
        }
    }

 /* take care of deletions -- mark gone if necessary */
    dsccir (cir);

    return;
}

/* here follow some routines that allocate memory for geometric structures,
   and store some information in them */

struct circle  *dsnwcir (double center[3], double radius, double axis[3])
{
    int     k;
    struct circle  *cir;

    cir = allocate_circle ();
    if (cir == NULL) {
        dserr (500, "dsnwcir");
        return (NULL);
    }
    cir -> radius = radius;
    cir -> width = 0.0;
    for (k = 0; k < 3; k++) {
        cir -> center[k] = center[k];
        cir -> axis[k] = axis[k];
    }
    cir -> atom = EOL;
    return (cir);
}

struct endpnt *dsnwept (int begin, double angle, int a)
{
    struct endpnt  *ept;

    atom = a;
    ept = (struct endpnt   *) allocate_object (ENDPNT);
    if (ept == NULL) {
        dserr (500, "dsnwept");
        return (NULL);
    }
    ept -> begin = (unsigned short) begin;
    ept -> angle = angle;
    ept -> atom = (unsigned short) a;
    return (ept);
}

struct probe *dsnwprb (int atom1, int atom2, int atom3, double center[3])
{
    int     k;
    struct probe   *prb;

    atom = atom1;
    prb = (struct probe *) allocate_object (PROBE);
    if (prb == NULL) {
        dserr (500, "dsnwprb");
        return (NULL);
    }
    /* link into list */
    if (hedprb == NULL)
        hedprb = prb;
    else
        talprb -> next = prb;
    talprb = prb;
    nprobe++;           /* increment counter */
    /* store info */
    for (k = 0; k < 3; k++)
        prb -> center[k] = center[k];
    prb -> atom[0] = (unsigned short) atom1;
    prb -> atom[1] = (unsigned short) atom2;
    prb -> atom[2] = (unsigned short) atom3;
    return (prb);
}

struct torus *dsnwtor (int atom1, int atom2, struct circle *cir)
{
    struct torus   *tor;

    atom = atom1;
    tor = (struct torus *) allocate_object (TORUS);
    if (tor == NULL) {
        dserr (500, "dsnwtor");
        return (NULL);
    }
 /* add to linked list of tori */
    if (hedtor == NULL)
        hedtor = tor;
    else
        taltor -> next = tor;
    taltor = tor;
    /* store info */
    tor -> atom[0] = (unsigned short) atom1;
    tor -> atom[1] = (unsigned short) atom2;
    tor -> central = cir;
    tor -> cusp = (cir -> radius < pradius);
    tor -> reverse = (atom2 < atom1);   /* used in dscyc */
    ntori++;            /* increment counter */
    return (tor);
}


int dsouch (struct dsdesc *dsd)                    /* output check of surface points */
{
    int     i, j, k, nmem;
    double   d, irad;
    double   pco[3];
    double  *normals, *snorm, *sco, *ico, *points;
    struct cluster *iren;

    dsstage ("surface point checking");

    if (dsd == NULL) {
        set_error1 ("null pointer passed to dsouch");
        return (550);
    }
    atom = EOL;
    dsdget (dsd);

    if (cvxsp == NULL) {
        dserr (555, "dsouch");
        dsdput (dsd);
        return (errflg);
    }
    if (rensp == NULL) {
        dserr (555, "dsouch");
        dsdput (dsd);
        return (errflg);
    }
    if (intern == NULL) {
        dserr (555, "dsouch");
        dsdput (dsd);
        return (errflg);
    }

    for (i = 0; i < natom; i++) {       /* atom loop */
        atom = i;
        iren = *(rensp + i);    /* only reentrant points checked */
        if (iren == NULL)
            continue;
        nmem = iren -> nmem;
        if (nmem <= 0)
            continue;
        points = iren -> points;
        normals = iren -> normals;
        irad = *(atmrad + i);
        ico = atmco + 3 * i;
        for (j = 0; j < nmem; j++) {    /* surface point loop */
            sco = points + 3 * j;
            /* check for surface point inside atom */
            d = dsdis (ico, sco);
            if (d - irad < -(10.0 * EPSILON)) {
                dserr (920, "dsouch");
                goto error;
            }
            /* check for probe sphere collision with atom */
            if (normals != NULL) {
                snorm = normals + 3 * j;
                for (k = 0; k < 3; k++)
                    pco[k] = *(sco + k) +
                        pradius * *(snorm + k);
                d = dsdis (ico, pco) - (irad + pradius);
                if (fabs (d) > 10.0 * EPSILON) {
					fprintf (stderr, "distance = %12.6f\n", dsdis (ico, pco));
					fprintf (stderr, "radii    = %12.6f\n", (irad + pradius));
                    dserr (925, "dsouch");
                    goto error;
                }
            }
        }                       /* end of cluster loop */
    }                           /* end of atom loop */

    if (verbose) inform("         all surface points pass check");

error:
    dsdput (dsd);
    return (errflg);
}


struct run *dspak (atomnum *nbrs)        /* pack neighbors of atom into runs */
{
    struct run *r, *p, *first;
    int     sflag;
    atomnum * ptr;

    if (nbrs == NULL)
        return (NULL);
    if (*nbrs == EOL)
        return (NULL);

 /* check whether nbrs array is sorted */

    if (!dsans (nbrs)) {
        dserr (685, "dspak");
        return (NULL);
    }

    sflag = 1;          /* turn on start run flag */
    r = NULL;           /* initialize run to null */
    /* go through the list, creating contiguous runs */
    for (ptr = nbrs; *ptr != EOL; ptr++) {
        p = r;
        /* the first time through the loop, sflag is always on;
           subsequently, it is on if we just passed a gap */
        if (ptr > nbrs)
            sflag = (*ptr > *(ptr - 1) + 1);
        /* start a new run */
        if (sflag) {
            r = (struct run *) allocate_object (RUN);
            if (r == NULL) {
                dserr (500, "dspak");
                return (NULL);
            }
            /* link into run list */
            if (p != NULL)
                p -> next = r;
            else
                first = r;
            r -> atom[0] = *ptr;        /* store first atom of run */
        }
        r -> atom[1] = *ptr;    /* store this atom of run as last (overwrite) */
                 atom = r -> atom[1];    /* for possible error message */
    }
    return (first);
}

void dspart (struct probe *prb)  /* partition dots of concave face: to atoms */
{
    int    i, k;
    short  *closest;
    int     i1, i2, i3, att1, att2, att3, a1, a2, a3, n1, n2, n3, nps;
    double   d1, d2, d3;
    double   pco[3];
    double  *f1, *f2, *f3, *f, *fn1, *fn2, *fn3;
    struct cluster *clu, *clu1, *clu2, *clu3;

    /* transfer info to local variables and arrays */
    for (k = 0; k < 3; k++)
        pco[k] = prb -> center[k];
    a1 = prb -> atom[0];
    a2 = prb -> atom[1];
    a3 = prb -> atom[2];
    atom = a1;
    att1 = *(atmatt + a1);
    att2 = *(atmatt + a2);
    att3 = *(atmatt + a3);
    /* initialization */
    n1 = 0;
    n2 = 0;
    n3 = 0;
    /* if all atoms ignored or blockers,
            there should be no probe position */
    if (att1 < AREA && att2 < AREA && att3 < AREA) {
        dserr (866, "dspart");
        return;
    }
    /* at least one atom must request surface if we are to continue */
    if (att1 <= AREA && att2 <= AREA && att3 <= AREA)
        return;

    clu = dsprb (prb);          /* call function to create concave face */
    if (errflg)
        return;
    if (clu == NULL)
        return;

    nps = clu -> nmem;
    if (nps <= 0) {
        dserr (630, "dspart");
        return;
    }

    /* allocate array to store cases saying
                which atom is closest to point */
    closest = allocate_shorts ((unsigned) nps);
    if (closest == NULL) {
        dserr (500, "dspart");
        return;
    }
    f = clu -> points;
    if (f == NULL) {
        dserr (625, "dspart");
        return;
    }
 /* for each point calculate distance to
             van der Waals sphere of each atom */
    for (i = 0; i < nps; i++, f += 3) {
        d1 = dsdis (f, atmco + 3 * a1) - *(atmrad + a1);
        d2 = dsdis (f, atmco + 3 * a2) - *(atmrad + a2);
        d3 = dsdis (f, atmco + 3 * a3) - *(atmrad + a3);
        /* store choice in temporary array, increment appropriate counter */
        if (d1 < d2 && d1 < d3) {
            *(closest + i) = 1;
            n1++;
        }
        else
            if (d2 < d3) {
                *(closest + i) = 2;
                n2++;
            }
            else {
                *(closest + i) = 3;
                n3++;
            }
    }
    /* if first atom has some points, allocate a cluster for it */
    if (n1 > 0 && att1 >= POINTS) {
        clu1 = (struct cluster *) allocate_object (CLUSTER);
        if (clu1 == NULL) {
            dserr (500, "dspart");
            return;
        }
        clu1 -> area = clu -> area * (double) n1 / (double) nps;
        fn1 = NULL;
        f1 = allocate_doubles ((unsigned) n1 * 3, 0, F1);
        if (f1 == NULL) {
            dserr (500, "dspart");
            return;
        }
        clu1 -> nmem = n1;
        if (att1 >= NORMALS) {
            fn1 = allocate_doubles ((unsigned) n1 * 3, 0, FN1);
            if (fn1 == NULL) {
                dserr (500, "dspart");
                return;
            }
        }
        clu1 -> points = f1;
        clu1 -> normals = fn1;
    }

    /* if second atom has some points, allocate a cluster for it */
    if (n2 > 0 && att2 >= POINTS) {
        clu2 = (struct cluster *) allocate_object (CLUSTER);
        if (clu2 == NULL) {
            dserr (500, "dspart");
            return;
        }
        clu2 -> area = clu -> area * (double) n2 / (double) nps;
        fn2 = NULL;
        f2 = allocate_doubles ((unsigned) n2 * 3, 0, F2);
        if (f2 == NULL) {
            dserr (500, "dspart");
            return;
        }
        clu2 -> nmem = n2;
        if (att2 >= NORMALS) {
            fn2 = allocate_doubles ((unsigned) n2 * 3, 0, FN2);
            if (fn2 == NULL) {
                dserr (500, "dspart");
                return;
            }
        }
        clu2 -> points = f2;
        clu2 -> normals = fn2;
    }

    /* if third atom has some points, allocate a cluster for it */
    if (n3 > 0 && att3 >= POINTS) {
        clu3 = (struct cluster *) allocate_object (CLUSTER);
        if (clu3 == NULL) {
            dserr (500, "dspart");
            return;
        }
        clu3 -> area = clu -> area * (double) n3 / (double) nps;
        fn3 = NULL;
        f3 = allocate_doubles ((unsigned) n3 * 3, 0, F3);
        if (f3 == NULL) {
            dserr (500, "dspart");
            return;
        }
        clu3 -> nmem = n3;
        if (att3 >= NORMALS) {
            fn3 = allocate_doubles ((unsigned) n3 * 3, 0, FN3);
            if (fn3 == NULL) {
                dserr (500, "dspart");
                return;
            }
        }
        clu3 -> points = f3;
        clu3 -> normals = fn3;
    }

    /* initialize indices */
    i1 = 0;
    i2 = 0;
    i3 = 0;
    /* transfer surface points and normals to appropriate atom */
    for (i = 0, f = clu -> points; i < nps; i++, f += 3) {
        switch ((long) * (closest + i)) {
            case 1:
                if (att1 < POINTS)
                    break;
                for (k = 0; k < 3; k++)
                    *(f1 + 3 * i1 + k) = *(f + k);
                if (att1 >= NORMALS)
                    for (k = 0; k < 3; k++)
                        *(fn1 + 3 * i1 + k) = (pco[k] - *(f + k)) / pradius;
                i1++;
                break;
            case 2:
                if (att2 < POINTS)
                    break;
                for (k = 0; k < 3; k++)
                    *(f2 + 3 * i2 + k) = *(f + k);
                if (att2 >= NORMALS)
                    for (k = 0; k < 3; k++)
                        *(fn2 + 3 * i2 + k) = (pco[k] - *(f + k)) / pradius;
                i2++;
                break;
            case 3:
                if (att3 < POINTS)
                    break;
                for (k = 0; k < 3; k++)
                    *(f3 + 3 * i3 + k) = *(f + k);
                if (att3 >= NORMALS)
                    for (k = 0; k < 3; k++)
                        *(fn3 + 3 * i3 + k) = (pco[k] - *(f + k)) / pradius;
                i3++;
                break;
            default:
                dserr (658, "dspart");
                return;
        }
    }

    /* free temporary memory */
    free_shorts (closest);

    /* link each new cluster into concave cluster list for that atom */

    if (n1 > 0 && att1 >= POINTS) {
        clu1 -> next = *(hedcav + a1);
        *(hedcav + a1) = clu1;
    }
    if (n2 > 0 && att2 >= POINTS) {
        clu2 -> next = *(hedcav + a2);
        *(hedcav + a2) = clu2;
    }
    if (n3 > 0 && att3 >= POINTS) {
        clu3 -> next = *(hedcav + a3);
        *(hedcav + a3) = clu3;
    }

    /* free memory of original cluster */
    dsfclu (clu);
    if (errflg)
        return;
    return;
}

void dsplace ()              /* place probes tangent to three atoms */
{
    int    k;
    int     att1, att2, att3, atom1, atom2, atom3;
    double   cosang, sinang;
    double   zenith[3], center[3], tp[3];
    struct torus   *tor;
    struct endpnt  *ept;
    struct circle  *cir;
    struct probe   *prb;

    if (natom <= 2)
        return;

    /* loop through tori,
            looking for endpoints of arcs on central circles */

    for (tor = hedtor; tor != NULL; tor = tor -> next) {
        cir = tor -> central;
        if (cir -> head == NULL)        /* no probe placement for free torus */
            continue;

        dscross (cir -> axis, cir -> base, zenith); /* base vector exists */

        /* one probe placement for each endpoint */
        for (ept = cir -> head; ept != NULL; ept = ept -> next) {
            atom1 = tor -> atom[0];     /* atom on one side of arc */
            atom2 = tor -> atom[1];     /* atom on other side of arc */
            atom3 = ept -> atom;                /* atom at end of arc */
            atom = atom3;
            att1 = *(atmatt + atom1);
            att2 = *(atmatt + atom2);
            att3 = *(atmatt + atom3);

            /* careful thought is required to avoid placement duplication */

            /* for connected rolling, we return after first probe is placed,
               so we don't need to worry about duplication */

            if (!connected) {
                if (att1 >= AREA && att2 >= AREA) {
                    if (att3 >= AREA) {
                    /* 1/3 probe place per endpnt */
                        if (atom3 < atom1)
                            continue;
                        if (atom3 < atom2)
                            continue;
                    }
                /* okay */
                }
                else
                    if (att3 >= AREA)
                        continue;
                    else
                        if (att1 >= AREA && atom3 < atom2)
                            continue;
                        else
                            if (att2 >= AREA && atom3 < atom1)
                                continue;
            }

        /* calculate probe center */
            cosang = cos (ept -> angle);
            sinang = sin (ept -> angle);
            for (k = 0; k < 3; k++)
                tp[k] = cir -> radius * (cosang * cir -> base[k]
                        + sinang * zenith[k]);
            for (k = 0; k < 3; k++)
                center[k] = cir -> center[k] + tp[k];
            /* determine order of three atoms correctly in order
                   to avoid creating a clockwise concave face */
            /* allocate new probe */
            if (ept -> begin)
                prb = dsnwprb (atom1, atom2, atom3, center);
            else
                prb = dsnwprb (atom2, atom1, atom3, center);
            if (errflg)
                return;
            prb -> low = (dshei (prb) < pradius);
            if (errflg)
                return;

            if (connected) {
            /* get rid of torus */
                dsfcir (cir);
                if (errflg)
                    return;
                free_object (TORUS, (short *) tor);
                if (ntori != 1) {
                    dserr (855, "dsplace");
                    return;
                }
                ntori = 0;
                hedtor = NULL;
                taltor = NULL;
                /* connected rolling will generate all tori */
                return;
            }
        }                       /* end of ept loop */
    }                           /* end of torus loop */
    return;
}

void dsplin ()                        /* probe list initialization */
{
    /* this routine is used for connected rolling only */
    int     a1, a2, a3, hash;
    struct probe   *prb;
    struct plist   *pl;

    /* set up hashing table */
    nplist = 3 * natom;
    hedpl = (struct plist **) allocate_pointers (PLIST, nplist);
    if (hedpl == NULL) {
        dserr (500, "dsplin");
        return;
    }
 /* put first probe in hashing table */
    prb = hedprb;
    a1 = prb -> atom[0];
    a2 = prb -> atom[1];
    a3 = prb -> atom[2];
    atom = a1;
    hash = a1 + a2 + a3;                /* sophisticated hashing key */
    if (hash < 0 || hash >= nplist) {
        dserr (665, "dsplin");
        return;
    }
    /* allocate head of probe list for this hashing key */
    pl = (struct plist *) allocate_object (PLIST);
    if (pl == NULL) {
        dserr (500, "dsplin");
        return;
    }
    pl -> prb = prb;                    /* store pointer to probe */
    pl -> next = NULL;                  /* mark as tail of list */
    *(hedpl + hash) = pl;               /* store head pointer */
}

struct probe *dsplook (int a1, int a2, int a3)    /* look up probe in table */
{
    int    k, hash, l1, l2, l3;
    double   center[3];
    struct probe   *prb;
    struct plist   *pl;


 /* hash key must be independent of atom order */
    atom = a1;
    hash = a1 + a2 + a3;
    if (hash < 0 || hash >= nplist) {
        dserr (665, "dsplook");
        return (NULL);
    }

    /* go down the linked list for this hash key */
    for (pl = *(hedpl + hash); pl != NULL; pl = pl -> next) {
        prb = pl -> prb;
        l1 = prb -> atom[0];
        l2 = prb -> atom[1];
        l3 = prb -> atom[2];
        /* check the three cyclic permutations */
        if (a1 == l1 && a2 == l2 && a3 == l3
                || a1 == l2 && a2 == l3 && a3 == l1
                || a1 == l3 && a2 == l1 && a3 == l2)
                                /* found probe in list */
            return (prb);
    }

 /* a probe with these three atoms and the same orientation
    has not been found, so create new probe */

    for (k = 0; k < 3; k++)
        center[k] = 0.0;
    prb = dsnwprb (a1, a2, a3, center);
    if (errflg)
        return (NULL);

 /* new plist entry */
    pl = (struct plist *) allocate_object (PLIST);
    if (pl == NULL) {
        dserr (500, "dsplook");
        return (NULL);
    }
    pl -> prb = prb;                    /* store pointer to probe */
    pl -> next = *(hedpl + hash);       /* insert at head of linked list */

 /* leave center, links, prb->low alone (caller sets them) */

    *(hedpl + hash) = pl;               /* store new head pointer */
    return (prb);
}

/* plop probe down at center and create new probe placement */
struct probe *dsplop (struct probe *frmprb, double center[3], int  a1, int a2, int a3)
{
    int    new, k;
    struct probe   *prb;

    atom = a1;
    prb = dsplook (a1, a2, a3);
    if (errflg)
        return (NULL);

    if (prb == NULL) {
        dserr(750, "dsplop");
        return(NULL);
    }

    /* mark whether we've been here before */
    new = (prb -> link[0] == NULL && prb -> link[1] == NULL &&
            prb -> link[2] == NULL);

    /* check whether dsplook gaves us what we asked for */
    if (new && (prb -> atom[0] != a1 || prb -> atom[1] != a2 ||
                prb -> atom[2] != a3)) {
        dserr (765, "dsplop");
        return (NULL);
    }

 /* store geometric info for new probe,
             also pointer to from where we came */
    if (new) {
        for (k = 0; k < 3; k++)
            prb -> center[k] = center[k];
        prb -> low = (dshei (prb) < pradius);
        if (errflg)
            return (NULL);
        prb -> link[0] = frmprb;
    }
    else {
        /* store pointer to probe we came from in appropriate link
           element, and make sure it has not already been used */
        if (a1 == prb -> atom[1] && a2 == prb -> atom[2])
            if (prb -> link[1] != NULL) {
                dserr (672, "dsplop");
                return (NULL);
            }
            else
                prb -> link[1] = frmprb;
        else
            if (a1 == prb -> atom[2] && a2 == prb -> atom[0])
                if (prb -> link[2] != NULL) {
                    dserr (672, "dsplop");
                    return (NULL);
                }
                else
                    prb -> link[2] = frmprb;
            else {
                dserr (622, "dsplop");
                return (NULL);
            }
    }
    return (prb);
}

struct cluster *dsprb (struct probe *prb)     /* compute surface points for probe */
{
    int    k, i, j, nlat;
    int     a1, a2, a3, kn, nnl;
    double   z, rad, dt1, dt2, dt3, mindot, angle, width, p2, d2;
    double   den1, den2, den3;
    double   axis[3], center[3];
    double   central[3], zenith[3], vectl[3];
    double   v1[3], v2[3], v3[3];
    double   norms[3][3];
    double  *f;
    struct circle  *cir;
    struct endpnt  *ept1, *ept2;
    struct cluster *merclu, *clu, *pclu, *fclu, *prbclu;
    /* structures and pointers for handling cusps */
    struct probe   *prb2;
    struct probe  **nearlow;
    struct midpln  *midptr;
    struct midpln  *nearmid;

    nnl = 0;
 /* if low probe, form list of near low probes */
    if (prb -> low) {
        nearlow = (struct probe **) allocate_pointers (PROBE, nlow);
        if (nearlow == NULL) {
            dserr (500, "dsprb");
            return (NULL);
        }
        p2 = 4 * pradius * pradius;
        for (i = 0; i < nlow; i++) {
            prb2 = *(lowp + i);
            if (prb == prb2)
                continue;
            d2 = dsdis2 (prb -> center, prb2 -> center);
            if (d2 < p2)
                *(nearlow + nnl++) = prb2;
        }
    }

    /* if this is a low probe and there are other nearby low probes
       calculate a bisecting plane (midplane) for each nearby low probe */
    if (nnl > 0) {
        nearmid = (struct midpln   *) allocate_objects (MIDPLN, nnl);
        if (nearmid == NULL) {
            dserr (500, "dsprb");
            return (NULL);
        }
        for (j = 0; j < nnl; j++) {
            midptr = nearmid + j;       /* pointer to storage space */
            prb2 = (*(nearlow + j));    /* pointer to nearby low probe */
            /* center and axis: */
            for (k = 0; k < 3; k++) {
                midptr -> center[k] = (prb2 -> center[k] +
                        prb -> center[k]) / 2;
                midptr -> axis[k] = prb2 -> center[k] -
                    prb -> center[k];
            }
            /* midpoint cannot be too far away from central probe center */
            if (dsdis (midptr -> center, prb -> center) > pradius) {
                dserr (869, "dsprb");
                return (NULL);
            }
            /* normalize plane axis vector */
            if (!dsnize (midptr -> axis)) {
                dserr (800, "dsprb");
                return (NULL);
            }
        }
    }


 /* set up coordinates, vectors for probe */
    a1 = prb -> atom[0];
    a2 = prb -> atom[1];
    a3 = prb -> atom[2];
    atom = a1;
    /* use average surface point density of three atoms */
    den1 = *(atmden + a1);
    den2 = *(atmden + a2);
    den3 = *(atmden + a3);
    density = (den1 + den2 + den3) / 3;
    /* vectors from probe center to three atom centers define triangle */
    for (k = 0; k < 3; k++) {
        v1[k] = *(atmco + 3 * a1 + k) - prb -> center[k];
        v2[k] = *(atmco + 3 * a2 + k) - prb -> center[k];
        v3[k] = *(atmco + 3 * a3 + k) - prb -> center[k];
    }
    if (!dsnize (v1)) {
        dserr (800, "dsprb");
        return (NULL);
    }
    if (!dsnize (v2)) {
        dserr (800, "dsprb");
        return (NULL);
    }
    if (!dsnize (v3)) {
        dserr (800, "dsprb");
        return (NULL);
    }

    /* compute normal vectors to the 3 planes that cut out the face */
    dscross (v1, v2, norms[0]);
    if (!dsnize (norms[0])) {
        dserr (800, "dsprb");
        return (NULL);
    }

    dscross (v2, v3, norms[1]);
    if (!dsnize (norms[1])) {
        dserr (800, "dsprb");
        return (NULL);
    }

    dscross (v3, v1, norms[2]);
    if (!dsnize (norms[2])) {
        dserr (800, "dsprb");
        return (NULL);
    }

    /* check that the concave face has counter-clockwise orientation */
    if (dsdot (norms[0], v3) > 0.0) {
        dserr (868, "dsprb");
        return (NULL);
    }


    /* choose average of vertex vectors
            to define northward central vector */
    for (k = 0; k < 3; k++)
        central[k] = v1[k] + v2[k] + v3[k];
    if (!dsnize (central)) {
        dserr (800, "dsprb");
        return (NULL);
    }

 /* cone enclosing concave face with central vector as its axis */
    mindot = 1.0;               /* initialize minimum dot product */
    dt1 = dsdot (v1, central);
    dt2 = dsdot (v2, central);
    dt3 = dsdot (v3, central);
    if (dt1 < mindot)
        mindot = dt1;
    if (dt2 < mindot)
        mindot = dt2;
    if (dt3 < mindot)
        mindot = dt3;
    if (mindot < -1.0 - (10.0 * EPSILON) || mindot > 1.0 + (10.0 * EPSILON)) {
        dserr (862, "dsprb");
        return (NULL);
    }
    if (mindot > 1.0)
        mindot = 1.0;
    else
        if (mindot < -1.0)
            mindot = -1.0;
    angle = acos (mindot);      /* cone angle */
    if (angle <= 0.0) {
        dserr (860, "dsprb");
        return (NULL);
    }

    /* create (meridian) arc between north pole and cone */
    ept1 = dsnwept (1, 0.0, EOL);
    if (errflg)
        return (NULL);
    ept2 = dsnwept (0, angle, EOL);
    if (errflg)
        return (NULL);
    ept1 -> next = ept2;

    dsarbp (central, axis);
    if (errflg)
        return (NULL);
    dscross (axis, central, zenith);
    cir = dsnwcir (prb -> center, pradius, axis);
    if (errflg)
        return (NULL);
    for (k = 0; k < 3; k++)
        cir -> base[k] = central[k];
    cir -> atom = EOL;
    cir -> head = ept1;

 /* subdivide meridian arc so we know where to put the latitudes */
    merclu = dsdivid (cir);
    if (errflg)
        return (NULL);
    nlat = merclu -> nmem;
    if (nlat <= 0) {
        dserr (865, "dsprb");
        return (NULL);
    }
    /* width adjusted for curvature of sphere */
    width = 2 * pradius * sin (angle / (2 * nlat));
    dsfcir (cir);       /* free circle */
    if (errflg)
        return (NULL);

 /* latitude center varies along central axis */

    f = merclu -> points;       /* these points are not surface points */
    if (f == NULL) {
        dserr (625, "dsprb");
        return (NULL);
    }

    pclu = NULL;        /* initialization for latitude cluster list */
    fclu = NULL;

    /* loop for each subdivision point of meridian arc */
    for (i = 0; i < nlat; i++, f += 3) {

    /* project meridian point onto central axis */
        for (k = 0; k < 3; k++)
            vectl[k] = *(f + k) - prb -> center[k];
        z = dsdot (central, vectl);
        /* calculate radius of latitude circle */
        rad = pradius * pradius - z * z;
        if (rad < 0.0) {
            dserr (861, "dsprb");
            return (NULL);
        }
        rad = sqrt (rad);
        if (rad <= 0.0) {
            dserr (859, "dsprb");
            return (NULL);
        }
    /* allocate latitude circle */
        for (k = 0; k < 3; k++)
            center[k] = prb -> center[k] + z * central[k];
        cir = dsnwcir (center, rad, central);
        if (errflg)
            return (NULL);
        cir -> width = width;

    /* cut with three planes */
        for (kn = 0; kn < 3; kn++) {
            dscpi (cir, prb -> center, norms[kn]);
            if (errflg)
                return (NULL);
            if (cir -> gone)
                break;
        }
        if (cir -> gone) {
            dsfcir (cir);
            if (errflg)
                return (NULL);
            continue;
        }

        if (nnl > 0) {
        /* cut with midlow planes */
            for (j = 0; j < nnl; j++) {
                midptr = nearmid + j;
                dscpi (cir, midptr -> center, midptr -> axis);
                if (errflg)
                    return (NULL);
                if (cir -> gone)
                    break;
            }
            if (cir -> gone) {
                dsfcir (cir);
                if (errflg)
                    return (NULL);
                continue;
            }
        }


        /* convert arc(s) of latitude circle into surface points */
        clu = dsdivid (cir);
        if (errflg)
            return (NULL);
        dsfcir (cir);           /* free latitude circle */
        if (errflg)
            return (NULL);
        if (clu == NULL)
            continue;
        if (clu -> nmem <= 0) {
            dsfclu (clu);
            if (errflg)
                return (NULL);
            continue;
        }

    /* link into list of latitude clusters */
        if (pclu == NULL)
            fclu = clu;
        else
            pclu -> next = clu;
        pclu = clu;
    }

    dsfclu (merclu);            /* free meridian cluster */
    if (errflg)
        return (NULL);
    if (fclu == NULL)           /* check for no latitude clusters */
        return (NULL);
    prbclu = dsmerge (fclu);    /* merge lat clusters into face cluster */
    if (errflg)
        return (NULL);
    if (prb -> low)
        free_pointers (PROBE, (void *) nearlow);        /* free memory */
    if (nnl > 0)
        free_objects (MIDPLN, (short *) nearmid);        /* free memory */

    return (prbclu);
}

void dsprbs ()       /* generate concave surface for all probe positions */
{
    int     ia;
    struct probe   *prb;
    struct probe  **pp;

 /* create array of head pointers for concave surface */
    hedcav = (struct cluster  **) allocate_pointers (CLUSTER, natom);
    if (hedcav == NULL) {
        dserr (500, "dsprbs");
        return;
    }

    /* check for van der Waals surface */
    if (pradius <= 0.0)
        return;

 /* create list of low probes */

    nlow = 0;
    for (prb = hedprb; prb != NULL; prb = prb -> next)
        if (prb -> low)
            nlow++;
    if (nlow <= 0)
        lowp = NULL;
    else {
        lowp = (struct probe  **) allocate_pointers (PROBE, nlow);
        if (lowp == NULL) {
            dserr (500, "dsprbs");
            return;
        }
        for (prb = hedprb, pp = lowp; prb != NULL; prb = prb -> next)
            if (prb -> low)
                *pp++ = prb;
    }

 /* generate and partition probe surface */
    for (prb = hedprb; prb != NULL; prb = prb -> next) {
        dspart (prb);
        if (errflg)
            return;
    }

 /* merge concave surface lists */

    for (ia = 0; ia < natom; ia++)
        if (*(hedcav + ia) != NULL) {
            atom = ia;
            *(hedcav + ia) = dsmerge (*(hedcav + ia));
            if (errflg)
                return;
        }

    return;
}

void dsputir (struct hidden *ir)                    /* store internal record info */
{
    ir -> ntori = ntori;
    ir -> nprobe = nprobe;
    ir -> nlow = nlow;
    ir -> ndead = ndead;
    ir -> hedtor = hedtor;
    ir -> taltor = taltor;
    ir -> hedprb = hedprb;
    ir -> talprb = talprb;
    ir -> lowp = lowp;
    ir -> hedrun = hedrun;
    ir -> hedcav = hedcav;
    ir -> hedpl = hedpl;
    ir -> problem = problem;
    ir -> nplist = nplist;
}




void dsstage (char *string)                /* store current stage in calculation */
{
    if (strlen (string) >= (unsigned long) 71)  /* be careful about smashing memory */
        *(string + 70) = 0;
	strcpy (stage, "         ");
    strcpy (stage+9, string);
    if (verbose)
        inform (stage);
}

void dsstep (struct probe *prb)     /* step to new probe positions from current */
{
    int     step;
    int     newlk[3];           /* 1 if new link to new probe positions */

    /* this function is recursive,
            but only the storage above is automatic */

    static int  k, ia, ja, ma, oa, ea, nmut, swapm;
    static  atomnum * mut, *m;
    static double    irad, jrad, dij, dij2, root1, root2;
    static double    cosang, sinang, radius, d1;
    static double    tp[3], zenith[3], ico[3], jco[3];
    static double    vij[3], axis[3], center[3];
    static double    ept1cen[3], ept2cen[3];
    static struct circle   *cir;
    static struct endpnt   *ept1, *ept2;
    static struct torus *tor;
    static struct probe *badprb;


    if (prb == NULL) {
        dserr (760, "dsstep");
        return;
    }

 /* check for no work */
    if (prb -> link[0] != NULL && prb -> link[1] != NULL
            && prb -> link[2] != NULL)
        return;

    /* initialize */
    for (step = 0; step < 3; step++)
        newlk[step] = 0;

    /* from each probe position there are three directions,
       corresponding to crevices between the three atoms */

    for (step = 0; step < 3; step++) {
        if (prb -> link[step] != NULL)  /* skip if link already present */
            continue;
        /* set up three atom numbers */
        switch (step) {
            case 0:
                ia = prb -> atom[0];
                ja = prb -> atom[1];
                oa = prb -> atom[2];
                break;
            case 1:
                ia = prb -> atom[1];
                ja = prb -> atom[2];
                oa = prb -> atom[0];
                break;
            case 2:
                ia = prb -> atom[2];
                ja = prb -> atom[0];
                oa = prb -> atom[1];
                break;
            default:
                dserr (658, "dsstep");
                return;
        }

        /* we roll along atoms ia (on right) and ja (on left),
           while we leave (step off) from atom oa */

        atom = oa;              /* for possible error message */

        /* if both atoms we would roll along are <= BLOCKERS
           then we mark a dead end */

        if (*(atmatt + ia) < AREA && *(atmatt + ja) < AREA) {
            /* point at pseudo probe */
            badprb = (struct probe *) allocate_object (PROBE);
            if (badprb == NULL) {
                dserr (500, "dsstep");
                return;
            }
            for (k = 0; k < 3; k++) badprb -> atom[k] = EOL;
            prb -> link[step] = badprb; /* mark as dead end */
            ndead++;                    /* increment number of dead ends */
            continue;
        }

    /* torus code: */
        irad = *(atmrad + ia) + pradius;
        for (k = 0; k < 3; k++)
            ico[k] = *(atmco + 3 * ia + k);
        jrad = *(atmrad + ja) + pradius;
        for (k = 0; k < 3; k++)
            jco[k] = *(atmco + 3 * ja + k);
        for (k = 0; k < 3; k++)
            vij[k] = jco[k] - ico[k];
        dij = dsnorm (vij);
        if (dij <= 0.0) {
            dserr (852, "dsstep");
            return;
        }
        dij2 = dij * dij;
        if (dij2 <= 0.0) {
            dserr (852, "dsstep");
            return;
        }
        for (k = 0; k < 3; k++)
            axis[k] = vij[k] / dij;
        for (k = 0; k < 3; k++)
            center[k] = 0.5 * (ico[k] + jco[k])
                + 0.5 * vij[k] * (irad * irad - jrad * jrad) / dij2;
        root1 = (irad + jrad) * (irad + jrad) - dij2;
        if (root1 < 0.0) {
            dserr (876, "dsstep");
            return;
        }
        root1 = sqrt (root1);
        root2 = dij2 - (irad - jrad) * (irad - jrad);
        if (root2 < 0.0) {
            dserr (876, "dsstep");
            return;
        }
        root2 = sqrt (root2);
        radius = 0.5 * root1 * root2 / dij;
        if (radius <= 0.0) {
            dserr (859, "dsstep");
            return;
        }
        cir = dsnwcir (center, radius, axis);
        if (errflg)
            return;

    /* get list of mutual neighbors */
        mut = dsmut (ia, ja);
        if (errflg)
            return;

        nmut = 0;
        if (mut != NULL)
            for (m = mut; *m != EOL; m++)
                nmut++;
        if (nmut < 1) {
            dserr (870, "dsstep");
            return;
        }

        /* make sure atom oa creates the base vector for circle */
        swapm = 0;
        for (m = mut; *m != EOL; m++)
            if (*m == oa) {
                *m = *mut;
                *mut = (unsigned short) oa;
                swapm = 1;
                break;
            }
        if (!swapm) {
            dserr (875, "dsstep");
            return;
        }
        /* atom oa is now the first atom in the list of mutual neighbors */

        /* determine the end of our step
           by cutting the torus central circle with the expanded sphere
           of each mutual neighbor of atoms ia and ja */

        for (m = mut; *m != EOL; m++) {
            ma = *m;
            dscsi (cir, ma);    /* cut out inaccessible region */
            if (errflg)
                return;
            if (cir -> gone)
                break;
            dslfa (cir);        /* our step involves only the first arc */
            if (errflg)
                return;
        }
        free_atomnums (mut);    /* free temporary memory */

        /* check for annoying inconsistency created by
           doubleing-point round off error */

        if (cir -> gone) {
            /* point at pseudo probe */
            badprb = (struct probe *) allocate_object (PROBE);
            if (badprb == NULL) {
                dserr (500, "dsstep");
                return;
            }
            for (k = 0; k < 3; k++) badprb -> atom[k] = EOL;
            prb -> link[step] = badprb; /* mark as dead end */
            ndead++;                    /* increment number of dead ends */
            /* mark all 3 as problem atoms */
            *(problem+ia) = 1;
            *(problem+ja) = 1;
            *(problem+oa) = 1;
            continue;
        }

        /* endpoints of step arc: */
        ept1 = cir -> head;
        ept2 = ept1 -> next;

        /* check for another kind of annoying inconsistency created by
           doubleing-point round off error */

        if (ept1 -> atom != oa) {
            /* point at pseudo probe */
            badprb = (struct probe *) allocate_object (PROBE);
            if (badprb == NULL) {
                dserr (500, "dsstep");
                return;
            }
            for (k = 0; k < 3; k++) badprb -> atom[k] = EOL;
            prb -> link[step] = badprb; /* mark as dead end */
            ndead++;                    /* increment number of dead ends */
            /* mark all 3 as problem atoms */
            *(problem+ia) = 1;
            *(problem+ja) = 1;
            *(problem+oa) = 1;
            continue;
        }

        /* this kind of inconsistency should not happen */
        if (ept1 -> angle != 0.0) {
            dserr (882, "dsstep");
            return;
        }

    /* we have a torus central circle that is partially accessible */
        tor = dsnwtor (ia, ja, cir);
        if (errflg)
            return;
        if (tor == NULL) {
            dserr (750, "dsstep");
            return;
        }

        /* check that the arc begins at the initial probe center */
        for (k = 0; k < 3; k++)
            tp[k] = cir -> radius * cir -> base[k];
        for (k = 0; k < 3; k++)
            ept1cen[k] = cir -> center[k] + tp[k];
        d1 = dsdis (prb -> center, ept1cen);
        if (fabs (d1) > (10.0 * EPSILON)) {
            dserr (882, "dsstep");
            return;
        }

        dscross (cir -> axis, cir -> base, zenith);

    /* calculate new probe center at new arc endpoint */
        cosang = cos (ept2 -> angle);
        sinang = sin (ept2 -> angle);
        for (k = 0; k < 3; k++)
            tp[k] = cir -> radius * (cosang * cir -> base[k]
                    + sinang * zenith[k]);
        for (k = 0; k < 3; k++)
            ept2cen[k] = cir -> center[k] + tp[k];

        ea = ept2 -> atom;      /* atom we bumped into */

        /* create new probe and store link (pointer)
                        from old probe to new probe */

        prb -> link[step] = dsplop (prb, ept2cen, ja, ia, ea);
        if (errflg)
            return;
        if (prb -> link[step] == NULL) {
            dserr (624, "dsstep");
            return;
        }
        newlk[step] = 1;        /* mark as new link */
    }                           /* end of step loop */

    /* for each new link to a new probe position, step off from there */
    for (step = 0; step < 3; step++) {
        if (newlk[step])
            dsstep (prb -> link[step]);         /* recursion */
        if (errflg)
            return;
    }
}



/* generate surface for toroidal face */
struct cluster *dstor (struct torus *tor, int which)
{
    int    k;
    int     ncen, a, ic, ao, n, i, aside;
    unsigned    size;
    double   np2, delta, factd, factt, dtq, dt, erad, erado, width;
    double   cen[3], av[3], pcen[3], pnt[3], thv[3];
    double   p1[3], p2[3], p3[3], q[3];
    double   midvec[3], pco[3], tco[3], vects[3];
    struct circle   arccir, bandcir;
    struct endpnt   ept1, ept2;
    double  *f, *fn, *fd;
    struct cluster *arcclu, *fclu, *clu, *pclu, *totclu;
    struct circle  *torcir;

    /* check for van der Waals surface */
    if (pradius <= 0.0)
        return (NULL);

    /* put info into local variables */
    torcir = tor -> central;
    for (k = 0; k < 3; k++)
        tco[k] = torcir -> center[k];
    a = tor -> atom[which];             /* atom we create surface for */
    ao = tor -> atom[1 - which];        /* other atom of torus */
    atom = a;
    density = *(atmden + a);
    erad = pradius + *(atmrad + a);
    if (erad <= 0.0) {
        dserr (859, "dstor");
        return(NULL);
    }
    erado = pradius + *(atmrad + ao);
    if (erado <= 0.0) {
        dserr (859, "dstor");
        return(NULL);
    }

    /* create arbitrary probe position along torus central circle */
    dsarbp (torcir -> axis, av);
    if (errflg)
        return (NULL);
    for (k = 0; k < 3; k++)
        pcen[k] = torcir -> center[k] + torcir -> radius * av[k];

    /* calculate unit vectors from probe to atoms */
    for (k = 0; k < 3; k++) {
        p2[k] = (*(atmco + 3 * a + k) - pcen[k]) / erad;
        p3[k] = (*(atmco + 3 * ao + k) - pcen[k]) / erado;
    }
    aside = (dsdot(p2,torcir->axis) * dsdot(p3,torcir->axis) > 0.0);
    /* check that it really is a unit vector */
    np2 = dsnorm (p2);
    if (fabs (np2 - 1.0) > 10.0 * EPSILON) {
        dserr (850, "dstor");
        return (NULL);
    }

    /* calculate vector midway between two contact circles */
    for (k = 0; k < 3; k++)
        midvec[k] = p2[k] + p3[k];
    if (!dsnize (midvec)) {
        dserr (800, "dstor");
        return (NULL);
    }
    /* make p2 vector reach to contact circle */
    for (k = 0; k < 3; k++)
        p2[k] = p2[k] * pradius;

 /* calculate p1 */
    if (!tor -> cusp || aside)           /* calculate vector to midline */
        for (k = 0; k < 3; k++)
            p1[k] = midvec[k] * pradius;

    else {                      /* calculate vector to cusp point */
        dtq = pradius * pradius - torcir -> radius * torcir -> radius;
        if (dtq < 0.0) {
            dserr (861, "dstor");
            return(NULL);
        }
        dtq = sqrt (dtq);
        for (k = 0; k < 3; k++)
            q[k] = torcir -> axis[k] * dtq * (2 * which - 1);
        for (k = 0; k < 3; k++)
            p1[k] = torcir -> center[k] + q[k] - pcen[k];
    }

 /* arc from p1 to p2 */
    for (k = 0; k < 3; k++)
        arccir.center[k] = pcen[k];
    dscross (torcir -> axis, av, arccir.axis);
    if (!which)
        for (k = 0; k < 3; k++)
            arccir.axis[k] = -arccir.axis[k];
    for (k = 0; k < 3; k++)
        arccir.base[k] = p1[k] / pradius;
    arccir.radius = pradius;
    arccir.width = 0.0;
    arccir.atom = (unsigned short) a;
    arccir.head = &ept1;
    arccir.gone = 0;
    ept1.angle = 0.0;
    ept1.begin = 1;
    ept1.atom = 0;
    ept1.next = &ept2;
    ept2.begin = 0;
    ept2.atom = (unsigned short) a;
    ept2.next = NULL;
    /* calculate angle between p1 and p2 */
    dt = dsdot (p1, p2) / (pradius * pradius);
    if (dt > 1.0 + 10.0 * EPSILON || dt < -1.0 - 10.0 * EPSILON) {
        dserr (862, "dstor");
        return (NULL);
    }
    if (dt > 1.0) dt = 1.0;
    if (dt < -1.0) dt = -1.0;
    ept2.angle = acos (dt);
    if (ept2.angle <= 0.0) {
        dserr (860, "dstor");
        return (NULL);
    }

 /* subdivide arc to generate points for new circles */
    arcclu = dsdivid (&arccir);
    if (errflg)
        return (NULL);
    if (arcclu == NULL)
        return (NULL);
    if (arcclu -> nmem <= 0) {
        dsfclu (arcclu);
        if (errflg)
            return (NULL);
        return (NULL);
    }

    /* half-width of toroidal bands */
    ncen = arcclu -> nmem;
    delta = ept2.angle / (ncen * 2);
    if (delta <= 0.0) {
        dserr (860, "dstor");
        return(NULL);
    }
    /* factor for area calculation */
    factd = 1.0 - sin (delta) / delta;

    fd = arcclu -> points;

    /* initialize linked list of clusters */
    pclu = NULL;
    fclu = NULL;

 /* loop for torus surface circles */
    for (ic = 0; ic < ncen; ic++) {
    /* coordinates relative to torus center */
        for (k = 0; k < 3; k++)
            pnt[k] = *(fd + 3 * ic + k) - torcir -> center[k];
        dt = dsdot (torcir -> axis, pnt);
        for (k = 0; k < 3; k++)
            cen[k] = torcir -> axis[k] * dt;
        bandcir.radius = dsdis (cen, pnt);
        if (bandcir.radius > torcir -> radius) {
            dserr (863, "dstor");
            return (NULL);
        }
        if (bandcir.radius <= 0.0) {
            dserr (859, "dstor");
            return (NULL);
        }

        /* there now follows some rather obscure calculations
           whose purpose it is to make the band area computation
           more accurate than simply multiplying width by arc length */

    /* get cosine theta */
        for (k = 0; k < 3; k++)
            thv[k] = pcen[k] - *(fd + 3 * ic + k);
        if (!dsnize (thv)) {
            dserr (800, "dstor");
            return (NULL);
        }
        dt = dsdot (av, thv);
        if (dt <= - (10.0 * EPSILON) || dt >= 1.0 + 10.0 * EPSILON) {
            dserr (862, "dstor");
            return (NULL);
        }
        /* correction factor */
        factt = 1.0 + (pradius / bandcir.radius) * dt * factd;
        width = 2 * delta * pradius * factt;
        bandcir.width = width;
        bandcir.atom = (unsigned short) a;
        for (k = 0; k < 3; k++) {
            bandcir.center[k] = torcir -> center[k] + cen[k];
            bandcir.axis[k] = torcir -> axis[k];
            bandcir.base[k] = torcir -> base[k];
        }
        bandcir.gone = 0;

     /* This is important: use same arcs as torus central circle.
        It works, because endpoints are specified by angle, not coordinates. */

        bandcir.head = torcir -> head;

        clu = dsdivid (&bandcir);       /* generate surface points */
        if (errflg)
            return (NULL);
        if (clu == NULL)
            continue;
        if (clu -> nmem <= 0) {
            dsfclu (clu);
            if (errflg)
                return (NULL);
            continue;
        }

    /* handle surface normals if requested */
        clu -> normals = NULL;
        n = clu -> nmem;
        if (*(atmatt + a) >= NORMALS && n > 0) {
            size = n * 3 * sizeof (double);
            clu -> normals = allocate_doubles (n * 3, 0, NORMALS);
            if (clu ->normals == NULL) {
                dserr (500, "dstor");
                return (NULL);
            }
            fn = clu -> normals;
            f = clu -> points;
            /* we calculate the probe center to get the normal vector */
            for (i = 0; i < n; i++)
                for (k = 0; k < 3; k++) {
                    vects[k] = (*(f + 3 * i + k) - bandcir.center[k])
                        / bandcir.radius;
                    pco[k] = torcir -> center[k] + torcir -> radius * vects[k];
                    *fn++ = (pco[k] - *(f + 3 * i + k)) / pradius;
                }
        }

    /* link into cluster list */
        if (pclu == NULL)
            fclu = clu;
        else
            pclu -> next = clu;
        pclu = clu;
    }

    dsfclu (arcclu);                    /* through with arc cluster */

    if (errflg)
        return (NULL);
    if (fclu == NULL)
        return (NULL);

    totclu = dsmerge (fclu);            /* merge clusters */
    if (errflg)
        return (NULL);
    return (totclu);
}

/* triple product of vectors */
double dstrip (double a[3], double b[3], double c[3])
{
    double   t;
    double   ab[3];

    dscross (a, b, ab);
    t = dsdot (ab, c);
    return (t);
}



/* unpack runs of atom numbers into list of atom numbers */
atomnum *dsupk (struct run *first)
{
    struct run *r;
    int     nn, i;
    atomnum * nbrs, *ptr;

    if (first == NULL)
        return (NULL);
    nn = dslen (first);         /* number of atoms in list of runs */
    if (errflg)
        return (NULL);
    if (nn <= 0)
        return (NULL);

    /* allocate memory for atom number array */
    nbrs = allocate_atomnums ((unsigned) (nn + 1));
    if (nbrs == NULL) {
        dserr (500, "dsupk");
        return (NULL);
    }


    for (ptr = nbrs, r = first; r != NULL; r = r -> next)
        for (i = r -> atom[0]; i <= r -> atom[1]; i++)
            *ptr++ = (unsigned short) i;         /* store atom number into array */

    *ptr = EOL;                 /* mark end of array */

 /* check whether nbrs array is sorted */
    if (!dsans (nbrs)) {
        dserr (685, "dsupk");
        return (NULL);
    }

    return (nbrs);
}

/* 

Copyright 1986 by Michael L. Connolly
All Rights Reserved

*/
