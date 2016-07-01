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

void dsbtree ()                      /* create box tree */
{
    int     a, k;
    double   v, maxw, mid;
    double   bounds[3][2];
    struct box *b;

    /* determine bounds of molecule */

    for (k = 0; k < 3; k++) {
        bounds[k][0] = 1000000.0;
        bounds[k][1] = -1000000.0;
    }
    for (a = 0; a < natom; a++)
        for (k = 0; k < 3; k++) {
            if (atmatt != NULL && *(atmatt+a) <= IGNORE) continue;
            v = *(atmco + 3 * a + k) - *(atmrad + a);
            if (bounds[k][0] > v)
                bounds[k][0] = v;
            v = *(atmco + 3 * a + k) + *(atmrad + a);
            if (bounds[k][1] < v)
                bounds[k][1] = v;
        }

    /* determine the maximum width of the bounds */
    maxw = 0.0;
    for (k = 0; k < 3; k++) {
        v = bounds[k][1] - bounds[k][0];
        if (v > maxw)
            maxw = v;
    }

    /* make the bounds cubical */
    for (k = 0; k < 3; k++) {
        mid = (bounds[k][0] + bounds[k][1]) / 2;
        bounds[k][0] = mid - maxw / 2;
        bounds[k][1] = mid + maxw / 2;
    }

    /* set up minimum width from defined parameter */
    minwid = 5.0;

    /* set up the root of the tree */
    rootbox = (struct box  *) allocate_object (BOX);
    if (rootbox == NULL) {
        dserr (500, "dsbtree");
        return ;
    }
    rootbox -> leaf = (maxw / 2 < minwid);
    for (k = 0; k < 3; k++) {
        rootbox -> bounds[k][0] = bounds[k][0];
        rootbox -> bounds[k][1] = bounds[k][1];
    }

    /* put each atom of the molecule into the box */
    for (a = 0; a < natom; a++) {
        if (atmatt != NULL && *(atmatt+a) <= IGNORE) continue;
        /* this recursive function returns the leaf accepting this atom */
        b = dspaib (a, rootbox);
        if (errflg)
            return;
        if (b == NULL) {
            dserr (750, "dsbtree");
            return;
        }
    }
}


void dsfbox (struct box *b)                      /* free box (recursive) */
{
    int    j;
    struct box *cb;

    if (b == NULL)
        return;
    if (b -> leaf) {
        free_atomnums (b -> atom);  /* free atom list */
    }
    else {
        for (j = 0; j < 8; j++) {
            cb = b -> child[j];
            if (cb != NULL)
                dsfbox (cb);    /* recursion */
            if (errflg)
                return;
        }
    }
    free_object (BOX, (short *) b);          /* free this box */
}

atomnum *dsgair(double rect[3][2], struct box *b)      /* get atoms in rectangular range */
{
    int    nin, nout, a, i, j;
    atomnum * lin, *lout, *lptr;
    atomnum * chlist[8];
    struct box *cb;

    if (b == NULL)
        return (NULL);

    /* check whether range intersects box */
    if (!dsrib (rect, b -> bounds))
        return (NULL);
    if (errflg)
        return (NULL);

    if (b -> leaf) {
        lin = b -> atom;
        if (lin == NULL)
            return (NULL);
        atom = *lin;
        /* count number of atoms in leaf list */
        for (nin = 0; *(lin + nin) != EOL; nin++);
        if (nin <= 0)
            return (NULL);
        /* allocate memory for output list to be the same size, initially */
        lout = allocate_atomnums ((unsigned) nin + 1);
        if (lout == NULL) {
            dserr (500, "dsgair");
            return (NULL);
        }
        lptr = lout;
        /* form the output list from the input list */
        for (i = 0; i < nin; i++) {
            a = *(lin + i);
            if (dsair (rect, a))        /* if atom in range, add to out list */
                *lptr++ = (unsigned short) a;
            if (errflg)
                return (NULL);
        }
        *lptr = EOL;                    /* mark End Of List */
        nout = lptr - lout;             /* number of atoms in output list */
        if (nout <= 0) {
            free_atomnums (lout);
            return (NULL);
        }
        else
            if (nout < nin) {
                /* reallocate less memory */
		lout = reallocate_atomnums (nout + 1, lout);
                if (lout == NULL) {
                    dserr (500, "dsgair");
                    return (NULL);
                }
            }
    }

    else {                      /* not a leaf */
        /* gather atom lists for each non-null child */
        for (j = 0; j < 8; j++) {
            cb = b -> child[j];
            if (cb == NULL)
                chlist[j] = NULL;
            else
                chlist[j] = dsgair (rect, cb);
            if (errflg)
                return (NULL);
        }
        /* merge the child lists into one */
        lout = dsmrga (chlist);
        if (errflg)
            return (NULL);
    }

    return (lout);
}

struct box *dsmkbox (double   bounds[3][2], int i)         /* make box for given octant */
{
    int     k;
    int     is[3];
    double   h;
    struct box *nb;

    /* convert octant number (0-7) to three numbers (each 0 or 1) */
    for (k = 0; k < 3; k++) {
        is[k] = i % 2;
        i /= 2;
    }
    h = (bounds[0][1] - bounds[0][0]) / 2;      /* half-width */
    /* allocate new box */
    nb = (struct box   *) allocate_object (BOX);
    if (nb == NULL) {
        dserr (500, "dsmkbox");
        return (NULL);
    }
    /* store new bounds */
    for (k = 0; k < 3; k++) {
        nb -> bounds[k][0] = bounds[k][0] + h * is[k];
        nb -> bounds[k][1] = bounds[k][1] + h * (is[k] - 1);
    }
    nb -> leaf = (h / 2 < minwid);      /* if small enough, mark as leaf */
    return (nb);
}

struct box *dspaib (int a, struct box *b)       /* put atom in box */
{
    int    k, i;
    double  h;
    int     powers[3];
    struct box *cb, *ib;

    atom = a;
    if (b == NULL) {
        dserr (760, "dspaib");
        return (NULL);
    }


 /* check whether atom is contained in box bounds */
    if (!dsair (b -> bounds, a)) {
        dserr (853, "dspaib");
        return (NULL);
    }

    if (b -> leaf) {
        /* add atom to list for this box */
        b -> atom = dsaa (b -> atom, a);
        if (errflg)
            return (NULL);
        return (b);             /* return pointer to this box */
    }

    /* box is not a leaf, atom must be added to one of 8 children */

    /* half-width of box */
    h = (b -> bounds[0][1] - b -> bounds[0][0]) / 2;

    /* determine octant */
    powers[0] = 1;
    powers[1] = 2;
    powers[2] = 4;
    i = 0;
    for (k = 0; k < 3; k++)
        if (*(atmco + 3 * a + k) > b -> bounds[k][0] + h)
            i += powers[k];

    cb = b -> child[i];         /* this is the child we want */
    /* if the box does not yet exist, create it */
    if (cb == NULL) {
        cb = dsmkbox (b -> bounds, i);
        if (errflg)
            return (NULL);
        if (cb == NULL) {
            dserr (750, "dspaib");
            return (NULL);
        }
        b -> child[i] = cb;     /* store pointer */
    }

    ib = dspaib (a, cb);        /* recursion */
    if (errflg)
        return (NULL);
    return (ib);                /* return pointer to eventual recipient box */
}

atomnum *dsaa (atomnum *list, int a)        /* add atom to list */
{
    int     i, n;
    atomnum * newl;

    atom = a;
    if (list == NULL) {         /* new list */
        newl = (atomnum *) allocate_atomnums ((unsigned) 2);
        if (newl == NULL) {
            dserr (500, "dsaa");
            return (NULL);
        }
        n = 0;
        i = -1;
    }
    else {                      /* old list */
        for (n = 0; *(list + n) != EOL; n++);
	newl = reallocate_atomnums (n + 2, list);
        if (newl == NULL) {
            dserr (500, "dsaa");
            return (NULL);
        }
        /* shift higher numbers upward */
        for (i = n - 1; i >= 0; i--)
            if (*(newl + i) > a)
                *(newl + i + 1) = *(newl + i);
            else
                break;
    }

    *(newl + i + 1) = (unsigned short) a;        /* insert new atom number */
    *(newl + n + 1) = EOL;      /* mark end of list */
    return (newl);              /* return pointer to moved list */
}

int dsair (double rect[3][2], int a)                 /* return whether atom is inside range */
{
    int     k;
    double  *f, v, d;

    atom = a;
    f = atmco + 3 * a;


    for (k = 0; k < 3; k++) {
        v = *(f + k);
		d = (v - rect[k][0]);
		if (fabs (d) > EPSILON * EPSILON && d < 0.0) return (0);
        /* if (v < rect[k][0]) return (0); */
		d = (v - rect[k][1]);
		if (fabs (d) > EPSILON * EPSILON && d > 0.0) return (0);
        /* if (v > rect[k][1]) return (0); */
    }
    return (1);

}


void dserr (int errno, char *routine)          /* fatal error function */
{
    int     lenrout, lenmsg, lentot, idx;
    char   *totstr, *msg;
	char un[100];


    /* search table for error number */
    for (idx = 0; idx < ntab; idx++)
        if (errtab[idx].number == errno)
            break;
    if (idx >= ntab) {
		sprintf (un, "unidentifiable error %d", errno);
		msg = un;
	}
    else
        msg = errtab[idx].msg;

    /* compose an error string */
    lenmsg = strlen (msg);
    lenrout = strlen (routine);
    lentot = lenmsg + lenrout + 10;
    totstr = allocate_chars (lentot);
    sprintf (totstr, "%s: %s", routine, msg);
    errstr = totstr;
    errflg = errno;
    /* clean up temporary memory (unless error occurred there) */
    if (strcmp (routine, "dsclean") != 0)
        dsclean ();
    return;
}


atomnum *dsmrga (atomnum *eight[8])        /* merge eight atom lists */
{
    int     emin, amin, value, e, l, totlen;
    int     len[8];
    unsigned    size;
    atomnum * tot, *tptr;
    atomnum * eptr[8];

    /* determine sum of eight lengths */
    totlen = 0;
    for (e = 0; e < 8; e++) {
        if (eight[e] == NULL)
            len[e] = 0;
        else {
            for (l = 0; *(eight[e] + l) != EOL; l++);
            len[e] = l;
        }
        totlen += len[e];
    }

    if (totlen > 0) {
        /* allocate memory to hold new, combined atom list */
        size = (totlen + 1) * sizeof (atomnum);
	tot = allocate_atomnums (totlen + 1);
        if (tot == NULL) {
            dserr (500, "dsmrga");
            return (NULL);
        }


        tptr = tot;                     /* initialize output list pointer */
        for (e = 0; e < 8; e++)
            eptr[e] = eight[e];         /* initialize input list pointers */

    /* transfer atom numbers */
        while (tptr < tot + totlen) {
            /* look for list whose current atom is lowest */
            amin = natom;
            emin = -1;
            for (e = 0; e < 8; e++) {
                if (eptr[e] == NULL)
                    continue;
                value = *eptr[e];
                if (value != EOL && value < amin) {
                    amin = value;
                    emin = e;
                }
            }
            if (emin < 0) {
                dserr (702, "dsmrga");
                return (NULL);
            }
            /* of the eight input atom lists, emin currently has lowest */
            *tptr++ = *eptr[emin]++;
        }
        *tptr = EOL;            /* mark End Of List */
    }
    else
        tot = NULL;

 /* free input atom number lists */
    for (e = 0; e < 8; e++)
        if (eight[e] != NULL)
            free_atomnums (eight[e]);

    return (tot);
}

/* return true if range intersects box */
int dsrib (double rect[3][2], double bounds[3][2])
{
    int     k;

    /* project onto each axis and check for interval intersection */
    for (k = 0; k < 3; k++)
        if (!dsii (rect[k], bounds[k]))
            return (0);

    /* intersection in all projections implies intersection in space */
    return (1);
}


void dsclean ()                      /* cleaup temporary memory */
{
    struct run *r;
    struct run **ptr;
    struct torus   *tor, *nxttor;
    struct probe   *prb, *nxtprb;
    struct plist   *pl, *npl;
    struct plist  **ltr;

    /* array of neighbor lists (run encoded) for each atom */
    if (hedrun != NULL) {
        for (ptr = hedrun; ptr < hedrun + natom; ptr++) {
            r = *ptr;
            if (r != NULL)
                dsfrun (r);
            if (errflg)
                return;
        }
        free_pointers (RUN, hedrun);
        hedrun = NULL;
    }
	

    /* array flagging problem atoms */
    if (problem != NULL) {
        free_shorts (problem);
        problem = NULL;
    }
	

    /* array of pointers to clusters of concave surface points */
    if (hedcav != NULL) {
        free_pointers (CLUSTER, hedcav);
        hedcav = NULL;
    }

    /* probe placement lists */
    if (hedprb != NULL) {
		for (prb = hedprb, nxtprb = NULL; prb != NULL; prb = nxtprb) {
			nxtprb = prb -> next;
			free_object (PROBE, (short *) prb);
		}
		hedprb = NULL;
    }

    /* tori list */
    for (tor = hedtor, nxttor = NULL; tor != NULL; tor = nxttor) {
        nxttor = tor -> next;
                /* free torus central circle */
                if (tor -> central != (struct circle *) NULL)
                        dsfcir (tor -> central);
        free_object (TORUS, (short *) tor);
    }
    if (hedtor != NULL) hedtor = NULL;

    /* probe hashing table for connected rolling */
    if (connected && hedpl != NULL) {
        for (ltr = hedpl; ltr < hedpl + nplist; ltr++)
            for (pl = *ltr, npl = NULL; pl != NULL; pl = npl) {
                npl = pl -> next;
                free_object (PLIST, (short *) pl);
            }
        free_pointers (PLIST, hedpl);
    }
}

int dsii (double i1[2], double i2[2])                   /* intervals intersect ? */
{
    if (i1[0] > i2[1])
        return (0);
    if (i1[1] < i2[0])
        return (0);
    return (1);
}

void dsfcir (struct circle *cir)                    /* free circle */
{
    struct endpnt  *ept1, *ept2;

    if (cir == NULL)
        return;
    for (ept1 = cir -> head, ept2 = NULL; ept1 != NULL; ept1 = ept2) {
        ept2 = ept1 -> next;
        free_object (ENDPNT, (short *) ept1);
    }
    free_object (CIRCLE, (short *) cir);
}


void dsfrun (struct run *first)                  /* free run length encoded atom list */
{
    struct run *r, *next;

    if (first == NULL)
        return;
    for (r = first, next = NULL; r != NULL; r = next) {
        next = r -> next;       /* store ptr before free */
        free_object (RUN, (short *) r);
    }
}


