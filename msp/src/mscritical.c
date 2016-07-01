/* Molecular Surface Package Copyright 1995 Michael L. Connolly */
/* February 4, 2000 */

#include "ms.h"
#include "msenum.h"
#include "mslimit.h"
#include "msdefault.h"
#include "msstruct.h"
#include "msglobal.h"
#include "msfunction.h"



/* mscritical */

long id_critical (struct surface *cpoly)
{
	long v, nc;
	struct phnvtx *cv;
	
	set_opposite (cpoly);
	set_link (cpoly);
	id_crits (cpoly);
	/* return critical point types */
	nc = 0;
	for (v = 0; v < cpoly -> n_phnvtx; v++) {
		cv = num2phnvtx (cpoly, v + 1);
		if (cv == NULL)  {
			set_error1 ("id_critical: invalid vertex number");
			return (0);
		}
		if (cv -> critter != 0) nc++;
	}
	return (nc);
}

void set_opposite (struct surface *cpoly)
{
	struct phntri *ctri;
	struct phnedg *ce, *ce1;
	long e, e1, ue, ue1, v, nt, t;
	int j, j1, orn, orn1;

	nt = cpoly -> n_phntri;
	for (t = 0; t < nt; t++) {
		ctri = num2phntri (cpoly, t + 1);
		if (ctri == NULL){
			set_error1 ("set_opposite: invalid triangle number");
			return;
		}
		for (j = 0; j < 3; j++) {
			e = ctri -> edgnum[j];
			j1 = (j < 2) ? j+1 : 0;
			e1 = ctri -> edgnum[j1];
			ue = abs(e); ue1 = abs(e1);
			orn = (e > 0) ? 0 : 1;
			orn1 = (e1 > 0) ? 0 : 1;
			ce1 = num2phnedg (cpoly, ue1);
			if (ce1 == NULL) {
				set_error1 ("set_opposite: invalid edge number");
				return;
			}
			ce = num2phnedg (cpoly, ue);
			if (ce == NULL) {
				set_error1 ("set_opposite: invalid edge number");
				return;
			}
			v = ce1 -> vtxnum[1-orn1];
			ce -> opposites[orn] = v;
		}
	}
}

void set_link (struct surface *cpoly)
{
	int j;
	long e, ne, v;
	struct phnedg *ce;
	struct phnvtx *cv;
	
	ne = cpoly -> n_phnedg;
	for (e = 0; e < ne; e++) {
		ce = num2phnedg (cpoly, e + 1);
		if (ce == NULL) {
			set_error1 ("set_link: invalid edge number");
			return;
		}
		for (j = 0; j < 2; j++) {
			v = ce -> opposites[j];
			cv = num2phnvtx (cpoly, v);
			if (cv == NULL) {
				set_error1 ("set_link: invalid vertex number");
				return;
			}
			add_link (cv, ce, j);
		}
	}
}

void add_link (struct phnvtx *cv, struct phnedg *ce, int orn)
{
	struct critlink *cl;
	
	cl = (struct critlink *) allocate_object (CRITLINK);
	if (cl == NULL) {
		set_error1 ("add_link: memory allocation fails");
		return;
	}
	cl -> next = cv -> head;
	cv -> head = cl;
	cl -> edg = ce;
	cl -> orn = orn;
}

void id_crits (struct surface *cpoly)
{
	struct phnvtx *cv;
	long nv, v;
	
	nv = cpoly -> n_phnvtx;
	for (v = 0; v < nv; v++) {
		cv = num2phnvtx (cpoly, v + 1);
		if (cv == NULL) {
			set_error1 ("id_crits: invalid vertex number");
			return ;
		}
		id_crit (cv, cpoly);
	}
}

void id_crit (struct phnvtx *cv, struct surface *cpoly)
{
	int nchanges, c, found, ismin, kiss;
	double valued[MAX_DEGREE], aval;
	int d, idx, degree, orn;
	long vtx0, vtx1, pvtx;
	struct critlink *cl;
	struct phnedg *ce;
	struct phnvtx *cp;
	long unsorted[MAX_DEGREE][2], sorted[MAX_DEGREE];
	int used[MAX_DEGREE];
	
	/* compute degree */

	degree = 0;
	for (cl = cv -> head; cl != NULL; cl = cl -> next)
		degree++;
	if (degree > MAX_DEGREE) {
		set_error1 ("id_crit: degree too large");
		return;
	}
	for (cl = cv -> head, idx = 0; cl != NULL; cl = cl -> next, idx++) {
		ce = cl -> edg; orn = cl -> orn;
		vtx0 = ce -> vtxnum[orn];
		vtx1 = ce -> vtxnum[1-orn];
		unsorted[idx][0] = vtx0;
		unsorted[idx][1] = vtx1;
	}
	for (idx = 0; idx < degree; idx++)
		used[idx] = 0;
	pvtx = 0;
	kiss = 0;
	for (d = 0; d < degree; d++) {
		/* find unused, connected place */
		found = 0;
		for (idx = 0; idx < degree; idx++) {
			if (used[idx]) continue;
			if (pvtx == 0) {
				found = 1;
				break;
			}
			if (pvtx == unsorted[idx][0]) {
				found = 1;
				break;
			}
		}
		if (!found) {
			inform ("surface(s) kiss at vertex");
			kiss = 1;
			break;
		}
		sorted[d] = unsorted[idx][0];
		cp = num2phnvtx(cpoly, sorted[d]);
		if (cp == NULL) {
			set_error1 ("id_crit: invalid vertex number");
			return;
		}

		valued[d] = cp -> values[0] - cv -> values[0];
		pvtx = unsorted[idx][1];
		used[idx] = 1;
	}
	if (kiss) c = 4;
	else {
		/* vertex numbers are now sorted around the link */
		nchanges = count_changes (degree, valued);
		aval = average (degree, valued);
		switch (nchanges) {
		case 0:
			ismin = (aval > 0.0);
			c = ismin ? 1 : 3;
			break;
		case 2:
			c = 0;	/* gradient */
			break;
		default:
			c = 2;	/* saddle */
		}
	}
	cv -> critter = (short) c;
}

int count_changes (int degree, double valued[])
{
	int nc, d, e;
	double val0, val1;
	
	nc = 0;
	/* later: deal with val = 0.0 case */
	for (d = 0; d < degree; d++) {
		e = (d < degree - 1) ? d + 1 : 0;
		val0 = valued[d]; val1 = valued[e];
		if (val0 * val1 < 0.0) nc++;
	}
	return (nc);
}

double average (int degree, double valued[])
{
	double aval;
	int d;
	
	aval = 0.0;
	for (d = 0; d < degree; d++)
		aval += valued[d];
	aval /= degree;
	return (aval);
}

long identify (struct surface *msphn)
{
	long nc;
	long v;
	char wmessage[MAXLINE];
	struct phnvtx *pv;
	
	nc = id_critical (msphn);
	
	sprintf (wmessage, "%8ld critical points", nc);
	inform(wmessage);

	for (v = 0; v < msphn -> n_phnvtx; v++) {
		pv = num2phnvtx (msphn, v + 1);
		if (pv == NULL) {
			set_error1 ("identify: invalid vertex number");
			return (0);
		}
		pv -> hue = pv -> critter + 1;
	}
	return (nc);
}


