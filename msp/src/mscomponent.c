/*
	Molecular Surface Package
	Copyright 1986, 1989 by Michael L. Connolly
	All rights reserved
	March 6, 2000
*/

#include "ms.h"
#include "msenum.h"
#include "mslimit.h"
#include "msdefault.h"
#include "msstruct.h"
#include "msglobal.h"
#include "msfunction.h"

/* CONNECTED COMPONENTS */

/* find connected components of pq molecular surface */

void find_components (struct surface *srf)	
{
	long i, k, idx, cidx, cmax, cnum, comp;
	long *cnums, *nfc;
	char message[MAXLINE];
	struct edge *edg;
	struct face *fac;
	struct cycle *cyc;
	struct component *cmp_ptr;
	struct sphere *atm;

	srf -> n_component = 0;
	/* initialize to disconnected faces */
	i = 0;
	for (fac = srf -> head_face; fac != NULL; fac = fac -> next, i++)
		fac -> comp = i;

	/* make the edges point at faces */
	for (fac = srf -> head_face; fac != NULL; fac = fac -> next) {
		if (fac -> problem) continue;
		for (cyc = fac -> first_cycle; cyc != NULL; cyc = cyc -> next)
			for (edg = cyc -> first_edge; edg != NULL; edg = edg -> next)
				edg -> fac = fac;
	}

	/* connect the faces by means of common arcs */
	for (fac = srf -> head_face; fac != NULL; fac = fac -> next) {
		if (fac -> problem) continue;
		for (cyc = fac -> first_cycle; cyc != NULL; cyc = cyc -> next)
			for (edg = cyc -> first_edge; edg != NULL; edg = edg -> next)
				if (edg -> orn == 0) {
					same_component (edg -> arcptr);
					if (error()) return;
				}
	}

	/* allocate memory preliminary to renumbering */
	cnums = allocate_longs (srf -> n_face, 0, CNUMS);
	if (cnums == NULL) {
		set_error1("find_components: ran out of memory");
		return;
	}
	nfc = allocate_longs (srf -> n_face, 0, NFC);
	if (nfc == NULL) {
		set_error1 ("find_components: ran out of memory");
		return;
	}

	/* new component number from old, count faces per component */
	for (fac = srf -> head_face; fac != NULL; fac = fac -> next) {
		cnum = fac -> comp;
		if (*(cnums + cnum) == 0)
			*(cnums + cnum) = ++(srf -> n_component);
		idx = *(cnums + cnum) - 1;
		(*(nfc + idx))++;
	}

	/* find component with the most faces */
	cmax = 0;
	cidx = -1;
	for (idx = 0; idx < srf -> n_component; idx++)
		if (*(nfc + idx) > cmax) {
			cidx = idx;
			cmax = (*(nfc + idx));
		}
	if (cidx < 0) {
		set_error1 ("(find_components): maximum component not found");
		return;
	}

	/* sometimes biggest is not first */
	if (cidx > 0) {  /* swap */
		for (i = 0; i < srf -> n_face; i++) {
			idx = *(cnums + i) - 1;
			if (idx == 0) *(cnums + i) = cidx + 1;
			else if (idx == cidx) *(cnums + i) = 1;
		}
	}

	/* renumber components */
	for (fac = srf -> head_face; fac != NULL; fac = fac -> next) {
		cnum = fac -> comp;
		fac -> comp = *(cnums + cnum);
	}

	/* free temporary memory */
	free_longs (cnums, 0, CNUMS);
	free_longs (nfc, 0, NFC);

	sprintf (message,"%8ld components", srf -> n_component); inform(message);

	/* allocate memory for component info */
	srf -> component_handles = (struct component **)
		allocate_pointers (COMPONENT, srf -> n_component);
	if (srf -> component_handles == (struct component **) NULL) {
		set_error1 ("(find_components): ran out of memory");
		return;
	}
	for (comp = 0; comp < srf -> n_component; comp++) {
		cmp_ptr = allocate_component ();
		if (error()) {
			return;
		}
		else if (cmp_ptr == NULL) {
			set_error1 ("find_components: no mem for component");
			return;
		}
		*(srf -> component_handles + comp) = cmp_ptr;
	}

	for (comp = 0; comp < srf -> n_component; comp++) {
		cmp_ptr = get_component_ptr (srf, comp+1);
		if (cmp_ptr == NULL) {
			set_error1("find_components: get_component_ptr returns null");
			return;
		}
		cmp_ptr -> volume = 0.0;
		cmp_ptr -> area = 0.0;
		cmp_ptr -> accessible = 0.0;
		for (k = 0; k < 3; k++)
			cmp_ptr -> center[k] = 0.0;
		cmp_ptr -> count = 0;
	}

	/* compute component centers */
	for (fac = srf -> head_face; fac != NULL; fac = fac -> next) {
		if (fac -> problem) {
			continue;
		}
		comp = fac -> comp;
		cmp_ptr = get_component_ptr (srf, comp);
		switch ((int) (fac -> shape)) {
		case CONVEX:
			atm = fac -> ptr.atm;
			for (k = 0; k < 3; k++)
				cmp_ptr -> center[k] += atm -> center[k];
			cmp_ptr -> count += 1;
			break;
		default:
			break;
		}
	}

	for (comp = 0; comp < srf -> n_component; comp++) {
		cmp_ptr = get_component_ptr (srf, comp+1);
		if (cmp_ptr -> count <= 0) {
			set_error1 ("find_components: no atoms in component");
			return;
		}
		for (k = 0; k < 3; k++)
			cmp_ptr -> center[k] /= cmp_ptr -> count;
	}
}

/* faces sharing arc are in same component */
void same_component (struct arc *a)
{
	long comp0, comp1;
	struct edge *edg0, *edg1;
	struct face *fac0, *fac1, *fac;
	struct surface *srf;

	if (a == (struct arc *) NULL) {
		set_error1 ("same_component: null arc argument received");
		return;
	}
	/* move info to local variables */
	edg0 = a -> edg[0];
	edg1 = a -> edg[1];
	if (edg0 == (struct edge *) NULL) {
		set_error1 ("same_component: null edge");
		return;
	}
	if (edg1 == (struct edge *) NULL) {
		set_error1 ("same_component: null edge");
		return;
	}
	fac0 = edg0 -> fac;
	fac1 = edg1 -> fac;
	if (fac0 == (struct face *) NULL) {
		set_error1 ("same_component: null or problem face");
		return;
	}
	if (fac1 == (struct face *) NULL) {
		set_error1 ("same_component: null or problem face");
		return;
	}
	comp0 = fac0 -> comp;
	comp1 = fac1 -> comp;

	/* if two faces currently with different component numbers
	   share an arc, all faces with the higher component number
	   will be changed to the lower component number */

	srf = fac0 -> srf;
	if (srf == NULL) {
		set_error1 ("same_component: null surface");
		return;
	}
	if (comp0 < comp1) {
		for (fac = srf -> head_face; fac != NULL; fac = fac -> next)
			if (fac -> comp == comp1)
				fac -> comp = comp0;
	}
	else if (comp1 < comp0) {
		for (fac = srf -> head_face; fac != NULL; fac = fac -> next)
			if (fac -> comp == comp0)
				fac -> comp = comp1;
	}
	else return;
}
