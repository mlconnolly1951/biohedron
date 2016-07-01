/*
 * MSP
 * Copyright 1996 by Michael L. Connolly
 * All Rights Reserved

 * February 3, 2000

 */

#include "ms.h"
#include "msenum.h"
#include "mslimit.h"
#include "msdefault.h"
#include "msstruct.h"
#include "msglobal.h"
#include "msfunction.h"


/* contours */

struct surface *contour_polyhedron (struct msscene *ms, double level)
{
	int j, k, n, orn1, orn2;
	long n_ctrvtx, n_ctredg, v, e;
	double val1, val2, range, fraction;
	char message[MAXLINE];
	struct surface *obj, *po;
	struct phnvtx *vtx1, *vtx2, *ctr_vtx, *vtx;
	struct phnvtx **ctr_vertices;
	struct phnedg *edg, *ctr_edg, *next_edg, *head_phnedg;
	struct phnedg **ctr_edges;
	struct polygon *tri;
	struct phnctr *ctr;
	
	sprintf (message,"%8.3f level", level);
	informd(message);;
	if (ms -> current_molecule == (struct molecule *) NULL) {
		set_error1 ("no molecule to contour");
		return ((struct surface *) NULL);
	}
	po = polyhedron_bunch (ms -> current_molecule);
	if (po == (struct surface *) NULL || po -> type != PHN_SURFACE) {
		set_error1 ("no polyhedron to contour");
		return ((struct surface *) NULL);
	}

	/* allocate contour bunch */
	obj = new_surface ();
	if (obj == NULL) {
		set_error1 ("(contour_molecule): mem alloc failure");
		return ((struct surface *) NULL);
	}
	/* initialize */
	obj -> type = CTR_SURFACE;
	obj -> n_polygon = 0L;
	obj -> n_phnvtx = 0L;
	obj -> n_phnedg = 0L;
	obj -> head_phnvtx = (struct phnvtx *) NULL;
	obj -> head_phnedg = (struct phnedg *) NULL;
	obj -> head_polygon = (struct polygon *) NULL;
	
	/* count edges that will have contour vertices */
	n_ctrvtx = 0L;
	for (edg = po -> head_phnedg; edg != NULL; edg = edg -> next) {
		/* clear old contour vertex */
		edg -> pvt[2] = (struct phnvtx *) NULL;
		edg -> upper = 0;
		vtx1 = edg -> pvt[0];
		if (vtx1 == NULL) {
			set_error1 ("contour_polyhedron: null polyhedron edge vertex pointer");
			return (NULL);
		}
		val1 = get_function (po, vtx1);
		vtx2 = edg -> pvt[1];
		if (vtx2 == NULL) {
			set_error1 ("contour_polyhedron: null polyhedron edge vertex pointer");
			return (NULL);
		}
		val2 = get_function (po, vtx2);
		if (val1 < level && level < val2) n_ctrvtx++;
		else if (val2 < level && level < val1) n_ctrvtx++;
	}
	sprintf (message,"%8ld edges needing contour vertices", n_ctrvtx);
	informd(message);
	
	if (n_ctrvtx <= 0) {
		/* empty bunch */
		sprintf (message, "%8.3f empty contour level (no vertices)", level);
		informd (message);
		return (obj);
	}
	/* allocate memory for contour vertices */
	ctr_vertices = (struct phnvtx **)
		allocate_pointers (PHNVTX, n_ctrvtx);
	if (ctr_vertices == NULL) {
		set_error1 ("(contour_molecule): memory full");
		return ((struct surface *) NULL);
	}
	
	/* create contour vertices */
	v = 0L;
	for (edg = po -> head_phnedg; edg != NULL; edg = edg -> next) {
		vtx1 = edg -> pvt[0];
		val1 = get_function (po, vtx1);
		vtx2 = edg -> pvt[1];
		val2 = get_function (po, vtx2);
		if (val1 < level && level < val2) {
			if (v >= n_ctrvtx) break;
			/* compute coordinates */
			range = val2 - val1;
			if (range <= 0.0) continue;
			fraction = (level - val1) / range;
			ctr_vtx = allocate_phnvtx ();
			if (ctr_vtx ==  NULL) {
				set_error1 ("(contour_polyhedron): vertex allocation failure");
				return ((struct surface *) NULL);
			}
			*(ctr_vertices + v) = ctr_vtx;
			ctr_vtx -> on = edg;
			ctr_vtx -> values[0] = level;
			for (k = 0; k < 3; k++)
				ctr_vtx -> center[k] = (1.0 - fraction) * vtx1 -> center[k] +
					fraction * vtx2 -> center[k];
			edg -> pvt[2] = ctr_vtx; edg -> upper = 1;
			v++;
			ctr_vtx -> number = v;
		}
		else if (val2 < level && level < val1) {
			if (v >= n_ctrvtx) break;
			/* compute coordinates */
			range = val1 - val2;
			if (range <= 0.0) continue;
			fraction = (level - val2) / range;
			ctr_vtx = allocate_phnvtx ();
			if (ctr_vtx ==  NULL) {
				set_error1 ("(contour_polyhedron): vertex allocation failure");
				return ((struct surface *) NULL);
			}
			*(ctr_vertices + v) = ctr_vtx;
			ctr_vtx -> on = edg;
			ctr_vtx -> values[0] = level;
			for (k = 0; k < 3; k++)
				ctr_vtx -> center[k] = (1.0 - fraction) * vtx2 -> center[k] +
					fraction * vtx1 -> center[k];
			edg -> pvt[2] = ctr_vtx; edg -> upper = 0;
			v++;
			ctr_vtx -> number = v;
		}
	}
	/* set up linked list */
	obj -> head_phnvtx = *ctr_vertices;
	for (v = 0; v < n_ctrvtx - 1; v++) {
		ctr_vtx = *(ctr_vertices + v);
		ctr_vtx -> next = *(ctr_vertices + v + 1);
	}
	
	/* count polygons that will have contour edges */
	n_ctredg = 0L;
	for (tri = po -> head_polygon; tri != NULL; tri = tri -> next) {
		n = 0;
		for (j = 0; j < 3; j++) {
			if (tri -> edg[j] -> pvt[2] != (struct phnvtx *) NULL) n++;
		}
		if (n == 0) continue;
		if (n != 2) {
			sprintf (message,
				"(contour_molecule) broken contour for parity = %d", n);
			set_error1(message);
			return ((struct surface *) NULL);
		}
		n_ctredg++;
	}
	if (n_ctredg <= 0) {
		/* empty bunch */
		sprintf (message, "%8.3f  empty contour level (no edges)", level);
		informd (message);
		return (obj);
	}
	
	/* allocate memory for contour edges */
	ctr_edges = (struct phnedg **)
		allocate_pointers (PHNEDG, n_ctredg);
	if (ctr_edges == NULL) {
		set_error1 ("(contour_molecule): memory full");
		return ((struct surface *) NULL);
	}
	
	/* create contour edges */
	e = 0L;
	for (tri = po -> head_polygon; tri != NULL; tri = tri -> next) {
		n = 0;
		for (j = 0; j < 3; j++) {
			if (tri -> edg[j] -> pvt[2] != (struct phnvtx *) NULL) {
				if (n == 0) {
					vtx1 = tri -> edg[j] -> pvt[2];
					orn1 = (tri -> orn[j] == tri -> edg[j] -> upper);
				}
				else if (n == 1) {
					vtx2 = tri -> edg[j] -> pvt[2];
					orn2 = (tri -> orn[j] == tri -> edg[j] -> upper);
				}
				n++;
			}
		}
		if (n == 0) continue;
		if (n != 2) {
			sprintf (message,
				"contour_polyhedron: broken contour for parity = %d", n);
			set_error1(message);
			return ((struct surface *) NULL);
		}
		/* order of vertices does matter */
		if (orn1 == orn2) {
			set_error1 ("contour_polyhedron: contour edge orientation ambiguity");
			sprintf (message, "orn1 = %d, orn2 = %d", orn1, orn2);
			set_error2 (message);
			return ((struct surface *) NULL);
		}
		if (orn1 == 1 && orn2 == 0) {
			vtx = vtx1;
			vtx1 = vtx2;
			vtx2 = vtx;
		}
		ctr_edg = allocate_phnedg ();
		if (ctr_edg ==  NULL) {
			set_error1 ("(contour_polyhedron): edge allocation failure");
			return ((struct surface *) NULL);
		}
		*(ctr_edges + e) = ctr_edg;
		ctr_edg -> on[0] = tri;
		ctr_edg -> pvt[0] = vtx1; ctr_edg -> vtxnum[0] = vtx1 -> number;
		ctr_edg -> pvt[1] = vtx2; ctr_edg -> vtxnum[1] = vtx2 -> number;
		e++;
		if (e >= n_ctredg) break;
	}

	/* set up linked list */
	obj -> head_phnedg = *ctr_edges;
	for (e = 0; e < n_ctredg - 1; e++) {
		ctr_edg = *(ctr_edges + e);
		ctr_edg -> next = *(ctr_edges + e + 1);
	}

	obj -> n_phnctr = 0;
	/* set up linked list of contours of linked lists of edges */
	for (;;) {
		/* look for unused edge to start contour */
		head_phnedg = NULL;
		for (e = 0; e < n_ctredg; e++) {
			ctr_edg = *(ctr_edges + e);
			if (ctr_edg -> used) continue;
			head_phnedg = ctr_edg;
			ctr_edg -> used = 1;
			break;
		}
		if (head_phnedg == NULL) break;
		ctr = allocate_phnctr ();
		if (ctr == NULL) {
			set_error1 ("contour_polyhedron: not enough memory");
			return (NULL);
		}
		if (obj -> head_phnctr == NULL)
			obj -> head_phnctr = ctr;
		else obj -> tail_phnctr -> next = ctr;
		obj -> tail_phnctr = ctr;
		ctr -> head_phnedg = head_phnedg;
		ctr -> tail_phnedg = head_phnedg;
		ctr -> n_phnedg = 1;
		obj -> n_phnctr++;
		for (;;) {
			/* look for unused edge to continue contour */
			next_edg = NULL;
			for (e = n_ctredg-1; e >= 0; e--) {
				ctr_edg = *(ctr_edges + e);
				if (ctr_edg -> used) continue;
				if (ctr -> tail_phnedg -> pvt[1] != ctr_edg -> pvt[0]) continue;
				ctr_edg -> used = 1;
				next_edg = ctr_edg;
				break;
			}
			if (next_edg == NULL) break;
			ctr -> tail_phnedg -> next_ctredg = next_edg;
			ctr -> tail_phnedg = next_edg;
			ctr -> n_phnedg++;
		}
	}
	
	/* transfer to structure variables */
	obj -> n_polygon = 0L;
	obj -> n_phnvtx = n_ctrvtx;
	obj -> n_phnedg = n_ctredg;
	obj -> head_polygon = (struct polygon *) NULL;
	obj -> phnvtx_handles = ctr_vertices;
	obj -> phnedg_handles = ctr_edges;

	return (obj);
}

struct phnctr *allocate_phnctr ()
{
	struct phnctr *ctr;

	ctr = (struct phnctr *) allocate_object (PHNCTR);
	if (ctr == NULL) {
		set_error1 ("(allocate_phnctr): mem alloc fails");
		return(NULL);
	}
	return (ctr);
}

void free_phnctr (struct phnctr *ctr)
{
	free_object (PHNCTR, (short *) ctr);
}


/*
	MSP
	Copyright 1993 by Michael L. Connolly
	All Rights Reserved
*/
