/* Molecular Surface Package Copyright 1995 Michael L. Connolly */
/* November 22, 2001 */

#include "ms.h"
#include "msenum.h"
#include "mslimit.h"
#include "msdefault.h"
#include "msstruct.h"
#include "msglobal.h"
#include "msfunction.h"


/* 
The purpose of hedron functions is to perform operations on a polyhedron:
(1) Creating a polyhedron by adding vertices, edges and triangles
(2) Making the vertex, edge and triangle elements contiguous
(3) Converting a surface-style polyhedron into a hedron style
(4) Converting a hedron into a surface-phnvtx-phnedg-phntri style
(5) Identifying vertices
(6) Removing degenerate elements (edge with two equal vertices)
(7) Removing duplicate elements (two edges pointing at same pair of vertices)
*/

/*
later: (a) worry about duplicate edges with reverse direction
       (b) check atom number match for identified vertices
       (c) use averaged coordinates for central vertex
*/

/* identify vertices closer than specified radius */
long purifyPolyhedron (struct surface *phn, double radius) {
	struct hedron *hed;
	char message[MAX_STRING];
	long removed = 0;
	long euler0 = 0;
	long euler1 = 0;
	sprintf (message, "%8.3f coalesce vertices distance cutoff", radius);
	inform (message);
	euler0 = phn -> n_phnvtx - phn -> n_phnedg + phn -> n_phntri;
	sprintf(message, "%8ld vertices", phn -> n_phnvtx);
	inform (message);
	sprintf(message, "%8ld edges", phn -> n_phnedg);
	inform (message);
	sprintf(message, "%8ld triangles", phn -> n_phntri);
	inform (message);
	hed = surfaceHedron (phn);
	removed = purifyHedron(hed, radius);
	sprintf(message, "%8ld elements removed from polyhedron", removed);
	informd (message);
	if (error()) return (0L);
	revisePolyhedron (phn, hed);
	euler1 = phn -> n_phnvtx - phn -> n_phnedg + phn -> n_phntri;
	if (euler0 != euler1) {
		sprintf(message, "%8ld Euler characteristic before", euler0);
		inform (message);
		sprintf(message, "%8ld Euler characteristic after", euler1);
		inform (message);
	}
	freeHedron(hed);
	return (removed);
}

/* identify vertices closer than specified radius */
long purifyHedron (struct hedron *hed, double radius) {
	char message[MAX_STRING];
	long total = 0;
	long removed = 0;
	long changed = 0;
	allVertexNeighbors (hed);
	if (error()) return (0L);
	allVertexGroups (hed, radius);
	if (error()) return (0L);
	sprintf(message, "%8ld vertex groups", hed -> nvtxgrp);
	if (hed -> nvtxgrp) inform(message);
	changed = collectChangedVertices (hed);
	if (error()) return (0L);
	sprintf(message, "%8ld changed vertices", changed);
	informd(message);
	changed = collectChangedEdges (hed);
	if (error()) return (0L);
	sprintf(message, "%8ld changed edges", changed);
	informd(message);
	changed = collectChangedTriangles (hed);
	if (error()) return (0L);
	sprintf(message, "%8ld changed triangles", changed);
	informd(message);
	updateEdges (hed);
	if (error()) return (0L);
	markDegenerateEdges (hed);
	if (error()) return (0L);
	hashEdges (hed);
	if (error()) return (0L);
	markDuplicateEdges (hed);
	if (error()) return (0L);
	updateTriangles (hed);
	if (error()) return (0L);
	markDegenerateTriangles (hed);
	if (error()) return (0L);
	hashTriangles (hed);
	if (error()) return (0L);
	markDuplicateTriangles (hed);
	if (error()) return (0L);
	removed = removeDegenerateTriangles (hed);
	if (error()) return (0L);
	total += removed;
	sprintf(message, "%8ld degenerate triangles removed", removed);
	if (removed > 0) inform(message);
	restoreDuplicateTriangles (hed);
	if (error()) return (0L);
	removed = removeDuplicateTriangles (hed);
	if (error()) return (0L);
	total += removed;
	sprintf(message, "%8ld duplicate triangles removed", removed);
	if (removed > 0) inform(message);
	removed = removeDegenerateEdges (hed);
	if (error()) return (0L);
	total += removed;
	sprintf(message, "%8ld degenerate edges removed", removed);
	if (removed > 0) inform(message);
	removed = removeDuplicateEdges (hed);
	if (error()) return (0L);
	total += removed;
	sprintf(message, "%8ld duplicate edges removed", removed);
	if (removed > 0) inform(message);
	removed = removeDuplicateVertices (hed);
	if (error()) return (0L);
	total += removed;
	sprintf(message, "%8ld duplicate vertices removed", removed);
	if (removed > 0) inform(message);
	markUnreferenced (hed);
	if (error()) return (0L);
	makeContiguous (hed);
	if (error()) return (0L);
	return (total);
}

/* convert a surface to a hedron */
struct hedron *surfaceHedron (struct surface *phn) {
	long v, e, t;
	struct hedron *hed = NULL;
	struct phnvtx *pvtx;
	struct phnedg *pedg;
	struct phntri *ptri;

	hed = newHedron (phn -> n_phnvtx, phn -> n_phnedg, phn -> n_phntri);
	if (hed == NULL) return (NULL);
	if (error()) return (NULL);
	hed -> van_der_Waals = phn -> van_der_Waals;
	for (v = 0; v < hed -> nvtx; v++) {
		pvtx = num2phnvtx (phn, v + 1);
		if (pvtx == NULL) return (NULL);
		addVertex (hed, v + 1, pvtx);
		if (error()) return (NULL);
	}
	for (e = 0; e < hed -> nedg; e++) {
		pedg = num2phnedg (phn, e + 1);
		if (pedg == NULL) return (NULL);
		addEdge (hed, e + 1, pedg);
		if (error()) return (NULL);
	}
	for (t = 0; t < hed -> ntri; t++) {
		ptri = num2phntri (phn, t + 1);
		if (ptri == NULL) return (NULL);
		addTriangle (hed, t + 1, ptri);
		if (error()) return (NULL);
	}
	return (hed);
}

/* copy polyhedron fields from second to first */
void revisePolyhedron (struct surface *phn, struct hedron *hed) {
	long v, e, t;
	struct hedvtx *hvtx;
	struct hededg *hedg;
	struct hedtri *htri;

	/* replace handle arrays */
	if (phn -> n_phnvtx == 0) return;
	if (phn -> n_phnedg == 0) return;
	if (phn -> n_phntri == 0) return;
	if (phn -> phnvtx_handles == NULL) return;
	if (phn -> phnedg_handles == NULL) return;
	if (phn -> phntri_handles == NULL) return;
	phn -> n_phnvtx = hed -> nvtx;
	phn -> n_phnedg = hed -> nedg;
	phn -> n_phntri = hed -> ntri;
	for (v = 0; v < hed -> nvtx; v++) {
		hvtx = *(hed -> vertices + v);
		*(phn -> phnvtx_handles + v) = hvtx -> vtx;
	}
	for (v = hed -> nvtx; v < hed -> maxvtx; v++) {
		*(phn -> phnvtx_handles + v) = NULL;
	}
	for (e = 0; e < hed -> nedg; e++) {
		hedg = *(hed -> edges + e);
		*(phn -> phnedg_handles + e) = hedg -> edg;
	}
	for (e = hed -> nedg; e < hed -> maxedg; e++) {
		*(phn -> phnedg_handles + e) = NULL;
	}
	for (t = 0; t < hed -> ntri; t++) {
		htri = *(hed -> triangles + t);
		*(phn -> phntri_handles + t) = htri -> tri;
	}
	for (t = hed -> ntri; t < hed -> maxtri; t++) {
		*(phn -> phntri_handles + t) = NULL;
	}

	/* redo links */
	link_polyhedron (phn);
}

void freeHedron(struct hedron *hed) {

	long ivtx;
	long iedg;
	long itri;
	struct hedvtx *vtx;
	struct hededg *edg;
	struct hedtri *tri;
	struct vtxgrp *grp = NULL;
	struct vtxgrp *nxtgrp = NULL;


	for (grp = hed -> head_vtxgrp; grp != NULL; grp = nxtgrp) {
		nxtgrp = grp -> next;
		free_vtxgrp (grp);
	}

	for (ivtx = 0; ivtx < hed -> nvtx; ivtx++) {
		vtx = *(hed -> vertices + ivtx);
		if (vtx != NULL) free_hedvtx (vtx);
	}
	free_pointers (HEDVTX, hed -> vertices);

	free_pointers (HEDVTX, hed -> changedVertices);

	for (iedg = 0; iedg < hed -> nedg; iedg++) {
		edg = *(hed -> edges + iedg);
		if (edg != NULL) free_hededg (edg);
	}
	free_pointers (HEDEDG, hed -> edges);

	free_pointers (HEDEDG, hed -> changedEdges);

	for (itri = 0; itri < hed -> ntri; itri++) {
		tri = *(hed -> triangles + itri);
		if (tri != NULL) free_hedtri (tri);
	}
	free_pointers (HEDTRI, hed -> triangles);

	free_pointers (HEDTRI, hed -> changedTriangles);

	free_longs(hed -> neighbors, 0, NEIGHBORS);

	free_pointers (HEDEDG, hed -> edghash);

	free_pointers (HEDTRI, hed -> trihash);

	free_object (HEDRON, (short *) hed);

	free_cache (VTXGRP);
	free_cache (HEDVTX);
	free_cache (HEDEDG);
	free_cache (HEDTRI);
	free_cache (HEDRON);
}

/* find all neighbors each vertex - that is, joined to it by an edge */
void allVertexNeighbors (struct hedron *hed) {
	long e, v, idx, vn0, vn1, vn, idx0, idx1;
	int j;
	struct hedvtx *hvtx, *hvtx0, *hvtx1;
	struct hededg *hedg;
	long *neighbors;
	char message[MAX_STRING];
    struct cept *ex;

	for (e = 0; e < hed -> nedg; e++) {
		hedg = *(hed -> edges + e);
		for (j = 0; j < 2; j++) {
			vn = hedg -> vtxnum[j];
			hvtx = *(hed -> vertices + vn - 1);
			hvtx -> degree++;
		}
	}
	idx = 0;
	for (v = 0; v < hed -> nvtx; v++) {
		hvtx = *(hed -> vertices + v);
		hvtx -> idx = idx;
		hvtx -> inc = 0;
		idx += hvtx -> degree;
	}
	neighbors = allocate_longs (2 * hed -> maxedg, 0, NEIGHBORS);
	if (neighbors == NULL) {
		ex = new_cept ( MEMORY_ERROR,  ALLOCATION,  FATAL_SEVERITY);
        add_variable (ex, NEIGHBORS, "neighbors");
        add_function (ex, "allVertexNeighbors");
        add_source (ex, "mshedron.c");
		return;
	}
	for (e = 0; e < hed -> nedg; e++) {
		hedg = *(hed -> edges + e);
		vn0 = hedg -> vtxnum[0];
		vn1 = hedg -> vtxnum[1];
		hvtx0 = *(hed -> vertices + vn0 - 1);
		hvtx1 = *(hed -> vertices + vn1 - 1);
		idx0 = hvtx0 -> idx + hvtx0 -> inc++;
		idx1 = hvtx1 -> idx + hvtx1 -> inc++;
		*(neighbors + idx0) = vn1;
		*(neighbors + idx1) = vn0;
		sprintf (message, "%6ld and %6ld are neighboring vertices",
			vn0, vn1);
		informd2 (message);
	}
	hed -> neighbors = neighbors;
}

/* find all groups of nearby vertices */
void allVertexGroups (struct hedron *hed, double radius) {
	int finished = 0;
	int g;
	struct hedvtx *gv;
	long vn = 0;
	struct vtxgrp *grp = NULL;
	char message[MAX_STRING];

	while (!finished) {
		vn = findUnusedVertex(hed);
		if (vn <= 0) {
			finished = 1;
			continue;
		}
		grp = findVertexGroup(hed, radius, vn);
		if (error()) return;
		if (grp == NULL) continue;
		sprintf(message, "vertex group with %ld members, central vertex %6ld",
			grp -> nvtx, grp -> central);
		informd2 (message);
		for (g = 0; g < grp -> nvtx; g++) {
			gv = grp -> vtxptrs[g];
			sprintf(message, "group member %6ld", gv -> number);
			informd2 (message);
		}
		if (hed -> head_vtxgrp == NULL)
			hed -> head_vtxgrp = grp;
		else hed -> tail_vtxgrp -> next = grp;
		hed -> tail_vtxgrp = grp;
		hed -> nvtxgrp++;
		identifyVertexGroup (hed, grp);
		if (error()) return;
	}
}

long findUnusedVertex (struct hedron *hed) {
	long v = 0;
	struct hedvtx *hv;
	for (v = 0; v < hed -> nvtx; v++) {
		hv = *(hed -> vertices + v);
		if (hv == NULL) continue;
		if (!hv -> used) {
			return (v+1);
		}
	}
	return (-1);
}

/* find the group, if it exists, of vertices too close to given */
struct vtxgrp *findVertexGroup (struct hedron *hed, double radius, long central) {
	int iv;
	long idx;
	long v = 0;
	long vn = 0;
	long nfound = 0;
	int enough = 0;
	struct hedvtx *hv, *cv;
	struct phnvtx *hpv, *cpv;
	struct vtxgrp *grp = NULL;
    struct cept *ex;

	cv = *(hed -> vertices + central - 1);
	cv -> used = 1;
	for (iv = 0; iv < cv -> degree; iv++) {
		idx = cv -> idx + iv;
		vn = *(hed -> neighbors + idx);
		v = vn - 1;
		hv = *(hed -> vertices + v);
		if (hv == NULL) continue;
		if (hv -> used) continue;
		cpv = cv -> vtx;
		hpv = hv -> vtx;
		enough = close_enough (cpv, hpv, radius, !hed -> van_der_Waals);
		if (enough) {
			if (nfound >= MAXVTXGRP) {
				ex = new_cept ( ARRAY_ERROR,  MSOVERFLOW,  FATAL_SEVERITY);
				add_function (ex, "findVertexGroup");
				add_source (ex, "mshedron.c");
                add_long (ex, "nfound", nfound);
                add_long (ex, "MAXVTXGRP", MAXVTXGRP);
				return (NULL);
			}
			if (grp == NULL) {
				grp = allocate_vtxgrp ();
			}
			grp -> vtxptrs[nfound] = hv;
			hv -> used = 1;
			hv -> changed = 1;
			nfound++;
		}
	}
	if (grp != NULL) {
		grp -> central = central;
		grp -> radius = radius;
		grp -> nvtx = nfound;
		cv -> iscentral = 1;
	}
	return (grp);
}

int close_enough (struct phnvtx *vtx0, struct phnvtx *vtx1, double radius, int same_atom) {
	double d2 = 0.0;
	double radius2 = 0.0;
	char message[MAX_STRING];

	/* must be same atom (except vdw) and component */
	if (same_atom && (vtx0 -> atm != vtx1 -> atm)) return (0);
	if (vtx0 -> comp != vtx1 -> comp) return (0);
	/* check whether two vertices are close enough */
	radius2 = radius * radius;
	d2 = distance_squared (vtx0 -> center, vtx1 -> center);
	sprintf (message, "%6.3f %6.3f", radius2, d2);
	informd2 (message);
	return (d2 <= radius2);
}

/* make all vertices point at central vertex as representative */
void identifyVertexGroup (struct hedron *hed, struct vtxgrp *grp) {
	int g = 0;
	long central = 0;
	struct hedvtx *hv;
	central = grp -> central;
	for (g = 0; g < grp -> nvtx; g++) {
		hv = grp -> vtxptrs[g];
		if (hv == NULL) continue;
		hv -> central = central;
	}
}

/* compute average position of group vertices, store in central vertex */
void averageVertexGroup (struct hedron *hed, struct vtxgrp *grp) {
	int g = 0;
	int k = 0;
	long central = 0;
	double sum[3];
	struct hedvtx *hv, *cv;
	struct phnvtx *phv, *pcv;

	central = grp -> central;
	cv = *(hed -> vertices + central - 1);
	pcv = cv -> vtx;
	for (k = 0; k < 3; k++)
		sum[k] = pcv -> center[k];
	for (g = 0; g < grp -> nvtx; g++) {
		hv = grp -> vtxptrs[g];
		if (hv == NULL) continue;
		phv = hv -> vtx;
		for (k = 0; k < 3; k++)
			sum[k] += phv -> center[k];
	}
	for (k = 0; k < 3; k++)
		pcv -> center[k] = sum[k] / (1 + grp -> nvtx);
}

/* central = original = representative */

/* if changed flag is true, store pointer in array */
long collectChangedVertices (struct hedron *hed) {
	long v = 0;
	struct hedvtx *hv;
	struct hedvtx **changedVertices;
    struct cept *ex;

	long nchgvtx = 0;
	for (v = 0; v < hed -> nvtx; v++) {
		hv = *(hed -> vertices + v);
		if (hv -> central != 0) {
			hv -> changed = 1;
			nchgvtx++;
		}
	}
	if (nchgvtx <= 0) {
		hed -> changedVertices = NULL;
		hed -> nchgvtx = 0;
		return (0L);
	}
	changedVertices = (struct hedvtx **)
		allocate_pointers (HEDVTX, nchgvtx);
	if (changedVertices == NULL) {
		ex = new_cept ( MEMORY_ERROR,  ALLOCATION,  FATAL_SEVERITY);
        add_object (ex, HEDVTX, "changedVertices");
        add_function (ex, "collectChangedVertices");
        add_source (ex, "mshedron.c");
		return (0L);
	}
	hed -> nchgvtx = nchgvtx;
	/* reset number of changed vertices to reuse as counter */
	nchgvtx = 0;
	for (v = 0; v < hed -> nvtx; v++) {
		hv = *(hed -> vertices + v);
		if (hv -> changed)
			*(changedVertices + nchgvtx++) = hv;
	}
	hed -> changedVertices = changedVertices;
	return (nchgvtx);
}

/* if changed flag is true, store pointer in array */
long collectChangedEdges (struct hedron *hed) {
	long e = 0;
	long vn = 0;
	int j = 0;
	struct hededg *he;
	struct hedvtx *hv;
	struct hededg **changedEdges;
    struct cept *ex;

	long nchgedg = 0;
	for (e = 0; e < hed -> nedg; e++) {
		he = *(hed -> edges + e);
		for (j = 0; j < 2; j++) {
			vn = he -> vtxnum[j];
			hv = *(hed -> vertices + vn - 1);
			if (hv -> changed) {
				he -> changed = 1;
				nchgedg++;
				break;
			}
		}
	}
	if (nchgedg <= 0) {
		hed -> changedEdges = NULL;
		hed -> nchgedg = 0;
		return (0L);
	}
	changedEdges = (struct hededg **)
		allocate_pointers (HEDEDG, nchgedg);
	if (changedEdges == NULL) {
		ex = new_cept ( MEMORY_ERROR,  ALLOCATION,  FATAL_SEVERITY);
        add_object (ex, HEDEDG, "changedEdges");
        add_function (ex, "collectChangedEdges");
        add_source (ex, "mshedron.c");
		return (0L);
	}
	hed -> nchgedg = nchgedg;
	/* reset number of changed edges to reuse as counter */
	nchgedg = 0;
	for (e = 0; e < hed -> nedg; e++) {
		he = *(hed -> edges + e);
		if (he -> changed)
			*(changedEdges + nchgedg++) = he;
	}
	hed -> changedEdges = changedEdges;
	return (nchgedg);
}

/* if changed flag is true, store pointer in array */
long collectChangedTriangles (struct hedron *hed) {
	long t = 0;
	long vn = 0;
	int j = 0;
	struct hedtri *ht;
	struct hedvtx *hv;
	struct hedtri **changedTriangles;
    struct cept *ex;

	long nchgtri = 0;
	for (t = 0; t < hed -> ntri; t++) {
		ht = *(hed -> triangles + t);
		for (j = 0; j < 3; j++) {
			vn = ht -> vtxnum[j];
			hv = *(hed -> vertices + vn - 1);
			if (hv -> changed) {
				ht -> changed = 1;
				nchgtri++;
				break;
			}
		}
	}
	if (nchgtri <= 0) {
		hed -> changedTriangles = NULL;
		hed -> nchgtri = 0;
		return (0L);
	}
	changedTriangles = (struct hedtri **)
		allocate_pointers (HEDTRI, nchgtri);
	if (changedTriangles == NULL) {
		ex = new_cept ( MEMORY_ERROR,  ALLOCATION,  FATAL_SEVERITY);
        add_object (ex, HEDTRI, "changedTriangles");
        add_function (ex, "collectChangedTriangles");
        add_source (ex, "mshedron.c");
		return (0L);
	}
	hed -> nchgtri = nchgtri;
	/* reset number of changed triangles to reuse as counter */
	nchgtri = 0;
	for (t = 0; t < hed -> ntri; t++) {
		ht = *(hed -> triangles + t);
		if (ht -> changed)
			*(changedTriangles + nchgtri++) = ht;
	}
	hed -> changedTriangles = changedTriangles;
	return (nchgtri);
}

/* change vertex pointers and indices to point at original vertex */
void updateEdges (struct hedron *hed) {
	long e = 0;
	long vn = 0;
	int j = 0;
	struct hededg *he;
	struct hedvtx *hv, *cv;
	for (e = 0; e < hed -> nchgedg; e++) {
		he = *(hed -> changedEdges + e);
		for (j = 0; j < 2; j++) {
			vn = he -> vtxnum[j];
			hv = *(hed -> vertices + vn - 1);
			if (hv -> central != 0) {
				cv = *(hed -> vertices + hv -> central - 1);
				he -> vtxnum[j] = hv -> central;
				he -> edg -> vtxnum[j] = hv -> central;
				/* is this a good idea? */
				he -> edg -> pvt[j] = cv -> vtx;
			}
		}
	}
}

/* if the endpoints are the same, mark as degenerate for later removal */
void markDegenerateEdges (struct hedron *hed) {
	long e;
	struct hededg *he;
	for (e = 0; e < hed -> nchgedg; e++) {
		he = *(hed -> changedEdges + e);
		if (he -> vtxnum[0] == he -> vtxnum[1])
			he -> degenerate = 1;
	}
}

void hashEdges (struct hedron *hed) {
	long e, idx;
	struct hededg *he;
	struct hededg **edghash;
    struct cept *ex;

	edghash = (struct hededg **)
		allocate_pointers (HEDEDG, 2 * hed -> maxvtx);
	if (edghash == NULL) {
		ex = new_cept (MEMORY_ERROR,  ALLOCATION,  FATAL_SEVERITY);
        add_variable (ex, HEDEDG, "edghash");
        add_function (ex, "hashEdges");
        add_source (ex, "mshedron.c");
		return;
	}
	hed -> edghash = edghash;
	for (e = 0; e < hed -> nedg; e++) {
		he = *(hed -> edges + e);
		if (he -> degenerate) continue;
		idx = he -> vtxnum[0] + he -> vtxnum[1] - 2;
		he -> next = *(edghash + idx);
		*(edghash + idx) = he;
	}
}

/* if two edges have the same endpoints, identify them */
void markDuplicateEdges (struct hedron *hed) {
	long idx = 0;
	long vn0 = 0;
	long vn1 = 0;
	char message[MAX_STRING];
	struct hededg *he0, *he1;
	struct hedvtx *hv0, *hv1;
 
	/* compare non-degenerate edges */
	for (idx = 0; idx < 2 * hed -> maxvtx; idx++) {
		for (he0 = *(hed -> edghash + idx); he0 != NULL; he0 = he0 -> next) {
			if (he0 -> duplicate) continue;
			vn0 = he0 -> vtxnum[0];
			vn1 = he0 -> vtxnum[1];
			hv0 = *(hed -> vertices + vn0 - 1);
			hv1 = *(hed -> vertices + vn1 - 1);
			if (!hv0 -> iscentral && !hv1 -> iscentral) continue;
			for (he1 = *(hed -> edghash + idx); he1 != NULL; he1 = he1 -> next) {
				if (he0 == he1) continue;
				if (he1 -> duplicate) continue;
				if (he1 -> isrepresentative) continue;
				if (he0 -> vtxnum[0] == he1 -> vtxnum[0] &&
					he0 -> vtxnum[1] == he1 -> vtxnum[1]) {
					he1 -> duplicate = 1;
					he1 -> representative = he0 -> number;
					he0 -> isrepresentative = 1;
					sprintf(message, "%6ld %6ld is duplicate edge",
						he0 -> vtxnum[0], he0 -> vtxnum[1]);
					informd2(message);
				}
				/* check for edge with opposite direction */
				else if (he0 -> vtxnum[0] == he1 -> vtxnum[1] &&
					he0 -> vtxnum[1] == he1 -> vtxnum[0]) {
					he1 -> duplicate = 1;
					he1 -> representative = - he0 -> number;
					he0 -> isrepresentative = 1;
					sprintf(message, "%6ld %6ld is duplicate edge",
						he0 -> vtxnum[0], he0 -> vtxnum[1]);
					informd2(message);
				}
			}
		}
	}
}

/* change vertex pointers and indices to point at central vertex */
/* change edge pointers and indices to point at representative edge */
void updateTriangles (struct hedron *hed) {
	char message[MAX_STRING];
	long vn = 0;
	long en = 0;
	long ren = 0;
	long aben = 0;
	long abren = 0;
	int j = 0;
	long t = 0;
	struct hedtri *ht;
	struct hededg *he, *her;
	struct hedvtx *hv;
	for (t = 0; t < hed -> ntri; t++) {
		ht = *(hed -> triangles + t);
		for (j = 0; j < 3; j++) {
			vn = ht -> vtxnum[j];
			hv = *(hed -> vertices + vn - 1);
			if (hv -> central != 0) {
				ht -> vtxnum[j] = hv -> central;
				ht -> tri -> vtxnum[j] = hv -> central;
			}
		}
		for (j = 0; j < 3; j++) {
			en = ht -> edgnum[j];
			aben = abs(en);
			he = *(hed -> edges + aben - 1);
			ren = he -> representative;
			abren = abs(ren);
			if (ren != 0) {
				ht -> edgnum[j] = ren;
				if (en < 0) ht -> edgnum[j] *= -1;
				her = *(hed -> edges + abren - 1);
				ht -> tri -> edg[j] = her -> edg;
			}
		}
	}
}

/* if two vertices the same, two edges the same, or degenerate edge,
   then the triangle is degenerate; mark for removal */
void markDegenerateTriangles (struct hedron *hed) {
	char message[MAX_STRING];
	int j;
	long en, aben, t;
	struct hedtri *ht;
	struct hededg *he;

	for (t = 0; t < hed -> nchgtri; t++) {
		ht = *(hed -> changedTriangles + t);
		if (ht -> vtxnum[0] == ht -> vtxnum[1])
			ht -> degenerate = 1;
		else if (ht -> vtxnum[1] == ht -> vtxnum[2])
			ht -> degenerate = 1;
		else if (ht -> vtxnum[2] == ht -> vtxnum[0])
			ht -> degenerate = 1;
		else if (ht -> edgnum[0] == ht -> edgnum[1])
			ht -> degenerate = 1;
		else if (ht -> edgnum[1] == ht -> edgnum[2])
			ht -> degenerate = 1;
		else if (ht -> edgnum[2] == ht -> edgnum[0])
			ht -> degenerate = 1;
		for (j = 0; j < 3; j++) {
			en = ht -> edgnum[j];
			aben = abs(en);
			he = *(hed -> edges + aben - 1);
			if (he -> degenerate) ht -> degenerate = 1;
		}
	}
}

void hashTriangles (struct hedron *hed) {
	long t, idx;
	struct hedtri *ht;
	struct hedtri **trihash;
    struct cept *ex;

	trihash = (struct hedtri **)
		allocate_pointers (HEDTRI, 3 * hed -> maxvtx);
	if (trihash == NULL) {
		ex = new_cept (MEMORY_ERROR,  ALLOCATION,  FATAL_SEVERITY);
        add_variable (ex, HEDTRI, "trihash");
        add_function (ex, "hashTriangles");
        add_source (ex, "mshedron.c");
		return;
	}
	hed -> trihash = trihash;
	for (t = 0; t < hed -> ntri; t++) {
		ht = *(hed -> triangles + t);
		if (ht -> degenerate) continue;
		idx = ht -> vtxnum[0] + ht -> vtxnum[1] + ht -> vtxnum[2] - 3;
		ht -> next = *(trihash + idx);
		*(trihash + idx) = ht;
	}
}


/* if two triangles have the same endpoints, error condition */
void markDuplicateTriangles (struct hedron *hed) {
	char message[MAX_STRING];
	long idx;
	struct hedtri *ht0, *ht1;
	long vn0 = 0;
	long vn1 = 0;
	long vn2 = 0;
	long vn3 = 0;
	long vn4 = 0;
	long vn5 = 0;
	long en0 = 0;
	long en1 = 0;
	long en2 = 0;
	long en3 = 0;
	long en4 = 0;
	long en5 = 0;
	struct hedvtx *hv0, *hv1, *hv2;
	struct hedvtx *hv3, *hv4, *hv5;

	/* compare non-degenerate triangles */
	for (idx = 0; idx < 3 * hed -> maxvtx; idx++) {
		for (ht0 = *(hed -> trihash + idx); ht0 != NULL; ht0 = ht0 -> next) {
			if (ht0 -> duplicate) continue;
			vn0 = ht0 -> vtxnum[0];
			vn1 = ht0 -> vtxnum[1];
			vn2 = ht0 -> vtxnum[2];
			hv0 = *(hed -> vertices + vn0 - 1);
			hv1 = *(hed -> vertices + vn1 - 1);
			hv2 = *(hed -> vertices + vn2 - 1);
			en0 = ht0 -> edgnum[0];
			en1 = ht0 -> edgnum[1];
			en2 = ht0 -> edgnum[2];
			if (!hv0 -> iscentral && !hv1 -> iscentral && !hv2 -> iscentral) continue;
			for (ht1 = *(hed -> trihash + idx); ht1 != NULL; ht1 = ht1 -> next) {
				if (ht1 -> duplicate) continue;
				if (ht1 -> isrepresentative) continue;
				if (ht0 == ht1) continue;
				vn3 = ht1 -> vtxnum[0];
				vn4 = ht1 -> vtxnum[1];
				vn5 = ht1 -> vtxnum[2];
				hv3 = *(hed -> vertices + vn3 - 1);
				hv4 = *(hed -> vertices + vn4 - 1);
				hv5 = *(hed -> vertices + vn5 - 1);
				en3 = ht1 -> edgnum[0];
				en4 = ht1 -> edgnum[1];
				en5 = ht1 -> edgnum[2];
				if (vn0 == vn3 && vn1 == vn4 && vn2 == vn5) {
					ht1 -> duplicate = 1;
					ht1 -> representative = ht0 -> number;
					ht0 -> isrepresentative = 1;
					sprintf (message, "%8ld %8ld %8ld duplicate triangle (standard)", vn0, vn1, vn2);
					informd2 (message);
					sprintf (message, "%8ld %8ld %8ld edges", en0, en1, en2);
					informd2 (message);
					sprintf (message, "%8ld %8ld %8ld edges", en3, en4, en5);
					informd2 (message);
				}
				else if (vn0 == vn4 && vn1 == vn5 && vn2 == vn3) {
					ht1 -> duplicate = 1;
					ht1 -> representative = ht0 -> number;
					ht0 -> isrepresentative = 1;
					sprintf (message, "%8ld %8ld %8ld duplicate triangle (turned)", vn0, vn1, vn2);
					informd2 (message);
					sprintf (message, "%8ld %8ld %8ld edges", en0, en1, en2);
					informd2 (message);
					sprintf (message, "%8ld %8ld %8ld edges", en3, en4, en5);
					informd2 (message);
				}
				else if (vn0 == vn5 && vn1 == vn3 && vn2 == vn4) {
					ht1 -> duplicate = 1;
					ht1 -> representative = ht0 -> number;
					ht0 -> isrepresentative = 1;
					sprintf (message, "%8ld %8ld %8ld duplicate triangle (turned)", vn0, vn1, vn2);
					informd2 (message);
					sprintf (message, "%8ld %8ld %8ld edges", en0, en1, en2);
					informd2 (message);
					sprintf (message, "%8ld %8ld %8ld edges", en3, en4, en5);
					informd2 (message);
				}
				/* check for other orientation */
				if (vn0 == vn5 && vn1 == vn4 && vn2 == vn3) {
					ht0 -> duplicate = 1;
					ht1 -> duplicate = 1;
					ht1 -> representative = ht0 -> number;
					ht0 -> isrepresentative = 1;
					ht1 -> isrepresentative = 1;
					sprintf (message, "%8ld %8ld %8ld duplicate triangle (reversed)", vn0, vn1, vn2);
					informd2 (message);
					sprintf (message, "%8ld %8ld %8ld edges", en0, en1, en2);
					informd2 (message);
					sprintf (message, "%8ld %8ld %8ld edges", en3, en4, en5);
					informd2 (message);
				}
				else if (vn0 == vn4 && vn1 == vn3 && vn2 == vn5) {
					ht0 -> duplicate = 1;
					ht1 -> duplicate = 1;
					ht1 -> representative = ht0 -> number;
					ht0 -> isrepresentative = 1;
					ht1 -> isrepresentative = 1;
					sprintf (message, "%8ld %8ld %8ld duplicate triangle (reversed)", vn0, vn1, vn2);
					informd2 (message);
					sprintf (message, "%8ld %8ld %8ld edges", en0, en1, en2);
					informd2 (message);
					sprintf (message, "%8ld %8ld %8ld edges", en3, en4, en5);
					informd2 (message);
				}
				else if (vn0 == vn3 && vn1 == vn5 && vn2 == vn4) {
					ht0 -> duplicate = 1;
					ht1 -> duplicate = 1;
					ht1 -> representative = ht0 -> number;
					ht0 -> isrepresentative = 1;
					ht1 -> isrepresentative = 1;
					sprintf (message, "%8ld %8ld %8ld duplicate triangle (reversed)", vn0, vn1, vn2);
					informd2 (message);
					sprintf (message, "%8ld %8ld %8ld edges", en0, en1, en2);
					informd2 (message);
					sprintf (message, "%8ld %8ld %8ld edges", en3, en4, en5);
					informd2 (message);
				}
			}
		}
	}
}

/* if triangle edges will become unreferenced, don't remove triangle */
void restoreDuplicateTriangles (struct hedron *hed) {
	int keep = 0;
	int j = 0;
	long nKept = 0;
	long t = 0;
	long en0, aben0;
	long nrefs = 0;
	struct hededg *he0;
	struct hedtri *ht;
	char message[MAX_STRING];
	for (t = 0; t < hed -> ntri; t++) {
		ht = *(hed -> triangles + t);
		if (ht == NULL) continue;
		if (ht -> degenerate) continue;
		if (!ht -> duplicate) continue;
		keep = 0;
		/* check each edge */
		for (j = 0; j < 3; j++) {
			en0 = ht -> edgnum[j];
			aben0 = abs (en0);
			he0 = *(hed -> edges + aben0 - 1);
			nrefs = countEdgeReferences (hed, aben0);
			if (error ()) return;
			if (nrefs <= 2) {
				/* could be 1 if already removed other duplicate */
				keep++;
			}
		}
		if (keep >= 3) {
			/* keep triangle by marking not duplicate */
			ht -> duplicate = 0;
			nKept++;
		}
	}
	sprintf(message, "%8ld duplicate triangles kept", nKept);
	if (nKept > 0) inform (message);
}

long countEdgeReferences (struct hedron *hed, long aben) {
	long t = 0;
	int j = 0;
	long en0, aben0;
	long nrefs = 0;
	struct hedtri *ht;

	for (t = 0; t < hed -> ntri; t++) {
		ht = *(hed -> triangles + t);
		if (ht == NULL) continue;
		if (ht -> degenerate) continue;
		/* check each edge */
		for (j = 0; j < 3; j++) {
			en0 = ht -> edgnum[j];
			aben0 = abs (en0);
			if (aben0 == aben) nrefs++;
		}
	}
	return (nrefs);
}

long removeDegenerateTriangles (struct hedron *hed) {
	long nDegenerates = 0;
	long t = 0;
	struct hedtri *ht;
	for (t = 0; t < hed -> ntri; t++) {
		ht = *(hed -> triangles + t);
		if (ht == NULL) continue;
		if (ht -> degenerate) {
			removeTriangle (hed, t + 1);
			nDegenerates++;
		}
	}
	return (nDegenerates);
}

/* if two triangles have the same vertices and edges,
   then they are duplicates; remove one of them */
long removeDuplicateTriangles (struct hedron *hed) {
	long nDuplicates = 0;
	long t = 0;
	struct hedtri *ht;
	for (t = 0; t < hed -> ntri; t++) {
		ht = *(hed -> triangles + t);
		if (ht == NULL) continue;
		if (ht -> duplicate) {
			removeTriangle (hed, t + 1);
			nDuplicates++;
		}
	}
	return (nDuplicates);
}

long removeDegenerateEdges (struct hedron *hed) {
	long nDegenerates = 0;
	long e = 0;
	struct hededg *he;
	for (e = 0; e < hed -> nedg; e++) {
		he = *(hed -> edges + e);
		if (he == NULL) continue;
		if (he -> degenerate) {
			removeEdge (hed, e + 1);
			nDegenerates++;
		}
	}
	return (nDegenerates);
}

long removeDuplicateEdges (struct hedron *hed) {
	long nDuplicates = 0;
	long e = 0;
	struct hededg *he;
	for (e = 0; e < hed -> nedg; e++) {
		he = *(hed -> edges + e);
		if (he == NULL) continue;
		if (he -> duplicate) {
			removeEdge (hed, e + 1);
			nDuplicates++;
		}
	}
	return (nDuplicates);
}

long removeDuplicateVertices (struct hedron *hed) {
	int j;
	long nDuplicates = 0;
	struct hedvtx *hv = NULL;
	struct vtxgrp *grp = NULL;
	for (grp = hed -> head_vtxgrp; grp != NULL; grp = grp -> next) {
		for (j = 0; j < grp -> nvtx; j++) {
			hv = grp -> vtxptrs[j];
			removeVertex (hed, hv -> number);
			nDuplicates++;
		}
	}
	return (nDuplicates);
}

void markUnreferenced (struct hedron *hed) {
	long e = 0;
	long v = 0;
	long t = 0;
	long vn0, vn1, vn2;
	long en0, en1, en2;
	long aben0, aben1, aben2;
	long unv = 0;
	long une = 0;
	struct hedvtx *hv0, *hv1, *hv2;
	struct hedvtx *hv;
	struct hededg *he;
	struct hededg *he0, *he1, *he2;
	struct hedtri *ht;
	char message[MAX_STRING];

	/* initialization */
	for (v = 0; v < hed -> nvtx; v++) {
		hv = *(hed -> vertices + v);
		if (hv == NULL) continue;
		hv -> unreferenced = 1;
	}
	for (e = 0; e < hed -> nedg; e++) {
		he = *(hed -> edges + e);
		if (he == NULL) continue;
		he -> unreferenced = 1;
	}
	for (e = 0; e < hed -> nedg; e++) {
		he = *(hed -> edges + e);
		if (he == NULL) continue;
		vn0 = he -> vtxnum[0];
		vn1 = he -> vtxnum[1];
		hv0 = *(hed -> vertices + vn0 - 1);
		hv1 = *(hed -> vertices + vn1 - 1);
		/* you are not alone */
		hv0 -> unreferenced = 0;
		hv1 -> unreferenced = 0;
	}
	for (t = 0; t < hed -> ntri; t++) {
		ht = *(hed -> triangles + t);
		if (ht == NULL) continue;
		vn0 = ht -> vtxnum[0];
		vn1 = ht -> vtxnum[1];
		vn2 = ht -> vtxnum[2];
		hv0 = *(hed -> vertices + vn0 - 1);
		hv1 = *(hed -> vertices + vn1 - 1);
		hv2 = *(hed -> vertices + vn2 - 1);
		/* you are not alone */
		hv0 -> unreferenced = 0;
		hv1 -> unreferenced = 0;
		hv2 -> unreferenced = 0;
		en0 = ht -> edgnum[0];
		en1 = ht -> edgnum[1];
		en2 = ht -> edgnum[2];
		aben0 = abs (en0);
		he0 = *(hed -> edges + aben0 - 1);
		aben1 = abs (en1);
		he1 = *(hed -> edges + aben1 - 1);
		aben2 = abs (en2);
		he2 = *(hed -> edges + aben2 - 1);
		/* you are not alone */
		he0 -> unreferenced = 0;
		he1 -> unreferenced = 0;
		he2 -> unreferenced = 0;
	}
	/* count */
	unv = 0;
	for (v = 0; v < hed -> nvtx; v++) {
		hv = *(hed -> vertices + v);
		if (hv == NULL) continue;
		if (hv -> unreferenced) {
			unv++;
			sprintf(message, "%8ld vertex unreferenced", v + 1);
			informd(message);
		}
	}
	une = 0;
	for (e = 0; e < hed -> nedg; e++) {
		he = *(hed -> edges + e);
		if (he == NULL) continue;
		if (he -> unreferenced) {
			une++;
			vn0 = he -> vtxnum[0];
			vn1 = he -> vtxnum[1];
			sprintf(message, "%8ld edge %8ld %8ld unreferenced", e + 1, vn0, vn1);
			informd(message);
		}
	}
	sprintf(message, "%8ld unreferenced vertices", unv);
	if (unv > 0) inform(message);
	sprintf(message, "%8ld unreferenced edges", une);
	if (une > 0) inform(message);
}


/* make vertices, edges and triangles contiguously numbered */

void makeContiguous (struct hedron *hed) {
	newNumberVertices (hed);
	if (error()) return;
	newNumberEdges (hed);
	if (error()) return;
	newNumberTriangles (hed);
	if (error()) return;
	updateEdgeToVertex (hed);
	if (error()) return;
	updateTriangleToVertex (hed);
	if (error()) return;
	updateTriangleToEdge (hed);
	if (error()) return;
	compactVertices (hed);
	if (error()) return;
	compactEdges (hed);
	if (error()) return;
	compactTriangles (hed);
	if (error()) return;
}

void newNumberVertices (struct hedron *hed) {
	long v = 0;
	long newber = 0;
	struct hedvtx *hv;
	for (v = 0; v < hed -> nvtx; v++) {
		hv = *(hed -> vertices + v);
		if (hv != NULL)
			hv -> newber = ++newber;
	}
}

void newNumberEdges (struct hedron *hed) {
	long e = 0;
	long newber = 0;
	struct hededg *he;
	for (e = 0; e < hed -> nedg; e++) {
		he = *(hed -> edges + e);
		if (he != NULL)
			he -> newber = ++newber;
	}
}

void newNumberTriangles (struct hedron *hed) {
	long t = 0;
	long newber = 0;
	struct hedtri *ht;
	for (t = 0; t < hed -> ntri; t++) {
		ht = *(hed -> triangles + t);
		if (ht != NULL)
			ht -> newber = ++newber;
	}
}

void updateEdgeToVertex (struct hedron *hed) {
	char message[MAX_STRING];
	long e = 0;
	long vn = 0;
	int j = 0;
	struct hededg *he;
	struct hedvtx *hv;
	for (e = 0; e < hed -> nedg; e++) {
		he = *(hed -> edges + e);
		if (he == NULL) continue;
		for (j = 0; j < 2; j++) {
			vn = he -> vtxnum[j];
			hv = *(hed -> vertices + vn - 1);
			he -> vtxnum[j] = hv -> newber;
			he -> edg -> vtxnum[j] = hv -> newber;
		}
	}
}

void updateTriangleToVertex (struct hedron *hed) {
	char message[MAX_STRING];
	long t = 0;
	long vn = 0;
	int j = 0;
	struct hedtri *ht;
	struct hedvtx *hv;
	for (t = 0; t < hed -> ntri; t++) {
		ht = *(hed -> triangles + t);
		if (ht == NULL) continue;
		for (j = 0; j < 3; j++) {
			vn = ht -> vtxnum[j];
			hv = *(hed -> vertices + vn - 1);
			ht -> vtxnum[j] = hv -> newber;
			ht -> tri -> vtxnum[j] = hv -> newber;
		}
	}
}

void updateTriangleToEdge (struct hedron *hed) {
	char message[MAX_STRING];
	long t = 0;
	long en = 0;
	long aben = 0;
	int j = 0;
	struct hedtri *ht;
	struct hededg *he;
	for (t = 0; t < hed -> ntri; t++) {
		ht = *(hed -> triangles + t);
		if (ht == NULL) continue;
		for (j = 0; j < 3; j++) {
			en = ht -> edgnum[j];
			aben = abs (en);
			he = *(hed -> edges + aben - 1);
			ht -> edgnum[j] = he -> newber;
			if (en < 0) ht -> edgnum[j] *= -1;
			ht -> tri -> edgnum[j] = he -> newber;
			if (en < 0) ht -> tri -> edgnum[j] *= -1;
		}
	}
}

void compactVertices (struct hedron *hed) {
	long v0 = 0;
	long v1 = 0;
	struct hedvtx *hv0 = NULL;
	struct hedvtx *hv1 = NULL;
    struct cept *ex;

	for (v1 = 0; v1 < hed -> nvtx; v1++) {
		hv0 = *(hed -> vertices + v0);
		hv1 = *(hed -> vertices + v1);
		if (hv0 == NULL && hv1 != NULL) {
			*(hed -> vertices + v0) = hv1;
			*(hed -> vertices + v1) = NULL;
		}
		if (*(hed -> vertices + v0) != NULL) v0++;
	}
	hed -> nvtx = v0;
	/* check */
	for (v0 = 0; v0 < hed -> nvtx; v0++) {
		hv0 = *(hed -> vertices + v0);
		if (hv0 == NULL) {
			ex = new_cept ( LOGIC_ERROR,  INCONSISTENCY,  FATAL_SEVERITY);
			add_function (ex, "compactVertices");
			add_source (ex, "mshedron.c");
			add_long (ex, "vertex number", v0 + 1);
			return;
		}
	}
}

void compactEdges (struct hedron *hed) {
	long e0 = 0;
	long e1 = 0;
	struct hededg *he0 = NULL;
	struct hededg *he1 = NULL;
    struct cept *ex;

	for (e1 = 0; e1 < hed -> nedg; e1++) {
		he0 = *(hed -> edges + e0);
		he1 = *(hed -> edges + e1);
		if (he0 == NULL && he1 != NULL) {
			*(hed -> edges + e0) = he1;
			*(hed -> edges + e1) = NULL;
		}
		if (*(hed -> edges + e0) != NULL) e0++;
	}
	hed -> nedg = e0;
	/* check */
	for (e0 = 0; e0 < hed -> nedg; e0++) {
		he0 = *(hed -> edges + e0);
		if (he0 == NULL) {
			ex = new_cept ( LOGIC_ERROR,  INCONSISTENCY,  FATAL_SEVERITY);
			add_function (ex, "compactEdges");
			add_source (ex, "mshedron.c");
			add_long (ex, "edge number", e0 + 1);
			return;
		}
	}
}

void compactTriangles (struct hedron *hed) {
	long t0 = 0;
	long t1 = 0;
	struct hedtri *ht0 = NULL;
	struct hedtri *ht1 = NULL;
    struct cept *ex;

	for (t1 = 0; t1 < hed -> ntri; t1++) {
		ht0 = *(hed -> triangles + t0);
		ht1 = *(hed -> triangles + t1);
		if (ht0 == NULL && ht1 != NULL) {
			*(hed -> triangles + t0) = ht1;
			*(hed -> triangles + t1) = NULL;
		}
		if (*(hed -> triangles + t0) != NULL) t0++;
	}
	hed -> ntri = t0;
	/* check */
	for (t0 = 0; t0 < hed -> ntri; t0++) {
		ht0 = *(hed -> triangles + t0);
		if (ht0 == NULL) {
			ex = new_cept ( LOGIC_ERROR,  INCONSISTENCY,  FATAL_SEVERITY);
			add_function (ex, "compactTriangles");
			add_source (ex, "mshedron.c");
			add_long (ex, "triangle number", t0 + 1);
			return;
		}
	}
}

struct hedvtx *getVertex (struct hedron *hed, long number) {
	struct hedvtx *vtx = NULL;
	long idx = 0;
	idx = number - 1;
	if (idx < 0) return (NULL);
	if (idx >= hed -> maxvtx) return (NULL);
	vtx = *(hed -> vertices + idx);
	return (vtx);
}

struct hededg *getEdge (struct hedron *hed, long number) {
	struct hededg *edg = NULL;
	long idx = 0;
	idx = number - 1;
	if (idx < 0) return (NULL);
	if (idx >= hed -> maxedg) return (NULL);
	edg = *(hed -> edges + idx);
	return (edg);
}

struct hedtri *getTriangle (struct hedron *hed, long number) {
	struct hedtri *tri = NULL;
	long idx = 0;
	idx = number - 1;
	if (idx < 0) return (NULL);
	if (idx >= hed -> maxtri) return (NULL);
	tri = *(hed -> triangles + idx);
	return (tri);
}

struct hedvtx *allocate_hedvtx ()
{
	struct hedvtx *vtx;
    struct cept *ex;

	vtx = (struct hedvtx *) allocate_object (HEDVTX);
	if (vtx == NULL) {
		ex = new_cept ( MEMORY_ERROR,  ALLOCATION,  FATAL_SEVERITY);
        add_object (ex, HEDVTX, "vtx");
        add_function (ex, "allocate_hedvtx");
        add_source (ex, "mshedron.c");
		return(NULL);
	}
	return (vtx);
}

void free_hedvtx (struct hedvtx *vtx)
{
	free_object (HEDVTX, (short *) vtx);
}

struct hededg *allocate_hededg ()
{
	struct hededg *edg;
    struct cept *ex;

	edg = (struct hededg *) allocate_object (HEDEDG);
	if (edg == NULL) {
		ex = new_cept ( MEMORY_ERROR,  ALLOCATION,  FATAL_SEVERITY);
        add_object (ex, HEDEDG, "edg");
        add_function (ex, "allocate_hededg");
        add_source (ex, "mshedron.c");
		return(NULL);
	}
	return (edg);
}

void free_hededg (struct hededg *edg)
{
	free_object (HEDEDG, (short *) edg);
}

struct hedtri *allocate_hedtri ()
{
	struct hedtri *tri;
    struct cept *ex;

	tri = (struct hedtri *) allocate_object (HEDTRI);
	if (tri == NULL) {
		ex = new_cept ( MEMORY_ERROR,  ALLOCATION,  FATAL_SEVERITY);
        add_object (ex, HEDTRI, "tri");
        add_function (ex, "allocate_hedtri");
        add_source (ex, "mshedron.c");
		return(NULL);
	}
	return (tri);
}

void free_hedtri (struct hedtri *tri)
{
	free_object (HEDTRI, (short *) tri);
}

struct vtxgrp *allocate_vtxgrp ()
{
	struct vtxgrp *grp;
    struct cept *ex;

	grp = (struct vtxgrp *) allocate_object (VTXGRP);
	if (grp == NULL) {
		ex = new_cept ( MEMORY_ERROR,  ALLOCATION,  FATAL_SEVERITY);
        add_object (ex, VTXGRP, "grp");
        add_function (ex, "allocate_vtxgrp");
        add_source (ex, "mshedron.c");
		return(NULL);
	}
	return (grp);
}

void free_vtxgrp (struct vtxgrp *grp)
{
	free_object (VTXGRP, (short *) grp);
}

void removeVertex (struct hedron *hed, long number) {
	long idx = 0;
	struct hedvtx *hvtx = NULL;
	struct phnvtx *pvtx = NULL;
	idx = number - 1;
	hvtx = getVertex (hed, number);
	if (hvtx == NULL) return;
	/* careful */
	pvtx = hvtx -> vtx;
	free_phnvtx (pvtx);
	free_hedvtx (hvtx);
	*(hed -> vertices + idx) = NULL;
}

void removeEdge (struct hedron *hed, long number) {
	long idx = 0;
	struct hededg *hedg = NULL;
	struct phnedg *pedg = NULL;
	idx = number - 1;
	hedg = getEdge (hed, number);
	if (hedg == NULL) return;
	pedg = hedg -> edg;
	free_phnedg (pedg);
	free_hededg (hedg);
	*(hed -> edges + idx) = NULL;
}

void removeTriangle (struct hedron *hed, long number) {
	long idx = 0;
	struct hedtri *htri = NULL;
	struct phntri *ptri = NULL;
	idx = number - 1;
	htri = getTriangle (hed, number);
	if (htri == NULL) return;
	ptri = htri -> tri;
	free_phntri (ptri);
	free_hedtri (htri);
	*(hed -> triangles + idx) = NULL;
}

/* allocate memory for new hedron */
struct hedron *newHedron (long maxvtx, long maxedg, long maxtri) {
	struct hedron *hed = NULL;
    struct cept *ex;

	hed = (struct hedron *) allocate_object (HEDRON);
	if (hed == NULL) {
		ex = new_cept ( MEMORY_ERROR,  ALLOCATION,  FATAL_SEVERITY);
        add_object (ex, HEDRON, "hed");
        add_function (ex, "newHedron");
        add_source (ex, "mshedron.c");
		return(NULL);
	}
	hed -> maxvtx = maxvtx;
	hed -> maxedg = maxedg;
	hed -> maxtri = maxtri;
	hed -> nvtx = hed -> maxvtx;
	hed -> nedg = hed -> maxedg;
	hed -> ntri = hed -> maxtri;
	hed -> vertices = (struct hedvtx **)
		allocate_pointers (HEDVTX, hed -> maxvtx);
	hed -> edges = (struct hededg **)
		allocate_pointers (HEDEDG, hed -> maxedg);
	hed -> triangles = (struct hedtri **)
		allocate_pointers (HEDTRI, hed -> maxtri);
	return (hed);
}

/* add a vertex to the hedron */
void addVertex (struct hedron *hed, long number, struct phnvtx *pvtx) {
	long idx = 0;
	struct hedvtx *hvtx = NULL;
	idx = number - 1;
	hvtx = allocate_hedvtx ();
	if (hvtx == NULL) return;
	*(hed -> vertices + idx) = hvtx;
	/* careful */
	hvtx -> vtx = pvtx;
	hvtx -> number = number;
}

/* add an edge to the hedron */
void addEdge (struct hedron *hed, long number, struct phnedg *pedg) {
	long idx = 0;
	int j;
	struct hededg *hedg = NULL;
	idx = number - 1;
	hedg = allocate_hededg ();
	if (hedg == NULL) return;
	*(hed -> edges + idx) = hedg;
	hedg -> edg = pedg;
	hedg -> number = number;
	for (j = 0; j < 2; j++)
		hedg -> vtxnum[j] = pedg -> vtxnum[j];
}

/* add a triangle to the hedron */
void addTriangle (struct hedron *hed, long number, struct phntri *ptri) {
	int j;
	long idx = 0;
	struct hedtri *htri = NULL;
	idx = number - 1;
	htri = allocate_hedtri ();
	if (htri == NULL) return;
	*(hed -> triangles + idx) = htri;
	htri -> tri = ptri;
	htri -> number = number;
	for (j = 0; j < 3; j++)
		htri -> vtxnum[j] = ptri -> vtxnum[j];
	for (j = 0; j < 3; j++)
		htri -> edgnum[j] = ptri -> edgnum[j];
}



