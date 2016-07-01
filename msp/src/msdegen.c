/*

	MSRoll

	Copyright 1986, 1989, 1996 by Michael L. Connolly
	All rights reserved

	Written by Michael L. Connolly.
	March 6, 2000

*/

#include "ms.h"
#include "msenum.h"
#include "mslimit.h"
#include "msdefault.h"
#include "msstruct.h"
#include "msglobal.h"
#include "msfunction.h"


void concave_degenerate (struct surface *this_srf, struct face *face_ptr)
{
	int i;
	struct edge *e;
	struct arc *a;
	struct cycle *cyc;
	char message[MAX_STRING];

	cyc = face_ptr -> first_cycle;
	sprintf(message, "concave_degenerate: ofn = %8ld", face_ptr -> ofn);
	informd2(message);
	/* kludge to prevent null original face number */
	for (i = 0, e = cyc -> first_edge; e != NULL; e = e -> next, i++) {
		a = e -> arcptr;
		if (a -> ofn <= 0) a -> ofn = face_ptr -> ofn;
	}
	if (!convert_face (this_srf, face_ptr)) return;
}

void saddle_degenerate (struct surface *this_srf, struct face *face_ptr)
{
	int orn0, orn1, orn2, orn3, ntie, ntoe, nbox;
	long n_e, n_e0, n_e1, i, j, m, n;
	double d = 0.0;
	struct edge *e, *already;
	struct edge *edge0, *edge1, *edge2, *edge3;
	struct arc *a;
	struct arc *arc0, *arc1;
	struct arc *arc2, *arc3;
	struct cycle *cyc0, *cyc1;
	struct vertex *vtx0, *vtx1;
	struct vertex *vtx2, *vtx3;
	struct face **face_list;
	struct vertex **vertex_list;
	struct edge **edge_list0, **tie_list;
	long *used;
	char message[MAX_STRING];

	cyc0 = face_ptr -> first_cycle;
	n_e0 = edges_in_cycle (cyc0);
	cyc1 = cyc0 -> next;
	if (cyc1 != NULL) {
		/* saddle hoop */
		n_e1 = edges_in_cycle(cyc1);
	}
	else n_e1 = 0;
	n_e = n_e0 + n_e1;

	sprintf(message, "saddle_degenerate: ofn = %8ld", face_ptr -> ofn);
	informd2(message);
	/* kludge to prevent null original face number */
	for (i = 0, e = cyc0 -> first_edge; e != NULL; e = e -> next, i++) {
		a = e -> arcptr;
		if (a -> ofn <= 0) a -> ofn = face_ptr -> ofn;
	}
	/* allocate local arrays */
	used = (long *) allocate_longs ((unsigned) n_e, 0, USED);
	if (used == NULL) {
		set_error1 ("(saddle_degenerate): mem allocation failure");
		return;
	}
	face_list = (struct face **)
		allocate_pointers (FACE, n_e);
	if (face_list == NULL) {
		set_error1 ("(saddle_degenerate): mem allocation failure");
		return;
	}
	edge_list0 = (struct edge **)
		allocate_pointers (EDGE, n_e);
	if (edge_list0 == NULL) {
		set_error1 ("(saddle_degenerate): mem allocation failure");
		return;
	}
	tie_list = (struct edge **)
		allocate_pointers (EDGE, 2 * n_e);
	if (tie_list == NULL) {
		set_error1 ("(saddle_degenerate): mem allocation failure");
		return;
	}
	vertex_list = (struct vertex **)
		allocate_pointers (VERTEX, n_e);
	if (vertex_list == NULL) {
		set_error1 ("(saddle_degenerate): mem allocation failure");
		return;
	}

	/* array of pre-existing edges */
	for (i = 0, edge0 = cyc0 -> first_edge; edge0 != NULL; i++, edge0 = edge0 -> next) {
		*(edge_list0+i) = edge0;
		if (edge0 == NULL) continue;
		arc0 = edge0 -> arcptr;
		orn0 = edge0 -> orn;
		vtx0 = arc0 -> vtx[orn0];
		vtx1 = arc0 -> vtx[1-orn0];
		sprintf (message, "pre edge %6ld %6ld (circumference)",
			vtx0 -> number, vtx1 -> number);
		informd2 (message);
	}
	if (cyc1 != NULL) {
		for (i = n_e0, edge0 = cyc1 -> first_edge; edge0 != NULL; i++, edge0 = edge0 -> next) {
			*(edge_list0+i) = edge0;
			if (edge0 == NULL) continue;
			arc0 = edge0 -> arcptr;
			orn0 = edge0 -> orn;
			vtx0 = arc0 -> vtx[orn0];
			vtx1 = arc0 -> vtx[1-orn0];
			sprintf (message, "pre edge %6ld %6ld (circumference2)",
				vtx0 -> number, vtx1 -> number);
			informd2 (message);
		}
	}

	/* create vertex array */
	for (i = 0; i < n_e; i++) {
		edge0 = *(edge_list0+i);
		arc0 = edge0 -> arcptr;
		orn0 = edge0 -> orn;
		vtx1 = arc0 -> vtx[orn0];
		*(vertex_list+i) = vtx1;
		if (!convert_vertex (this_srf, vtx1)) {
			informd("saddle_degenerate: failure to convert vertex");
			return;
		}
	}

	/* create tie edges */
	/* ending edges already exist for non-hoop */
	ntie = 0;
	ntoe = 0;
	for (i = 0; i < n_e-1; i++) {
		vtx0 = *(vertex_list+i);
		for (j = i+1; j < n_e; j++) {
			vtx1 = *(vertex_list+j);
			d = distance(vtx0 -> center, vtx1 -> center);
			sprintf (message, "distance %12.6f from %6ld to %6ld",
				d, vtx0 -> number, vtx1 -> number);
			informd2(message);
			if (d > 0.000001) continue;
			if (ntie >= n_e) {
				set_error1("saddle_degenerate: too many tie edges");
				return;
			}
			/* check whether it already exists */
			already = NULL;
			for (m = 0; m < n_e; m++) {
				edge0 = *(edge_list0+m);
				arc0 = edge0 -> arcptr;
				orn0 = edge0 -> orn;
				vtx2 = arc0 -> vtx[orn0];
				vtx3 = arc0 -> vtx[1-orn0];
				if (vtx0 == vtx2 && vtx1 == vtx3) already = edge0;
				if (vtx0 == vtx3 && vtx1 == vtx2) already = edge0;
			}
			if (already) {
				sprintf (message, "tie edge %6ld %6ld (already)",
					vtx0 -> number, vtx1 -> number);
				informd2 (message);
				*(tie_list+ntoe++) = already;
			}
			else {
				edge1 = make_straight_edge (this_srf, vtx0, vtx1);
				if (error()) return;
				sprintf (message, "tie edge %6ld %6ld (new)",
					vtx0 -> number, vtx1 -> number);
				informd2 (message);
				arc1 = edge1 -> arcptr;
				if (!convert_edge (this_srf, arc1)) {
					inform("saddle_degenerate: failure to convert edge");
				}
				orn1 = edge1 -> orn;
				edge2 = new_edge (edge1 -> arcptr, 1, (struct face *) NULL, NULL);
				*(tie_list+ntoe++) = edge1;
				*(tie_list+ntoe++) = edge2;
				ntie++;
			}
		}
	}
	sprintf (message, "%6d ties and %6d toes from %d edges", ntie, ntoe, n_e);
	informd (message);

	nbox = 0;
	for (i = 0; i < n_e; i++) {
		if (*(used + i)) continue;
		edge0 = *(edge_list0+i);
		arc0 = edge0 -> arcptr;
		orn0 = edge0 -> orn;
		vtx0 = arc0 -> vtx[orn0];
		vtx1 = arc0 -> vtx[1-orn0];
		for (j = 0; j < n_e; j++) {
			if (*(used + j)) continue;
			edge1 = *(edge_list0+j);
			arc1 = edge1 -> arcptr;
			orn1 = edge1 -> orn;
			vtx2 = arc1 -> vtx[orn1];
			vtx3 = arc1 -> vtx[1-orn1];
			if (vtx0 == vtx2 || vtx0 == vtx3) continue;
			if (vtx1 == vtx2 || vtx1 == vtx3) continue;
			for (m = 0; m < ntoe; m++) {
				edge2 = *(tie_list+m);
				arc2 = edge2 -> arcptr;
				orn2 = edge2 -> orn;
				if (arc2 -> vtx[orn2] != vtx1) continue;
				if (arc2 -> vtx[1-orn2] != vtx2) continue;
				for (n = 0; n < ntoe; n++) {
					if (m == n) continue;
					edge3 = *(tie_list+n);
					arc3 = edge3 -> arcptr;
					orn3 = edge3 -> orn;
					if (arc3 -> vtx[orn3] != vtx3) continue;
					if (arc3 -> vtx[1-orn3] != vtx0) continue;
					sprintf(message, "box: %6ld %6ld %6ld %6ld",
						vtx0 -> number, vtx1 -> number, vtx2 -> number, vtx3 -> number);
					informd2(message);
					make_two(this_srf, face_ptr, edge0, edge2, edge1, edge3, (nbox++ == 0));
					if (error()) return;
					*(used + i) = 1;
					*(used + j) = 1;
				}
			}
		}
	}

	free_pointers (FACE, face_list);
	free_pointers (EDGE, edge_list0);
	free_pointers (EDGE, tie_list);
	free_pointers (VERTEX, vertex_list);
	free_longs (used, 0, USED);

	/* simple_tent (face_ptr, 1); */
}

int make_two (struct surface *srf, struct face *fac, struct edge *e0, struct edge *e1, struct edge *e2, struct edge *e3, int use_given) {
	int orn0, orn2;
	struct face *f0 = NULL;
	struct face *f1 = NULL;
	struct edge *e4, *e5;
	struct arc *diagonal;
	struct vertex *vtx0, *vtx2;
	struct cycle *cyc = NULL;
	char message[MAX_STRING];

	orn0 = e0 -> orn;
	vtx0 = e0 -> arcptr -> vtx[orn0];
	orn2 = e2 -> orn;
	vtx2 = e2 -> arcptr -> vtx[orn2];
	e4 = make_straight_edge (srf, vtx0, vtx2);
	if (error()) return (0);
	sprintf (message, "tie edge %6ld %6ld (diagonal)",
		vtx0 -> number, vtx2 -> number);
	informd2 (message);
	diagonal = e4 -> arcptr;
	e5 = new_edge (diagonal, 1, (struct face *) NULL, NULL);
	if (use_given) {
		f0 = fac;
		cyc = fac -> first_cycle;
		cyc -> first_edge = e0;
		cyc -> next = NULL;
	}
	else {
		f0 = duplicate_face (fac);
		if (error()) return (0);
		cyc = new_cycle (NULL, e0);
		f0 -> first_cycle = cyc;
		if (error()) return (0);
		ink_cycle (srf);
	}
	/* modify links */
	e0 -> next = e1;
	e1 -> next = e5;
	e5 -> next = NULL;
	/* make edges point at triangle face */
	e0 -> fac = f0;
	e1 -> fac = f0;
	e5 -> fac = f0;
	f1 = duplicate_face (f0);
	if (error()) return (0);
	cyc = new_cycle (NULL, e2);
	f1 -> first_cycle = cyc;
	if (error()) return (0);
	ink_cycle (srf);
	/* modify links */
	e2 -> next = e3;
	e3 -> next = e4;
	e4 -> next = NULL;
	/* make edges point at triangle face */
	e2 -> fac = f1;
	e3 -> fac = f1;
	e4 -> fac = f1;
	if (!convert_face (srf, f0)) return(0);
	if (!convert_face (srf, f1)) return(0);
	return (1);
}


/*
   MSRoll
   Copyright 1996 by Michael L. Connolly
   All Rights Reserved

*/
