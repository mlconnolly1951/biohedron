/*

	MSRoll

	Copyright 1986, 1989, 1996 by Michael L. Connolly
	All rights reserved

	Written by Michael L. Connolly.
	March 7, 2000

*/

#include "ms.h"
#include "msenum.h"
#include "mslimit.h"
#include "msdefault.h"
#include "msstruct.h"
#include "msglobal.h"
#include "msfunction.h"


/* TENT routines */

int fac_centroid (struct face *fac, double center[3])
{
	int i, n, k;
	double sum[3], cyc_center[3];
	struct cycle *cyc;
	char message[MAXLINE];

	n = 0;
	for (k = 0; k < 3; k++)
		sum[k] = 0.0;
	for (cyc = fac -> first_cycle; cyc != NULL; cyc = cyc -> next) {
		n++;
		i = edges_in_cycle (cyc);
		if (i <= 0)  {
			set_error1 ("(fac_centroid): no edges for cycle");
			return (0);
		}
		if (!cyc_centroid (cyc, cyc_center)) return (0);
		for (k = 0; k < 3; k++)
			sum[k] += cyc_center[k];
	}
	if (n <= 0) {
		set_error1 ("(fac_centroid): no cycles");
		return (0);
	}
	for (k = 0; k < 3; k++)
		center[k] = sum[k] / n;
	return (1);
}

int cyc_centroid (struct cycle *cyc, double center[3])
{
	int n, k, orn;
	double sum[3];
	struct edge *e;
	struct arc *a;
	struct vertex *vtx;
	struct circle *cir;
	char message[MAXLINE];

	n = edges_in_cycle (cyc);
	if (n <= 0) {
		set_error1 ("(cyc_centroid): no edges");
		return (0);
	}
	for (k = 0; k < 3; k++)
		sum[k] = 0.0;

	for (e = cyc -> first_edge; e != NULL; e = e -> next) {
		a = e -> arcptr;
		cir = a -> cir;
		orn = e -> orn;
		vtx = a -> vtx[orn];
		if (vtx != (struct vertex *) NULL)
			for (k = 0; k < 3; k++)
				sum[k] += vtx -> center[k];
		else
			for (k = 0; k < 3; k++)
				sum[k] += cir -> center[k];
	}
	for (k = 0; k < 3; k++)
		center[k] = sum[k] / n;
	return (1);
}

int multiple_tent (struct face *given_face)
{
	double center[3];
	struct cycle *cyc;
	struct vertex *vtx0;
	struct surface *this_srf;

	this_srf = given_face -> srf;
	if (this_srf == NULL) return (0);

	cyc = given_face -> first_cycle;
	if (cyc == NULL)  {
		set_error1 ("(multiple_tent): face has no boundary");
		return (0);
	}
	if (!fac_centroid (given_face, center)) return (0);
	if (error()) return (0);
	vtx0 = new_vertex (center, NULL, NULL, (struct arc *) NULL, given_face);
	if (error()) return (0);
	link_vertex (this_srf, vtx0);
	if (!tent (given_face, vtx0)) return (0);
	if (error()) return (0);
	return (1);
}


/* subdivide the given face in a kludgy, but robust fashion */
/* the face is assumed to be simply connected */

int simple_tent (struct face *given_face, int troid)
{
	double center[3];
	struct cycle *cyc;
	struct vertex *vtx0;
	struct surface *this_srf;

	this_srf = given_face -> srf;
	if (this_srf == NULL) return (0);

	cyc = given_face -> first_cycle;
	if (cyc == NULL)  {
		set_error1 ("(simple_tent): face has no boundary");
		return (0);
	}
	if (cyc -> next != NULL) {
		set_error1 ("(simple_tent): > 1 cycle");
		return (0);
	}
	if (!troid && given_face -> shape == CONVEX) {
		if (!face_center (given_face, center)) {
			if (!fac_centroid (given_face, center)) return (0);
		}
	}
	else {
		if (!fac_centroid (given_face, center)) return (0);
	}
	if (error()) return (0);


	vtx0 = new_vertex (center, NULL, NULL, (struct arc *) NULL, given_face);
	if (error()) return (0);
	link_vertex (this_srf, vtx0);

	if (!tent (given_face, vtx0)) return (0);
	if (error()) return (0);
	return (1);
}

int tent (struct face *given_face, struct vertex *vtx0)
{
	int i, n;
	struct cycle *cyc;
	struct surface *this_srf;

	this_srf = given_face -> srf;
	if (this_srf == NULL) return (0);

	n = 0;
	for (cyc = given_face -> first_cycle; cyc != NULL; cyc = cyc -> next)
		n++;
	if (n <= 0)  {
		set_error1 ("(tent): face has no boundary");
		return (0);
	}
	for (i = 0, cyc = given_face -> first_cycle; cyc != NULL;
		i++, cyc = cyc -> next) {
		if (!cycle_tent (given_face, cyc, vtx0, (i == n - 1))) return (0);
		if (error()) return (0);
	}
	return (1);
}

int cycle_tent (struct face *given_face, struct cycle *given_cyc, struct vertex *vtx0, int use_given)
{
	int orn, orn0, orn1;
	long n_e, i, j;
	char message[MAXLINE];
	struct vertex *vtx1, *vtx2;
	struct arc *arc0, *arc1;
	struct edge *edge0, *edge1, *edge2;
	struct face *triangle_face;
	struct face **face_list;
	struct vertex **vertex_list;
	struct edge **edge_list0, **edge_list1, **edge_list2;
	struct surface *this_srf;

	this_srf = given_face -> srf;
	if (this_srf == NULL) return (0);

	if (given_cyc == NULL)  {
		set_error1 ("(cycle_tent): face has no boundary");
		return (0);
	}
	n_e = edges_in_cycle (given_cyc);
	if (error()) return (0);
	if (n_e < 3) {
		sprintf (message, "(cycle_tent): face with %d sides", n_e);
		set_error1 (message);
		return (0);
	}

	/* allocate local arrays */
	edge_list0 = (struct edge **)
		allocate_pointers (EDGE, n_e);
	if (edge_list0 == NULL) {
		set_error1 ("(cycle_tent): mem allocation failure");
		return (0);
	}
	edge_list1 = (struct edge **)
		allocate_pointers (EDGE, n_e);
	if (edge_list1 == NULL) {
		set_error1 ("(cycle_tent): mem allocation failure");
		return (0);
	}
	edge_list2 = (struct edge **)
		allocate_pointers (EDGE, n_e);
	if (edge_list2 == NULL) {
		set_error1 ("(cycle_tent): mem allocation failure");
		return (0);
	}
	face_list = (struct face **)
		allocate_pointers (FACE, n_e);
	if (face_list == NULL) {
		set_error1 ("(cycle_tent): mem allocation failure");
		return (0);
	}
	vertex_list = (struct vertex **)
		allocate_pointers (VERTEX, n_e);
	if (vertex_list == NULL) {
		set_error1 ("(cycle_tent): mem allocation failure");
		return (0);
	}

	for (i = 0, edge0 = given_cyc -> first_edge; edge0 != NULL; i++, edge0 = edge0 -> next) {
		*(edge_list0+i) = edge0;
		if (edge0 == NULL) continue;
		arc0 = edge0 -> arcptr;
		if (arc0 == NULL) {
			set_error1 ("(cycle_tent): null arcptr");
			return (0);
		}
		orn0 = edge0 -> orn;
	}

	/* create vertex array */
	for (i = 0; i < n_e; i++) {
		edge0 = *(edge_list0+i);
		arc0 = edge0 -> arcptr;
		orn = edge0 -> orn;
		vtx1 = arc0 -> vtx[orn];
		*(vertex_list+i) = vtx1;
	}
	for (i = 0; i < n_e-1; i++) {
		for (j = i + 1; j < n_e; j++) {
			vtx1 = *(vertex_list+i);
			vtx2 = *(vertex_list+j);
			if (vtx1 == vtx2 && ! this_srf -> van_der_Waals &&
				(j == i + 1 || i == 0 && j == n_e-1)) {
				sprintf (message, "cycle_tent: repeated vertex: %8ld", vtx1 -> number);
				set_error1(message);
				sprintf (message, "cycle_tent: i = %d, j = %d",
					i, j);
				set_error2(message);
				return (0);
			}
		}
	}

	/* create spoke edges */
	for (i = 0; i < n_e; i++) {
		vtx1 = *(vertex_list+i);
		edge1 = make_straight_edge (this_srf, vtx0, vtx1);
		if (error()) return (0);
		arc1 = edge1 -> arcptr;
		orn1 = edge1 -> orn;
		*(edge_list1+i) = edge1;
		*(edge_list2+i) = new_edge (edge1 -> arcptr, 1, (struct face *) NULL, NULL);
	}

	for (i = 0; i < n_e; i++) {
		edge0 = *(edge_list0+i);
		j = ((i < n_e - 1) ? i + 1 : 0);
		edge1 = *(edge_list1+i);
		edge2 = *(edge_list2+j);
		/* use given or duplicate face */
		if (use_given && i == 0) {
			triangle_face = given_face;
			given_cyc -> first_edge = edge1;
			triangle_face -> first_cycle = given_cyc;
		}
		else {
			triangle_face = duplicate_face (given_face);
			if (error()) return (0);
			triangle_face -> first_cycle = new_cycle (NULL, edge1);
			if (error()) return (0);
			ink_cycle (this_srf);
		}
		/* later: triangle_face -> shape = STRAIGHT; */
		*(face_list + i) = triangle_face;
		edge1 -> next = edge0;
		edge0 -> next = edge2;
		edge2 -> next = NULL;
	}

	/* edge -> face pointers */
	for (i = 0; i < n_e; i++) {
		edge0 = *(edge_list0+i);
		j = ((i < n_e - 1) ? i + 1 : 0);
		edge1 = *(edge_list1+i);
		edge2 = *(edge_list2+j);
		triangle_face = *(face_list+i);
		/* pointer from old edge to face */
		edge0 -> fac = triangle_face;
		/* pointers from new edges to faces */
		edge1 -> fac = triangle_face;
		edge2 -> fac = triangle_face;
	}
	
	for (i = 0; i < n_e; i++) {
		triangle_face = *(face_list+i);
		if (!convert_face (this_srf, triangle_face)) {
			inform("cycle_tent: failure to convert face to triangle");
			return(0);
		}
		if (error()) return (0);
	}

	free_pointers (EDGE, edge_list0);
	free_pointers (EDGE, edge_list1);
	free_pointers (EDGE, edge_list2);
	free_pointers (FACE, face_list);
	free_pointers (VERTEX, vertex_list);
	return (1);
}

/*
   MSRoll
   Copyright 1996 by Michael L. Connolly
   All Rights Reserved

*/
