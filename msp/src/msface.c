#include "ms.h"
#include "msenum.h"
#include "mslimit.h"
#include "msdefault.h"
#include "msstruct.h"
#include "msglobal.h"
#include "msfunction.h"

/* Copyright 1995 by Michael L. Connolly */
/* Last revised: December 16, 2001 */


/* CONCAVE FACE ROUTINES */

/* new concave face */
void new_concave_face (struct surface *this_srf, struct probe *prb)
{
    int k;
    double probe_center[3];
	double vertex1_coor[3], vertex2_coor[3], vertex3_coor[3], vertex4_coor[3];
	double circle1_axis[3], circle2_axis[3], circle3_axis[3], circle4_axis[3];
	struct circle *circle1, *circle2, *circle3, *circle4;
	struct vertex *vertex1, *vertex2, *vertex3, *vertex4;
	struct arc *arc1, *arc2, *arc3, *arc4;
    struct sphere *atom1, *atom2, *atom3, *atom4;
	struct face *concave_fac;
    struct pair *torus1, *torus2, *torus3, *torus4;
	double dot1, dot2, dot3, dot4;
	char message[MAX_STRING];
    struct cept *ex;

	if (prb -> natom > 3) {
		sprintf (message, "%8ld probe square concave face", prb -> number);
		informd (message);
	}
	else {
		sprintf (message, "%8ld probe triangular concave face", prb -> number);
		informd2 (message);
	}
    atom1 = prb -> atm[0];
    atom2 = prb -> atm[1];
    atom3 = prb -> atm[2];
    atom4 = prb -> atm[3];
    torus1 = prb -> pairs[0];
    torus2 = prb -> pairs[1];
    torus3 = prb -> pairs[2];
    torus4 = prb -> pairs[3];
    for (k = 0; k < 3; k++)
		probe_center[k] = prb -> center[k];

	/* compute vertex coordinates */
	for (k = 0; k < 3; k++) {
		vertex1_coor[k] = (atom1 -> radius * probe_center[k] +
			this_srf -> probe_radius * atom1 -> center[k]) /
			(atom1 -> radius + this_srf -> probe_radius);
		vertex2_coor[k] = (atom2 -> radius * probe_center[k] +
			this_srf -> probe_radius * atom2 -> center[k]) /
			(atom2 -> radius + this_srf -> probe_radius);
		vertex3_coor[k] = (atom3 -> radius * probe_center[k] +
			this_srf -> probe_radius * atom3 -> center[k]) /
				(atom3 -> radius + this_srf -> probe_radius);
		if (atom4 != NULL) {
			vertex4_coor[k] = (atom4 -> radius * probe_center[k] +
				this_srf -> probe_radius * atom4 -> center[k]) /
					(atom4 -> radius + this_srf -> probe_radius);
		}
	}

	/* set up vertices */
	vertex1 = new_vertex (vertex1_coor, (struct sphere *) atom1, prb, NULL, NULL);
	if (vertex1 == NULL) return;
	link_vertex (this_srf, vertex1);
	vertex2 = new_vertex (vertex2_coor, (struct sphere *) atom2, prb, NULL, NULL);
	if (vertex2 == NULL) return;
	link_vertex (this_srf, vertex2);
	vertex3 = new_vertex (vertex3_coor, (struct sphere *) atom3, prb, NULL, NULL);
	if (vertex3 == NULL) return;
	link_vertex (this_srf, vertex3);
	if (atom4 != NULL) {
		vertex4 = new_vertex (vertex4_coor, (struct sphere *) atom4, prb, NULL, NULL);
		if (vertex4 == NULL) return;
		link_vertex (this_srf, vertex4);
	}
	else vertex4 = NULL;

	/* calculate axes and set up circles */
	setup_axis (probe_center, prb -> atm[0] -> center,
		prb -> atm[1] -> center, circle1_axis);
	setup_axis (probe_center, prb -> atm[1] -> center,
		prb -> atm[2] -> center, circle2_axis);
	if (atom4 == NULL) {
		setup_axis (probe_center, prb -> atm[2] -> center,
			prb -> atm[0] -> center, circle3_axis);
	}
	else if (atom4 != NULL) {
		setup_axis (probe_center, prb -> atm[2] -> center,
			prb -> atm[3] -> center, circle3_axis);
		setup_axis (probe_center, prb -> atm[3] -> center,
			prb -> atm[0] -> center, circle4_axis);
	}

	circle1 = new_circle (probe_center, this_srf -> probe_radius, circle1_axis);
	if (circle1 == NULL) return;
	link_circle (this_srf, circle1);
	circle1 -> theta = 0.0;
	circle1 -> subtype = GREAT_SUBTYPE;
	circle2 = new_circle (probe_center, this_srf -> probe_radius, circle2_axis);
	if (circle2 == NULL) return;
	link_circle (this_srf, circle2);
	circle2 -> theta = 0.0;
	circle2 -> subtype = GREAT_SUBTYPE;
	circle3 = new_circle (probe_center, this_srf -> probe_radius, circle3_axis);
	if (circle3 == NULL) return;
	link_circle (this_srf, circle3);
	circle3 -> theta = 0.0;
	circle3 -> subtype = GREAT_SUBTYPE;
	if (atom4 != NULL) {
		circle4 = new_circle (probe_center, this_srf -> probe_radius, circle4_axis);
		if (circle4 == NULL) return;
		link_circle (this_srf, circle4);
		circle4 -> theta = 0.0;
		circle4 -> subtype = GREAT_SUBTYPE;
	}
	else circle4 = NULL;

	/* set up arcs */
	arc1 = new_arc (circle1, vertex1, vertex2, CONCAVE, 0, (double) 0.0, 0L, 0L, 0L);
	if (arc1 == NULL) return;
	this_srf -> n_arc++;
	arc2 = new_arc (circle2, vertex2, vertex3, CONCAVE, 0, (double) 0.0, 0L, 0L, 0L);
	if (arc2 == NULL) return;
	this_srf -> n_arc++;
	if (circle4 == NULL) {
		arc3 = new_arc (circle3, vertex3, vertex1, CONCAVE, 0, (double) 0.0, 0L, 0L, 0L);
		if (arc3 == NULL) return;
		this_srf -> n_arc++;
		arc4 = NULL;
	}
	else if (circle4 != NULL) {
		arc3 = new_arc (circle3, vertex3, vertex4, CONCAVE, 0, (double) 0.0, 0L, 0L, 0L);
		if (arc3 == NULL) return;
		this_srf -> n_arc++;
		arc4 = new_arc (circle4, vertex4, vertex1, CONCAVE, 0, (double) 0.0, 0L, 0L, 0L);
		if (arc4 == NULL) return;
		this_srf -> n_arc++;
	}

	/* add arcs to tori */
	arc1 -> next = torus1 -> first_arc;
	torus1 -> first_arc = arc1;
	sprintf (message, "%8ld %8ld torus: add arc", 
		torus1 -> sph[0] -> number, torus1 -> sph[1] -> number);
	informd2(message);
	arc2 -> next = torus2 -> first_arc;
	torus2 -> first_arc = arc2;
	sprintf (message, "%8ld %8ld torus: add arc", 
		torus2 -> sph[0] -> number, torus2 -> sph[1] -> number);
	informd2(message);
	arc3 -> next = torus3 -> first_arc;
	torus3 -> first_arc = arc3;
	sprintf (message, "%8ld %8ld torus: add arc", 
		torus3 -> sph[0] -> number, torus3 -> sph[1] -> number);
	informd2(message);
	if (torus4 != NULL) {
		arc4 -> next = torus4 -> first_arc;
		torus4 -> first_arc = arc4;
		sprintf (message, "%8ld %8ld torus: add arc", 
			torus4 -> sph[0] -> number, torus4 -> sph[1] -> number);
		informd(message);
	}

	/* new concave face */
	concave_fac = new_face (NULL, CONCAVE);
	if (concave_fac == NULL) return;
	concave_fac -> ptr.prb = prb;
	concave_fac -> chi = 1;
	link_face (this_srf, concave_fac);
	add2face (concave_fac, arc1, arc2, arc3, arc4);

	/* back pointer for cusp corrections */
	prb -> fac = concave_fac;
	arc1 -> fac = concave_fac;
	arc2 -> fac = concave_fac;
	arc3 -> fac = concave_fac;
	if (arc4 != NULL) {
		arc4 -> fac = concave_fac;
	}
}




/* CONVEX SURFACE ROUTINES: */

void create_convex (struct surface *this_srf)
{
	int nc = 0;
	long number1;
	char message[MAX_STRING];
	struct sphere *atm;
	struct face *convex_fac;

	/* atoms with convex arcs are not buried */
	for (atm = this_srf -> head_atom; atm != NULL; atm = atm -> next)
		if (atm -> first_arc != NULL)
			atm -> buried = FALSE;

	/* non-buried atoms have accessible areas */
	/* to compute the Euler characteristic, we must count the cycles */
	for (atm = this_srf -> head_atom; atm != NULL; atm = atm -> next)
		if (!atm -> buried) {
			if (atm -> problem) continue;
			convex_fac = pre_form (this_srf, atm); if (error()) return;
			if (atm -> problem) continue;
			nc = form_cycles (convex_fac, NULL);
			if (nc < 0) {
				number1 = atm -> number;
				sprintf (message,
					"(create_convex) cycle does not close for atom %5ld",
					number1);
				inform(message);
				/* mark as a problem face */
				atm -> problem = TRUE;
				convex_fac -> problem = TRUE;
			}
			this_srf -> n_cycle += nc;
			if (error()) return;

			if (atm -> problem) continue;
			group_cycles (this_srf, convex_fac, NULL); if (error()) return;
		}
}


/* SADDLE SURFACE ROUTINES: */

/* create saddle surfaces */
void create_saddles (struct surface *this_srf)
{
	long n_free_saddle;
	char message[MAXLINE];
	struct torus *torus_ptr;

	n_free_saddle = 0;
	for (torus_ptr = this_srf -> head_torus; torus_ptr != NULL;
		torus_ptr = torus_ptr -> next)
		if (torus_ptr -> free) {
			free_saddle (this_srf, torus_ptr);			/* no concave edges */
			n_free_saddle++;
		}
		else
			regular_saddle (this_srf, torus_ptr);			/* regular saddles */
	sprintf (message, "%8ld free saddle faces", n_free_saddle);
	if (n_free_saddle > 0) inform (message);
}

/* free (no collision) torus saddle surface */
void free_saddle (struct surface *this_srf, struct torus *torus_ptr)
{
	struct face *saddle_fac;
	struct circle *circle1, *circle2;
	struct arc *arc1, *arc2;
	struct sphere *atom1, *atom2;

	atom1 = (torus_ptr -> atm[0]);
	atom2 = (torus_ptr -> atm[1]);

	/* set up circles */

	circle1 = new_contact_circle (this_srf, torus_ptr, 0); if (error()) return;
	circle2 = new_contact_circle (this_srf, torus_ptr, 1); if (error()) return;

	/* set up arcs */
	arc1 = new_arc (circle2, NULL, NULL, CONVEX, 0, (double) 0.0, 0L, 0L, 0L);
	if (error()) return;
	this_srf -> n_arc++;
	arc2 = new_arc (circle1, NULL, NULL, CONVEX, 0, (double) 0.0, 0L, 0L, 0L);
	if (error()) return;
	this_srf -> n_arc++;
	arc1 -> phi = 2 * PI;
	arc2 -> phi = 2 * PI;

	/* add convex arcs to atom list */
	arc1 -> next = atom2 -> first_arc;
	atom2 -> first_arc = arc1;
	arc2 -> next = atom1 -> first_arc;
	atom1 -> first_arc = arc2;

	/* allocate saddle face */
	saddle_fac = new_face (NULL, SADDLE); if (error()) return;
	saddle_fac -> ptr.tor = torus_ptr;
	link_face (this_srf, saddle_fac); if (error()) return;
	add2face (saddle_fac, arc1, arc2, NULL, NULL); if (error()) return;
}

/* one new saddle face with 4 edges */
void newsad (struct surface *this_srf, struct arc *arc1, struct arc *arc2, struct circle *circle1, struct circle *circle2, struct torus *torus_ptr, double wrap_angle)
{
	struct vertex *vertex1, *vertex2, *vertex3, *vertex4;
	struct arc *arc3, *arc4;
	struct face *saddle_fac;
	struct sphere *atom1, *atom2;

	atom1 = circle1 -> atm;
	atom2 = circle2 -> atm;

	/* gather vertices */
	vertex1 = arc1 -> vtx[0];
	vertex2 = arc1 -> vtx[1];
	vertex3 = arc2 -> vtx[0];
	vertex4 = arc2 -> vtx[1];

	/* set up arcs */
	arc3 = new_arc (circle2, vertex2, vertex3, CONVEX, 0, (double) 0.0, 0L, 0L, 0L);
	if (error()) return;
	this_srf -> n_arc++;
	arc4 = new_arc (circle1, vertex4, vertex1, CONVEX, 0, (double) 0.0, 0L, 0L, 0L);
	if (error()) return;
	this_srf -> n_arc++;
	arc3 -> phi = wrap_angle;
	arc4 -> phi = wrap_angle;

	/* add convex arcs to atom list */
	arc3 -> next = atom2 -> first_arc;
	atom2 -> first_arc = arc3;
	arc4 -> next = atom1 -> first_arc;
	atom1 -> first_arc = arc4;

	/* set up saddle face */
	saddle_fac = new_face (NULL, SADDLE); if (error()) return;
	saddle_fac -> ptr.tor = torus_ptr;
	link_face (this_srf, saddle_fac); if (error()) return;
	add2face (saddle_fac, arc1, arc3, arc2, arc4); if (error()) return;
}

/* regular saddle surfaces for this torus */
void regular_saddle (struct surface *this_srf, struct torus *torus_ptr)
{
	int k, n_torus_arcs, arc_index;
	int index0, index1;
	struct arc *arc1, *arc2;
	struct circle *circle1, *circle2;
	double circle_points[MAX_SORT][3];
	double point_vector1[3];
	double point_vector2[3];
	short arc_orn[MAX_SORT];
	short indices[MAX_SORT];
	struct arc **arc_list, **arc_hdl;
	double wrap_angle;
    struct cept *ex;

	/* set up circles */

	circle1 = new_contact_circle (this_srf, torus_ptr, 0); if (error()) return;
	circle2 = new_contact_circle (this_srf, torus_ptr, 1); if (error()) return;

	/* count arcs belonging to torus */
	n_torus_arcs = 0;
	for (arc1 = torus_ptr -> first_arc; arc1 != NULL;
		arc1 = arc1 -> next) n_torus_arcs++;
	if (n_torus_arcs <= 0) return;

	/* must be even */
	if (0 != n_torus_arcs % 2) {
        ex = new_cept (GEOMETRY_ERROR, INCONSISTENCY, FATAL_SEVERITY);
        add_function (ex, "regular_saddle");
        add_source (ex, "msface.c");
        add_long (ex, "n_torus_arcs", n_torus_arcs);
        add_atom (ex, torus_ptr -> atm[0]);
        add_atom (ex, torus_ptr -> atm[1]);
		add_message (ex, "odd number of torus arcs");
		return;
	}
	if (n_torus_arcs > MAX_SORT) {
        ex = new_cept (GEOMETRY_ERROR, MSOVERFLOW, FATAL_SEVERITY);
        add_function (ex, "regular_saddle");
        add_source (ex, "msface.c");
        add_long (ex, "n_torus_arcs", n_torus_arcs);
        add_long (ex, "MAX_SORT", MAX_SORT);
        add_atom (ex, torus_ptr -> atm[0]);
        add_atom (ex, torus_ptr -> atm[1]);
		return;
	}

	/* allocate memory */
	arc_list = (struct arc **)
		allocate_pointers (ARC, n_torus_arcs);
	if (arc_list == NULL) {
        ex = new_cept (MEMORY_ERROR, ALLOCATION, FATAL_SEVERITY);
        add_function (ex, "regular_saddle");
        add_source (ex, "msface.c");
        add_atom (ex, torus_ptr -> atm[0]);
        add_atom (ex, torus_ptr -> atm[1]);
		return;
	}

	/* set up torus arc pointer list */
	for (arc1 = torus_ptr -> first_arc, arc_index = 0;
		arc_index < n_torus_arcs;
		arc1 = arc1 -> next, arc_index++) {
		arc_hdl = arc_list + arc_index;
		*arc_hdl = arc1;			/* pointer to arc */
	}

	setup_torus_arcs (torus_ptr -> axis, n_torus_arcs, arc_list,
		circle_points, arc_orn);
	if (error()) return;

	sort_points (torus_ptr -> center, torus_ptr -> radius,
		torus_ptr -> axis, n_torus_arcs, circle_points,
		arc_orn, indices);
	if (error()) return;

	for (arc_index = 0; arc_index < n_torus_arcs; arc_index += 2) {
		index0 = indices[arc_index];
		arc_hdl = arc_list + index0;
		arc1 = *arc_hdl;
		index1 = indices[arc_index+1];
		arc_hdl = arc_list + index1;
		arc2 = *arc_hdl;

		/* compute saddle wrap angle */
		for (k = 0; k < 3; k++) {
			point_vector1[k] =
				(circle_points[index0][k] - torus_ptr -> center[k]) /
					torus_ptr -> radius;
			point_vector2[k] =
				(circle_points[index1][k] - torus_ptr -> center[k]) /
					torus_ptr -> radius;
		}
		wrap_angle = positive_angle (point_vector1, point_vector2,
			torus_ptr -> axis);
		/* create saddle face */
		newsad (this_srf, arc1, arc2, circle1, circle2, torus_ptr, wrap_angle);
		if (error()) return;
	}

	/* free temporary memory */
	free_pointers (ARC, arc_list);
}

void setup_torus_arcs (double axis[3], int n_torus_arcs, struct arc **arc_list, double circle_points[][3], short arc_orn[])
{
	int arc_index, k;
	struct arc **arc_hdl;
	struct arc *arc1;
	double vertex_vector[3];

	for (arc_index = 0; arc_index < n_torus_arcs; arc_index++) {
		arc_hdl = arc_list + arc_index;
		arc1 = *arc_hdl;			/* pointer to arc */
		for (k = 0; k < 3; k++) {
			circle_points[arc_index][k] = arc1 -> cir -> center[k];
			vertex_vector[k] = arc1 -> vtx[1] -> atm -> center[k] -
				arc1 -> vtx[0] -> atm -> center[k];
		}
		arc_orn[arc_index] = negdot (vertex_vector, axis);
	}
}


/* contact circle routine: */

struct circle *new_contact_circle (struct surface *this_srf, struct torus *torus_ptr, int which_side) 
{
	int k;
	double signed_distance, circle_radius;
	double circle_center[3], circle_atom_vector[3], circle_axis[3];
	struct circle *circle_ptr;
	struct sphere *atm_ptr;

	atm_ptr = torus_ptr -> atm[which_side];

	/* computations for circle */
	circle_radius =
		torus_ptr -> radius * atm_ptr -> radius /
		(atm_ptr -> radius + this_srf -> probe_radius);
	for (k = 0; k < 3; k++) {
		circle_center[k] =
			(atm_ptr -> radius * torus_ptr -> center[k] +
			this_srf -> probe_radius * atm_ptr -> center[k])
			/ (atm_ptr -> radius + this_srf -> probe_radius);
		circle_atom_vector[k] =
			circle_center[k] - atm_ptr -> center[k];
		circle_axis[k] = (2 * which_side - 1) * torus_ptr -> axis[k];
	}

	/* allocate memory, setup fields */
	circle_ptr = new_circle (circle_center, circle_radius, circle_axis);
	if (error()) return(NULL);
	link_circle (this_srf, circle_ptr);
	circle_ptr -> subtype = CONTACT_SUBTYPE;
	signed_distance = (-dot_product (circle_axis, circle_atom_vector));
	circle_ptr -> theta = atan2 (signed_distance, circle_radius);
	circle_ptr -> atm = atm_ptr;
	torus_ptr -> cir[which_side] = circle_ptr;
	return (circle_ptr);
}

int add2face (struct face *fac, struct arc *arc0, struct arc *arc1, struct arc *arc2, struct arc *arc3)
{
	fac -> n_cycle = 0;
	if (arc0 == NULL) return (0);
	fac -> arcsp[0] = arc0;
	fac -> n_arc++;
	if (arc1 == NULL) return (1);
	fac -> arcsp[1] = arc1;
	fac -> n_arc++;
	if (arc2 == NULL) return (2);
	fac -> arcsp[2] = arc2;
	fac -> n_arc++;
	if (arc3 == NULL) return (3);
	fac -> arcsp[3] = arc3;
	fac -> n_arc++;
	return (4);
}
