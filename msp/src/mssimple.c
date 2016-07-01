/*

	MSRoll

	Copyright 1986, 1989, 1996 by Michael L. Connolly
	All rights reserved

	Written by Michael L. Connolly.
	January 8, 2002

*/

#include "ms.h"
#include "msenum.h"
#include "mslimit.h"
#include "msdefault.h"
#include "msstruct.h"
#include "msglobal.h"
#include "msfunction.h"


/* FACE */


/* make face simply connected */
int simplify (struct face *given_face)
{
	int result, n_point;
    long atmnum0, atmnum1, atmnum2, atmnum3;
	struct sphere *atm;
	struct probe *prb;
	struct arc *connecting_arc;
	struct cycle *cyc1, *cyc2;
	struct variety *vty;
	double omega_before, omega_after;
	struct surface *this_srf;
    struct cept *ex;

	vty = given_face -> vty;
    atmnum0 = vty -> atmnum[0];
    atmnum1 = vty -> atmnum[1];
    atmnum2 = vty -> atmnum[2];
    atmnum3 = vty -> atmnum[3];
	this_srf = given_face -> srf;
	if (this_srf == NULL) return (0);
	/* first cycle of face */
	cyc1 = given_face -> first_cycle;
	if (cyc1 == NULL) return (1);

	/* second cycle of face */
	cyc2 = cyc1 -> next;
	if (cyc2 == NULL) return (1);

	if (given_face -> shape == CONVEX) {
		atm = given_face -> ptr.atm;
	}
	else if (given_face -> shape == CONCAVE) {
		prb = given_face -> ptr.prb;
	}
	omega_before = compute_omega (vty -> radii[0], given_face);
	if (error()) return(0);

	/* find closest approach of two cycles */
	connecting_arc = closest_pair (given_face, cyc1, cyc2);
	if (error()) return (0);

	/* connect the two cycles by an arc */
	result = connect_pair (cyc1, cyc2, connecting_arc);
	if (error()) return (0);
	if (!result) return (0);

	omega_after = compute_omega (vty -> radii[0], given_face);
	if (error()) return (0);
	if (! equal (omega_after , omega_before)) {
		ex = new_cept (GEOMETRY_ERROR,  INCONSISTENCY,  FATAL_SEVERITY);
		add_function (ex, "simplify");
		add_source (ex, "mssimple.c");
        add_double (ex, "omega_before", omega_before);
        add_double (ex, "omega_after", omega_after);
        add_long (ex, "first  atom number", atmnum0);
        add_long (ex, "second atom number", atmnum1);
        add_long (ex, "third  atom number", atmnum2);
		if (atmnum3 > 0)
			add_long (ex, "fourth atom number", atmnum3);
		return (0);
	}

	/* Note: do not free cycle 2 memory */
	/* it was not independently allocated */
	/* Well, maybe that is not true anymore */
	free_cycle (cyc2);
	this_srf -> n_simplified_face++;

	n_point = subdivide_arc (this_srf, connecting_arc);		/* subdivide new arc */
	if (error()) return (0);
	if (n_point <= 0) return (0);

	/* recursion */

	if (!simplify (given_face)) return (0);
	if (error()) return (0);

	given_face -> simplified = 1;
	return (1);
}

/* connect cycle 1 and cycle 2 by an arc to make big cycle */
int connect_pair (struct cycle *cyc1, struct cycle *cyc2, struct arc *connecting_arc)
{
	int orn, after_count, cyc1_edge_count, cyc2_edge_count;
    long atmnum0, atmnum1, atmnum2, atmnum3;
	struct edge *e1, *e2, *e3, *e4, *e5, *e6;
	struct edge *e, *cyc1_end_edge, *cyc2_end_edge;
	struct vertex *vtx1, *vtx2;
	struct arc *a;
	struct face *fac;
    struct cept *ex;
    struct variety *vty;

	cyc1_edge_count = edges_in_cycle (cyc1);
	cyc2_edge_count = edges_in_cycle (cyc2);

	/* copy vertex pointers to local variables */
	vtx1 = connecting_arc -> vtx[0];
	vtx2 = connecting_arc -> vtx[1];

	/* get face */
	fac = cyc1 -> first_edge -> fac;
    vty = fac -> vty;
    atmnum0 = vty -> atmnum[0];
    atmnum1 = vty -> atmnum[1];
    atmnum2 = vty -> atmnum[2];
    atmnum3 = vty -> atmnum[3];

	/* initialization */
	e1 = e2 = e3 = e4 = NULL;

	/* look for edges adjacent to vertex 1 */

	for (e = cyc1 -> first_edge; e != NULL; e = e -> next) {
		a = e -> arcptr;
		orn = e -> orn;
		if (a -> vtx[orn] == vtx1) e2 = e;
		if (a -> vtx[1-orn] == vtx1) e1 = e;
	}

	/* look for edges in cycle 2 adjacent to connecting arc */
	for (e = cyc2 -> first_edge; e != NULL; e = e -> next) {
		a = e -> arcptr;
		orn = e -> orn;
		if (a -> vtx[orn] == vtx2) e4 = e;
		if (a -> vtx[1-orn] == vtx2) e3 = e;
	}

	/* allocate two new edges for connecting arc */
	e5 = new_edge (connecting_arc, 0, fac, NULL);
	if (e5 == NULL) return (0);
	e6 = new_edge (connecting_arc, 1, fac, NULL);
	if (e6 == NULL) return (0);

	/* find last edges of the two cycles */

	for (e = cyc1 -> first_edge; e != NULL; e = e -> next)
		if (e -> next == NULL) cyc1_end_edge = e;
	for (e = cyc2 -> first_edge; e != NULL; e = e -> next)
		if (e -> next == NULL) cyc2_end_edge = e;

	/* cyclize */

	cyc1_end_edge -> next = cyc1 -> first_edge;
	cyc2_end_edge -> next = cyc2 -> first_edge;

	/* change links to there will be one big cycle */
	e1 -> next = NULL;
	e5 -> next = e4;
	e3 -> next = e6;
	e6 -> next = e2;

	/* make first cycle point at new first edge */
	cyc1 -> first_edge = e5;

	cyc1 -> next = cyc2 -> next;		/* delink 2nd cycle of face */

	after_count = edges_in_cycle (cyc1);

	if (after_count != cyc1_edge_count + cyc2_edge_count + 2) {
		ex = new_cept (GEOMETRY_ERROR,  INCONSISTENCY,  FATAL_SEVERITY);
		add_function (ex, "connect_pair");
		add_source (ex, "mssimple.c");
        add_long (ex, "edge count after making two cycles", after_count);
        add_long (ex, "first cycle count", cyc1_edge_count);
        add_long (ex, "second cycle count", cyc2_edge_count);
        add_long (ex, "first  atom number", atmnum0);
        add_long (ex, "second atom number", atmnum1);
        add_long (ex, "third  atom number", atmnum2);
		if (atmnum3 > 0)
			add_long (ex, "fourth atom number", atmnum3);
		return (0);
	}
	return (1);
}

/* find closest pair of vertices for two cycles */
struct arc *closest_pair (struct face *given_face, struct cycle *cyc1, struct cycle *cyc2)
{
	int k, orn1, orn2, reverse_axis, arc_shape;
	int n_edge1, n_edge2, is_contained, orn;
	long lfn, ofn;
    long atmnum0, atmnum1, atmnum2, atmnum3;
	double minimum_distance, vertex_separation, circle_radius, dot_axis;
	double alpha, ang, delta1, delta2, delta3, delta4;
	double circle_center[3], circle_axis[3];
	double vector1[3], vector2[3], midpoint[3];
	double vtx1_vector[3], vtx2_vector[3];
	struct variety *vty;
	struct vertex *vtx1, *vtx2, *min_vtx1, *min_vtx2;
	struct edge *edg1, *edg2;
	struct edge *e, *e1, *e2, *e3, *e4;
	struct arc *arc1, *arc2, *a;
	struct arc *connecting_arc;
	struct circle *cir;
	struct circle temporary_circle;
	struct arc temporary_arc;
	struct edge temporary_edge1, temporary_edge2;
	struct surface *this_srf;
    struct cept *ex;

	/* store variety pointer in local variable */
	vty = given_face -> vty;
	this_srf = given_face -> srf;
	if (this_srf == NULL) return (0);
	alpha = given_face -> alpha;
	lfn = given_face -> lfn;
	ofn = given_face -> ofn;
    atmnum0 = vty -> atmnum[0];
    atmnum1 = vty -> atmnum[1];
    atmnum2 = vty -> atmnum[2];
    atmnum3 = vty -> atmnum[3];

	n_edge1 = edges_in_cycle (cyc1);
	n_edge2 = edges_in_cycle (cyc2);
	if (n_edge1 < 2 || n_edge2 < 2) {
		ex = new_cept (GEOMETRY_ERROR,  MSUNDERFLOW,  FATAL_SEVERITY);
		add_function (ex, "closest_pair");
		add_source (ex, "mssimple.c");
        add_long (ex, "n_edge1", n_edge1);
        add_long (ex, "n_edge2", n_edge2);
        add_long (ex, "first  atom number", atmnum0);
        add_long (ex, "second atom number", atmnum1);
        add_long (ex, "third  atom number", atmnum2);
		if (atmnum3 > 0)
			add_long (ex, "fourth atom number", atmnum3);
		return (0);
	}

	/* initialization */
	minimum_distance = 1000000.0;
	min_vtx1 = NULL;
	min_vtx2 = NULL;

	/* nexted loops check every pair of vertices */
	for (edg1 = cyc1 -> first_edge; edg1 != NULL; edg1 = edg1 -> next) {
		/* transfer to local variables */
		arc1 = edg1 -> arcptr;
		orn1 = edg1 -> orn;
		vtx1 = arc1 -> vtx[orn1];
		for (edg2 = cyc2 -> first_edge; edg2 != NULL; edg2 = edg2 -> next) {
			/* transfer to local variables */
			arc2 = edg2 -> arcptr;
			orn2 = edg2 -> orn;
			vtx2 = arc2 -> vtx[orn2];
			vertex_separation = distance (vtx1 -> center, vtx2 -> center);
			/* check whether arc intersects well */
			if (vertex_separation < minimum_distance &&
				given_face -> shape == CONVEX) {
				for (k = 0; k < 3; k++)
					circle_center[k] = vty -> center[k];
				circle_radius = vty -> radii[0];
				for (k = 0; k < 3; k++) {
					vector1[k] = vtx1 -> center[k] -
						circle_center[k];
					vector2[k] = vtx2 -> center[k] -
						circle_center[k];
				}
				cross (vector1, vector2, circle_axis);
				if (norm (circle_axis) <= 0.0) continue;
				normalize (circle_axis);
				for (k = 0; k < 3; k++) {
					temporary_circle.center[k] = vty -> center[k];
					temporary_circle.axis[k] = circle_axis[k];
				}
				temporary_circle.radius = vty -> radii[0];
				temporary_arc.cir = &temporary_circle;
				temporary_arc.vtx[0] = vtx1;
				temporary_arc.vtx[1] = vtx2;
				temporary_arc.shape = CONVEX;
				temporary_arc.small = 0;
				/* check for arc midpoint outside face */
				middle (&temporary_arc, midpoint);
				if (error()) return (NULL);
				is_contained = point_in_face (midpoint, given_face, 0);
				if (error()) return (NULL);
				if (is_contained != 1) continue;	/* no good */
				/* temporary edges for edge-edge angle checks */
				temporary_edge1.arcptr = &temporary_arc;
				temporary_edge1.orn = 0;
				temporary_edge1.fac = given_face;
				temporary_edge2.arcptr = &temporary_arc;
				temporary_edge2.orn = 1;
				temporary_edge2.fac = given_face;
				for (k = 0; k < 3; k++) {
					vtx1_vector[k] = (vtx1 -> center[k] -
						vty -> center[k]) / vty -> radii[0];
					vtx2_vector[k] = (vtx2 -> center[k] -
						vty -> center[k]) / vty -> radii[0];
				}
				/* initialization */
				e1 = e2 = e3 = e4 = NULL;
				/* look for edges adjacent to vertex 1 */
				for (e = cyc1 -> first_edge; e != NULL; e = e -> next) {
					a = e -> arcptr;
					orn = e -> orn;
					if (a -> vtx[orn] == vtx1) e2 = e;
					if (a -> vtx[1-orn] == vtx1) e1 = e;
				}
				/* look for edges in cycle 2
					adjacent to connecting arc */
				for (e = cyc2 -> first_edge; e != NULL; e = e -> next) {
					a = e -> arcptr;
					orn = e -> orn;
					if (a -> vtx[orn] == vtx2) e4 = e;
					if (a -> vtx[1-orn] == vtx2) e3 = e;
				}
				/* check for failure to find these 4 edges */
				if (e1 == NULL) continue;
				if (e2 == NULL) continue;
				if (e3 == NULL) continue;
				if (e4 == NULL) continue;
				delta1 = edge_delta (e1, &temporary_edge1, vtx1_vector);
				delta2 = edge_delta (&temporary_edge2, e2, vtx1_vector);
				delta3 = edge_delta (e3, &temporary_edge2, vtx2_vector);
				delta4 = edge_delta (&temporary_edge1, e4, vtx2_vector);
				if (delta1 <= 0.0 || delta1 >= PI) continue;
				if (delta2 <= 0.0 || delta2 >= PI) continue;
				if (delta3 <= 0.0 || delta3 >= PI) continue;
				if (delta4 <= 0.0 || delta4 >= PI) continue;
				/* end of meeting angle checks */
				ang = arc_ang (&temporary_arc);
				if (ang <= 0.0 || ang >= PI) continue;
				vertex_separation = circle_radius * ang;
			}

			/* check for closer pair */
			if (vertex_separation < minimum_distance) {
				/* store new minimum */
				minimum_distance = vertex_separation;
				min_vtx1 = vtx1;
				min_vtx2 = vtx2;
			}
		}
	}

	/* check for failure */
	if (min_vtx1 == NULL || min_vtx2 == NULL) {
		ex = new_cept (GEOMETRY_ERROR,  DEGENERACY,  FATAL_SEVERITY);
		add_function (ex, "closest_pair");
		add_source (ex, "mssimple.c");
        add_message (ex, "no minimum-distance pair of vertices found");
        add_long (ex, "first  atom number", atmnum0);
        add_long (ex, "second atom number", atmnum1);
        add_long (ex, "third  atom number", atmnum2);
		if (atmnum3 > 0)
			add_long (ex, "fourth atom number", atmnum3);
		return (0);
	}

	/* calculate center, radius, axis of circle */

	if (vty -> type == SPHERE) {
		for (k = 0; k < 3; k++) circle_center[k] = vty -> center[k];
		circle_radius = vty -> radii[0];
		for (k = 0; k < 3; k++) {
			vector1[k] = min_vtx1 -> center[k] - circle_center[k];
			vector2[k] = min_vtx2 -> center[k] - circle_center[k];
		}
		cross (vector1, vector2, circle_axis);
		normalize (circle_axis);
	}
	else if (vty -> type == TORUS) {
		for (k = 0; k < 3; k++) {
			midpoint[k] =
				(min_vtx1 -> center[k] + min_vtx2 -> center[k]) / 2;
			vector2[k] = min_vtx2 -> center[k] - min_vtx1 -> center[k];
		}
		reverse_axis = (dot_product (vector2, vty -> axis) < 0.0);
		for (k = 0; k < 3; k++)
			vector1[k] = midpoint[k] - vty -> center[k];
		dot_axis = dot_product (vector1, vty -> axis);
		for (k = 0; k < 3; k++)
			vector1[k] -= (dot_axis * vty -> axis[k]);
		normalize (vector1);
		for (k = 0; k < 3; k++)
			circle_center[k] =
				vty -> center[k] + vty -> radii[0] * vector1[k];
		circle_radius = this_srf -> probe_radius;
		cross (vty -> axis, vector1, circle_axis);
		normalize (circle_axis);
		if (reverse_axis)
			for (k = 0; k < 3; k++)
				circle_axis[k] = (-circle_axis[k]);
	}
	else if (vty -> type == CYLINDER) {
		for (k = 0; k < 3; k++) {
			circle_center[k] = (min_vtx1 -> center[k] +
				min_vtx2 -> center[k]) / 2;
			circle_axis[k] = 0.0;
		}
		circle_radius = 0.0;
	}
	else {
		set_error1 ("(closest_pair): invalid variety type");
		return (0);
	}

	/* new circle and arc */
	cir = new_circle (circle_center, circle_radius, circle_axis);
	if (cir == NULL) return (0);
	link_circle (this_srf, cir);
	if (given_face -> shape == CONVEX)
		arc_shape = CONVEX;
	else if (given_face -> shape == SADDLE)
		arc_shape = CONCAVE;
	else if (given_face -> shape == CONCAVE)
		arc_shape = CONCAVE;
	else if (given_face -> shape == CYLINDRICAL)
		arc_shape = STRAIGHT;

	connecting_arc = new_arc (cir, min_vtx1, min_vtx2, arc_shape, 0, alpha, lfn, ofn, ofn);
	if (error()) return (NULL);
	link_arc (this_srf, connecting_arc);
	return (connecting_arc);
}

/* special split routine for entire sphere convex faces */

struct face *entire_sphere (struct face *given_face)
{
	int result, n_point;
	long lfn, ofn;
    long atmnum0, atmnum1, atmnum2, atmnum3;
	double alpha, omega, omega1, omega2, omega3, omega4;
	double hemi_axis[3];
	struct face *second_face;
	struct cycle *cyc1, *cyc2;
	struct arc *hemi_arc;
	struct edge *edg1, *edg2;
	struct circle *hemi_circle;
	struct variety *vty;
	struct surface *this_srf;
    struct cept *ex;

	this_srf = given_face -> srf;
	if (this_srf == NULL) return (0);
	lfn = given_face -> lfn;
	ofn = given_face -> ofn;

	/* transfer info to local variables */

	vty = given_face -> vty;
	alpha = given_face -> alpha;
    atmnum0 = vty -> atmnum[0];
    atmnum1 = vty -> atmnum[1];
    atmnum2 = vty -> atmnum[2];
    atmnum3 = vty -> atmnum[3];

	omega = compute_omega (vty -> radii[0], given_face);
	if (error()) return (0);

	this_srf -> n_entire_face++;
	hemi_axis[0] = 0.0;
	hemi_axis[1] = 0.0;
	hemi_axis[2] = 1.0;

	/* allocate new circle and arc */
	hemi_circle = new_circle (vty -> center, vty -> radii[0], hemi_axis);
	if (error()) return (0);
	link_circle (this_srf, hemi_circle);
	hemi_arc = new_arc (hemi_circle, (struct vertex *) NULL,
		(struct vertex *) NULL, CONVEX, 0, alpha, lfn, ofn, ofn);
	if (error()) return (0);
	link_arc (this_srf, hemi_arc);

	/* need face pointer for newedg call */
	second_face = duplicate_face (given_face);
	if (error()) return (0);

	/* set okay to use add-circle routine (kludge) */
	given_face -> largeok = 1;
	second_face -> largeok = 1;

	/* two new edges */
	edg1 = new_edge (hemi_arc, 0, second_face, NULL);
	if (error()) return (0);
	edg2 = new_edge (hemi_arc, 1, given_face, NULL);
	if (error()) return (0);

	/* two new cycles */
	cyc1 = new_cycle (NULL, edg1);
	if (error()) return (0);
	ink_cycle (this_srf);
	cyc2 = new_cycle (NULL, edg2);
	if (error()) return (0);
	ink_cycle (this_srf);

	/* new face has just one cycle */
	second_face -> first_cycle = cyc1;

	/* old face has 1 cycle now */
	given_face -> first_cycle = cyc2;

	omega1 = compute_omega (vty -> radii[0], given_face);
	if (error()) return (0);
	omega2 = compute_omega (vty -> radii[0], second_face);
	if (error()) return (0);

	if (!equal (omega1 + omega2, omega)) {
		ex = new_cept (GEOMETRY_ERROR,  INCONSISTENCY,  FATAL_SEVERITY);
		add_function (ex, "connect_pair");
		add_source (ex, "mssimple.c");
        add_double (ex, "omega1", omega1);
        add_double (ex, "omega2", omega2);
        add_double (ex, "omega1+omega2", omega1+omega2);
        add_double (ex, "omega", omega);
        add_long (ex, "first  atom number", atmnum0);
        add_long (ex, "second atom number", atmnum1);
        add_long (ex, "third  atom number", atmnum2);
		if (atmnum3 > 0)
			add_long (ex, "fourth atom number", atmnum3);
		return(0);
	}


	n_point = subdivide_arc (this_srf, hemi_arc);		/* subdivide new arc */
	if (error()) return (0);
	if (n_point <= 0) return (0);

	omega3 = compute_omega (vty -> radii[0], given_face);
	if (error()) return (0);
	omega4 = compute_omega (vty -> radii[0], second_face);
	if (error()) return (0);
	if (!equal (omega1, omega3)) {
		ex = new_cept (GEOMETRY_ERROR,  INCONSISTENCY,  FATAL_SEVERITY);
		add_function (ex, "connect_pair");
		add_source (ex, "mssimple.c");
        add_double (ex, "omega1", omega1);
        add_double (ex, "omega3", omega3);
        add_double (ex, "omega1+omega3", omega1+omega3);
        add_double (ex, "omega", omega);
        add_long (ex, "first  atom number", atmnum0);
        add_long (ex, "second atom number", atmnum1);
        add_long (ex, "third  atom number", atmnum2);
		if (atmnum3 > 0)
			add_long (ex, "fourth atom number", atmnum3);
 		return (0);
	}
	if (!equal (omega2, omega4)) {
		ex = new_cept (GEOMETRY_ERROR,  INCONSISTENCY,  FATAL_SEVERITY);
		add_function (ex, "connect_pair");
		add_source (ex, "mssimple.c");
        add_double (ex, "omega2", omega2);
        add_double (ex, "omega4", omega4);
        add_double (ex, "omega2+omega4", omega2+omega4);
        add_double (ex, "omega", omega);
        add_long (ex, "first  atom number", atmnum0);
        add_long (ex, "second atom number", atmnum1);
        add_long (ex, "third  atom number", atmnum2);
		if (atmnum3 > 0)
			add_long (ex, "fourth atom number", atmnum3);
 		return (0);
	}

	return (second_face);	/* return pointer to new face */
}

/* special split routine for large or troublesome convex faces */

struct face *circle_face (struct face *given_face)
{
	int k, is_contained, result, n_point;
	long lfn, ofn;
    long atmnum0, atmnum1, atmnum2, atmnum3;
	double alpha, circle_radius, cone_angle;
	double omega, omega1, omega2, omega3, omega4;
	double cone_axis[3], circle_center[3], cone_base_center[3];
	char message[MAXLINE];
	struct face *second_face;
	struct cycle *cyc, *cyc1, *cyc2;
	struct arc *circle_arc;
	struct edge *edg1, *edg2;
	struct circle *cone_circle;
	struct variety *vty;
	struct surface *this_srf;
    struct cept *ex;

	this_srf = given_face -> srf;
	if (this_srf == NULL) return (NULL);
	lfn = given_face -> lfn;
	ofn = given_face -> ofn;
	/* transfer info to local variables */
	cyc = given_face -> first_cycle;
	vty = given_face -> vty;
    atmnum0 = vty -> atmnum[0];
    atmnum1 = vty -> atmnum[1];
    atmnum2 = vty -> atmnum[2];
    atmnum3 = vty -> atmnum[3];
	alpha = given_face -> alpha;
	omega = compute_omega (vty -> radii[0], given_face);
	if (error()) return (0);
	if (omega < 0.0) {
		ex = new_cept (GEOMETRY_ERROR,  DEGENERACY,  FATAL_SEVERITY);
		add_function (ex, "circle_face");
		add_source (ex, "mssimple.c");
        add_double (ex, "omega", omega);
        add_long (ex, "atom number", atmnum0);
        add_double (ex, "current fineness angle", given_face -> alpha);
		add_remedy (ex, "Decrease fineness angle for atom");
		return (NULL);
	}

	/* get a vector pointing away from face */
	if (!away (given_face, cone_axis, &cone_angle)) return (NULL);
	if (error()) return (NULL);

	/* error checking */
	if (cone_angle <= 0.0 || cone_angle >= PI) {
		ex = new_cept (GEOMETRY_ERROR,  DEGENERACY,  FATAL_SEVERITY);
		add_function (ex, "circle_face");
		add_source (ex, "mssimple.c");
        add_double (ex, "cone_angle", cone_angle);
        add_long (ex, "first  atom number", atmnum0);
        add_long (ex, "second atom number", atmnum1);
        add_long (ex, "third  atom number", atmnum2);
		if (atmnum3 > 0)
			add_long (ex, "fourth atom number", atmnum3);
		return (0);
	}

	if (cone_angle > PI/3) cone_angle = PI/3;

	for (k = 0; k < 3; k++)
		cone_base_center[k] = vty -> center[k] +
			vty -> radii[0] * cone_axis[k];

	is_contained = point_in_face (cone_base_center, given_face, 0);
	if (error()) return (NULL);

	if (is_contained != 1) return (NULL);

	this_srf -> n_large_face++;

	/* compute center, radius of new circle */
	for (k = 0; k < 3; k++)
		circle_center[k] = vty -> center[k] +
			vty -> radii[0] * cos (cone_angle) * cone_axis[k];
	circle_radius = vty -> radii[0] * sin (cone_angle);
	/* allocate new circle and arc */
	cone_circle = new_circle (circle_center, circle_radius, cone_axis);
	if (error()) return (NULL);
	link_circle (this_srf, cone_circle);
	circle_arc = new_arc (cone_circle, (struct vertex *) NULL,
		(struct vertex *) NULL, CONVEX, 0, alpha, lfn, ofn, ofn);
	if (error()) return (NULL);
	link_arc (this_srf, circle_arc);

	/* need face pointer for newedg call */
	second_face = duplicate_face (given_face);
	if (error()) return (NULL);

	/* two new edges */
	edg1 = new_edge (circle_arc, 0, second_face, NULL);
	if (error()) return (NULL);
	edg2 = new_edge (circle_arc, 1, given_face, NULL);
	if (error()) return (NULL);

	/* two new cycles */
	cyc1 = new_cycle (NULL, edg1);
	if (error()) return (NULL);
	ink_cycle (this_srf);
	cyc2 = new_cycle (NULL, edg2);
	if (error()) return (NULL);
	ink_cycle (this_srf);

	/* new face has just one cycle */
	second_face -> first_cycle = cyc1;

	/* old face has 1 or 2 cycles */
	if (cyc == NULL) given_face -> first_cycle = cyc2;
	else cyc -> next = cyc2;

	omega1 = compute_omega (vty -> radii[0], given_face);
	if (error()) return(0);
	if (omega1 < 0.0) {
		ex = new_cept (GEOMETRY_ERROR,  DEGENERACY,  WARNING_SEVERITY);
		add_function (ex, "circle_face");
		add_source (ex, "mssimple.c");
        add_double (ex, "omega1", omega1);
        add_long (ex, "atom number", atmnum0);
        add_double (ex, "current fineness angle", given_face -> alpha);
		add_remedy (ex, "Decrease fineness angle for atom");
		simple_tent (second_face, 1);
		return (NULL);
	}
	omega2 = compute_omega (vty -> radii[0], second_face);
	if (error()) return(0);
	if (omega2 < 0.0) {
		ex = new_cept (GEOMETRY_ERROR,  DEGENERACY,  WARNING_SEVERITY);
		add_function (ex, "circle_face");
		add_source (ex, "mssimple.c");
        add_double (ex, "omega2", omega2);
        add_long (ex, "atom number", atmnum0);
        add_double (ex, "current fineness angle", given_face -> alpha);
		add_remedy (ex, "Decrease fineness angle for atom");
		simple_tent (second_face, 1);
		return (NULL);
	}

	if (!equal (omega1 + omega2, omega)) {
		ex = new_cept (GEOMETRY_ERROR,  INCONSISTENCY,  FATAL_SEVERITY);
		add_function (ex, "circle_face");
		add_source (ex, "mssimple.c");
        add_double (ex, "omega1", omega1);
        add_double (ex, "omega2", omega2);
        add_double (ex, "omega1+omega2", omega1+omega2);
        add_double (ex, "omega", omega);
        add_long (ex, "first  atom number", atmnum0);
        add_long (ex, "second atom number", atmnum1);
        add_long (ex, "third  atom number", atmnum2);
		if (atmnum3 > 0)
			add_long (ex, "fourth atom number", atmnum3);
		return (NULL);
	}


	/* don't create problems for the simplify function */
	if (circle_arc -> alpha < MIN_ALPHA)
		circle_arc -> alpha = MIN_ALPHA;

	n_point = subdivide_arc (this_srf, circle_arc);		/* subdivide new arc */
	if (n_point <= 0) return (0);

	omega3 = compute_omega (vty -> radii[0], given_face);
	if (error()) return(0);
	if (omega3 < 0.0) {
		sprintf (message, "(circle_face): Warning: decrease angle for atom %d to < %6.3f",
			(int) given_face -> vty -> atmnum[0], given_face -> alpha);
		inform(message);
		simple_tent (second_face, 1);
		return (NULL);
	}
	omega4 = compute_omega (vty -> radii[0], second_face);
	if (error()) return(0);
	if (omega4 < 0.0) {
		sprintf (message, "(circle_face): Warning: decrease angle for atom %d to < %6.3f",
			(int) given_face -> vty -> atmnum[0], given_face -> alpha);
		inform(message);
		simple_tent (second_face, 1);
		return (NULL);
	}

	if (!equal (omega1, omega3)) {
		ex = new_cept (GEOMETRY_ERROR,  INCONSISTENCY,  FATAL_SEVERITY);
		add_function (ex, "circle_face");
		add_source (ex, "mssimple.c");
		add_message (ex, "omega changes with arc subdivision");
        add_double (ex, "omega1", omega1);
        add_double (ex, "omega3", omega3);
        add_long (ex, "first  atom number", atmnum0);
        add_long (ex, "second atom number", atmnum1);
        add_long (ex, "third  atom number", atmnum2);
		if (atmnum3 > 0)
			add_long (ex, "fourth atom number", atmnum3);
 		return (0);
	}
	if (!equal (omega2, omega4)) {
		ex = new_cept (GEOMETRY_ERROR,  INCONSISTENCY,  FATAL_SEVERITY);
		add_function (ex, "circle_face");
		add_source (ex, "mssimple.c");
		add_message (ex, "omega changes with arc subdivision");
        add_double (ex, "omega2", omega2);
        add_double (ex, "omega4", omega4);
        add_long (ex, "first  atom number", atmnum0);
        add_long (ex, "second atom number", atmnum1);
        add_long (ex, "third  atom number", atmnum2);
		if (atmnum3 > 0)
			add_long (ex, "fourth atom number", atmnum3);
 		return (0);
	}

	return (second_face);	/* return pointer to new face */
}

/* return TRUE if numbers are approximately equal */
int equal (double double1, double double2)
{
	double diff;

	diff = fabs (double1 - double2);
	return (diff < EPSILON);
}

/* compute vector away from face,
and also angle of cone not intersecting face */

int away (struct face *spherical_face, double vector[3], double *cone_angle)
{
	int k, orn;
	double maximum_dot_product, trial_angle, this_dot_product;
	double alpha;
	double vertex_vector[3], cone_base_center[3];
	struct edge *edg;
	struct arc *a;
	struct cycle *cyc;
	struct vertex *vtx;
	struct variety *vty;
	struct surface *this_srf;
	
	this_srf = spherical_face -> srf;
	if (this_srf == NULL) return (0);
	vty = spherical_face -> vty;
	cyc = spherical_face -> first_cycle;
	alpha = spherical_face -> alpha;

	/* check for face = sphere */
	if (cyc == NULL) {
		/* arbitrary vector okay */
		for (k = 0; k < 3; k++)
			vector[k] = ((k == 2 ) ? 1.0 : 0.0);
		*cone_angle = PI / 3;	/* a nice angle */
		return (1);
	}

	if (cyc -> next != NULL) {
		set_error1 ("(away): > 1 cycle");
		return (0);
	}

	/* pick a point near the face center */
	if (!face_center (spherical_face, cone_base_center)) return (0);
	for (k = 0; k < 3; k++)
		vector[k] = cone_base_center[k] - vty -> center[k];

	normalize (vector);

	/* compute angle of cone not intersecting face */
	/* initialize to maximum dot product */
	maximum_dot_product = (- vty -> radii[0]);
	for (edg = cyc -> first_edge; edg != NULL; edg = edg -> next) {
		a = edg -> arcptr;
		orn = edg -> orn;
		vtx = a -> vtx[orn];
		if (vtx == NULL) {
			set_error1 ("(away): edge without vertices");
			return(0) ;
		}
		for (k = 0; k < 3; k++)
			vertex_vector[k] = (vtx -> center[k] - vty -> center[k]);
		this_dot_product = dot_product (vertex_vector, vector);
		if (this_dot_product > maximum_dot_product)
			maximum_dot_product = this_dot_product;
	}

	/* compute cone angle from maximum dot product */
	maximum_dot_product /= vty -> radii[0];
	if (maximum_dot_product < -1.0) maximum_dot_product = -1.0;
	else if (maximum_dot_product > 1.0) maximum_dot_product = 1.0;
	trial_angle = acos (maximum_dot_product);
	if (trial_angle < 0.15) return (0);

	/* move away from cycle a bit */
	if (trial_angle >  alpha + MIN_ALPHA) trial_angle -= alpha/2;
	else trial_angle /= 2.0;

	if (trial_angle < 0.05) return (0);

	/* store angle in calling function */
	*cone_angle = trial_angle;
	return (1);
}

int face_center (struct face *given_face, double returned_point[3])
{
	int i, j, k, is_contained, point_okay, nlat, nlong;
	double sphere_radius, best_vertex_distance, this_distance;
	double phi, theta, alpha, sin_phi, cos_phi, x, y, z;
	double test_point[3], sphere_center[3], best_point[3];

	if (given_face -> shape != CONVEX) {
		set_error1 ("(face_center): face not convex");
		return (0);
	}
	alpha = given_face -> alpha;
	if (alpha <= MIN_ALPHA) alpha = MIN_ALPHA;
	if (alpha > 1.0) alpha = 1.0;

	sphere_radius = given_face -> vty -> radii[0];
	for (k = 0; k < 3; k++)
		sphere_center[k] = given_face -> vty -> center[k];

	nlat = N_LATITUDE / alpha + 0.5;
	nlong = N_LONGITUDE / alpha + 0.5;

	best_vertex_distance = 0.0;
	point_okay = 0;
	for (i = 0; i < nlat; i++) {
		phi = PI * (i + 0.5) / nlat;
		cos_phi = cos (phi);
		sin_phi = sin (phi);
		for (j = 0; j < nlong; j++) {
			theta = 2 * PI * (j + 0.5) / nlong;
			x = sphere_radius * sin_phi * cos (theta);
			y = sphere_radius * sin_phi * sin (theta);
			z = sphere_radius * cos_phi;
			test_point[0] = sphere_center[0] + x;
			test_point[1] = sphere_center[1] + y;
			test_point[2] = sphere_center[2] + z;
			this_distance = closest_vertex (test_point, given_face);
			if (this_distance >= best_vertex_distance) {
				is_contained = point_in_face (test_point, given_face, 0);
				if (error()) return (0);
				if (is_contained == 1) {
					point_okay = 1;
					best_vertex_distance = this_distance;
					for (k = 0; k < 3; k++)
						best_point[k] = test_point[k];
				}
			}
		}
	}
	if (!point_okay) return (0);

	for (k = 0; k < 3; k++)
		returned_point[k] = best_point[k];
	return (1);
}

double closest_vertex (double pnt[3], struct face *fac)
{
	int k, orn;
	double maximum_dot_product, trial_angle, this_dot_product;
	double vertex_vector[3], vector[3];
	struct edge *edg;
	struct arc *a;
	struct cycle *cyc;
	struct vertex *vtx;
	struct variety *vty;

	vty = fac -> vty;
	if (vty -> type != SPHERE) {
		set_error1 ("(closest_vertex): variety not sphere");
		return (0.0);
	}
	cyc = fac -> first_cycle;
	for (k = 0; k < 3; k++)
		vector[k] = pnt[k] - vty -> center[k];
	normalize (vector);

	/* compute angle of cone not intersecting face */
	/* initialize to maximum dot product */
	maximum_dot_product = (- vty -> radii[0]);
	for (edg = cyc -> first_edge; edg != NULL; edg = edg -> next) {
		a = edg -> arcptr;
		orn = edg -> orn;
		vtx = a -> vtx[orn];
		if (vtx == NULL) {
			set_error1 ("(closest_vertex): edge without vertices");
			return (0.0);
		}
		for (k = 0; k < 3; k++)
			vertex_vector[k] = (vtx -> center[k] - vty -> center[k]);
			this_dot_product = dot_product (vertex_vector, vector);
			if (this_dot_product > maximum_dot_product)
				maximum_dot_product = this_dot_product;
	}

	/* compute cone angle from maximum dot product */
	maximum_dot_product /= vty -> radii[0];
	if (maximum_dot_product < -1.0) maximum_dot_product = -1.0;
	else if (maximum_dot_product > 1.0) maximum_dot_product = 1.0;
	trial_angle = acos (maximum_dot_product);
	if (trial_angle <= 0.0) trial_angle = 0.0;

	return (trial_angle * vty -> radii[0]);
}


int is_saddle_cone (struct face *given_face)
{
	int n_concave, n_straight, n_different, one_circle;
	struct edge *e, *e1, *e2;
	struct cycle *cyc;

	if (given_face -> shape != SADDLE) return (0);
	cyc = given_face -> first_cycle;
	if (cyc == NULL) return (0);
	if (cyc -> next != NULL) return (0);

	n_concave = 0;
	/* count number concave edges */
	for (e = cyc -> first_edge; e != NULL; e = e -> next)
		if (e -> arcptr -> shape == CONCAVE) n_concave++;
	if (n_concave != 0) return (0);

	n_straight = 0;
	/* count number straight edges */
	for (e = cyc -> first_edge; e != NULL; e = e -> next)
		if (e -> arcptr -> shape == CONCAVE) n_straight++;
	if (n_straight != 0) return (0);

	/* all edges are convex */

	n_different = 0;
	for (e1 = cyc -> first_edge; e1 != NULL; e1 = e1 -> next)
		for (e2 = cyc -> first_edge; e2 != NULL; e2 = e2 -> next)
			if (e1 -> arcptr -> cir != e2 -> arcptr -> cir)
				n_different++;
	one_circle = (n_different == 0);
	return (one_circle);
}

struct vertex *find_cone_vertex (struct face *given_face)
{
	double d_to_vtx, best;
	double center[3];
	struct vertex *vtx, *best_vtx;
	struct cycle *cyc;
	struct surface *this_srf;

	this_srf = given_face -> srf;
	if (this_srf == NULL) return (NULL);

	if (given_face -> shape != SADDLE)  {
		set_error1 ("(find_cone_vertex): face not saddle");
		return (0);
	}
	cyc = given_face -> first_cycle;
	if (cyc == NULL)  {
		set_error1 ("(find_cone_vertex): face has no boundary");
		return (0);
	}
	if (!cyc_centroid (cyc, center)) return (NULL);
	if (error()) return (0);

	/* find closest point cusp vertex */
	best_vtx = NULL;
	best = 1000000.0;
	for (vtx = this_srf -> head_vertex; vtx != NULL; vtx = vtx -> next) {
		if (!vtx -> cusp) continue;
		d_to_vtx = distance (center, vtx -> center);
		if (d_to_vtx < best) {
			best = d_to_vtx;
			best_vtx = vtx;
		}
	}
	return (best_vtx);
}


int do_saddle_cone (struct face *given_face)
{
	char message[MAXLINE];
	struct vertex *best_vtx;
	struct cycle *cyc;
	struct surface *this_srf;

	this_srf = given_face -> srf;
	if (this_srf == NULL) return (0);
	if (given_face -> shape != SADDLE)  {
		set_error1 ("(do_saddle_cone): face not saddle");
		return (0);
	}
	cyc = given_face -> first_cycle;
	if (cyc == NULL)  {
		set_error1 ("(do_saddle_cone): face has no boundary");
		return (0);
	}
	best_vtx = find_cone_vertex (given_face);
	if (best_vtx == NULL) {
		set_error1 ("(do_saddle_cone): vertex not found");
		return (0);
	}
	sprintf (message, "%8ld vertex number for saddle cone", best_vtx -> number);
	informd (message);
	if (!tent (given_face, best_vtx)) return (0);
	if (error()) return(0);
	return (1);
}


/*
   MSRoll
   Copyright 1996 by Michael L. Connolly
   All Rights Reserved

*/
