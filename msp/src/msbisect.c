/*

	MSRoll

	Copyright 1986, 1989, 1996 by Michael L. Connolly
	All rights reserved

	Written by Michael L. Connolly.
	November 27, 2001

*/

#include "ms.h"
#include "msenum.h"
#include "mslimit.h"
#include "msdefault.h"
#include "msstruct.h"
#include "msglobal.h"
#include "msfunction.h"


/* FACE SUBDIVISION */

/* bisect face recursively */
int bisect_face (struct face *first_face)
{
	int n_point;
	int n_edge_per_face, result, ne1, ne2;
	double omega, circumference;
	char message[MAXLINE];
	struct arc *bisecting_arc;
	struct face *second_face;
	struct edge *e1, *e2, *e3, *e4;
	struct surface *this_srf;
	struct variety *vty;

	/* check argument validity */
	if (first_face -> problem) {
		inform ("problem face skipped in bisection");
		return(1);
	}
	this_srf = first_face -> srf;
	if (this_srf == NULL) return(0);
	n_edge_per_face = 0;
	vty = first_face -> vty;
	if (first_face -> first_cycle == NULL) {
		/* sphere */
		second_face = entire_sphere (first_face);
		if (error()) return(0);
		result = bisect_face (first_face);
		if (!result) return (0);
		if (error()) return(0);
		result = bisect_face (second_face);
		if (!result) return (0);
		if (error()) return(0);
		return(1);
	}
	if (is_saddle_cone (first_face)) {
		do_saddle_cone (first_face);
		return(1);
	}
	/* make face simply connected */
	if (!simplify (first_face)) {	/* routine failed */
		informd("simplify failed, call multiple tent routine");
		multiple_tent (first_face);
		return(1);
	}
	/* count edges in face and compute face circumference */
	if (first_face -> first_cycle == NULL) {
		n_edge_per_face = 0;
		circumference = 0.0;
	}
	else {
		n_edge_per_face = edges_in_cycle (first_face -> first_cycle);
		circumference = circum (first_face -> first_cycle);
	}
	if (error()) return (0);
	n_edge_per_face = edges_in_cycle (first_face -> first_cycle);
	if (n_edge_per_face == 0)
		return(0);
	if (n_edge_per_face <= 3) {		/* already a triangle */
		if (!convert_face (this_srf, first_face)) return(0);
		return(1);
	}
	n_edge_per_face = edges_in_cycle (first_face -> first_cycle);
	if (first_face -> shape == CONVEX) {
		n_edge_per_face = edges_in_cycle (first_face -> first_cycle);
		omega = compute_omega (vty -> radii[0], first_face);
		if (error()) return(0);
		if (omega < 0.0) {
			sprintf (message, "(bisect_face): omega = %8.3f; call simple tent (1)", omega);
			informd (message);
			simple_tent (first_face, 1);
			if (error()) return(0);
			return(1);
		}
	}
	else omega = 0.0;
	n_edge_per_face = edges_in_cycle (first_face -> first_cycle);

	/* later: add depth integer to face structure */

	if (first_face -> largeok &&
		first_face -> shape == CONVEX &&
		omega > this_srf -> large_omega &&
		circumference < 4 * PI * first_face -> vty -> radii[0]) {
		n_edge_per_face = edges_in_cycle (first_face -> first_cycle);
		/* much or most of sphere but not all */
		second_face = circle_face (first_face);
		result = bisect_face (first_face);
		if (!result) return (0);
		if (error()) return(0);
		result = bisect_face (second_face);
		if (!result) return (0);
		if (error()) return(0);
		return(1);
	}
	n_edge_per_face = edges_in_cycle (first_face -> first_cycle);
	bisecting_arc = select_pair (first_face, &e1, &e2, &e3, &e4);
	if (error()) return(0);

	if (bisecting_arc == NULL) {	/* routine failed */
		simple_tent (first_face, 0);
		return(1);
	}

	/* split face into two faces */
	second_face = split_face (first_face, bisecting_arc, e1, e2, e3, e4);
	if (second_face == NULL) return(1);	/* unable to split, did tent */

	n_point = subdivide_arc (this_srf, bisecting_arc);			/* subdivide new arc */
	if (error()) return(0);
	if (n_point <= 0) return (0);

	ne1 = edges_in_cycle (first_face -> first_cycle);
	ne2 = edges_in_cycle (second_face -> first_cycle);
	result = bisect_face (first_face);			/* bisect first subface */
	if (error()) return(0);
	if (!result) return (0);
	result = bisect_face (second_face);		/* bisect second subface */
	if (error()) return(0);
	if (!result) return (0);
	return(1);
}

/* split face into two */

struct face *split_face (struct face *given_face, struct arc *bisecting_arc, struct edge *e1, struct edge *e2, struct edge *e3, struct edge *e4)
{
	int original_n_edge, cyc1_n_edge, cyc2_n_edge;
    int later_n_edge;
	double omega, omega1, omega2;
	char message[MAXLINE];
	struct edge *e5, *e6;
	struct edge *edg, *ending_edge;
	struct cycle *cyc1, *cyc2;
	struct face *second_face;
	struct variety *vty;
	struct surface *this_srf;
    struct cept *ex;

	this_srf = given_face -> srf;
	if (this_srf == NULL) return (NULL);
	/* pointer to cycle */
	cyc1 = given_face -> first_cycle;
	vty = given_face -> vty;

	original_n_edge = edges_in_cycle (cyc1);
	if (given_face -> shape == CONVEX) {
		omega = compute_omega (vty -> radii[0], given_face);
		if (error()) return(0);
		if (omega < 0.0) {
			/* face not correct */
			sprintf (message, "(split_face): omega < 0.0; call simple tent (3) %8.3f", omega);
			informd (message);
			simple_tent (given_face, 1);
			return (NULL);
		}
	}
	else omega = 0.0;

	/* two new edges for bisecting arc */
	e5 = new_edge (bisecting_arc, 0, given_face, NULL);
	if (error()) return(NULL);
	e6 = new_edge (bisecting_arc, 1, given_face, NULL);
	if (error()) return(NULL);

	/* find last arc of cycle */
	for (edg = cyc1 -> first_edge; edg != NULL; edg = edg -> next)
		if (edg -> next == NULL) ending_edge = edg;

	/* cyclize */
	ending_edge -> next = cyc1 -> first_edge;

	/* change links so there will be two cycles */

	e1 -> next = NULL;
	e5 -> next = e4;
	e3 -> next = NULL;
	e6 -> next = e2;

	/* change beginning edge of first cycle */
	cyc1 -> first_edge = e5;

	cyc1_n_edge = edges_in_cycle (cyc1);

	/* allocate new cycle */

	cyc2 = new_cycle (NULL, e6);
	ink_cycle (this_srf);

	cyc2_n_edge = edges_in_cycle (cyc2);

	later_n_edge = cyc1_n_edge + cyc2_n_edge - 2;
	if (original_n_edge != later_n_edge) {
		ex = new_cept (LOGIC_ERROR,  INCONSISTENCY,  FATAL_SEVERITY);
		add_function (ex, "split_face");
		add_source (ex, "msbisect.c");
		add_long (ex, "original_n_edge", original_n_edge);
		add_long (ex, "later_n_edge", later_n_edge);
		return (NULL);
	}
	if (cyc1_n_edge <= 2) {
		ex = new_cept (LOGIC_ERROR,  INVALID_VALUE,  FATAL_SEVERITY);
		add_function (ex, "split_face");
		add_source (ex, "msbisect.c");
		add_long (ex, "cyc1_n_edge", cyc1_n_edge);
		return (NULL);
	}
	if (cyc2_n_edge <= 2) {
		ex = new_cept (LOGIC_ERROR,  INVALID_VALUE,  FATAL_SEVERITY);
		add_function (ex, "split_face");
		add_source (ex, "msbisect.c");
		add_long (ex, "cyc2_n_edge", cyc2_n_edge);
		return (NULL);
	}

	/* duplicate face */
	second_face = duplicate_face (given_face);	
	if (error()) return(NULL);
	second_face -> first_cycle = cyc2;

	/* pointers from new edges to faces */
	e5 -> fac = given_face;
	e6 -> fac = second_face;

	/* check whether area split */

	if (given_face -> shape == CONVEX && original_n_edge > 6) {
		omega1 = compute_omega (vty -> radii[0], given_face);
		if (error()) return(0);
		if (omega1 < 0.0) {
			informd ("(split_face): omega < 0.0; call simple tent (4)");
			simple_tent (given_face, 1);
			if (error()) return(NULL);
			informd ("(split_face): omega < 0.0; call simple tent (5)");
			simple_tent (second_face, 1);
			return (NULL);
		}
		omega2 = compute_omega (vty -> radii[0], second_face);
		if (error()) return(0);
		if (omega2 < 0.0) {
			informd ("(split_face): omega < 0.0; call simple tent (6)");
			simple_tent (given_face, 1);
			if (error()) return(NULL);
			informd ("(split_face): omega < 0.0; call simple tent (7)");
			simple_tent (second_face, 1);
			return (NULL);
		}
		if (!equal (omega1+omega2, omega)) {
			/* area not split correctly */
			informd ("(split_face): omega not equal; call simple tent (8)");
			simple_tent (given_face, 1);
			if (error()) return(NULL);
			informd ("(split_face): omega not equal; call simple tent (9)");
			simple_tent (second_face, 1);
			return (NULL);
		}
	}
	return (second_face);		/* return pointer to 2nd face */
}

/* select a pair of vertices to connect for face bisection */
/* this function is a mess */

struct arc *select_pair (struct face *given_face, struct edge **e1, struct edge **e2, struct edge **e3, struct edge **e4)
{
	int n_edge_f, n_concave, n_straight;
	int n_edge_per_face;
	int repeat_found;
	int n_contain_fail, n_delta_fail, n_pair;
	int one_three, orn, i, j, k, ibest, jbest;
	int orn1, orn2, orn3, orn4;
	int is_contained, weird;
	int i_before, j_before, small, shape, sign, around_axis_first;
	int one_circle, n_different, inserted;
	long lfn, ofn;
	double alpha, around, vtx_vtx_dist, doti, dotj;
	double ij_distance, ij_value, ji_value, this_value;
	double i_perim, j_perim, circle_radius;
	double prev_forward, prev_backward, difference;
	double delta1, delta2, delta3, delta4;
	double i_vector[3], j_vector[3], ij_vector[3], mid_point[3];
	double circle_center[3], circle_axis[3];
	double vtx1_vector[3], vtx2_vector[3];
	char message[MAXLINE];
	struct circle temporary_circle;
	struct variety *vty;
	struct edge *e;
	struct edge *edg1, *edg2, *edg3, *edg4;
	struct edge temporary_edge1, temporary_edge2;
	struct arc *arc_new, *arc_ptr;
	struct circle *cir, *prev_cir;
	struct vertex *v, *vtx1, *vtx2, *vtx3, *vtx4;
	struct vertex *first_vtx, *second_vtx;
	struct cycle *cyc;
	struct arc temporary_arc;
	struct surface *this_srf;
	double *perimeter = NULL;
	double *arc_direction = NULL;
	double *median = NULL;
	double *distance_forward = NULL;
	double *distance_backward = NULL;
	struct vertex **vertex_list = NULL;
	struct arc **arc_list = NULL;
	struct edge **edge_list = NULL;
	long *bad_vertex = NULL;
	long *bad_neighbor = NULL;
    struct cept *ex;

	/* the linked list of pairs is in order of decreasing goodness */
	struct vertex_pair *head_pair = NULL;
	struct vertex_pair *tail_pair = NULL;
	struct vertex_pair *new_pair = NULL;
	struct vertex_pair *pair_ptr = NULL;
	struct vertex_pair *prev_pair_ptr = NULL;
	struct vertex_pair *next_pair_ptr = NULL;

	/* store pointers to variety and (only) cycle */
	vty = given_face -> vty;
	cyc = given_face -> first_cycle;
	alpha = given_face -> alpha;
	this_srf = given_face -> srf;
	if (this_srf == NULL) return (NULL);

	/* error checking */
	if (cyc == NULL) {
		ex = new_cept (PARAMETER_ERROR,  NULL_VALUE,  FATAL_SEVERITY);
        add_object (ex, CYCLE, "cyc");
		add_function (ex, "select_pair");
		add_source (ex, "msbisect.c");
		return (NULL);
	}

	n_edge_per_face = edges_in_cycle (given_face -> first_cycle);
	if (n_edge_per_face == 0) {
		ex = new_cept (GEOMETRY_ERROR,  NULL_VALUE,  FATAL_SEVERITY);
		add_function (ex, "select_pair");
		add_source (ex, "msbisect.c");
		add_message (ex, "cycle has no edges");
		return (NULL);
	}
	if (cyc -> next != NULL) {
		ex = new_cept (GEOMETRY_ERROR,  NULL_VALUE,  FATAL_SEVERITY);
		add_function (ex, "select_pair");
		add_source (ex, "msbisect.c");
		add_message (ex, "cycle is not simply connected");
		return (NULL);
	}

	/* initialization */
	head_pair = NULL;
	tail_pair = NULL;
	n_edge_f = 0;
	n_concave = 0;
	n_straight = 0;
	lfn = given_face -> lfn;
	ofn = given_face -> ofn;

	/* count */
	for (e = cyc -> first_edge; e != NULL; e = e -> next) {
		n_edge_f++;
		if (e -> arcptr -> shape == CONCAVE) n_concave++;
		if (e -> arcptr -> shape == STRAIGHT) n_straight++;
	}
	around_axis_first =
		((given_face -> shape == SADDLE && n_concave > 2) ||
		(given_face -> shape == CYLINDRICAL && n_straight > 2));
	if (n_edge_f <= 3) {
		ex = new_cept (LOGIC_ERROR,  INCONSISTENCY,  FATAL_SEVERITY);
		add_function (ex, "select_pair");
		add_source (ex, "msbisect.c");
		add_long (ex, "edges in face", n_edge_f);
		add_message (ex, "select_pair was mistakenly called for triangular face");
		return (NULL);	/* don't bisect triangle */
	}

	if (!around_axis_first && n_edge_f == 4 &&
		(given_face -> shape != CONVEX)) {
		edg1 = cyc -> first_edge;
		edg2 = edg1 -> next;
		edg3 = edg2 -> next;
		edg4 = edg3 -> next;
		orn1 = edg1 -> orn;
		orn2 = edg2 -> orn;
		orn3 = edg3 -> orn;
		orn4 = edg4 -> orn;
		vtx1 = edg1 -> arcptr -> vtx[orn1];
		vtx2 = edg2 -> arcptr -> vtx[orn2];
		vtx3 = edg3 -> arcptr -> vtx[orn3];
		vtx4 = edg4 -> arcptr -> vtx[orn4];
		one_three = first_and_third (vtx1 -> center, vtx2 -> center,
			vtx3 -> center, vtx4 -> center);
		if (error()) return(NULL);
		if (one_three) {
			first_vtx = vtx1;
			second_vtx = vtx3;
			/* store edge pointers in calling routine */
			*e1 = edg4;
			*e2 = edg1;
			*e3 = edg2;
			*e4 = edg3;
		}
		else {
			first_vtx = vtx2;
			second_vtx = vtx4;
			/* store edge pointers in calling routine */
			*e1 = edg1;
			*e2 = edg2;
			*e3 = edg3;
			*e4 = edg4;
		}
		/* short and straight */
		arc_new =  make_straight_arc (this_srf, first_vtx, second_vtx);
		if (error()) return(NULL);
		return (arc_new);
	}

	/* allocate arrays */

	/* pointers to vertices */
	vertex_list = (struct vertex **) 
		allocate_pointers (VERTEX, n_edge_f);
	if (vertex_list == NULL) {
		print_counts();
		ex = new_cept (MEMORY_ERROR,  ALLOCATION,  FATAL_SEVERITY);
        add_object (ex, VERTEX, "vertex_list");
        add_function (ex, "select_pair");
		add_source (ex, "msbisect.c");
		add_long (ex, "n_edge_f", n_edge_f);
		return (NULL);
	}

	/* bad vertices (on 4 edges ) */
	bad_vertex = allocate_longs (n_edge_f, 0, BAD_VERTEX);
	if (bad_vertex == NULL) {
		print_counts();
		ex = new_cept (MEMORY_ERROR,  ALLOCATION,  FATAL_SEVERITY);
        add_object (ex, VERTEX, "bad_vertex");
        add_function (ex, "select_pair");
		add_source (ex, "msbisect.c");
		add_long (ex, "n_edge_f", n_edge_f);
		return (NULL);
	}

	/* pointers to arcs */
	arc_list = (struct arc **)
		allocate_pointers (ARC, n_edge_f);
	if (arc_list == NULL) {
		print_counts();
		ex = new_cept (MEMORY_ERROR,  ALLOCATION,  FATAL_SEVERITY);
        add_object (ex, ARC, "arc_list");
        add_function (ex, "select_pair");
		add_source (ex, "msbisect.c");
		add_long (ex, "n_edge_f", n_edge_f);
		return (NULL);
	}

	/* pointers to edges */
	edge_list = (struct edge **)
		allocate_pointers (EDGE, n_edge_f);
	if (edge_list == NULL) {
		print_counts();
		ex = new_cept (MEMORY_ERROR,  ALLOCATION,  FATAL_SEVERITY);
        add_object (ex, EDGE, "edge_list");
        add_function (ex, "select_pair");
		add_source (ex, "msbisect.c");
		add_long (ex, "n_edge_f", n_edge_f);
		return (NULL);
	}

	/* distance around perimeter */
	perimeter = allocate_doubles (n_edge_f, 0, PERIMETER);
	if (perimeter == NULL) {
		print_counts();
		ex = new_cept (MEMORY_ERROR,  ALLOCATION,  FATAL_SEVERITY);
        add_variable (ex, PERIMETER, "perimeter");
        add_function (ex, "select_pair");
		add_source (ex, "msbisect.c");
		add_long (ex, "n_edge_f", n_edge_f);
		return (NULL);
	}

	/* direction of arc relative to torus or cylinder axis */
	arc_direction = allocate_doubles (n_edge_f, 0, ARC_DIRECTION);
	if (arc_direction == NULL) {
		print_counts();
		ex = new_cept (MEMORY_ERROR,  ALLOCATION,  FATAL_SEVERITY);
        add_variable (ex, ARC_DIRECTION, "arc_direction");
        add_function (ex, "select_pair");
		add_source (ex, "msbisect.c");
		add_long (ex, "n_edge_f", n_edge_f);
		return (NULL);
	}

	distance_forward = allocate_doubles (n_edge_f, 0, DISTANCE_FORWARD);
	if (distance_forward == NULL) {
		print_counts();
		ex = new_cept (MEMORY_ERROR,  ALLOCATION,  FATAL_SEVERITY);
        add_variable (ex, DISTANCE_FORWARD, "distance_forward");
        add_function (ex, "select_pair");
		add_source (ex, "msbisect.c");
		add_long (ex, "n_edge_f", n_edge_f);
		return (NULL);
	}

	distance_backward = allocate_doubles (n_edge_f, 0, DISTANCE_BACKWARD);
	if (distance_backward == NULL) {
		print_counts();
		ex = new_cept (MEMORY_ERROR,  ALLOCATION,  FATAL_SEVERITY);
        add_variable (ex, DISTANCE_BACKWARD, "distance_backward");
        add_function (ex, "select_pair");
		add_source (ex, "msbisect.c");
		add_long (ex, "n_edge_f", n_edge_f);
		return (NULL);
	}

	/* distance around perimeter */
	median = allocate_doubles (n_edge_f, 0, MEDIAN);
	if (median == NULL) {
		print_counts();
		ex = new_cept (MEMORY_ERROR,  ALLOCATION,  FATAL_SEVERITY);
        add_variable (ex, MEDIAN, "median");
        add_function (ex, "select_pair");
		add_source (ex, "msbisect.c");
		add_long (ex, "n_edge_f", n_edge_f);
		return (NULL);
	}

	/* set up arrays */

	for (i = 0, e = cyc -> first_edge; e != NULL; e = e -> next, i++) {
		arc_ptr = e -> arcptr;
		orn = e -> orn;
		v = arc_ptr -> vtx[orn];
		if (v == NULL) {
			ex = new_cept (POINTER_ERROR,  NULL_VALUE,  FATAL_SEVERITY);
			add_object (ex, VERTEX, "v");
			add_function (ex, "select_pair");
			add_source (ex, "msbisect.c");
			add_message (ex, "no vertex for edge");
			return (NULL);
		}
		*(vertex_list + i) = v;	/* store vertex */
		*(bad_vertex + i) = 0;	/* initialize to okay */
		*(arc_list + i) = arc_ptr;	/* store arc */
		*(edge_list + i) = e;		/* store edge */
	}

	/* check for one circle (or almost one circle ) */
	n_different = 0;
	for (i = 0; i < n_edge_f - 1; i++)
		for (j = i + 1; j < n_edge_f; j++) {
			if ((*(arc_list+i)) -> cir != (*(arc_list+j)) -> cir)
				n_different++;
		}
	repeat_found = 0;
	for (i = 0; i < n_edge_f - 1; i++) {
		for (j = i + 1; j < n_edge_f; j++) {
			if (*(vertex_list+i)  == *(vertex_list+j)) {
				repeat_found = 1;
				break;
			}
		}
		if (repeat_found) break;
	}
	if (repeat_found && !given_face -> simplified) {
		print_face (given_face);
		ex = new_cept (LOGIC_ERROR,  INVALID_VALUE,  FATAL_SEVERITY);
		add_function (ex, "select_pair");
		add_source (ex, "msbisect.c");
		add_message (ex, "repeated vertex");
		return (NULL);
	}
	else if (!repeat_found && given_face -> simplified) {
		given_face -> simplified = 0;
	}
	one_circle = (n_different < n_edge_f);

	around = 0.0;		/* initialize */

	/* distances around perimeter */
	for (i = 0; i < n_edge_f; i++) {
		*(perimeter + i) = around;
		arc_ptr = *(arc_list + i);
		j = ((i < n_edge_f - 1) ? i + 1 : 0);
		vtx1 = *(vertex_list + i);
		vtx2 = *(vertex_list + j);
		if (arc_ptr -> shape == STRAIGHT) /* through space */
			vtx_vtx_dist = distance (vtx1 -> center, vtx2 -> center);
		else /* along arc */
			vtx_vtx_dist = arc_ptr -> cir -> radius *
				arc_ang (arc_ptr);
		around += vtx_vtx_dist;	/* increment cummulative distance */
		for (k = 0; k < 3; k++)
			ij_vector[k] = vtx2 -> center[k] - vtx1 -> center[k];
		*(arc_direction+i) = dot_product (ij_vector, vty -> axis);
	}

	/* set up middle-of-arc-ness arrays */
	prev_cir = NULL;
	for (i = 0; i < n_edge_f; i++) {
		arc_ptr = *(arc_list + i);
		if (arc_ptr -> cir != prev_cir)
			prev_forward = (*(perimeter + i));
		*(distance_forward + i) = (*(perimeter + i)) - prev_forward;
		prev_cir = arc_ptr -> cir;
	}

	prev_cir = NULL;
	for (i = n_edge_f - 1; i >= 0; i--) {
		arc_ptr = *(arc_list+i);
		if (arc_ptr -> cir != prev_cir)
			prev_backward = ((i == n_edge_f - 1)
				? around : (*(perimeter + i + 1)));
		*(distance_backward + i) = prev_backward - (*(perimeter + i ));
		prev_cir = arc_ptr -> cir;
	}

	for (i = 0; i < n_edge_f; i++) {
		arc_ptr = *(arc_list + i);
		difference = (*(distance_forward + i ))
			- (*(distance_backward + i));
		*(median + i) = fabs (difference);
	}


	/* intialization */
	weird = 0;

	/* check for bad vertices (if any, face is weird) */

	for (i = 0; i < n_edge_f - 1; i++)
		for (j = i + 1; j < n_edge_f; j++) {
			vtx1 = *(vertex_list + i);
			vtx2 = *(vertex_list + j);
			if (vtx1 == vtx2 && !around_axis_first) {
				weird = 1;
				*(bad_vertex + i) = 1;
				*(bad_vertex + j) = 1;
			}
		}

	/* forbid neighborhood of bad vertices */
	bad_neighbor = allocate_longs (n_edge_f, 0, BAD_NEIGHBOR);
	if (bad_neighbor == NULL) {
		print_counts();
		ex = new_cept (MEMORY_ERROR,  ALLOCATION,  FATAL_SEVERITY);
        add_variable (ex, BAD_NEIGHBOR, "bad_neighbor");
        add_function (ex, "select_pair");
		add_source (ex, "msbisect.c");
		add_long (ex, "n_edge_f", n_edge_f);
		return (NULL);
	}

	if (given_face -> shape != CYLINDRICAL) { 
		if (weird)
			for (i = 0; i < n_edge_f; i++) {
				if (! *(bad_vertex + i)) continue;
				j = i - 1;
				if (j < 0) j = n_edge_f - 1;
				*(bad_neighbor + j) = 1;
				j = i + 1;
				if (j > n_edge_f - 1) j = 0;
				*(bad_neighbor + j) = 1;
			}

		if (n_edge_f > 10)
			for (i = 0; i < n_edge_f; i++)
				if (*(bad_neighbor + i))
					*(bad_vertex + i) = 1;
	}

	free_longs (bad_neighbor, 0, BAD_NEIGHBOR);

	/* check each pair of vertices for possible connecting arc */

	for (i = 0; i < n_edge_f - 1; i++) {
		if (*(bad_vertex + i)) continue;
		i_before = ((i == 0) ? n_edge_f - 1 : i - 1);
		/* connect concave arcs first */
		if (around_axis_first) {
			if ((*(arc_list+i)) -> shape == CONVEX) continue;
			if ((*(arc_list+i_before)) -> shape == CONVEX) continue;
		}

		for (j = i + 1; j < n_edge_f; j++) {
			shape = 0; /* init */
			if (*(bad_vertex + j)) continue;
			j_before = ((j == 0) ? n_edge_f - 1 : j - 1);
			if (around_axis_first) {
				if ((*(arc_list+j)) -> shape == CONVEX) continue;
				if ((*(arc_list+j_before)) -> shape == CONVEX) continue;
			}
			/* store vertex pointers in local variables */
			vtx1 = *(vertex_list + i);
			vtx2 = *(vertex_list + j);

			if ( ! (around_axis_first && vtx1 == vtx2)) {
				/* check for not on same arc */
				if (*(arc_list + i) == *(arc_list + j)) continue;
				if (*(arc_list + i_before) == *(arc_list + j)) continue;
				if (*(arc_list + i) == *(arc_list + j_before)) continue;
				if (*(arc_list + i_before) == *(arc_list + j_before))
					continue;
				if (!one_circle) {
					/* check for not on same circle */
					if ((*(arc_list + i)) -> cir ==
						(*(arc_list + j)) -> cir) continue;
					if ((*(arc_list + i_before)) -> cir ==
						(*(arc_list + j)) -> cir) continue;
					if ((*(arc_list + i)) -> cir ==
						(*(arc_list + j_before)) -> cir) continue;
					if ((*(arc_list + i_before)) -> cir ==
						(*(arc_list + j_before)) -> cir) continue;
				}
			}

			if (vtx1 == vtx2 && ! around_axis_first) continue;

			for (k = 0; k < 3; k++) {
				i_vector[k] = vtx1 -> center[k] - vty -> center[k];
				j_vector[k] = vtx2 -> center[k] - vty -> center[k];
			}

			/* connect concave arcs of saddle face by convex arc */
			/* or connect straight arcs of cylinder by convex arc */
			if (around_axis_first) {
				/* project vertices onto torus axis */
				doti = dot_product (i_vector, vty -> axis);
				dotj = dot_product (j_vector, vty -> axis);
				/* check for good opposition */
				if (fabs (doti - dotj) > EPSILON) continue;
				shape = CONVEX;
			}

			/* basic idea: minimize connecting arc length and
			   new face perimeter inequality */

			if (given_face -> shape != CONVEX)
				ij_distance = distance (vtx1 -> center, vtx2 -> center);
			else {
				cross (i_vector, j_vector, circle_axis);
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
				ij_distance = temporary_circle.radius *
					arc_ang (&temporary_arc);
			}

			/* choosing among possibilities by scoring */

			i_perim = *(perimeter + i);
			j_perim = *(perimeter + j);
			ij_value = j_perim - i_perim;
			ji_value = (i_perim + around) - j_perim;
			this_value = (double) fabs (ij_value - ji_value) / 2;
			if (given_face -> shape == SADDLE && n_concave <= 2)
				this_value += ((*(median + i)) + (*(median+j))) / 4;
			if (given_face -> shape == CYLINDRICAL)
				this_value = ij_distance;
			else
				this_value =
					(1-this_srf->weight) * this_value + this_srf -> weight * ij_distance;

			/* allocate pair structure */
			new_pair = (struct vertex_pair *) allocate_object (VERTEX_PAIR);
			if (new_pair == NULL) {
				print_counts();
				ex = new_cept (MEMORY_ERROR,  ALLOCATION,  FATAL_SEVERITY);
				add_object (ex, VERTEX_PAIR, "new_pair");
				add_function (ex, "select_pair");
				add_source (ex, "msbisect.c");
				return (NULL);
			}

			/* link into list */
			if (head_pair == NULL) {
				/* start linked list */
				head_pair = new_pair;
				tail_pair = new_pair;
			}
			else {
				prev_pair_ptr = NULL;
				inserted = 0;
				for (pair_ptr = head_pair; pair_ptr != NULL;
					pair_ptr = pair_ptr -> next) {
					if (this_value <= pair_ptr -> value) {
						if (prev_pair_ptr == NULL) {
							/* put at beginning */
							new_pair -> next = head_pair;
							head_pair = new_pair;
						}
						else {
							/* put into middle */
							new_pair -> next = prev_pair_ptr -> next;
							prev_pair_ptr -> next = new_pair;
						}
						inserted = 1;
						break;
					}
					prev_pair_ptr = pair_ptr;
				}
				if (!inserted) {
					/* put at end */
					tail_pair -> next = new_pair;
					tail_pair = new_pair;
				}
			}

			/* put in info */
			new_pair -> i = i;
			new_pair -> j = j;
			new_pair -> vtx1 = vtx1;
			new_pair -> vtx2 = vtx2;
			new_pair -> value = this_value;
			new_pair -> dij = ij_distance;
			new_pair -> shape = shape;
			if (given_face -> shape == CONVEX) {
				for (k = 0; k < 3; k++) {
					new_pair -> c.center[k] = temporary_circle.center[k];
					new_pair -> c.axis[k] = temporary_circle.axis[k];
				}
				new_pair -> c.radius = temporary_circle.radius;
				new_pair -> a.cir = &(new_pair -> c);
				new_pair -> a.vtx[0] = vtx1;
				new_pair -> a.vtx[1] = vtx2;
				new_pair -> a.shape = CONVEX;
				new_pair -> a.small = 0;
			}
		}
	}
	if (head_pair == NULL)
		informd("select_pair: no pair found");

	ibest = -1;
	jbest = -1;
	n_contain_fail = 0;
	n_delta_fail = 0;
	n_pair = 0;
	for (pair_ptr = head_pair; pair_ptr != NULL;
		pair_ptr = pair_ptr -> next) n_pair++;
	for (pair_ptr = head_pair; pair_ptr != NULL;
		pair_ptr = pair_ptr -> next) {
		i = pair_ptr -> i;
		j = pair_ptr -> j;
		ij_distance = pair_ptr -> dij;
		vtx1 = pair_ptr -> vtx1;
		vtx2 = pair_ptr -> vtx2;
		shape = pair_ptr -> shape;
		if (given_face -> shape != CONVEX) {
			ibest = i;
			jbest = j;
			break;
		}
		/* check for arc midpoint outside face */
		middle (&(pair_ptr -> a), mid_point);
		if (error()) return(NULL);
		is_contained = point_in_face (mid_point, given_face, 0);
		if (error()) return(NULL);
		if (is_contained != 1) n_contain_fail++;
		if (is_contained != 1) continue;	/* no good */
		/* temporary edges for edge-edge angle checks */
		temporary_edge1.arcptr = &(pair_ptr -> a);
		temporary_edge1.orn = 0;
		temporary_edge1.fac = given_face;
		temporary_edge2.arcptr = &(pair_ptr -> a);
		temporary_edge2.orn = 1;
		temporary_edge2.fac = given_face;
		for (k = 0; k < 3; k++) {
			vtx1_vector[k] = (vtx1 -> center[k] - vty -> center[k]) /
				vty -> radii[0];
			vtx2_vector[k] = (vtx2 -> center[k] - vty -> center[k]) /
				vty -> radii[0];
		}
		edg1 = ((i > 0) ? *(edge_list + i - 1) :
			*(edge_list + n_edge_f - 1));
		edg2 = *(edge_list + i);
		edg3 = ((j > 0) ? *(edge_list + j - 1) :
			*(edge_list + n_edge_f - 1));
		edg4 = *(edge_list + j);
		delta1 = edge_delta (edg1, &temporary_edge1, vtx1_vector);
		delta2 = edge_delta (&temporary_edge2, edg2, vtx1_vector);
		delta3 = edge_delta (edg3, &temporary_edge2, vtx2_vector);
		delta4 = edge_delta (&temporary_edge1, edg4, vtx2_vector);
		n_delta_fail++;
		if (n_edge_f > 6) {
			if (delta1 <= 0.0 || delta1 >= PI) continue;
			if (delta2 <= 0.0 || delta2 >= PI) continue;
			if (delta3 <= 0.0 || delta3 >= PI) continue;
			if (delta4 <= 0.0 || delta4 >= PI) continue;
		}
		n_delta_fail--;
		ibest = i;
		jbest = j;
		break;
	}

	/* check for failure */

	if (ibest < 0 || jbest < 0) {
		free_pointers (VERTEX, vertex_list);
		free_longs (bad_vertex, 0, BAD_VERTEX);
		free_pointers (ARC, arc_list);
		free_pointers (EDGE, edge_list);
		free_doubles (perimeter, 0, PERIMETER);
		free_doubles (arc_direction, 0, ARC_DIRECTION);
		free_doubles (median, 0, MEDIAN);
		free_doubles (distance_backward, 0, DISTANCE_BACKWARD);
		free_doubles (distance_forward, 0, DISTANCE_FORWARD);
		next_pair_ptr = NULL;
		for (pair_ptr = head_pair; pair_ptr != NULL;
			pair_ptr = next_pair_ptr) {
			next_pair_ptr = pair_ptr -> next;
			free_object (VERTEX_PAIR, (short *) pair_ptr);
		}
		return (NULL);
	}

	small = 0;
	/* vertices to connect */
	vtx1 = *(vertex_list + ibest);
	vtx2 = *(vertex_list + jbest);

	/* store edge pointers in calling routine */

	*e1 = ((ibest > 0) ? *(edge_list + ibest - 1) :
		*(edge_list + n_edge_f - 1));
	*e2 = *(edge_list + ibest);
	*e3 = ((jbest > 0) ? *(edge_list + jbest - 1) :
		*(edge_list + n_edge_f - 1));
	*e4 = *(edge_list + jbest);

	/* vectors from center to endpoints of new arc */
	for (k = 0; k < 3; k++) {
		i_vector[k] = vtx1 -> center[k] - vty -> center[k];
		j_vector[k] = vtx2 -> center[k] - vty -> center[k];
	}

	/* center and axis */
	if (vty -> type == TORUS || vty -> type == CYLINDER) {
		/* if shape of proposed arc not already set (to convex),
			then short and straight */
		if (shape == 0) {
			shape = STRAIGHT;
			small = 1;
			for (k = 0; k < 3; k++) {
				circle_center[k] =
					(vtx1 -> center[k] + vtx2 -> center[k]) / 2;
				circle_axis[k] = 0.0;
			}
		}
		else {
			/* convex arc joining concave (or straight) arcs */
			doti = dot_product (i_vector, vty -> axis);
			for (k = 0; k < 3; k++)
				circle_center[k] =
					vty -> center[k] + doti * vty -> axis[k];
			sign = ((*(arc_direction + jbest) > 0.0) ? 1 : -1);
			for (k = 0; k < 3; k++)
				circle_axis[k] = sign * vty -> axis[k];
		}
	}
	else if (vty -> type == SPHERE) {
		shape = given_face -> shape;		/* convex or concave */
		small = 0;
		for (k = 0; k < 3; k++)
			circle_center[k] = vty -> center[k];
		cross (i_vector, j_vector, circle_axis);
		normalize (circle_axis);
	}
	else {
		ex = new_cept (PARAMETER_ERROR,  INVALID_VALUE,  FATAL_SEVERITY);
        add_object (ex, VARIETY, "vty");
		add_function (ex, "select_pair");
		add_source (ex, "msbisect.c");
        add_long (ex, "vty -> type", (long) vty -> type);
		return (NULL);
	}

	/* radius */
	circle_radius = distance (circle_center, vtx1 -> center);

	/* new circle and arc */

	cir = new_circle (circle_center, circle_radius, circle_axis);
	if (error()) return(NULL);
	link_circle (this_srf, cir);
	arc_new = new_arc (cir, vtx1, vtx2, shape, small, alpha, lfn, ofn, ofn);
	if (error()) return(NULL);
	link_arc (this_srf, arc_new);

	/* check for arc meeting cycle at bad angles */
	/* later: add code */

	free_pointers (VERTEX, vertex_list);
	free_longs (bad_vertex, 0, BAD_VERTEX);
	free_pointers (ARC, arc_list);
	free_pointers (EDGE, edge_list);
	free_doubles (perimeter, 0, PERIMETER);
	free_doubles (arc_direction, 0, ARC_DIRECTION);
	free_doubles (median, 0, MEDIAN);
	free_doubles (distance_backward, 0, DISTANCE_BACKWARD);
	free_doubles (distance_forward, 0, DISTANCE_FORWARD);
	for (pair_ptr = head_pair; pair_ptr != NULL;
		pair_ptr = next_pair_ptr) {
		next_pair_ptr = pair_ptr -> next;
		free_object (VERTEX_PAIR, (short *) pair_ptr);
	}

	return (arc_new);
}



int first_and_third (double v1[3], double v2[3], double v3[3], double v4[3])
{
	double min_ang;
	double angles[4];
	int i, min_v;
    struct cept *ex;

	angles[0] = quad_angle (v1, v2, v3, v4);
	angles[1] = quad_angle (v2, v3, v4, v1);
	angles[2] = quad_angle (v3, v4, v1, v2);
	angles[3] = quad_angle (v4, v1, v2, v3);
	min_ang = PI;
	min_v = -1;
	for (i = 0; i < 4; i++)
		if (angles[i] < min_ang) {
			min_ang = angles[i];
			min_v = i;
		}
	if (min_v < 0) {
		ex = new_cept (LOGIC_ERROR,  NOT_FOUND,  FATAL_SEVERITY);
		add_function (ex, "first_and_third");
		add_source (ex, "msbisect.c");
		add_message (ex, "no minimum angle");
		return (0);
	}

	return ((min_v % 2) == 0);
}


/*
   MSRoll
   Copyright 1996 by Michael L. Connolly
   All Rights Reserved

*/
