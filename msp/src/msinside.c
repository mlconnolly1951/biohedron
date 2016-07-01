/*
 * Molecular Surface Package
 * Surface Rendering by Foliation
 * Copyright 1986 by Michael L. Connolly
 * All Rights Reserved

 * January 5, 2002
 */


#include "ms.h"
#include "msenum.h"
#include "mslimit.h"
#include "msdefault.h"
#include "msstruct.h"
#include "msglobal.h"
#include "msfunction.h"


double compute_omega (double probe_radius, struct face *fac) 	/* calculate solid angle of face */
{
	int orn, k, prev_orn, chi, ne, result;
	long cycle_number, edge_number;
	double fsgn, radius, phi, theta, beta, omega, cir_dist;
	double outward[3], tangent1[3], tangent2[3], center[3];
	char shape_name[40];
	struct variety *vty;
	struct arc *arc_ptr, *prev_arc;
	struct edge *edg;
	struct vertex *vtx;
	struct circle *cir;
	struct cycle *cyc;
	struct sphere *atm;
	struct probe *prb;
    struct cept *ex;

	atm = NULL;
	prb = NULL;
	vty = NULL;
	prev_orn = -1;
	cycle_number = 0;
	edge_number = 0;
	if (fac -> shape == CONVEX) {
		fsgn = 1.0;
		vty = fac -> vty;
		atm = fac -> ptr.atm;
		if (vty != NULL) {
			radius = vty -> radii[0];
			for (k = 0; k < 3; k++)
				center[k] = vty -> center[k];
		}
		else if (atm != NULL) {
			radius = atm -> radius;
			for (k = 0; k < 3; k++)
				center[k] = atm -> center[k];
		}
		else return (0.0);
	}
	else if (fac -> shape == CONCAVE) {
		fsgn = -1.0;
		vty = fac -> vty;
		prb = fac -> ptr.prb;
		if (vty != NULL) {
			radius = vty -> radii[0];
			for (k = 0; k < 3; k++)
				center[k] = vty -> center[k];
		}
		else if (prb != NULL) {
			radius = probe_radius;
			for (k = 0; k < 3; k++)
				center[k] = prb -> center[k];
		}
	}
	else return (0.0);

	if (radius <= 0.0) return (0.0);

	/* All the way with Gauss-Bonnet! */

	chi = 2;
	for (cyc = fac -> first_cycle; cyc != NULL; cyc = cyc -> next) chi--;
	fac -> chi = (short) chi;
	omega = 2 * PI * chi;

	cycle_number = 0;
	for (cyc = fac -> first_cycle; cyc != NULL; cyc = cyc -> next) {
		/* find last edge of linked list */
		edge_number = 0;
		for (edg = cyc -> first_edge; edg != NULL; edg = edg -> next) {
			prev_arc = edg -> arcptr;
			prev_orn = edg -> orn;
			if (prev_orn < 0 || prev_orn > 1) {
				ex = new_cept (PARAMETER_ERROR, INVALID_VALUE, FATAL_SEVERITY);
				add_function (ex, "compute_omega");
				add_source (ex, "msinside.c");
				add_long (ex, "prev_orn", (long) prev_orn);
				add_long (ex, "cycle_number", cycle_number);
				add_long (ex, "edge_number", edge_number);
				add_long (ex, "face shape", (long) fac -> shape);
				add_long (ex, "probe pointer", (long) prb);
				add_message (ex, "unable to find last edge of linked list");
				return(0.0);
			}
			edge_number++;
		}
		edge_number = 0;
		for (edg = cyc -> first_edge; edg != NULL; edg = edg -> next) {
			arc_ptr = edg -> arcptr;
			orn = edg -> orn;
			arc_ptr -> phi = arc_ang (arc_ptr);
			phi = arc_ptr -> phi;
			vtx = arc_ptr -> vtx[orn];
			if (orn < 0 || orn > 1) {
				ex = new_cept (PARAMETER_ERROR, INVALID_VALUE, FATAL_SEVERITY);
				add_function (ex, "compute_omega");
				add_source (ex, "msinside.c");
				add_long (ex, "orn", (long) orn);
				return(0.0);
			}
			if (prev_orn < 0 || prev_orn > 1) {
				ex = new_cept (PARAMETER_ERROR, INVALID_VALUE, FATAL_SEVERITY);
				add_function (ex, "compute_omega");
				add_source (ex, "msinside.c");
				add_long (ex, "prev_orn", (long) prev_orn);
				return(0.0);
			}
			if (vtx == NULL) beta = 0.0;
			else if (arc_ptr -> cir == prev_arc -> cir) beta = 0.0;
			else {
				/* computation of beta */
				result = get_tangent (prev_arc, prev_orn, 1 - prev_orn, tangent1);
				if (!result) {
					add_function(tail_cept, "compute_omega");
					return(0.0);
				}
				result = get_tangent (arc_ptr, orn, orn, tangent2);
				if (!result) {
					add_function(tail_cept, "compute_omega");
					return(0.0);
				}
				for (k = 0; k < 3; k++)
					outward[k] = vtx -> center[k] - center[k];
				normalize (outward);
				beta = odd_angle (tangent1, tangent2, outward, fsgn);
				if (error()) return(0.0);
			}
			cir = arc_ptr -> cir;
			cir_dist = distance (cir -> center, center);
			if (cir_dist > radius) {
                ex = new_cept (GEOMETRY_ERROR, INCONSISTENCY, FATAL_SEVERITY);
                add_source (ex, "msinside.c");
                add_function (ex, "compute_omega");
				if (fac -> shape == CONVEX)
					add_message (ex, "convex face");
				else if (fac -> shape == CONCAVE)
					add_message (ex, "concave face");
				add_message (ex, "badd circle pointer");
                add_double (ex, "cir_dist", cir_dist);
                add_double (ex, "radius", radius);
				return (0.0);
			}
			if (vty != NULL)
				theta = theta_circle (cir, orn, vty, fsgn);
			else theta = cir -> theta;
			omega = omega + phi * sin (theta);
			if (beta < - 5 * PI / 6) beta += 2 * PI;
			omega -= beta;
			prev_arc = edg -> arcptr;
			prev_orn = edg -> orn;
			edge_number++;
		}
		cycle_number++;
	}

	if (omega > 4 * PI + EPSILON) {
		if (chi == 1) { /* assume tiny face */
			omega = EPSILON;
			return (omega);
		}
		ex = new_cept (GEOMETRY_ERROR, INVALID_VALUE, FATAL_SEVERITY);
		add_source (ex, "msinside.c");
		add_function (ex, "compute_omega");
		add_message (ex, "invalid solid angle");
        add_double (ex, "omega", omega);
        add_long (ex, "chi", (long) chi);
		return (- omega);
	}
	if (omega > 4 * PI) omega = 4 * PI;

	if (omega < -EPSILON) {
		ne = 0;
		for (cyc = fac -> first_cycle; cyc != NULL; cyc = cyc -> next) {
			ne += edges_in_cycle (cyc);
		}
		return (omega);
	}
	if (omega < 0.0) omega = 0.0;
	return (omega);
}


double numerical_omega (struct surface *msphn, struct variety *vty)
{
	int nvert, nhor, nequat, nu, k, i, j, inside, n_inside, n_fail;
	double nomega;
	double ux, uy, uz, uxy, fi, fj;
	double pnt[3];
	double units[MAX_UNIT][3];

	/* generate n unit vectors equally spaced over evaluation sphere */
	nvert = 0.5 * sqrt(MAX_UNIT * PI);
	if (nvert <= 0) return (-1.0);
	nequat = 2 * nvert;
	nu = 0;
	for (i = 0; i < nvert; i++) {
		fi = (PI * i) / nvert;
		uz = cos (fi);
		uxy = sin (fi);
		nhor = nequat * uxy;
		if (nhor <= 0) nhor = 1;
		for (j = 0; j < nhor - 1; j++) {
			fj = (2 * PI * j) / nhor;
			ux = cos(fj) * uxy;
			uy = sin(fj) * uxy;
			if (nu >= MAX_UNIT) break;
			units[nu][0] = ux;
			units[nu][1] = uy;
			units[nu][2] = uz;
			nu = nu + 1;
		}
		if (nu >= MAX_UNIT) break;
	}
	if (nu <= 0) return (-1.0);
	
	/* count points on evaluation sphere inside polyhedron */
	n_inside = 0; n_fail = 0;
	for (i = 0; i < nu; i++) {
		for (k = 0; k < 3; k++)
			pnt[k] = vty -> center[k] + vty -> radii[0] * units[i][k];
		/* check each point for being inside polyhedron */
		inside = point_in_polyhedron (msphn, pnt);
		if (inside < 0) {
			n_fail++;
			if (n_fail > 10) {
				inform ("numerical_omega: too many failures, give up");
				return (0.0);
			}
		}
		if (inside) n_inside++;
	}
	nomega = 4 * PI * (double) n_inside / (double) nu;
	return (nomega);
}



/* point in polyhedron? */

int point_in_polyhedron (struct surface *msphn, double pnt[3])
{
	int is_inside;
	double solid;

	solid = polyhedron_subtends (msphn, pnt);
	is_inside = (solid > 2 * PI);
	return (is_inside);
}


double polyhedron_subtends (struct surface *msphn, double pnt[3])
{
	int i, j, orn;
	double solid, delta;
	struct phnedg *edg;
	struct phntri *tri;
	struct phnvtx *vs[3];

	/* initialize solid angle subtended by polyhedron as seen by point */
	solid = 0.0;

	for (i = 0; i < msphn -> n_phntri; i++) {
		tri = num2phntri (msphn, i + 1);
		if (tri == NULL) return (0.0);
		for (j = 0; j < 3; j++) {
			edg = tri -> edg[j];
			orn = tri -> orn[j];
			vs[j] = edg -> pvt[orn];
		}
		/* compute solid angle of triangles as seen by point */
		delta = tetra_solid_angle (pnt, vs[0] -> center, vs[1] -> center, vs[2] -> center);
		solid += delta;
	}
	return (solid);
}

/* is point in triangle? */

int point_in_triangle (double pnt[3], struct phntri *tri)
{
	int j, k, orn;
	int is_inside;
	double winding;
	double vects[3][3];
	struct phnedg *e;
	struct phnvtx *vtx;


	/* compute vectors from center to vertices */

	for (j = 0; j < 3; j++) {
		e = tri -> edg[j];
		orn = tri -> orn[j];
		vtx = e -> pvt[orn];
		for (k = 0; k < 3; k++)
			vects[j][k] = vtx -> center[k] - pnt[k];
		if (norm (vects[j]) <= 0.0) {
			is_inside = 0;
			return (is_inside);
		}
		if (!normalize (vects[j])) {
			is_inside = 0;
			return (is_inside);
		}
	}
	/* compute winding number */
	winding = 0.0;
	for (j = 0; j < 3; j++) {
		k = ((j < 2) ? j + 1 : 0);
		winding += odd_angle (vects[j], vects[k], tri -> axis, (double) 1.0);
	}
	is_inside = (winding > PI);
	return (is_inside);
}


/* does face fac contain point pnt ? */

int point_in_face (double pnt[3], struct face *given_face, int rendering)
{
	int k, result, one_cycle;
	double sphere_radius, face_sign;
	double sphere_center[3];
	struct variety *vty;
	struct cycle *cyc;

	vty = given_face -> vty;
	face_sign = ((given_face -> shape == CONVEX) ? 1.0 : -1.0);
	for (k = 0; k < 3; k++)
		sphere_center[k] = vty -> center[k];
	sphere_radius = vty -> radii[0];
	one_cycle = 0;
	if (given_face -> first_cycle != NULL)
		if (given_face -> first_cycle -> next == NULL) one_cycle = 1;
	for (cyc = given_face -> first_cycle; cyc != NULL; cyc = cyc -> next) {
		result = cycle_contain_point (cyc, pnt, sphere_center, sphere_radius, face_sign, one_cycle, rendering);
		if (error()) return(0);
		if (result != 1) return (result);
	}
	return (1);
}

/* does cycle given_cycle contain point given_point ? */

int cycle_contain_point (struct cycle *given_cycle, double given_point[3], double center[3], double radius, double face_sign, int one_cycle, int rendering)
{
	int n_edge_cyc, k, ivtx, jvtx, orn;
	double rotation_angle, dot_prod, extension, side_length;
	double vector[3], south[3];
	double *projected_vertices, *polygon_side, *exterior_angle;
	struct arc *this_arc;
	struct edge *edg;
	struct vertex *vtx;
    struct cept *ex;

	/* count number of arcs in cycle */
	n_edge_cyc = 0;
	for (edg = given_cycle -> first_edge; edg != NULL;
		edg = edg -> next) n_edge_cyc++;

	/* quick containment check */
	if (rendering) {
		if (!quick(given_cycle, given_point, face_sign)) return (0);
		if (error()) return(0);

		/* if < 3 arcs, quick containment sufficient */
		if (n_edge_cyc < 3) return (1);

		/* if <= 3 arcs, one cycle for face,
		   and face is concave, then quick containment sufficient */
		if (one_cycle && face_sign < 0.0 && n_edge_cyc == 3) return (1);
	}

	/* stereographic projection */
	for (k = 0; k < 3; k++)
		south[k] = (center[k] - given_point[k]) / radius;
	/* allocate memory for projected vertices of a */
	projected_vertices = allocate_doubles (n_edge_cyc * 3, 0, PROJECTED_VERTICES);
	if (projected_vertices == NULL) {
		ex = new_cept (MEMORY_ERROR,  ALLOCATION,  FATAL_SEVERITY);
        add_variable (ex, PROJECTED_VERTICES, "projected vertices");
        add_function (ex, "cycle_contain_point");
        add_source (ex, "msinside.c");
		return (0);
	}
	polygon_side = allocate_doubles (n_edge_cyc * 3, 0, POLYGON_SIDE);
	if (polygon_side == NULL) {
		ex = new_cept (MEMORY_ERROR,  ALLOCATION,  FATAL_SEVERITY);
        add_variable (ex, POLYGON_SIDE, "polygon side");
        add_function (ex, "cycle_contain_point");
        add_source (ex, "msinside.c");
		return (0);
	}
	exterior_angle = allocate_doubles (n_edge_cyc, 0, EXTERIOR_ANGLE);
	if (exterior_angle == NULL) {
		ex = new_cept (MEMORY_ERROR,  ALLOCATION,  FATAL_SEVERITY);
        add_variable (ex, EXTERIOR_ANGLE, "exterior angle");
        add_function (ex, "cycle_contain_point");
        add_source (ex, "msinside.c");
		return (0);
	}
	/* project vertices */
	ivtx = 0;
	for (edg = given_cycle -> first_edge; edg != NULL; edg = edg -> next) {
		this_arc = edg -> arcptr;
		orn = edg -> orn;
		vtx = this_arc -> vtx[orn];
		for (k = 0; k < 3; k++)
			vector[k] = vtx -> center[k] - given_point[k];
		dot_prod = dot_product (vector, south);
		if (dot_prod <= 0.0) {
			/* bad projection */
			free_doubles (projected_vertices, 0, PROJECTED_VERTICES);
			free_doubles (polygon_side, 0, POLYGON_SIDE);
			free_doubles (exterior_angle, 0, EXTERIOR_ANGLE);
			return (-1);
		}
		extension = 2 * radius / dot_prod;
		for (k = 0; k < 3; k++)
			*(projected_vertices + 3 * ivtx + k) =
				given_point[k] + extension * vector[k];
		ivtx++;
	}
	/* later: add code for null edges */
	/* compute polygon_side vectors */
	for (ivtx = 0; ivtx < n_edge_cyc; ivtx++) {
		jvtx = ivtx + 1;
		if (jvtx >= n_edge_cyc) jvtx = 0;
		for (k = 0; k < 3; k++)
			*(polygon_side+3*ivtx+k) =
				*(projected_vertices+3*jvtx+k) -
					*(projected_vertices+3*ivtx+k);
		side_length = norm (polygon_side + 3 * ivtx);
		if (side_length > 0.0) {
			for (k = 0; k < 3; k++)
				*(polygon_side + 3 * ivtx + k) /= side_length;
        }
	}

	/* compute angles */
	for (ivtx = 0; ivtx < n_edge_cyc; ivtx++) {
		jvtx = ivtx - 1;
		if (jvtx < 0) jvtx = n_edge_cyc - 1;
		side_length = norm (polygon_side + 3 * ivtx);
		if (side_length <= 0.0) {
			*(exterior_angle+ivtx) = 0.0;
			continue;
		}
		side_length = norm (polygon_side + 3 * jvtx);
		if (side_length <= 0.0) {
			*(exterior_angle+ivtx) = 0.0;
			continue;
		}
		*(exterior_angle+ivtx) =
			odd_angle (polygon_side + 3 * jvtx, polygon_side+3*ivtx,
				south, -1.0);
	}
	/* calculate rotation angle */
	rotation_angle = 0.0;
	for (ivtx = 0; ivtx < n_edge_cyc; ivtx++)
		rotation_angle += *(exterior_angle + ivtx);
	/* free temporary memory */
	free_doubles (projected_vertices, 0, PROJECTED_VERTICES);
	free_doubles (polygon_side, 0, POLYGON_SIDE);
	free_doubles (exterior_angle, 0, EXTERIOR_ANGLE);
	return (face_sign * rotation_angle > 0.0);
}
	

/* quick containment check */

int quick (struct cycle *cyc, double pnt[3], double fsgn)
{
	int k, orn, sgn;
	double dt;
	double vect[3];
	struct edge *edg;
	struct arc *a;
	struct circle *cir;

	/* check whether point is on accessible side of each circle plane */

	for (edg = cyc -> first_edge; edg != NULL; edg = edg -> next) {
		a = edg -> arcptr;
		orn = edg -> orn;
		cir = a -> cir;
		sgn = 1 - 2 * orn;
		for (k = 0; k < 3; k++)
			vect[k] = pnt[k] - cir -> center[k];
		dt = dot_product (vect, cir -> axis);
		if (dt * sgn * fsgn < 0.0) return (0);
	}
	return (1);
}

