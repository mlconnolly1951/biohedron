#include "ms.h"
#include "msenum.h"
#include "mslimit.h"
#include "msdefault.h"
#include "msstruct.h"
#include "msglobal.h"
#include "msfunction.h"

/* Copyright 1995 by Michael L. Connolly */
/* Last revised: January 5, 2002 */

/* SURFACE GENERATION ROUTINES: */


/* initialize the surface */
void initialize_surface (struct surface *this_srf)
{
	struct sphere *atm_ptr;
	
	this_srf -> n_circle = this_srf -> n_arc = 0;
	this_srf -> n_vertex = this_srf -> n_face = 0;
	this_srf -> n_edge = this_srf -> n_cycle = 0;
	this_srf -> head_circle = this_srf -> tail_circle = NULL;
	this_srf -> head_vertex = this_srf -> tail_vertex = NULL;
	this_srf -> head_face = this_srf -> tail_face = NULL;
	
	for (atm_ptr = (this_srf -> head_atom);
		atm_ptr != NULL;
		atm_ptr = atm_ptr -> next) {
		atm_ptr -> problem = FALSE;
		atm_ptr -> first_arc = NULL;
		atm_ptr -> first_face = NULL;
	}

	this_srf -> n_component = 0;
	this_srf -> n_problem_face = 0;
	this_srf -> n_problem_atom = 0;
	this_srf -> molecule_area = 0.0;
	this_srf -> total_contact_area = 0.0;
	this_srf -> total_reentrant_area = 0.0;
	this_srf -> total_accessible_area = 0.0;
	this_srf -> surface_volume = 0.0;
	this_srf -> polyhedron_volume = 0.0;
	this_srf -> old_volume = 0.0;
	this_srf -> mol_vol = 0.0;
	this_srf -> surface_completed = 0;
	this_srf -> volume_computed = 0;
	this_srf -> old_volume_computed = 0;
	
}

/* calculate an analytical molecular surface */
void calculate_surface (struct surface *this_srf)
{
	char message[MAXLINE];
    struct cept *ex;

	sprintf (message,"%8.4f probe radius", this_srf -> probe_radius);
	inform (message);
	initialize_surface (this_srf);	/* initialize surface */
	probe_rolling (this_srf); if (error()) return;
	create_permanent (this_srf); if (error()) return;
	free_pairs (this_srf);
	free_centrals (this_srf);
	create_saddles (this_srf); if (error()) return;
	create_convex (this_srf); if (error()) return;
	if (this_srf -> probe_radius > 0.0) {
		cusp_trimming (this_srf); if (error()) return;
	}
	finish_up (this_srf); if (error()) return;
	count_problem_faces2 (this_srf); if (error()) return;
	number_arcs (this_srf); if (error()) return;
	clean_circles (this_srf); if (error()) return;
	count_edges (this_srf); if (error()) return;

	sprintf (message,"%8ld faces", this_srf -> n_face); inform(message);
	sprintf (message,"%8ld problem faces, %5ld problem atoms",
		this_srf -> n_problem_face, this_srf -> n_problem_atom);
	if (this_srf -> n_problem_face > 0)
		inform(message);
    if (this_srf -> n_face <= 0) {
		ex = new_cept (GEOMETRY_ERROR, MSUNDERFLOW, FATAL_SEVERITY);
        add_function (ex, "calculate_surface");
        add_source (ex, "msroll.c");
        add_long (ex, "number of faces", (long) this_srf -> n_face);
    }
}

/* PART A: ANALYTICAL SURFACE CALCULATION */


void probe_rolling (struct surface *this_srf)
{
	this_srf -> head_central = NULL;
	this_srf -> tail_central = NULL;
	this_srf -> n_central = 0;
	
	if (this_srf -> use_grid) setup_sgrid (this_srf); if (error()) return;
	create_neighbors (this_srf); if (error()) return;
	if (this_srf -> use_grid) free_cgrid (this_srf -> sg);
	create_tori (this_srf); if (error()) return;
	if (this_srf -> use_grid) {
		init_bgrid (this_srf -> sg, this_srf -> bit_width);
		if (error()) return;
	}
	create_probes (this_srf); if (error()) return;
	if (this_srf -> use_grid) {
		free_sgrid (this_srf -> sg);
		free_bgrid (this_srf -> sg);
	}
	find_quartets (this_srf);
	create_concave_faces (this_srf); if (error()) return;
	mark_not_buried (this_srf); if (error()) return;
	free_nbr (this_srf); if (error()) return;
	create_centrals (this_srf); if (error()) return;
}

void setup_sgrid (struct surface *this_srf)
{
	long i; int k;
	long nspheres;
	double *centers;
	double *radii;
    struct cept *ex;

	struct sphere *sphere_ptr;
	centers = allocate_doubles (this_srf -> n_atom * 3, 0, CENTERS);
	if (centers == NULL) {
		ex = new_cept ( MEMORY_ERROR,  ALLOCATION,  FATAL_SEVERITY);
        add_variable (ex, CENTERS, "sphere grid centers");
        add_function (ex, "setup_sgrid");
        add_source (ex, "msroll.c");
		return;
	}
	radii = allocate_doubles (this_srf -> n_atom, 0, RADII);
	if (radii == NULL) {
		ex = new_cept ( MEMORY_ERROR,  ALLOCATION,  FATAL_SEVERITY);
        add_variable (ex, RADII, "sphere grid radii");
        add_function (ex, "setup_sgrid");
        add_source (ex, "msroll.c");
		return;
	}

	for (sphere_ptr = (struct sphere *) (this_srf -> head_atom), i = 0; sphere_ptr != NULL;
		sphere_ptr = sphere_ptr -> next, i++) {
		for (k = 0; k < 3; k++)
			*(centers + 3 * i + k) = sphere_ptr -> center[k];
		*(radii + i) = sphere_ptr -> radius + this_srf -> probe_radius;
	}
	nspheres = this_srf -> n_atom;
	this_srf -> sg = init_sgrid (nspheres, centers, radii);
	init_cgrid (this_srf -> sg, this_srf -> cube_width);
	free_doubles (centers, 0, CENTERS);
	free_doubles (radii, 0, RADII);
}

/* create neighbor lists */
void create_neighbors (struct surface *this_srf)
{
	int i, count, c;
    int are;
	long n_of_neighbors, total_neighbors, sn, sn2;
	unsigned size;
	unsigned long lsize, l, lnsize;
	long *numbers;
	double expanded_radius;
	char message[MAXLINE];
	struct sphere *sphere1, *sphere2;
	struct sphere *sphere_ptr;
	struct sphere **sphere_handles;
	struct neighbor *neighbor_hdl, *nbr_hdl;
    struct cept *ex;

	
	/* allocate temporary memory to store sphere handles */
	sphere_handles = (struct sphere **)
		allocate_pointers (SPHERE, this_srf -> n_atom);
	if (sphere_handles == NULL) {
		ex = new_cept ( MEMORY_ERROR,  ALLOCATION,  FATAL_SEVERITY);
        add_object (ex, SPHERE, "sphere handles");
        add_function (ex, "create_neighbors");
        add_source (ex, "msroll.c");
		return;
	}
	sn = 0;
	for (sphere_ptr = (struct sphere *) (this_srf -> head_atom); sphere_ptr != NULL;
		sphere_ptr = sphere_ptr -> next) {
		sphere_ptr -> buried = TRUE;
		sphere_ptr -> first_neighbor = NULL;
		sphere_ptr -> last_neighbor = NULL;
		expanded_radius = sphere_ptr -> radius + this_srf -> probe_radius;
		sphere_ptr -> ers = expanded_radius * expanded_radius;
		*(sphere_handles + sn++) = sphere_ptr;
		if (sn != sphere_ptr -> number) {
			ex = new_cept (LOGIC_ERROR,  INCONSISTENCY,  FATAL_SEVERITY);
			add_object (ex, SPHERE, "neighbors");
			add_function (ex, "create_neighbors");
			add_source (ex, "msroll.c");
			add_long (ex, "sphere number", sn);
			add_atom (ex, sphere_ptr);
			return;
		}
	}
	this_srf -> pair_array = NULL;
	this_srf -> n_pair = 0;
	this_srf -> n_tori = this_srf -> n_probe = 0;
	this_srf -> head_torus = this_srf -> tail_torus = NULL;
	this_srf -> head_probe = this_srf -> tail_probe = NULL;

	/* allocate temporary memory to store neighbors */
	neighbor_hdl = (struct neighbor *) allocate_objects (NEIGHBOR, this_srf -> n_atom);
	if (neighbor_hdl == NULL) {
		ex = new_cept (MEMORY_ERROR,  ALLOCATION,  FATAL_SEVERITY);
        add_object (ex, NEIGHBOR, "neighbor handles");
        add_function (ex, "create_neighbors");
		add_source (ex, "msroll.c");
		return;
	}
	total_neighbors = 0;
	lnsize = 0L;
	numbers = allocate_longs (this_srf -> n_atom, 0, NUMBERS);
	if (numbers == NULL) {
		ex = new_cept (MEMORY_ERROR,  ALLOCATION,  FATAL_SEVERITY);
        add_object (ex, NEIGHBOR, "neighbor numbers");
        add_function (ex, "create_neighbors");
			add_source (ex, "msroll.c");
		return;
	}
	/* loop through spheres */
	for (sphere1 = (struct sphere *) (this_srf -> head_atom); sphere1 != NULL;
		sphere1 = sphere1 -> next) {
		n_of_neighbors = 0;
		if (this_srf -> use_grid) {
			sn = sphere1 -> number;
			count = sphere_inquiry (this_srf -> sg, sn,  this_srf -> n_atom, numbers);
			if (error()) {
				add_function (tail_cept, "create_neighbors");
				add_source (ex, "msroll.c");
				return;
			}
			/* loop through spheres in box enclosing sphere */
			nbr_hdl = neighbor_hdl;
			for (c = 0; c < count; c++) {
				sn2 = *(numbers + c);
				if (sn2 <= 0 || sn2 > this_srf -> n_atom) {
					ex = new_cept (LOGIC_ERROR,  BOUNDS,  FATAL_SEVERITY);
					add_object (ex, SPHERE, "neighbors");
					add_function (ex, "create_neighbors");
					add_source (ex, "msroll.c");
					add_long (ex, "sphere number", sn2);
					add_long (ex, "n_atom", this_srf -> n_atom);
					return;
				}
				sphere2 = *(sphere_handles + sn2 - 1);
				if (sphere1 == sphere2) continue;
                are = are_neighbors (this_srf, sphere1, sphere2);
                if (error()) {
                    if (tail_cept != NULL) {
                        add_function (tail_cept, "create_neighbors");
						add_source (ex, "msroll.c");
                        return;
                    }
                }
				if (are) {
					nbr_hdl -> sphptr = sphere2; /* store pointer to nbr */
					nbr_hdl++;
				}
			}
		}
		else {
			/* loop through spheres in a nested loop */
			for (sphere2 = (struct sphere *) (this_srf -> head_atom),
				nbr_hdl = neighbor_hdl;
				sphere2 != NULL; sphere2 = sphere2 -> next)
				if (sphere1 != sphere2) {
					are = are_neighbors (this_srf, sphere1, sphere2);
					if (error()) {
						if (tail_cept != NULL) {
							add_function (tail_cept, "create_neighbors");
							add_source (ex, "msroll.c");
							return;
						}
					}
                    if (are) {
						nbr_hdl -> sphptr = sphere2; /* store pointer to nbr */
						nbr_hdl++;
					}
				}
		}
		/* number of neighbors */
		n_of_neighbors = nbr_hdl - neighbor_hdl;
		if (n_of_neighbors <= 0) {
			/* no neighbors */
			sprintf (message, "%8ld atom has no neighbors", sphere1 -> number);
			informd(message);
			sphere1 -> first_neighbor = sphere1 -> last_neighbor = NULL;
			sphere1 -> buried = FALSE;	/* sphere w/o nbrs ! buried */
		}
		else {
			total_neighbors += n_of_neighbors;
			/* allocate neighbor list */
			size = n_of_neighbors * sizeof (struct neighbor);
			lnsize += size;
			sphere1 -> first_neighbor = (struct neighbor *)
				allocate_objects (NEIGHBOR, n_of_neighbors);
			if (sphere1 -> first_neighbor == (struct neighbor *) NULL) {
				ex = new_cept (MEMORY_ERROR,  ALLOCATION,  FATAL_SEVERITY);
				add_object (ex, SPHERE, "neighbors");
				add_function (ex, "create_neighbors");
				add_source (ex, "msroll.c");
				add_atom (ex, sphere1);
				add_long (ex, "lnsize", lnsize);
				return;
			}
			for (i = 0, nbr_hdl = neighbor_hdl; i < n_of_neighbors;
				i++, nbr_hdl++)
				(sphere1 -> first_neighbor + i) -> sphptr = nbr_hdl -> sphptr;
			sphere1 -> last_neighbor = sphere1 -> first_neighbor + n_of_neighbors - 1;
		}
	}
	if (total_neighbors % 2 != 0) {
		ex = new_cept (LOGIC_ERROR,  INVALID_VALUE,  FATAL_SEVERITY);
		add_function (ex, "create_neighbors");
		add_source (ex, "msroll.c");
        add_message (ex, "total number of neighbors not even");
		add_long (ex, "total_neighbors", total_neighbors);
		return;
	}
	if (total_neighbors <= 0) {
		this_srf -> n_pair = 0;
		this_srf -> pair_array = NULL;
		return;
	}
	free_objects (NEIGHBOR, (short *) neighbor_hdl);
	free_longs (numbers, 0, NUMBERS);
	free_pointers (SPHERE, sphere_handles);
	/* allocate temporary tori array */
	this_srf -> n_pair = total_neighbors / 2;
	lsize = this_srf -> n_pair * sizeof (struct pair);
	this_srf -> pair_array = (struct pair *) allocate_objects (PAIR, this_srf -> n_pair);
	if (this_srf -> pair_array == NULL) {
		ex = new_cept (MEMORY_ERROR,  ALLOCATION,  FATAL_SEVERITY);
		add_function (ex, "create_neighbors");
		add_source (ex, "msroll.c");
        add_object (ex, PAIR, "pair_array");
		return;
	}
	/* do the zeroing ourselves - historical reasons, not remembered */
	for (l = 0; l < lsize; l++)
		*((char *) (this_srf -> pair_array) + l) = 0;
}


/* place probe tangent to three spheres */
void create_probes (struct surface *this_srf)
{
	int k, n_mutual, should_break;
	char message[MAXLINE];
	struct pair *torus12, *torus13, *torus23, *pair_ptr;
	struct sphere *sphere1, *sphere2, *sphere3;
	struct mutual *mutuals, *last_mutual;
	struct mutual *neighbor3;
    struct cept *ex;

	if (this_srf -> n_pair <= 0) return;
	/* allocate temporary memory to store mutual neighbors */
	mutuals = (struct mutual *) allocate_objects (MUTUAL, this_srf -> n_atom);
	if (mutuals == NULL) {
		ex = new_cept (MEMORY_ERROR, ALLOCATION, FATAL_SEVERITY);
		add_function (ex, "create_probes");
		add_source (ex, "msroll.c");
        add_object (ex, MUTUAL, "mutuals");
		return;
	}

	/* loop through tmp tori */
	for (torus12 = this_srf -> pair_array; torus12 < this_srf -> pair_array +
		this_srf -> n_pair; torus12++) {
		sphere1 = torus12 -> sph[0];
		sphere2 = torus12 -> sph[1];
		/* get mutual neighbors of the two spheres */
		n_mutual = form_mutual_list (sphere1,sphere2, mutuals);
		if (n_mutual <= 0) {
			/* no neighbors implies not buried */
			torus12 -> buried = FALSE;
			/* probably not necessary: */
			torus12 -> free = TRUE;
			continue;
		}
		last_mutual = mutuals + n_mutual - 1;
		/* loop through mutual neighbors */
		for (neighbor3 = mutuals; neighbor3 <= last_mutual; neighbor3++) {
			sphere3 = neighbor3 -> sphptr; /* retrieve pointer to sphere */
			torus13 = neighbor3 -> tor13;
			torus23 = neighbor3 -> tor23;
			should_break = place_probe (this_srf, mutuals, last_mutual, 
				sphere1, sphere2, sphere3, torus12, torus13, torus23);
            if (error()) return;
            if (should_break) break;
		}	/* end of mutual neighbor loop */
	}
	/* free list of mutual neighbors */
	free_objects (MUTUAL, (short *) mutuals);
	/* mark free tori as not buried */
	for (pair_ptr = this_srf -> pair_array; pair_ptr < this_srf -> pair_array + this_srf -> n_pair; pair_ptr++)
		if (pair_ptr -> free) pair_ptr-> buried = FALSE;
}

/* place probe tangent to three spheres */
/* return true if the caller should break out of the mutual neighbor loop */
int place_probe (struct surface *this_srf, 
struct mutual *mutuals, struct mutual *last_mutual, 
struct sphere *sphere1, struct sphere *sphere2, struct sphere *sphere3, 
struct pair *torus12, struct pair *torus13, struct pair *torus23)
{
	int k, side, hit;
    double probe_radius;
	double distance13, dot1213, angle123, sin_angle123, radius12;
	double height, special_radius, special_distance;
	double sphere1_center[3], sphere2_center[3], sphere3_center[3];
	double axis13[3], axis12[3];
	double torus12_center[3], torus13_center[3];
	double vector_torus_base[3], base_point[3];
	double probe_position[3], altitude[3];
    struct probe *prb;
    struct cept *ex;

	/* computations for probe placement */
	/* transfer data to local variables */
    probe_radius = this_srf -> probe_radius;
	for (k = 0; k < 3; k++) {
		sphere1_center[k] = sphere1 -> center[k];
		sphere2_center[k] = sphere2 -> center[k];
		sphere3_center[k] = sphere3 -> center[k];
	}
	setup_torus_fields (probe_radius, torus12, 
		torus12_center, &radius12, axis12);
	for (k = 0; k < 3; k++)
		axis13[k] = sphere3_center[k] - sphere1_center[k];
	distance13 = norm (axis13);
	if (distance13 <= 0.0) {
		ex = new_cept (GEOMETRY_ERROR,  DEGENERACY,  FATAL_SEVERITY);
		add_function (ex, "place_probe");
		add_source (ex, "msroll.c");
		add_message (ex, "coincident atoms");
		add_atom (ex, sphere1);
		add_atom (ex, sphere3);
	}
	for (k = 0; k < 3; k++)
		axis13[k] /= distance13;
	dot1213 = dot_product (axis12, axis13);
	if (dot1213 < -1.0) dot1213 = -1.0;
	else if (dot1213 > 1.0) dot1213 = 1.0;
	angle123 = acos (dot1213);
	cross (axis12, axis13, altitude);
	if (norm(altitude) <= 0.0) {
		/* colinear atoms */
		ex = new_cept (GEOMETRY_ERROR,  DEGENERACY,  FATAL_SEVERITY);
		add_function (ex, "place_probe");
		add_source (ex, "msroll.c");
		add_message (ex, "colinear atoms");
		add_atom (ex, sphere1);
		add_atom (ex, sphere2);
		add_atom (ex, sphere3);
		return (0);
	}
	sin_angle123 = sin (angle123);
	if (sin_angle123 <= 0.0) {
		/* check from Appendix I of Connolly's 1983 JAC article */
		special_distance = distance_squared (torus12_center, sphere3_center);
		special_radius = (sphere3 -> radius + probe_radius) *
				(sphere3 -> radius + probe_radius) - radius12 * radius12;
		if (special_distance < special_radius) {
			/* torus completely buried */
			/* next statement probably not necessary */
			torus12 -> buried = TRUE;
			torus12 -> free = FALSE;
			return (1);
		}
		else return (0);		/* far away */
	}
	for (k = 0; k < 3; k++)
		altitude[k] /= sin_angle123;
	cross (altitude, axis12, vector_torus_base);
	compute_torus_center (probe_radius, sphere1, sphere3, torus13_center);
	dot1213 = 0.0;
	for (k = 0; k < 3; k++)
		dot1213 += axis13[k] * (torus13_center[k] - torus12_center[k]);
	for (k = 0; k < 3; k++)
		base_point[k] = torus12_center[k] +
			vector_torus_base[k] * dot1213 / sin_angle123;
	height = sphere1 -> ers -
				distance_squared (base_point, sphere1 -> center);
	if (height < 0.0) {
		/* check from Appendix I
			of Connolly's 1983 JAC article */
		special_distance = distance_squared (torus12_center, sphere3_center);
		special_radius =
			(sphere3 -> radius + probe_radius) *
				(sphere3 -> radius + probe_radius) - radius12 * radius12;
		if (special_distance < special_radius) {
			/* torus completely buried */
			/* next assignment is probably not necessary */
			torus12 -> buried = TRUE;
			torus12 -> free = FALSE;
			return (1);
		}
		else return (0);	/* far away */
	}
	height = sqrt (height);

	if (torus13 == NULL) {
		ex = new_cept (LOGIC_ERROR,  NULL_VALUE,  FATAL_SEVERITY);
		add_function (ex, "place_probe");
		add_source (ex, "msroll.c");
		add_message (ex, "missing torus between atoms 1 and 3");
		return (0);
	}
	if (torus23 == NULL) {
		ex = new_cept (LOGIC_ERROR,  NULL_VALUE,  FATAL_SEVERITY);
		add_function (ex, "place_probe");
		add_source (ex, "msroll.c");
		add_message (ex, "missing torus between atoms 2 and 3");
		return (0);
	}
	torus12 -> free = torus23 -> free = torus13 -> free = FALSE;

	/* no probe placement triplication */
	if (sphere3 < sphere1 || sphere3 < sphere2) return (0);

	/* try above and below */
	for (side = -1; side <= 1; side += 2) {
		for (k = 0; k < 3; k++)
			probe_position[k] = base_point[k] +
				side * height * altitude[k];
		hit = check_collision (this_srf, probe_position, mutuals, last_mutual, sphere3);
        if(error()) return (0);
        if (hit) continue;
		/* good probe position, tori not buried */
		torus12 -> buried = FALSE;
		torus23 -> buried = FALSE;
		torus13 -> buried = FALSE;
		/* new probe */
		prb = new_probe (this_srf, probe_position, height,
			altitude, side, sphere1, sphere2, sphere3, torus12, torus23, torus13);
        if (error()) return (0);
	}	/* end of above-below loop*/
    return (0);
}

/* return 1 if collision */
int check_collision (struct surface *this_srf, double probe_position[3],
struct mutual *mutuals, struct mutual *last_mutual, struct sphere *sphere3) 
{
	int k, collide;
    int blocked, empty, where;
	double px, py, pz, dto4, dx, dy, dz;
	double epsilon;
	struct sphere *sphere4;
	struct mutual *neighbor4;
    struct cept *ex;

	blocked = 0; empty = 0;
	epsilon = EPSILON * EPSILON;
	if (this_srf -> use_grid) {
		/* sphere grid check */
		where = point_inquiry (this_srf -> sg, probe_position);
		if (error()) {
			if (tail_cept != NULL) {
				add_function (tail_cept, "check_collision");
				add_source (tail_cept, "msroll.c");
				add_message (tail_cept, "sphere_grid error");
			}
			this_srf -> use_grid = 0;
			free_sgrid (this_srf -> sg);
			return (0); /* die for now */
		}
		/* if where is 1, the cube is partially blocked */
		if (where == 0) empty = 1;
		else if (where == 2) blocked = 1;
	}
	if (blocked) return (1);
	if (!empty) {
		/* collision check */
		collide = FALSE;
		px = probe_position[0];
		py = probe_position[1];
		pz = probe_position[2];
		for (neighbor4 = mutuals; neighbor4 <= last_mutual;
			neighbor4++) {
			sphere4 = neighbor4 -> sphptr;
			if (sphere4 == sphere3) continue;
			dx = px - sphere4 -> center[0];
			dy = py - sphere4 -> center[1];
			dz = pz - sphere4 -> center[2];
			dto4 = dx * dx + dy * dy + dz * dz;
			/* just testing */
			if (dto4 < sphere4 -> ers - epsilon) {
				collide = TRUE;
				break;
			}
		}
		if (collide) return (1);
	}
    return (0);
}

void create_concave_faces (struct surface *this_srf)
{
    struct probe *prb;

    for (prb = this_srf -> head_probe; prb != NULL; prb = prb -> next) {
        /* new concave face */
        new_concave_face (this_srf, prb);
	}
}

void mark_not_buried (struct surface *this_srf)
{
	int nt;
	struct sphere *sphere_ptr;
	struct neighbor *nbr, *first_neighbor, *last_neighbor;

	for (sphere_ptr = (struct sphere *) (this_srf -> head_atom); sphere_ptr != NULL;
		sphere_ptr = sphere_ptr -> next) {
		nt = 0; /* count tori involving this sphere */
		first_neighbor = sphere_ptr -> first_neighbor;
		last_neighbor = sphere_ptr -> last_neighbor;
		if (first_neighbor != NULL)
			for (nbr = first_neighbor; nbr <= last_neighbor; nbr++)
				if (nbr -> torptr != NULL) nt++;
		/* if there are no tori, the sphere is not buried */
		if (nt <= 0)
			sphere_ptr -> buried = FALSE;
	}

}

void free_nbr (struct surface *this_srf)
{
	struct sphere *sphere_ptr;
	
	for (sphere_ptr = (struct sphere *) (this_srf -> head_atom); sphere_ptr != NULL;
		sphere_ptr = sphere_ptr -> next) {
		if (sphere_ptr -> first_neighbor != NULL) {
			free_objects (NEIGHBOR, (short *) (sphere_ptr -> first_neighbor));
			sphere_ptr -> first_neighbor = sphere_ptr -> last_neighbor = NULL;
		}
	}
}

/* LOW LEVEL SUBROUTINES: */

/* return true if atoms sph1 and sph2 are neighbors */
int are_neighbors (struct surface *this_srf, struct sphere *sph1, struct sphere *sph2)
{
	int k;
	double radius_sum, radius_sum_squared;
	double distance12, distance12_squared;
	double *sph1_center, *sph2_center;
	char message[MAXLINE];
    struct cept *ex;

	sph1_center = sph1 -> center;
	sph2_center = sph2 -> center;
	radius_sum = sph1 -> radius + sph2 -> radius + 2 * this_srf -> probe_radius;
	radius_sum_squared = radius_sum * radius_sum;

	distance12_squared = 0.0;
	for (k = 0; k < 3; k++) {
		distance12 = *sph2_center++ - *sph1_center++;
		if (distance12 >= radius_sum) return (0);
		if (distance12 <=  (-radius_sum)) return (0);
		distance12_squared += distance12 * distance12;
	}
	if (sph1 != sph2 && distance12_squared <= 0.0) {
		ex = new_cept (GEOMETRY_ERROR,  DEGENERACY,  FATAL_SEVERITY);
		add_function (ex, "are_neighbors");
		add_source (ex, "msroll.c");
		add_message (ex, "coincident atoms");
		add_atom (ex, sph1);
		add_atom (ex, sph2);
		return(0);
	}
	/* check whether atoms can be bridged by probe sphere */
	return ((int) (distance12_squared < radius_sum_squared));
}

/* form a list of mutual neighbors of sphere1 and sphere2 */
int form_mutual_list (struct sphere *sphere1, struct sphere *sphere2, struct mutual *mutuals)
{
	int n_mutual;
	struct sphere *sphere1_ptr, *sphere2_ptr;
	struct neighbor *sphere1_list, *sphere1_last;
	struct neighbor *sphere2_list, *sphere2_last;
	struct neighbor *sphere1_hdl, *sphere2_hdl;
	struct mutual *mutual_ptr;

	sphere1_list = sphere1 -> first_neighbor;
	sphere2_list = sphere2 -> first_neighbor;
	sphere1_last = sphere1 -> last_neighbor;
	sphere2_last = sphere2 -> last_neighbor;
	if (sphere1_list == NULL || sphere2_list == NULL) return (0);
	mutual_ptr = mutuals;
	sphere1_hdl = sphere1_list;
	sphere2_hdl = sphere2_list;
	/* only one loop required because neighbors are in atom-number order */
	while (sphere1_hdl <= sphere1_last &&  sphere2_hdl <= sphere2_last) {
		sphere1_ptr = sphere1_hdl -> sphptr;
		sphere2_ptr = sphere2_hdl -> sphptr;
		if (sphere1_ptr -> number < sphere2_ptr -> number) {
			sphere1_hdl++;
			continue;
		}
		else if (sphere1_ptr -> number > sphere2_ptr -> number) {
			sphere2_hdl++;
			continue;
		}
		/* we have a mutual neighbor */
		mutual_ptr -> sphptr = sphere1_ptr;
		mutual_ptr -> tor13 = sphere1_hdl -> torptr;
		mutual_ptr -> tor23 = sphere2_hdl -> torptr;
		mutual_ptr++;
		sphere1_hdl++;
		sphere2_hdl++;
	}
	n_mutual = mutual_ptr - mutuals;
	return (n_mutual);
}

/* low level tori routines: */


void finish_up (struct surface *this_srf)
{
	int n_arc, iarc;
	struct arc *arcptr;
	struct edge *edg, *prev_edg;
	struct cycle *cyc1, *cyc2, *cyc;
	struct face *fac;
	struct probe *prb;
	char message[MAX_STRING];
    struct torus *tor;
    struct cept *ex;


	informd ("finish_up");
	for (fac = this_srf -> head_face; fac != NULL; fac = fac -> next) {
		switch ((int) fac -> shape) {
		case CONVEX:
			/* already done */
			break;
		case SADDLE:
			/* form edges and cycles for saddle face */
			n_arc = fac -> n_arc;
			cyc1 = new_cycle (fac, NULL);
			if (cyc1 == NULL) return;
			this_srf -> n_cycle++;
			prev_edg = NULL;
			for (iarc = n_arc - 1; iarc >= 0; iarc--) {
				arcptr = fac -> arcsp[iarc];
				if (arcptr == NULL) continue;
				if (arcptr -> eaten > 0) {
                    tor = fac -> ptr.tor;
					ex = new_cept (GEOMETRY_ERROR,  INCONSISTENCY,  FATAL_SEVERITY);
					add_function (ex, "finish_up");
					add_source (ex, "msroll.c");
					add_message(ex, "eaten arc in saddle face");
					add_atom (ex, tor -> atm[0]);
					add_atom (ex, tor -> atm[1]);
					return;
				}
				edg = new_edge (arcptr, 1, fac, NULL);
				if (edg == NULL) {
					inform("null edge pointer returned from new_edge");
					return;
				}
				if (prev_edg == NULL) cyc1 -> first_edge = edg;
				else prev_edg -> next = edg;
				prev_edg = edg;
				if (n_arc == 2) break;	/* hoop has two cycles */
			}
			/* 2nd cycle for hoop */
			if (n_arc == 2) {
				arcptr = fac -> arcsp[0];
				if (arcptr == NULL) {
					ex = new_cept (POINTER_ERROR,  NULL_VALUE,  FATAL_SEVERITY);
					add_function (ex, "finish_up");
					add_source (ex, "msroll.c");
                    tor = fac -> ptr.tor;
					add_atom (ex, tor -> atm[0]);
					add_atom (ex, tor -> atm[1]);
					return;
				}
				if (arcptr -> eaten > 0) {
					ex = new_cept (GEOMETRY_ERROR,  INCONSISTENCY,  FATAL_SEVERITY);
					add_function (ex, "finish_up");
					add_source (ex, "msroll.c");
					add_message(ex, "eaten arc in hoop face");
                    tor = fac -> ptr.tor;
					add_atom (ex, tor -> atm[0]);
					add_atom (ex, tor -> atm[1]);
					return;
				}
				edg = new_edge (arcptr, 1, fac, NULL);
				if (edg == NULL) return;
				cyc2 = new_cycle (fac, edg);
				if (cyc2 == NULL) return;
				this_srf -> n_cycle++;
			}
			break;
		case CONCAVE:
			prb = fac -> ptr.prb;
			if (prb -> low) break;	 /* already done in cusp_trimming */
			/* form edges and cycles for concave triangle */
			cyc = new_cycle (fac, NULL);
			if (cyc == NULL) return;
			this_srf -> n_cycle++;
			prev_edg = NULL;
			for (iarc = 0; iarc < MAXPA; iarc++) {
				arcptr = fac -> arcsp[iarc];
				if (arcptr == NULL) continue;
				if (arcptr -> eaten > 0) {
					ex = new_cept (GEOMETRY_ERROR,  INCONSISTENCY,  FATAL_SEVERITY);
					add_function (ex, "finish_up");
					add_source (ex, "msroll.c");
					add_message(ex, "eaten arc in concave face");
					add_atom (ex, prb -> atm[0]);
					add_atom (ex, prb -> atm[1]);
					add_atom (ex, prb -> atm[2]);
					if (prb -> atm[3] != NULL)
						add_atom(ex, prb -> atm[3]);
					return;
				}
				edg = new_edge (arcptr, 0, fac, NULL);
				if (edg == NULL) return;
				if (prev_edg == NULL) cyc -> first_edge = edg;
				else prev_edg -> next = edg;
				prev_edg = edg;
			}
			break;
		default:
			ex = new_cept (GEOMETRY_ERROR,  INCONSISTENCY,  FATAL_SEVERITY);
			add_function (ex, "finish_up");
			add_source (ex, "msroll.c");
			add_message(ex, "invalid face shape");
			return;
		}
	}
}

/* count the number of edges */
void count_edges (struct surface *this_srf)
{
	char message[MAXLINE];
	struct face *fac;

	informd ("count_edges");
	this_srf -> n_edge = 0;
	for (fac = this_srf -> head_face; fac != NULL; fac = fac -> next)
		this_srf -> n_edge += fac -> n_arc;

	if (this_srf -> n_edge % 2 != 0) {
		inform ("(count_edges): odd number of edges");
		return;
	}
	if (this_srf -> n_edge / 2 != this_srf -> n_arc) {
		sprintf (message,
		"(count_edges): arc-edg inconsistency: n_edge = %6ld; n_arc = %6ld",
			this_srf -> n_edge, this_srf -> n_arc);
		inform(message); return;
	}
}

void count_problem_faces2 (struct surface *this_srf)
{
	struct face *f;
	struct sphere *a;

	informd ("count_problem_faces2");
	this_srf -> n_problem_face = 0;
	for (f = this_srf -> head_face; f != NULL; f = f -> next)
		if (f -> problem) this_srf -> n_problem_face++;

	this_srf -> n_problem_atom = 0;
	for (a = this_srf -> head_atom; a != NULL; a = a -> next)
		if (a -> problem) this_srf -> n_problem_atom++;
}

/* number the arcs */
void number_arcs (struct surface *this_srf)
{
	int iarc, m;
	struct arc *a;
	struct face *fac;
	struct cycle *cyc;
	struct edge *edg;
	struct edge *e;
    struct cept *ex;

	iarc = 0;

	informd ("number_arcs");
	for (fac = this_srf -> head_face; fac != NULL; fac = fac -> next) {
		if (fac -> problem) continue;
		switch ((int) fac -> shape) {
		case CONVEX:
			if (fac -> n_cycle <= 0) break;
			for (cyc = fac -> first_cycle; cyc != NULL; cyc = cyc -> next)
				for (e = cyc -> first_edge; e != NULL; e = e -> next)
					e -> arcptr -> number = ++iarc;
			break;
		case CONCAVE:
			if (fac -> ptr.prb -> low)
				for (cyc = fac -> first_cycle; cyc != NULL;
					cyc = cyc -> next)
					for (edg = cyc -> first_edge; edg != NULL;
						edg = edg -> next) {
						if (edg -> orn) continue;
						a = edg -> arcptr;
						a -> number = ++iarc;
					}
			else
				for (m = 0; m < fac -> n_arc; m++) {
					a = fac -> arcsp[m];
                    if (a == NULL) continue;
					a -> number = ++iarc;
				}	/* end of arc loop */
			break;
		default:
			break;
		}
	}

	if (iarc > this_srf -> n_arc) {
		ex = new_cept (LOGIC_ERROR,  MSOVERFLOW,  FATAL_SEVERITY);
        add_function (ex, "number_arcs");
		add_source (ex, "msroll.c");
		add_message (ex, "arc count inconsistency");
		return;
	}
	if (iarc < this_srf -> n_arc) this_srf -> n_arc = iarc;
}

/* clean circles */
void clean_circles (struct surface *this_srf)
{
	int m, icir;
	long ncleaned = 0;
	struct circle *cir, *pcir, *nxtcir;
	struct face *fac;
	struct cycle *cyc;
	struct edge *edg;
	struct edge *e;
    struct cept *ex;
	char message[MAXLINE];

	informd ("clean_circles");
	/* initialize to unused in any arc */
	for (cir = this_srf -> head_circle; cir != NULL; cir = cir -> next)
		cir -> unused = 1;

	/* mark circles used by arcs */
	for (fac = this_srf -> head_face; fac != NULL; fac = fac -> next) {
		if (fac -> problem) continue;
		for (cyc = fac -> first_cycle; cyc != NULL; cyc = cyc -> next)
			for (e = cyc -> first_edge; e != NULL; e = e -> next)
				e -> arcptr -> cir -> unused = 0;
	}

	/* unlink and return unused circles */

	pcir = NULL;
	icir = 0;
	nxtcir = NULL;

	for (cir = this_srf -> head_circle; cir != NULL; cir = nxtcir) {
		nxtcir = cir -> next;
		if (cir -> unused) {
			if (pcir == NULL) this_srf -> head_circle = nxtcir;
			else pcir -> next = nxtcir;
			if (cir == this_srf -> tail_circle)
				this_srf -> tail_circle = pcir;
			free_circle (cir);
			ncleaned++;
		}
		else {
			pcir = cir;
			cir -> number = ++icir;
		}
	}

	if (icir > this_srf -> n_circle) {
		ex = new_cept (LOGIC_ERROR,  MSOVERFLOW,  FATAL_SEVERITY);
        add_function (ex, "clean_circles");
		add_source (ex, "msroll.c");
		add_message (ex, "circle count inconsistency");
		return;
	}
	if (icir < this_srf -> n_circle) this_srf -> n_circle = icir;
	sprintf(message, "%8ld unused circles returned", ncleaned);
	informd (message);
}

void roll_mem ()
{

	unsigned long size;
	int type;
	char type_name[32];
	
	init_mem ((unsigned long) N_OBJECTS);

	type = MSSCENE;
	size = sizeof (struct msscene);
	strcpy (type_name, "msscene");
	define_type (type, size, type_name);
	if (error()) return;
	
	type = SPHERE;
	size = sizeof (struct sphere);
	strcpy (type_name, "sphere");
	define_type (type, size, type_name);
	if (error()) return;
	
	type = TORUS;
	size = sizeof (struct torus);
	strcpy (type_name, "torus");
	define_type (type, size, type_name);
	if (error()) return;
	
	type = PROBE;
	size = sizeof (struct probe);
	strcpy (type_name, "probe");
	define_type (type, size, type_name);
	if (error()) return;
	
	type = CENTRAL;
	size = sizeof (struct central);
	strcpy (type_name, "central");
	define_type (type, size, type_name);
	if (error()) return;
	
	type = VERTEX;
	size = sizeof (struct vertex);
	strcpy (type_name, "vertex");
	define_type (type, size, type_name);
	if (error()) return;
	
	type = CIRCLE;
	size = sizeof (struct circle);
	strcpy (type_name, "circle");
	define_type (type, size, type_name);
	if (error()) return;
	
	type = ARC;
	size = sizeof (struct arc);
	strcpy (type_name, "arc");
	define_type (type, size, type_name);
	if (error()) return;
	
	type = FACE;
	size = sizeof (struct face);
	strcpy (type_name, "face");
	define_type (type, size, type_name);
	if (error()) return;
	
	type = CYCLE;
	size = sizeof (struct cycle);
	strcpy (type_name, "cycle");
	define_type (type, size, type_name);
	if (error()) return;
	
	type = EDGE;
	size = sizeof (struct edge);
	strcpy (type_name, "edge");
	define_type (type, size, type_name);
	if (error()) return;
	
	type = CUSP;
	size = sizeof (struct cusp);
	strcpy (type_name, "cusp");
	define_type (type, size, type_name);
	if (error()) return;
	
	type = CUSP_LINK;
	size = sizeof (struct cusp_link);
	strcpy (type_name, "cusp_link");
	define_type (type, size, type_name);
	if (error()) return;
	
	type = COMPONENT;
	size = sizeof (struct component);
	strcpy (type_name, "component");
	define_type (type, size, type_name);
	if (error()) return;
	
	type = CHUNK;
	size = sizeof (struct chunk);
	strcpy (type_name, "chunk");
	define_type (type, size, type_name);
	if (error()) return;

	type = VARIETY;
	size = sizeof (struct variety);
	strcpy (type_name, "variety");
	define_type (type, size, type_name);
	if (error()) return;

	type = PHNVTX;
	size = phnvtx_size ();
	strcpy (type_name, "phnvtx");
	define_type (type, size, type_name);
	if (error()) return;

	type = PHNEDG;
	size = phnedg_size ();
	strcpy (type_name, "phnedg");
	define_type (type, size, type_name);
	if (error()) return;

	type = PHNTRI;
	size = phntri_size ();
	strcpy (type_name, "phntri");
	define_type (type, size, type_name);
	if (error()) return;

	type = POLYGON;
	size = sizeof (struct polygon);
	strcpy (type_name, "polygon");
	define_type (type, size, type_name);
	if (error()) return;

	type = HEDVTX;
	size = sizeof (struct hedvtx);
	strcpy (type_name, "hedvtx");
	define_type (type, size, type_name);
	if (error()) return;

	type = HEDEDG;
	size = sizeof (struct hededg);
	strcpy (type_name, "hededg");
	define_type (type, size, type_name);
	if (error()) return;

	type = HEDTRI;
	size = sizeof (struct hedtri);
	strcpy (type_name, "hedtri");
	define_type (type, size, type_name);
	if (error()) return;

	type = VTXGRP;
	size = sizeof (struct vtxgrp);
	strcpy (type_name, "vtxgrp");
	define_type (type, size, type_name);
	if (error()) return;

	type = HEDRON;
	size = sizeof (struct hedron);
	strcpy (type_name, "hedron");
	define_type (type, size, type_name);
	if (error()) return;

	type = OBJECT_SCHEME;
	size = sizeof (struct object_scheme);
	strcpy (type_name, "object_scheme");
	define_type (type, size, type_name);
	if (error()) return;

	type = COLOR_RAMP;
	size = sizeof (struct color_ramp);
	strcpy (type_name, "color_ramp");
	define_type (type, size, type_name);
	if (error()) return;

	type = MOLECULE;
	size = sizeof (struct molecule);
	strcpy (type_name, "molecule");
	define_type (type, size, type_name);
	if (error()) return;

	type = SURFACE;
	size = sizeof (struct surface);
	strcpy (type_name, "surface");
	define_type (type, size, type_name);
	if (error()) return;

	type = NEIGHBOR;
	size = sizeof (struct neighbor);
	strcpy (type_name, "neighbor");
	define_type (type, size, type_name);
	if (error()) return;

	type = MUTUAL;
	size = sizeof (struct mutual);
	strcpy (type_name, "mutual");
	define_type (type, size, type_name);
	if (error()) return;

	type = PAIR;
	size = sizeof (struct pair);
	strcpy (type_name, "pair");
	define_type (type, size, type_name);
	if (error()) return;

	type = SCUBE;
	size = sizeof (struct scube);
	strcpy (type_name, "scube");
	define_type (type, size, type_name);
	if (error()) return;
	
	type = SNUMBER;
	size = sizeof (struct snumber);
	strcpy (type_name, "snumber");
	define_type (type, size, type_name);
	if (error()) return;
	
	type = SNUMBERBLOCK;
	size = sizeof (struct sNumberBlock);
	strcpy (type_name, "sNumberBlock");
	define_type (type, size, type_name);
	if (error()) return;
	
	type = SPHEREGRID;
	size = sizeof (struct spheregrid);
	strcpy (type_name, "spheregrid");
	define_type (type, size, type_name);
	if (error()) return;
	
	type = PATTERNRECORD;
	size = sizeof (struct patternRecord);
	strcpy (type_name, "patternRecord");
	define_type (type, size, type_name);
	if (error()) return;
	
	type = RADIUSRECORD;
	size = sizeof (struct radiusRecord);
	strcpy (type_name, "radiusRecord");
	define_type (type, size, type_name);
	if (error()) return;
	
	type = RESIDUEINDEX;
	size = sizeof (struct residueIndex);
	strcpy (type_name, "residueIndex");
	define_type (type, size, type_name);
	if (error()) return;
	
	type = VERTEX_PAIR;
	size = sizeof (struct vertex_pair);
	strcpy (type_name, "vertex_pair");
	define_type (type, size, type_name);
	if (error()) return;
	
	type = CRITLINK;
	size = sizeof (struct critlink);
	strcpy (type_name, "critlink");
	define_type (type, size, type_name);
	if (error()) return;
	
	type = CUSP_EXTENSION;
	size = sizeof (struct cusp_extension);
	strcpy (type_name, "cusp_extension");
	define_type (type, size, type_name);
	if (error()) return;
	
	type = ARRAY;
	size = sizeof (struct array);
	strcpy (type_name, "array");
	define_type (type, size, type_name);
	if (error()) return;
	
	type = ATOM;
	size = sizeof (struct atom);
	strcpy (type_name, "atom");
	define_type (type, size, type_name);
	if (error()) return;
	
	type = BOND;
	size = sizeof (struct bond);
	strcpy (type_name, "bond");
	define_type (type, size, type_name);
	if (error()) return;
	
	type = OBJECT_HEADER;
	size = sizeof (struct object_header);
	strcpy (type_name, "object_header");
	define_type (type, size, type_name);
	if (error()) return;
	
	type = OBJECT_BLOCK;
	size = sizeof (struct object_block);
	strcpy (type_name, "object_block");
	define_type (type, size, type_name);
	if (error()) return;
	
	type = ELEMENT;
	size = sizeof (struct element);
	strcpy (type_name, "element");
	define_type (type, size, type_name);
	if (error()) return;
	
	type = ELEMENT_BLOCK;
	size = sizeof (struct element_block);
	strcpy (type_name, "element_block");
	define_type (type, size, type_name);
	if (error()) return;
	
	type = SET;
	size = sizeof (struct set);
	strcpy (type_name, "set");
	define_type (type, size, type_name);
	if (error()) return;
	
	type = SYMBOL;
	size = sizeof (struct symbol);
	strcpy (type_name, "symbol");
	define_type (type, size, type_name);
	if (error()) return;
	
	type = BOOLEAN;
	size = sizeof (struct boolean);
	strcpy (type_name, "boolean");
	define_type (type, size, type_name);
	if (error()) return;
	
	type = INTEGER;
	size = sizeof (struct integer);
	strcpy (type_name, "integer");
	define_type (type, size, type_name);
	if (error()) return;
	
	type = REAL;
	size = sizeof (struct real);
	strcpy (type_name, "real");
	define_type (type, size, type_name);
	if (error()) return;
	
	type = STRING;
	size = sizeof (struct string);
	strcpy (type_name, "string");
	define_type (type, size, type_name);
	if (error()) return;
	
	type = REGION;
	size = sizeof (struct region);
	strcpy (type_name, "region");
	define_type (type, size, type_name);
	if (error()) return;
	
	type = DSDESC;
	size = sizeof (struct dsdesc);
	strcpy (type_name, "dsdesc");
	define_type (type, size, type_name);
	if (error()) return;
	
	type = HIDDEN;
	size = sizeof (struct hidden);
	strcpy (type_name, "hidden");
	define_type (type, size, type_name);
	if (error()) return;
	
	type = BOX;
	size = sizeof (struct box);
	strcpy (type_name, "box");
	define_type (type, size, type_name);
	if (error()) return;
	
	type = RUN;
	size = sizeof (struct run);
	strcpy (type_name, "run");
	define_type (type, size, type_name);
	if (error()) return;
	
	type = TLIST;
	size = sizeof (struct tlist);
	strcpy (type_name, "tlist");
	define_type (type, size, type_name);
	if (error()) return;
	
	type = PLIST;
	size = sizeof (struct plist);
	strcpy (type_name, "plist");
	define_type (type, size, type_name);
	if (error()) return;
	
	type = EDG;
	size = sizeof (struct edg);
	strcpy (type_name, "edg");
	define_type (type, size, type_name);
	if (error()) return;
	
	type = ENDPNT;
	size = sizeof (struct endpnt);
	strcpy (type_name, "endpnt");
	define_type (type, size, type_name);
	if (error()) return;
	
	type = MIDPLN;
	size = sizeof (struct midpln);
	strcpy (type_name, "midpln");
	define_type (type, size, type_name);
	if (error()) return;
	
	type = CLUSTER;
	size = sizeof (struct cluster);
	strcpy (type_name, "cluster");
	define_type (type, size, type_name);
	if (error()) return;
	
	type = CEPT;
	size = sizeof (struct cept);
	strcpy (type_name, "cept");
	define_type (type, size, type_name);
	if (error()) return;
	
	type = PHOLDER;
	size = sizeof (struct pholder);
	strcpy (type_name, "pholder");
	define_type (type, size, type_name);
	if (error()) return;
	
	type = QUARTET;
	size = sizeof (struct quartet);
	strcpy (type_name, "quartet");
	define_type (type, size, type_name);
	if (error()) return;
	
}

