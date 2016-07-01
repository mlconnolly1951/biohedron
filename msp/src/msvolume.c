/*
	Molecular Surface Package
	Copyright 1986, 1989 by Michael L. Connolly
	All rights reserved

	December 16, 2001
*/

/* PART V: AREAS AND VOLUMES */

#include "ms.h"
#include "msenum.h"
#include "mslimit.h"
#include "msdefault.h"
#include "msstruct.h"
#include "msglobal.h"
#include "msfunction.h"

/* calculate solvent-excluded volume and also areas */


int compute_volume (struct surface *this_srf)
{
	int comp, k, atom_number;
	double partial_area, accessible_area, area_ratio;
	double radius, expanded_radius, radius_ratio;
	struct sphere *atm;
	struct face *fac;
	struct component *cmp_ptr;
    struct cept *ex;

	informd("compute_volume");
	if (! this_srf -> surface_completed) {
		inform ("premature volume command");
		return(0);
	}
	/* initialization */


	this_srf -> molecule_area = this_srf -> total_contact_area = this_srf -> total_reentrant_area = 0.0;
	this_srf -> total_accessible_area = 0.0;
	this_srf -> mol_vol = 0.0;

	/* compute surface center */
	for (k = 0; k < 3; k++)
		this_srf -> center[k] = 0.0;
	for (atm = (this_srf -> head_atom);
		atm != NULL; atm = (atm -> next))
		for (k = 0; k < 3; k++)
			this_srf -> center[k] += atm -> center[k];
	for (k = 0; k < 3; k++)
		this_srf -> center[k] /= this_srf -> n_atom;

	/* accumulate various area and volume sums */
	for (fac = this_srf -> head_face; fac != NULL; fac = fac -> next) {
		if (fac -> problem) {
			inform ("problem face skipped in area/volume computation");
			continue;
		}
		accessible_area = 0.0;
		comp = fac -> comp;
		cmp_ptr = get_component_ptr (this_srf, comp);
		switch ((int) fac -> shape) {
		case CONVEX:
			atm = fac -> ptr.atm;
			radius = atm -> radius;
			expanded_radius = radius + this_srf -> probe_radius;
			radius_ratio = expanded_radius / radius;
			area_ratio = radius_ratio * radius_ratio;
			convex_face_area (fac);
			accessible_area = fac -> area * area_ratio;
			this_srf -> total_accessible_area += accessible_area;
			this_srf -> total_contact_area += fac -> area;
			add_chunk ((struct atom *) (fac -> ptr.atm), (int) fac -> comp,
				 (double) fac -> area, (double) 0.0,
				 accessible_area);
			convex_piece_volume (fac);
			if (error()) return(0);
			break;
		case SADDLE:
			saddle_face_area (fac);
			this_srf -> total_reentrant_area += fac -> area;
			for (k = 0; k < 2; k++) {
				atom_number = fac -> ptr.tor -> atm[k] -> number;
				partial_area = aliquot (fac, atom_number);
				add_chunk ((struct atom *) fac -> ptr.tor -> atm[k],
					(int) fac -> comp, (double) 0.0,
					partial_area, (double) 0.0);
			}
			saddle_piece_volume (fac);
			if (error()) return(0);
			break;
		case CONCAVE:
			concave_face_area (fac);
			this_srf -> total_reentrant_area += fac -> area;
			for (k = 0; k < MAXPA; k++) {
				atm = fac -> ptr.prb -> atm[k];
				if (atm == NULL) continue;
				atom_number = atm -> number;
				partial_area = aliquot (fac, atom_number);
				add_chunk ((struct atom *) fac -> ptr.prb -> atm[k],
					(int) fac -> comp, (double) 0.0,
					partial_area, (double) 0.0);
			}
			if (error()) return(0);
			break;
		default:
			ex = new_cept (ENUM_ERROR,  INVALID_VALUE,  FATAL_SEVERITY);
			add_function (ex, "compute_volume");
			add_source (ex, "msvolume.c");
            add_long (ex, "fac -> shape", (long) fac -> shape);
			return(0);
		}
		this_srf -> molecule_area += fac -> area;
		cmp_ptr -> area  += fac -> area;
		cmp_ptr -> accessible += accessible_area;
		face_volume (fac);
		this_srf -> mol_vol += fac -> vol2;
		cmp_ptr -> volume += fac -> vol2;
	}

	setup_second_links (this_srf);
	if (error()) return(0);
	informd("         volume computed");
	return (1);
}

int compute_old_volume (struct surface *this_srf)
{
	double prbvol;
	char message[MAXLINE];
	struct face *fac;
	struct probe *prb;
    struct cept *ex;

	if (! this_srf -> surface_completed) {
		inform ("premature old_volume command");
		return(0);
	}
	if (this_srf -> old_volume_computed) {
		sprintf (message,
			"volume (old method)  =   %10.3f", this_srf -> old_volume);
		inform(message);
		return(1);
	}
	/* initialization */


	this_srf -> surface_volume = 0.0;

	/* accumulate various area and volume sums */
	for (fac = this_srf -> head_face; fac != NULL; fac = fac -> next) {
		if (fac -> problem) {
			inform ("problem face skipped in area/volume computation");
			continue;
		}
		switch ((int) fac -> shape) {
		case CONVEX:
			convex_piece_volume (fac);
			if (error()) return(0);
			this_srf -> surface_volume += fac -> vol1;
			break;
		case SADDLE:
			saddle_piece_volume (fac);
			if (error()) return(0);
			this_srf -> surface_volume += fac -> vol1;
			break;
		case CONCAVE:
			break;
		default:
			ex = new_cept (ENUM_ERROR,  INVALID_VALUE,  FATAL_SEVERITY);
			add_function (ex, "compute_old_volume");
			add_source (ex, "msvolume.c");
            add_long (ex, "face shape", (long) fac -> shape);
			return(0);
		}
	}

	for (prb = this_srf -> head_probe; prb != NULL; prb = prb -> next) {
		prbvol = probe_volume (prb);
		if (error()) return(0);
		this_srf -> surface_volume += prbvol;
	}

	this_srf -> polyhedron_volume = hedron (this_srf);
	this_srf -> old_volume = this_srf -> surface_volume + this_srf -> polyhedron_volume;
	
	compute_old_area (this_srf);	/* numerical areas for low probes */
	if (error()) return(0);
	this_srf -> old_volume_computed = 1;
	sprintf (message,
		"volume (old method)  =   %10.3f", this_srf -> old_volume);
	inform(message);
	return(1);
}


/* AREA routines: */

void convex_face_area (struct face *fac) 	/* calculate area of convex face */
{
	int problem;
	double omega;
	struct sphere *atm;
	struct surface *this_srf;

	atm = fac -> ptr.atm;
	this_srf = fac -> srf;
	omega = compute_omega (this_srf -> probe_radius, fac);
	problem = check_omega (fac, omega);
	if (problem) {
		fac -> problem = TRUE;
	}

	if (omega > 4 * PI) omega = 4 * PI;
	if (omega < 0.0) omega = 0.0;

	fac -> area = atm -> radius * atm -> radius * omega;
}

void saddle_face_area (struct face *fac)			/* calculate area of saddle face */
{
	int n_fac_arc, arc_index;
	double ratio, phi, atm1_theta, atm2_theta;
	double area, cusp_angle, cusp_area;
	struct arc *arc1, *arc2, *ark;
	struct circle *cir1, *cir2;
	struct torus *tor;
	struct surface *this_srf;
    struct cept *ex;

	this_srf = fac -> srf;
	if (this_srf -> probe_radius <= 0.0) return;

	n_fac_arc = fac -> n_arc;

	if (n_fac_arc == 1) {		/* conical saddle face */
		arc1 = fac -> arcsp[0];
		if (arc1 == NULL) {
			ex = new_cept (POINTER_ERROR,  NULL_VALUE,  FATAL_SEVERITY);
			add_object (ex, SURFACE, "this_srf");
			add_function (ex, "saddle_face_area");
			add_source (ex, "msvolume.c");
			return;
		}
		arc2 = NULL;
	}
	else if (n_fac_arc == 2) {	/* hoop */
		arc2 = fac -> arcsp[0];
		if (arc2 == NULL) {
			ex = new_cept (POINTER_ERROR,  NULL_VALUE,  FATAL_SEVERITY);
			add_object (ex, SURFACE, "this_srf");
			add_function (ex, "saddle_face_area");
			add_source (ex, "msvolume.c");
			return;
		}
		arc1 = fac -> arcsp[1];
		if (arc1 == NULL) {
			ex = new_cept (POINTER_ERROR,  NULL_VALUE,  FATAL_SEVERITY);
			add_object (ex, SURFACE, "this_srf");
			add_function (ex, "saddle_face_area");
			add_source (ex, "msvolume.c");
			return;
		}
	}
	else if (n_fac_arc == 3) {	/* conical */
		arc1 = NULL;
		for (arc_index = 0; arc_index < n_fac_arc; arc_index++) {
            ark = fac -> arcsp[arc_index];
			if (ark == NULL) {
				ex = new_cept (POINTER_ERROR,  NULL_VALUE,  FATAL_SEVERITY);
				add_object (ex, SURFACE, "this_srf");
				add_function (ex, "saddle_face_area");
				add_source (ex, "msvolume.c");
				return;
			}
			if (ark -> shape == CONVEX) {
				arc1 = ark;
				break;
			}
		}
		if (arc1 == NULL) {
			ex = new_cept (POINTER_ERROR,  NULL_VALUE,  FATAL_SEVERITY);
			add_function (ex, "saddle_face_area");
			add_source (ex, "msvolume.c");
            add_long (ex, "number of arcs in face", (long) n_fac_arc);
			add_message (ex, "bad saddle cone");
			return;
		}
		arc2 = NULL;
	}
	else if (n_fac_arc == 4) {		/* regular saddle face */
		arc2 = fac -> arcsp[1];
		if (arc2 == NULL) {
			ex = new_cept (POINTER_ERROR,  NULL_VALUE,  FATAL_SEVERITY);
			add_object (ex, SURFACE, "this_srf");
			add_function (ex, "saddle_face_area");
			add_source (ex, "msvolume.c");
			return;
		}
		arc1 = fac -> arcsp[3];
		if (arc1 == NULL) {
			ex = new_cept (POINTER_ERROR,  NULL_VALUE,  FATAL_SEVERITY);
			add_object (ex, SURFACE, "this_srf");
			add_function (ex, "saddle_face_area");
			add_source (ex, "msvolume.c");
			return;
		}
	}
	else {
		ex = new_cept (GEOMETRY_ERROR,  INVALID_VALUE,  FATAL_SEVERITY);
		add_function (ex, "saddle_face_area");
		add_source (ex, "msvolume.c");
        add_long (ex, "number of arcs in face", (long) n_fac_arc);
		return;
	}

	phi = arc1 -> phi;
	cir1 = arc1 -> cir;
	atm1_theta = cir1 -> theta;

	if (arc2 != NULL) {
		cir2 = arc2 -> cir;
		atm2_theta = cir2 -> theta;
	}
	else atm2_theta = 0.0;

	tor = fac -> ptr.tor;
	area = tor -> radius * this_srf -> probe_radius * (atm1_theta + atm2_theta);
	area = area - this_srf -> probe_radius * this_srf -> probe_radius *
		(sin (atm1_theta) + sin (atm2_theta));

	/* cusp correction */

	if (arc2 == NULL) {
		ratio = tor -> radius / this_srf -> probe_radius;
		if (ratio > 1.0) ratio = 1.0;
		if (ratio < -1.0) ratio = -1.0;
		cusp_angle = acos (ratio);
		cusp_area = tor -> radius * this_srf -> probe_radius * cusp_angle
			- this_srf -> probe_radius * this_srf -> probe_radius * sin (cusp_angle);
		area -= cusp_area;
	}

	area *= phi;		/* multiply by saddle wrap angle */

	fac -> area = area;
}

/* calculate area of negatively curved face */
void concave_face_area (struct face *fac)
{
	int problem;
	double omega;
	struct surface *this_srf;

	this_srf = fac -> srf;
	if (this_srf -> probe_radius <= 0.0) return;

	omega = compute_omega (this_srf -> probe_radius, fac);
	problem = check_omega (fac, omega);
	if (problem) {
		fac -> problem = TRUE;
	}

	if (omega > 4 * PI) omega = 4 * PI;
	if (omega < 0.0) omega = 0.0;

	fac -> area =  omega * (this_srf -> probe_radius * this_srf -> probe_radius);
}


int check_omega (struct face *fac, double omega)
{
	int problem;

	if (fac == NULL) return (0);
	problem = (omega > 4 * PI + 0.001 || omega < -0.001);
	return (problem);
}

void compute_old_area (struct surface *this_srf)
{
	int nnear, k, i, j, okay, l1;
	double darea, probe_area;
	double dco[3];
	double axes[3][3];
	double *dotptr;
	char message[MAXLINE];
	struct probe *prb, *prb2;
	struct probe **near, **lowptr;
	struct face *fac;
    struct cept *ex;

	if (this_srf -> n_prbdot <= 0) return;
	if (this_srf -> probe_radius <= 0.0) return;
	
	/* accumulate analytical face areas for low probes */
	for (fac = this_srf -> head_face; fac != NULL; fac = fac -> next) {
		if (fac -> problem) {
			inform ("problem face skipped in area/volume computation");
			continue;
		}
		if (fac -> shape != CONCAVE) continue;
		prb = fac -> ptr.prb;
		if (!prb -> low) continue;
		prb -> area += fac -> area;
	}
	
	/* allocate memory for list of low neighboring probes */
	near = (struct probe **)
		allocate_pointers (PROBE, this_srf -> n_low_probe);
	if (near == NULL) {
		ex = new_cept (MEMORY_ERROR,  ALLOCATION,  FATAL_SEVERITY);
		add_function (ex, "compute_old_area");
		add_source (ex, "msvolume.c");
        add_variable (ex, PROBE, "low neighboring probes");
		return;
	}
	
	
	darea = 4 * PI * this_srf -> probe_radius * this_srf -> probe_radius  / this_srf -> n_prbdot;
	
	for (l1 = 0; l1 < this_srf -> n_low_probe; l1++) {
		/* retrieve pointer to probe */
		prb = *(this_srf -> low_probe_hdl + l1);
	
		/* calculate axes and set up circles */
		setup_axis (prb -> center, prb -> atm[0] -> center,
			prb -> atm[1] -> center, axes[0]);
			if (error()) return;
		setup_axis (prb -> center, prb -> atm[1] -> center,
			prb -> atm[2] -> center, axes[1]);
			if (error()) return;
		setup_axis (prb -> center, prb -> atm[2] -> center,
			prb -> atm[0] -> center, axes[2]);
			if (error()) return;
			
		/* form list of low neighbors */
		nnear = 0;
		for (lowptr = this_srf -> low_probe_hdl; lowptr - this_srf -> low_probe_hdl < this_srf -> n_low_probe;
			lowptr++) {
			prb2 = *lowptr;
			if (prb == prb2) continue;
			if (distance (prb -> center, prb2 -> center) < 2 * this_srf -> probe_radius)
				*(near + nnear++) = prb2;
		}
		if (nnear <= 0) continue;
		
		/* calculate area for this probe */
		
		probe_area = 0.0;
		
		for (i = 0; i < this_srf -> n_prbdot; i++) {
			dotptr = this_srf -> prbdot + 3 * i;
			/* check against 3 planes */
			okay = 1;
			for (k = 0; k < 3; k++)
				if (dot_product (dotptr, axes[k]) >= 0.0) {
					okay = 0;
					break;
				}
			if (!okay) continue;
			/* absolute dot coordinates */
			for (k = 0; k < 3; k++)
				dco[k] = prb -> center[k] + *(dotptr + k);
			/* check against near low probes */
			okay = 1;
			for (j = 0; j < nnear; j++) {
				prb2 = *(near + j);
				if (distance (dco, prb2 -> center) < this_srf -> probe_radius) {
					okay = 0;
					break;
				}
			}
			if (okay) probe_area += darea;
		}
		prb -> numerical_area = probe_area;
		if (debug) {
			sprintf (message,
			"low probe analytical area = %10.4f, numerical area = %10.4f",
				prb -> area, prb -> numerical_area);
			inform(message);
		}
	}
	if (near != (struct probe **) NULL)
		free_pointers (PROBE, near);
}


void makdot (struct surface *this_srf)
{					/* make the array of dots on the probe surface */
	int nuse;
	unsigned size;
	double delang, dellat, lat, lon, x, y, xy, z;
	double *dotptr;
	char message[MAXLINE];
    struct cept *ex;
	
	if (this_srf -> probe_radius <= 0.0) {
		this_srf -> n_prbdot = 0;
		this_srf -> prbdot = NULL;
		return;
	}
	/* angle between dots */
	delang = 0.1;
	/* make the number a little higher so we don't overflow */
	this_srf -> n_prbdot = (1.2) * 4.0 * PI / (delang * delang);
	this_srf -> n_prbdot += 12;
	/* allocate memory */
	size = 3 * this_srf -> n_prbdot * sizeof (double);
	this_srf -> prbdot = allocate_doubles (3 * this_srf -> n_prbdot, 0, PRBDOT);
	if (this_srf -> prbdot == NULL) {
		ex = new_cept (MEMORY_ERROR,  ALLOCATION,  FATAL_SEVERITY);
		add_function (ex, "makdot");
        add_variable (ex, PRBDOT, "probe dots");
		add_source (ex, "msvolume.c");
		return;
	}
	
	/* initialization */
	dotptr = this_srf -> prbdot;
	nuse = 0;
	
	/* compute the coordinates of the dots */
	for (lat = delang / 2; lat < PI; lat += delang) {
		z = this_srf -> probe_radius * cos (lat);
		xy = this_srf -> probe_radius * sin (lat);
		if (xy <= 0.0) {
			ex = new_cept (GEOMETRY_ERROR,  NEGATIVE_VALUE,  FATAL_SEVERITY);
			add_function (ex, "makdot");
			add_source (ex, "msvolume.c");
            add_double (ex, "latitude radius", xy);
			return;
		}
		/* angle between dots on latitude */
		dellat = (delang * this_srf -> probe_radius) / xy;
		if (dellat > PI)
			dellat = PI;
		for (lon = dellat / 2; lon <= 2 * PI; lon += dellat) {
			x = xy * cos (lon);
			y = xy * sin (lon);
			if (nuse >= this_srf -> n_prbdot) {
				ex = new_cept (ARRAY_ERROR,  MSOVERFLOW,  FATAL_SEVERITY);
				add_function (ex, "makdot");
				add_source (ex, "msvolume.c");
				add_long (ex, "number of probe dots", (long) nuse);
				return;
			}
			/* store coordinates */
			*dotptr++ = x;
			*dotptr++ = y;
			*dotptr++ = z;
			nuse++;
		}
	}
	this_srf -> n_prbdot = nuse;
	sprintf (message,"%8ld dots made for numerical cusp area", this_srf -> n_prbdot);
	informd(message);
	if (this_srf -> n_prbdot <= 0) {
		ex = new_cept (ARRAY_ERROR,  MSUNDERFLOW,  FATAL_SEVERITY);
		add_function (ex, "makdot");
		add_source (ex, "msvolume.c");
		add_long (ex, "number of probe dots", (long) this_srf -> n_prbdot);
		return;
	}
}



/* VOLUME ROUTINES: */

/* first method: */

void convex_piece_volume (struct face *fac)
{
	struct sphere *atm;

	atm = fac -> ptr.atm;
	fac -> vol1 = (1.0 / 3.0) * atm -> radius * fac -> area;
	fac -> vol2 = 0.0;
}

void saddle_piece_volume (struct face *fac)
{
	int n_fac_arc, k, arc_index;
	double ratio, phi, atm1_theta, atm2_theta;
	double cir_radius_squared;
	double volume_cone1, volume_cone2, volume_saddle1, volume_saddle2;
	double volume_subtract1, volume_subtract2;
	double atom1_radius, atom2_radius, torus_radius;
	double cos1, cos2, sin1, sin2;
	double term1, term2, term3, cusp_angle, cusp_volume;
	double perp_dist, poly_area;
	double vector[3];
	struct arc *arc1, *arc2, *ark;
	struct circle *cir1, *cir2, *cir;
	struct torus *tor;
	struct sphere *atm1, *atm2;
	struct surface *this_srf;
    struct cept *ex;
 
 	this_srf = fac -> srf;
	n_fac_arc = fac -> n_arc;

	if (n_fac_arc == 1) {		/* conical saddle face */
		arc1 = fac -> arcsp[0];
		if (arc1 == NULL) {
			ex = new_cept (POINTER_ERROR,  NULL_VALUE,  FATAL_SEVERITY);
			add_object (ex, SURFACE, "this_srf");
			add_function (ex, "saddle_piece_volume");
			add_source (ex, "msvolume.c");
			return;
		}
		arc2 = NULL;
	}
	else if (n_fac_arc == 2) {	/* hoop */
		arc2 = fac -> arcsp[0];
		if (arc2 == NULL) {
			ex = new_cept (POINTER_ERROR,  NULL_VALUE,  FATAL_SEVERITY);
			add_object (ex, SURFACE, "this_srf");
			add_function (ex, "saddle_piece_volume");
			add_source (ex, "msvolume.c");
			return;
		}
		arc1 = fac -> arcsp[1];
		if (arc1 == NULL) {
			ex = new_cept (POINTER_ERROR,  NULL_VALUE,  FATAL_SEVERITY);
			add_object (ex, SURFACE, "this_srf");
			add_function (ex, "saddle_piece_volume");
			add_source (ex, "msvolume.c");
			return;
		}
	}
	else if (n_fac_arc == 3) {	/* conical */
		arc1 = NULL;
		for (arc_index = 0; arc_index < n_fac_arc; arc_index++) {
            ark = fac -> arcsp[arc_index];
			if (ark == NULL) {
				ex = new_cept (POINTER_ERROR,  NULL_VALUE,  FATAL_SEVERITY);
				add_object (ex, SURFACE, "this_srf");
				add_function (ex, "saddle_piece_volume");
				add_source (ex, "msvolume.c");
				return;
			}
			if (ark -> shape == CONVEX) {
				arc1 = ark;
				break;
			}
		}
		if (arc1 == NULL) {
			ex = new_cept (POINTER_ERROR,  NULL_VALUE,  FATAL_SEVERITY);
			add_function (ex, "saddle_piece_volume");
			add_source (ex, "msvolume.c");
            add_long (ex, "number of arcs in face", (long) n_fac_arc);
			add_message (ex, "bad saddle cone");
			return;
		}
		arc2 = NULL;
	}
	else if (n_fac_arc == 4) {		/* regular saddle face */
		arc2 = fac -> arcsp[1];
		if (arc2 == NULL) {
			ex = new_cept (POINTER_ERROR,  NULL_VALUE,  FATAL_SEVERITY);
			add_object (ex, SURFACE, "this_srf");
			add_function (ex, "saddle_piece_volume");
			add_source (ex, "msvolume.c");
			return;
		}
		arc1 = fac -> arcsp[3];
		if (arc1 == NULL) {
			ex = new_cept (POINTER_ERROR,  NULL_VALUE,  FATAL_SEVERITY);
			add_object (ex, SURFACE, "this_srf");
			add_function (ex, "saddle_piece_volume");
			add_source (ex, "msvolume.c");
			return;
		}
	}
	else {
		ex = new_cept (POINTER_ERROR,  NULL_VALUE,  FATAL_SEVERITY);
		add_function (ex, "saddle_piece_volume");
		add_source (ex, "msvolume.c");
        add_long (ex, "number of arcs in face", (long) n_fac_arc);
		add_message (ex, "bad saddle face");
		return;
	}

	phi = arc1 -> phi;
	cir1 = arc1 -> cir;
	atm1_theta = cir1 -> theta;
	atm1 = cir1 -> atm;
	atom1_radius = atm1 -> radius;

	if (arc2 != NULL) {
		cir2 = arc2 -> cir;
		atm2_theta = cir2 -> theta;
		atm2 = cir2 -> atm;
		atom2_radius = atm2 -> radius;
	}
	else {
		cir2 = NULL;
		atm2_theta = 0.0;
		atm2 = NULL;
		atom2_radius = 0.0;
	}

	tor = fac -> ptr.tor;
	torus_radius = tor -> radius;

	cos1 = cos (atm1_theta);
	cos2 = cos (atm2_theta);
	sin1 = sin (atm1_theta);
	sin2 = sin (atm2_theta);

	/* cone for atom 1 */
	volume_cone1 = (phi / 6.0) * (atom1_radius * atom1_radius *
		atom1_radius) * sin1 * (cos1 * cos1);

	/* inside torus */
	term1 = (torus_radius * torus_radius) * this_srf -> probe_radius * sin1;
	term2 = torus_radius * (this_srf -> probe_radius * this_srf -> probe_radius) *
		(sin1 * cos1 + atm1_theta);
	term3 = (this_srf -> probe_radius * this_srf -> probe_radius * this_srf -> probe_radius / 3.0) *
		(sin1 * (cos1 * cos1) + 2 * sin1);
	volume_saddle1 = (phi / 2.0) * (term1 - term2 + term3);
	cir = arc1 -> cir;
	poly_area =  0.5 * cir -> radius * cir -> radius * arc1 -> phi;
	for (k = 0; k < 3; k++)
		vector[k] = tor -> center[k] - cir -> center[k];
	perp_dist = dot_product (vector, cir -> axis);
	volume_subtract1 = perp_dist * poly_area / 3.0;

	if (arc2 != NULL) {

		/* cone in atom 2 */
		volume_cone2 = (phi / 6.0) * (atom2_radius * atom2_radius *
			atom2_radius) * sin2 * (cos2 * cos2);

		/* inside torus */
		term1 = (torus_radius * torus_radius) * this_srf -> probe_radius * sin2;
		term2 = torus_radius * (this_srf -> probe_radius * this_srf -> probe_radius) *
			(sin2 * cos2 + atm2_theta);
		term3 = (this_srf -> probe_radius * this_srf -> probe_radius * this_srf -> probe_radius / 3.0) *
			(sin2 * (cos2 * cos2) + 2 * sin2);
		volume_saddle2 = (phi / 2.0) * (term1 - term2 + term3);
		cir = arc2 -> cir;
		cir_radius_squared = cir -> radius * cir -> radius;
		poly_area =  0.5 * cir_radius_squared * arc2 -> phi;
		for (k = 0; k < 3; k++)
			vector[k] = tor -> center[k] - cir -> center[k];
		perp_dist = dot_product (vector, cir -> axis);
		volume_subtract2 = perp_dist * poly_area / 3.0;
	}
	else {
		volume_cone2 = 0.0;
		volume_saddle2 = 0.0;
		volume_subtract2 = 0.0;
	}

	/* cusp correction */

	if (arc2 == NULL && this_srf -> probe_radius > 0.0) {
		ratio = tor -> radius / this_srf -> probe_radius;
		if (ratio > 1.0) ratio = 1.0;
		if (ratio < -1.0) ratio = -1.0;
		cusp_angle = acos (ratio);
		cos1 = cos (cusp_angle);
		sin1 = sin (cusp_angle);
		term1 = (torus_radius * torus_radius) * this_srf -> probe_radius * sin1;
		term2 = torus_radius * (this_srf -> probe_radius * this_srf -> probe_radius) *
			(sin1 * cos1 + cusp_angle);
		term3 = (this_srf -> probe_radius * this_srf -> probe_radius * this_srf -> probe_radius / 3.0) *
			(sin1 * (cos1 * cos1) + 2 * sin1);
		cusp_volume = (phi / 2.0) * (term1 - term2 + term3);
	}
	else cusp_volume = 0.0;

	fac -> vol1 = volume_cone1 + volume_saddle1 + volume_saddle2 +
		volume_cone2 - cusp_volume;
	fac -> vol2 = volume_subtract1 + volume_saddle1 + volume_saddle2 +
		volume_subtract2 - cusp_volume;
}

double probe_volume (struct probe *prb)
{
	struct cusp *csp;
	struct face *fac;
	struct cycle *cyc;
	struct edge *edg;
	struct cusp_link *clk;
	int n_cusp_edge;
	double pyramid_volume, sector_volume, prbvol;
	struct edge *cusp_edges[MAX_POLY_EDGE];
	struct surface *this_srf;
	double base_area;
	struct sphere *atm1, *atm2, *atm3;
    struct cept *ex;

	this_srf = prb -> srf;
	atm1 = prb -> atm[0];
	atm2 = prb -> atm[1];
	atm3 = prb -> atm[2];
	base_area = triangle_area (atm1 -> center, atm2 -> center, atm3 -> center);
	pyramid_volume = (1.0 / 3.0) * (prb -> height * base_area);
	prbvol = pyramid_volume;

	/* the usual */
	for (fac = this_srf -> head_face; fac != NULL; fac = fac -> next) {
		if (fac -> shape != CONCAVE) continue;
		if (prb != fac -> ptr.prb) continue;
		if (fac -> problem) {
			inform ("problem face skipped in probe_volume computation");
			continue;
		}
		sector_volume = (1.0 / 3.0) * (this_srf -> probe_radius * fac -> area);
		prbvol -= sector_volume;
	}

	if (!prb -> low) return (prbvol);

	/* cusp circle cones */

	for (clk = prb -> first_cusp; clk != NULL; clk = clk -> next) {
		csp = clk -> csp;
		n_cusp_edge = 0;
		for (fac = this_srf -> head_face; fac != NULL; fac = fac -> next) {
			if (fac -> shape != CONCAVE) continue;
			if (prb != fac -> ptr.prb) continue;
			if (fac -> problem) {
				inform ("problem face skipped in probe_volume computation");
				continue;
			}
			for (cyc = fac -> first_cycle; cyc != NULL; cyc = cyc -> next)
				for (edg = cyc -> first_edge; edg != NULL; edg = edg -> next)
					if (csp == edg -> arcptr -> csp) {
						if (n_cusp_edge >= MAX_POLY_EDGE/2) {
							ex = new_cept (ARRAY_ERROR,  MSOVERFLOW,  FATAL_SEVERITY);
							add_function (ex, "probe_volume");
							add_source (ex, "msvolume.c");
							add_long (ex, "number of cusp edges", (long) n_cusp_edge);
							add_long (ex, "MAX_POLY_EDGE/2", (long) MAX_POLY_EDGE/2);
							return (0.0);
						}
						cusp_edges[n_cusp_edge] = edg;
						++n_cusp_edge;
					}
		}	/* end of face loop */
		if (n_cusp_edge <= 0) continue;		/* eaten ? */
		prbvol -= cusp_cone_volume (prb, csp, n_cusp_edge, cusp_edges);
		if (error()) return(0.0);
	}	/* end of cusp loop */

	return (prbvol);
}

double cusp_cone_volume (struct probe *prb, struct cusp *csp, int n_cusp_edge, struct edge *cusp_edges[MAX_POLY_EDGE])
{
	struct edge polyedges[MAX_POLY_EDGE];
	struct arc polyarcs[MAX_POLY_EDGE];
	struct polygon poly;
	struct vertex *vtx0, *vtx1, *vtx2;
	struct circle *cir;
	struct torus *tor;
	struct cusp_extension *vce, *tce;
	int n_cone_edge, cusp_index, k, orn0, orn1, orn2;
	double vol, sign;
    struct cept *ex;

	sign = ((cusp_edges[0] -> orn) ? 1.0 : -1.0);

	cir = csp -> cir;
	for (k = 0; k < 3; k++) {
		poly.axis[k] = sign * cir -> axis[k];
		poly.center[k] = cir -> center[k];
	}

	poly.first_edge = &(polyedges[0]);

	/* set up first edge of polygon */
	polyarcs[0].vtx[0] = cusp_edges[0] -> arcptr -> vtx[0];
	polyarcs[0].vtx[1] = cusp_edges[0] -> arcptr -> vtx[1];
	polyarcs[0].cir = cusp_edges[0] -> arcptr -> cir;
	polyarcs[0].phi = cusp_edges[0] -> arcptr -> phi;
	polyarcs[0].shape = CONVEX;	/* not concave anymore */
	polyedges[0].arcptr = &(polyarcs[0]);
	polyedges[0].orn = 1 - cusp_edges[0] -> orn;
	polyedges[0].next = NULL;

	/* store starting vertex for cycle close test */
	orn0 = polyedges[0].orn;
	vtx0 = polyedges[0].arcptr -> vtx[orn0];

	/* check for free circle */
	if (vtx0 == NULL) {
		if (n_cusp_edge > 1) {
			ex = new_cept (GEOMETRY_ERROR,  INVALID_VALUE,  FATAL_SEVERITY);
			add_function (ex, "cusp_cone_volume");
			add_source (ex, "msvolume.c");
			add_long (ex, "number of cusp edges", (long) n_cusp_edge);
			add_message (ex, "free circle not only edge");
			return(0.0);
		}
		polygon_area (&poly);
		vol = cone_volume (prb -> center, &poly);
		if (error()) return(0.0);
		return (vol);
	}

	n_cone_edge = 0;

	while (n_cone_edge < MAX_POLY_EDGE) {
		orn1 = polyedges[n_cone_edge].orn;
		vtx1 = polyedges[n_cone_edge].arcptr -> vtx[1 - orn1];
		if (vtx1 == NULL) {
			ex = new_cept (POINTER_ERROR,  NULL_VALUE,  FATAL_SEVERITY);
			add_function (ex, "cusp_cone_volume");
			add_source (ex, "msvolume.c");
			add_long (ex, "number of cone edges", (long) n_cone_edge);
			add_message (ex, "null vertex for curved arc");
			return(0.0);
		}

		/* get opposing vertex */
		vce = vtx1 -> ce;
		if (vce == (struct cusp_extension *) NULL) continue;
		vtx2 = NULL;
		tor = vce -> tor;
		if (tor == NULL) {
			vtx2 = vce -> partner;
		}
		else {
			tce = tor -> ce;
			if (tce == (struct cusp_extension *) NULL) continue;
			if (tce -> vtx[0] == vtx1) vtx2 = tce -> vtx[1];
			else if (tce -> vtx[1] == vtx1) vtx2 = tce -> vtx[0];
		}

		if (vtx2 == NULL) {
			ex = new_cept (LOGIC_ERROR,  NOT_FOUND,  FATAL_SEVERITY);
			add_function (ex, "cusp_cone_volume");
			add_source (ex, "msvolume.c");
			add_message (ex, "cannot find opposing vertex");
			return(0.0);
		}

		/* create straight edge */
		n_cone_edge++;
		if (n_cone_edge >= MAX_POLY_EDGE) {
			ex = new_cept (ARRAY_ERROR,  MSOVERFLOW,  FATAL_SEVERITY);
			add_function (ex, "cusp_cone_volume");
			add_source (ex, "msvolume.c");
			add_long (ex, "number of cone edges", (long) n_cone_edge);
			return(0.0);
		}
		polyarcs[n_cone_edge].vtx[0] = vtx1;
		polyarcs[n_cone_edge].vtx[1] = vtx2;
		polyarcs[n_cone_edge].shape = STRAIGHT;
		polyarcs[n_cone_edge].cir = NULL;
		polyedges[n_cone_edge].arcptr = &(polyarcs[n_cone_edge]);
		polyedges[n_cone_edge].orn = 0;
		polyedges[n_cone_edge-1].next = &(polyedges[n_cone_edge]);
		polyedges[n_cone_edge].next = NULL;

		/* check for end of cycle */
		if (vtx2 == vtx0) {
			if (n_cone_edge+1 != 2 * n_cusp_edge) {
				inform ("(cusp_cone_volume): edge count inconsistency");
				inform ("Will attempt to recover from error.");
				inform ("If no further error messages,");
				inform ("assume error recovery successful,");
				inform ("except for a small area and volume inaccuracy.");
				return (0.0);
			}
			break;
		}

		/* look for following arc */
		for (cusp_index = 0; cusp_index < n_cusp_edge; cusp_index++) {
			/* reverse orientation */
			orn2 = 1 - cusp_edges[cusp_index] -> orn;
			if (cusp_edges[cusp_index] -> arcptr -> vtx[orn2] == vtx2)
				break;
		}
		if (cusp_index >= n_cusp_edge) {
			ex = new_cept (LOGIC_ERROR,  NOT_FOUND,  FATAL_SEVERITY);
			add_function (ex, "cusp_cone_volume");
			add_source (ex, "msvolume.c");
            add_long (ex, "n_cusp_edge", n_cusp_edge);
			add_message (ex, "following arc not found");
			return(0.0);
		}

		/* found */
		/* create curved edge with opposite orientation */
		n_cone_edge++;
		if (n_cone_edge >= MAX_POLY_EDGE) {
			ex = new_cept (ARRAY_ERROR,  MSOVERFLOW,  FATAL_SEVERITY);
			add_function (ex, "cusp_cone_volume");
			add_source (ex, "msvolume.c");
			add_long (ex, "number of cone edges", (long) n_cone_edge);
			return(0.0);
		}
		polyarcs[n_cone_edge].vtx[0] =
			cusp_edges[cusp_index] -> arcptr -> vtx[0];
		polyarcs[n_cone_edge].vtx[1] =
			cusp_edges[cusp_index] -> arcptr -> vtx[1];
		polyarcs[n_cone_edge].cir =
			cusp_edges[cusp_index] -> arcptr -> cir;
		polyarcs[n_cone_edge].phi =
			cusp_edges[cusp_index] -> arcptr -> phi;
		/* not concave anymore */
		polyarcs[n_cone_edge].shape = CONVEX;
		polyedges[n_cone_edge].arcptr = &(polyarcs[n_cone_edge]);
		polyedges[n_cone_edge].orn = 1 - cusp_edges[cusp_index] -> orn;
		polyedges[n_cone_edge-1].next = &(polyedges[n_cone_edge]);
		polyedges[n_cone_edge].next = NULL;
	}

	/* polygon finished */
	polygon_area (&poly);
	vol = cone_volume (prb -> center, &poly);
	if (error()) return(0.0);
	return (vol);
}

/* interior polyhedron volume */
double hedron (struct surface *this_srf)
{
	double polyhedron_volume, prism, average_height, base_area;
	struct face *fac;
	struct probe *prb;
	struct sphere *atm1, *atm2, *atm3, *atm4;

	polyhedron_volume = 0.0;

	/* one triangular prism for each concave face */

	for (fac = this_srf -> head_face; fac != NULL; fac = fac -> next) {
		if (fac -> shape != CONCAVE) continue;
		prb = fac -> ptr.prb;
		atm1 = prb -> atm[0];
		atm2 = prb -> atm[1];
		atm3 = prb -> atm[2];
		atm4 = prb -> atm[3];
		/* later: handle non-planar case */
		if (atm4 == NULL)
			average_height = (1.0 / 3.0) *
				(atm1 -> center[2] + atm2 -> center[2] + atm3 -> center[2]);
		else average_height = (1.0 / 4.0) *
				(atm1 -> center[2] + atm2 -> center[2] +
				 atm3 -> center[2] + atm4 -> center[2]);
		base_area = triangle_area (atm1 -> center, atm2 -> center, atm3 -> center);
		if (atm4 != NULL)
			base_area += triangle_area (atm4 -> center, atm1 -> center, atm3 -> center);
		prism = average_height * prb -> unit_altitude_z * base_area;
		polyhedron_volume += prism;
	}
	return (polyhedron_volume);
}

/* volume routines (2nd method) */

void face_volume (struct face *fac)
{
	int k;
	double edge_join_vol;
	double center[3];
	struct sphere *i;
	struct torus *torus_ptr;
	struct probe *prb;
	struct cycle *cyc;
	struct edge *edg;
	struct surface *this_srf;

	this_srf = fac -> srf;
	switch ((int) fac -> shape) {
	case CONVEX:
		i = fac -> ptr.atm;
		fac -> vol2 = (1.0 / 3.0) * i -> radius * fac -> area;
		for (k = 0; k < 3; k++)
			center[k] = i -> center[k];
		break;
	case SADDLE:
		torus_ptr = fac -> ptr.tor;
		for (k = 0; k < 3; k++)
			center[k] = torus_ptr -> center[k];
		break;
	case CONCAVE:
		prb = fac -> ptr.prb;
		fac -> vol2 = (-1.0 / 3.0) * (this_srf -> probe_radius * fac -> area);
		for (k = 0; k < 3; k++)
			center[k] = prb -> center[k];
		break;
	}

	/* add joins */

	edge_join_vol = 0.0;

	for (cyc = fac -> first_cycle; cyc != NULL; cyc = cyc -> next)
		for (edg = cyc -> first_edge; edg != NULL; edg = edg -> next)
			edge_join_vol += join_volume (this_srf -> center, center, edg);

	fac -> vol2 += edge_join_vol;
}

/* low level measuring routines */

void polygon_center (struct polygon *poly)
{
	int orn, n_poly_side, k;
	double center[3];
	struct edge *edg;
	struct arc *arcptr;
	struct vertex *vtx;
	struct circle *cir;
    struct cept *ex;

	n_poly_side = 0;
	for (k = 0; k < 3; k++)
		center[k] = 0.0;

	for (edg = poly -> first_edge; edg != NULL; edg = edg -> next) {
		arcptr = edg -> arcptr;
		orn = edg -> orn;
		cir = arcptr -> cir;
		vtx = arcptr -> vtx[orn];
		if (cir != NULL && arcptr -> shape == CONVEX) {
			for (k = 0; k < 3; k++)
				poly -> center[k] = cir -> center[k];
			return;
		}
		if (vtx == NULL) {
			ex = new_cept (POINTER_ERROR,  NULL_VALUE,  FATAL_SEVERITY);
			add_function (ex, "polygon_center");
			add_source (ex, "msvolume.c");
			add_message (ex, "polygon vertex");
			return;
		}

		for (k = 0; k < 3; k++)
			center[k] += vtx -> center[k];
		n_poly_side++;
	}

	if (n_poly_side <= 0) {
		ex = new_cept (GEOMETRY_ERROR,  NULL_VALUE,  FATAL_SEVERITY);
		add_function (ex, "polygon_center");
		add_source (ex, "msvolume.c");
        add_long (ex, "n_poly_side", (long) n_poly_side);
		return;
	}

	for (k = 0; k < 3; k++)
		center[k] /= n_poly_side;
	for (k = 0; k < 3; k++)
		poly -> center[k] = center[k];
}

void polygon_area (struct polygon *poly)
{
	int orn;
	double poly_area;
	double area_triangle,  area_segment;
	struct edge *edg;
	struct arc *arcptr;
	struct vertex *vtx0, *vtx1;
	struct circle *cir;

	poly_area = 0.0;
	for (edg = poly -> first_edge; edg != NULL; edg = edg -> next) {
		arcptr = edg -> arcptr;
		orn = edg -> orn;
		vtx0 = arcptr -> vtx[orn];
		vtx1 = arcptr -> vtx[1-orn];
		cir = arcptr -> cir;
		if (vtx0 != NULL && vtx1 != NULL)
			area_triangle = signed_triangle_area (poly -> center,
				vtx0 -> center, vtx1 -> center, poly -> axis);
		else area_triangle = 0.0;
		if (cir != NULL) area_segment = segment_area (arcptr);
		else area_segment = 0.0;
		poly_area += (area_triangle + area_segment);
	}
	poly -> area = poly_area;
}

double segment_area (struct arc *arcptr)
{
	double area_sector, area_triangle, area_segment;
	struct vertex *vtx0, *vtx1;
	struct circle *cir;

	vtx0 = arcptr -> vtx[0];
	vtx1 = arcptr -> vtx[1];
	cir = arcptr -> cir;
	if (cir == NULL) return (0.0);

	area_sector = sector_area (arcptr);

	if (vtx0 != NULL && vtx1 != NULL)
		area_triangle = triangle_area (cir -> center,
			vtx0 -> center, vtx1 -> center);
	else area_triangle = 0.0;

	if (arcptr -> phi > PI) area_triangle *= (-1.0);

	area_segment = area_sector - area_triangle;

	if (arcptr -> shape == CONCAVE) area_segment *= (-1.0);

	return (area_segment);
}

double sector_area (struct arc *arcptr)
{
	double area_sector;
	struct circle *cir;

	cir = arcptr -> cir;
	if (cir == NULL) return (0.0);
	area_sector =  0.5 * cir -> radius * cir -> radius * arcptr -> phi;
	return (area_sector);
}

double cone_volume (double pnt[3], struct polygon *poly)
{
	int k;
	double perp_dist, vol;
	double vector[3];

	for (k = 0; k < 3; k++)
		vector[k] = pnt[k] - poly -> center[k];
	perp_dist = dot_product (vector, poly -> axis);
	vol = perp_dist * poly -> area / 3.0;
	return (vol);
}

double join_volume (double pnt0[3], double pnt1[3], struct edge *edg)
{
	int k, orn;
	double perp_dist, join_vol, tetra_volume1, tetra_volume2;
	double sector_vol;
	double vector[3];
	struct vertex *vtx0, *vtx1;
	struct circle *cir;
	struct arc *arcptr;
    struct cept *ex;

	arcptr = edg -> arcptr;
	orn = edg -> orn;
	if (arcptr -> shape == STRAIGHT) {
		ex = new_cept (ENUM_ERROR,  INVALID_VALUE,  FATAL_SEVERITY);
		add_function (ex, "join_volume");
		add_source (ex, "msvolume.c");
		add_message(ex, "straight arc");
		return(0.0);
	}
	cir = arcptr -> cir;
	if (cir == NULL) {
		ex = new_cept (POINTER_ERROR,  NULL_VALUE,  FATAL_SEVERITY);
		add_function (ex, "join_volume");
        add_object (ex, CIRCLE, "cir");
		add_source (ex, "msvolume.c");
		return(0.0);
	}
	if (arcptr -> phi <= 0.0) return (0.0);
	if (cir -> radius <= 0.0) return (0.0);

	for (k = 0; k < 3; k++)
		vector[k] = pnt1[k] - pnt0[k];
	perp_dist = dot_product (vector, cir -> axis);
	if (orn) perp_dist *= (-1);
	sector_vol = (1.0 / 6.0) * perp_dist * cir -> radius *
		cir -> radius * arcptr -> phi;

	vtx0 = arcptr -> vtx[orn];
	vtx1 = arcptr -> vtx[1-orn];
	if (vtx0 == NULL || vtx1 == NULL) {
		tetra_volume1 = 0.0; tetra_volume2 = 0.0;
	}
	else {
		tetra_volume1 = tetrahedron_volume (pnt0, pnt1, cir -> center,
			vtx0 -> center);
		tetra_volume2 = tetrahedron_volume (pnt0, pnt1, vtx1 -> center,
			cir -> center);
	}
	join_vol = sector_vol + tetra_volume1 + tetra_volume2;
	return (join_vol);
}



/*
	PQMS

	PIECEWISE QUARTIC MOLECULAR SURFACE Computer Program

	Copyright 1986, 1989 by Michael L. Connolly
	All rights reserved
*/
