/*

	MSRoll

	Copyright 1986, 1989 by Michael L. Connolly
	All rights reserved

	Written by Michael L. Connolly.
	December 16, 2001

*/

#include "ms.h"
#include "msenum.h"
#include "mslimit.h"
#include "msdefault.h"
#include "msstruct.h"
#include "msglobal.h"
#include "msfunction.h"


int ds_stream (struct msscene *ms, int atoms, FILE *fpe, FILE *fpj, FILE *fpk, double pr, int connected)
{
	/* local variables, arrays and structures for ds_stream */
    int     i, j, okay;
    int     maxatom;            /* maximum number of atoms */
    int     errflg;             /* error flag returned from ds */
    double   carea, rarea;       /* contact and reentrant areas */
	double  *atmco;                  /* pointer to atomic coordinate array */
	double  *atmrad;                 /* pointer to atomic radii array */
	double  *atmden;                 /* pointer to array of surface densities */
	short  *atmatt;                 /* pointer to array of attention numbers */
    char    sform[MAXLINE];     /* surface point format */
    char    snform[MAXLINE];    /* surface point with normal format */
    char message[MAXLINE];
    double dot_area;
    struct dsdesc  *curdsd;     /* current surface descriptor */
	/* pointers for writing out atom's area and surface points */
    struct cluster *clu1, *clu2;
    struct surface *dsnmls;
    
    int k, result;
    int atom_set, atom_number, atom_srn;
    double atom_radius, atom_angle, atom_density, average_density, factor;
    double atom_center[3];
    double length, red, green, blue;
	char atom_name[MAX_ATNAME], residue[MAX_ATNAME], sequence[MAX_ATNAME];

 	/* initialization */
    curdsd = NULL;
	/* control strings for surface points and surface normals */
    /* strcpy (sform, "%5d %c %8.3f %8.3f %8.3f\n");
    strcpy (snform, "%5d %c %8.3f %8.3f %8.3f %6.3f %6.3f %6.3f\n"); */
    strcpy (sform, "%5d%5d%5d%2d%9.3f%9.3f%9.3f\n");
    strcpy (snform, "%5d%5d%5d%2d%9.3f%9.3f%9.3f%7.3f%7.3f%7.3f%7.3f\n");
    
    atom_set = atoms;
    maxatom = count_set (atom_set);
	/* check value of maximum number of atoms */
	if (maxatom < 1) {
		sprintf (message, "ds: maxatom = %5d too small", maxatom);
		set_error1 (message);
		return (0);
	}
	if (maxatom > MAXATOM) {
		sprintf (message, "ds: maxatom = %5d too large", maxatom);
		set_error1 (message);
		return (0);
	}
	/* allocate space for:
		coordinates, radii, densities, attention numbers */
	atmco = allocate_doubles (3 * maxatom, 0, CENTERS);
	if (atmco == NULL) {
		sprintf (message, "ds: cannot allocate memory for coordinates");
		set_error1 (message);
		return (0);
	}
	atmrad = allocate_doubles (maxatom, 0, RADII);
	if (atmrad == NULL) {
		sprintf (message, "ds: cannot allocate memory for radii");
		set_error1 (message);
		return (0);
	}
	atmden = allocate_doubles (maxatom, 0, ATMDEN);
	if (atmden == NULL) {
		sprintf (message, "ds: cannot allocate memory for densities");
		set_error1 (message);
		return (0);
	}
	atmatt = allocate_shorts (maxatom);
	if (atmatt == NULL) {
		sprintf (message, "ds: cannot allocate memory for attention");
		set_error1 (message);
		return (0);
	}
	/* allocate new surface descriptor */
	curdsd = (struct dsdesc *) allocate_object (DSDESC);
	if (curdsd == NULL) {
		sprintf (message, "ds: mem alloc fail for: descriptor");
		set_error1 (message);
		return (0);
	}
	curdsd -> intern = NULL;
	curdsd -> maxatom = maxatom;
	curdsd -> atmco = atmco;
	curdsd -> atmrad = atmrad;
	curdsd -> atmden = atmden;
	curdsd -> atmatt = atmatt;
	i = 0;              /* initialize atom number */
	average_density = 0.0;
	for (atom_number = init_for (atom_set); atom_number != 0;
		atom_number = next_for (atom_set)) {
		if (i >= maxatom) break;
		get_atom_center (atom_number, atom_center);
		if (error()) {
			sprintf (message, "ds: invalid atom number %d", atom_number);
			set_error1 (message);
			return (0);
		}
		atom_radius = get_atom_radius (atom_number);
		atom_srn = get_atom_srn (atom_number);
		atom_angle = get_atom_angle (atom_number);
		atom_density = 1.0 / (atom_radius * atom_radius * atom_angle * atom_angle);
		for (k = 0; k < 3; k++)
			*(curdsd -> atmco + 3 * i + k) = atom_center[k];
		*(curdsd -> atmrad + i) = atom_radius;
		*(curdsd -> atmatt + i) = (short) atom_srn;
		*(curdsd -> atmden + i) = atom_density;
		average_density += atom_density;
		/* if attention >= 2, increment for consistency with ms */
		if (*(curdsd -> atmatt + i) >= 2)
			*(curdsd -> atmatt + i) += 1;
		/* require normals for inventor output */
		*(curdsd -> atmatt + i) = 4;
		i++;
	}
	if (i <= 0) {
		set_error1 ("ds_stream: no atoms");
		return (0);
	}
	curdsd -> natom = i;        /* store number of atoms read */
	average_density /= curdsd -> natom;
	/* informative messages */
	sprintf (message, "%8ld atoms", curdsd -> natom);
	inform (message);
	sprintf (message, "%8.3f average surface point density", average_density);
	inform (message);

	/* defaults: */
	curdsd -> pradius = DEFAULT_PROBE;
	curdsd -> connected = DEFAULT_CONNECTED;
	curdsd -> verbose = 1;
	if (pr >= 0.0) curdsd -> pradius = pr;
	if (connected >= 0) curdsd -> connected = connected;

	/* calculate a dot surface for the molecule */
	/* informative messages */
	sprintf (message, "%8.3f probe radius", curdsd -> pradius);
	inform (message);
	if (curdsd -> connected)
		inform ("         connected surface");
	else
		inform ("         complete surface");
	inform ("         start calculation");

	/* CALCULATE SURFACE */
	if (errflg = ds (curdsd)) {
		/* error flag non-zero */
		if (errflg == 550)
			return (0);
		sprintf (message, "ds: ERROR in %s", curdsd -> errstr);
		set_error1 (message);
		sprintf (message, "stage: %s; atom: %5ld",
			curdsd -> stage, curdsd -> atom + 1);
		set_error2 (message);
		return (0);
	}
	/* summary messages */
	/* check whether any atom's att >= POINTS */
	/* if not, there's no point in saying 0 points were computed */
	okay = 0;
	for (i = 0; i < curdsd -> natom; i++)
		if (*(curdsd -> atmatt + i) >= POINTS)
			okay = 1;
	if (okay) {
		sprintf (message, "%8d surface points", curdsd -> nsrfpnt);
		inform (message);
	}
	sprintf (message, "%8.2f contact   surface area", curdsd -> cvxarea);
	inform (message);

	/* check whether any atom's att = AREA */
	/* if so, reentrant area is underreported */
	okay = 1;
	for (i = 0; i < curdsd -> natom; i++)
		if (*(curdsd -> atmatt + i) == AREA)
			okay = 0;
	if (okay) {
		sprintf (message, "%8.2f reentrant surface area", curdsd -> renarea);
		inform(message);
		sprintf (message, "%8.2f total     surface area",
			curdsd -> cvxarea + curdsd -> renarea);
		inform(message);
	}

	if (errflg = dsouch (curdsd)) {
		/* error flag non-zero */
		if (errflg == 550) {
			set_error1 ("ds fails");
			return (0);
		}
		sprintf (message, "dsstream: ERROR in %s", curdsd -> errstr);
		set_error1 (message);
		sprintf (message, "stage: %s   atom: %5d",
			curdsd -> stage, curdsd -> atom + 1);
		set_error2 (message);
		return (0);
	}

	/* write atoms surface areas or points to file */
	if (fpj != NULL || fpk != NULL) {
		/* write output file */
		/* loop through all the atoms of the molecule */
		for (atom_number = init_for (atom_set), i = 0; i < curdsd -> natom && atom_number != 0;
			atom_number = next_for (atom_set), i++) {
			/* each atom has two cluster pointers (possibly null)
				one is for contact (convex) surface
				the other is for reentrant surface
				each cluster has an area, and optionally,
				surface points, and optionally, surface normals */
			clu1 = *(curdsd -> cvxsp + i);  /* convex */
			clu2 = *(curdsd -> rensp + i);  /* reentrant */
			if (fpj != NULL && clu1 != NULL)   {       /* write cluster */
				if (clu1 -> nmem > 0)
					dot_area = clu1 -> area / clu1 -> nmem;
				else dot_area = 0.0;
				for (j = 0; j < clu1 -> nmem; j++)
					if (clu1 -> normals == NULL)
						fprintf (fpj, sform, i + 1, 0, 0, 1,
							*(clu1 -> points + 3 * j),
							*(clu1 -> points + 3 * j + 1),
							*(clu1 -> points + 3 * j + 2));
					else
						fprintf (fpj, snform, i + 1, 0, 0, 1,
							*(clu1 -> points + 3 * j),
							*(clu1 -> points + 3 * j + 1),
							*(clu1 -> points + 3 * j + 2),
							dot_area,
							*(clu1 -> normals + 3 * j),
							*(clu1 -> normals + 3 * j + 1),
							*(clu1 -> normals + 3 * j + 2));
			}
			if (fpj != NULL && clu2 != NULL) {       /* write cluster */
				if (clu2 -> nmem > 0)
					dot_area = clu2 -> area / clu2 -> nmem;
				else dot_area = 0.0;
				for (j = 0; j < clu2 -> nmem; j++)
					if (clu2 -> normals == NULL)
						fprintf (fpj, sform, i + 1, i + 1, 0, 2,
							*(clu2 -> points + 3 * j),
							*(clu2 -> points + 3 * j + 1),
							*(clu2 -> points + 3 * j + 2));
					else
						fprintf (fpj, snform, i + 1, i + 1, 0, 2,
							*(clu2 -> points + 3 * j),
							*(clu2 -> points + 3 * j + 1),
							*(clu2 -> points + 3 * j + 2),
							dot_area,
							*(clu2 -> normals + 3 * j),
							*(clu2 -> normals + 3 * j + 1),
							*(clu2 -> normals + 3 * j + 2));
			}

			if (fpk != NULL) {              /* write cluster areas */
				if (clu1 == NULL) carea = 0.0;
				else carea = clu1 -> area;
				if (clu2 == NULL) rarea = 0.0;
				else rarea = clu2 -> area;
				atom_radius = *(curdsd -> atmrad + i);
				factor = (atom_radius + curdsd -> pradius) / atom_radius;
				factor = factor * factor;
				get_atom_name (atom_number, atom_name);
				get_atom_group (atom_number, residue);
				get_atom_sequence (atom_number, sequence);
				fprintf (fpk, "%-5s %-5s %-5s %8.3f %8.3f %8.3f %8.3f\n",
						residue, sequence, atom_name, carea, rarea, carea + rarea, carea * factor);
			}
		}/* end of cluster loop */
	}
	if (fpj != NULL) fclose (fpj);
	if (fpk != NULL) fclose (fpk);
	length = 0.5; 
	free_doubles (atmco, 0, CENTERS);
	free_doubles (atmrad, 0, RADII);
	free_doubles (atmden, 0, ATMDEN);
	free_shorts (atmatt);
	return (1);
}


/* normals */

struct surface *ds_normals (struct dsdesc  *curdsd, double length)
{
	long j, k;
	long i, v;
	double vtx_center[3], vtx_outward[3];
	char message[MAXLINE];
	struct surface *obj;
	struct phnvtx *vtx1, *vtx2, *nml_vtx;
	struct phnvtx **nml_vertices;
	struct phnedg *nml_edg;
	struct phnedg **nml_edges;
    struct cluster *clu1, *clu2;
	
	sprintf (message,"%8.3f length for polyhedron normals", length);
	inform (message);

	/* allocate normals bunch */
	obj = new_surface ();
	if (obj == NULL) {
		set_error1 ("msroll (ds_normals): mem alloc failure");
		return(NULL);
	}
	obj -> type = NML_SURFACE;

	/* allocate memory for normal vertices */
	nml_vertices = (struct phnvtx **)
		allocate_pointers (PHNVTX, 2 * curdsd -> nsrfpnt);
	if (nml_vertices == NULL) {
		set_error1 ("msroll (ds_normals): memory full");
		return(NULL);
	}

	/* allocate memory for normal edges */
	nml_edges = (struct phnedg **)
		allocate_pointers (PHNEDG, curdsd -> nsrfpnt);
	if (nml_edges == NULL) {
		set_error1 ("msroll (ds_normals): memory full");
		return(NULL);
	}

	/* loop through all the atoms of the molecule */
	for (i = 0, v = 0; i < curdsd -> natom; i++) {
		/* each atom has two cluster pointers (possibly null)
			one is for contact (convex) surface
			the other is for reentrant surface
			each cluster has an area, and optionally,
			surface points, and optionally, surface normals */
		clu1 = *(curdsd -> cvxsp + i);  /* convex */
		if (clu1 != NULL) {
			if (clu1 -> normals == NULL) {
				set_error1 ("ds_normals: missing normals");
				return (NULL);
			}
			for (j = 0; j < clu1 -> nmem; j++, v++) {
				for (k = 0; k < 3; k++) {
					vtx_center[k] =	*(clu1 -> points + 3 * j + k);
					vtx_outward[k] =	*(clu1 -> normals + 3 * j + k);
				}
				/* create normal vertices and edges */
				vtx1 = allocate_phnvtx ();
				if (vtx1 == (struct phnvtx *) NULL) {
					set_error1 ("msplot (polyhedron_normals): vertex allocation failure");
					return ((struct surface *) NULL);
				}
				vtx1 -> number = 2 * v + 1;
				*(nml_vertices + 2 * v) = vtx1;
				for (k = 0; k < 3; k++)
					vtx1 -> center[k] = vtx_center[k];
				vtx2 = allocate_phnvtx ();
				if (vtx2 == (struct phnvtx *) NULL) {
					set_error1 ("msroll (ds_normals): vertex allocation failure");
					return ((struct surface *) NULL);
				}
				vtx2 -> number = 2 * v + 2;
				*(nml_vertices + 2 * v + 1) = vtx2;
				for (k = 0; k < 3; k++)
					vtx2 -> center[k] = vtx_center[k] + length * vtx_outward[k];
				nml_edg = allocate_phnedg ();
				if (nml_edg == (struct phnedg *) NULL) {
					set_error1 ("msroll (ds_normals): edge allocation failure");
					return ((struct surface *) NULL);
				}
				*(nml_edges + v) = nml_edg;
				nml_edg -> pvt[0] = vtx1;
				nml_edg -> pvt[1] = vtx2;
				nml_edg -> vtxnum[0] = vtx1 -> number;
				nml_edg -> vtxnum[1] = vtx2 -> number;
			} /* end of surface points of first cluster loop */
		}

		clu2 = *(curdsd -> rensp + i);  /* reentrant */
		if (clu2 != NULL) {
			if (clu2 -> normals == NULL) {
				set_error1 ("ds_normals: missing normals");
				return (NULL);
			}
			for (j = 0; j < clu2 -> nmem; j++, v++) {
				for (k = 0; k < 3; k++) {
					vtx_center[k] =	*(clu2 -> points + 3 * j + k);
					vtx_outward[k] =	*(clu2 -> normals + 3 * j + k);
				}
				/* create normal vertices and edges */
				vtx1 = allocate_phnvtx ();
				if (vtx1 == (struct phnvtx *) NULL) {
					set_error1 ("msplot (polyhedron_normals): vertex allocation failure");
					return ((struct surface *) NULL);
				}
				vtx1 -> number = 2 * v + 1;
				*(nml_vertices + 2 * v) = vtx1;
				for (k = 0; k < 3; k++)
					vtx1 -> center[k] = vtx_center[k];
				vtx2 = allocate_phnvtx ();
				if (vtx2 == (struct phnvtx *) NULL) {
					set_error1 ("msroll (ds_normals): vertex allocation failure");
					return ((struct surface *) NULL);
				}
				vtx2 -> number = 2 * v + 2;
				*(nml_vertices + 2 * v + 1) = vtx2;
				for (k = 0; k < 3; k++)
					vtx2 -> center[k] = vtx_center[k] + length * vtx_outward[k];
				nml_edg = allocate_phnedg ();
				if (nml_edg == (struct phnedg *) NULL) {
					set_error1 ("msroll (ds_normals): edge allocation failure");
					return ((struct surface *) NULL);
				}
				*(nml_edges + v) = nml_edg;
				nml_edg -> pvt[0] = vtx1;
				nml_edg -> pvt[1] = vtx2;
				nml_edg -> vtxnum[0] = vtx1 -> number;
				nml_edg -> vtxnum[1] = vtx2 -> number;
			} /* end of surface points of second cluster loop */
		} /* end of surface points of second cluster loop */
	}/* end of cluster loop */


	/* set up linked list */
	obj -> head_phnvtx = *nml_vertices;
	for (v = 0; v < 2 * curdsd -> nsrfpnt - 1; v++) {
		nml_vtx = *(nml_vertices + v);
		if (nml_vtx == NULL) {
			sprintf (message, "ds normal vertex %8ld null", v + 1);
			set_error1 (message);
			return (NULL);
		}
		nml_vtx -> next = *(nml_vertices + v + 1);
	}
	obj -> head_phnedg = *nml_edges;
	for (v = 0; v < curdsd -> nsrfpnt - 1; v++) {
		nml_edg = *(nml_edges + v);
		if (nml_edg == NULL) {
			sprintf (message, "ds normal edge %8ld null", v + 1);
			set_error1 (message);
			return (NULL);
		}
		nml_edg -> next = *(nml_edges + v + 1);
	}
		
	sprintf (message,"%8.3f normal length;  %8ld normal vectors",
		length, curdsd -> nsrfpnt);
	inform(message);
	obj -> scheme = NULL;
	/* transfer to structure variables */
	obj -> n_polygon = 0L;
	obj -> n_phnvtx = 2 * curdsd -> nsrfpnt;
	obj -> n_phnedg = curdsd -> nsrfpnt;
	obj -> head_phnvtx = *nml_vertices;
	obj -> head_phnedg = *nml_edges;
	obj -> head_polygon = (struct polygon *) NULL;
	
	obj -> phnvtx_handles = nml_vertices;
	obj -> phnedg_handles = nml_edges;

	return (obj);
}


int pqms_stream (struct msscene *ms, struct molecule *mol, FILE *fpe, FILE *fpq, FILE *fpa, FILE *fpv, double pr, double grid, int dox, char *order, double surface_center[3])
{
	int result;
	struct surface *this_srf;

	if (fpe == NULL) return (0);
	this_srf = new_surface ();
	if (this_srf == NULL) return (0);
	this_srf -> n_atom = mol -> n_atom;
	this_srf -> head_atom = (struct sphere *) (mol -> head_atom);
	this_srf -> tail_atom = (struct sphere *) (mol -> tail_atom);
	this_srf -> mol = mol;
	this_srf -> do_cusp_intersections = dox;
	ms -> this_srf = this_srf;
	ms -> purge_frequency = DEFAULT_PURGE_FREQUENCY;
	this_srf -> probe_radius = DEFAULT_PROBE;
	this_srf -> bit_width = DEFAULT_BIT_WIDTH;
	this_srf -> cube_width = DEFAULT_CUBE_WIDTH;
	if (pr >= 0.0) this_srf -> probe_radius = pr; /* otherwise, use default */
	if (grid > 0.0) {
		this_srf -> bit_width = grid;
		this_srf -> cube_width = 2 * this_srf -> bit_width;
		this_srf -> use_grid = 1;
	}
	calculate_surface (ms -> this_srf);
	if (error())  {
		return (0);
	}
	if (this_srf -> n_problem_atom > 0) {
		inform ("msroll: some problem atoms encountered, quit");
		return (0);
	}
	find_components (this_srf); if (error()) return (0);
	this_srf -> surface_completed = 1;

	/* always compute the volume */
	compute_volume (this_srf);
	if (error())  {
		return (0);
	}
	this_srf -> volume_computed = 1;

	result = cvt_surface (this_srf, surface_center);
	if (error())  {
		return (0);
	}
	
	if (fpa != NULL) {
		result = write_area (this_srf, fpa, order);
		if (error()) {
			return (0);
		}
		if (!result) {
			inform("msroll: problem writing analytica area file");
			return (0);
		}
		if (fpa != stdout && fpa != stderr) fclose (fpa);
	}
	if (fpv != NULL) {
		fprintf (fpv, "molecule name           = %s\n", this_srf -> mol -> name);
		result = write_volume1 (this_srf, fpv);
		if (error()) {
			return (0);
		}
		if (!result) {
			inform("msroll: problem writing analytical volume file");
			return (0);
		}
		if (fpv != stdout && fpv != stderr) fclose (fpv);
	}
	if (fpq != NULL) {
		result = write_surface (this_srf, fpq);
		if (error()) {
			return (0);
		}
		if (!result) {
			inform("msroll: problem writing analytical surface file");
			return (0);
		}
	}
	return (1);
}


int trb_stream (struct msscene *ms, FILE *fpe, FILE *fpt, FILE *fpc, struct molecule *mol)
{
	int result;
	int selected_comp;
	struct surface *tphn;
	struct surface *cphn;


	if (ms == NULL) {
		set_error1 ("trb_stream: no scene specified");
		return (0);
	}
	if (ms -> this_srf == NULL) {
		set_error1 ("trb_stream: no surface specified");
		return (0);
	}
	ms -> this_srf -> scheme = NULL;
	sort_faces (ms -> this_srf);
	count_problem_faces (ms -> this_srf);
	if (!original_lfn (ms -> this_srf)) return(0);
	if (error()) return(0);
	if (!original_face_number (ms -> this_srf)) return(0);
	if (error()) return(0);
	ms -> this_mol = mol;
	setup_alpha (ms -> this_mol);
	if (error()) {
		return (0);
	}
	face_alpha (ms -> this_mol, ms -> this_srf);
	if (error()) {
		return (0);
	}

	result = triangulate(ms);
	if (error()) {
		return (0);
	}
	if (!result) {
		inform("msroll: problem triangulating analytical surface");
		return (0);
	}
	clear_pqms (ms -> this_srf);
	if (error()) return (0);
	result = write_volume2 (ms -> this_srf);
	if (error()) {
		return (0);
	}
	if (!result) {
		inform("msroll: problem writing polyhedral volume file");
		return (0);
	}
	if (fpc != NULL && ms -> this_srf -> n_component >= 2) {
		selected_comp = 1;
		tphn = filter_phn (ms -> this_srf);
		cphn = tphn -> next;
		free_phn (ms -> this_srf);
		ms -> this_srf = tphn;
		if (cphn != NULL) {
			/* write file of triangles */
			result = write_vet (cphn, fpc);
			if (error()) {
				return (0);
			}
			if (!result) {
				inform("msroll: problem writing cavity file");
				return (0);
			}
		}
	}
	if (fpt != NULL) {
		/* write file of triangles */
		result = write_vet (ms -> this_srf, fpt);
		if (error()) {
			return (0);
		}
		if (!result) {
			inform("msroll: problem writing polyhedron file");
			return (0);
		}
	}

	return (1);
}



/*
   MSRoll
   Copyright 1996 by Michael L. Connolly
   All Rights Reserved

*/
