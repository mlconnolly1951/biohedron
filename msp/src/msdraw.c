/*
 * MSDraw
 * Copyright 1996 by Michael L. Connolly
 * All Rights Reserved

 * February 16, 2006

 * This computer program was written by Michael L. Connolly
 * while self-employed as owner of Biohedron.  It is based
 * partially upon the algorithm in this publication:

 * M. L. Connolly (1985)
 * "Depth Buffer Algorithms for Molecular Modeling",
 * Journal of Molecular Graphics, Volume 3, pages 19-24.
 * Some new tricks have been added.
 */

#include "ms.h"
#include "msenum.h"
#include "mslimit.h"
#include "msdefault.h"
#include "msstruct.h"
#include "msglobal.h"
#include "msfunction.h"

void do_draw (struct msscene *ms, FILE *fpi)
{
	double angle;
	char message[MAXLINE];
	
	ms -> current_molecule = NULL;
	read_draw_commands (ms, fpi);
	if (error()) return;
	if (ms -> raster_format == 0 && ms -> plot_format == 0 && ms -> vector_format == 0) return;

	if (error()) return;
	if (ms -> stereo) {
		angle = ms -> stereo_angle * (3.14159 / 180.0);
		make_matrix (angle, ms -> stereomat);
	}
	init_bounds (ms);		/* to find rotation center */
	if (error()) return;
	init_center (ms);
	if (error()) return;
	transform_all_surface (ms);
	if (error()) return;
	/* to determining scaling */
	init_bounds (ms);
	if (error()) return;
	init_window (ms);
	if (error()) return;
	if (ms -> stereo) {
		rotate_lsource (ms);
		if (error()) return;
		rotate_clip (ms);
		if (error()) return;
	}

	if (ms -> raster_format != 0) {
		init_zbuffer (ms);
		sprintf (message,"%8ld times %8ld = %8ld pixels",
			ms -> horizontal, ms -> vertical, ms -> size);
		inform(message);
	}
	if (ms -> plot_format || ms -> vector_format) {
		if (ms -> clipping) {
			do_clipping (ms);
			if (error()) return;
		}
	}
	if (ms -> plot_format != 0) {
		setup_polygons (ms);
		if (error()) return;
		setup_squares (ms);
		if (error()) return;
		/* use same clipping scaling for left and right images */
		create_linsegs (ms);
		if (error()) return;
		/* mark edges that are backfacing, frontfacing, silhouette or boundary */
		mark_molecules (ms);
		if (error()) return;
	}

	draw_scene (ms);
	if (error()) return;

	if (ms -> raster_format != 0) {
		sprintf (message,"%8ld arc-boundary intersection one-pixel errors", ms -> n_bdy_parity);
		if (ms -> n_bdy_parity > 0) inform(message);
		sprintf (message,"%8ld arc-clipping intersection one-pixel errors", ms -> n_clip_parity);
		if (ms -> n_clip_parity > 0) inform(message);
		sprintf (message,"%8ld containment bad projection one-pixel errors", ms -> n_bad_projection);
		if (ms -> n_bad_projection > 0) inform(message);
		sprintf (message,"%8ld missing surface point one-pixel errors", ms -> n_missing_points);
		if (ms -> n_missing_points > 0) inform(message);
		sprintf (message,"%8ld missing leaf multi-pixel errors", ms -> n_missing_leaves);
		if (ms -> n_missing_leaves > 0) inform(message);
		setup_color_table (ms);
		if (error()) return;
		write_image (ms);
		if (error()) return;
		fclose (ms -> fp_raster);
		free_buffers (ms);
		if (error()) return;
	}
	if (ms -> plot_format != 0) {
		write_plot (ms);
		if (error()) return;
		free_squares (ms);
		if (error()) return;
	}
}


/* INITIALIZATION */


void init_parameters (struct msscene *ms)
{
	int k, j;

	ms -> fineness = FINENESS;
	ms -> border = 0.0;
	ms -> alignment = DEFAULT_ALIGNMENT;
	ms -> interpolate = 1;
	ms -> screendoor = 1;
	ms -> stereo_angle = DEFAULT_STEREO;
	make_matrix (0.0, ms -> stereomat);
	ms -> screen_size = DEFAULT_SCREEN;
	ms -> tiltx = DEFAULT_TILTX;
	ms -> tilty = DEFAULT_TILTY;
	ms -> clip_fraction = 0.0;
	ms -> clipping = 0;
	ms -> n_row = 8;
	ms -> n_column = 8;
	for (j = 0; j < 2; j++)
		for (k = 0; k < 2; k++)
			ms -> viewport[j][k] = ((j == 0) ? 0 : ms -> screen_size);
	ms -> lsource[0] = 0.5;
	ms -> lsource[1] = 0.5;
	ms -> lsource[2] = 1.0;
	normalize (ms -> lsource);
	ms -> model[0] = 0.0;
	ms -> model[1] = 0.0;
	ms -> model[2] = 1.0;
	ms -> model[3] = 0.0;
	ms -> model[4] = 1.0;
	strcpy (ms -> title, "plot");	/* default title */
}

void init_counters (struct msscene *ms)
{
	ms -> n_bad_projection = 0;
	ms -> n_bdy_parity = 0;
	ms -> n_clip_parity = 0;
	ms -> n_missing_points = 0;
	ms -> n_missing_leaves = 0;
}

/* called before and after transformations */
void init_bounds (struct msscene *ms)
{
	int k;
	double low, up;

	for (k = 0; k < 3; k++) {
		ms -> bounds[0][k] = 1000000.0;
		ms -> bounds[1][k] = (-1000000.0);
	}
	if (ms -> head_molecule == NULL) {
		for (k = 0; k < 3; k++) {
			ms -> bounds[0][k] = -100.0;
			ms -> bounds[1][k] = 100.0;
		}
		return;
	}

	for (ms -> current_molecule = ms -> head_molecule; ms -> current_molecule != NULL;
		ms -> current_molecule = ms -> current_molecule -> next) {
		setup_mol_bounds (ms -> current_molecule);
		for (k = 0; k < 3; k++) {
			low = ms -> current_molecule -> bounds[0][k];
			up = ms -> current_molecule -> bounds[1][k];
			if (low < ms -> bounds[0][k]) ms -> bounds[0][k] = low;
			if (up > ms -> bounds[1][k]) ms -> bounds[1][k] = up;
		}
	}
}

void init_center (struct msscene *ms)
{
	int k;

	for (k = 0; k < 3; k++)
		ms -> center[k] = (ms -> bounds[0][k] + ms -> bounds[1][k]) / 2;
}

void init_window (struct msscene *ms)
{
	int k;
	double remainder;
	long quotient;
	double center[3];
	double max_width, horizontal_width, vertical_width;
	char message[MAXLINE];

	/* compute center */
	for (k = 0; k < 3; k++)
		center[k] = (ms -> bounds[0][k] + ms -> bounds[1][k])/ 2;
	horizontal_width = ms -> bounds[1][0] - ms -> bounds[0][0];
	vertical_width = ms -> bounds[1][1] - ms -> bounds[0][1];
	max_width = ((horizontal_width > vertical_width) ? horizontal_width :
		vertical_width);
	if (max_width <= 0.0) {
		set_error1 ("init_window: zero-width molecules");
		return;
	}
	max_width += 2 * ms -> border;
	if (max_width <= 0.0) {
		set_error1 ("init_window: zero-width window");
		return;
	}
	/* if window not specified, figure it out */
	if (ms -> window[0][0] == ms -> window[1][0] ||
		ms -> window[0][1] == ms -> window[1][1]) {
		ms -> window[0][0] = center[0] - max_width / 2;
		ms -> window[1][0] = center[0] + max_width / 2;
		ms -> window[0][1] = center[1] - max_width / 2;
		ms -> window[1][1] = center[1] + max_width / 2;
	}
	ms -> window[0][2] = ms -> bounds[0][2] - 1.0;
	ms -> window[1][2] = ms -> bounds[1][2] + 1.0;
	/* handle alignment */
	if (ms -> alignment > 0.0) {
		for (k = 0; k < 3; k++) {
			quotient = (long) floor ((ms -> window[0][k] / ms -> alignment));
			remainder = ms -> window[0][k] - (ms -> alignment * quotient);
			if (remainder != 0.0)
				ms -> window[0][k] = ms -> alignment * (quotient-1);
			quotient = (long) floor ((ms -> window[1][k] / ms -> alignment));
			remainder = ms -> window[1][k] - (ms -> alignment * quotient);
			if (remainder != 0.0)
				ms -> window[1][k] = ms -> alignment * (quotient+1);
		}
		while ((ms -> window[1][0] - ms -> window[0][0]) >
			(ms -> window[1][1] - ms -> window[0][1]))
			ms -> window[1][1] += ms -> alignment;
		while ((ms -> window[1][0] - ms -> window[0][0]) <
			(ms -> window[1][1] - ms -> window[0][1]))
			ms -> window[1][0] += ms -> alignment;
	}
	horizontal_width = ms -> window[1][0] - ms -> window[0][0];
	vertical_width = ms -> window[1][1] - ms -> window[0][1];
	ms -> plot_width = ((horizontal_width > vertical_width) ?
		horizontal_width : vertical_width);
	sprintf (message,"%8.2f %8.2f x  %8.2f %8.2f y  %8.2f %8.2f z",
		ms -> window[0][0], ms -> window[1][0],
		ms -> window[0][1],ms ->  window[1][1],
		ms -> window[0][2], ms -> window[1][2]);
	inform(message);
	/* initialize some clipping planes */
	ms -> zrange = ms -> window[1][2] - ms -> window[0][2];
	if (ms -> clipping) {
		ms -> clip_center[0] = ms -> center[0];
		ms -> clip_center[1] = ms -> center[1];
		ms -> clip_center[2] = ms -> center[2] + ms -> clip_fraction * ms -> zrange/2;
		ms -> clip_axis[0] = ms -> tiltx;
		ms -> clip_axis[1] = ms -> tilty;
		ms -> clip_axis[2] = 1.0;
		normalize (ms -> clip_axis);
		sprintf (message, "%8.3f clipping plane z-level", ms -> clip_center[2]);
		inform(message);
	}
}

void setup_mol_bounds (struct molecule *mol)
{
	int k;
	double low, up;
	struct surface *grb;

	for (k = 0; k < 3; k++) {
		mol -> bounds[0][k] = 1000000.0;
		mol -> bounds[1][k] = (-1000000.0);
	}
	for (grb = mol -> head_surface; grb != NULL; grb = grb -> next) {
		for (k = 0; k < 3; k++) {
			low = grb -> bounds[0][k];
			up = grb -> bounds[1][k];
			if (low < mol -> bounds[0][k]) mol -> bounds[0][k] = low;
			if (up > mol -> bounds[1][k]) mol -> bounds[1][k] = up;
		}
	}
	for (k = 0; k < 3; k++)
		mol -> centroid[k] =
			(mol -> bounds[0][k] + mol -> bounds[1][k])/2;
}

void setup_surface_bounds (struct surface *obj)
{
	int k;
	long i;
	double low, up;
	struct variety *vty;
	struct phnvtx *pvtx;

	for (k = 0; k < 3; k++) {
		obj -> bounds[0][k] = 1000000.0;
		obj -> bounds[1][k] = (-1000000.0);
	}
	if (obj -> type == PHN_SURFACE ||
		obj -> type == CTR_SURFACE || obj -> type == NML_SURFACE) {
		for (pvtx = obj -> head_phnvtx; pvtx != NULL; pvtx = pvtx -> next) {
			for (k = 0; k < 3; k++) {
				low = pvtx -> center[k];
				up = pvtx -> center[k];
				if (low < obj -> bounds[0][k]) obj -> bounds[0][k] = low;
				if (up > obj -> bounds[1][k]) obj -> bounds[1][k] = up;
			}
		}
	}
	else if (obj -> type == DEN_SURFACE) {
		for (k = 0; k < 3; k++) {
			obj -> bounds[0][k] = obj -> cube_width * obj -> origin[k];
			obj -> bounds[1][k] = obj -> cube_width * obj -> limit[k];
		}
	}
	else {
		for (i = 0; i < obj -> n_atom; i++) {
			vty = *(obj -> variety_handles + i);
			for (k = 0; k < 3; k++) {
				low = vty -> center[k] - vty -> radii[0];
				up = vty -> center[k] + vty -> radii[0];
				if (low < obj -> bounds[0][k]) obj -> bounds[0][k] = low;
				if (up > obj -> bounds[1][k]) obj -> bounds[1][k] = up;
			}
		}
	}
	for (k = 0; k < 3; k++)
		obj -> centroid[k] =
			(obj -> bounds[0][k] + obj -> bounds[1][k])/2;
}

/* TRANSFORMATION */

void transform_all_surface (struct msscene *ms)
{
	struct surface *head_surface, *current_surface;

	for (ms -> current_molecule = ms -> head_molecule;
		ms -> current_molecule != NULL;
		ms -> current_molecule = ms -> current_molecule -> next) {
		head_surface = ms -> current_molecule -> head_surface;
		for (current_surface = head_surface; current_surface != NULL;
			current_surface = current_surface -> next) {
			transform_surface (ms, current_surface);
			setup_surface_bounds (current_surface);
		}
	}
}

/* transform pqms data or bas by rotation and translation */

void transform_surface (struct msscene *ms, struct surface *obj)
{
	int n_side, j, k;
	long i, t;
	double tmpmat[3][3];
	double circle_center[3];
	double center[3], translate[3], rotation[3][3], stereomat[3][3];
	struct variety *vty, **vtyh;
	struct circle *cir, **cirh;
	struct vertex **vtxh;
	struct phnvtx *pvtx;
	struct phntri *tri;
	struct polygon *pgn;
	char message[MAXLINE];

	for (k = 0; k < 3; k++) {
		center[k] = ms -> center[k];
		translate[k] = 0.0;
	}
	for (j = 0; j < 3; j++)
		for (k = 0; k < 3; k++) {
			rotation[j][k] = ms -> rotation[j][k];
			stereomat[j][k] = ms -> stereomat[j][k];
	}
	if (ms -> stereo) {
		for (j = 0; j < 3; j++)
			for (k = 0; k < 3; k++)
				tmpmat[j][k] =
					stereomat[j][0] * rotation[0][k] +
					stereomat[j][1] * rotation[1][k] +
					stereomat[j][2] * rotation[2][k];
		for (j = 0; j < 3; j++)
			for (k = 0; k < 3; k++)
				rotation[j][k] = tmpmat[j][k];
		
	}
	for (j = 0; j < 3; j++)	{
		sprintf (message, "%8.3f %8.3f %8.3f",
			rotation[j][0], rotation[j][1], rotation[j][2]);
		informd (message);
	}
	/* polyhedron of trianglular faces */
	/* transform vertices */
	for (pvtx = obj -> head_phnvtx; pvtx != NULL; pvtx = pvtx -> next) {
		trpnt (rotation, center, translate, pvtx -> center);
		trvec (rotation, pvtx -> outward);
	}
	/* transform polyhedron triangles */
	if (obj -> phntri_handles != NULL) {
		for (t = 0; t < obj -> n_phntri; t++) {
			tri = num2phntri (obj, t+1);
			if (tri == NULL) {
				set_error1 ("transform_surface: null triangle pointer");
				return;
			}
			trpnt (rotation, center, translate, tri -> center);
			trvec (rotation, tri -> axis);
		}
	}
	/* transform polyhedron polygons */
	for (pgn = obj -> head_polygon; pgn != NULL; pgn = pgn -> next) {
		trpnt (rotation, center, translate, pgn -> center);
		n_side = pgn -> n_side;
		for (j = 0; j < n_side; j++) {
			trpnt (rotation, center, translate, pgn -> vc[j]);
		}
		trvec (rotation, pgn -> axis);
	}
	/* transform atom centers */
	if (obj -> atom_centers != NULL) {
		for (i = 0; i < obj -> n_atom; i++)
			trpnt (rotation, center, translate, obj -> atom_centers + 3 * i);
	}
	if (obj -> variety_handles != NULL) {
		for (vtyh = obj -> variety_handles; vtyh - obj -> variety_handles < obj -> n_variety; vtyh++) {
			vty = *vtyh;
			trpnt (rotation, center, translate, vty -> center);
			trvec (rotation, vty -> axis);
			if (vty -> tube) {	/* kludge */
				for (j = 0; j < 2; j++) {
					for (k = 0; k < 3; k++)
						circle_center[k] = vty -> ccens[j][k];
					trpnt (rotation, center, translate, circle_center);
					for (k = 0; k < 3; k++)
						vty -> ccens[j][k] = circle_center[k];
				}
			}
		}
	}
	/* if (obj -> type != PQMS_SURFACE) return; */
	if (obj -> circle_handles != NULL) {
		for (cirh = obj -> circle_handles; cirh - obj -> circle_handles < obj -> n_circle; cirh++) {
			cir = *cirh;
			trpnt (rotation, center, translate, cir -> center);
			trvec (rotation, cir -> axis);
		}
	}
	if (obj -> vertex_handles != NULL) {
		for (vtxh = obj -> vertex_handles; vtxh - obj -> vertex_handles < obj -> n_vertex; vtxh++)
			trpnt (rotation, center, translate, (*vtxh) -> center);
	}
}

void rotate_lsource (struct msscene *ms)
{
	int j;
	double translate[3];

	for (j = 0; j < 3; j++)
		translate[j] = 0.0;
	trvec (ms -> stereomat, ms -> lsource);
}

void rotate_clip (struct msscene *ms)
{
	int j;
	double translate[3];

	for (j = 0; j < 3; j++)
		translate[j] = 0.0;
	trpnt (ms -> stereomat, ms -> center, translate, ms -> clip_center);
	trvec (ms -> stereomat, ms -> clip_axis);
}


void setup_polygons (struct msscene *ms)
{
	int j;
	double xmax, xmin;
	double ymax, ymin;
	double zmax, zmin;
	struct polygon *tri;
	struct molecule *mol;
	struct surface *po;
	
	
	for (mol = ms -> head_molecule; mol != NULL; mol = mol -> next) {
		po = polyhedron_bunch (mol);
		if (po == NULL) continue;
		if (po -> type != PHN_SURFACE && po -> type != DEN_SURFACE) continue;
		for (tri = po -> head_polygon; tri != NULL; tri = tri -> next) {
			xmax = -MS_INFINITY;
			ymax = -MS_INFINITY;
			zmax = -MS_INFINITY;
			xmin =  MS_INFINITY;
			ymin =  MS_INFINITY;
			zmin =  MS_INFINITY;
			for (j = 0; j < tri -> n_side; j++) {
				if (xmax < tri ->vc[j][0])
					xmax = tri -> vc[j][0];
				if (xmin > tri ->vc[j][0])
					xmin = tri -> vc[j][0];
				if (ymax < tri ->vc[j][1])
					ymax = tri -> vc[j][1];
				if (ymin > tri ->vc[j][1])
					ymin = tri -> vc[j][1];
				if (zmax < tri ->vc[j][2])
					zmax = tri -> vc[j][2];
				if (zmin > tri ->vc[j][2])
					zmin = tri -> vc[j][2];
			}
			tri -> xmax = xmax;
			tri -> ymax = ymax;
			tri -> zmax = zmax;
			tri -> xmin = xmin;
			tri -> ymin = ymin;
			tri -> zmin = zmin;
		}
	}
	
}


void min_max (struct molecule *mol)
{
	long v;
	double minv, maxv, f;
	char message[MAXLINE];
	struct phnvtx *vtx;
	struct surface *phn;

	minv = MS_INFINITY;
	maxv = -MS_INFINITY;
	phn = polyhedron_bunch (mol);
	if (phn == NULL) return;
	for (v = 0; v < phn -> n_phnvtx; v++) {
		vtx = *(phn -> phnvtx_handles + v);
		if (vtx == NULL) {
			set_error1 ("min_max: null polyhedron vertex pointer");
			return;
		}
		f = get_function (phn, vtx);
		if (f < minv) minv = f;
		if (f > maxv) maxv = f;
	}
	sprintf (message, "%8.3f min; %8.3f max %2s values", minv, maxv, phn -> function);
	inform(message);
}



void draw_scene (struct msscene *ms)
{
	for (ms -> current_molecule = ms -> head_molecule; 
		ms -> current_molecule != NULL;
		ms -> current_molecule = ms -> current_molecule -> next) {
		draw_molecule (ms);
		if (error()) return;
	}
	if (ms -> raster_format != 0) {
		merge_buffers (ms);
	}
}

void draw_molecule (struct msscene *ms)
{
	int result;
	double fine_pixel;
	char message[MAXLINE];
	struct surface *head_surface;
	struct surface *current_surface, *poly_surface;

	fine_pixel = ms -> pixel_width / ms -> fineness;
	head_surface = ms -> current_molecule -> head_surface;
	sprintf (message,"%8ld graphical objects for molecule: %s",
		ms -> current_molecule -> n_surface, ms -> current_molecule -> name);
	inform(message);

	poly_surface = polyhedron_bunch (ms -> current_molecule);
	for (current_surface = head_surface; current_surface != NULL;
		current_surface = current_surface -> next) {
		if (ms -> current_molecule -> blank && current_surface == poly_surface) continue;
		if (ms -> raster_format != 0) {
			if (current_surface -> solid_shade > 0) {
				current_surface -> db = allocate_buffer (ms, 0);
				if (error()) {
					add_function (tail_cept, "draw_molecule");
					add_source (tail_cept, "msdraw.c");
					return;
				}
			}
			if (current_surface -> type == PQMS_SURFACE) {
				sprintf (message,"%8ld atoms to render for molecular surface",
					current_surface -> n_atom);
				inform(message);
				render_surface (ms, current_surface);
			}
			else if (current_surface -> type == BAS_SURFACE) {
				sprintf (message,"%8ld atoms to render ball_and_stick model",
					current_surface -> n_atom);
				inform(message);
				render_surface (ms, current_surface);
			}
			else if (current_surface -> type == PHN_SURFACE) {
				sprintf (message,"%8ld triangles to render for polyhedron",
					current_surface -> n_phntri);
				inform(message);
				render_polyhedron (ms, current_surface, fine_pixel, ms -> interpolate);
			}
			if (error()) return;
		}
		if (ms -> plot_format != 0) {
			if (current_surface -> n_polygon != 0)
				sprintf (message,"%8ld polygons to plot for molecular surface",
					current_surface -> n_polygon);
			else if (current_surface -> n_phntri > 0)
				sprintf (message,"%8ld triangles to plot for molecular surface",
					current_surface -> n_phntri);
			else if (current_surface -> n_phnedg > 0)
				sprintf (message,"%8ld edges to plot for molecular surface",
					current_surface -> n_phnedg);
			else strcpy (message, "         unknown object being plotted");
			inform (message);
			if (current_surface -> type != PQMS_SURFACE) {
				hle_bunch (ms, current_surface);
				if (error()) return;
			}
		}
		if (ms -> vector_format != 0) {
			write_header_iv (ms, ms -> fp_vector);
			if (current_surface -> type == PQMS_SURFACE) {
				sprintf (message,"%8ld arcs for 3D display list", current_surface -> n_arc);
				inform(message);
				result = write_arc_iv (ms, current_surface, ms -> fp_vector, 18L);
				if (result == 0 || error()) return;
			}
			else if (current_surface -> type == BAS_SURFACE) {
				sprintf (message,"%8ld edges for ball_and_stick 3D display list",
						current_surface -> n_phnedg);
				inform(message);
				result = write_edg_iv (ms, current_surface , ms -> fp_vector);
			}
			else if (current_surface -> type == PHN_SURFACE) {
				if (current_surface -> n_polygon != 0) {
					sprintf (message,"%8ld polygons for 3D display list", current_surface -> n_polygon);
					inform (message);
					result = write_phn_iv (ms, current_surface, ms -> fp_vector);
				}
				else if (current_surface -> n_phntri > 0) {
					sprintf (message,"%8ld triangles for 3D display list", current_surface -> n_phntri);
					inform (message);
					result = write_tri_iv (ms, current_surface, ms -> fp_vector, NULL);
				}
				else if (current_surface -> n_phnedg > 0) {
					sprintf (message,"%8ld edges for 3D display list", current_surface -> n_phnedg);
					inform (message);
					result = write_edg_iv (ms, current_surface , ms -> fp_vector);
				}
				else strcpy (message, "         unknown object being 3D");
				if (error()) return;
			}
			else if (current_surface -> type == CTR_SURFACE) {
				if (current_surface -> n_phnctr > 0) {
					result = write_ctr_iv (ms, current_surface , ms -> fp_vector);
				}
			}
			else if (current_surface -> type == NML_SURFACE) {
				if (current_surface -> n_phnedg > 0) {
					result = write_edg_iv (ms, current_surface , ms -> fp_vector);
				}
			}
		}
	}
}

void draw_mem ()
{
	unsigned long size;
	int type;
	char type_name[32];
	
	init_mem ((unsigned long) N_OBJECTS);
	if (error()) return;

	type = MSSCENE;
	size = sizeof (struct msscene);
	strcpy (type_name, "msscene");
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

	type = VARIETY;
	size = sizeof (struct variety);
	strcpy (type_name, "variety");
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

	type = LEAF;
	size = sizeof (struct leaf);
	strcpy (type_name, "leaf");
	define_type (type, size, type_name);
	if (error()) return;

	type = MOLECULE;
	size = sizeof (struct molecule);
	strcpy (type_name, "molecule");
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

	type = LINSEG;
	size = sizeof (struct linseg);
	strcpy (type_name, "linseg");
	define_type (type, size, type_name);
	if (error()) return;

	type = PHNVTX;
	size = sizeof (struct phnvtx);
	strcpy (type_name, "phnvtx");
	define_type (type, size, type_name);
	if (error()) return;

	type = PHNEDG;
	size = sizeof (struct phnedg);
	strcpy (type_name, "phnedg");
	define_type (type, size, type_name);
	if (error()) return;

	type = PHNTRI;
	size = sizeof (struct phntri);
	strcpy (type_name, "phntri");
	define_type (type, size, type_name);
	if (error()) return;
	
	type = POLYGON;
	size = sizeof (struct polygon);
	strcpy (type_name, "polygon");
	define_type (type, size, type_name);
	if (error()) return;

	type = PHNCTR;
	size = sizeof (struct phnctr);
	strcpy (type_name, "phnctr");
	define_type (type, size, type_name);
	if (error()) return;
	
	type = EVALPNT;
	size = sizeof (struct evalpnt);
	strcpy (type_name, "evalpnt");
	define_type (type, size, type_name);
	if (error()) return;
	
	type = SOLID_ANGLE;
	size = sizeof (struct solid_angle);
	strcpy (type_name, "solid_angle");
	define_type (type, size, type_name);
	if (error()) return;
	
	type = COMPONENT;
	size = sizeof (struct component);
	strcpy (type_name, "component");
	define_type (type, size, type_name);
	if (error()) return;
	
	type = SURFACE;
	size = sizeof (struct surface);
	strcpy (type_name, "surface");
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

	type = EDGER;
	size = sizeof (struct edger);
	strcpy (type_name, "edger");
	define_type (type, size, type_name);
	if (error()) return;

	type = SQUARE;
	size = sizeof (struct square);
	strcpy (type_name, "square");
	define_type (type, size, type_name);
	if (error()) return;

	type = LAX;
	size = sizeof (struct lax);
	strcpy (type_name, "lax");
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
	
	type = MATERIAL_TABLE;
	size = sizeof (struct material_table);
	strcpy (type_name, "material_table");
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
	
	type = CEPT;
	size = sizeof (struct cept);
	strcpy (type_name, "cept");
	define_type (type, size, type_name);
	if (error()) return;
	
	type = DEPTH_BUFFER;
	size = sizeof (struct depth_buffer);
	strcpy (type_name, "depth_buffer");
	define_type (type, size, type_name);
	if (error()) return;
	
}

/*
	MSDraw
	Copyright 1993 by Michael L. Connolly
	All Rights Reserved
*/
