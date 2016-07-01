/* Molecular Surface Package Copyright 1993 by Michael L. Connolly */
/* February 8, 2000 */

#include "ms.h"
#include "msenum.h"
#include "mslimit.h"
#include "msdefault.h"
#include "msstruct.h"
#include "msglobal.h"
#include "msfunction.h"

int write_header_iv (struct msscene *ms, FILE *fpiv)
{
	if (ms != NULL && ms -> vrml)
		fprintf (fpiv, "#VRML V1.0 ascii\n\n");
	else
		fprintf (fpiv, "#Inventor V2.0 ascii\n\n");
	write_material_iv (ms -> table, fpiv, 0);
	return (1);
}

int write_material_iv (struct material_table *table, FILE *fpiv, int emissive)
{
	long i, nmaterial;
	
	nmaterial = table -> nmaterial;
	fprintf (fpiv, "Material { \n");
	if (emissive) fprintf (fpiv, "emissiveColor [ \n");
	else fprintf (fpiv, "diffuseColor [ \n");
	for (i = 0; i < nmaterial; i++) {
		fprintf (fpiv, "%6.3f %6.3f %6.3f",
			table -> red[i], table -> green[i], table -> blue[i]);
		if (i != nmaterial - 1) fprintf (fpiv, ",");
		fprintf (fpiv, "\n");
	}
	fprintf (fpiv, "] \n");
	fprintf (fpiv, "}\n");
	return (1);
}

int write_edg_iv (struct msscene *ms, struct surface *phn, FILE *fpiv)
{
	int material0, material1, which_value, color_type;
	int inner, ind, hue, ccomp;
	long v, e, v0, v1;
	char message[MAXLINE];
	double value0, value1;
	struct phnvtx *pv, *pv0, *pv1;
	struct phnedg *pe;
	struct object_scheme *scheme;

	inner = 0;
	scheme = phn -> scheme;
	if (scheme == NULL) {
		set_error1 ("write_edg_iv: missing color scheme");
		return (0);
	}
	which_value = scheme -> which_value;
	color_type = scheme -> color_type;
	fprintf (fpiv, "Separator {\n");
	fprintf (fpiv, "    MaterialBinding {\n");
	fprintf (fpiv, "	value	PER_VERTEX_INDEXED\n");
	fprintf (fpiv, "    }\n");
	fprintf (fpiv, "    Coordinate3 {\n");
	fprintf (fpiv, "	point	[ \n");

	/* write vertex coordinates */
	for (v = 0; v < phn -> n_phnvtx; v++) {
		pv = num2phnvtx (phn, v + 1);
		if (pv == NULL) return (0);
		fprintf (fpiv, "%9.5f %9.5f %9.5f,\n",
			pv -> center[0], pv -> center[1], pv -> center[2]);
	}
	fprintf (fpiv, "    ] \n");
	fprintf (fpiv, "    }\n");


	sprintf (message,"%8ld vertices written to Inventor file", phn -> n_phnvtx);
	inform(message);

	fprintf (fpiv, "    IndexedLineSet {\n");
	fprintf (fpiv, "	coordIndex	[ \n");
	for (e = 0; e < phn -> n_phnedg; e++) {
		pe = num2phnedg (phn, e + 1);
		if (pe == NULL) return (0);
		v0 = pe -> vtxnum[0]; v1 = pe -> vtxnum[1];
		fprintf (fpiv, "%6ld, %6ld,",
			v0-1, v1-1);
		fprintf (fpiv, " -1,\n");
	}
	fprintf (fpiv, "    ] \n");

	fprintf (fpiv, "	materialIndex	[ \n");
	for (e = 0; e < phn -> n_phnedg; e++) {
		pe = num2phnedg (phn, e + 1);
		if (pe == NULL) return (0);
		v0 = pe -> vtxnum[0]; v1 = pe -> vtxnum[1];
		pv0 = num2phnvtx (phn, v0);
		if (pv0 == NULL) return (0);
		pv1 = num2phnvtx (phn, v1);
		if (pv1 == NULL) return (0);
		if (color_type == VERTEX_COLORING && which_value >= 0 && scheme -> ramp != NULL) {
			value0 = pv0 -> values[which_value];
			material0 = value2material(scheme, value0);
			value1 = pv1 -> values[which_value];
			material1 = value2material(scheme, value1);
		}
		else if (color_type == ATOM_COLORING) {
			material0 = pe -> hue;
			material1 = pe -> hue;
		}
		else if (color_type == UNIFORM_COLORING) {
			ind = inner;
			hue = scheme -> uniform_colors[ind];
			material0 = hue;
			material1 = hue;
		}
		else if (color_type == COMPONENT_COLORING) {
			ccomp = (pe -> comp >= 2) ? 2 : 1;
			ind = 2 * (ccomp - 1) + inner;
			hue = scheme -> component_colors[ind];
			material0 = hue;
			material1 = hue;
		}
		else {
			material0 = 7;
			material1 = 7;
		}
		fprintf (fpiv, "%6d, %6d, ",
			material0, material1);
		fprintf (fpiv, " -1,\n");
	}
	fprintf (fpiv, "    ] \n");
	fprintf (fpiv, "    }\n");
	
	fprintf (fpiv, "}\n");
	sprintf (message,"%8ld edges written to Inventor file", phn -> n_phnedg);
	inform(message);
	return (1);
}


int write_arc_iv (struct msscene *ms, struct surface *srf, FILE *fpiv, long max_point)
{
	char message[MAXLINE];
	double center[3];
	double axis[3]; double radius; double ends[2][3];
	int noends; int m; long l;
	struct arc *a;
	int i,  k;
	long hue;
	double theta;
	double cos_angle, sin_angle, angle;
	double base[3], altitude[3];
	long vtxnum[2];
	struct circle *cir;
	struct vertex *vtx[2];
	double points[MAX_ARC_POINT][3];
	long n_point;
	double vectors[2][3];


	if (max_point <= 1) max_point = 2;
	if (max_point > MAX_ARC_POINT) {
		sprintf (message, "write_arc_iv: max_point = %ld > limit (%ld)",
			max_point, (long) MAX_ARC_POINT);
		set_error1 (message);
		return (0);
	}
	fprintf (fpiv, "Separator {\n");
	for (l = 0; l < srf -> n_arc; l++) {
		a = *(srf -> arc_handles + l);
		cir = a -> cir;
		for (k = 0; k < 3; k++) {
			axis[k] = cir -> axis[k];
			center[k] = cir -> center[k];
		}
		radius = cir -> radius;
		for (m = 0; m < 2; m++) {
			vtx[m] = a -> vtx[m];
			vtxnum[m] = ((vtx[m] == NULL) ? 0 : vtx[m] -> number);
		}
		noends = (vtxnum[0] == vtxnum[1]);

		/* set up theta, base and altitude */
		if (noends) {
			theta = 2 * PI;
			arbprp (axis, base);
		}
		else {
			for (m = 0; m < 2; m++) {
				for (k = 0; k < 3; k++)
					ends[m][k] = vtx[m] -> center[k];
			}
			for (k = 0; k < 3; k++) {
				vectors[0][k] = ends[0][k] - center[k];
				vectors[1][k] = ends[1][k] - center[k];
			}
			theta = positive_angle (vectors[0], vectors[1], axis);
			for (k = 0; k < 3; k++)
				base[k] = (ends[0][k] - center[k]) / radius;
		}
		cross (axis, base, altitude);
		n_point = theta * (max_point / (2 * PI));
		if (n_point < 2) n_point = 2;
		if (n_point >= max_point-1) n_point = max_point-1;
		/* write a polygonal curve */
		for (i = 0; i <= n_point; i++) {
			angle = i * theta / n_point;
			cos_angle = (double) cos ((double) angle);
			sin_angle = (double) sin ((double) angle);
			for (k = 0; k < 3; k++) {
				points[i][k] = center[k] + radius *
					(cos_angle * base[k] + sin_angle * altitude[k]);
			}
		}
		hue = a -> cir -> subtype + 1;
		write_one_arc (ms, points, n_point + 1, hue, fpiv);
	}
	fprintf (fpiv, "}\n");

	sprintf (message,"%8ld arcs written to Inventor file", srf -> n_arc);
	inform(message);
	return (1);
}

int write_one_arc (struct msscene *ms, double points[MAX_ARC_POINT][3], long n_point, long hue, FILE *fpiv)
{
	long i;

	fprintf (fpiv, "Separator {\n");
	fprintf (fpiv, "    MaterialBinding {\n");
	fprintf (fpiv, "    value  PER_VERTEX_INDEXED\n");
	fprintf (fpiv, "    }\n");
	fprintf (fpiv, "    Coordinate3 {\n");
	fprintf (fpiv, "    point  [ \n");
	for (i = 0; i < n_point; i++) {
		fprintf (fpiv, "%9.5f %9.5f %9.5f,\n",
			points[i][0], points[i][1], points[i][2]);
	}
	fprintf (fpiv, "    ] \n");
	fprintf (fpiv, "    }\n");
	fprintf (fpiv, "    IndexedLineSet {\n");
	fprintf (fpiv, "	coordIndex	[ \n");
	for (i = 0; i < n_point; i++) {
		fprintf (fpiv, "%6ld,", i);
	}
	fprintf (fpiv, " -1,\n");
	fprintf (fpiv, "    ] \n");

	fprintf (fpiv, "	materialIndex	[ \n");
	for (i = 0; i < n_point; i++) {
		fprintf (fpiv, "%6ld,", hue);
	}
	fprintf (fpiv, " -1,\n");
	fprintf (fpiv, "    ] \n");
	fprintf (fpiv, "    }\n");
	fprintf (fpiv, "}\n");
	return (1);
}


int write_tri_iv (struct msscene *ms, struct surface *phn, FILE *fpiv, struct color_ramp *ramp)
{
	int which_value;
	int result;
	long n_phnvtx;
	long n_phntri;
	struct object_scheme *ivms;

	ivms = phn -> scheme;
	if (ivms == NULL) {
		set_error1 ("write_tri_iv: null scheme");
		return (0);
	}
	which_value = ivms -> which_value;
	/* use default yellow-green-blue color ramp if none specified */
	if (ramp == NULL)
		ivms -> ramp = lookup_ramp (ms, "green");
	else ivms -> ramp = ramp;
	if (ivms -> ramp == NULL) {
		set_error1 ("write_tri_iv: null color ramp ramp");
		return(0);
	}
	n_phnvtx = material_phn_vtx (phn);
	if (n_phnvtx <= 0) {
		inform("(write_tri_iv): no vertices for polyhedron");
		return(0);
	}
	n_phntri = material_phn_tri (phn);
	if (n_phntri <= 0) {
		inform("(write_tri_iv): no triangles for polyhedron");
		return(0);
	}
	result = write_pgn_iv (ms, phn, n_phnvtx, 0, fpiv);
	if (result <= 0) return (0);
	return (1);
}


int write_mol_iv (struct msscene *ms, struct molecule *mol, FILE *fp_tri)
{
	return (1);
}

int write_ctr_iv (struct msscene *ms, struct surface *phn, FILE *fpiv)
{
	int material0, material1, which_value;
	long v, v0, v1;
	char message[MAXLINE];
	double value0, value1;
	struct phnvtx *pv, *pv0, *pv1;
	struct phnedg *pe;
	struct phnctr *pc;
	struct object_scheme *scheme;

	scheme = phn -> scheme;
	which_value = scheme -> which_value;

	fprintf (fpiv, "Separator {\n");
	fprintf (fpiv, "    MaterialBinding {\n");
	fprintf (fpiv, "	value	PER_VERTEX_INDEXED\n");
	fprintf (fpiv, "    }\n");
	fprintf (fpiv, "    Coordinate3 {\n");
	fprintf (fpiv, "	point	[ \n");

	/* write vertex coordinates */
	for (v = 0; v < phn -> n_phnvtx; v++) {
		pv = num2phnvtx (phn, v + 1);
		if (pv == NULL) return (0);
		fprintf (fpiv, "%9.5f %9.5f %9.5f,\n",
			pv -> center[0], pv -> center[1], pv -> center[2]);
	}
	fprintf (fpiv, "    ] \n");
	fprintf (fpiv, "    }\n");


	sprintf (message,"%8ld vertices written to Inventor file", phn -> n_phnvtx);
	inform(message);

	fprintf (fpiv, "    IndexedLineSet {\n");
	fprintf (fpiv, "	coordIndex	[ \n");
	for (pc = phn -> head_phnctr; pc != NULL; pc = pc -> next) {
		v0 = pc -> head_phnedg -> vtxnum[0];
		fprintf (fpiv, "%6ld,", v0-1);
		for (pe = pc -> head_phnedg; pe != NULL; pe = pe -> next_ctredg) {
			v1 = pe -> vtxnum[1];
			fprintf (fpiv, "%6ld,", v1-1);
		}
		fprintf (fpiv, " -1,\n");
	}
	fprintf (fpiv, "    ] \n");

	fprintf (fpiv, "	materialIndex	[ \n");
	for (pc = phn -> head_phnctr; pc != NULL; pc = pc -> next) {
		v0 = pc -> head_phnedg -> vtxnum[0];
		pv0 = num2phnvtx (phn, v0);
		if (pv0 == NULL) return (0);
		value0 = pv0 -> values[which_value];
		if (scheme -> color_type == VERTEX_COLORING && scheme -> ramp != NULL) {
			value0 = pv0 -> values[which_value];
			material0 = value2material(scheme, value0);
		}
		else {
			material0 = scheme -> uniform_colors[0];
		}
		fprintf (fpiv, "%6d,", material0);
		for (pe = pc -> head_phnedg; pe != NULL; pe = pe -> next_ctredg) {
			v1 = pe -> vtxnum[1];
			pv1 = num2phnvtx (phn, v1);
			if (pv1 == NULL) return (0);
			if (scheme -> color_type == VERTEX_COLORING && scheme -> ramp != NULL) {
				value1 = pv1 -> values[which_value];
				material1 = value2material(scheme, value1);
			}
			else {
				material1 = scheme -> uniform_colors[0];
			}
			fprintf (fpiv, "%6d,", material1);
			if (error()) return (0);
		}
		fprintf (fpiv, " -1,\n");
	}
	fprintf (fpiv, "    ] \n");
	fprintf (fpiv, "    }\n");
	
	fprintf (fpiv, "}\n");

	sprintf (message,"%8ld edges written to Inventor file", phn -> n_phnedg);
	inform(message);
	sprintf (message,"%8ld contours written to Inventor file", phn -> n_phnctr);
	inform(message);
	return (1);
}

int write_pgn_iv (struct msscene *ms, struct surface *phn, long nnormal, int outline, FILE *fpiv)
{
	int j;
	char message[MAXLINE];
	long v, p, m, n, n_side;
	long v_written, e_written, p_written;
	struct phnvtx **phnvtx_handles;
	long n_phnvtx;
	struct polygon **polygon_handles;
	long n_polygon;
	struct phntri **phntri_handles;
	long n_phntri;
	struct phnvtx *pv;
	struct phntri *pt;
	struct polygon *pp;
	struct object_scheme *ivms;

	v_written = 0; e_written = 0; p_written = 0;
	ivms = phn -> scheme;
	n_phnvtx = phn -> n_phnvtx;
	n_polygon = phn -> n_polygon;
	n_phntri = phn -> n_phntri;
	phnvtx_handles = phn -> phnvtx_handles;
	phntri_handles = phn -> phntri_handles;
	polygon_handles = phn -> polygon_handles;

	
	fprintf (fpiv, "Separator {\n");
	if (outline) {
		fprintf (fpiv, "    DrawStyle {\n");
		fprintf (fpiv, "	style	LINES\n");
		fprintf (fpiv, "    }\n");
	}
	if (ivms -> color_type != VERTEX_COLORING) {
		fprintf (fpiv, "    MaterialBinding {\n");
		fprintf (fpiv, "	value	PER_FACE_INDEXED\n");
		fprintf (fpiv, "    }\n");
	}
	else if (ivms -> color_type == VERTEX_COLORING) {
		fprintf (fpiv, "    MaterialBinding {\n");
		fprintf (fpiv, "	value	PER_VERTEX_INDEXED\n");
		fprintf (fpiv, "    }\n");
	}
	if (outline) {
		write_material_iv (ms -> table, fpiv, 1);
		fprintf (fpiv, "    Material { \n");
		fprintf (fpiv, "    diffuseColor [ \n");
		fprintf (fpiv, "    %6.2f %6.2f %6.2f, \n", 0.0, 0.0, 0.0);
		fprintf (fpiv, "     ] \n");
		fprintf (fpiv, "    emissiveColor [ \n");
		for (m = 0; m < ms -> table -> nmaterial; m++)
			fprintf (fpiv, "    %6.2f %6.2f %6.2f, \n",
				ms -> table -> red[m],
				ms -> table -> green[m],
				ms -> table -> blue[m]);
		fprintf (fpiv, "     ] \n");
		fprintf (fpiv, "    }\n");
	}
	fprintf (fpiv, "    ShapeHints {\n");
	fprintf (fpiv, "	vertexOrdering	COUNTERCLOCKWISE\n");
	fprintf (fpiv, "    }\n");
	fprintf (fpiv, "    ShapeHints {\n");
	fprintf (fpiv, "	creaseAngle	0\n");
	fprintf (fpiv, "    }\n");
	fprintf (fpiv, "    Coordinate3 {\n");
	fprintf (fpiv, "	point	[ \n");

	/* write vertex coordinates */
	for (v = 0; v < n_phnvtx; v++) {
		pv = *(phnvtx_handles + v);
		fprintf (fpiv, "%9.5f %9.5f %9.5f,\n",
			pv -> center[0], pv -> center[1], pv -> center[2]);
		v_written++;
	}
	fprintf (fpiv, "    ] \n");
	fprintf (fpiv, "    }\n");
	fprintf (fpiv, "    ShapeHints {\n");
	fprintf (fpiv, "	creaseAngle	0\n");
	fprintf (fpiv, "    }\n");

	/*  NORMALS     NORMALS    NORMALS */
	if (nnormal > 0) {
		fprintf (fpiv, "    Normal {\n");
		fprintf (fpiv, "	vector	[ \n");
		/* write vertex normals */
		for (v = 0; v < nnormal; v++) {
			pv = *(phnvtx_handles + v);
			fprintf (fpiv, "%9.5f %9.5f %9.5f,\n",
				pv -> outward[0], pv -> outward[1], pv -> outward[2]);
		}
		fprintf (fpiv, "    ] \n");
		fprintf (fpiv, "    }\n");
	}
	fprintf (fpiv, "    NormalBinding {\n");
	fprintf (fpiv, "	value	DEFAULT\n");
	fprintf (fpiv, "    }\n");

	sprintf (message,"%8ld vertices written to Inventor file",
		n_phnvtx);
	inform(message);

	fprintf (fpiv, "    IndexedFaceSet {\n");
	fprintf (fpiv, "	coordIndex	[ \n");
	n = (n_polygon > 0) ? n_polygon : n_phntri;
	for (p = 0; p < n; p++) {
		if (polygon_handles != NULL) {
			pp = *(polygon_handles + p);
			if (pp -> clipped) continue;
		}
		else pp = NULL;
		if (phntri_handles != NULL)
			pt = *(phntri_handles + p);
		else pt = NULL;
		if (pp != NULL) n_side = pp -> n_side;
		else n_side = 3;
		for (j = 0; j < n_side; j++) {
			if (pp != NULL)
				v = pp -> vertex_index[j];
			else v = pt -> vtxnum[j] - 1;
			if (v < 0 || v >= n_phnvtx) {
				sprintf(message, "invalid vertex index: %8ld, n_phnvtx = %8ld",
					v, n_phnvtx);
				set_error1 (message);
				return (0);
			}
			fprintf (fpiv, "%6ld, ", v);
		}
		p_written++;
		fprintf (fpiv, " -1,\n");
	}
	fprintf (fpiv, "    ] \n");
	/* normals */
	if (nnormal > 0) {
		fprintf (fpiv, "	normalIndex	[ \n");
		n = (n_polygon > 0) ? n_polygon : n_phntri;
		for (p = 0; p < n; p++) {
			if (polygon_handles != NULL) {
				pp = *(polygon_handles + p);
				if (pp -> clipped) continue;
			}
			else pp = NULL;
			if (phntri_handles != NULL)
				pt = *(phntri_handles + p);
			else pt = NULL;
			if (pp != NULL) n_side = pp -> n_side;
			else n_side = 3;
			for (j = 0; j < n_side; j++) {
				if (pp != NULL)
					v = pp -> vertex_index[j];
				else v = pt -> vtxnum[j] - 1;
				if (v < 0 || v >= n_phnvtx) {
					sprintf(message, "invalid vertex index: %8ld, n_phnvtx = %8ld",
						v, n_phnvtx);
					set_error1 (message);
					return (0);
				}
				fprintf (fpiv, "%6ld, ", v);
			}
			fprintf (fpiv, " -1,\n");
		}
		fprintf (fpiv, "    ] \n");
	}
	if (ivms -> color_type == VERTEX_COLORING) {
		fprintf (fpiv, "	materialIndex	[ \n");
		n = (n_polygon > 0) ? n_polygon : n_phntri;
		for (p = 0; p < n; p++) {
			if (polygon_handles != NULL) {
				pp = *(polygon_handles + p);
				if (pp -> clipped) continue;
			}
			else pp = NULL;
			if (phntri_handles != NULL)
				pt = *(phntri_handles + p);
			else pt = NULL;
			if (pp != NULL) n_side = pp -> n_side;
			else n_side = 3;
			for (j = 0; j < n_side; j++) {
				if (pp != NULL)
					v = pp -> vertex_index[j];
				else v = pt -> vtxnum[j] - 1;
				if (v < 0 || v >= n_phnvtx) {
					sprintf(message, "invalid vertex index: %8ld, n_phnvtx = %8ld",
						v, n_phnvtx);
					set_error1 (message);
					return (0);
				}
				pv = *(phnvtx_handles + v);
				fprintf (fpiv, "%6ld, ", pv -> material_index);
			}
			fprintf (fpiv, " -1,\n");
		}
		fprintf (fpiv, "    ] \n");
	}
	else if (ivms -> color_type == UNIFORM_COLORING) {
		fprintf (fpiv, "	materialIndex	[ \n");
		n = (n_polygon > 0) ? n_polygon : n_phntri;
		for (p = 0; p < n; p++) {
			if (polygon_handles != NULL) {
				pp = *(polygon_handles + p);
				if (pp -> clipped) continue;
			}
			else pp = NULL;
			if (phntri_handles != NULL)
				pt = *(phntri_handles + p);
			else pt = NULL;
			fprintf (fpiv, "%6ld, \n", ivms -> uniform_colors[0]);
		}
		fprintf (fpiv, "    ] \n");
	}
	else {
		fprintf (fpiv, "	materialIndex	[ \n");
		n = (n_polygon > 0) ? n_polygon : n_phntri;
		for (p = 0; p < n; p++) {
			if (polygon_handles != NULL) {
				pp = *(polygon_handles + p);
				if (pp -> clipped) continue;
			}
			else pp = NULL;
			if (phntri_handles != NULL)
				pt = *(phntri_handles + p);
			else pt = NULL;
			if (pp != NULL)
				fprintf (fpiv, "%6ld, \n", pp -> material_index);
			else
				fprintf (fpiv, "%6ld, \n", pt -> material_index);
		}
		fprintf (fpiv, "    ] \n");
	}
	fprintf (fpiv, "    }\n");
	
	fprintf (fpiv, "}\n");
	if (n_polygon > 0)
		sprintf (message,"%8ld polygons written to SGI Inventor ascii file",
			p_written);
	else sprintf (message,"%8ld triangles written to SGI Inventor ascii file",
			n_phntri);
	inform(message);
	return (1);
}

int write_phn_iv (struct msscene *ms, struct surface *phn, FILE *fpiv)
{
	int result;
	long n_polygon, n_phnvtx, n_edge;
	char message[MAXLINE];

	n_edge = collect_edges (phn);
	if (error()) return (0);
	if (n_edge <= 0) {
		inform("(write_phn_iv); no edges for polyhedron");
		return(0);
	}
	sprintf(message, "%8ld edges collected", n_edge);
	inform (message);
	n_phnvtx = material_phn_vtx (phn);
	if (n_phnvtx <= 0) {
		inform("(write_phn_iv): no vertices for polyhedron");
		return(0);
	}
	n_polygon = material_phn_pgn (phn);
	if (n_polygon <= 0) {
		inform("(write_phn_iv): no polygons for polyhedron");
		return(0);
	}
	result = write_pgn_iv (ms, phn, n_phnvtx, 0, fpiv);
	if (result <= 0) return (0);
	return (1);
}


int write_den_iv (struct msscene *ms, struct surface *den, double ctrlev, FILE *fpiv)
{
	int result;
	long n_phnvtx, n_polygon;

	n_polygon = material_den_pgn (den,  ctrlev);
	if (n_polygon <= 0) return (0);
	n_phnvtx = den -> n_phnvtx;
	result = write_pgn_iv (ms, den, n_phnvtx, 0, fpiv);
	if (result <= 0) return (0);
	return (1);
}


long material_phn_vtx (struct surface *phn)
{
	long num;
	long n_phnvtx;
	double val;
	struct phnvtx *pv;
	struct object_scheme *ivms;

	ivms = phn -> scheme;
	if (ivms == NULL) {
		set_error1 ("material_phn_vtx: null object scheme for polyhedron");
		return (0L);
	}
	n_phnvtx = phn -> n_phnvtx;
	/* store vertex materials */
	for (num = 1; num <= n_phnvtx; num++) {
		pv = num2phnvtx (phn, num);
		/* under user choice */
		if (ivms -> which_value >= 0 && ivms -> ramp != NULL) {
			val = pv -> values[ivms -> which_value];
			pv -> material_index = value2material(ivms, val);
		}
		else {
			pv -> material_index = pv -> hue;
		}
	}
	return (n_phnvtx);
}

long material_phn_tri (struct surface *phn)
{
	long idx;
	long n_phntri;
	struct phntri *pt;
	struct phntri **phntri_handles;

	n_phntri = phn -> n_phntri;
	phntri_handles = phn -> phntri_handles;
	for (idx = 0; idx < n_phntri; idx++) {
		pt = *(phntri_handles + idx);
		if (pt == NULL) return (0);
		pt -> material_index = pt -> hue;
	}
	return (n_phntri);
}

long material_phn_pgn (struct surface *phn)
{
	long idx;
	long n_polygon;
	struct polygon *poly;
	struct polygon **polygon_handles;

	n_polygon = phn -> n_polygon;
	polygon_handles = phn -> polygon_handles;
	for (idx = 0; idx < n_polygon; idx++) {
		poly = *(polygon_handles + idx);
		if (poly == NULL) return (0);
		poly -> material_index = poly -> hue;
	}
	return (n_polygon);
}


long material_den_pgn (struct surface *den, double ctrlev)
{
	long ix, ix0, ix1, iy, iy0, iy1, iz, iz0, iz1;
	long idx011;
	long idx101, idx110, idx111;
	double den011;
	double den101;
	double den110;
	double den111;
	long size1, size2;
	long n_polygon;
	struct polygon **polygon_handles, **ip;
	struct polygon *poly;

	size1 = den -> width[0];
	size2 = size1 * den -> width[1];
	polygon_handles = den -> polygon_handles;
	n_polygon = den -> n_polygon;
	ip = polygon_handles;
	for (iz = den -> origin[2]+1; iz < den -> limit[2]; iz++) {
		for (iy = den -> origin[1]+1; iy < den -> limit[1]; iy++) {
			for (ix = den -> origin[0]+1; ix < den -> limit[0]; ix++) {
				ix0 = ix - den -> origin[0] - 1; ix1 = ix0 + 1;
				iy0 = iy - den -> origin[1] - 1; iy1 = iy0 + 1;
				iz0 = iz - den -> origin[2] - 1; iz1 = iz0 + 1;
				/* idx000 = iz0 * size2 + iy0 * size1 + ix0; */
				/* idx001 = iz0 * size2 + iy0 * size1 + ix1; */
				/* idx010 = iz0 * size2 + iy1 * size1 + ix0; */
				idx011 = iz0 * size2 + iy1 * size1 + ix1;
				/* idx100 = iz1 * size2 + iy0 * size1 + ix0; */
				idx101 = iz1 * size2 + iy0 * size1 + ix1;
				idx110 = iz1 * size2 + iy1 * size1 + ix0;
				idx111 = iz1 * size2 + iy1 * size1 + ix1;
				den011 = *(den -> densities + idx011);
				den101 = *(den -> densities + idx101);
				den110 = *(den -> densities + idx110);
				den111 = *(den -> densities + idx111);
				if (ip - polygon_handles >= n_polygon) break;
				poly = *ip;
				if (poly == NULL) return (0);
				if      (den110 < ctrlev && ctrlev < den111) {
					poly -> material_index = 0;
					ip++;
				}
				else if (den110 > ctrlev && ctrlev > den111) {
					poly -> material_index = 0;
					ip++;
				}
				if (ip - polygon_handles >= n_polygon) break;
				poly = *ip;
				if (poly == NULL) return (0);
				if (den101 < ctrlev && ctrlev < den111) {
					poly -> material_index = 0;
					ip++;
				}
				else if (den101 > ctrlev && ctrlev > den111) {
					poly -> material_index = 0;
					ip++;
				}
				if (ip - polygon_handles >= n_polygon) break;
				poly = *ip;
				if (poly == NULL) return (0);
				if (den011 < ctrlev && ctrlev < den111) {
					poly -> material_index = 0;
					ip++;
				}
				else if (den011 > ctrlev && ctrlev > den111) {
					poly -> material_index = 0;
					ip++;
				}
			}
		}
	}
	if (ip - polygon_handles != n_polygon) {
		set_error1 ("inconsistent polygon count");
		return (0);
	}
	return (n_polygon);
}







