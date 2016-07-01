#include "ms.h"
#include "msenum.h"
#include "mslimit.h"
#include "msdefault.h"
#include "msstruct.h"
#include "msglobal.h"
#include "msfunction.h"

/* last revised February 23, 2006 */

/* DRAW USER INTERFACE */


void read_draw_commands (struct msscene *ms, FILE *fp_command)
{
	int n_words, result, w;
	char command_line[MAXLINE+1];
	char words[MAX_WORD][MAX_NAME];
	char verb[MAX_NAME];

	for (;;) {
		fgets (command_line, MAXLINE, fp_command);
		if (feof (fp_command)) break;
		for (w = 0; w < MAX_WORD; w++)
			strcpy (words[w], " ");
		null_terminate (command_line, '#');
		n_words = sscanf (command_line, "%s %s %s %s %s %s %s %s %s %s",
			words[0], words[1], words[2], words[3], words[4],
			words[5], words[6], words[7], words[8], words[9]);
		if (n_words <= 0) continue;
		if (words[0][0] == '!') continue;
		if (words[0][0] == '#') continue;
		strcpy (verb, words[0]);
		if (strcmp (verb, "source") == 0) {
			result = parse_source (ms, n_words, words);
		}
		else if (strcmp (verb, "shading_model") == 0) {
			result = parse_shading_model (ms, n_words, words);
		}
		else if (strcmp (verb, "light_source") == 0) {
			result = parse_light_source (ms, n_words, words);
		}
		else if (strcmp (verb, "rotation") == 0) {
			result = parse_rotation (ms, n_words, words);
		}
		else if (strcmp (verb, "elbow") == 0) {
			result = parse_elbow (ms, n_words, words);
		}
		else if (strcmp (verb, "outer_line_width") == 0) {
			result = parse_outer_line_width (ms, n_words, words);
		}
		else if (strcmp (verb, "inner_line_width") == 0) {
			result = parse_inner_line_width (ms, n_words, words);
		}
		else if (strcmp (verb, "cavity_line_width") == 0) {
			result = parse_cavity_line_width (ms, n_words, words);
		}
		else if (strcmp (verb, "bond_line_width") == 0) {
			result = parse_bond_line_width (ms, n_words, words);
		}
		else if (strcmp (verb, "molecule") == 0) {
			result = parse_molecule (ms, n_words, words);
		}
		else if (strcmp (verb, "surface") == 0) {
			result = parse_surface (ms, n_words, words);
		}
		else if (strcmp (verb, "read_surface") == 0) {
			result = parse_surface (ms, n_words, words);
		}
		else if (strcmp (verb, "ball_and_stick") == 0) {
			result = parse_ball_and_stick (ms, n_words, words);
		}
		else if (strcmp (verb, "polyhedron") == 0) {
			result = parse_read_polyhedron (ms, n_words, words);
		}
		else if (strcmp (verb, "read_polyhedron") == 0) {
			result = parse_read_polyhedron (ms, n_words, words);
		}
		else if (strcmp (verb, "blank") == 0) {
			result = parse_blank (ms);
		}
		else if (strcmp (verb, "contour") == 0 ||
			strcmp (verb, "contours") == 0) {
			result = parse_contour (ms, n_words, words);
		}
		else if (strcmp (verb, "normals") == 0) {
			result = parse_normals (ms, n_words, words);
		}
		else if (strcmp (verb, "no_clipping") == 0) {
			result = parse_no_clipping (ms);
		}
		else if (strcmp (verb, "solid_shade") == 0) {
			result = parse_solid_shade (ms, n_words, words);
		}
		else if (strcmp (verb, "overlap_hue") == 0) {
			result = parse_overlap_color (ms, n_words, words);
		}
		else if (strcmp (verb, "surface_thickness") == 0) {
			result = parse_surface_thickness (ms, n_words, words);
		}
		else if (strcmp (verb, "define_color") == 0) {
			result = parse_define_color (ms, n_words, words);
		}
		else if (strcmp (verb, "define_ramp") == 0) {
			result = parse_define_ramp (ms, n_words, words);
		}
		else if (strcmp (verb, "atom_coloring") == 0) {
			result = parse_atom_coloring (ms);
		}
		else if (strcmp (verb, "uniform_coloring") == 0) {
			result = parse_uniform_coloring (ms, n_words, words);
		}
		else if (strcmp (verb, "shape_coloring") == 0) {
			result = parse_shape_coloring (ms, n_words, words);
		}
		else if (strcmp (verb, "component_coloring") == 0) {
			result = parse_component_coloring (ms, n_words, words);
		}
		else if (strcmp (verb, "input_coloring") == 0) {
			result =  parse_input_coloring (ms);
		}
		else if (strcmp (verb, "atom_opacity") == 0) {
			result = parse_atom_opacity(ms);
		}
		else if (strcmp (verb, "uniform_opacity") == 0) {
			result = parse_uniform_opacity (ms, n_words, words);
		}
		else if (strcmp (verb, "shape_opacity") == 0) {
			result = parse_shape_opacity (ms, n_words, words);
		}
		else if (strcmp (verb, "component_opacity") == 0) {
			result = parse_component_opacity (ms, n_words, words);
		}
		else if (strcmp (verb, "input_opacity") == 0) {
			result = parse_input_opacity (ms);
		}
		else if (strcmp (verb, "sphere") == 0) {
			result = parse_sphere (ms, n_words, words);
		}
		else if (strcmp (verb, "plane") == 0) {
			result = parse_plane (ms, n_words, words);
		}
		else if (strcmp (verb, "clear") == 0) {
			result = parse_clear (ms, n_words, words);
		}
		else if (strcmp (verb, "tolerance") == 0) {
			result = parse_tolerance (ms, n_words, words);
		}
		else if (strcmp (verb, "connect") == 0) {
			result = parse_connect (ms, n_words, words);
		}
		else if (strcmp (verb, "disconnect") == 0) {
			result = parse_disconnect (ms, n_words, words);
		}
		else if (strcmp (verb, "bond_radius") == 0) {
			result = parse_bond_radius (ms, n_words, words);
		}
		else if (strcmp (verb, "ball_radius") == 0) {
			result = parse_ball_radius (ms, n_words, words);
		}
		else {
			result = pdb_command (command_line);
		}
		if (result == 0 && !error ()) {
			set_error1("msdraw: unrecognizable command");
			set_error2(command_line);
			return;
		}
		if (!result) break;
	}
}

int parse_source (struct msscene *ms, int n_words, char words[MAX_WORD][MAX_NAME])
{
	char scriptfile[MAX_NAME];
	char message[MAXLINE];
	FILE *fp_command;

	strcpy(scriptfile,"");
	if (n_words < 2) {
		set_error1("incomplete source");
		return (0);
	}
	/* read commands from another file */
	strcpy(scriptfile,words[1]);
	fp_command = fopen (scriptfile, "r");
	if (fp_command == NULL) {
		sprintf (message,"msdraw: cannot open nested script file: %s",
			scriptfile);
		set_error1(message);
		return(0);
	}
	read_draw_commands (ms, fp_command);
	fclose (fp_command);
	return (1);
}

int parse_shading_model (struct msscene *ms, int n_words, char words[MAX_WORD][MAX_NAME])
{
	int j;
	if (ms -> current_molecule != NULL) {
		set_error1("too late to declare shading_model");
		return (0);
	}
	if (n_words < 6) {
		set_error1("incomplete shading_model");
		return (0);
	}
	/* depth, ambient, diffuse, specular, exponent */
	for (j = 0; j < 5; j++) {
		ms -> model[j] = atof (words[j + 1]);
	}
	return (1);
}


int parse_light_source (struct msscene *ms, int n_words, char words[MAX_WORD][MAX_NAME])
{
	int j;
	if (ms -> current_molecule != NULL) {
		set_error1("too late to declare light_source");
		return (0);
	}
	if (n_words < 4) {
		set_error1("incomplete light_source");
		return (0);
	}
	/* light source */
	for (j = 0; j < 3; j++) {
		ms -> lsource[j] = atof (words[j + 1]);
	}
	/* error checking */
	if (norm (ms -> lsource) <= 0.0) {
		set_error1("invalid light_source");
		return (0);
	}
	normalize (ms -> lsource);
	return (1);
}

int parse_rotation (struct msscene *ms, int n_words, char words[MAX_WORD][MAX_NAME])
{
	char message[MAXLINE];
	int result, j, k;
	double fangle, t;
	double axis[3];
	if (ms -> current_molecule != NULL) {
		set_error1("too late to declare rotation matrix");
		return (0);
	}
	if (n_words == 5) {
		fangle = atof (words[1]);
		for (k = 0; k < 3; k++)
			axis[k] = atof (words[k+2]);
		/* simple form */
		result = normalize (axis);
		if (!result) {
			set_error1("cannot normalize axis for rotation matrix");
			sprintf (message, "angle = %12.6f; axis = %12.6f %12.6f %12.6f",
				fangle, axis[0], axis[1], axis[2]);
			set_error2(message);
			return (0);
		}
		fangle *= (PI/180.0);
		result = generate_rotation (axis, fangle, ms -> rotation);
		if (!result) {
			set_error1("cannot compute rotation matrix");
			sprintf (message, "angle = %12.6f; axis = %12.6f %12.6f %12.6f",
				fangle, axis[0], axis[1], axis[2]);
			set_error2(message);
			return (0);
		}
		for (j = 0; j < 3; j++) {
			sprintf (message, "%8.3f %8.3f %8.3f",
				ms -> rotation[j][0],
				ms -> rotation[j][1],
				ms -> rotation[j][2]);
			informd(message);
		}

	}
	else if (n_words == 10) {
		/* rotation matrix */
		for (j = 0; j < 3; j++) {
			for (k = 0; k < 3; k++)
				ms -> rotation[j][k] = atof (words[3*j+k+1]);
		}
	}
	else {
		set_error1 ("wrong number of fields for rotation");
		return (0);
	}
	t = triple_product (ms -> rotation[0],
		ms -> rotation[1],
		ms -> rotation[2]);
	if (t <= 0.0) {
		set_error1("rotation matrix not right-handed");
		return (0);
	}
	if (fabs (1.0 - t) > 0.001) {
		set_error1("rotation matrix not orthogonal");
		return (0);
	}
	return (1);
}


int parse_elbow (struct msscene *ms, int n_words, char words[MAX_WORD][MAX_NAME])
{
	char message[MAXLINE];
	if (ms -> current_molecule == NULL) return (0);
	if (n_words < 2) {
		ms -> current_molecule -> elbow = DEFAULT_ELBOW;
	}
	else ms -> current_molecule -> elbow = atof (words[1]);
	/* error checking */
	if (ms -> current_molecule -> elbow < 0.0 ||
		ms -> current_molecule -> elbow > 10.0) {
		sprintf(message,"invalid elbow (%8.3f)",
			ms -> current_molecule -> elbow);
		set_error1 (message);
		return (0);
	}
	return (1);
}


int parse_outer_line_width (struct msscene *ms, int n_words, char words[MAX_WORD][MAX_NAME])
{
	if (n_words < 2) {
		set_error1 ("msdraw: incomplete outer_line_width");
		return (0);
	}
	if (ms -> current_molecule == NULL) return (0);
	ms -> current_molecule -> outer_width = atof (words[1]);
	if (ms -> current_molecule -> outer_width <= 0.0)
		ms -> current_molecule -> outer_width = DEFAULT_OUTER_WIDTH;
	return (1);
}


int parse_inner_line_width (struct msscene *ms, int n_words, char words[MAX_WORD][MAX_NAME])
{
	if (n_words < 2) {
		set_error1 ("msdraw: incomplete inner_line_width");
		return (0);
	}
	if (ms -> current_molecule == NULL) return (0);
	ms -> current_molecule -> inner_width = atof (words[1]);
	if (ms -> current_molecule -> inner_width <= 0.0)
		ms -> current_molecule -> inner_width = DEFAULT_INNER_WIDTH;
	return (1);
}


int parse_cavity_line_width (struct msscene *ms, int n_words, char words[MAX_WORD][MAX_NAME])
{
	if (n_words < 2) {
		set_error1 ("msdraw: incomplete cavity_line_width");
		return (0);
	}
	if (ms -> current_molecule == NULL) return (0);
	ms -> current_molecule -> cavity_width = atof (words[1]);
	if (ms -> current_molecule -> cavity_width <= 0.0)
		ms -> current_molecule -> cavity_width = DEFAULT_CAVITY_WIDTH;
	return (1);
}


int parse_bond_line_width (struct msscene *ms, int n_words, char words[MAX_WORD][MAX_NAME])
{
	if (n_words < 2) {
		set_error1 ("msdraw: incomplete bond_line_width");
		return (0);
	}
	if (ms -> current_molecule == NULL) return (0);
	ms -> current_molecule -> bond_width = atof (words[1]);
	if (ms -> current_molecule -> bond_width <= 0.0)
		ms -> current_molecule -> bond_width = DEFAULT_BOND_WIDTH;
	return (1);
}

int parse_ball_radius (struct msscene *ms, int n_words, char words[MAX_WORD][MAX_NAME])
{
	int atom_set;
	double ball_radius;
	/* ball_radius  real */
	if (n_words < 2) {
		set_error1 ("invalid ball_radius command");
		return (0);
	}
	if (ms -> current_molecule != NULL) {
		ball_radius = atof (words[1]);
		ms -> current_molecule -> ball_radius = ball_radius;
		atom_set = ms -> current_molecule -> atom_set;
		setup_ball_radii (ball_radius, atom_set);
	}
	return (1);
}


int parse_bond_radius (struct msscene *ms, int n_words, char words[MAX_WORD][MAX_NAME])
{
	/* bond_radius  real */
	if (n_words < 2) {
		set_error1 ("invalid bond_radius command");
		return (0);
	}
	if (ms -> current_molecule != NULL)
	ms -> current_molecule -> bond_radius = atof (words[1]);
	return (1);
}

int parse_tolerance (struct msscene *ms, int n_words, char words[MAX_WORD][MAX_NAME])
{
	/* tolerance  real */
	if (n_words < 2) {
		set_error1 ("invalid tolerance command");
		return (0);
	}
	if (ms -> current_molecule != NULL)
	ms -> current_molecule -> tolerance = atof (words[1]);
	return (1);
}

int parse_molecule (struct msscene *ms, int n_words, char words[MAX_WORD][MAX_NAME])
{
	int k, result;
	long atom_set;
	double ball_radius;
	char message[MAXLINE];
	FILE *fpm, *fpr;
	struct molecule *mol;
	struct string_array *pdb_lines;
	struct string_array *radius_lines;

	if (n_words < 2) {
		set_error1 ("missing molecule name");
		return (0);
	}
	/* check for uniqueness */
	for (mol = ms -> head_molecule; mol != NULL; mol = mol -> next) {
		if (strcmp (mol -> name, words[1]) == 0) {
			set_error1 ("parse_molecule: duplicate molecule name");
			return (0);
		}
	}
	fpm = fpr = NULL;
	if (n_words >= 3) {
		/* read pdb file */
		fpm = fopen (words[2], "r");
		if (fpm == NULL) {
			sprintf (message, "msdraw: cannot open molecule file: %s", words[2]);
			set_error1 (message);
			return (0);
		}
		if (n_words >= 4) {
			fpr = fopen (words[3], "r");
			if (fpr == NULL) {
				sprintf (message, "msdraw: cannot open radii file: %s", words[3]);
				set_error1 (message);
				return (0);
			}
		}
		else fpr = NULL;
		mol = read_molecule (ms, words[1], words[2], DEFAULT_ALPHA, fpm, fpr, NULL, NULL);
		if (error()) {
			return (0);
		}
	}
	else { /* set up molecule with no atoms */
		mol = new_molecule (0L, 0L, NULL, NULL);
		if (mol == NULL) {
			set_error1("memory allocation fails");
			return (0);
		}
		if (ms -> head_molecule == NULL)
			ms -> head_molecule = mol;
		else ms -> tail_molecule -> next = mol;
		ms -> tail_molecule = mol;
		ms -> current_molecule = mol;
		strcpy (mol -> name, words[1]);
		mol -> head_surface = NULL;
		mol -> tail_surface = NULL;
		for (k = 0; k < 3; k++) {
			mol -> center[k] = 0.0;
		}
		atom_set = named_set (mol -> name, ATOM);
		if (error()) return (0);
		if (atom_set == 0) {
			set_error1 ("no set for molecule name");
			return (0);
		}
		mol -> atom_set = atom_set;
	}
	return (1);
}

int parse_connect (struct msscene *ms, int n_words, char words[MAX_WORD][MAX_NAME])
{
	/* connect setname */
	if (n_words < 1) {
		set_error1 ("invalid connect command");
		return (0);
	}
	if (ms -> current_molecule == NULL) return (0);
	connect_command (words[1], ms -> current_molecule -> tolerance);
	if (error ()) return (0);
	return (1);
}

int parse_disconnect (struct msscene *ms, int n_words, char words[MAX_WORD][MAX_NAME])
{
	/* disconnect set1 set2 */
	if (n_words != 2) {
		set_error1 ("invalid disconnect command");
		return (0);
	}
	if (ms -> current_molecule == NULL) return (0);
	disconnect_command (words[1], words[2]);
	if (error ()) return (0);
	return (1);
}


int parse_surface (struct msscene *ms, int n_words, char words[MAX_WORD][MAX_NAME])
{
	char message[MAXLINE];
	FILE *fpsurf;					/* surface file */
	if (ms -> current_molecule == NULL) {
		set_error1("surface for undeclared molecule");
		return (0);
	}
	if (n_words < 2) {
		set_error1 ("msdraw: incomplete parse_surface");
		return (0);
	}
	/* open surface file */
	fpsurf = fopen (words[1], "rb");
	if (fpsurf == NULL) {
		sprintf(message, "cannot open pqms file: %s", words[1]);
		set_error1 (message);
		return (0);
	}
	/* read pqms file */
	ms -> current_molecule -> current_surface = surf_input (fpsurf);
	if (ms -> current_molecule -> current_surface == NULL) return (0);
	/* close input file */
	fclose (fpsurf);
	if (error()) return (0);

	if (ms -> current_molecule -> head_surface == NULL)
		ms -> current_molecule -> head_surface = ms -> current_molecule -> current_surface;
	else
		ms -> current_molecule -> tail_surface -> next = ms -> current_molecule -> current_surface;
	ms -> current_molecule -> tail_surface = ms -> current_molecule -> current_surface;
	ms -> current_molecule -> n_surface++;

	setup_surface_bounds (ms -> current_molecule -> current_surface);
	if (error()) return (0);
	return (1);
}


int parse_ball_and_stick (struct msscene *ms, int n_words, char words[MAX_WORD][MAX_NAME])
{
	int h;
	long bond_set;
	double red;
	double green;
	double blue;
	char color_name[MAX_NAME];
	if (ms -> current_molecule == NULL) {
		set_error1("ball_and_stick for undeclared molecule");
		return (0);
	}
	bond_set = select_bonds (ms -> current_molecule -> atom_set);
	ms -> current_molecule -> bond_set = bond_set;
	ms -> current_molecule -> current_surface = make_bas (ms -> current_molecule);
	if (error()) return (0);
	setup_surface_bounds (ms -> current_molecule -> current_surface);
	if (error()) return (0);
	if (ms -> current_molecule -> head_surface == NULL)
		ms -> current_molecule -> head_surface = ms -> current_molecule -> current_surface;
	else
		ms -> current_molecule -> tail_surface -> next = ms -> current_molecule -> current_surface;
	ms -> current_molecule -> tail_surface = ms -> current_molecule -> current_surface;
	ms -> current_molecule -> n_surface++;
	setup_linewidth (ms -> current_molecule -> current_surface,
		ms -> current_molecule -> outer_width, ms -> current_molecule -> inner_width,
		ms -> current_molecule -> cavity_width, ms -> current_molecule -> bond_width);
	if (error()) return (0);
	return (1);
}


int parse_read_polyhedron (struct msscene *ms, int n_words, char words[MAX_WORD][MAX_NAME])
{
	int h;
	char color_name[MAX_NAME];
	char message[MAXLINE];
	double red;
	double green;
	double blue;
	FILE *fp_polyhedron;			/* polyhedron file */
	if (ms -> current_molecule == NULL) {
		set_error1 ("parse_polyhedron: no molecule declared");
		return (0);
	}
	/* open polyhedron file */
	fp_polyhedron = fopen (words[1], "r");
	if (fp_polyhedron == NULL) {
		sprintf (message, "cannot open vet file: %s", words[1]);
		set_error1(message);
		return (0);
	}
	strcpy (ms -> current_molecule -> name, words[1]);
	/* read polyhedron file */
	ms -> current_molecule -> current_surface = read_phn_surface (fp_polyhedron);
	/* close input file */
	fclose (fp_polyhedron);
	if (error()) return (0);
	if (ms -> current_molecule -> head_surface == NULL)
		ms -> current_molecule -> head_surface = ms -> current_molecule -> current_surface;
	else
		ms -> current_molecule -> tail_surface -> next = ms -> current_molecule -> current_surface;
	ms -> current_molecule -> tail_surface = ms -> current_molecule -> current_surface;
	ms -> current_molecule -> n_surface++;
	do_bounds (ms -> current_molecule -> current_surface);
	if (error()) return (0);
	do_axes (ms -> current_molecule -> current_surface);
	if (error()) return (0);
	setup_surface_bounds (ms -> current_molecule -> current_surface);
	if (error()) return (0);
	ms -> current_molecule -> current_surface -> function[0] = 'z';
	/* min_max (ms -> current_molecule); */
	setup_linewidth (ms -> current_molecule -> current_surface,
		ms -> current_molecule -> outer_width, ms -> current_molecule -> inner_width,
		ms -> current_molecule -> cavity_width, ms -> current_molecule -> bond_width);
	if (error()) return (0);
	return (1);
}



/* later:  to we really need the edges ? */

struct surface *read_phn_surface (FILE *fp_polyhedron)
{
	long n_vertex, n_edge, n_triangle;
	struct surface *phn;

	phn = read_vet (fp_polyhedron);
	if (phn == NULL) {
		set_error2 ("(read_phn_surface): read_vet returns null");
		return (NULL);
	}
	if (phn == NULL) {
		set_error2 ("(read_phn_bunch): no polyhedron returned");
		return(NULL);
	}
	fclose (fp_polyhedron);
	phn -> type = PHN_SURFACE;
	/* set to input coloring (set by trb) */
	phn -> scheme = define_scheme (INPUT_COLORING, 1, 0, 0.0, 0.0);
	/* store counts of vertices and triangles */
	n_vertex = phn -> n_phnvtx;
	n_edge = phn -> n_phnedg;
	n_triangle = phn -> n_phntri;
	if (n_vertex <= 0) {
		set_error1 ("(read_phn_surface): bad n_vertex");
		return(NULL);
	}
	if (n_edge <= 0) {
		set_error1 ("(read_phn_surface): bad n_edge");
		return(NULL);
	}
	if (n_triangle <= 0) {
		set_error1 ("(read_phn_surface): bad n_triangle");
		return(NULL);
	}
	phn -> n_component = 1;
	phn -> clipping = 1;
	setup_triangles (phn);
	if (error()) return (NULL);;
	polygonize_polyhedron (phn);
	if (error()) return(NULL);
	do_polyhedron_bounds (phn);
	if (error()) return(NULL);
	return (phn);
}


/* setup each triangle */

void setup_triangles (struct surface *phn)
{
	int j, k, orn;
	long t;
	double dvc;
	double pt1[3], pt2[3];
	double vc[3][3];
	struct phntri *tri;
	struct phnedg *e;
	struct phnvtx *vtx;

	/* n_triangle = phn -> n_face; */
	if (phn -> phntri_handles == NULL) return;
	for (t = 0; t < phn -> n_phntri; t++) {
		tri = num2phntri (phn, t+1);
		if (tri == NULL) {
			set_error1 ("setup_triangles: null triangle pointer");
			return;
		}
		/* compute triangle center */
		for (k = 0; k < 3; k++)
			tri -> center[k] = 0.0;
		for (j = 0; j < 3; j++) {
			e = tri -> edg[j];
			orn = tri -> orn[j];
			vtx = e -> pvt[orn];
			for (k = 0; k < 3; k++) {
				tri -> center[k] += vtx -> center[k];
				vc[j][k] = vtx -> center[k];
			}
		}
		for (k = 0; k < 3; k++)
			tri -> center[k] /= 3;
		/* compute radius of circumscribing circle */
		tri -> radius = 0.0;
		for (j = 0; j < 3; j++) {
			for (k = 0; k < 3; k++) {
				pt1[k] = vc[j][k];
				pt2[k] = tri -> center[k];
			}
			dvc = distance (pt1, pt2);
			if (dvc > tri -> radius) tri -> radius = dvc;
		}
	}
}

int parse_blank (struct msscene *ms)
{
	if (ms -> current_molecule == NULL) {
		inform ("msdraw: premature blank command ignored");
		return (1);
	}
	ms -> current_molecule -> blank = 1;
	return (1);
}


int parse_no_clipping (struct msscene *ms)
{
	if (ms -> current_molecule -> current_surface == NULL) {
		set_error1("no_clipping for undeclared surface");
		return (0);
	}
	ms -> current_molecule -> current_surface -> clipping = 0;
	return (1);
}

int parse_solid_shade (struct msscene *ms, int n_words, char words[MAX_WORD][MAX_NAME])
{
	int grey_level;
	if (ms -> current_molecule -> current_surface == NULL) {
		set_error1("solid_shade for undeclared surface");
		return (0);
	}
	if (n_words < 2) grey_level = 127;
	else grey_level = atoi (words[1]);
	if (grey_level < 0) grey_level = 0;
	if (grey_level > 255) grey_level = 255;
	ms -> current_molecule -> current_surface -> solid_shade = grey_level;
	return (1);
}

int parse_overlap_color (struct msscene *ms, int n_words, char words[MAX_WORD][MAX_NAME])
{
	int hue;
	if (n_words > 1)
		hue = get_color_number (ms -> table, words[1]);
	else hue = 1;
	ms -> overlap_hue = hue;
	return (1);
}

int parse_surface_thickness (struct msscene *ms, int n_words, char words[MAX_WORD][MAX_NAME])
{
	double surface_thickness;
	if (ms -> current_molecule -> current_surface == NULL) {
		set_error1("surface_thickness for undeclared surface");
		return (0);
	}
	if (n_words < 2) surface_thickness = DEFAULT_THICKNESS;
	else surface_thickness = atof (words[1]);
	if (surface_thickness < MINIMUM_THICKNESS) surface_thickness = MINIMUM_THICKNESS;
	if (surface_thickness > MAXIMUM_THICKNESS) surface_thickness = MAXIMUM_THICKNESS;
	ms -> current_molecule -> current_surface -> surface_thickness = surface_thickness;
	return (1);
}


int parse_normals (struct msscene *ms, int n_words, char words[MAX_WORD][MAX_NAME])
{
	double length;
	double red;
	double green;
	double blue;
	char color_name[MAX_NAME];
	if (ms -> current_molecule == NULL) {
		set_error1 ("msdraw: normals for undeclared molecule");
		return (0);
	}
	length = atof (words[1]);
	if (n_words == 3)
		strcpy (color_name, words[2]);
	else strcpy (color_name, "black");
	name_to_rgb (ms -> table, color_name, &red, &green, &blue);
	/* make fuzz bunch */
	ms -> current_molecule -> current_surface =  polyhedron_normals (ms, length);
	if (error()) return (0);
	if (ms -> current_molecule -> current_surface == NULL) {
		set_error1 ("msdraw: polyhedron_normals returns null");
		return (0);
	}

	ms -> current_molecule -> tail_surface -> next = ms -> current_molecule -> current_surface;
	ms -> current_molecule -> tail_surface = ms -> current_molecule -> current_surface;
	ms -> current_molecule -> n_surface++;

	setup_surface_bounds (ms -> current_molecule -> current_surface);
	if (error()) return (0);
	setup_linewidth (ms -> current_molecule -> current_surface,
		ms -> current_molecule -> outer_width, ms -> current_molecule -> inner_width,
		ms -> current_molecule -> cavity_width, ms -> current_molecule -> bond_width);
	if (error()) return (0);
	return (1);
}


/* contour function from to ramp */

int parse_contour (struct msscene *ms, int n_words, char words[MAX_WORD][MAX_NAME])
{
	long w, h, hue;
	long n_levels, nmaterial, n_ramp;
	double level;
	double red;
	double green;
	double blue;
	double from, to;
	char message[MAXLINE];
	struct object_scheme *cs;
	struct surface *obj, *ctr, *phn;
	double floats[MAX_LEVELS];
	long colors[MAX_LEVELS];
	struct color_ramp *ramp;
	struct material_table *table;
	struct molecule *mol;

	mol = ms -> current_molecule;
	if (mol == NULL) {
		set_error1 ("msdraw: contours for undeclared molecule");
		return (0);
	}
	if (n_words < 5) {
		set_error1 ("msdraw: contour command has missing arguments");
		return (0);
	}
	table = ms -> table;
	if (table == NULL) {
		set_error1 ("msdraw: null material table");
		return (0);
	}
	nmaterial = table -> nmaterial;
	obj = mol -> current_surface;
	phn = polyhedron_bunch (mol);
	if (phn == NULL) {
		set_error1 ("msdraw: no polyhedron to contour for this molecule");
		return (0);
	}
	if (strcmp (words[1], "u") == 0 ||
		strcmp (words[1], "v") == 0 ||
		strcmp (words[1], "w") == 0 ||
		strcmp (words[1], "x") == 0 ||
		strcmp (words[1], "y") == 0 ||
		strcmp (words[1], "z") == 0) {
		strcpy (phn -> function, words[1]);
		min_max (mol);
	}
	else {
		set_error1 ("msdraw: invalid function definition");
		return (0);
	}
	from = atof (words[2]);
	to = atof (words[3]);
	ramp = lookup_ramp (ms, words[4]);
	if (error ()) return (0);
	if (ramp == NULL) {
		set_error1 ("msdraw: invalid ramp name");
		return (0);
	}
	n_ramp = ramp -> n_ramp;
	if (n_ramp < 2) {
		set_error1 ("msdraw: invalid ramp");
		return (0);
	}
	n_levels = n_ramp;
	if (n_levels >= MAX_LEVELS) {
		set_error1 ("msdraw: contour command has too many arguments");
		return (0);
	}
	/* rendered image contours = vertex coloring */
	cs = phn -> scheme;
	cs -> color_type = VERTEX_COLORING;
	cs -> minval = from;
	cs -> maxval = to;
	cs -> ramp = ramp;

	for (h = 0; h < n_levels; h++)
		floats[h] = from + h * (to - from)/(n_levels-1) ;
	for (h = 0; h < n_levels; h++) {
		level = floats[h];
		hue = ramp -> colors[h];
		ramp -> vertex_colors[h][0] = hue;
		ramp -> vertex_colors[h][1] = hue;
		ramp -> vertex_ranges[h][0] = level;
		ramp -> vertex_ranges[h][1] = MS_INFINITY;
		if (h > 0) ramp -> vertex_ranges[h-1][1] = level;
	}
	if (ms -> plot_format != 0 || ms -> vector_format != 0) {
		/* plotting */
		for (h = 0; h < n_levels; h++) {
			level = floats[h];
			hue = ramp -> colors[h];
			sprintf (message, "%8.3f contour hue %2d", level, hue);
			informd (message);
			number_to_rgb (table, hue, &red, &green, &blue);
			/* make contour bunch */
			ctr =  contour_polyhedron (ms, level);
			if (error()) return (0);
			if (ctr == NULL) {
				set_error1 ("msdraw: contour_polyhedron returns null");
				return (0);
			}
			mol -> current_surface =  ctr;
			mol -> tail_surface -> next = ctr;
			mol -> tail_surface = ctr;
			mol -> n_surface++;
			setup_surface_bounds (ctr);
			if (error()) return (0);
			setup_linewidth (ctr, mol -> outer_width, mol -> inner_width,
				mol -> cavity_width, mol -> bond_width);
			if (error()) return (0);
			ctr -> scheme = define_scheme (UNIFORM_COLORING, hue, -1, 0.0, 0.0);
			sprintf (message,"%8.3f = %2s %15s %8ld vertices; %8ld edges; %8ld contours",
				level, words[1], table -> name[hue],
				ctr -> n_phnvtx, ctr -> n_phnedg, ctr -> n_phnctr);
			inform(message);
		}
	}
	return (1);
}

int parse_define_color (struct msscene *ms, int n_words, char words[MAX_WORD][MAX_NAME])
{
	int h, nmaterial;
	double red, green, blue;
	struct material_table *table;

	if (n_words < 5) {
		set_error1 ("msdraw: define_color command has missing arguments");
		return (0);
	}
	table = ms -> table;
	if (table == NULL) {
		set_error1 ("msdraw: null material table");
		return (0);
	}
	red = atof (words[2]); green = atof (words[3]); blue = atof (words[4]);
	nmaterial = table -> nmaterial;
	for (h = 0; h < nmaterial; h++) {
		if (strcmp (words[1], table -> name[h]) == 0) {
			/* found, change */
			table -> red[h] = red;
			table -> green[h] = green;
			table -> blue[h] = blue;
			return (1);
		}
	}
	/* not found, add new color */
	add_material (table, red, green, blue, words[1]);
	if (error ()) return (0);
	return (1);
}

int parse_define_ramp (struct msscene *ms, int n_words, char words[MAX_WORD][MAX_NAME])
{
	long w, h, hue;
	long n_color, nmaterial;
	double level;
	char message[MAXLINE];
	long divisions[MAX_LEVELS];
	long colors[MAX_LEVELS];
	struct color_ramp *ramp;
	struct material_table *table;

	if (n_words < 5) {
		set_error1 ("msdraw: define_ramp command has missing arguments");
		return (0);
	}
	table = ms -> table;
	if (table == NULL) {
		set_error1 ("msdraw: null material table");
		return (0);
	}
	nmaterial = table -> nmaterial;
	for (h = 0; h < MAX_LEVELS; h++) {
		divisions[h] = 0;
		colors[h] = 0;
	}
	n_color = 0;
	for (h = 0; h < MAX_LEVELS; h++) {
		w = 2 * h + 2;
		if (w >= n_words) break;
		colors[h] = get_color_number (table, words[w]);
		if (w < n_words) divisions[h] = atoi (words[w+1]);
		n_color++;
	}
	if (n_color >= MAX_LEVELS) {
		set_error1 ("msdraw: define_ramp command has too many arguments");
		return (0);
	}

	/* check values */
	for (h = 0; h < n_color; h++) {
		if (colors[h] < 0 || colors[h] >= nmaterial) {
			set_error1 ("msdraw: define_ramp command has invalid color");
			return (0);
		}
		if (h < n_color - 1) {
			if (divisions[h] < 1 || divisions[h] > MAX_LEVELS - 3) {
				set_error1 ("msdraw: define_ramp has invalid # intermediate colors");
				return (0);
			}
		}
	}
	/* create a ramp */
	ramp = create_ramp (ms, n_color, colors, divisions, words[1]);
	return (1);
}

int parse_atom_coloring (struct msscene *ms)
{
	struct surface *surface;
	if (ms -> current_molecule == NULL) {
		set_error1("atom_coloring for undeclared molecule");
		return (0);
	}
	/* color current graphical object of the current molecule */
	for (surface = ms -> current_molecule -> head_surface; surface != NULL; surface = surface -> next) {
		if (surface != ms -> current_molecule -> current_surface) continue;
		set_atom_hues (ms, ms -> current_molecule, surface);
		if (error()) return (0);
	}
	return (1);
}


int parse_uniform_coloring (struct msscene *ms, int n_words, char words[MAX_WORD][MAX_NAME])
{
	int i;
	int twohues[2];
	if (ms -> current_molecule -> current_surface == NULL) {
		set_error1("uniform_coloring for undeclared surface");
		return (0);
	}
	if (n_words < 2) {
		set_error1("(uniform_coloring): incomplete line");
		set_error2("format: uniform_coloring   outside_hue   [inside_hue]");
		return (0);
	}
	twohues[0] = get_color_number (ms -> table, words[1]);
	if (n_words >= 3) twohues[1] = get_color_number (ms -> table, words[2]);
	else twohues[1] = twohues[0];

	/* error checking */
	for (i = 0; i < 2; i++) {
		if (twohues[i] <= 0) {
			set_error1("invalid uniform hue");
			return (0);
		}
		if (twohues[i] >= ms -> table -> nmaterial) {
			set_error1("invalid uniform hue");
			return (0);
		}
	}
	ms -> current_molecule -> current_surface -> scheme -> color_type = UNIFORM_COLORING;
	for (i = 0; i < 2; i++)
		ms -> current_molecule -> current_surface -> scheme -> uniform_colors[i] =
			twohues[i];
	return (1);
}


int parse_shape_coloring (struct msscene *ms, int n_words, char words[MAX_WORD][MAX_NAME])
{
	int i;
	int sixhues[6];
	if (ms -> current_molecule -> current_surface == NULL) {
		set_error1("shape_coloring for undeclared surface");
		return (0);
	}
	switch (n_words-1) {
	case 2:
		sixhues[0] = get_color_number (ms -> table, words[1]);
		sixhues[1] = get_color_number (ms -> table, words[2]);
		sixhues[2] = get_color_number (ms -> table, words[2]);
		sixhues[3] = get_color_number (ms -> table, words[1]);
		sixhues[4] = get_color_number (ms -> table, words[2]);
		sixhues[5] = get_color_number (ms -> table, words[2]);
		break;
	case 3:
		sixhues[0] = get_color_number (ms -> table, words[1]);
		sixhues[1] = get_color_number (ms -> table, words[2]);
		sixhues[2] = get_color_number (ms -> table, words[3]);
		sixhues[3] = get_color_number (ms -> table, words[1]);
		sixhues[4] = get_color_number (ms -> table, words[2]);
		sixhues[5] = get_color_number (ms -> table, words[3]);
		break;
	case 6:
		sixhues[0] = get_color_number (ms -> table, words[1]);
		sixhues[1] = get_color_number (ms -> table, words[2]);
		sixhues[2] = get_color_number (ms -> table, words[3]);
		sixhues[3] = get_color_number (ms -> table, words[4]);
		sixhues[4] = get_color_number (ms -> table, words[5]);
		sixhues[5] = get_color_number (ms -> table, words[6]);
		break;
	default:
		set_error1("shape_coloring: not 2, 3 or 6 arguments");
		return (0);
	}
	/* error checking */
	for (i = 0; i < 6; i++) {
		if (sixhues[i] <= 0) {
			set_error1("invalid shape hue");
			return (0);
		}
		if (sixhues[i] >= ms -> table -> nmaterial) {
			set_error1("invalid shape hue");
			return (0);
		}
	}
	ms -> current_molecule -> current_surface -> scheme -> color_type = SHAPE_COLORING;
	for (i = 0; i < 6; i++)
		ms -> current_molecule -> current_surface -> scheme -> shape_colors[i] =
			sixhues[i];
	return (1);
}


int parse_component_coloring (struct msscene *ms, int n_words, char words[MAX_WORD][MAX_NAME])
{
	int i, j;
	int fourhues[4];
	if (ms -> current_molecule -> current_surface == NULL) {
		set_error1("component_coloring for undeclared surface");
		return (0);
	}
	if (n_words < 4) {
		set_error1("component_coloring: incomplete line");
		set_error2("format: component_coloring   surface_outside   surface_inside  cavity");
		return (0);
	}
	for (j = 0; j < 3; j++)
		fourhues[j] = get_color_number (ms -> table, words[j+1]);
	fourhues[3] = fourhues[2];
	/* error checking */
	for (i = 0; i < 4; i++) {
		if (fourhues[i] <= 0) {
			set_error1("invalid component hue");
			return (0);
		}
		if (fourhues[i] >= ms -> table -> nmaterial) {
			set_error1("invalid component hue");
			return (0);
		}
	}
	ms -> current_molecule -> current_surface -> scheme -> color_type = COMPONENT_COLORING;
	for (i = 0; i < 4; i++)
		ms -> current_molecule -> current_surface -> scheme -> component_colors[i] =
			fourhues[i];
	return (1);
}

int parse_input_coloring (struct msscene *ms)
{
	struct surface *surface;
	if (ms -> current_molecule == NULL) {
		set_error1("input_coloring for undeclared molecule");
		return (0);
	}
	/* color current graphical object of the current molecule */
	for (surface = ms -> current_molecule -> head_surface; surface != NULL; surface = surface -> next) {
		if (surface != ms -> current_molecule -> current_surface) continue;
		ms -> current_molecule -> current_surface -> scheme -> color_type = INPUT_COLORING;
	}
	return (1);
}


int parse_input_opacity (struct msscene *ms)
{
	int i;
	struct surface *surface;
	if (ms -> current_molecule == NULL) {
		set_error1("input_opacity for undeclared molecule");
		return (0);
	}
	/* opacify current graphical object of the current molecule */
	for (surface = ms -> current_molecule -> head_surface; surface != NULL; surface = surface -> next) {
		if (surface != ms -> current_molecule -> current_surface) continue;
		ms -> current_molecule -> current_surface -> scheme -> opacity_type = INPUT_OPACITY;
	}
	for (i = 0; i < 2; i++)
		ms -> current_molecule -> current_surface -> scheme -> uniform_opacities[i] = (1-i);
	ms -> translucency = 1;
	return (1);
}

int parse_atom_opacity (struct msscene *ms)
{
	if (ms -> current_molecule -> current_surface == NULL) {
		set_error1("atom_opacity for undeclared surface");
		return (0);
	}
	set_atom_opacities (ms -> current_molecule, ms -> current_molecule -> current_surface);
	if (error()) return (0);
	ms -> translucency = 1;
	return (1);
}


int parse_uniform_opacity (struct msscene *ms, int n_words, char words[MAX_WORD][MAX_NAME])
{
	int i, j;
	double twoopacities[2];
	if (ms -> current_molecule -> current_surface == NULL) {
		set_error1("uniform_opacity for undeclared surface");
		return (0);
	}
	if (n_words < 2) {
		set_error1("(uniform opacities): incomplete line");
		set_error2("format: uniform_opacity   real_number");
		return (0);
	}
	for (j = 0; j < 2; j++)
		twoopacities[j] = get_opacity_number (words[1+j]);
	/* error checking */
	for (i = 0; i < 2; i++) {
		if (twoopacities[i] < 0.0) {
			set_error1("invalid uniform opacity");
			return (0);
		}
		if (twoopacities[i] > 1.0) {
			set_error1("invalid uniform opacity");
			return (0);
		}
	}
	if (ms -> current_molecule -> current_surface -> scheme == NULL) {
		set_error1("uniform_opacity for null opacity scheme");
		return (0);
	}
	ms -> current_molecule -> current_surface -> scheme -> opacity_type = UNIFORM_OPACITY;
	for (i = 0; i < 2; i++)
		ms -> current_molecule -> current_surface -> scheme -> uniform_opacities[i] =
			twoopacities[i];
	ms -> translucency = 1;
	return (1);
}


int parse_shape_opacity (struct msscene *ms, int n_words, char words[MAX_WORD][MAX_NAME])
{
	int i;
	double sixopacities[6];
	if (ms -> current_molecule -> current_surface == NULL) {
		set_error1("shape_opacity for undeclared surface");
		return (0);
	}
	if (n_words < 3) {
		set_error1("shape_opacity: not 2, 3 or 6 arguments");
		return (0);
	}
	switch (n_words-1) {
	case 2:
		sixopacities[0] = get_opacity_number (words[1]);
		sixopacities[1] = get_opacity_number (words[2]);
		sixopacities[2] = get_opacity_number (words[2]);
		sixopacities[3] = get_opacity_number (words[1]);
		sixopacities[4] = get_opacity_number (words[2]);
		sixopacities[5] = get_opacity_number (words[2]);
		break;
	case 3:
		sixopacities[0] = get_opacity_number (words[1]);
		sixopacities[1] = get_opacity_number (words[2]);
		sixopacities[2] = get_opacity_number (words[3]);
		sixopacities[3] = get_opacity_number (words[1]);
		sixopacities[4] = get_opacity_number (words[2]);
		sixopacities[5] = get_opacity_number (words[3]);
		break;
	case 6:
		sixopacities[0] = get_opacity_number (words[1]);
		sixopacities[1] = get_opacity_number (words[2]);
		sixopacities[2] = get_opacity_number (words[3]);
		sixopacities[3] = get_opacity_number (words[4]);
		sixopacities[4] = get_opacity_number (words[5]);
		sixopacities[5] = get_opacity_number (words[6]);
		break;
	default:
		set_error1("shape_opacity: not 2, 3 or 6 arguments");
		return (0);
	}
	/* error checking */
	for (i = 0; i < 6; i++) {
		if (sixopacities[i] < 0.0) {
			set_error1("invalid shape opacity");
			return (0);
		}
		if (sixopacities[i] > 1.0) {
			set_error1("invalid shape opacity");
			return (0);
		}
	}
	ms -> current_molecule -> current_surface -> scheme -> opacity_type = SHAPE_OPACITY;
	for (i = 0; i < 6; i++)
		ms -> current_molecule -> current_surface -> scheme -> shape_opacities[i] =
			sixopacities[i];
	ms -> translucency = 1;
	return (1);
}


int parse_component_opacity (struct msscene *ms, int n_words, char words[MAX_WORD][MAX_NAME])
{
	int i;
	double fouropacities[4];
	if (ms -> current_molecule -> current_surface == NULL) {
		set_error1("component_opacity for undeclared surface");
		return (0);
	}
	if (n_words < 3) {
		set_error1("component_opacity: incomplete line");
		set_error2("format: component_opacity surface_outside   surface_inside  cavity");
		return (0);
	}
	fouropacities[0] = get_opacity_number (words[1]);
	fouropacities[1] = get_opacity_number (words[2]);
	fouropacities[2] = get_opacity_number (words[3]);
	fouropacities[3] = get_opacity_number (words[3]);
	/* error checking */
	for (i = 0; i < 4; i++) {
		if (fouropacities[i] < 0.0) {
			set_error1("invalid shape opacity");
			return (0);
		}
		if (fouropacities[i] > 1.0) {
			set_error1("invalid shape opacity");
			return (0);
		}
	}
	ms -> current_molecule -> current_surface -> scheme -> opacity_type = COMPONENT_OPACITY;
	for (i = 0; i < 4; i++)
		ms -> current_molecule -> current_surface -> scheme -> component_opacities[i] =
			fouropacities[i];
	ms -> translucency = 1;
	return (1);
}

/* set atom colors (hues) */

void set_atom_hues (struct msscene *ms, struct molecule *mol, struct surface *obj)
{
	int j;
	long i, b, atom1;
	char message[MAXLINE];
	long *atom_colors;
	int hue;
	long n_atom;
	long atom_set, atom_number;
	struct phnedg *edg;

	obj -> scheme -> color_type = ATOM_COLORING;
	n_atom = mol -> n_atom;
	if (n_atom <= 0) {
		n_atom = obj -> n_atom;
	}
	if (n_atom == 0) {
		sprintf (message, "(set_atom_hues): invalid number of atoms %ld", n_atom);
		set_error1 (message);
		return;
	}
	/* allocate memory for atom colors */
	atom_colors = allocate_longs (n_atom * 4, 0, ATOM_COLORS);
	if (atom_colors == NULL) {
		set_error1 ("(set_atom_hues): memory allocation failure");
		return;
	}
	obj -> scheme -> atom_colors = atom_colors;
	sprintf (message, "%8d atoms have had their coloring assigned", n_atom);
	informd (message);

	atom_set = mol -> atom_set;
	for (atom_number = init_for (atom_set), i = 0; atom_number != 0;
		atom_number = next_for (atom_set), i++) {
		if (i >= obj -> n_atom) break;
		hue = get_atom_color (atom_number);
		*(atom_colors+4*i) = hue;
		*(atom_colors+4*i+1) = hue;
		*(atom_colors+4*i+2) = hue;
		*(atom_colors+4*i+3) = hue;
		sprintf (message, "atom %4d colors: %2d %2d %2d %2d", i + 1, hue, hue, hue, hue);
		informd (message);
		/* contact, reentrant */
		/* error checking */
		for (j = 0; j < 4; j++) {
			if (*(atom_colors+4*i+j) < 0) {
				set_error1 ("(set_atom_hues): invalid atom hue");
				return;
			}
			if (*(atom_colors+4*i+j) >= ms -> table -> nmaterial) {
				set_error1 ("(set_atom_hues): invalid atom hue");
				return;
			}
		}
	}
	if (obj -> type == BAS_SURFACE) {
		/* bond colors for stick model */
		for (b = 0; b < obj -> n_phnedg; b++) {
			edg = *(obj -> phnedg_handles + b);
			atom1 = edg -> atm;
			hue = determine_hue (table, obj -> scheme, atom1, 1, 1, 0, 0, 0.0);
			edg -> hue = hue;
		}
	}
}

