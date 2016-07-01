/*
 * Molecular Surface Package

 * Copyright 1986 by Michael L. Connolly
 * All Rights Reserved

 * January 1, 2002

 */

#include "ms.h"
#include "msenum.h"
#include "mslimit.h"
#include "msdefault.h"
#include "msstruct.h"
#include "msglobal.h"
#include "msfunction.h"

double detopac (int atm, int inner, int comp, int shape, struct object_scheme *os, double input_opacity)
{
	int ind;
	int ccomp;
	char message[MAXLINE];
	double opacity;

	if (os == NULL) {
		sprintf (message, "detopac: render encounters null opacity");
		set_error1 (message);
		return(0.0);
	}

	opacity = 1.0; /* default */
	if (os -> opacity_type == ATOM_OPACITY) {
		if (shape > 2) shape = 2;
		ind = 4 * (atm - 1) + 2 * inner + shape - 1;
		opacity = *(os -> atom_opacities + ind);
	}
	else if (os -> opacity_type == SHAPE_OPACITY) {
		ind = 3 * inner + shape - 1;
		opacity = os -> shape_opacities[ind];
	}
	else if (os -> opacity_type == UNIFORM_OPACITY) {
		ind = inner;
		opacity = os -> uniform_opacities[ind];
	}
	else if (os -> opacity_type == INPUT_OPACITY) {
		opacity = input_opacity;
	}
	else if (os -> opacity_type == COMPONENT_OPACITY) {
		ccomp = (comp >= 2) ? 2 : 1;
		ind = 2 * (ccomp - 1) + inner;
		opacity = os -> component_opacities[ind];
	}
	else {
		set_error1 ("detopac: render encounters invalid opacity scheme");
		return(0.0);
	}

	if (opacity < 0.0 || opacity > 1.0) {
		sprintf (message, "detopac: render encounters invalid opacity = %8.3f", opacity);
		return(0.0);
	}
	return (opacity);
}


int opaque (int x, int y, double opacity)
{	
	/* screen-door transparency */
	if (opacity <= 0.125) return (0);
	else if (opacity <= 0.375) return (!((x+2*y) % 4));
	else if (opacity <= 0.625) return ((x+y) % 2);
	else if (opacity <= 0.875) return ((x+2*y) % 4);
	else return (1);
}


void setup_color_table (struct msscene *ms)
{
	long h, s, idx, i, j, this_hue, pidx, this_alpha;
	unsigned char ucr, ucg, ucb;
	double r, g, b;
	double fnshades, fraction;
	long nmaterial, nonzero, nbackground;
	char message[MAXLINE];
	struct material_table *table;
	
	table = ms -> table;
	nmaterial = table -> nmaterial;
	for (i = 0; i < MAX_MATERIAL; i++)
		table -> hue_count[i] = 0;
	nbackground = 0;
	for (i = 0; i < ms -> vertical; i++) 
		for (j = 0; j < ms -> horizontal; j++) {
			/* compute index into depth buffer arrays */
			pidx = i * ms -> horizontal + j;
			this_alpha = *(ms -> db -> alphas+pidx);
			if (this_alpha == 0) {
				nbackground++;
				continue;
			}
			this_hue = *(ms -> db -> hues+pidx);
			if (this_hue < 0) this_hue = 0;
			if (this_hue >= table -> nmaterial) this_hue = ms -> table -> nmaterial - 1;
			table -> hue_count[this_hue]++;
		}
	nonzero = 0;
	for (h = 0; h < ms -> table -> nmaterial; h++) {
		if (table -> hue_count[h] == 0) continue;
		nonzero++;
	}
	if (nonzero == 0) nonzero = 1;
	fnshades = (double) MAX_GRAY_LEVELS / (nonzero+1);
	/* fnshades = ms -> n_shades; */
	/* compute color table */
	/* background color: black */
	table -> ucrs[0] = 0;
	table -> ucgs[0] = 0;
	table -> ucbs[0] = 0;
	idx = 1;
	for (h = 0; h < nmaterial; h++) {
		if (idx > 255) break;
		if (table -> hue_count[h] == 0) continue;
		table -> start[h] = idx;
		for (s = 0; s < fnshades; s++, idx++) {
			if (idx > 255) break;
			table -> nshades[h] = s + 1;
			fraction = (s + 0.5) / fnshades;
			if (fraction >= 1.0) fraction = 1.0;
			if (fraction <= 0.0) fraction = 0.0;
			r = fraction * table -> red[h];
			g = fraction * table -> green[h];
			b = fraction * table -> blue[h];
			ucr = (256 - MAX_GRAY_LEVELS) + MAX_GRAY_LEVELS * r;
			ucg = (256 - MAX_GRAY_LEVELS) + MAX_GRAY_LEVELS * g;
			ucb = (256 - MAX_GRAY_LEVELS) + MAX_GRAY_LEVELS * b;
			table -> ucrs[idx] = ucr;
			table -> ucgs[idx] = ucg;
			table -> ucbs[idx] = ucb;
		}
	}
	sprintf (message, "%8ld background pixels", nbackground);
	inform (message);
	for (h = 0; h < ms -> table -> nmaterial; h++) {
		if (ms -> table -> hue_count[h] == 0) continue;
		sprintf (message,"%8ld %4.2f %4.2f %4.2f pixels; table range: %3d-%3d",
			table -> hue_count[h], table -> red[h], table -> green[h], table -> blue[h],
			table -> start[h], table -> start[h] + table -> nshades[h] - 1);
		inform(message);
	}
}

int hsa2value (struct material_table *table, int this_hue, int this_shade, int this_alpha)
{
	int value;
	int s;
	
	if (this_alpha == 0) return (0); /* background pixel */
	if (this_shade == 0) return (0); /* background pixel */
	s = (this_shade / 256.0) * table -> nshades[this_hue];
	if (s <= 1) s = 1;
	value = table -> start[this_hue] + s;
	return (value);
}

int hsa2red (struct material_table *table, int this_hue, int this_shade, int this_alpha)
{
	int ucr;
	/* value = hsa2value (table, this_hue, this_shade, this_alpha); ucr = ms -> ucrs[value]; */
	if (this_alpha == 0) return (0); /* background pixel */
	ucr = this_shade * table -> red[this_hue];
	return (ucr);
}

int hsa2green (struct material_table *table, int this_hue, int this_shade, int this_alpha)
{
	int ucg;
	/* value = hsa2value (table, this_hue, this_shade, this_alpha); ucg = ms -> ucgs[value]; */
	if (this_alpha == 0) return (0); /* background pixel */
	ucg = this_shade * table -> green[this_hue];
	return (ucg);
}

int hsa2blue (struct material_table *table, int this_hue, int this_shade, int this_alpha)
{
	int ucb;
	/* value = hsa2value (table, this_hue, this_shade, this_alpha); ucb = ms -> ucbs[value]; */
	if (this_alpha == 0) return (0); /* background pixel */
	ucb = this_shade * table -> blue[this_hue];
	return (ucb);
}

int code_color (double red, double green, double blue, char name[MAX_NAME])
{
	int ired, igreen, iblue;
	char cn[MAX_NAME];
	ired = 0.5 + 100.0 * red;
	igreen = 0.5 + 100.0 * green;
	iblue = 0.5 + 100.0 * blue;
	sprintf (cn, "r%3dg%3db%3d", ired, igreen, iblue);
	/* leading zeroes */
	if (cn[1] == ' ') cn[1] = '0';
	if (cn[2] == ' ') cn[2] = '0';
	if (cn[5] == ' ') cn[5] = '0';
	if (cn[6] == ' ') cn[6] = '0';
	if (cn[9] == ' ') cn[9] = '0';
	if (cn[10] == ' ') cn[10] = '0';
	strcpy (name, cn);
	return (1);
}

int rgb_to_color_number (struct material_table *table, double red, double green, double blue)
{
	int hue, m, nmaterial;

	hue = -1; /* illegal value */
	nmaterial = table -> nmaterial;
	for (m = 0; m < nmaterial; m++) {
		if (red == table -> red[m] && green == table -> green[m] && blue == table -> blue[m]) {
			hue = m;
			break;
		}
	}
	return (hue);
}

double rgb_to_wheel (double red, double green, double blue)
{
	/*
	red     120
	green   240
	blue      0
	cyan    300
	magenta  60
	yellow  180
	*/
	double wheel;
	if (red == green && green == blue) wheel = 0.0;
	else if (red <= green && red <= blue)
		wheel = 240.0 * green / (green + blue) + 360.0 * blue / (green + blue);
	else if (green <= red && green <= blue)
		wheel = 120.0 * red / (red + blue) + 0.0 * blue / (red + blue);
	else wheel = 240.0 * green / (green + red) + 120.0 * red / (green + red);
	return (wheel);
}

void name_to_rgb (struct material_table *table, char *color_name, double *red, double *green, double *blue)
{
	int hue;
	hue = get_color_number (table, color_name);
	if (error()) return;
	number_to_rgb (table, hue, red, green, blue);
}

void number_to_rgb (struct material_table *table, int hue, double *red, double *green, double *blue)
{
	*red = table -> red[hue];
	*green = table -> green[hue];
	*blue = table -> blue[hue];
	return;
}

double get_opacity_number (char *opacity_name)
{
	double opacity;

	opacity = 0;
	if (isdigit (*opacity_name)) {
		opacity = atof (opacity_name);
	}
	else if (strcmp(opacity_name, "transparent") == 0) {
		opacity = 0.0;
	}
	else if (strcmp(opacity_name, "opaque") == 0) {
		opacity = 1.0;
	}
	else if (strcmp(opacity_name, "translucent") == 0) {
		opacity = 0.5;
	}
	else {	/* default: opaque */
		opacity = 1.0;
	}
	return (opacity);
}

struct object_scheme *define_scheme (int color_type, int uniform_hue, int which_value, double minval, double maxval)
{
	struct object_scheme *cs;

	cs = allocate_object_scheme ();
	if (cs == NULL) {
		set_error1 ("(define_scheme): memory allocation failure");
		return (NULL);
	}
	cs -> color_type = color_type;
	cs -> which_value = which_value;
	cs -> minval = minval;
	cs -> maxval = maxval;
	cs -> uniform_colors[0] = uniform_hue;
	cs -> uniform_colors[1] = uniform_hue;
	cs -> opacity_type = UNIFORM_OPACITY;
	cs -> uniform_opacities[0] = 1.0;
	cs -> uniform_opacities[1] = 1.0;
	return (cs);
}

struct object_scheme *allocate_object_scheme ()
{
	struct object_scheme *sch;

	sch = (struct object_scheme *) allocate_object (OBJECT_SCHEME);
	if (sch == NULL) {
		set_error1 ("allocate_object_scheme: mem alloc fails");
		return(NULL);
	}
	return (sch);
}

int determine_hue (struct material_table *mt, struct object_scheme *cs, int atm, int comp, int shape, int inner, int input_hue, double value)
{
	int i, hue, ind, comp2, shape2, shape3, nmaterial;
	char coloring[MAX_NAME];
	char message[MAXLINE];
	struct color_ramp *ramp;

	if (mt == NULL) {
		sprintf (message, "determine_hue: null material table");
		set_error1 (message);
		return(0);
	}
	nmaterial = mt -> nmaterial;
	if (cs == NULL) {
		sprintf (message, "determine_hue: null color scheme");
		set_error1 (message);
		return(0);
	}
	if (cs -> color_type == ATOM_COLORING)
		strcpy (coloring, "atom coloring");
	else if (cs -> color_type == COMPONENT_COLORING)
		strcpy (coloring, "component coloring");
	else if (cs -> color_type == UNIFORM_COLORING)
		strcpy (coloring, "uniform coloring");
	else if (cs -> color_type == INPUT_COLORING)
		strcpy (coloring, "input coloring");
	else if (cs -> color_type == VERTEX_COLORING)
		strcpy (coloring, "vertex coloring");
	else if (cs -> color_type == SHAPE_COLORING)
		strcpy (coloring, "shape coloring");
	shape3 = shape;
	if (shape3 < 1) shape3 = 1;
	if (shape3 > 3) shape3 = 3;
	shape2 = shape3;
	if (shape2 > 2) shape2 = 2;
	comp2 = comp;
	if (comp2 < 1) comp2 = 1;
	if (comp2 > 2) comp2 = 2;

	if (cs -> color_type == UNIFORM_COLORING) {
		ind = inner;
		hue = cs -> uniform_colors[ind];
	}
	else if (cs -> color_type == INPUT_COLORING) {
		if (inner) hue = 9 - input_hue;
		else hue = input_hue;
	}
	else if (cs -> color_type == ATOM_COLORING) {
		ind = 4 * (atm - 1) + 2 * inner + shape2 - 1;
		hue = *(cs -> atom_colors + ind);
	}
	else if (cs -> color_type == SHAPE_COLORING) {
		ind = 3 * inner + shape3 - 1;
		hue = cs -> shape_colors[ind];
		}
	else if (cs -> color_type == COMPONENT_COLORING) {
		ind = 2 * (comp2 - 1) + inner;
		hue = cs -> component_colors[ind];
	}
	else if (cs -> color_type == VERTEX_COLORING) {
		ramp = cs -> ramp;
		if (ramp == NULL) {
			set_error1 ("determine_hue: null ramp");
			return (0);
		}
		for (i = 0; i < ramp -> n_ramp; i++) {
			if ((i == 0 || value >= ramp -> vertex_ranges[i][0]) &&
				(i == ramp -> n_ramp-1 || value <= ramp -> vertex_ranges[i][1])) {
				hue = ramp -> vertex_colors[i][inner];
				break;
			}
		}
	}
	else {
		sprintf (message, "(determine_hue): invalid color scheme type: %s", coloring);
		set_error1 (message);
		return (0);
	}

	if (hue < 0 || hue >= nmaterial) {
		sprintf (message, "(determine_hue): invalid hue = %8d",
				hue);
		set_error1 (message);
		sprintf (message, "%s atom = %4d, nmaterial = %3d", coloring, atm, nmaterial);
		set_error2 (message);
		return(0);
	}
	return (hue);
}

int get_color_number (struct material_table *table, char *color_name)
{
	int hue, m, nmaterial;

	hue = -1; /* illegal value */
	nmaterial = table -> nmaterial;
	if (isdigit (*color_name)) {
		hue = atoi (color_name);
	}
	else {
		for (m = 0; m < nmaterial; m++) {
			if (strcmp (table -> name[m], color_name) == 0) {
				hue = m;
				break;
			}
		}
	}
	return (hue);
}


struct material_table *new_material_table (struct msscene *ms)
{
	struct material_table *mt;
	long n_ramp;
	long colors[MAX_LEVELS];
	long divisions[MAX_LEVELS];
	struct color_ramp *ramp;
	
	mt = (struct material_table *) allocate_object (MATERIAL_TABLE);
	if (mt == NULL) {
		set_error1 ("new_material_table fails");
		return (NULL);
	}
	/* set up materials */
	add_material (mt, 0.25, 0.25, 0.25, "black");
	if (error ()) return (NULL);
	add_material (mt, 1.0, 1.0, 1.0, "white");
	if (error ()) return (NULL);
	add_material (mt, 1.0, 0.0, 0.0, "red");
	if (error ()) return (NULL);
	add_material (mt, 0.0, 1.0, 0.0, "green");
	if (error ()) return (NULL);
	add_material (mt, 0.0, 0.0, 1.0, "blue");
	if (error ()) return (NULL);
	add_material (mt, 0.0, 1.0, 1.0, "cyan");
	if (error ()) return (NULL);
	add_material (mt, 1.0, 0.0, 1.0, "magenta");
	if (error ()) return (NULL);
	add_material (mt, 1.0, 1.0, 0.0, "yellow");
	if (error ()) return (NULL);
	add_material (mt, 0.50, 0.50, 0.50, "grey");
	if (error ()) return (NULL);
	add_material (mt, 0.75, 0.75, 0.75, "gray");
	if (error ()) return (NULL);
	add_material (mt, 0.4, 0.2, 0.1, "brown");
	if (error ()) return (NULL);
	add_material (mt, 1.0, 0.4, 0.0, "orange");
	if (error ()) return (NULL);
	add_material (mt, 0.5, 0.0, 0.3, "purple");
	if (error ()) return (NULL);
	add_material (mt, 0.0, 0.0, 0.3, "navy");
	if (error ()) return (NULL);
	add_material (mt, 0.4, 0.4, 1.0, "sky");
	if (error ()) return (NULL);
	add_material (mt, 1.0, 0.5, 0.5, "pink");
	if (error ()) return (NULL);
	add_material (mt, 0.8, 1.0, 0.0, "yellow_green");
	if (error ()) return (NULL);
	add_material (mt, 0.8, 0.4, 0.1, "tan");
	if (error ()) return (NULL);

	/* kludge for PSL and color name-number conversion */
	table = mt;
	ms -> table = mt;


	/* set up short ramps */
	n_ramp = 3;
	colors[0] = get_color_number (mt, "yellow");
	colors[1] = get_color_number (mt, "green");
	colors[2] = get_color_number (mt, "blue");
	ramp = define_ramp (ms, n_ramp, colors, "green");
	if (error ()) return (NULL);

	colors[0] = get_color_number (mt, "yellow");
	colors[1] = get_color_number (mt, "red");
	colors[2] = get_color_number (mt, "blue");
	ramp = define_ramp (ms, n_ramp, colors, "red");
	if (error ()) return (NULL);

	colors[0] = get_color_number (mt, "sky");
	colors[1] = get_color_number (mt, "blue");
	colors[2] = get_color_number (mt, "navy");
	ramp = define_ramp (ms, n_ramp, colors, "blue");
	if (error ()) return (NULL);

	/* set up medium ramps */
	/* subdivisions between given colors */
	divisions[0] = 2;
	divisions[1] = 2;
	divisions[2] = 2;

	colors[0] = get_color_number (mt, "yellow");
	colors[1] = get_color_number (mt, "tan");
	colors[2] = get_color_number (mt, "brown");
	ramp = create_ramp (ms, n_ramp, colors, divisions, "tan");
	if (error ()) return (NULL);

	colors[0] = get_color_number (mt, "yellow");
	colors[1] = get_color_number (mt, "orange");
	colors[2] = get_color_number (mt, "red");
	ramp = create_ramp (ms, n_ramp, colors, divisions, "orange");
	if (error ()) return (NULL);

	colors[0] = get_color_number (mt, "white");
	colors[1] = get_color_number (mt, "magenta");
	colors[2] = get_color_number (mt, "purple");
	ramp = create_ramp (ms, n_ramp, colors, divisions, "magenta");
	if (error ()) return (NULL);

	colors[0] = get_color_number (mt, "white");
	colors[1] = get_color_number (mt, "pink");
	colors[2] = get_color_number (mt, "red");
	ramp = create_ramp (ms, n_ramp, colors, divisions, "pink");
	if (error ()) return (NULL);

	return (mt);
}


struct color_ramp *lookup_ramp (struct msscene *ms, char *name)
{
	struct color_ramp *ramp;
	for (ramp = ms -> head_ramp; ramp != NULL; ramp = ramp -> next) {
		if (strcmp (ramp -> name, name) == 0) return (ramp);
	}
	return (NULL);
}


int value2material (struct object_scheme *scheme, double value)
{
	double fraction;
	double range;
	double minval, maxval;
	int r;
	int material;
	struct color_ramp *ramp;

	ramp = scheme -> ramp;
	if (ramp == NULL) {
		set_error1 ("value2material: null ramp");
		return (0);
	}
	minval = scheme -> minval;
	maxval = scheme -> maxval;
	range = maxval - minval;
	if (range <= 0.0) return (-1);
	fraction = (value - minval) / range;
	r = fraction * ramp -> n_ramp;
	if (r < 0) r = 0;
	else if (r > ramp -> n_ramp - 1) r = ramp -> n_ramp - 1;
	material = ramp -> colors[r];
	return (material);
}


struct color_ramp *define_ramp (struct msscene *ms, long n_ramp, long colors[], char *name)
{
	int n;
	struct color_ramp *ramp;
	ramp = (struct color_ramp *) allocate_object (COLOR_RAMP);
	if (ramp == NULL) {
		set_error1 ("define_ramp: no memory");
		return (NULL);
	}
	if (n_ramp < 1 || n_ramp > MAX_LEVELS) {
		set_error1 ("define_ramp: invalid number of levels for color ramp");
		return (0);
	}
	if (strlen (name) > MAX_NAME) {
		set_error1 ("define_ramp: ramp name too long");
		return (0);
	}
	for (n = 0; n < n_ramp; n++)
		ramp -> colors[n] = colors[n];
	strcpy (ramp -> name, name);
	ramp -> n_ramp = n_ramp;
	if (ms -> head_ramp == NULL) ms -> head_ramp = ramp;
	else ms -> tail_ramp -> next = ramp;
	ms -> tail_ramp = ramp;
	
	return (ramp);
}

int add_material (struct material_table *table, double red, double green, double blue, char *name)
{
	int nmaterial, n;
	if (strlen (name) > 15) {
		set_error1 ("add_material: color name too long (>15 chars)");
		return (0);
	}
	nmaterial = table -> nmaterial;
	if (nmaterial >= MAX_MATERIAL) {
		set_error1 ("add_material: material table overflow");
		return (0);
	}
	for (n = 0; n < nmaterial; n++) {
		if (strcmp (table -> name[n], name) == 0) {
			set_error1 ("add_material: duplicate color name");
			return (0);
		}
	}
	table -> red[nmaterial] = red;
	table -> green[nmaterial] = green;
	table -> blue[nmaterial] = blue;
	strcpy (table -> name[nmaterial], name);
	table -> nmaterial++;
	return (table -> nmaterial-1);
}


struct color_ramp *create_ramp (struct msscene *ms, long n_color, long colors[], long divisions[], char *name)
{
	long total, div, n, i;
	long finer[MAX_LEVELS];
	long some[MAX_LEVELS];
	struct color_ramp *ramp;
	struct material_table *mt;
	
	mt = ms -> table;
	if (mt == NULL) {
		set_error1 ("create_ramp: null material table");
		return(NULL);
	}

	if (n_color < 2) {
		set_error1 ("create_ramp: too few colors");
		return (NULL);
	}
	total = 0;
	for (n = 0; n < n_color - 1; n++) {
		div = divisions[n];
		if (div < 0 || div > MAX_LEVELS - 3) {
			set_error1 ("create_ramp: too few/many divisions");
			return (NULL);
		}
		if (total + div > MAX_LEVELS - 3) {
			set_error1 ("create_ramp: too few/many divisions");
			return (NULL);
		}
		finer[total] = colors[n];
		total++;
		if (div == 0) continue;
		create_part (mt, colors[n], colors[n+1], div, some);
		if (error ()) return (NULL);
		for (i = 0; i < div; i++) {
			finer[total] = some[i];
			total++;
		}
	}
	finer[total] = colors[n];
	total++;
	ramp = define_ramp (ms, total, finer, name);
	return (ramp);
}

int create_part (struct material_table *mt, long color0, long color1, long div, long some [])
{
	double reds[MAX_LEVELS];
	double greens[MAX_LEVELS];
	double blues[MAX_LEVELS];
	double red0, green0, blue0, red1, green1, blue1;
	long i, hue;
	char names[MAX_LEVELS][MAX_NAME];
	
	red0 = mt -> red[color0]; green0 = mt -> green[color0]; blue0 = mt -> blue[color0];
	red1 = mt -> red[color1]; green1 = mt -> green[color1]; blue1 = mt -> blue[color1];
	for (i = 0; i < div; i++) {
		reds[i] = red0 + (i+1) * (red1 - red0) / (div+1);
		greens[i] = green0 + (i+1) * (green1 - green0) / (div+1);
		blues[i] = blue0 + (i+1) * (blue1 - blue0) / (div+1);
	}
	/* round down to hundredths */
	for (i = 0; i < div; i++) {
		reds[i] = floor (100.0 * reds[i]) / 100.0;
		greens[i] = floor (100.0 * greens[i]) / 100.0;
		blues[i] = floor (100.0 * blues[i]) / 100.0;
	}
	for (i = 0; i < div; i++) {
		code_color (reds[i], greens[i], blues[i], names[i]);
		/* look for same name */
		hue = get_color_number (mt, names[i]);
		if (error ()) return (0);
		if (hue >= 0) {
			if (mt -> red[hue] != reds[i] || mt -> green[hue] != greens[i]
				|| mt -> blue[hue] != blues[i]) {
				set_error1 ("create_part: interpolated color name already defined differently");
				return (0);
			}
			/* we are okay */
			some[i] = hue;
			continue;
		}
		/* look for same rgb */
		hue = rgb_to_color_number (mt, reds[i], greens[i], blues[i]);
		if (error ()) return (0);
		if (hue >= 0) {
			some[i] = hue;
			continue;
		}
		/* define new color */
		hue = add_material (mt, reds[i], greens[i], blues[i], names[i]);
		if (error ()) return (0);
		some[i] = hue;
	}
	return (1);
}

