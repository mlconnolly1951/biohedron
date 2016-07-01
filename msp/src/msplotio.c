#include "ms.h"
#include "msenum.h"
#include "mslimit.h"
#include "msdefault.h"
#include "msstruct.h"
#include "msglobal.h"
#include "msfunction.h"

/* Molecular Surface Package Copyright 1993 by Michael L. Connolly */
/* February 9, 2000 */


/* Plotting */

void write_plot (struct msscene *ms)
{
	long total_ls, n_ls;
	char file_type[40], typename[25];
	char message[MAXLINE];
	struct surface *head_surface, *b;

	start_plot (ms);

	total_ls = 0;

	for (ms -> current_molecule = ms -> head_molecule; ms -> current_molecule != NULL;
		ms -> current_molecule = ms -> current_molecule -> next) {
		/* copy to global variables */
		head_surface = ms -> current_molecule -> head_surface;
		sprintf (message,"%8ld objects to plot for molecule %s",
			ms -> current_molecule -> n_surface, ms -> current_molecule -> name);
		inform(message);

		for (b = head_surface; b != NULL; b = b -> next) {
			if (b -> type == PHN_SURFACE &&
				ms -> current_molecule -> blank) {
				informd ("polyhedron blanked");
				continue;
			}
			if (b -> head_linseg == (struct linseg *) NULL) continue;
			n_ls = plot_bunch (ms, b);
			if (error()) return;
			total_ls += n_ls;
			get_bunch_name (b -> type, typename);
			sprintf (message,"%8ld line segments plotted for %s object",
				n_ls, typename);
			inform(message);
		}
	}
	terminate_plot (ms);

	fclose (ms -> fp_plot);
	
	if (ms -> plot_format == HPGL)
		strcpy (file_type, "hpgl");
	else if (ms -> plot_format == POSTSCRIPT)
		strcpy (file_type, "PostScript");
	else strcpy (file_type, "unknown-format");

	sprintf (message,"%8ld line segments plotted in %s format",
		total_ls, file_type);
	inform(message);
}

/* OUTPUT */


void start_plot (struct msscene *ms)
{
	double factor, hfactor, width;
	
	factor = (11.0/8.5);
	hfactor = (1 + factor) / 2;
	width = ms -> window[1][0] - ms -> window[0][0];
	ms -> psfactor = 72 * 8.5 / width;

	if (ms -> plot_format == HPGL) {
		fprintf (ms -> fp_plot, "IN;\n");
		fprintf (ms -> fp_plot, "SP1;\n");
		fprintf (ms -> fp_plot, "NP8;\n");
		fprintf (ms -> fp_plot, "SC%10.5f,%10.5f,%10.5f,%10.5f;\n",
			ms -> window[1][0] - hfactor * width,
			ms -> window[0][0] + hfactor * width,
			ms -> window[0][1],
			ms -> window[1][1]);
		fprintf (ms -> fp_plot, "PU;\n");
		fprintf (ms -> fp_plot, "PA0.0,0.0\n");
	}
	else if (ms -> plot_format == POSTSCRIPT) {
		fprintf (ms -> fp_plot, "%%!PS-Adobe-2.0\n");
		fprintf (ms -> fp_plot, "%%%%Title: %s\n", ms -> title);
		fprintf (ms -> fp_plot, "%%%%Creator:  MSDraw %d.%d\n",
			(int) VERSION, (int) SUBVERSION);
		fprintf (ms -> fp_plot, "%%%%BoundingBox: 0 0 %4d %4d\n",
			(int) (72 * 8.5), (int) (72 * 11.0));
		fprintf (ms -> fp_plot, "%4.0f %4.0f translate\n",
			- ms -> psfactor * ms -> window[0][0], - ms -> psfactor * ms -> window[0][1]);
	}
}

/* plot bunch */

long plot_bunch (struct msscene *ms, struct surface *b)
{
	int hue, h, inner, comp, atm, shape, input_hue;
	long n_ls;
	double red, green, blue;
	double linewidth, value;
	struct linseg *ls;
	struct phnedg *edg;
	struct object_scheme *cs;
	struct material_table *mt;

	n_ls = 0;
	cs = b -> scheme;
	mt = ms -> table;
	if (mt == NULL) {
		set_error1 ("plot_bunch: null material table");
		return(0L);
	}
	for (ls = b -> head_linseg; ls != NULL; ls = ls -> next) {
		if (ls -> hidden) continue;
		n_ls++;
		shape = 1;
		edg = ls -> edg;
		inner = ls -> inner;
		atm = ls -> atm;
		comp = ls -> comp;
		value = ls -> value;
		input_hue = ls -> hue;
		hue = determine_hue (mt, cs, atm, comp, shape, inner, input_hue, value);
		number_to_rgb (ms -> table, hue, &red, &green, &blue);
		set_rgb (ms, red, green, blue);
		/* linewidth */
		h = 0;
		if (edg != NULL && edg -> comp > 1) h += 2;
		if (ls -> inner) h += 1;
		linewidth = b -> linewidth[h];
		set_linewidth (ms, linewidth);
		output_linseg (ms, ls);
	}
	return (n_ls);
}

/* the pen choice is based upon HP's documentation
   on Hewlett-Packard Graphics Language (HP-GL) as described in:
   The HP-GL/2 Reference Guide,
   Addison-Wesley, Menlo Park, CA (1990) */

void print_rgb (struct msscene *ms, double red, double green, double blue)
{
	int pen_number;

	if (ms -> plot_format == HPGL) {
		if (red > 0.5 && green > 0.5 && blue > 0.5)
			pen_number = 0; /* white */
		else if (red > 0.5 && green <= 0.5 && blue <= 0.5)
			pen_number = 2; /* red */
		else if (red <= 0.5 && green > 0.5 && blue <= 0.5)
			pen_number = 3; /* green */
		else if (red <= 0.5 && green <= 0.5 && blue > 0.5)
			pen_number = 5; /* blue */
		else if (red <= 0.5 && green > 0.5 && blue > 0.5)
			pen_number = 7; /* cyan */
		else if (red > 0.5 && green <= 0.5 && blue > 0.5)
			pen_number = 6; /* magenta */
		else if (red > 0.5 && green > 0.5 && blue <= 0.5)
			pen_number = 4; /* yellow */
		else 
			pen_number = 1; /* black */
		fprintf (ms -> fp_plot, "SP%d;\n", pen_number);
	}
	else if (ms -> plot_format == POSTSCRIPT) {
		fprintf (ms -> fp_plot, "%4.2f %4.2f %4.2f setrgbcolor\n",
			red, green, blue);
	}
}

void output_linewidth (struct msscene *ms, double factor)
{
	if (ms -> plot_format == HPGL)
		fprintf (ms -> fp_plot, "PW%5.2f;\n",0.35 * factor);
	else if (ms -> plot_format == POSTSCRIPT)
		fprintf (ms -> fp_plot,"%5.2f setlinewidth\n",factor);
}

void output_linseg (struct msscene *ms, struct linseg *ls)
{
	if (ms -> plot_format == HPGL)
		fprintf (ms -> fp_plot, "PU;PA%10.5f,%10.5f;PD;PA%10.5f,%10.5f;\n",
			ls -> ends[0][0], ls -> ends[0][1],
			ls -> ends[1][0], ls -> ends[1][1]);
	else if (ms -> plot_format == POSTSCRIPT)
		fprintf (ms -> fp_plot,
			"%4.0f %4.0f moveto %4.0f %4.0f lineto stroke\n",
			ms -> psfactor * ls -> ends[0][0], ms -> psfactor * ls -> ends[0][1],
			ms -> psfactor * ls -> ends[1][0], ms -> psfactor * ls -> ends[1][1]);
}

void terminate_plot (struct msscene *ms)
{
	if (ms -> plot_format == HPGL)
		/* set pen 0, page */
		fprintf (ms -> fp_plot, "SP0;\nPG;\n");
	else if (ms -> plot_format == POSTSCRIPT)
		fprintf (ms -> fp_plot, "showpage\n");
}


void set_rgb (struct msscene *ms, double red, double green, double blue)
{
	if (red == ms -> red_amount && green == ms -> green_amount && blue == ms -> blue_amount)
		return;
	ms -> red_amount = red;
	ms -> green_amount = green;
	ms -> blue_amount = blue;
	print_rgb (ms, ms -> red_amount, ms -> green_amount, ms -> blue_amount);
}

void set_linewidth (struct msscene *ms, double linewidth)
{
	if (linewidth == ms -> current_linewidth)
		return;
	ms -> current_linewidth = linewidth;
	output_linewidth (ms, ms -> current_linewidth);
}

void setup_linewidth (struct surface *obj, double outer_width, double inner_width, double cavity_width, double bond_width)
{
	switch (obj -> type) {
	case BAS_SURFACE:
		obj -> linewidth[0] = bond_width;
		obj -> linewidth[1] = bond_width;
		obj -> linewidth[2] = bond_width;
		obj -> linewidth[3] = bond_width;
		break;
	default:
		obj -> linewidth[0] = outer_width;
		obj -> linewidth[1] = inner_width;
		obj -> linewidth[2] = cavity_width;
		obj -> linewidth[3] = cavity_width;
		break;
	}
}

