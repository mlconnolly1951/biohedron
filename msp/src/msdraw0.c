/* MSDraw 
   Copyright 2006 by Michael L. Connolly
   All rights reserved
   Last revised: February 17, 2006
*/

#include "ms.h"
#include "msenum.h"
#include "mslimit.h"
#include "msdefault.h"
#include "msstruct.h"
#include "msglobal.h"
#include "msfunction.h"

char *intprop = "Copyright 2006 by Michael L. Connolly";
char *southam = "All Rights Reserved";
struct argrec argarray[] = {
	{ 'a', 0, 0, 0, 0.0, ""},
	{ 'b', 0, 0, 0, 0.0, ""},
	{ 'd', 0, 0, 0, 0.0, ""},
	{ 'e', 0, 0, 0, 0.0, ""},
	{ 'f', 0, 0, 0, 0.0, ""},
	{ 'g', 0, 0, 0, 0.0, ""},
	{ 'h', 0, 0, 0, 0.0, ""},
	{ 'i', 0, 1, 0, 0.0, ""},
	{ 'n', 0, 1, 0, 0.0, ""},
	{ 'p', 0, 2, 0, 0.0, ""},
	{ 'r', 0, 2, 0, 0.0, ""},
	{ 's', 0, 0, 0, 0.0, ""},
	{ 't', 0, 0, 0, 0.0, ""},
	{ 'v', 0, 0, 0, 0.0, ""},
	{ 'x', 0, 0, 0, 0.0, ""},
	{ 'y', 0, 0, 0, 0.0, ""},
	{ 'z', 0, 0, 0, 0.0, ""}
};


main (int argc, char *argv[])
{
	struct msscene *ms;
	int nitems, ngotten, nw, k, clipping;
	unsigned long size1, size2;
	int arg, vrml;
	int raster_format, plot_format, format_number, vector_format;
	int screen_size;
	double alignment, border, clip_fraction, stereo_angle, tiltx, tilty, fineness;
	char word1[64], word2[64], word3[64], title[80];
	char shading[64];
	char raster_file[64], raster_right[64], rformat[64], vformat[64];
	char plot_file[64], plot_right[64], pformat[64], vector_file[64];
	char scriptfile[MAXLINE];
	char message[MAXLINE];
	FILE *fpi;				/* command file */
	FILE *fpp, *fpp2;
	FILE *fpr, *fpr2;
	FILE *fpe;
	FILE *fpv;

	fperror = stderr;
	fpinform = stderr;
	fpdebug = stderr;

	sprintf (message, "msdraw   version  %d.%d.%d", (int) VERSION, (int) SUBVERSION, (int) STUBVERSION);
	strcat (message,"   Copyright  2006  by  Michael L. Connolly");
	inform(message);
	if (argc < 2) {
		fprintf (stderr, "usage: msdraw -i command_file -r image [format] -p plot [format] -v vectors\n");
		fprintf (stderr, "[-a stereo_angle] [-b border] [-f fineness] [-g shading]  [-h header]\n");
		fprintf (stderr, "[-n alignment] [-s image_size] [-t title] [-x tiltx] [-y tilty] [-z zclip]\n");
		exit (1);
	}
	/* initialize */
	fpi = fpp = fpp2 = fpr = fpr2 = fpv = (FILE *) NULL;	/* no default */
	fpe = stderr;			/* default */
	debug = 0; screen_size = 0; vrml = 0;
	alignment = 0.0; border = 0.0; clip_fraction = 0.0; stereo_angle = 0.0; clipping = 0;
	tiltx = 0.0; tilty = 0.0; fineness = 0.0;
	/* set up default output format */
	strcpy(pformat,"ps");
	strcpy(rformat,"sgi");
	strcpy(vformat,"sgi");
	strcpy(scriptfile,"");
	strcpy(title, "");
	strcpy(shading, "");
	
	/* streams and flags */
	/* get flags from command line */
	size1 = sizeof(argarray);
	size2 = sizeof(struct argrec);
	nitems = size1/size2;
	ngotten = get_arguments (argc, argv, argarray, nitems);
	if (ngotten <= 0) {
		if (ngotten < -1)
			fprintf (stderr, "msdraw: invalid flag: %c\n", -ngotten);
		else fprintf (stderr, "msdraw: problem reading command-line arguments\n");
		exit (1);
	}
	for (arg = 0; arg < nitems; arg++) {
		if (!argarray[arg].on) continue;
		switch (argarray[arg].flag) {
		case 'a':
			stereo_angle = argarray[arg].doublearg;
			break;
		case 'b':
			border = argarray[arg].doublearg;
			break;
		case 'd':
			debug = argarray[arg].longarg;
			break;
		case 'e':
			if (strcmp(argarray[arg].stringarg,"stderr") == 0) fpe = stderr;
			else if (strcmp(argarray[arg].stringarg,"stdout") == 0) fpe = stdout;
			else {
				fpe = fopen (argarray[arg].stringarg, "w");
				if (fpe == NULL) {
					fprintf (stderr,
						"msdraw: cannot open error file: %s\n",
						argarray[arg].stringarg);
					exit (1);
				}
			}
			break;
		case 'f':
			fineness = argarray[arg].doublearg;
			break;
		case 'g':
			nw = sscanf(argarray[arg].stringarg, "%s %s %s", word1, word2, word3);
			if (nw < 1) continue;
			strcpy (shading, word1);
			break;
		case 'h':
			nw = sscanf(argarray[arg].stringarg, "%s %s %s", word1, word2, word3);
			if (nw < 1) continue;
			if (strcmp(word1, "vrml") == 0) vrml = 1;
			if (strcmp(word1, "VRML") == 0) vrml = 1;
			if (strcmp(word1, "inventor") == 0) vrml = 0;
			if (strcmp(word1, "INVENTOR") == 0) vrml = 0;
			if (strcmp(word1, "Inventor") == 0) vrml = 0;
			break;
		case 'i':
			strcpy(scriptfile,argarray[arg].stringarg);
			fpi = fopen (scriptfile, "r");
			if (fpi == NULL) {
				fprintf (stderr,"msdraw: cannot open script file\n");
				exit (1);
			}
			break;
		case 'n':
			alignment = argarray[arg].doublearg;
			break;
		case 'p':
			nw = sscanf(argarray[arg].stringarg, "%s %s %s", word1, word2, word3);
			if (nw < 1) continue;
			if (nw == 2) {
				format_number = get_format_number(word2);
				if (format_number > 0) strcpy (pformat, word2);
			}
			else format_number = 0;
			if (nw == 1 || nw == 2 && format_number > 0) {
				strcpy (plot_file, word1);
				fpp = fopen (plot_file, "w");
				if (fpp == NULL) {
					fprintf (stderr,
						"msdraw: cannot open plot output file: %s\n", plot_file);
					exit (1);
				}
				break;
			}
			/* later stuff for both left and right in one run */
			if (nw >= 3) {
				format_number = get_format_number(word3);
				if (format_number > 0) strcpy (pformat, word3);
			}
			strcpy (plot_file, word1);
			fpp = fopen (plot_file, "w");
			if (fpp == NULL) {
				fprintf (stderr,
					"msdraw: cannot open left plot output file: %s\n", plot_file);
				exit (1);
			}
			strcpy (plot_right, word2);
			fpp2 = fopen (plot_right, "w");
			if (fpp2 == NULL) {
				fprintf (stderr,
					"msdraw: cannot open right plot output file: %s\n", plot_right);
				exit (1);
			}
			break;
		case 'r':
			nw = sscanf(argarray[arg].stringarg, "%s %s %s", word1, word2, word3);
			if (nw < 1) continue;
			if (nw == 2) {
				format_number = get_format_number(word2);
				if (format_number > 0) strcpy (rformat, word2);
			}
			else format_number = 0;
			if (nw == 1 || nw == 2 && format_number > 0) {
				strcpy (raster_file, word1);
				fpr = fopen (raster_file, "w");
				if (fpr == NULL) {
					fprintf (stderr,
						"msdraw: cannot open raster output file: %s\n", raster_file);
					exit (1);
				}
				break;
			}
			/* later stuff for both left and right in one run */
			if (nw >= 3) {
				format_number = get_format_number(word3);
				if (format_number > 0) strcpy (rformat, word3);
			}
			strcpy (raster_file, word1);
			fpr = fopen (raster_file, "w");
			if (fpr == NULL) {
				fprintf (stderr,
					"msdraw: cannot open left raster output file: %s\n", raster_file);
				exit (1);
			}
			strcpy (raster_right, word2);
			fpr2 = fopen (raster_right, "w");
			if (fpr2 == NULL) {
				fprintf (stderr,
					"msdraw: cannot open right raster output file: %s\n", raster_right);
				exit (1);
			}
			break;
		case 's':
			screen_size = argarray[arg].longarg;
			break;
		case 't':
			strcpy (title, argarray[arg].stringarg);
			break;
		case 'v':
			nw = sscanf(argarray[arg].stringarg, "%s %s %s", word1, word2, word3);
			if (nw < 1) continue;
			if (nw == 2) {
				format_number = get_format_number(word2);
				if (format_number > 0) strcpy (vformat, word2);
			}
			else format_number = 0;
			if (nw == 1 || nw == 2 && format_number > 0) {
				strcpy (vector_file, word1);
				fpv = fopen (vector_file, "w");
				if (fpv == NULL) {
					fprintf (stderr,
						"msdraw: cannot open vector output file: %s\n", vector_file);
					exit (1);
				}
			}
			break;
		case 'x':
			tiltx = argarray[arg].doublearg;
			break;
		case 'y':
			tilty = argarray[arg].doublearg;
			break;
		case 'z':
			clip_fraction = argarray[arg].doublearg;
			clipping = 1;
			break;
		default:
			fprintf (stderr,
				"msdraw: invalid flag: %s\n", argarray[arg].stringarg);
			exit (1);
		}
	}


	init_routing (fpe, fpe, fpe);


	/* some error checking */

	if (fpi == NULL) {
		fprintf (stderr, "msdraw: no script file specified\n");
		exit (1);
	}
	if (strcmp(pformat,"hpgl") == 0)
		plot_format = HPGL;
	else if (strcmp(pformat,"ps") == 0)
		plot_format = POSTSCRIPT;
	else {
		fprintf (stderr, "msdraw: invalid plot format, try hpgl or ps\n");
		exit (1);
	}
	if (strcmp(rformat,"sun") == 0)
		raster_format = SUNRASTER;
	else if (strcmp(rformat,"sgi") == 0)
		raster_format = SGIIMAGE;
	else if (strcmp(rformat,"avsimage") == 0)
		raster_format = AVSIMAGE;
	else if (strcmp(rformat,"avsfield") == 0)
		raster_format = AVSFIELD;
	else if (strcmp(rformat,"bmp") == 0)
		raster_format = BMP;
	else {
		fprintf (stderr, "msdraw: invalid raster format, try bmp, sun or sgi\n");
		exit (1);
	}
	if (strcmp(vformat,"inventor") == 0)
		vector_format = INVENTOR;
	else if (strcmp(vformat,"iv") == 0)
		vector_format = INVENTOR;
	else if (strcmp(vformat,"sgi") == 0)
		vector_format = INVENTOR;
	else {
		fprintf (stderr, "msdraw: invalid vector format, try sgi, iv or inventor\n");
		exit (1);
	}
	if (fpr == NULL) raster_format = 0;
	if (fpp == NULL) plot_format = 0;
	if (fpv == NULL) vector_format = 0;

	/* initialize memory, scene and material table */
	draw_mem ();
	if (error()) return(0);
	ms = new_msscene ();
	if (ms == NULL || error ()) {
		set_error2 ("msdraw: new_msscene fails");
		print_error ();
		exit (1);
	}

	/* initialization of Protein Shape Library */
	psl_init ();
	if (error()) {
		print_error ();
		exit (1);
	}
	init_mol();
	if (error()) {
		print_error ();
		exit (1);
	}
	init_region();
	if (error()) {
		print_error ();
		exit (1);
	}
	init_pdb ();
	if (error()) {
		print_error ();
		exit (1);
	}

	init_counters (ms);
	if (error()) return(0);
	ms -> table = new_material_table (ms);
	if (error ()) {
		set_error2 ("msdraw: new_material_table fails");
		print_error ();
		exit (1);
	}

	/* initialize parameters */
	init_parameters (ms);
	if (error()) {
		fprintf (stderr, "msdraw: parameter initialization failed\n");
		exit (1);
	}

	/* set parameters based upon command line flags */
	if (debug) ms -> debug = debug;
	ms -> raster_format = raster_format;
	ms -> plot_format = plot_format;
	ms -> vector_format = vector_format;
	ms -> fp_raster = fpr;
	ms -> fp_raster2 = fpr2;
	ms -> fp_plot = fpp;
	ms -> fp_plot2 = fpp2;
	ms -> fp_vector = fpv;
	if (strlen (title) > (unsigned long) 0L && strlen (title) < (unsigned long) 80L)
		strcpy (ms -> title, title);
	if (vrml != 0) ms -> vrml = vrml;
	if (strlen(shading) > (unsigned long) 0L) {
		if (strcmp(shading, "flat") == 0) ms -> interpolate = 0;
		else if (strcmp(shading, "facet") == 0) ms -> interpolate = 0;
		else if (strcmp(shading, "smooth") == 0) ms -> interpolate = 1;
		else if (strcmp(shading, "gouraud") == 0) ms -> interpolate = 1;
		else if (strcmp(shading, "phong") == 0) ms -> interpolate = 1;
	}
	if (screen_size > 0) {
		/* error checking */
		for (k = 0; k < 2; k++) {
			if (screen_size < MIN_SCREEN || screen_size > MAX_SCREEN) {
				set_error1("invalid screen size");
				return (0);
			}
		}
		ms -> screen_size = screen_size; 
		ms -> viewport[0][0] = 0;
		ms -> viewport[1][0] = ms -> screen_size;
		ms -> viewport[0][1] = 0;
		ms -> viewport[1][1] = ms -> screen_size;
	}
	if (border > 0.0) {
		ms -> border = -border;
	}
	if (alignment > 0.0) {
		ms -> alignment = alignment;
	}
	if (fineness > 0.0) {
		ms -> fineness = fineness;
	}
	if (clipping) {
		if (clip_fraction > 1.0) {
			set_error1("z clip > 1");
			print_error ();
			return (0);
		}
		if (clip_fraction < -1.0) {
			set_error1("z clip < -1");
			print_error ();
			return (0);
		}
		ms -> clipping = 1;
		ms -> clip_fraction = clip_fraction;
	}
	if (stereo_angle != 0.0) {
		ms -> stereo_angle = stereo_angle;
		if (ms -> stereo_angle > 0.0) ms -> stereo = 1;
		else if (ms -> stereo_angle < 0.0) ms -> stereo = 1;
	}
	if (tiltx != 0.0) {
		ms -> tiltx = tiltx;
		/* error checking */
		if (ms -> tiltx < MIN_TILTX || ms -> tiltx > MAX_TILTX) {
			sprintf (message, "minimum tiltx = %d", (int) MIN_TILTX);
			set_error1(message);
			sprintf (message, "maximum tiltx = %d", (int) MAX_TILTX);
			set_error2(message);
			print_error ();
			return (0);
		}
	}
	if (tilty != 0.0) {
		ms -> tilty = tilty;
		/* error checking */
		if (ms -> tilty < MIN_TILTY || ms -> tilty > MAX_TILTY) {
			sprintf (message, "minimum tilty = %d", (int) MIN_TILTY);
			set_error1(message);
			sprintf (message, "maximum tilty = %d", (int) MAX_TILTY);
			set_error2(message);
			print_error ();
			return (0);
		}
	}


	do_draw (ms, fpi);
	free_scene (ms);
	if (error()) {
		print_error();
		exit (1);
	}
	free_all_pdb ();
	if (error()) {
		print_error();
		exit (1);
	}
	free_all_psl ();
	if (error()) {
		print_error();
		exit (1);
	}
	print_counts();
	if (error()) {
		print_error();
		exit (1);
	}
	return (0);
}



/*
 * MSDraw
 * Copyright 2006 by Michael L. Connolly
 * All Rights Reserved
 */
