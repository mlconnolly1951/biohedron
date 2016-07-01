/*
   MSTran

   MS Translate to other formats

   Copyright 2006 by Michael L. Connolly
   All Rights Reserved

   February 16, 2006
*/

#include "ms.h"
#include "msenum.h"
#include "mslimit.h"
#include "msdefault.h"
#include "msstruct.h"
#include "msglobal.h"
#include "msfunction.h"

struct argrec argarray[] = {
	{ 't', 0, 0, 0, 0.0, ""},
	{ 'd', 0, 0, 0, 0.0, ""},
	{ 'g', 0, 0, 0, 0.0, ""},
	{ 'h', 0, 0, 0, 0.0, ""},
	{ 'l', 0, 0, 0, 0.0, ""},
	{ 'u', 0, 0, 0, 0.0, ""},
	{ 'c', 0, 0, 0, 0.0, ""},
	{ 'm', 0, 0, 0, 0.0, ""},
	{ 'f', 0, 0, 0, 0.0, ""},
	{ 'e', 0, 0, 0, 0.0, ""},
	{ 'o', 0, 0, 0, 0.0, ""},
};

main (int argc, char *argv[])
{
	int arg, result, nitems, ngotten;
	int hue = 0;
	int nw;
	int use_function;
	int which_value;
	int fn;
	double mingiven, maxgiven, ctrlev;
	char fname[64], format[64];
	char message[MAX_STRING];
	char color_name[32];
	unsigned long size1, size2;
	struct msscene *ms;
	struct color_ramp *ramp;
	FILE *fpt;
	FILE *fpd;
	FILE *fpg;
	FILE *fpm;
	FILE *fpe;
	FILE *fpo;

	fperror = stderr;
	fpinform = stderr;
	fpdebug = stderr;
	
	sprintf (message,
		"mstran   version %d.%d.%d    Copyright 2006 by Michael L. Connolly",
		(int) VERSION, (int) SUBVERSION, (int) STUBVERSION);
	inform(message);

	/* initialize */
	strcpy(format,"inventor");	/* default format: SGI Inventor */
	fn = get_format_number (format);
	use_function = ' ';
	which_value = -1;
	fpt = fpd = fpm = fpg = fpo = (FILE *) NULL;	/* no default */
	fpe = stderr;			/* default */
	mingiven = 0.0;
	maxgiven = 1.0;
	ctrlev = 0.5;
	strcpy (color_name, "gray");
	
	/* streams and flags */
	
	if (argc <= 2) {
		fprintf (stderr,
			"usage: mstran -t polyhedron -f u|v|w -l lower -u upper -h color_ramp -o inventor_file\n");
		fprintf (stderr, "usage: mstran -d density -c contour_level -h hue -o inventor_file\n");
		exit(1);
	}
	/* get flags from command line */
	size1 = sizeof(argarray);
	size2 = sizeof(struct argrec);
	nitems = size1/size2;
	ngotten = get_arguments (argc, argv, argarray, nitems);
	if (ngotten < 0) {
		if (ngotten < -1)
			fprintf (stderr, "invalid flag: %c\n", -ngotten);
		else fprintf (stderr, "problem reading command-line arguments\n");
		exit (1);
	}
	for (arg = 0; arg < nitems; arg++) {
		if (!argarray[arg].on) continue;
		switch (argarray[arg].flag) {
		case 't':
			if (strcmp(argarray[arg].stringarg,"stdin") == 0) fpt = stdin;
			else {
				fpt = fopen (argarray[arg].stringarg, "r");
				if (fpt == NULL) {
					fprintf (stderr, "mstran: cannot open polyhedron file: %s\n",
							argarray[arg].stringarg);
					exit (1);
				}
			}
			break;
		case 'd':
			if (strcmp(argarray[arg].stringarg,"stdin") == 0) fpd = stdin;
			else {
				fpd = fopen (argarray[arg].stringarg, "rb");
				if (fpd == NULL) {
					fprintf (stderr, "mstran: cannot open density file: %s\n",
							argarray[arg].stringarg);
					exit (1);
				}
			}
			break;
		case 'g':
			if (strcmp(argarray[arg].stringarg,"stdin") == 0) fpg = stdin;
			else {
				fpg = fopen (argarray[arg].stringarg, "r");
				if (fpg == NULL) {
					fprintf (stderr, "mstran: cannot open oringe file: %s\n",
							argarray[arg].stringarg);
					exit (1);
				}
			}
			break;
		case 'f':
			if (argarray[arg].stringarg != NULL)
				use_function = argarray[arg].stringarg[0];
			break;
		case 'h':
			if (argarray[arg].stringarg != NULL)
				strcpy (color_name, argarray[arg].stringarg);
			break;
		case 'l':
			if (argarray[arg].doublearg != 0.0)
				mingiven = argarray[arg].doublearg;
			break;
		case 'u':
			if (argarray[arg].doublearg != 0.0)
				maxgiven = argarray[arg].doublearg;
			break;
		case 'c':
			if (argarray[arg].doublearg != 0.0)
				ctrlev = argarray[arg].doublearg;
			break;
		case 'm':
			if (strcmp(argarray[arg].stringarg,"stdin") == 0) fpm = stdin;
			else {
				fpm = fopen (argarray[arg].stringarg, "r");
				if (fpm == NULL) {
					fprintf (stderr, "mstran: cannot open atom file: %s\n",
							argarray[arg].stringarg);
					exit (1);
				}
			}
			break;
		case 'o':
			nw = sscanf(argarray[arg].stringarg, "%s %s", fname, format);
			if (nw < 1) continue;
			if (nw >= 2) fn = get_format_number (format);
			if (fn <= 0) {
				fprintf (stderr, "mstran: invalid output format: %s\n", format);
				exit (1);
			}
			if (strcmp(fname,"stdout") == 0) fpo = stdout;
			else {
				fpo = fopen (fname, "w");
				if (fpo == NULL) {
					fprintf (stderr,
						"mstran: cannot open output file: %s\n", fname);
					exit (1);
				}
			}
			break;
		case 'e':
			if (strcmp(argarray[arg].stringarg,"stderr") == 0) fpe = stderr;
			else if (strcmp(argarray[arg].stringarg,"stdout") == 0) fpe = stdout;
			else {
				fpe = fopen (argarray[arg].stringarg, "w");
				if (fpe == NULL) {
					fprintf (stderr, "mstran: cannot open error file: %s\n",
						argarray[arg].stringarg);
					exit (1);
				}
			}
			break;
		}
	}
	if ((fpt == NULL && fpd == NULL && fpg == NULL && fpm == NULL) || fpo == NULL) {
		fprintf (stderr, "missing polyhedron, density, oringe, atom or inventor file name\n");
		exit(1);
	}
	if (use_function == 'u') which_value = 0;
	else if (use_function == 'v') which_value = 1;
	else if (use_function == 'w') which_value = 2;
	if (ctrlev <= 0.0 || ctrlev >= 1.0) ctrlev = 0.5;


	/* initialize scene, memory and material table */
	init_routing (fpe, fpe, fpe);

	tran_mem ();
	if (error()) {
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
	ms = new_msscene ();
	if (ms == NULL || error ()) {
		set_error2 ("mstran: new_msscene fails");
		print_error ();
		exit (1);
	}
	if (debug) ms -> debug = debug;
	ms -> table = new_material_table (ms);
	if (error ()) {
		set_error2 ("mstran: new_material_table fails");
		print_error ();
		exit (1);
	}
	/* global variable kludge */
	table = ms -> table;
	/* hue will be -1 if not valid color name */
	hue = get_color_number (ms -> table, color_name);
	/* ramp will be NULL if not valid ramp name */
	ramp = lookup_ramp (ms, color_name);

	result = mstran (ms, fpt,fpd,fpm,fpg,fpo, fn, which_value,mingiven,maxgiven, ramp, ctrlev, hue);
	fclose(fpo);
	if (!result || error()) {
		print_error();
		exit(1);
	}
	return(0);
}

