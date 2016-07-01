/* 
   MSForm 
   Copyright 2006 by Michael L. Connolly
   All rights reserved
   Last revised: February 16, 2006
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
	{ 'c', 0, 0, 0, 0.0, ""},
	{ 'd', 0, 0, 0, 0.0, ""},
	{ 'e', 0, 0, 0, 0.0, ""},
	{ 'f', 0, 0, 0, 0.0, ""},
	{ 'g', 0, 0, 0, 0.0, ""},
	{ 'h', 0, 0, 0, 0.0, ""},
	{ 'k', 0, 0, 0, 0.0, ""},
	{ 'l', 0, 0, 0, 0.0, ""},
	{ 'm', 0, 0, 0, 0.0, ""},
	{ 'n', 0, 0, 0, 0.0, ""},
	{ 'o', 0, 0, 0, 0.0, ""},
	{ 'p', 0, 0, 0, 0.0, ""},
	{ 'r', 0, 0, 0, 0.0, ""},
	{ 't', 0, 0, 0, 0.0, ""},
	{ 'v', 0, 0, 0, 0.0, ""},
	{ 'w', 0, 0, 0, 0.0, ""},
	{ 'x', 0, 0, 0, 0.0, ""},
	{ 'y', 0, 0, 0, 0.0, ""}
};


main (int argc, char *argv[])
{
	int arg, result, nitems, ngotten, nw;
	int smooth, minus;
	long nfrequency, nsample, nbin1, nsector, nbin;
	unsigned long size1, size2;
	double cube_width, ctrlev, sphere_radius, x_radius;
	char fname[64], format[64];
	char message[128];
	char partition[16];
	FILE *fpa;
	FILE *fpd;
	FILE *fpe;
	FILE *fpf;
	FILE *fpg;
	FILE *fph;
	FILE *fpn;
	FILE *fpo;
	FILE *fpt;
	FILE *fpy;
	struct msscene *ms;

	fperror = stderr;
	fpinform = stderr;
	fpdebug = stderr;

	sprintf (message,
		"msform   version %d.%d.%d    Copyright 2006 by Michael L. Connolly",
		(int) VERSION, (int) SUBVERSION, (int) STUBVERSION);
	inform(message);

	/* initialize */
	nfrequency = DEFAULT_FREQUENCY;
	nbin1 = DEFAULT_BIN;
	nbin = nbin1 * nbin1;
	strcpy (partition, "maze");
	cube_width = 0.0;
	ctrlev = 0.0;
	sphere_radius = 0.0;
	x_radius = 0.0;
	smooth = 0;
	minus = 0;
	nsample = 0;
	nsector = 0;
	strcpy(format, "vet");
	fpa = fpf = (FILE *) NULL;
	fpt = fpn = fph = (FILE *) NULL;
	fpo = (FILE *) NULL;
	fpe = stderr;
	fpt = fpd = (FILE *) NULL;
	fpy = (FILE *) NULL;
	fpg = (FILE *) NULL;
	
	/* streams and flags */
	
	if (argc <= 2) {
		fprintf (stderr, "usage: msform -t polyhedron -r radius -o polyhedron\n");
		fprintf (stderr, "usage: msform -t polyhedron -r radius -x expansion -o polyhedron\n");
		fprintf (stderr, "usage: msform -t polyhedron -n polyhedron -r radius -o polyhedron\n");
		fprintf (stderr, "usage: msform -t polyhedron -n polyhedron -r radius -x extension -o polyhedron\n");
		fprintf (stderr, "usage: msform -t polyhedron -n polyhedron -h polyhedron -r radius -o polyhedron\n");
		fprintf (stderr, "usage: msform -t polyhedron -o density -w width\n");
		fprintf (stderr, "usage: msform -d density -r radius -c level [-m] -o polyhedron \n");
		fprintf (stderr, "usage: msform -d density -y density -r radius -c level [-m] -o polyhedron\n");
		fprintf (stderr, "usage: msform -t polyhedron -d density -r radius -o polyhedron\n");
		fprintf (stderr, "usage: msform -t polyhedron -d density -y density -r radius -o polyhedron \n");
		fprintf (stderr, "usage: msform -t polyhedron -d density -g oringe -o polyhedron\n");
		fprintf (stderr, "usage: msform -t polyhedron -d density -r radius -k nsector -o polyhedron\n");
		fprintf (stderr, "usage: msform -g standard_oringe(s) -o fourier_oringe(s)\n");
		fprintf (stderr, "usage: msform -f fourier_oringe(s) -k nsector -o standard_oringe(s)\n");
		fprintf (stderr, "usage: msform -g standard_oringe(s) [-p partition] [-b nbin] -l nsample -o abstract_oringe(s)\n");
		return (0);
	}

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
			nw = sscanf(argarray[arg].stringarg, "%s %s", fname, format);
			if (nw < 1) continue;
			fpa = fopen (fname, "r");
			if (fpa == NULL) {
				fprintf (stderr, "msform: cannot open abstract input file: %s\n",
						fname);
				return (0);
			}
			break;
		case 'b':
			if (argarray[arg].longarg != 0) {
				nbin = argarray[arg].longarg;
				if (nbin > 256) {
					fprintf (stderr, "msform: nbin (%d) > 256\n", nbin);
					return (0);
				}
				nbin1 = sqrt (nbin);
				nbin = nbin1 * nbin1;
			}
			break;
		case 'c':
			ctrlev = argarray[arg].doublearg;
			break;
		case 'd':
			nw = sscanf(argarray[arg].stringarg, "%s %s", fname, format);
			if (nw < 1) continue;
			fpd = fopen (fname, "rb");
			if (fpd == NULL) {
				fprintf (stderr, "msform: cannot open density file: %s\n",
						fname);
				return (0);
			}
			break;
		case 'e':
			if (strcmp(argarray[arg].stringarg,"stderr") == 0) fpe = stderr;
			else if (strcmp(argarray[arg].stringarg,"stdout") == 0) fpe = stdout;
			else {
				fpe = fopen (argarray[arg].stringarg, "w");
				if (fpe == NULL) {
					fprintf (stderr, "msform: cannot open error file: %s",
						argarray[arg].stringarg);
					exit (1);
				}
			}
			break;
		case 'f':
			nw = sscanf(argarray[arg].stringarg, "%s %s", fname, format);
			if (nw < 1) continue;
			fpf = fopen (fname, "r");
			if (fpf == NULL) {
				fprintf (stderr, "msform: cannot open fourier input file: %s\n",
						fname);
				return (0);
			}
			break;
		case 'g':
			nw = sscanf(argarray[arg].stringarg, "%s %s", fname, format);
			if (nw < 1) continue;
			fpg = fopen (fname, "r");
			if (fpg == NULL) {
				fprintf (stderr, "msform: cannot open oringe file: %s\n",
						fname);
				return (0);
			}
			break;
		case 'h':
			if (strcmp(argarray[arg].stringarg,"stdin") == 0) fph = stdin;
			else {
				fph = fopen (argarray[arg].stringarg, "r");
				if (fph == NULL) {
					fprintf (stderr, "msform: cannot open 3rd polyhedron file: %s",
							argarray[arg].stringarg);
					exit (1);
				}
			}
			break;
		case 'k':
			nsector = argarray[arg].longarg;
			if (nsector > MAX_SECTOR) {
				fprintf (stderr,
					"msform: nsector (%ld) > maximum allowed (%ld)\n",
					nsector, MAX_SECTOR);
				exit (1);
			}
			break;
		case 'l':
			nsample = argarray[arg].longarg;
			break;
		case 'm':
			smooth = 1;
			minus = 1;
			break;
		case 'n':
			if (strcmp(argarray[arg].stringarg,"stdin") == 0) fpn = stdin;
			else {
				fpn = fopen (argarray[arg].stringarg, "r");
				if (fpn == NULL) {
					fprintf (stderr, "msform: cannot open 2nd polyhedron file: %s",
							argarray[arg].stringarg);
					exit (1);
				}
			}
			break;
		case 'o':
			nw = sscanf(argarray[arg].stringarg, "%s %s", fname, format);
			if (nw < 1) continue;
			if (strcmp(fname,"stdout") == 0) fpo = stdout;
			else {
				if (strcmp(format,"avs") == 0)
					fpo = fopen (fname, "wb");
				else
					fpo = fopen (fname, "w");
				if (fpo == NULL) {
					fprintf (stderr,
						"msform: cannot open output file: %s\n", fname);
					return (0);
				}
			}
			break;
		case 'p':
			strcpy (partition, argarray[arg].stringarg);
			if (strlen(partition) < 1) strcpy (partition, "maze");
			if (strcmp (partition, "maze") == 0) {
			}
			else if (strcmp (partition, "dart") == 0) {
			}
			else if (strcmp (partition, "stagger") == 0) {
			}
			else if (strcmp (partition, "spiral") == 0) {
			}
			else if (strcmp (partition, "chess") == 0) {
			}
			else if (strcmp (partition, "old") == 0) {
			}
			else {
				fprintf (stderr, "msform: invalid partition type: %s\n", partition);
				return (0);
			}
			break;
		case 'r':
			sphere_radius = argarray[arg].doublearg;
			if (sphere_radius < MINIMUM_RADIUS) sphere_radius = MINIMUM_RADIUS;
			if (sphere_radius > MAXIMUM_RADIUS) sphere_radius = MAXIMUM_RADIUS;
			break;
		case 't':
			if (strcmp(argarray[arg].stringarg,"stdin") == 0) fpt = stdin;
			else {
				fpt = fopen (argarray[arg].stringarg, "r");
				if (fpt == NULL) {
					fprintf (stderr, "msform: cannot open polyhedron file: %s\n",
							argarray[arg].stringarg);
					return (0);
				}
			}
			break;
		case 'w':
			cube_width = argarray[arg].doublearg;
			break;
		case 'x':
			x_radius = argarray[arg].doublearg;
			if (x_radius < 0.0) x_radius = 0.0;
			if (x_radius > sphere_radius) x_radius = sphere_radius;
			break;
		case 'y':
			nw = sscanf(argarray[arg].stringarg, "%s %s", fname, format);
			if (nw < 1) continue;
			if (strcmp(fname,"stdin") == 0) fpy = stdin;
			else {
				if (strcmp(format,"avs") == 0)
					fpy = fopen (fname, "rb");
				else
				fpy = fopen (fname, "r");
				if (fpy == NULL) {
					fprintf (stderr, "msform: cannot open density file: %s\n",
							fname);
					return (0);
				}
			}
			break;
		default:
			fprintf (stderr,
				"invalid flag: %s\n", argarray[arg].stringarg);
			exit (1);
		}
	}


	/* initialize scene, memory and material table */
	init_routing (fpe, fpe, fpe);
	form_mem ();
	if (error()) {
		print_error ();
		exit (1);
	}
	ms = new_msscene ();
	if (ms == NULL || error ()) {
		set_error2 ("msform: new_msscene fails");
		print_error ();
		exit (1);
	}
	if (debug) ms -> debug = debug;
	ms -> table = NULL;

	if (sphere_radius > 0.0) {
		sprintf (message, "%8.3f sphere radius", sphere_radius);
		inform (message);
	}
	if (ctrlev > 0.0) {
		sprintf (message, "%8.3f contour level", ctrlev);
		inform (message);
	}
	if (nsector > 0.0) {
		sprintf (message, "%8ld sectors", nsector);
		inform (message);
	}
	/* parameters based upon command line flags */
	
	if (fpt != NULL && fpd == NULL && cube_width == 0.0)
		result = omega_stream (ms, fpt, fpn, fph, fpe, fpo, format, sphere_radius, x_radius);
	else if (fpt != NULL || fpd != NULL) {
		result = density_stream (ms, fpt, fpd, fpy, fpg, fpe, fpo, format, sphere_radius, cube_width, ctrlev, smooth, nsector);
	}
	else {
		sprintf(message, "%8ld bins in %s partition", nbin, partition);
		inform (message);
		sprintf(message, "%8ld frequencies ", nfrequency);
		inform (message);
		result = oringe_stream (nsector, nfrequency, nbin1, nsample, fpg, fpf, fpa, fpo, partition, minus);
	}

	if (error()) {
		print_error ();
		exit (1);
	}
	free_scene (ms);
	if (error()) exit (1);
	print_counts();
	if (error()) {
		print_error();
		exit (1);
	}
	return (0);
}


/*
 * MSForm
 * Copyright 2006 by Michael L. Connolly
 * All Rights Reserved
 */
