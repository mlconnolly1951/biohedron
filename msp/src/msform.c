/* MSForm Copyright 1996 by Michael L. Connolly */
/* Last revised: December 20, 2001 */

#include "ms.h"
#include "msenum.h"
#include "mslimit.h"
#include "msdefault.h"
#include "msstruct.h"
#include "msglobal.h"
#include "msfunction.h"


int omega_stream (struct msscene *ms, FILE *fpt, FILE *fpn, FILE *fph, FILE *fpe, FILE *fpo, char *format, double sphere_radius, double x_radius)
{
	double expansion_radius, extension_radius;
	char message[128];
	struct surface *msphn1, *msphn2, *msphn3;

	if (sphere_radius <= 0.0) {
		sprintf (message, "%8.3f negative sphere radius", sphere_radius);
		set_error1 (message);
		return (0);
	}
	
	expansion_radius = 0.0;
	extension_radius = 0.0;
	if (fpt != NULL) {
		/* first polyhedron */
		msphn1 = read_vet (fpt);
		if (error()) return (0);
		fclose (fpt);
		do_bounds (msphn1);
		if (error()) return (0);
		do_axes(msphn1);
		if (error()) return (0);
	}
	else msphn1 = NULL;
	
	if (fpn != NULL) {
		/* second polyhedron */
		msphn2 = read_vet (fpn);
		if (error()) return (0);
		fclose (fpn);
		do_bounds (msphn2);
		if (error()) return (0);
		do_axes(msphn2);
		if (error()) return (0);
	}
	else msphn2 = NULL;

	if (fph != NULL) {
		/* third polyhedron */
		msphn3 = read_vet (fph);
		if (error()) return (0);
		fclose (fph);
		do_bounds (msphn3);
		if (error()) return (0);
		do_axes(msphn3);
		if (error()) return (0);
	}
	else msphn3 = NULL;

	/* evaluation loci */
	sprintf (message, "%8.3f evaluation radius", sphere_radius);
	inform(message);
	if (msphn1 != NULL && msphn2 == NULL && msphn3 == NULL) {
		expansion_radius = x_radius;
		sprintf (message, "%8.3f expansion radius", expansion_radius);
		inform(message);
	}
	if (msphn1 != NULL && msphn2 != NULL && msphn3 == NULL) {
		extension_radius = x_radius;
		sprintf (message, "%8.3f extension radius", extension_radius);
		inform(message);
	}
	set_omega_radii (msphn1, sphere_radius, expansion_radius, extension_radius);
	do_evaluation (msphn1, msphn2, msphn3, (int) (extension_radius > 0.0));
	if (error()) return (0);
	if (msphn2 == NULL && msphn3 == NULL) {
		identify (msphn1);
		if (error()) return (0);
	}
	msphn1 -> scheme = NULL;
	if (format == NULL || strlen(format) == (unsigned) 0 || strcmp(format,"vet") == 0)
		write_vet (msphn1, fpo);
	if (error()) return (0);

	if (msphn1 != NULL) {
		free_phn (msphn1);
		if (error()) return(0);
		free_object (SURFACE, (short *) msphn1);
		if (error()) return(0);
	}
	if (msphn2 != NULL) {
		free_phn (msphn2);
		if (error()) return(0);
		free_object (SURFACE, (short *) msphn2);
		if (error()) return(0);
	}
	if (msphn3 != NULL) {
		free_phn (msphn3);
		if (error()) return(0);
		free_object (SURFACE, (short *) msphn3);
		if (error()) return(0);
	}
	return (1);
}


int density_stream (struct msscene *ms, FILE *fpt, FILE *fpd, FILE *fpy, FILE *fpg, FILE *fpe, FILE *fpo, char *format, double sphere_radius, double cube_width, double ctrlev, int smooth, long nsector)
{
	int result;
	char message[MAX_STRING];

	if (fpt != NULL && fpd != NULL && fpg != NULL) {
		if (fpo == NULL) {
			set_error1 ("msform: missing output file name");
			return (0);
		};
		result = msorange (ms, fpt,fpd,fpg,fpe,fpo);
		if (error()) return (0);
		if (!result) return (0);
	}
	else if (fpt != NULL && fpd != NULL && fpy != NULL) {
		if (fpo == NULL) {
			set_error1 ("msform: missing output file name");
			return (0);
		};
		if (sphere_radius <= 0.0) {
			sprintf (message, "%8.3f negative sphere radius", sphere_radius);
			set_error1 (message);
			return (0);
		}
		result = mspolden2 (ms, fpt,fpd,fpy,fpe,fpo, sphere_radius, nsector);
		if (error()) return (0);
		if (!result) return (0);
	}
	else if (fpt != NULL && fpd != NULL) {
		if (fpo == NULL) {
			set_error1 ("msform: missing output file name");
			return (0);
		};
		if (sphere_radius <= 0.0) {
			sprintf (message, "%8.3f negative sphere radius", sphere_radius);
			set_error1 (message);
			return (0);
		}
		result = mspolden (ms, fpt,fpd,fpe,fpo, sphere_radius, nsector, smooth);
		if (error()) return (0);
		if (!result) return (0);
	}
	else if (fpd != NULL && fpy != NULL) {
		if (fpo == NULL) {
			set_error1 ("msform: missing output file name");
			return (0);
		};
		if (sphere_radius <= 0.0) {
			sprintf (message, "%8.3f negative sphere radius", sphere_radius);
			set_error1 (message);
			return (0);
		}
		result = mspol2 (ms, fpd,fpy,fpe,fpo,format, ctrlev, sphere_radius, smooth);
		if (error()) return (0);
		if (!result) return (0);
	}
	else if (fpt != NULL) {
		if (fpo == NULL) {
			set_error1 ("msform: missing output file name");
			return (0);
		};
		result = msden (ms,fpt,fpe,fpo,cube_width,format);
		if (error()) return (0);
		if (!result) return (0);
	}
	else if (fpd != NULL) {
		if (fpo == NULL) {
			set_error1 ("msform: missing output file name");
			return (0);
		};
		if (sphere_radius <= 0.0) {
			sprintf (message, "%8.3f negative sphere radius", sphere_radius);
			set_error1 (message);
			return (0);
		}
		result = mspol (ms, fpd,fpe,fpo,format, ctrlev, sphere_radius, smooth);
		if (error()) return (0);
		if (!result) return (0);
	}
	if (error()) return (0);

	return (1);
}

int msden (struct msscene *ms, FILE *fpt, FILE *fpe, FILE *fpo, double cube_width, char *format)
{
	int result;
	char message[MAXLINE];
	struct surface *den;
	struct surface *phn;

	sprintf (message, "%8.3f cube width", cube_width);
	inform (message);
	if (format == NULL) informd ("msform: no output format specified");
	phn = read_polyhedron (fpt);
	if (error ()) return (0);
	if (phn == NULL) {
		set_error2 ("msform: problem reading polyhedron");
		return (0);
	}
	ms -> this_srf = phn;
	den = polyhedron_to_density (phn, cube_width);
	if (error ()) return (0);
	if (den == NULL) return (0);
	phn -> next = den;
	result = write_density (den, fpo);
	if (error ()) return (0);
	if (result == 0) return (0);
	fclose(fpo);
	return(1);
}

int mspol (struct msscene *ms, FILE *fpd, FILE *fpe, FILE *fpo, char *format, double ctrlev, double sphere_radius, int smooth)
{
	int result;
	char message[MAXLINE];
	struct surface *den;
	struct surface *phn;

	sprintf (message, "%8.3f sphere radius", sphere_radius);
	inform (message);
	sprintf (message, "%8.3f contour level", ctrlev);
	inform (message);
	if (format == NULL) informd ("msform: no output format specified");
	if (smooth) inform ("         smoothed density contoured");
	den = read_density (fpd);
	if (den == NULL) {
		set_error1 ("mspol: problem reading density");
		return (0);
	}
	phn = density_to_polyhedron (den, ctrlev, sphere_radius, smooth);
	if (phn == NULL) return (0);
	result = write_vet (phn, fpo);
	if (result == 0) return (0);
	fclose(fpo);
	return(1);
}

int mspol2 (struct msscene *ms, FILE *fpd, FILE *fpy, FILE *fpe, FILE *fpo, char *format, double ctrlev, double sphere_radius, int smooth)
{
	int result;
	char message[MAXLINE];
	struct surface *den, *den1, *den2;
	struct surface *phn;

	if (format == NULL) informd ("msform: no output format specified");
	sprintf (message, "%8.3f sphere radius", sphere_radius);
	inform (message);
	sprintf (message, "%8.3f contour level", ctrlev);
	inform (message);
	if (smooth) inform ("         smoothed density contoured");
	den1 = read_density (fpd);
	if (den1 == NULL) {
		set_error1 ("mspol2: problem reading density");
		return (0);
	}
	den2 = read_density (fpy);
	if (den2 == NULL) {
		set_error1 ("mspol2: problem reading density");
		return (0);
	}
	den = subtract_densities (den1, den2);
	if (den == NULL) return (0);
	phn = density_to_polyhedron (den, ctrlev, sphere_radius, smooth);
	if (phn == NULL) return (0);
	result = write_vet (phn, fpo);
	if (result == 0) return (0);
	fclose(fpo);
	return(1);
}

int mspolden (struct msscene *ms, FILE *fpt, FILE *fpd, FILE *fpe, FILE *fpo, double sphere_radius, long nsector, int smooth)
{
	int result, action;
	char message[MAXLINE];
	struct surface *den;
	struct surface *phn;

	sprintf (message, "%8.3f sphere radius", sphere_radius);
	inform (message);
	phn = read_outer (fpt);
	if (phn == NULL) {
		set_error1 ("mspolden: problem reading polyhedron");
		return (0);
	}
	den = read_density (fpd);
	if (den == NULL) {
		set_error1 ("mspolden: problem reading density");
		return (0);
	}
	if (smooth) {
		action = GET_NORMAL;
		result = polyhedron_and_density (phn, den, sphere_radius, 0, action);
		if (result == 0) return (0);
	}
	else {
		outward_to_normal (phn);
		if (error()) return (0);
	}
	action = GET_FRAME;
	result = polyhedron_and_density (phn, den, sphere_radius, 0, action);
	if (result == 0) return (0);
	if (nsector == 0) {
		result = write_vet (phn, fpo);
		if (result == 0) return (0);
	}
	else {
		action = GET_ORINGE;
		result = polyhedron_and_density (phn, den, sphere_radius, nsector, action);
		if (result == 0) return (0);
		result = write_phnorg (phn, sphere_radius, nsector, fpo);
		if (result == 0) return (0);
	}
	fclose(fpo);
	return(1);
}

int mspolden2 (struct msscene *ms, FILE *fpt, FILE *fpd, FILE *fpy, FILE *fpe, FILE *fpo, double sphere_radius, long nsector)
{
	int result, action;
	char message[MAXLINE];
	struct surface *densum, *dendiff, *den1, *den2;
	struct surface *phn;

	sprintf (message, "%8.3f sphere radius", sphere_radius);
	inform (message);
	phn = read_polyhedron (fpt);
	if (phn == NULL) {
		set_error1 ("mspolden2: problem reading polyhedron");
		return (0);
	}
	den1 = read_density (fpd);
	if (den1 == NULL) {
		set_error1 ("mspolden2: problem reading density");
		return (0);
	}
	den2 = read_density (fpy);
	if (den2 == NULL) {
		set_error1 ("mspolden2: problem reading density");
		return (0);
	}
	dendiff = subtract_densities (den1, den2);
	if (dendiff == NULL) return (0);
	action = GET_NORMAL;
	result = polyhedron_and_density (phn, dendiff, sphere_radius, nsector, action);
	if (result == 0) return (0);
	densum = add_densities (den1, den2);
	if (densum == NULL) return (0);
	action = GET_FOURIER0;
	result = polyhedron_and_density (phn, den1, sphere_radius, nsector, action);
	if (result == 0) return (0);
	action = GET_FOURIER1;
	result = polyhedron_and_density (phn, den2, sphere_radius, nsector, action);
	if (result == 0) return (0);
	action = GET_LOCAL;
	result = polyhedron_and_density (phn, densum, sphere_radius, nsector, action);
	if (result == 0) return (0);
	result = write_vet (phn, fpo);
	if (result == 0) return (0);
	fclose(fpo);
	return(1);
}


struct surface *density_to_polyhedron (struct surface *den, double ctrlev, double sphere_radius, int smooth)
{
	int j, k;
	long i, dimension;
	long vn0, vn1, vn2, en0, en1, en2, uen0, uen1, uen2;
	long dimensions[3], origin[3], width[3], limit[3];
	long n_vertex, n_edge, n_triangle;
	long *dedg, *dtri;
	float *densities, *dvtx, *dnml;
	double bounds[2][3];
	double cube_width;
	struct phnvtx *vtx;
	struct phnedg *edg;
	struct phntri *tri;
	char message[128];
	struct msdata *msd;
	struct surface *phn;
	struct vanity *vanities;

	vanities = NULL;
	cube_width = den -> cube_width;
	for (j = 0; j < 3; j++) {
		origin[j] = den -> origin[j];
		width[j] = den -> width[j];
		limit[j] = den -> limit[j];
	}
	densities = den -> densities;
	for (k = 0; k < 3; k++) {
		bounds[0][k] = origin[k] * cube_width;
		bounds[1][k] = limit[k] * cube_width;
		dimensions[k] = width[k];
	}
	dimension = 3;
	vanities = smooth_density (dimension, dimensions, densities, (float) sphere_radius, (float) cube_width, smooth);
	msd = doXelp (dimension, dimensions, bounds, ctrlev, vanities);
	if (msd == NULL) {
		set_error1 ("density_to_polyhedron: doXelp returns null");
		return (NULL);
	}
	if (msd -> narray != 4) {
		set_error1 ("density_to_polyhedron: doXelp returns wrong number of arrays");
		return (NULL);
	}
	/* store msd data into polyhedron */
	n_vertex = msd -> counts[1];
	n_edge = msd -> counts[2];
	n_triangle = msd -> counts[3];
	if (n_vertex <= 0) {
		set_error1 ("(density_to_polyhedron): bad n_vertex");
		return (NULL);
	}
	if (n_edge <= 0) {
		set_error1 ("(density_to_polyhedron): bad n_edge");
		return (NULL);
	}
	if (n_triangle <= 0) {
		set_error1 ("(density_to_polyhedron): bad n_triangle");
		return (NULL);
	}
	phn = init_phn (n_vertex, n_edge, n_triangle);
	if (phn == NULL) {
		set_error1 ("density_to_polyhedron: no memory for polyhedron");
		return (NULL);
	}
	/* skip over plexlines */
	phn -> format = 1;

	/* initialize bounds */
	for (k = 0; k < 3; k++) {
		phn -> bounds[0][k] =  1000000.0;
		phn -> bounds[1][k] = -1000000.0;
	}
	phn -> minvals[0] =  ctrlev;
	phn -> maxvals[0] =  ctrlev;

	for (i = 0; i < phn -> n_phnvtx; i++) {
		dvtx = msd -> floats[1] + 6 * i;
		dnml = dvtx + 3;
		vtx = *(phn -> phnvtx_handles + i);
		for (k = 0; k < 3; k++) {
			vtx -> center[k] = *(dvtx+k);
			vtx -> outward[k] = *(dnml+k);
		}
		normalize (vtx -> outward);
		for (k = 0; k < 3; k++) {
			if (vtx -> center[k] < phn -> bounds[0][k])
				phn -> bounds[0][k] = vtx -> center[k];
			if (vtx -> center[k] > phn -> bounds[1][k])
				phn -> bounds[1][k] = vtx -> center[k];
		}
	}
	for (k = 0; k < 3; k++)
		phn -> center[k] = (phn -> bounds[0][k] + phn -> bounds[1][k]) / 2.0;

	for (i = 0; i < phn -> n_phnedg; i++) {
		dedg = msd -> longs[2] + 2 * i;
		edg = *(phn -> phnedg_handles + i);
		vn0 = *dedg;
		vn1 = *(dedg+1);
		if (vn0 < 1 || vn0 > phn -> n_phnvtx) {
			sprintf (message, "(density_to_polyhedron): bad vertex number: %6ld", vn0);
			set_error1(message);
			return (NULL);
		}
		if (vn1 < 1 || vn1 > phn -> n_phnvtx) {
			sprintf (message, "(density_to_polyhedron): bad vertex number: %6ld", vn1);
			set_error1(message);
			return (NULL);
		}
		edg -> vtxnum[0] = vn0;
		edg -> vtxnum[1] = vn1;
	}

	for (i = 0; i < phn -> n_phntri; i++) {
		tri = *(phn -> phntri_handles + i);
		dtri = msd -> longs[3] + 6 * i;
		en0 = *dtri;
		en1 = *(dtri+1);
		en2 = *(dtri+2);
		vn0 = *(dtri+3);
		vn1 = *(dtri+4);
		vn2 = *(dtri+5);
		uen0 = abs (en0);
		uen1 = abs (en1);
		uen2 = abs (en2);
		if (uen0 < 1 || uen0 > phn -> n_phnedg) {
			set_error1 ("(density_to_polyhedron): bad edge number");
			return (NULL);
		}
		if (uen1 < 1 || uen1 > phn -> n_phnedg) {
			set_error1 ("(density_to_polyhedron): bad edge number");
			return (NULL);
		}
		if (uen2 < 1 || uen2 > phn -> n_phnedg) {
			set_error1 ("(density_to_polyhedron): bad edge number");
			return (NULL);
		}
		if (vn0 < 1 || vn0 > phn -> n_phnvtx) {
			sprintf (message, "(density_to_polyhedron): bad vertex number: %6ld", vn0);
			set_error1(message);
			return (NULL);
		}
		if (vn1 < 1 || vn1 > phn -> n_phnvtx) {
			sprintf (message, "(density_to_polyhedron): bad vertex number: %6ld", vn1);
			set_error1(message);
			return (NULL);
		}
		if (vn2 < 1 || vn2 > phn -> n_phnvtx) {
			sprintf (message, "(density_to_polyhedron): bad vertex number: %6ld", vn2);
			set_error1(message);
			return (NULL);
		}
		tri -> edgnum[0] = en0;
		tri -> edgnum[1] = en1;
		tri -> edgnum[2] = en2;
		tri -> vtxnum[0] = vn0;
		tri -> vtxnum[1] = vn1;
		tri -> vtxnum[2] = vn2;
	}

	if (!freemsdata (msd)) return (NULL);
	free_objects (VANITY, (short *) vanities);

	return (phn);
}

int polyhedron_and_density (struct surface *phn, struct surface *den, double sphere_radius, long nsector, int action)
{
	int do_fourier;
	int j, k;
	long v, n_vertex, dimension;
	long dimensions[3], origin[3], width[3], limit[3];
	float *densities, *dfloat;
	unsigned char *dsb;
	double cube_width;
	double *vertices, *normals, *bases, *dnml, *dvtx, *dbas;
	double bounds[2][3];
	unsigned char *byte_density;
	struct msdata *msd;
	struct phnvtx *vtx;

	cube_width = den -> cube_width;
	n_vertex = phn -> n_phnvtx;
	vertices = allocate_doubles (n_vertex * 3, 0, VERTICES);
	if (vertices == NULL) {
		set_error1 ("polyhedron_and_density: not enough memory for vertices");
		return (0);
	}
	normals = allocate_doubles (n_vertex * 3, 0, NORMALSV);
	if (normals == NULL) {
		set_error1 ("polyhedron_and_density: not enough memory and normals");
		return (0);
	}
	bases = allocate_doubles (n_vertex * 3, 0, BASES);
	if (bases == NULL) {
		set_error1 ("polyhedron_and_density: not enough memory for bases");
		return (0);
	}
	for (j = 0; j < 3; j++) {
		origin[j] = den -> origin[j];
		width[j] = den -> width[j];
		limit[j] = den -> limit[j];
	}
	densities = den -> densities;
	for (k = 0; k < 3; k++) {
		bounds[0][k] = origin[k] * cube_width;
		bounds[1][k] = limit[k] * cube_width;
		dimensions[k] = width[k];
	}
	dimension = 3;
	/* make vertex coordinates relative to origin */
	for (v = 0; v < n_vertex; v++) {
		vtx = *(phn -> phnvtx_handles + v);
		dvtx = vertices + 3 * v;
		dnml = normals + 3 * v;
		dbas = bases + 3 * v;
		for (k = 0; k < 3; k++) {
			*(dvtx + k) = vtx -> center[k] - bounds[0][k];
			*(dnml + k) = vtx -> normal[k];	/* use smoothed normal */
			*(dbas + k) = vtx -> base[k];	/* use two-fold base */
		}
	}
	switch (action) {
	case GET_NORMAL:
	case GET_LOCAL:
		do_fourier = 0;
		break;
	case GET_FOURIER0:
	case GET_FOURIER1:
	case GET_FRAME:
	case GET_ORINGE:
		do_fourier = 1;
		break;
	default:
		do_fourier = 0;
		break;
	}
	msd = eval_density (dimension, dimensions, vertices, normals, bases, n_vertex, densities, (float) sphere_radius, (float) cube_width, do_fourier, nsector);
	if (msd == NULL) return (0);
	if (msd == NULL) {
		set_error2 ("polyhedron_and_density: eval_density returns null msd");
		return (0);
	}
	
	if (msd -> counts[0] != n_vertex ||
		msd -> widths[0] != 6 || msd -> floats[0] == NULL) {
		set_error2 ("polyhedron_and_density: eval_density returns invalid msd (no floats)");
		return (0);
	}
	if (action == GET_ORINGE) {
		if (nsector == 0 || msd -> narray < 2) {
			set_error2 ("polyhedron_and_density: eval_density returns invalid msd (no bytes)");
			return (0);
		}
		if (msd -> counts[1] != n_vertex ||
			msd -> widths[1] != nsector || msd -> bytes[1] == NULL) {
			set_error2 ("polyhedron_and_density: eval_density returns invalid msd (no sector bytes)");
			return (0);
		}
	}
	for (v = 0; v < phn -> n_phnvtx; v++) {
		vtx = *(phn -> phnvtx_handles + v);
		dfloat = msd -> floats[0] + 6 * v;
		switch (action) {
		case GET_NORMAL:
			for (k = 0; k < 3; k++)
				vtx -> normal[k] = *(dfloat+3+k);
			normalize (vtx -> normal);
			break;
		case GET_FOURIER0:
			vtx -> values[0] = *(dfloat + 1);
			break;
		case GET_FOURIER1:
			vtx -> values[1] = *(dfloat + 1);
			break;
		case GET_LOCAL:
			vtx -> values[2] = *(dfloat);
			break;
		case GET_FRAME:
			vtx -> values[0] = *(dfloat);
			vtx -> values[1] = *(dfloat + 1);
			vtx -> values[2] = *(dfloat + 2);
			/* normal vector already gotten */
			for (k = 0; k < 3; k++)
				vtx -> base[k] = *(dfloat+3+k);
			normalize (vtx -> base);
			cross (vtx -> normal, vtx -> base, vtx -> zenith);
			break;
		case GET_ORINGE:
			dsb = msd -> bytes[1] + nsector * v;
			byte_density = allocate_bytes (nsector);
			if (byte_density == NULL) return (0);
			for (j = 0; j < nsector; j++)
				*(byte_density+j) = *(dsb+j);
			vtx -> byte_density = byte_density;
			break;
		default:
			break;
		}
	}
	if (!freemsdata (msd)) return (0);
	free_doubles (vertices, 0, VERTICES);
	free_doubles (normals, 0, NORMALSV);
	free_doubles (bases, 0, BASES);
	return (1);
}


void form_mem ()
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
	
	type = SPHERE;
	size = sizeof (struct sphere);
	strcpy (type_name, "sphere");
	define_type (type, size, type_name);
	if (error()) return;
	
	type = TORUS;
	size = sizeof (struct torus);
	strcpy (type_name, "torus");
	define_type (type, size, type_name);
	if (error()) return;
	
	type = PROBE;
	size = sizeof (struct probe);
	strcpy (type_name, "probe");
	define_type (type, size, type_name);
	if (error()) return;
	
	type = CENTRAL;
	size = sizeof (struct central);
	strcpy (type_name, "central");
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
		
	type = COMPONENT;
	size = sizeof (struct component);
	strcpy (type_name, "component");
	define_type (type, size, type_name);
	if (error()) return;
	
	type = CHUNK;
	size = sizeof (struct chunk);
	strcpy (type_name, "chunk");
	define_type (type, size, type_name);
	if (error()) return;

	type = SURFACE;
	size = sizeof (struct surface);
	strcpy (type_name, "surface");
	define_type (type, size, type_name);
	if (error()) return;

	type = VARIETY;
	size = sizeof (struct variety);
	strcpy (type_name, "variety");
	define_type (type, size, type_name);
	if (error()) return;

	type = PHNVTX;
	size = phnvtx_size ();
	strcpy (type_name, "phnvtx");
	define_type (type, size, type_name);

	type = PHNEDG;
	size = phnedg_size ();
	strcpy (type_name, "phnedg");
	define_type (type, size, type_name);

	type = PHNTRI;
	size = phntri_size ();
	strcpy (type_name, "phntri");
	define_type (type, size, type_name);

	type = POLYGON;
	size = sizeof (struct polygon);
	strcpy (type_name, "polygon");
	define_type (type, size, type_name);

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

	type = OBJECT_SCHEME;
	size = sizeof (struct object_scheme);
	strcpy (type_name, "object_scheme");
	define_type (type, size, type_name);

	type = MOLECULE;
	size = sizeof (struct molecule);
	strcpy (type_name, "molecule");
	define_type (type, size, type_name);

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

	type = COLOR_RAMP;
	size = sizeof (struct color_ramp);
	strcpy (type_name, "color_ramp");
	define_type (type, size, type_name);
	if (error()) return;

	type = PLEXLINK;
	size = sizeof (struct PlexLink);
	strcpy (type_name, "PlexLink");
	define_type (type, size, type_name);
	if (error()) return;

	type = SUBEDGE;
	size = sizeof (struct subedge);
	strcpy (type_name, "subedge");
	define_type (type, size, type_name);
	if (error()) return;

	type = SUBPOLYGON;
	size = sizeof (struct subpolygon);
	strcpy (type_name, "subpolygon");
	define_type (type, size, type_name);
	if (error()) return;

	type = ORINGE;
	size = sizeof (struct oringe);
	strcpy (type_name, "oringe");
	define_type (type, size, type_name);
	if (error()) return;

	type = SECTOR;
	size = sizeof (struct sector);
	strcpy (type_name, "sector");
	define_type (type, size, type_name);
	if (error()) return;

	type = TERM;
	size = sizeof (struct term);
	strcpy (type_name, "term");
	define_type (type, size, type_name);
	if (error()) return;

	type = BAR;
	size = sizeof (struct bar);
	strcpy (type_name, "bar");
	define_type (type, size, type_name);
	if (error()) return;

	type = RECORD;
	size = sizeof (struct record);
	strcpy (type_name, "record");
	define_type (type, size, type_name);
	if (error()) return;

	type = TOKEN;
	size = sizeof (struct token);
	strcpy (type_name, "token");
	define_type (type, size, type_name);
	if (error()) return;

	type = EDGER;
	size = sizeof (struct edger);
	strcpy (type_name, "edger");
	define_type (type, size, type_name);
	if (error()) return;
	
	type = PLEX;
	size = sizeof (struct Plex);
	strcpy (type_name, "Plex");
	define_type (type, size, type_name);
	if (error()) return;
	
	type = PLEXI;
	size = sizeof (struct Plexi);
	strcpy (type_name, "Plexi");
	define_type (type, size, type_name);
	if (error()) return;
	
	type = PLEXVERTEX;
	size = sizeof (struct PlexVertex);
	strcpy (type_name, "PlexVertex");
	define_type (type, size, type_name);
	if (error()) return;
	
	type = PLEXEDGE;
	size = sizeof (struct PlexEdge);
	strcpy (type_name, "PlexEdge");
	define_type (type, size, type_name);
	if (error()) return;
	
	type = PLEXTRIANGLE;
	size = sizeof (struct PlexTriangle);
	strcpy (type_name, "PlexTriangle");
	define_type (type, size, type_name);
	if (error()) return;
	
	type = PLEXCUBE;
	size = sizeof (struct PlexCube);
	strcpy (type_name, "PlexCube");
	define_type (type, size, type_name);
	if (error()) return;
	
	type = PLEXJOIN;
	size = sizeof (struct PlexJoin);
	strcpy (type_name, "PlexJoin");
	define_type (type, size, type_name);
	if (error()) return;
	
	type = PROVINCE;
	size = sizeof (struct Province);
	strcpy (type_name, "Province");
	define_type (type, size, type_name);
	if (error()) return;
	
	type = GLASS;
	size = sizeof (struct Glass);
	strcpy (type_name, "Glass");
	define_type (type, size, type_name);
	if (error()) return;
	
	type = GLASS_BLOCK;
	size = sizeof (struct GlassBlock);
	strcpy (type_name, "GlassBlock");
	define_type (type, size, type_name);
	if (error()) return;
	
	type = VANITY;
	size = sizeof (struct vanity);
	strcpy (type_name, "vanity");
	define_type (type, size, type_name);
	if (error()) return;
	
	type = DENVTX;
	size = sizeof (struct denvtx);
	strcpy (type_name, "denvtx");
	define_type (type, size, type_name);
	if (error()) return;
	
	type = MSDATA;
	size = sizeof (struct msdata);
	strcpy (type_name, "msdata");
	define_type (type, size, type_name);
	if (error()) return;
	
	type = CRITLINK;
	size = sizeof (struct critlink);
	strcpy (type_name, "critlink");
	define_type (type, size, type_name);
	if (error()) return;
	
	type = CEPT;
	size = sizeof (struct cept);
	strcpy (type_name, "cept");
	define_type (type, size, type_name);
	if (error()) return;
	
}

void outward_to_normal (struct surface *phn)
{
	int k;
	long v;
	char message[MAX_STRING];
	struct phnvtx *cv;
	
	for (v = 0; v < phn -> n_phnvtx; v++) {
		cv = num2phnvtx (phn, v + 1);
		if (cv == NULL)  {
			set_error1 ("outward_to_normal: invalid vertex number");
			return;
		}
		for (k = 0; k < 3; k++)
			cv -> normal[k] = cv -> outward[k];
	}
	sprintf(message,
		"%8ld vertices using input outward vector as normal vector",
		phn -> n_phnvtx);
	inform (message);
}


/*
 * MSForm
 * Copyright 1986 by Michael L. Connolly
 * All Rights Reserved
 */
