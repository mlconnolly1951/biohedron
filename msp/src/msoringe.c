/* Molecular Surface Package */
/* Copyright 1995 by Michael L. Connolly */
/* February 8, 2000 */

#include "ms.h"
#include "msenum.h"
#include "mslimit.h"
#include "msdefault.h"
#include "msstruct.h"
#include "msglobal.h"
#include "msfunction.h"

/* Steps: (1) compute normal vector and frame,
          (2) compute sector densities,
          (3) set vertex values,
          (4) rotate polyhedron according to frame.
*/

int oringe_stream (long nsector, long nfrequency, long nbin1, long nsample, FILE *fpg, FILE *fpf, FILE *fpa, FILE *fpo, char *partition, int minus)
{
	long counter = 0;
	char message[128];
	struct oringe *iring, *jring, *kring;
	struct oringe *head_iring, *tail_iring, *next_iring;
	struct oringe *head_jring, *tail_jring, *next_jring;

	if (fpg != NULL && fpo != NULL && nsample > 0) {
		sprintf (message, "%8ld samples for abstract oringes", nsample);
		inform (message);
		while (1) {
			iring = read_oringe (fpg);
			if (error()) return (0);
			if (iring == NULL) break;
			if (minus) complement_standard (iring);
			if (error()) return (0);
			allocate_fourier (iring, nfrequency);
			if (error ()) return (0);
			standard_to_abstract (iring, nsample, nbin1, partition);
			if (error()) return (0);
			write_oringe (iring, fpo, "abstract", partition);
			if (error()) return (0);
			deep_free_oringe (iring);
			counter++;
		}
		fclose (fpg);
		fclose (fpo);
	}
	else if (fpg != NULL && fpo != NULL) {
		sprintf (message, "%8ld frequencies for fourier oringes", nfrequency);
		inform (message);
		while (1) {
			iring = read_oringe (fpg);
			if (error()) return (0);
			if (iring == NULL) break;
			allocate_fourier (iring, nfrequency);
			if (error ()) return (0);
			standard_to_fourier (iring, 0L);
			if (error()) return (0);
			write_oringe (iring, fpo, "fourier", partition);
			if (error()) return (0);
			deep_free_oringe (iring);
			counter++;
		}
		fclose (fpg);
		fclose (fpo);
	}
	else if (fpf != NULL && fpo != NULL && nsector > 0) {
		sprintf (message, "%8ld sectors for back transformed oringes", nsector);
		inform (message);
		while (1) {
			iring = read_oringe (fpf);
			if (error()) return (0);
			if (iring == NULL) break;
			allocate_standard (iring, nsector);
			if (error ()) return (0);
			fourier_to_standard (iring);
			if (error()) return (0);
			write_oringe (iring, fpo, "standard", partition);
			if (error()) return (0);
			deep_free_oringe (iring);
			counter++;
		}
		fclose (fpf);
		fclose (fpo);
	}
	sprintf (message, "%8ld oringes written", counter);
	inform (message);
	return (1);
}


int msorange (struct msscene *ms, FILE *fpt, FILE *fpd, FILE *fpg, FILE *fpe, FILE *fpo)
{
	int result, whichval;
	char message[MAXLINE];
	struct surface *den;
	struct surface *phn;
	struct oringe *org;

	phn = read_outer (fpt);
	if (phn == NULL) {
		set_error2 ("msorange: problem reading polyhedron");
		return (0);
	}
	den = read_density (fpd);
	if (den == NULL) {
		set_error2 ("msorange: problem reading density");
		return (0);
	}
	org = read_oringe (fpg);
	if (org == NULL) {
		set_error2 ("msorange: problem reading oringe");
		return (0);
	}
	whichval = 0;
	make_phn_simple (phn);
	if (error ()) return (0);
	result = compute_orange_density (phn, den, org);
	if (result == 0) return (0);
	rotate_polyhedron (org, phn);
	if (error ()) return (0);
	result = write_vet (phn, fpo);
	if (result == 0) return (0);
	fclose(fpo);
	return(1);
}

int compute_orange_density (struct surface *phn, struct surface *den, struct oringe *org) {
	int k;
	long n_vertex, v;
	float *values;
	double *vertices;
	char message[MAXLINE];

	values = (float *) NULL;
	n_vertex = phn -> n_phnvtx;
	vertices = phn -> vertex_centers;
	if (vertices == NULL) {
		set_error2 ("compute_orange_density: null polyhedron vertex_centers");
		return (0);
	}
	compute_orange_frame (den, org);
	if (error ()) return (0);
	compute_sector_densities (den, org);
	if (error ()) return (0);
	values = set_orange_values (org, phn);
	if (error ()) return (0);
	if (values == NULL) {
		set_error2 ("compute_orange_density: no vertex values returned by set_orange_values");
		return (0);
	}
	store_orange_values (values, org, phn);
	if (error ()) return (0);
	return (1);
}

int compute_orange_frame (struct surface *denptr, struct oringe *org) {
	int  k;
	long jdx, nx, ny, nz, x, y, z, nden;
	long xx, yy, zz, dx, dy, dz, iradius, v;
	double r2, d2, factor, d1, delta, hdelta, hwidth;
	double gdx, gdy, gdz;
	double den, numerator, denominator;
	double xden, yden, zden;
	double xnumerator, ynumerator, znumerator;
	double concavity, sphere_radius, cube_width;
	double fx, fy, fz, fdx, fdy, fdz;
	double axis[3], base[3], zenith[3];
	char message[MAXLINE];
	float *densities;

	nx = denptr -> width[0];
	ny = denptr -> width[1];
	nz = denptr -> width[2];
	nden = nx * ny * nz;
	if (nden <= 0) return (0);
	sphere_radius = org -> radius;
	cube_width = denptr -> cube_width;
	iradius = floor (sphere_radius / cube_width) + 2.5;
	r2 = (sphere_radius * sphere_radius);
	delta = DELTA_WIDTH;
	hdelta = delta/2;
	hwidth = cube_width/2;
	densities = denptr -> densities;
	fx = org -> center[0] - denptr -> origin[0] * denptr -> cube_width;
	fy = org -> center[1] - denptr -> origin[1] * denptr -> cube_width;
	fz = org -> center[2] - denptr -> origin[2] * denptr -> cube_width;
	x = floor (fx / cube_width);
	y = floor (fy / cube_width);
	z = floor (fz / cube_width);
	/* compute normal vector and local density */
	concavity = 0.0;
	numerator = 0.0; denominator = 0.0;
	xnumerator = 0.0; ynumerator = 0.0; znumerator = 0.0;
	for (k = 0; k < 3; k++) {
		base[k] = 0.0;
		zenith[k] = 0.0;
		axis[k] = 0.0;
	}
	for (dz = -iradius; dz <= iradius; dz++) {
		zz = z + dz;
		if (zz < 0 || zz >= nz) continue;
		fdz = (zz+0.5) * cube_width - fz;
		for (dy = -iradius; dy <= iradius; dy++) {
			yy = y + dy;
			if (yy < 0 || yy >= ny) continue;
			fdy = (yy+0.5) * cube_width - fy;
			for (dx = -iradius; dx <= iradius; dx++) {
				xx = x + dx;
				if (xx < 0 || xx >= nx) continue;
				jdx = zz * ny * nx + yy * nx + xx;
				den = *(densities+jdx);
				fdx = (xx+0.5) * cube_width - fx;
				for (gdz = fdz-hwidth + hdelta; gdz < fdz+hwidth; gdz += delta)
					for (gdy = fdy-hwidth + hdelta; gdy < fdy+hwidth; gdy += delta)
						for (gdx = fdx-hwidth + hdelta; gdx < fdx+hwidth; gdx += delta) {
							d2 = gdx * gdx + gdy * gdy + gdz * gdz;
							factor = 1.0 - d2/r2;
							if (factor <= 0.0) continue;
							factor = 1.0;
							d1 = sqrt (d2);
							if (d1 <= 0.0) continue;
							numerator += factor * den;
							denominator += factor;
							xden = gdx * den;
							xnumerator += factor * xden;
							yden = gdy * den;
							ynumerator += factor * yden;
							zden = gdz * den;
							znumerator += factor * zden;
						}
			}
		}
	}
	if (denominator > 0.0) {
		/* use negative of gradient to get outward-pointing vector */
		axis[0] = -(xnumerator/denominator);
		axis[1] = -(ynumerator/denominator);
		axis[2] = -(znumerator/denominator);
		concavity = (numerator/denominator);
	}
	normalize (axis);
	arbprp (axis, base);
	cross (axis, base, zenith);
	org -> concavity = concavity;
	for (k = 0; k < 3; k++) {
		org -> frame[0][k] = base[k];
		org -> frame[1][k] = zenith[k];
		org -> frame[2][k] = axis[k];
	}
	return (1);
}

int compute_sector_densities (struct surface *denptr, struct oringe *org) {
	int  i, k;
	long jdx, nx, ny, nz, x, y, z, nden;
	long isector, nsector;
	long xx, yy, zz, dx, dy, dz, iradius, v;
	double r2, d2, factor, d1, delta, hdelta, hwidth;
	double gdx, gdy, gdz;
	double den;
	double xden, yden, zden;
	double sphere_radius, cube_width;
	double c, s;
	double fx, fy, fz, fdx, fdy, fdz;
	double angle;
	double axis[3], base[3], zenith[3];
	char message[MAXLINE];
	float *densities;
	double sector_numerator[MAX_SECTOR];
	double sector_denominator[MAX_SECTOR];

	nx = denptr -> width[0];
	ny = denptr -> width[1];
	nz = denptr -> width[2];
	nden = nx * ny * nz;
	if (nden <= 0) return (0);
	sphere_radius = org -> radius;
	cube_width = denptr -> cube_width;
	iradius = floor (sphere_radius / cube_width) + 2.5;
	r2 = (sphere_radius * sphere_radius);
	delta = DELTA_WIDTH;
	hdelta = delta/2;
	hwidth = cube_width/2;
	densities = denptr -> densities;
	nsector = org -> nsector;
	fx = org -> center[0] - denptr -> origin[0] * denptr -> cube_width;
	fy = org -> center[1] - denptr -> origin[1] * denptr -> cube_width;
	fz = org -> center[2] - denptr -> origin[2] * denptr -> cube_width;
	x = floor (fx / cube_width);
	y = floor (fy / cube_width);
	z = floor (fz / cube_width);
	for (i = 0; i < org -> nsector; i++) {
		org -> sector_density[i] = 0.0;
		sector_numerator[i] = 0.0;
		sector_denominator[i] = 0.0;
	}
	/* compute sector densities */
	for (k = 0; k < 3; k++) {
		base[k] = org -> frame[0][k];
		zenith[k] = org -> frame[1][k];
		axis[k] = org -> frame[2][k];
	}
	for (dz = -iradius; dz <= iradius; dz++) {
		zz = z + dz;
		if (zz < 0 || zz >= nz) continue;
		fdz = (zz+0.5) * cube_width - fz;
		for (dy = -iradius; dy <= iradius; dy++) {
			yy = y + dy;
			if (yy < 0 || yy >= ny) continue;
			fdy = (yy+0.5) * cube_width - fy;
			for (dx = -iradius; dx <= iradius; dx++) {
				xx = x + dx;
				if (xx < 0 || xx >= nx) continue;
				jdx = zz * ny * nx + yy * nx + xx;
				den = *(densities+jdx);
				fdx = (xx+0.5) * cube_width - fx;
				for (gdz = fdz-hwidth + hdelta; gdz < fdz+hwidth; gdz += delta)
					for (gdy = fdy-hwidth + hdelta; gdy < fdy+hwidth; gdy += delta)
						for (gdx = fdx-hwidth + hdelta; gdx < fdx+hwidth; gdx += delta) {
							d2 = gdx * gdx + gdy * gdy + gdz * gdz;
							factor = 1.0 - d2/r2;
							if (factor <= 0.0) continue;
							factor = 1.0;
							d1 = sqrt (d2);
							if (d1 <= 0.0) continue;
							c = gdx * base[0] + gdy * base[1] + gdz * base[2];
							c /= d1;
							s = gdx * zenith[0] + gdy * zenith[1] + gdz * zenith[2];
							s /= d1;
							angle = atan2 (s, c);
							isector = floor (nsector * (angle+PI) / (2 * PI));
							if (isector < 0) isector = 0;
							if (isector >= nsector) isector = nsector - 1;
							sector_numerator[isector] += factor * den;
							sector_denominator[isector] += factor;
						}
			}
		}
	}
	for (i = 0; i < org -> nsector; i++) {
		if (sector_denominator[i] > 0.0)
			org -> sector_density[i] = sector_numerator[i]/sector_denominator[i];
	}
	return (1);
}

float *set_orange_values (struct oringe *org, struct surface *phn) {
	int  k;
	long n_vertex, n_edge, n_triangle;
	long isector, nsector;
	long xx, yy, zz, iradius, v;
	float *values, *val;
	double r2, d2, factor, d1;
	double c, s;
	double fx, fy, fz, fdx, fdy, fdz;
	double sphere_radius, angle;
	double *vertices;
	double *vtx;
	double axis[3], base[3], zenith[3];
	double center[3];
	char message[MAXLINE];

	n_vertex = phn -> n_phnvtx;
	n_edge = phn -> n_phnedg;
	n_triangle = phn -> n_phntri;
	nsector = org -> nsector;
	vertices = phn -> vertex_centers;
	if (vertices == NULL) {
		set_error1 ("set_orange_values: null polyhedron vertex_centers");
		return (NULL);
	}
	values = allocate_floats (n_vertex);
	if (values == NULL) {
		set_error1 ("set_orange_values: not enough memory");
		return (NULL);
	}
	for (k = 0; k < 3; k++) {
		center[k] = org -> center[k];
		base[k] = org -> frame[0][k];
		zenith[k] = org -> frame[1][k];
		axis[k] = org -> frame[2][k];
	}
	sphere_radius = org -> radius;
	r2 = (sphere_radius * sphere_radius);
	for (v = 0; v < n_vertex; v++) {
		vtx = vertices + 3 * v;
		val = values + v;
		*val = 500.0;
		fx = *(vtx + 0);
		fy = *(vtx + 1);
		fz = *(vtx + 2);
		fdx = fx - center[0];
		fdy = fy - center[1];
		fdz = fz - center[2];
		d2 = fdx * fdx + fdy * fdy + fdz * fdz;
		factor = 1.0 - d2/r2;
		if (factor <= 0.0) continue;
		factor = 1.0;
		d1 = sqrt (d2);
		if (d1 <= 0.0) continue;
		c = fdx * base[0] + fdy * base[1] + fdz * base[2];
		c /= d1;
		s = fdx * zenith[0] + fdy * zenith[1] + fdz * zenith[2];
		s /= d1;
		angle = atan2 (s, c);
		isector = floor (nsector * (angle+PI) / (2 * PI));
		if (isector < 0) isector = 0;
		if (isector >= nsector) isector = nsector - 1;
		*val = org -> sector_density[isector];
	}
	return (values);
}

void store_orange_values (float *values, struct oringe *org, struct surface *phn)
{
	long v, nv;
	struct phnvtx *pv;
	double average_value;
	char message[MAXLINE];

	average_value = 0.0;
	nv = phn -> n_phnvtx;
	if (nv <= 0) {
		set_error1 ("store_orange_values: no vertices in polyhedron");
		return;
	}
	for (v = 0; v < nv; v++) {
		pv = num2phnvtx (phn, v + 1);
		if (pv == NULL) {
			set_error1 ("store_orange_values: invalid vertex number");
			return;
		}
		pv -> values[0] = *(values + v);
		pv -> values[1] = 0.0;
		pv -> values[2] = 0.0;
		average_value += pv -> values[0];
	}
	average_value /= nv;
	sprintf (message, "%8.3f average polyhedron vertex value", average_value);
	inform (message);
}

void rotate_polyhedron (struct oringe *org, struct surface *phn)
{
	int n_side, j, k;
	long i, t, v;
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
		center[k] = org -> center[k];
		translate[k] = 0.0;
	}
	for (j = 0; j < 3; j++)
		for (k = 0; k < 3; k++) {
			rotation[j][k] = org -> frame[j][k];
	}
	/* polyhedron of trianglular faces */
	/* transform vertices */
	for (v = 0; v < phn -> n_phnvtx; v++) {
		pvtx = num2phnvtx (phn, v + 1);
		if (pvtx == NULL) return;
		trpnt (rotation, center, translate, pvtx -> center);
		trvec (rotation, pvtx -> outward);
	}
	/* transform polyhedron triangles */
	if (phn -> phntri_handles != NULL) {
		for (t = 0; t < phn -> n_phntri; t++) {
			tri = num2phntri (phn, t+1);
			if (tri == NULL) {
				set_error1 ("rotate_polyhedron: null triangle pointer");
				return;
			}
			trpnt (rotation, center, translate, tri -> center);
			trvec (rotation, tri -> axis);
		}
	}
}

int allocate_standard (struct oringe *rng, long nsector)
{
	int s;
	struct sector *sec;
	for (s = 0; s < nsector; s++) {
		sec = allocate_sector ();
		if (sec == NULL) {
			set_error1 ("allocate_standard: ran out of memory");
			return (0);
		}
		rng -> sectors[s] = sec;
	}
	rng -> nsector = nsector;
	return (1);
}

int allocate_fourier (struct oringe *rng, long nfrequency)
{
	int t;
	struct term *trm;
	for (t = 0; t < nfrequency; t++) {
		trm = allocate_term ();
		if (trm == NULL) {
			set_error1 ("allocate_fourier: ran out of memory");
			return (0);
		}
		rng -> terms[t] = trm;
	}
	rng -> nfrequency = nfrequency;
	return (1);
}

int complement_standard (struct oringe *iring)
{
	int k;
	long i;
	long nsector;
	double temp[3];
	char message[MAXLINE];
	struct sector *sori;

	nsector = iring -> nsector;
	for (i = 0; i < nsector; i++) {
		sori = iring -> sectors[i];
		sori -> value = 1.0 - sori -> value;
	}
	/* complement concavity */
	iring -> concavity = 1.0 - iring -> concavity;
	/* do not reverse normal vector */
	/* for (k = 0; k < 3; k++) iring -> frame[2][k] = (-1.0) * iring -> frame[2][k]; */
	/* do rotate base and zenith vector 90 degrees */
	for (k = 0; k < 3; k++)
		temp[k] = iring -> frame[0][k];
	for (k = 0; k < 3; k++)
		iring -> frame[0][k] = iring -> frame[1][k];
	for (k = 0; k < 3; k++)
		iring -> frame[1][k] = (-temp[k]);
	return (1);
}


int standard_to_fourier (struct oringe *iring, double offset)
{
	long i, j, k, idx;
	long nsector, nfrequency;
	double value, angle, cos_angle, sin_angle;
	double cos_coefficient;
	double sin_coefficient;
	double modulus;
	double phase;
	char message[MAXLINE];
	struct term *trm;
	struct sector *sori;

	nsector = iring -> nsector;
	nfrequency = iring -> nfrequency;
	for (j = 0; j < nfrequency; j++) {
		trm = iring -> terms[j];
		trm -> order = j;
		cos_coefficient = 0.0;
		sin_coefficient = 0.0;
		for (i = 0; i < nsector; i++) {
			sori = iring -> sectors[i];
			angle = sori -> angle + offset;
			sin_angle = sin (j * angle);
			cos_angle = cos (j * angle);
			cos_coefficient += cos_angle * sori -> value * sori -> width;
			sin_coefficient += sin_angle * sori -> value * sori -> width;
		}
		cos_coefficient /= (2 * PI);
		sin_coefficient /= (2 * PI);
		modulus = cos_coefficient * cos_coefficient + sin_coefficient * sin_coefficient;
		if (modulus <= 0.0) modulus = 0.0;
		else modulus = sqrt (modulus);
		phase = atan2 (sin_coefficient, cos_coefficient);
		trm -> cos_coefficient = cos_coefficient;
		trm -> sin_coefficient = sin_coefficient;
		trm -> modulus = modulus;
		trm -> phase = phase;
	}
	return (1);
}

int fourier_to_standard (struct oringe *iring)
{
	long i, j, k, idx;
	long nfrequency, nsector;
	double value, angle, cos_angle, sin_angle, width;
	char message[MAXLINE];
	struct term *trm;
	struct sector *sor;

	nfrequency = iring -> nfrequency;
	nsector = iring -> nsector;
	for (j = 0; j < nsector; j++) {
		sor = iring -> sectors[j];
		sor -> value = 0.0;
		sor -> angle = (j + 0.5) * (2 * PI / nsector);
		sor -> width = 2 * PI / nsector;
	}
	for (j = 0; j < nsector; j++) {
		sor = iring -> sectors[j];
		angle = sor -> angle;
		width = sor -> width;
		for (i = 0; i < nfrequency; i++) {
			sin_angle = sin (i * angle);
			cos_angle = cos (i * angle);
			trm = iring -> terms[i];
			sor -> value += trm -> cos_coefficient * cos_angle;
			sor -> value += trm -> sin_coefficient * sin_angle;
		}
	}
	return (1);
}

int standard_to_abstract (struct oringe *iring, long nsample, long nbin1, char *partition)
{
	long nsector, nfrequency, i, n;
	long digitized[DEFAULT_FREQUENCY];
	double offset;
	double fraction;
	unsigned long packed;

	iring -> nbin1 = nbin1;
	nsector = iring -> nsector;
	nfrequency = iring -> nfrequency;
	fraction = 1.0 / nsample;
	for (i = 0; i < nsample; i++) {
		offset = i * (2 * PI / nsample);
		standard_to_fourier (iring, offset);
		if (error ()) return (0);
		digitize_fourier (nfrequency, iring -> terms, nbin1, digitized, partition);
		if (error ()) return (0);
		packed = pack_one (nfrequency, digitized);
		if (error ()) return (0);
		n = insert_bar (iring, packed, fraction);
		if (error ()) return (0);
		iring -> nbar += n;
	}
	iring -> nsample = nsample;
	return (1);
}

int digitize_fourier (long nfrequency, struct term *terms[], long nbin1, long digitized[], char *partition)
{
	long nbin2, ibin1, ibin2, idx, one;
	double c, s;
	struct term *trm;
	
	for (idx = 0; idx < nfrequency; idx++) {
		trm = terms[idx];
		c = trm -> cos_coefficient;
		s = trm -> sin_coefficient;
		one = digitize_one (c, s, nbin1, idx, partition);
		digitized[idx] = one;
	}
	/* kludge to ignore low-order terms */
	digitized[0] = 0;
	digitized[1] = 0;
	return (1);
}

long digitize_one (double c, double s, long nbin1, long idx, char *partition)
{
	long is, ic, ir, ncol, nrow, nannuli, nsec;
	long one, nbin, nbin2, ibin1, ibin2;
	double angle, modulus, scaled, fs, fc, offset;

	nbin = nbin1 * nbin1;
	one = 0;
	if (idx == 0) {
		modulus = (c * c + s * s);
		modulus = sqrt (modulus);
		scaled = ((modulus + 1.0)/2.0) * nbin;
		one = floor (scaled);
		if (one >= nbin) one = nbin - 1;
		else if (one < 0) one = 0;
	}
	else {
		if (strcmp(partition, "old") == 0) {
			modulus = (c * c + s * s);
			if (modulus < 0.0) modulus = 0.0;
			else if (modulus > 1.0) modulus = 1.0;
			modulus = sqrt (modulus);
			modulus = sqrt (modulus);
			modulus = sqrt (modulus);
			ibin1 = modulus * nbin1;
			if (ibin1 >= nbin1) ibin1 = nbin1 - 1;
			else if (ibin1 < 0) ibin1 = 0;
			nbin2 = 2 * ibin1 + 1;
			angle = PI + atan2 (s, c);
			ibin2 = nbin2 * (angle / ( 2 * PI));
			if (ibin2 >= nbin2) ibin2 = nbin2 - 1;
			else if (ibin2 < 0) ibin2 = 0;
			one = ibin1 * ibin1 + ibin2;
		}
		else if (strcmp(partition, "maze") == 0) {
			nannuli = sqrt ((double) nbin);
			modulus = (c * c + s * s);
			modulus = sqrt (modulus);
			modulus = sqrt (modulus);
			ir = (modulus * nannuli);
			if (ir >= nannuli) ir = nannuli - 1;
			else if (ir < 0) ir = 0;
			nsec = 2 * ir + 1;
			angle = atan2 (s, c);
			if (angle < 0.0) angle += 2.0 * PI;
			is = nsec * (angle / ( 2.0 * PI));
			if (is >= nsec) is = nsec - 1;
			else if (is < 0) is = 0;
			one = ir * ir + is;
		}
		else if (strcmp(partition, "dart") == 0) {
			nannuli = (long) sqrt ((double) nbin);
			nsec = nannuli;
			modulus = (c * c + s * s);
			modulus = sqrt (modulus);
			ir = (long) floor ((modulus * nannuli));
			if (ir >= nannuli) ir = nannuli - 1;
			else if (ir < 0) ir = 0;
			angle = atan2 (s, c);
			if (angle < 0.0) angle += 2.0 * PI;
			is = (long)  floor ((double) nsec * (angle / ( 2.0 * PI)));
			if (is >= nsec) is = nsec - 1;
			else if (is < 0) is = 0;
			one = ir * nsec + is;
		}
		else if (strcmp(partition, "stagger") == 0) {
			nannuli = (long) sqrt ((double) nbin);
			nsec = nannuli;
			modulus = (c * c + s * s);
			modulus = sqrt (modulus);
			angle = atan2 (s, c);
			if (angle < 0.0) angle += 2.0 * PI;
			is = (long)  floor ((double) nsec * (angle / ( 2.0 * PI)));
			if (is >= nsec) is = nsec - 1;
			else if (is < 0) is = 0;
			if (is % 2) ir = (long) floor (modulus * nannuli);
			else ir = (long) floor (modulus * nannuli + 0.5);
			if (ir >= nannuli) ir = nannuli - 1;
			else if (ir < 0) ir = 0;
			one = ir * nsec + is;
		}
		else if (strcmp(partition, "spiral") == 0) {
			nannuli = (long) sqrt ((double) nbin);
			nsec = nannuli;
			modulus = (c * c + s * s);
			modulus = sqrt (modulus);
			angle = atan2 (s, c);
			if (angle < 0.0) angle += 2.0 * PI;
			is = (long)  floor ((double) nsec * (angle / ( 2.0 * PI)));
			if (is >= nsec) is = nsec - 1;
			else if (is < 0) is = 0;
			offset = ((double) is) / ((double) nsec);
			ir = (long) floor (modulus * nannuli + offset);
			if (ir >= nannuli) ir = nannuli - 1;
			else if (ir < 0) ir = 0;
			one = ir * nsec + is;
		}
		else if (strcmp(partition, "chess") == 0) {
			nrow = (long) floor(sqrt ((double) nbin));
			ncol = (long) floor(sqrt ((double) nbin));
			modulus = (c * c + s * s);
			modulus = sqrt (modulus);
			/* assume inside radius 0.5 circle */
			/* fs = ((s + 1.0)/2.0) * nrow; */
			fs = 4.0 * (s + 0.125) * nrow;
			is = (long) floor (fs);
			/* fc = ((c + 1.0)/2.0) * ncol; */
			fc = 4.0 * (c + 0.125) * ncol;
			ic = (long) floor (fc);
			if (is >= nrow) is = nrow - 1;
			else if (is < 0) is = 0;
			if (ic >= ncol) ic = ncol - 1;
			else if (ic < 0) ic = 0;
			one = is * ncol + ic;
		}
		else return (0);
	}
	return (one);
}

unsigned long pack_one (long nfrequency, long digitized[])
{
	long i;
	unsigned long packed;
	packed = 0;
	for (i = 2; i <= 5; i++) {
		packed += digitized[i];
		if (i < 5) packed *= 256;
	}
	return (packed);
}

void unpack_one (long nfrequency, unsigned long packed, long digitized[])
{
	long i;
	for (i = 0; i < 8; i++)
		digitized[i] = 0;
	for (i = 5; i >= 2; i--) {
		digitized[i] = packed & 0xFF;
		packed /= 256;
	}
}

/* return number of new bars added */
int insert_bar (struct oringe *iring, unsigned long packed, double fraction)
{
	struct bar *b, *p, *n;

	if (iring -> nbar == 0) {
		n = allocate_bar ();
		if (n == NULL) {
			set_error1 ("insert_bar: memory failure");
			return (0);
		}
		iring -> head_bar = n;
		iring -> tail_bar = n;
		n -> fraction = fraction;
		n -> shape = packed;
		return (1);
	}
	p = NULL;
	for (b = iring -> head_bar; b != NULL; b = b -> next) {
		if (b -> shape == packed) {
			/* add to pre-existing bar */
			b -> fraction += fraction;
			return (0);
		}
		else if (packed > b -> shape) {
			/* keep going */
			p = b;
			continue;
		}
		else {
			/* new bar */
			n = allocate_bar ();
			if (n == NULL) {
				set_error1 ("insert_bar: memory failure");
				return (0);
			}
			if (p == NULL) {
				/* insert at beginning */
				n -> next = iring -> head_bar;
				iring -> head_bar = n;
				n -> fraction = fraction;
				n -> shape = packed;
				return (1);
			}
			/* insert in middle */
			n -> next = b;
			p -> next = n;
			n -> fraction = fraction;
			n -> shape = packed;
			return (1);
		}
	}
	/* not found - insert at end */
	n = allocate_bar ();
	if (n == NULL) {
		set_error1 ("insert_bar: memory failure");
		return (0);
	}
	iring -> tail_bar -> next = n;
	iring -> tail_bar = n;
	n -> fraction = fraction;
	n -> shape = packed;
	return (1);
}

/* oringes */

struct surface *polyhedron_oringes (struct oringe *head, double org_length, int srf_type)
{
	int j, k;
	long v, e, nsector, noringe, nvtx, nedg;
	double angle, c, s, d0, d1, d2;
	char message[MAXLINE];
	struct sector *sor;
	struct surface *obj;
	struct phnvtx *vtx, *vtx1, *vtx2, *org_vtx;
	struct phnvtx **org_vertices;
	struct phnedg *org_edg;
	struct phnedg **org_edges;
	struct oringe *org;

	noringe = 0; nsector = 0;
	for (org = head; org != NULL; org = org -> next) {
		if (nsector == 0) nsector = org -> nsector;
		else if (nsector != org -> nsector) {
			set_error1 ("(polyhedron_oringes): inconsistent number of sectors");
			return(NULL);
		}
		noringe++;
	}
	nvtx = noringe * nsector;
	nedg = nvtx; 

	/* allocate oringes bunch */
	obj = new_surface ();
	if (obj == NULL) {
		set_error1 ("(polyhedron_oringes): mem alloc failure");
		return(NULL);
	}
	obj -> type = srf_type;
	
	/* allocate memory for oringe vertices */
	org_vertices = (struct phnvtx **)
		allocate_pointers (PHNVTX, nvtx);
	if (org_vertices == NULL) {
		set_error1 ("(polyhedron_oringes): memory full");
		return(NULL);
	}
	
	/* allocate memory for oringe edges */
	org_edges = (struct phnedg **)
		allocate_pointers (PHNEDG, nedg);
	if (org_edges == NULL) {
		set_error1 ("(polyhedron_oringes): memory full");
		return(NULL);
	}
	
	/* create oringe vertices and edges */
	for (v = 0, org = head; org != NULL; v++, org = org -> next) {
		for (j = 0; j < nsector; j++) {
			sor = org -> sectors[j];
			angle = ((j + 0.5) / nsector) * 2.0 * PI;
			c = cos (angle);
			s = sin (angle);
			vtx1 = allocate_phnvtx ();
			if (vtx1 ==  NULL) {
				set_error1 ("(polyhedron_oringes): vertex allocation failure");
				return ((struct surface *) NULL);
			}
			*(org_vertices + nsector * v + j) = vtx1;
			for (k = 0; k < 3; k++) {
				d0 = c * org_length * org -> frame[0][k];
				d1 = s * org_length * org -> frame[1][k];
				d2 = (sor -> value - 0.5) * org_length * org -> frame[2][k];
				vtx1 -> center[k] = org -> center[k] + d0 + d1 + d2;
			}
			vtx1 -> number = nsector * v + j + 1;
		}
		for (j = 0; j < nsector; j++) {
			vtx1 = *(org_vertices + nsector * v + j);
			if (j < nsector - 1) vtx2 = *(org_vertices + nsector * v + j + 1);
			else vtx2 = *(org_vertices + nsector * v);
			org_edg = allocate_phnedg ();
			if (org_edg ==  NULL) {
				set_error1 ("(polyhedron_oringes): edge allocation failure");
				return ((struct surface *) NULL);
			}
			*(org_edges + nsector * v + j) = org_edg;
			org_edg -> pvt[0] = vtx1;
			org_edg -> pvt[1] = vtx2;
			org_edg -> vtxnum[0] = vtx1 -> number;
			org_edg -> vtxnum[1] = vtx2 -> number;
		}
	}
	/* set up linked list */
	obj -> head_phnvtx = *org_vertices;
	for (v = 0; v < nvtx - 1; v++) {
		org_vtx = *(org_vertices + v);
		org_vtx -> next = *(org_vertices + v + 1);
	}
	obj -> head_phnedg = *org_edges;
	for (e = 0; e < nedg - 1; e++) {
		org_edg = *(org_edges + e);
		org_edg -> next = *(org_edges + e + 1);
	}
		
	sprintf (message,"%8.3f oringe length", org_length);
	inform(message);
	sprintf (message,"%8ld oringe vectors", nedg);
	inform(message);
	/* transfer to structure variables */
	obj -> n_polygon = 0L;
	obj -> n_phnvtx = nvtx;
	obj -> n_phnedg = nedg;
	obj -> head_phnvtx = *org_vertices;
	obj -> head_phnedg = *org_edges;
	obj -> head_polygon = (struct polygon *) NULL;
	
	obj -> phnvtx_handles = org_vertices;
	obj -> phnedg_handles = org_edges;

	return (obj);
}


struct polygon *read_polygon (FILE *fp)
{
	int nscan, j, k, nsector;
	long count, i, n, v, nv, is, ns;
	double x, y, z, r, fraction, factor;
	double vx, vy;
	char message[MAXLINE];
	char input_line[MAXLINE+1];
	struct polygon *poly;
	struct record *polygon_record;
	struct record *vertex_record;
	struct token *tok;

	polygon_record = get_record (POLYGON_RECORD, fp);
	if (error ()) return (NULL);
	if (polygon_record -> nfield != 1) {
		sprintf(message, "(read_polygon): only %d fields for polygon record, should be %d",
			polygon_record -> nfield, 1);
		set_error1 (message);
		return (NULL);
	}
	if (polygon_record -> type != POLYGON_RECORD) {
		set_error1("(read_polygon): wrong type for polygon record");
		return (NULL);
	}
	poly = allocate_polygon ();
	if (error()) return (NULL);
	if (poly == NULL) {
		set_error2("(read_polygon): no memory");
		return (NULL);
	}
	if (polygon_record -> field_type[2] != VERTEX_RECORD) {
		set_error1("(read_polygon): bad format for polygon record");
		return (NULL);
	}
	nv = polygon_record -> field_count[2];
	for (v = 0; v < nv; v++) {
		if (v == 0) vertex_record = polygon_record -> head[2];
		else vertex_record = vertex_record -> next;
		if (vertex_record == NULL) {
			set_error1("(read_polygon): missing vertex for polygon record");
			return (NULL);
		}
		tok = vertex_record -> first_token;
		vx = tok -> fvalue;
		tok = tok -> next;
		vy = tok -> fvalue;
		poly -> vc[v][0] = vx;
		poly -> vc[v][1] = vy;
		poly -> vc[v][2] = 0.0;
	}
	poly -> n_side = nv;
	deep_free_record (polygon_record);
	return (poly);
}

struct oringe *read_oringe (FILE *fp)
{
	int nscan, j, k, nsector, nfrequency;
	int f;
	long count, i, n, v, nv, is, ns;
	unsigned long packed;
	long digitized[DEFAULT_FREQUENCY];
	double x, y, z, r, fraction, factor;
	double vx, vy;
	double value, width, angle, c, s, delta;
	double rx0, ry0, rz0, rx1, ry1, rz1, rx2, ry2, rz2;
	char message[MAXLINE];
	char input_line[MAXLINE+1];
	struct oringe *iring;
	struct sector *sor;
	struct term *trm;
	struct record *oringe_record, *center_record, *radius_record, *nsector_record;
	struct record *vertex_record, *sector_record, *frequency_record;
	struct record *abstract_record, *rotation_record, *concavity_record;
	struct token *tok;

	oringe_record = get_record (ORINGE_RECORD, fp);
	if (error ()) return (NULL);
	if (oringe_record == NULL) return (NULL);
	if (oringe_record -> type != ORINGE_RECORD) {
		set_error1("(read_oringe): wrong type");
		return (NULL);
	}

	iring = allocate_oringe ();
	if (error()) return (NULL);
	if (iring == NULL) {
		set_error2("(read_oringe): no memory");
		return (NULL);
	}
	iring -> nsector = 0;
	iring -> center[0] = 0.0;
	iring -> center[1] = 0.0;
	iring -> center[2] = 0.0;
	iring -> radius = 0.0;
	for (j = 0; j < 3; j++)
		for (k = 0; k < 3; k++)
			iring -> frame[j][k] = (float) (j == k);
	iring -> concavity = 0.0;
	for (i = 0; i < iring -> nsector; i++)
		iring -> sector_density[i] = 0.0;
	angle = 0.0;
	width = 0.0;
	value = 0.0;

	f = 0;
	if (f >= oringe_record -> nfield) {
		deep_free_record (oringe_record);
		return (iring);
	}
	if (oringe_record -> field_type[f] == CENTER_RECORD) {
		center_record = oringe_record -> head[f];
		if (center_record == NULL) {
			set_error1("(read_oringe): missing center for standard oringe record");
			return (NULL);
		}
		tok = center_record -> first_token;
		iring -> center[0] = tok -> fvalue;
		tok = tok -> next;
		iring -> center[1] = tok -> fvalue;
		tok = tok -> next;
		iring -> center[2] = tok -> fvalue;
		f++;
	}
	if (f >= oringe_record -> nfield) {
		deep_free_record (oringe_record);
		return (iring);
	}
	if (oringe_record -> field_type[f] == RADIUS_RECORD) {
		radius_record = oringe_record -> head[f];
		if (radius_record == NULL) {
			set_error1("(read_oringe): missing radius for standard oringe record");
			return (NULL);
		}
		tok = radius_record -> first_token;
		iring -> radius = tok -> fvalue;
		f++;
	}
	if (f >= oringe_record -> nfield) {
		deep_free_record (oringe_record);
		return (iring);
	}
	if (oringe_record -> field_type[f] == NSECTOR_RECORD) {
		nsector_record = oringe_record -> head[f];
		if (nsector_record == NULL) {
			set_error1("(read_oringe): missing nsector simple for oringe record");
			return (NULL);
		}
		tok = nsector_record -> first_token;
		nsector = tok -> ivalue;
		if (nsector > MAX_SECTOR) nsector = MAX_SECTOR;
		if (nsector < 1) nsector = 1;
		iring -> nsector = nsector;
		f++;
	}
	if (f >= oringe_record -> nfield) {
		deep_free_record (oringe_record);
		return (iring);
	}
	if (oringe_record -> field_type[f] == ROTATION_RECORD) {
		rotation_record = oringe_record -> head[f];
		if (rotation_record == NULL) {
			set_error1("(read_oringe): missing rotation for standard oringe record");
			return (NULL);
		}
		tok = rotation_record -> first_token;
		iring -> frame[0][0] = tok -> fvalue;
		tok = tok -> next;
		iring -> frame[0][1] = tok -> fvalue;
		tok = tok -> next;
		iring -> frame[0][2] = tok -> fvalue;
		tok = tok -> next;
		iring -> frame[1][0] = tok -> fvalue;
		tok = tok -> next;
		iring -> frame[1][1] = tok -> fvalue;
		tok = tok -> next;
		iring -> frame[1][2] = tok -> fvalue;
		tok = tok -> next;
		iring -> frame[2][0] = tok -> fvalue;
		tok = tok -> next;
		iring -> frame[2][1] = tok -> fvalue;
		tok = tok -> next;
		iring -> frame[2][2] = tok -> fvalue;
		f++;
	}
	if (f >= oringe_record -> nfield) {
		deep_free_record (oringe_record);
		return (iring);
	}
	if (oringe_record -> field_type[f] == CONCAVITY_RECORD) {
		concavity_record = oringe_record -> head[f];
		if (concavity_record == NULL) {
			set_error1("(read_oringe): missing concavity for standard oringe record");
			return (NULL);
		}
		tok = concavity_record -> first_token;
		iring -> concavity = tok -> fvalue;
		f++;
	}
	if (f >= oringe_record -> nfield) {
		deep_free_record (oringe_record);
		return (iring);
	}
	if (oringe_record -> field_type[f] == SECTOR_RECORD) {
		ns = oringe_record -> field_count[f];
		allocate_standard (iring, ns);
		if (error ()) return (NULL);
		angle = 0.0;
		/* read values, two per line */
		for (is = 0; is < ns; is++) {
			if (is == 0) sector_record = oringe_record -> head[f];
			else sector_record = sector_record -> next;
			if (sector_record == NULL) {
				set_error1("(read_oringe): missing sector for standard oringe record");
				return (NULL);
			}
			tok = sector_record -> first_token;
			if (tok -> type != FLOAT_FIELD && tok -> type != INT_FIELD) {
				sprintf(message, "(read_oringe): missing float field for sector");
				set_error1 (message);
				return (NULL);
			}
			width = tok -> fvalue;
			tok = tok -> next;
			if (tok -> type != FLOAT_FIELD && tok -> type != INT_FIELD) {
				sprintf(message, "(read_oringe): missing float field for sector");
				set_error1 (message);
				return (NULL);
			}
			value = tok -> fvalue;
			sor = iring -> sectors[is];
			angle += width/2;
			sor -> angle = angle;
			angle += width/2;
			sor -> width = width;
			sor -> value = value;
		}
		factor = (2 * PI) / angle;
		for (is = 0; is < ns; is++) {
			sor = iring -> sectors[is];
			sor -> angle *= factor;
			sor -> width *= factor;
		}
		iring -> nsector = ns;
		f++;
	}
	if (f >= oringe_record -> nfield) {
		deep_free_record (oringe_record);
		return (iring);
	}
	if (oringe_record -> field_type[f] == FREQUENCY_RECORD) {
		ns = oringe_record -> field_count[f];
		allocate_fourier (iring, (long) ns);
		if (error ()) return (NULL);
		/* read values, two per line */
		for (is = 0; is < ns; is++) {
			if (is == 0) frequency_record = oringe_record -> head[f];
			else frequency_record = frequency_record -> next;
			if (frequency_record == NULL) {
				set_error1("(read_oringe): null frequency for fourier oringe record");
				return (NULL);
			}
			tok = frequency_record -> first_token;
			if (tok -> type != INT_FIELD) {
				sprintf(message, "(read_oringe): missing int field for frequency");
				set_error1 (message);
				return (NULL);
			}
			if (tok -> ivalue != is) {
				sprintf(message, "(read_oringe): frequency inconsistency, %ld != %ld",
					tok -> ivalue, is);
				set_error1 (message);
				return (NULL);
			}
			tok = tok -> next;
			if (tok -> type != FLOAT_FIELD && tok -> type != INT_FIELD) {
				sprintf(message, "(read_oringe): missing float field for frequency");
				set_error1 (message);
				return (NULL);
			}
			c = tok -> fvalue;
			tok = tok -> next;
			if (tok -> type != FLOAT_FIELD && tok -> type != INT_FIELD) {
				sprintf(message, "(read_oringe): missing float field for frequency");
				set_error1 (message);
				return (NULL);
			}
			s = tok -> fvalue;
			trm = iring -> terms[is];
			trm -> cos_coefficient = c;
			trm -> sin_coefficient = s;
		}
		iring -> nfrequency = ns;
		f++;
	}
	if (f >= oringe_record -> nfield) {
		deep_free_record (oringe_record);
		return (iring);
	}
	if (oringe_record -> field_type[f] == BAR_RECORD) {
		ns = oringe_record -> field_count[f];
		iring -> nbar = 0;
		/* read histogram */
		for (is = 0; is < ns; is++) {
			if (is == 0) abstract_record = oringe_record -> head[f];
			else abstract_record = abstract_record -> next;
			if (abstract_record == NULL) {
				set_error1("(read_oringe): null bar for abstract oringe record");
				return (NULL);
			}
			nfrequency = abstract_record -> ntoken - 1;
			tok = abstract_record -> first_token;
			if (tok -> type != FLOAT_FIELD && tok -> type != INT_FIELD) {
				sprintf(message, "(read_oringe): missing float field for bar");
				set_error1 (message);
				return (NULL);
			}
			fraction = tok -> fvalue;
			for (j = 0; j < nfrequency; j++) {
				tok = tok -> next;
				if (tok == NULL) break;
				if (tok -> type != INT_FIELD) {
					sprintf(message, "(read_oringe): missing int field for bar");
					set_error1 (message);
					return (NULL);
				}
				digitized[j] = tok -> ivalue;
			}
			packed = pack_one (nfrequency, digitized);
			if (error ()) return (0);
			n = insert_bar (iring, packed, fraction);
			if (error ()) return (0);
			iring -> nbar += n;
		}
	}
	deep_free_record (oringe_record);
	return (iring);
}

int write_oringe (struct oringe *oring, FILE *fp, char *fmt, char *partition)
{
	long i, j, nbin1, nfrequency, nsector, nbar;
	long unpacked[DEFAULT_FREQUENCY];
	int byte0, byte1, byte2, byte3;
	int byte4, byte5, byte6, byte7;
	struct sector *sor;
	struct term *trm;
	char message[MAXLINE];
	struct bar *b;

	nfrequency = oring -> nfrequency;
	nsector = oring -> nsector;
	nbin1 = oring -> nbin1;
	nbar = oring -> nbar;
	if (strcmp (fmt, "standard") == 0) {
		fprintf (fp, "oringe: 1 center, 1 radius, 1 rotation, 1 concavity, %ld sectors {\n", nsector);
		fprintf (fp, "center: 3 floats {  %8.3f, %8.3f, %8.3f }\n",
			oring -> center[0], oring -> center[1], oring -> center[2]);
		fprintf (fp, "radius: 1 float {  %7.3f } \n", oring -> radius);
		fprintf (fp,
		"rotation: 9 floats { %7.4f, %7.4f, %7.4f,  %7.4f, %7.4f, %7.4f,  %7.4f, %7.4f, %7.4f}\n",
			oring -> frame[0][0], oring -> frame[0][1], oring -> frame[0][2],
			oring -> frame[1][0], oring -> frame[1][1], oring -> frame[1][2],
			oring -> frame[2][0], oring -> frame[2][1], oring -> frame[2][2]);
		fprintf (fp, "concavity: 1 float { %7.3f }\n",  oring -> concavity);
		for (i = 0; i < nsector; i++) {
			sor = oring -> sectors[i];
			fprintf (fp, "sector: 2 floats { %12.6f %12.6f}\n",
				sor -> width, sor -> value);
		}
		fprintf (fp, "}\n");
	}
	else if (strcmp (fmt, "fourier") == 0) {
		fprintf (fp, "oringe: 1 center, 1 radius, 1 rotation, 1 concavity, %ld frequencies {\n", nfrequency);
		fprintf (fp, "center: 3 floats {  %8.3f, %8.3f, %8.3f }\n",
			oring -> center[0], oring -> center[1], oring -> center[2]);
		fprintf (fp, "radius: 1 float {  %7.3f } \n", oring -> radius);
		fprintf (fp,
		"rotation: 9 floats { %7.4f, %7.4f, %7.4f,  %7.4f, %7.4f, %7.4f,  %7.4f, %7.4f, %7.4f}\n",
			oring -> frame[0][0], oring -> frame[0][1], oring -> frame[0][2],
			oring -> frame[1][0], oring -> frame[1][1], oring -> frame[1][2],
			oring -> frame[2][0], oring -> frame[2][1], oring -> frame[2][2]);
		fprintf (fp, "concavity: 1 float { %7.3f }\n",  oring -> concavity);
		for (i = 0; i < nfrequency; i++) {
			trm = oring -> terms[i];
			fprintf (fp, "frequency: 1 int, 2 floats { %3ld %12.6f %12.6f}\n",
				i, trm -> cos_coefficient, trm -> sin_coefficient);
		}
		fprintf (fp, "}\n");
	}
	else if (strcmp (fmt, "abstract") == 0) {
		fprintf (fp, "oringe: 1 center, 1 radius, 1 rotation, 1 concavity, 1 partition, %d amplitudes, %ld bars {\n",
			nfrequency, nbar);
		fprintf (fp, "center: 3 floats {  %8.3f, %8.3f, %8.3f }\n",
			oring -> center[0], oring -> center[1], oring -> center[2]);
		fprintf (fp, "radius: 1 float {  %7.3f } \n", oring -> radius);
		fprintf (fp,
		"rotation: 9 floats { %7.4f, %7.4f, %7.4f,  %7.4f, %7.4f, %7.4f,  %7.4f, %7.4f, %7.4f}\n",
			oring -> frame[0][0], oring -> frame[0][1], oring -> frame[0][2],
			oring -> frame[1][0], oring -> frame[1][1], oring -> frame[1][2],
			oring -> frame[2][0], oring -> frame[2][1], oring -> frame[2][2]);
		fprintf (fp, "concavity: 1 float { %7.3f }\n",  oring -> concavity);
		fprintf (fp, "partition: 1 int, 1 enum, 1 int { %d, %s, %2ld }\n",
			nbin1 * nbin1, partition, nfrequency);
		for (i = 0; i < nfrequency; i++) {
			trm = oring -> terms[i];
			fprintf (fp, "amplitude: 1 int, 1 float { %3ld %12.6f}\n",
				i, trm -> modulus);
		}
		for (b = oring -> head_bar; b != NULL; b = b -> next) {
			unpack_one (nfrequency, b -> shape, unpacked);
			byte0 = unpacked[0];
			byte1 = unpacked[1];
			byte2 = unpacked[2];
			byte3 = unpacked[3];
			if (nfrequency >= 8) {
				byte4 = unpacked[4];
				byte5 = unpacked[5];
				byte6 = unpacked[6];
				byte7 = unpacked[7];
			}
			else {
				byte4 = 0;
				byte5 = 0;
				byte6 = 0;
				byte7 = 0;
			}
			fprintf (fp, "bar: 1 float, 8 ints { %7.4f %3d %3d %3d %3d %3d %3d %3d %3d}\n",
				b -> fraction, byte0, byte1, byte2, byte3, byte4, byte5, byte6, byte7);
		}
		fprintf (fp, "}\n");
	}
	else {
		set_error1("(write_oringe): unknown output format");
		return (0);
	}
	return (1);
}


int write_phnorg (struct surface *phn, double sphere_radius, long nsector, FILE *fp_org)
{
	int b, i;
	long v;
	float float_density[MAX_SECTOR];
	char message[MAXLINE];
	struct phnvtx *pv;

	if (phn -> phnvtx_handles == NULL) return (0);

	/* write vertex coordinates */
	for (v = 0; v < phn -> n_phnvtx; v++) {
		pv = num2phnvtx (phn, v + 1);
		if (pv == NULL) return (0);
		/* compress from float to unsigned char */
		for (i = 0; i < nsector; i++) {
			b = *(pv -> byte_density + i);
			float_density[i] = (-1.0) + b/128.0;
		}
		fprintf (fp_org, "oringe: 1 center, 1 radius, 1 rotation, 1 concavity, %ld sectors {\n", nsector);
		fprintf (fp_org, "center: 3 floats { %8.3f, %8.3f, %8.3f }\n",
			pv -> center[0], pv -> center[1], pv -> center[2]);
		fprintf (fp_org, "radius: 1 float { %7.3f }\n", sphere_radius);
		fprintf (fp_org,
		"rotation: 9 floats { %7.4f, %7.4f, %7.4f,  %7.4f, %7.4f, %7.4f,  %7.4f, %7.4f, %7.4f}\n",
			pv -> base[0], pv -> base[1], pv -> base[2],
			pv -> zenith[0], pv -> zenith[1], pv -> zenith[2],
			pv -> normal[0], pv -> normal[1], pv -> normal[2]);
		fprintf (fp_org, "concavity: 1 float { %7.3f }\n",  pv -> values[0]);
		for (i = 0; i < nsector; i++) {
			fprintf (fp_org, "sector: 2 floats { %12.6f %12.6f}\n",
				2 * PI / nsector, float_density[i]);
		}
		fprintf (fp_org, "}\n");
	}
	return (1);
}


struct oringe *allocate_oringe ()
{
	struct oringe *rng;

	/* allocate memory */
	rng = (struct oringe *) allocate_object (ORINGE);
	if (rng == NULL) {
		set_error1 ("allocate_oringe: ran out of memory");
		return(NULL);
	}
	return (rng);
}

void free_oringe (struct oringe *rng)
{
	free_object (ORINGE, (short *) rng);
	if (error()) return;
}

void deep_free_oringe (struct oringe *rng) 
{
	int s, t;
	struct sector *sec;
	struct term *trm;
	struct bar *b, *n;
	if (rng -> nsector > 0) {
		for (s = 0; s < rng -> nsector; s++) {
			sec = rng -> sectors[s];
			if (sec != NULL) free_sector (sec);
		}
	}
	if (rng -> nfrequency > 0) {
		for (t = 0; t < rng -> nfrequency; t++) {
			trm = rng -> terms[t];
			if (trm != NULL) free_term (trm);
		}
	}
	for (b = rng -> head_bar; b != NULL; b = n) {
		n = b -> next;
		free_bar (b);
	}
	free_oringe (rng);
}

struct sector *allocate_sector ()
{
	struct sector *sec;

	/* allocate memory */
	sec = (struct sector *) allocate_object (SECTOR);
	if (sec == NULL) {
		set_error1 ("allocate_sector: ran out of memory");
		return(NULL);
	}
	return (sec);
}

void free_sector (struct sector *sec)
{
	free_object (SECTOR, (short *) sec);
	if (error()) return;
}


struct term *allocate_term ()
{
	struct term *trm;

	/* allocate memory */
	trm = (struct term *) allocate_object (TERM);
	if (trm == NULL) {
		set_error1 ("allocate_term: ran out of memory");
		return(NULL);
	}
	return (trm);
}

void free_term (struct term *trm)
{
	free_object (TERM, (short *) trm);
	if (error()) return;
}


struct bar *allocate_bar ()
{
	struct bar *b;

	/* allocate memory */
	b = (struct bar *) allocate_object (BAR);
	if (b == NULL) {
		set_error1 ("allocate_bar: ran out of memory");
		return(NULL);
	}
	return (b);
}

void free_bar (struct bar *b)
{
	free_object (BAR, (short *) b);
	if (error()) return;
}


/* MSP: Copyright 1995 by Michael L. Connolly */
