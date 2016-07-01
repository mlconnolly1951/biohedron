/* Molecular Surface Package */
/* Copyright 1995 by Michael L. Connolly */
/* February 3, 2000 */

#include "ms.h"
#include "msenum.h"
#include "mslimit.h"
#include "msdefault.h"
#include "msstruct.h"
#include "msglobal.h"
#include "msfunction.h"


struct vanity *smooth_density (long dimension, long dimensions[3], float *densities, float sphere_radius, float cube_width, int smooth)
{
	long idx, jdx, nx, ny, nz, x, y, z, nden;
	long xx, yy, zz, dx, dy, dz, iradius;
	double r2, d2, factor;
	double den, numerator, denominator, den0;
	double xden, yden, zden, xnumerator, ynumerator, znumerator;
	double fdx, fdy, fdz, gdx, gdy, gdz;
	double delta, hdelta, hwidth, fr2;
	struct vanity *smoothies;
	
	nx = dimensions[0];
	ny = dimensions[1];
	nz = dimensions[2];
	nden = nx * ny * nz;
	if (nden <= 0) return (NULL);
	if (cube_width <= 0.0) return (NULL);
	hwidth = cube_width/2;
	delta = DELTA_WIDTH;
	if (hwidth < delta) delta = hwidth;
	hdelta = delta/2;
	smoothies = (struct vanity *) allocate_objects (VANITY, nden);
	if (smoothies == NULL) return (NULL);
	if (sphere_radius <= hwidth) {
		/* just transfer the data; no normals */
		for (idx = 0; idx < nden; idx++)
			(smoothies + idx) -> scalar = *(densities + idx);
		return (smoothies);
	}
	iradius = (sphere_radius / cube_width) + 1.5;
	r2 = (sphere_radius * sphere_radius) / (cube_width * cube_width);
	fr2 = (sphere_radius * sphere_radius);
	for (z = 0; z < nz; z++) {
		for (y = 0; y < ny; y++) {
			for (x = 0; x < nx; x++) {
				idx = z * ny * nx + y * nx + x;
				numerator = 0.0;
				xnumerator = 0.0;
				ynumerator = 0.0;
				znumerator = 0.0;
				denominator = 0.0;
				for (dz = -iradius; dz <= iradius; dz++) {
					zz = z + dz;
					if (zz < 0 || zz >= nz) continue;
					fdz = (dz+0.5) * cube_width;
					for (dy = -iradius; dy <= iradius; dy++) {
						yy = y + dy;
						if (yy < 0 || yy >= ny) continue;
						fdy = (dy+0.5) * cube_width;
						for (dx = -iradius; dx <= iradius; dx++) {
							xx = x + dx;
							if (xx < 0 || xx >= nx) continue;
							fdx = (dx+0.5) * cube_width;
							if (dimension == 2) {
								d2 = dx * dx + dy * dy + dz * dz;
								factor = 1.0  - d2/r2;
								if (factor <= 0.0) continue;
								idx = y * nx + x;
								jdx = yy * nx + xx;
								den0 = *(densities+idx);
								den = *(densities+jdx);
								numerator += factor * den;
								denominator += factor;
								xden = cube_width * dx * den;
								xnumerator += factor * xden;
								yden = cube_width * dy * den;
								ynumerator += factor * yden;
								zden = cube_width * dz * den;
								znumerator += factor * zden;
							}
							else if (dimension == 3) {
								idx = z * ny * nx + y * nx + x;
								jdx = zz * ny * nx + yy * nx + xx;
								den0 = *(densities+idx);
								den = *(densities+jdx);
								for (gdz = fdz-hwidth + hdelta; gdz < fdz+hwidth; gdz += delta)
									for (gdy = fdy-hwidth + hdelta; gdy < fdy+hwidth; gdy += delta)
										for (gdx = fdx-hwidth + hdelta; gdx < fdx+hwidth; gdx += delta) {
											d2 = gdx * gdx + gdy * gdy + gdz * gdz;
											factor = 1.0 - d2/fr2;
											if (factor <= 0.0) continue;
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
				}
				if (denominator <= 0.0) {
					(smoothies+idx) -> scalar = den0;
				}
				else {
					if (!smooth)
						(smoothies+idx) -> scalar = den0;
					else (smoothies+idx) -> scalar = (numerator/denominator);
					/* use negative of gradient to get outward-pointing vector */
					(smoothies+idx) -> vector[0] = -(xnumerator/denominator);
					(smoothies+idx) -> vector[1] = -(ynumerator/denominator);
					(smoothies+idx) -> vector[2] = -(znumerator/denominator);
				}
			}
		}
	}
	return (smoothies);
}

struct msdata *eval_density (long dimension, long dimensions[3], double *vertices, double *normals, double *bases, long n_vertex, float *densities, float sphere_radius, float cube_width, int do_fourier, long nsector)
{
	int dig;
	int  i, k, normal_exists;
	long isector;
	long jdx, nx, ny, nz, x, y, z, nden, npos;
	long xx, yy, zz, dx, dy, dz, iradius, v;
	unsigned char *sectbytes, *sb;
	float *values, *val;
	double r2, d2, factor, d1, delta, hdelta, hwidth;
	double gdx, gdy, gdz;
	double den, numerator, denominator;
	double xden, yden, zden;
	double xnumerator, ynumerator, znumerator;
	double c, s;
	double angle;
	double c2, s2, c2den, s2den, c2numerator, s2numerator, numerator2;
	double c3, s3, c3den, s3den, c3numerator, s3numerator, numerator3;
	double fx, fy, fz, fdx, fdy, fdz;
	double *vtx, *nml, *bas;
	double axis[3], base[3], zenith[3], newbase[3];
	unsigned char byte_density[MAX_SECTOR];
	double sector_density[MAX_SECTOR];
	double sector_numerator[MAX_SECTOR];
	double sector_denominator[MAX_SECTOR];
	char message[MAXLINE];
	struct msdata *msd;
	
	if (dimension != 3) {
		set_error1 ("eval_density: dimension must be 3");
		return (NULL);
	}
	nx = dimensions[0];
	ny = dimensions[1];
	nz = dimensions[2];
	nden = nx * ny * nz;
	if (nden <= 0) return (NULL);
	iradius = floor (sphere_radius / cube_width) + 2.5;
	r2 = (sphere_radius * sphere_radius);
	delta = DELTA_WIDTH;
	hdelta = delta/2;
	hwidth = cube_width/2;
	values = allocate_floats (n_vertex * 6);
	if (values == NULL) {
		set_error1 ("eval_density: not enough memory for values");
		return (NULL);
	}
	if (nsector > 0) {
		sectbytes = allocate_bytes (n_vertex * nsector);
		if (sectbytes == NULL) {
			set_error1 ("eval_density: not enough memory for sectbytes");
			return (NULL);
		}
	}
	else sectbytes = (unsigned char *) NULL;
	npos = 0;
	for (v = 0; v < n_vertex; v++) {
		vtx = vertices + 3 * v;
		nml = normals + 3 * v;
		bas = bases + 3 * v;
		val = values + 6 * v;
		if (nsector > 0) sb = sectbytes + nsector * v;
		else sb = NULL;
		fx = *(vtx + 0);
		fy = *(vtx + 1);
		fz = *(vtx + 2);
		x = floor (fx / cube_width);
		y = floor (fy / cube_width);
		z = floor (fz / cube_width);
		numerator = 0.0; denominator = 0.0;
		xnumerator = 0.0; ynumerator = 0.0; znumerator = 0.0;
		c2numerator = 0.0; s2numerator = 0.0;
		c3numerator = 0.0; s3numerator = 0.0;
		normal_exists = 1;
		for (i = 0; i < nsector; i++) {
			sector_density[i] = 0.0;
			sector_numerator[i] = 0.0;
			sector_denominator[i] = 0.0;
		}
		for (k = 0; k < 3; k++) {
			axis[k] = 0.0;
			base[k] = 0.0;
			zenith[k] = 0.0;
		}
		for (k = 0; k < 3; k++)
			axis[k] = *(nml + k);
		if (norm(axis) <= 0.0)
			normal_exists = 0;
		if (normal_exists) {
			if (nsector > 0) {
				/* use base vector from 2-fold call done earlier */
				for (k = 0; k < 3; k++)
					base[k] = *(bas + k);
			}
			else arbprp (axis, base);
			cross (axis, base, zenith);
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
								if (do_fourier && normal_exists) {
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
									/* two-theta fourier term */
									c2 = c * c - s * s;
									s2 = 2 * s * c;
									c2den = c2 * den;
									s2den = s2 * den;
									c2numerator += factor * c2den;
									s2numerator += factor * s2den;
									/* three-fold symmetry */
									c3 = c2 * c - s2 * s;
									s3 = s * c2 + s2 * c;
									c3den = c3 * den;
									s3den = s3 * den;
									c3numerator += factor * c3den;
									s3numerator += factor * s3den;
								}
							}
				}
			}
		}
		if (denominator > 0.001) {
			if (numerator > 0.0) npos++;
			/* use negative of gradient to get outward-pointing vector */
			*(val+0) = (numerator/denominator);
			*(val+1) = 0.0;
			*(val+2) = 0.0;
			*(val+3) = -(xnumerator/denominator);
			*(val+4) = -(ynumerator/denominator);
			*(val+5) = -(znumerator/denominator);
			if (do_fourier) {
				numerator2 = sqrt (c2numerator * c2numerator + s2numerator * s2numerator);
				numerator3 = sqrt (c3numerator * c3numerator + s3numerator * s3numerator);
				*(val+1) = (numerator2/denominator);
				*(val+2) = (numerator3/denominator);
				if (nsector == 0) {
					angle = atan2 (s2numerator, c2numerator);
					c = cos (angle/2.0);
					s = sin (angle/2.0);
					for (k = 0; k < 3; k++)
						newbase[k] = c * base[k] + s * zenith[k];
					*(val+3) = newbase[0];
					*(val+4) = newbase[1];
					*(val+5) = newbase[2];
				}
			}
		}
		if (nsector > 0) {
			for (i = 0; i < nsector; i++) {
				if (sector_denominator[i] > 0.001)
					sector_density[i] = sector_numerator[i]/sector_denominator[i];
				else sector_density[i] = 0.0;
			}
			/* compress from float to unsigned char */
			for (i = 0; i < nsector; i++) {
				dig = 128 * (1.0 + sector_density[i]);
				if (dig < 0) dig = 0;
				if (dig > 255) dig = 255;
				byte_density[i] = (unsigned char) dig;
			}
			for (i = 0; i < nsector; i++) {
				*(sb + i) = byte_density[i];
			}
		}
	}
	sprintf (message, "%8ld positive numerators, %8.3f delta", npos, delta);
	inform (message);
	/* transfer data to msdata */
	msd = newmsdata ();
	if (msd == NULL) {
		set_error1 ("eval_density: no memory for msd");
		return (NULL);
	}
	if (!addfloats (msd, 6L, n_vertex, values)) {
		set_error1 ("eval_density: addfloats fails for vertex values");
		return (NULL);
	}
	if (nsector > 0) {
		if (!addbytes (msd, (long) nsector, n_vertex, sectbytes)) {
			set_error1 ("eval_density: addbytes fails");
		}
	}
	return (msd);
}

/* Molecular Surface Package */
/* Copyright 1995 by Michael L. Connolly */
