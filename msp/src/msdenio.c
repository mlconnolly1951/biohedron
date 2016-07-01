#include "ms.h"
#include "msenum.h"
#include "mslimit.h"
#include "msdefault.h"
#include "msstruct.h"
#include "msglobal.h"
#include "msfunction.h"

/* msdensity.c Copyright 1995 by Michael L. Connolly */

/* February 3, 2000 */

struct surface *create_density (double cube_width, int origin[3], int width[3])
{
	int k;
	long n_cube;
	float *densities;
	struct surface *den;
	
	den = new_surface ();
	if (den == NULL) return (NULL);
	den -> type = DEN_SURFACE;
	n_cube = width[0] * width[1] * width[2];
	densities = allocate_floats (n_cube);
	if (densities == NULL) return (NULL);
	den -> densities = densities;
	den -> cube_width = cube_width;
	den -> n_cube = n_cube;
	for (k = 0; k < 3; k++) {
		den -> origin[k] = origin[k];
		den -> limit[k] = den -> origin[k] + width[k];
		den -> width[k] = width[k];
	}
	return (den);
}


struct surface *read_density (FILE *fpd)
{
	int j, k, eof,  nscan;
	int origin[3], width[3], limit[3];
	long nlines, ndim, nspace, veclen;
	long x, y, z;
	long dim1, dim2, dim3;
	long ncubes;
	float *ptr;
	float value;
	double widthx, widthy, widthz, cube_width;
	double average_value, minimum_value, maximum_value;
	double bounds[2][3];
	char first_word[32], second_word[32];
	char input_line[256];
	char message[MAXLINE];
	struct surface *den;

	/* set up defaults */
	ndim = 3;
	dim1 = 1;
	dim2 = 1;
	dim3 = 1;
	nspace = 3;
	veclen = 1;
	for (k = 0; k < 3; k++) {
		bounds[0][k] = 0.0;
		bounds[1][k] = 100.0;
	}

	nlines = 0;
	eof = 0;
	while (!eof) {
		eof = get_line (input_line, MAXLINE, fpd);
		if (eof) break;
		if (input_line[0] == '#') continue;	/* later: make smarter */
		nlines++;
		nscan = sscanf (input_line, "%[a-zA-Z_01-9]%*[ =]%s", first_word, second_word);
		if (nscan < 2) {
			sprintf (message,"nscan = %d, cannot interpret: %s", nscan, input_line);
			inform (message);
			continue;
		}
		if (strcmp(first_word, "ndim") == 0) {
			ndim = atoi(second_word);
			if (ndim != 3) {
				sprintf (message, "ndim = %ld", ndim);
				set_error1 (message);
				return (NULL);
			}
		}
		else if (strcmp(first_word, "dim1") == 0) {
			dim1 = atoi(second_word);
		}
		else if (strcmp(first_word, "dim2") == 0) {
			dim2 = atoi(second_word);
		}
		else if (strcmp(first_word, "dim3") == 0) {
			dim3 = atoi(second_word);
		}
		else if (strcmp(first_word, "nspace") == 0) {
			nspace = atoi(second_word);
			if (nspace != 3) {
				sprintf (message, "nspace = %ld", nspace);
				set_error1 (message);
				return (NULL);
			}
		}
		else if (strcmp(first_word, "veclen") == 0) {
			veclen = atoi(second_word);
			if (veclen != 1) {
				sprintf (message, "veclen = %ld", veclen);
				set_error1 (message);
				return (NULL);
			}
		}
		else if (strcmp(first_word, "data") == 0) {
			if (strcmp(second_word, "float") != 0) {
				sprintf(message, "invalid data type");
				set_error1 (message);
				return (NULL);
			}
		}
		else if (strcmp(first_word, "field") == 0) {
			if (strcmp(second_word, "uniform") != 0) {
				sprintf(message, "invalid field type");
				set_error1 (message);
				return (NULL);
			}
		}
		else if (strcmp(first_word, "min_ext") == 0) {
			nscan = sscanf (input_line, "%*[a-zA-Z_01-9]%*[ =]%lf %lf %lf",
				&bounds[0][0], &bounds[0][1], &bounds[0][2]);
			if (nscan < 3) {
				sprintf(message, "invalid min_ext");
				set_error1 (message);
				return (NULL);
			}
		}
		else if (strcmp(first_word, "max_ext") == 0) {
			nscan = sscanf (input_line, "%*[a-zA-Z_01-9]%*[ =]%lf %lf %lf",
				&bounds[1][0], &bounds[1][1], &bounds[1][2]);
			if (nscan < 3) {
				sprintf(message, "invalid min_ext");
				set_error1 (message);
				return (NULL);
			}
		}
	}
	sprintf (message, "%8ld lines read from AVS-format density file", nlines);
	inform (message);
	if (nlines < 2) {
		sprintf (message, "nlines = %ld", nlines);
		set_error1 (message);
		return (NULL);
	}
	if (dim1 > DIM_MAXIMUM) {
		sprintf (message, "dim1 = %ld is too large", dim1);
		set_error1 (message);
		return (NULL);
	}
	if (dim2 > DIM_MAXIMUM) {
		sprintf (message, "dim2 = %ld is too large", dim2);
		set_error1 (message);
		return (NULL);
	}
	if (dim3 > DIM_MAXIMUM) {
		sprintf (message, "dim3 = %ld is too large", dim3);
		set_error1 (message);
		return (NULL);
	}
	width[0] = dim1;
	width[1] = dim2;
	width[2] = dim3;
	ncubes = width[0] * width[1] * width[2];
	sprintf (message, "%8d by %8d by %8d  cubes",
		width[0], width[1], width[2]);
	inform (message);
	for (j = 0; j < 2; j++)
		for (k = 0; k < 3; k++)
			if (fabs(bounds[j][k]) > 1000.0) {
				set_error1 ("invalid min_ext");
				return (NULL);
			}
	/* figure some things out */
	widthx = (bounds[1][0] - bounds[0][0]) / dim1;
	widthy = (bounds[1][1] - bounds[0][1]) / dim2;
	widthz = (bounds[1][2] - bounds[0][2]) / dim3;
	cube_width = widthx;
	sprintf (message, "%8.3f cube_width", cube_width);
	inform (message);
	if (cube_width <= 0.0 || cube_width > 10.0) return (NULL);
	if (fabs(cube_width - widthy) > 0.001) {
		set_error1 ("non-cubical cube");
		return (0);
	}
	if (fabs(cube_width - widthz) > 0.001) {
		set_error1 ("non-cubical cube");
		return (0);
	}
	for (k = 0; k < 3; k++) {
		origin[k] = floor (0.5 + bounds[0][k] / cube_width);
		limit[k] =  floor (0.5 + bounds[1][k] / cube_width);
	}

	den = create_density (cube_width, origin, width);
	ptr = den -> densities;
	average_value = 0.0;
	maximum_value = -1000.0;
	minimum_value =  1000.0;
	/* in AVS, first array index varies most quickly */
	for (z = 0; z < den -> width[2]; z++) {
		for (y = 0; y < den -> width[1]; y++) {
			for (x = 0; x < den -> width[0]; x++) {
				fread ((char *) &value, sizeof (float), 1, fpd);
				if (feof(fpd)) {
					set_error1 ("density file: premature EOF");
					sprintf (message, "%12ld floats read",
						(long) (ptr - den -> densities));
					set_error2 (message);
					return (0);
				}
				average_value += value;
				if (value < minimum_value) minimum_value = value;
				if (value > maximum_value) maximum_value = value;
				*ptr++ = value;
			}
		}
	}
	average_value /= ncubes;
	sprintf (message, "%8.4f density minimum", minimum_value);
	inform (message);
	sprintf (message, "%8.4f density average", average_value);
	inform (message);
	sprintf (message, "%8.4f density maximum", maximum_value);
	inform (message);
	return (den);
}

int get_line (char input_line[], int max_len, FILE *fpi)
{
	int nff, one;
	long idx;
	/* look for start of binary data */
	nff = 0;
	idx = 0;
	one = 1;
	while (one != EOF) {
		one = getc(fpi);
		if (one == EOF) break;
		if (one == '\014') nff++;
		if (nff >= 2) break;
		if (idx >= max_len - 1) return (1);
		input_line[idx++] = (char) one;
		if (one == '\n') break;
		if (one == '\r') break;
	}
	input_line[idx] = 0;
	if (one == EOF) return (1);
	else if (nff >= 2) return (1);
	else return(0);
}

int write_density (struct surface *den, FILE *fpd)
{
	long x, y, z;
	int k;
	float *ptr;
	float value;
	double bounds[2][3];

	for (k = 0; k < 3; k++) {
		bounds[0][k] = den -> cube_width * den -> origin[k];
		bounds[1][k] = den -> cube_width * den -> limit[k];
	}
	/* in AVS, first array index varies most quickly */
	fprintf(fpd, "# AVS field file\n");
	fprintf(fpd, "ndim=3\n");
	fprintf(fpd, "dim1=%ld\n", den -> width[0]);
	fprintf(fpd, "dim2=%ld\n", den -> width[1]);
	fprintf(fpd, "dim3=%ld\n", den -> width[2]);
	fprintf(fpd, "nspace=3\n");
	fprintf(fpd, "veclen=1\n");
	fprintf(fpd, "data=float\n");
	fprintf(fpd, "field=uniform\n");
	fprintf(fpd, "min_ext=%8.3f %8.3f %8.3f\n",
		bounds[0][0], bounds[0][1], bounds[0][2]);
	fprintf(fpd, "max_ext=%8.3f %8.3f %8.3f\n",
		bounds[1][0], bounds[1][1], bounds[1][2]);
	fprintf(fpd, "# end of header\n");
	fprintf(fpd, "\014\014");
	ptr = den -> densities;
	for (z = 0; z < den -> width[2]; z++) {
		for (y = 0; y < den -> width[1]; y++) {
			for (x = 0; x < den -> width[0]; x++) {
				value = *ptr++;
				fwrite ((char *) &value, sizeof (float), 1, fpd);
			}
		}
	}
	return (1);
}


