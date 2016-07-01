#include "ms.h"
#include "msenum.h"
#include "mslimit.h"
#include "msdefault.h"
#include "msstruct.h"
#include "msglobal.h"
#include "msfunction.h"

/* Molecular Surface Package Copyright 1993 by Michael L. Connolly */
/* February 3, 2006 */


/* IMAGE OUTPUT */


struct {
	union {
		unsigned short short_type;                 /* Magic identifier            */
		char bytes[2];
	} type; 
    unsigned long size;                       /* File size in bytes          */    
    unsigned short int reserved1, reserved2;    
    unsigned long offset;                     /* Offset to image data, bytes */ 
} BMP_Header; 

struct {
    unsigned long size;               /* Header size in bytes      */    
	long width,height;                /* Width and height of image */    
	unsigned short int planes;       /* Number of colour planes   */    
	unsigned short int bits;         /* Bits per pixel            */    
	unsigned long compression;        /* Compression type          */    
	unsigned long imagesize;          /* Image size in bytes       */    
	long xresolution,yresolution;     /* Pixels per meter          */    
    unsigned long ncolours;           /* Number of colours         */    
    unsigned long importantcolours;   /* Important colours         */ 
} Info_BMP_Header; 

int big_endian ();
unsigned short btols (unsigned short big);
unsigned long btoll(unsigned long big);
unsigned short ltobs (unsigned short little);
unsigned long ltobl(unsigned long little);
int write_bmp_headers (FILE *fpimage, long width, long height);

/* write output frame */
void write_image (struct msscene *ms)
{
	write_image_header(ms);
	if (error()) return;
	write_color_table (ms);
	if (error()) return;
	write_pixels (ms);
	if (error()) return;
}

void write_image_header (struct msscene *ms)
{
	long i;
	char zerobyte;
	struct sgi_header sgiheader;
	struct sun_header sunheader;
	struct avs_header avsheader;
	FILE *fp_image;
	int output_format;

	fp_image = ms -> fp_raster;
	output_format = ms -> raster_format;
	zerobyte = 0;

	/* no header or color table for rgb and RGB formats */
	if (output_format == SGIIMAGE) {
		inform ("         SGI RGB Image Format Header, 24 bit color");
		sgiheader.magic = (unsigned short) 0732;
		sgiheader.storage = (unsigned char) 0;
		sgiheader.bpc = (unsigned char) 1;
		sgiheader.dimension = (unsigned short) 3;
		sgiheader.xsize = (unsigned short) (ms -> horizontal);
		sgiheader.ysize = (unsigned short) (ms -> vertical);
		sgiheader.zsize = (unsigned short) 3;
		sgiheader.min = (long) 0;
		sgiheader.max = (long) 255;
		sgiheader.dummy = (long) 0;
		strcpy(sgiheader.name, "");
		sgiheader.colormap = (long) 0;
		if (ms -> head_molecule != NULL)
			strncpy(sgiheader.name, ms -> head_molecule -> name, 79);
		sgiheader.name[79] = (char) 0;
		for (i = 0; i < 404; i++)
			sgiheader.unused[i] = (unsigned char) 0;
		if (!big_endian()) {
			sgiheader.magic = ltobs (sgiheader.magic);
			sgiheader.dimension = ltobs (sgiheader.dimension);
			sgiheader.xsize = ltobs (sgiheader.xsize);
			sgiheader.ysize = ltobs (sgiheader.ysize);
			sgiheader.zsize = ltobs (sgiheader.zsize);
			sgiheader.min = ltobl (sgiheader.min);
			sgiheader.max = ltobl (sgiheader.max);
			sgiheader.dummy = ltobl (sgiheader.dummy);
			sgiheader.colormap = ltobl (sgiheader.colormap);
		}
		fwrite ((char *) &sgiheader, sizeof (sgiheader), 1, fp_image);
	}
	else if (output_format == SUNRASTER) {
		inform ("         SUN Raster Format Header, 8 bit color");
		sunheader.magic = 0x59a66a95;
		sunheader.width = ms -> horizontal;
		sunheader.height = ms -> vertical;
		sunheader.depth = 8;
		sunheader.length = ms -> horizontal * ms -> vertical;
		sunheader.type = 1;
		sunheader.maptype = 1;
		sunheader.maplength = 3 * 256;
		if (!big_endian()) {
			sunheader.magic = ltobl (sunheader.magic);
			sunheader.width = ltobl (sunheader.width);
			sunheader.height = ltobl (sunheader.height);
			sunheader.depth = ltobl (sunheader.depth);
			sunheader.length = ltobl (sunheader.length);
			sunheader.type = ltobl (sunheader.type);
			sunheader.maptype = ltobl (sunheader.maptype);
			sunheader.maplength = ltobl (sunheader.maplength);
		}
		fwrite ((char *) &sunheader, sizeof (sunheader), 1, fp_image);
	}
	else if (output_format == AVSIMAGE) {
		inform ("         AVS Image Format Header");
		avsheader.width = ms -> horizontal;
		avsheader.height = ms -> vertical;
		if (!big_endian()) {
			avsheader.width = ltobl (avsheader.width);
			avsheader.height = ltobl (avsheader.height);
		}
		fwrite ((char *) &avsheader, sizeof (avsheader), 1, fp_image);
	}
	else if (output_format == AVSFIELD) {
		inform ("         AVS Field Format Header");
		fprintf(fp_image, "# AVS field file\n");
		fprintf(fp_image, "ndim=2\n");
		fprintf(fp_image, "dim1=%ld\n", ms -> horizontal);
		fprintf(fp_image, "dim2=%ld\n", ms -> vertical);
		fprintf(fp_image, "nspace=2\n");
		fprintf(fp_image, "veclen=4\n");
		fprintf(fp_image, "data=byte\n");
		fprintf(fp_image, "field=uniform\n");
		fprintf(fp_image, "\014\014");
	}
	else if (output_format == AVSDENSITY) {
		inform ("         AVS Density Format Header");
		/* in AVS, first array index varies most quickly */
		fprintf(fp_image, "# AVS field file\n");
		fprintf(fp_image, "ndim=3\n");
		fprintf(fp_image, "dim1=%ld\n", ms -> horizontal);
		fprintf(fp_image, "dim2=%ld\n", ms -> vertical);
		fprintf(fp_image, "dim3=%ld\n", (long) (ms -> zrange));
		fprintf(fp_image, "nspace=3\n");
		fprintf(fp_image, "veclen=1\n");
		fprintf(fp_image, "data=double\n");
		fprintf(fp_image, "field=uniform\n");
		fprintf(fp_image, "min_ext=%8.3f %8.3f %8.3f\n",
			ms -> window[0][0], ms -> window[0][1], ms -> window[0][2]);
		fprintf(fp_image, "max_ext=%8.3f %8.3f %8.3f\n",
			ms -> window[1][0], ms -> window[1][1], ms -> window[1][2]);
		fprintf(fp_image, "# end of header\n");
		fprintf(fp_image, "\014\014");
	}
	else if (output_format == BMP) {
		inform ("         Windows BMP Format Header");
		write_bmp_headers (fp_image, (long) ms -> horizontal, (long) ms -> vertical);
	}
}

void write_color_table (struct msscene *ms)
{
	FILE *fp_image;
	fp_image = ms -> fp_raster;
	/* write color table */
	if (ms -> raster_format == SUNRASTER) {
		fwrite ((char *) (ms -> table -> ucrs), sizeof (unsigned char), 256L, fp_image);
		fwrite ((char *) (ms -> table -> ucgs), sizeof (unsigned char), 256L, fp_image);
		fwrite ((char *) (ms -> table -> ucbs), sizeof (unsigned char), 256L, fp_image);
	}
}

void write_pixels (struct msscene *ms)
{
	long i, j, pidx;
	int value, this_hue, this_shade, this_alpha;
	int pass, npass, yreverse;
	int background_pixel;
	char format_name[MAXLINE];
	char message[MAX_STRING];
	double float_den;
	struct pixel_24 pixel24;
	struct pixel_32 pixel32;
	unsigned char ucr, ucg, ucb;
	unsigned char char_value;
	FILE *fp_image;
	int output_format;
	struct material_table *table;

	fp_image = ms -> fp_raster;
	output_format = ms -> raster_format;
	table = ms -> table;

	/* write pixel data */
	if (output_format == SGIIMAGE) npass = 3;
	else if (output_format == AVSDENSITY) npass = ms -> zrange;
	else npass = 1;
	yreverse = 1;
	if (output_format == AVSDENSITY) yreverse = 0;
	if (output_format == SGIIMAGE) yreverse = 0;
	if (output_format == BMP) yreverse = 0;
	for (pass = 0; pass < npass; pass++)
	for (i = 0; i < ms -> vertical; i++) 
		for (j = 0; j < ms -> horizontal; j++) {
			/* compute index into depth buffer arrays */
			if (yreverse) pidx = (ms -> vertical - i - 1) * ms -> horizontal + j;
			else pidx = i * ms -> horizontal + j;
			background_pixel = 0;
			this_hue = *(ms -> db -> hues+pidx);
			if (this_hue < 0) this_hue = 0;
			if (this_hue >= table -> nmaterial) this_hue = table -> nmaterial - 1;
			this_shade = *(ms -> db -> shades+pidx);
			if (this_shade <= 0) background_pixel = 1;
			this_alpha = *(ms -> db -> alphas+pidx);
			if (this_alpha < 0) this_alpha = 0;
			if (this_alpha >= 256) this_alpha = 255;
			if (background_pixel) {
					ucr = 0;
					ucg = 0;
					ucb = 0;
					char_value = (unsigned char) 0;
					pixel24.ucr = 0;
					pixel24.ucg = 0;
					pixel24.ucb = 0;
					pixel32.ucr = 0;
					pixel32.ucg = 0;
					pixel32.ucb = 0;
					pixel32.uct = 0;
			}
			else {
					value = hsa2value (table, this_hue, this_shade, this_alpha);
					ucr = hsa2red (table, this_hue, this_shade, this_alpha);
					if (ucr < 1) ucr = 1;
					ucg = hsa2green (table, this_hue, this_shade, this_alpha);
					if (ucg < 1) ucg = 1;
					ucb = hsa2blue (table, this_hue, this_shade, this_alpha);
					if (ucb < 1) ucb = 1;
					if (value < 0) value = 0;
					if (value >= 256) value = 255;
					if (this_alpha == 0) float_den = 0.0;
					else float_den = this_shade / 256.0;
					char_value = (unsigned char) value;
					pixel24.ucr = ucr;
					pixel24.ucg = ucg;
					pixel24.ucb = ucb;
					pixel32.ucr = ucr;
					pixel32.ucg = ucg;
					pixel32.ucb = ucb;
					pixel32.uct = (unsigned char) this_alpha;
			}
			if (output_format == SUNRASTER) {
				fwrite ((char *) &char_value, sizeof (unsigned char), 1, fp_image);
			}
			else if (output_format == AVSIMAGE || output_format == AVSFIELD) {
				fwrite ((char *) &pixel32, sizeof (pixel32), 1, fp_image);
			}
			else if (output_format == AVSDENSITY) {
				fwrite ((char *) &float_den, sizeof (float), 1, fp_image);
			}
			else if (output_format == SGIIMAGE) {
				switch (pass) {
				case 0:
					fwrite ((char *) &ucr, sizeof (unsigned char), 1, fp_image);
					break;
				case 1:
					fwrite ((char *) &ucg, sizeof (unsigned char), 1, fp_image);
					break;
				case 2:
					fwrite ((char *) &ucb, sizeof (unsigned char), 1, fp_image);
					break;
				default:
					break;
				}
			}
			else if (output_format == BMP) {
				fwrite ((char *) &ucb, sizeof (unsigned char), 1, fp_image);
				fwrite ((char *) &ucg, sizeof (unsigned char), 1, fp_image);
				fwrite ((char *) &ucr, sizeof (unsigned char), 1, fp_image);
			}
	}
	get_format_name (output_format, format_name);
	sprintf (message, "         image in %s format written to disk", format_name);
	inform (message);

}

int write_bmp_headers (FILE *fpimage, long width, long height)
{
	BMP_Header.type.bytes[0] = 'B';
	BMP_Header.type.bytes[1] = 'M';
	BMP_Header.size = 14 + 40 + 0 + 512L * 512L * 4;
	BMP_Header.reserved1 = 0;
	BMP_Header.reserved2 = 0;
	BMP_Header.offset = 14 + 40 + 0;
	
	Info_BMP_Header.size = 40;
	Info_BMP_Header.width = width;
	Info_BMP_Header.height = height;
	Info_BMP_Header.planes = 1;
	Info_BMP_Header.bits = 24;
	Info_BMP_Header.compression = 0;
	Info_BMP_Header.xresolution = 0;
	Info_BMP_Header.yresolution = 0;
	Info_BMP_Header.ncolours = 0;
	Info_BMP_Header.importantcolours = 0;
	/* swap bytes from big-endian (PowerPC) to little-endian (Intel BMP) */
	if (big_endian()) {
		BMP_Header.size = btoll (BMP_Header.size);
		BMP_Header.offset = btoll (BMP_Header.offset);
	}
	if (big_endian()) {
		Info_BMP_Header.size = btoll (Info_BMP_Header.size);
		Info_BMP_Header.width = btoll (Info_BMP_Header.width);
		Info_BMP_Header.height = btoll (Info_BMP_Header.height);
		Info_BMP_Header.planes = btols (Info_BMP_Header.planes);
		Info_BMP_Header.bits = btols (Info_BMP_Header.bits);
	}
	fwrite ((char*) &BMP_Header, 2, 1, fpimage);
	fwrite ((char*) &BMP_Header.size, sizeof (BMP_Header) - 4, 1, fpimage);
	fwrite ((char*) &Info_BMP_Header, sizeof (Info_BMP_Header), 1, fpimage);
}

unsigned short btols (unsigned short big)
{
	int c = 0;
	union {
			unsigned short big;
			char bytes[sizeof(unsigned short)];
	} swapper_big;
	union {
			char bytes[sizeof(unsigned short)];
			unsigned short little;
	} swapper_little;
	swapper_big.big = big;
	for (c = 0; c < sizeof(unsigned short); c++)
		swapper_little.bytes[c] = swapper_big.bytes[sizeof(unsigned short) - c - 1];
	return (swapper_little.little);
}

unsigned long btoll (unsigned long big)
{
	int c = 0;
	union {
			unsigned long big;
			char bytes[sizeof(unsigned long)];
	} swapper_big;
	union {
			char bytes[sizeof(unsigned long)];
			unsigned long little;
	} swapper_little;
	swapper_big.big = big;
	for (c = 0; c < sizeof(unsigned long); c++)
		swapper_little.bytes[c] = swapper_big.bytes[sizeof(unsigned long) - c - 1];
	return (swapper_little.little);
}
unsigned short ltobs (unsigned short little)
{
	int c = 0;
	union {
			unsigned short big;
			char bytes[sizeof(unsigned short)];
	} swapper_big;
	union {
			char bytes[sizeof(unsigned short)];
			unsigned short little;
	} swapper_little;
	swapper_little.little = little;
	for (c = 0; c < sizeof(unsigned short); c++)
		swapper_big.bytes[c] = swapper_little.bytes[sizeof(unsigned short) - c - 1];
	return (swapper_big.big);
}

unsigned long ltobl (unsigned long little)
{
	int c = 0;
	union {
			unsigned long big;
			char bytes[sizeof(unsigned long)];
	} swapper_big;
	union {
			char bytes[sizeof(unsigned long)];
			unsigned long little;
	} swapper_little;
	swapper_little.little = little;
	for (c = 0; c < sizeof(unsigned long); c++)
		swapper_big.bytes[c] = swapper_little.bytes[sizeof(unsigned long) - c - 1];
	return (swapper_big.big);
}

int big_endian () 
{
	union {
			int indian;
			char bytes[sizeof(int)];
	} endian;
	endian.indian = 1;
	if(endian.bytes[0] == 1) 
			return (0);
	else    return (1);
}
