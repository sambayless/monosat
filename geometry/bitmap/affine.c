/************************************************************************/
/*																		*/
/* Author			Philippe Thevenaz									*/
/* Affiliation		Ecole polytechnique federale de Lausanne			*/
/*					DMT/IOA/BIG											*/
/*					P.O. Box 127										*/
/*					CH-1015 Lausanne									*/
/*					Switzerland											*/
/* Telephone		+41(21)693.51.89									*/
/* Telefax			+41(21)693.37.01									*/
/* Email			philippe.thevenaz@epfl.ch							*/
/*																		*/
/************************************************************************/

#include		<stddef.h>
#include		<stdlib.h>
#include 		<stdio.h>
#include		<string.h>
#include		<math.h>
#include		<float.h>

#include		"affine.h"

static void		errorReport		(char				*errorMessage){
	fprintf(stderr,errorMessage);
}

/************************************************************************/
static	double*		allocLinD		(long				n);
static	float*		allocLinF		(long				n);
static	float*		allocVolF		(long				nx,
									 long				ny,
									 long				nz);
static	double		BsplnWght0		(long				i,
									 double				x);
static	double		BsplnWght1		(long				i,
									 double				x);
static	double		BsplnWght2		(long				i,
									 double				x);
static	double		BsplnWght3		(long				i,
									 double				x);
static	double		BsplnWght4		(long				i,
									 double				x);
static	double		BsplnWght5		(long				i,
									 double				x);
static	double		BsplnWght6		(long				i,
									 double				x);
static	double		BsplnWght7		(long				i,
									 double				x);
static	int			directBsplineMirror
									(float				*inPtr,
									 float				*outPtr,
									 long				nx,
									 long				ny,
									 long				nz,
									 int				degree);
static	long		fold			(long				i,
									 long				n);
static	void		freeLinD		(double				*line);
static	void		freeLinF		(float				*line);
static	void		freeVolF		(float				*volume);
static	int			gaussj			(double				a[][3],
									 double				b[][3]);
static	void		getxF2D			(float				*volumeSrc,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny,
									 double				*rowDst,
									 long				n);
static	void		getxF2F			(float				*volumeSrc,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny,
									 float				*rowDst,
									 long				n);
static	void		getyF2D			(float				*volumeSrc,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny,
									 double				*columnDst,
									 long				n);
static	void		getzF2D			(float				*volumeSrc,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny,
									 double				*pillarDst,
									 long				n);
static	void		iirConvolveMirror
									(double				*inPtr,
									 double				*outPtr,
									 long				n,
									 double				gain,
									 double				*realPoles,
									 long				np);
static	int			invertTrsf		(struct trsfStruct	*dirTrsf,
									 struct trsfStruct	*invTrsf);
static	void		putxD2F			(float				*volumeDst,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny,
									 double				*rowSrc,
									 long				n);
static	void		putxF2F			(float				*volumeDst,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny,
									 float				*rowSrc,
									 long				n);
static	void		putyD2F			(float				*volumeDst,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny,
									 double				*columnSrc,
									 long				n);
static	void		putzD2F			(float				*volumeDst,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny,
									 double				*pillarSrc,
									 long				n);

/************************************************************************/
static	double*		allocLinD		(long				n) {

		return((double *)malloc((size_t)n * sizeof(double)));
} /* End of allocLinD */

/************************************************************************/
static	float*		allocLinF		(long				n) {

		return((float *)malloc((size_t)n * sizeof(float)));
} /* End of allocLinF */

/************************************************************************/
static	float*		allocVolF		(long				nx,
									 long				ny,
									 long				nz) {

		return((float *)malloc((size_t)(nx * ny * nz) * sizeof(float)));
} /* End of allocVolF */

/************************************************************************/
static	double		BsplnWght0		(long				i,
									 double				x) {

		x -= (double)i;
		return((x >= 0.5) ? (0.0) : ((x < -0.5) ? (0.0) : (1.0)));
} /* End of BsplnWght0 */

/************************************************************************/
static	double		BsplnWght1		(long				i,
									 double				x) {

		x = fabs(x - (double)i);
		return((x > 1.0) ? (0.0) : (1.0 - x));
} /* End of BsplnWght1 */

/************************************************************************/
static	double		BsplnWght2		(long				i,
									 double				x) {

		x = fabs(x - (double)i);
		if (x < 0.5)
		  return(0.75 - x * x);
		if (x < 1.5) {
		  x = 1.5 - x;
		  return(0.5 * x * x);
		}
		return(0.0);
} /* End of BsplnWght2 */

/************************************************************************/
static	double		BsplnWght3		(long				i,
									 double				x) {

		x = fabs(x - (double)i);
		if (x < 1.0)
		  return((x * x * (x - 2.0) * 3.0 + 4.0) * (1.0 / 6.0));
		if (x < 2.0) {
		  x = 2.0 - x;
		  return(x * x * x * (1.0 / 6.0));
		}
		return(0.0);
} /* End of BsplnWght3 */

/************************************************************************/
static	double		BsplnWght4		(long				i,
									 double				x) {

		x = fabs(x - (double)i);
		if (x < 0.5) {
		  x *= x;
		  return(x * (x * 0.25 - 0.625) + 115.0 / 192.0);
		}
		if (x < 1.5)
		  return(x * (x * (x * (5.0 / 6.0 - x * (1.0 / 6.0)) - 1.25) + 5.0 / 24.0) + 55.0
			/ 96.0);
		if (x < 2.5) {
		  x -= 2.5;
		  x *= x;
		  return(x * x * (1.0 / 24.0));
		}
		return(0.0);
} /* End of BsplnWght4 */

/************************************************************************/
static	double		BsplnWght5		(long				i,
									 double				x) {

		double				f;

		x = fabs(x - (double)i);
		if (x < 1.0) {
		  f = x * x;
		  return(f * (f * (0.25 - x * (1.0 / 12.0)) - 0.5) + 0.55);
		}
		if (x < 2.0)
		  return(x * (x * (x * (x * (x * (1.0 / 24.0) - 0.375) + 1.25) - 1.75) + 0.625)
			+ 0.425);
		if (x < 3.0) {
		  f = 3.0 - x;
		  x = f * f;
		  return(f * x * x * (1.0 / 120.0));
		}
		return(0.0);
} /* End of BsplnWght5 */

/************************************************************************/
static	double		BsplnWght6		(long				i,
									 double				x) {

		x = fabs(x - (double)i);
		if (x < 0.5) {
		  x *= x;
		  return(x * (x * (7.0 / 48.0 - x * (1.0 / 36.0)) - 77.0 / 192.0) + 5887.0 / 11520.0);
		}
		if (x < 1.5)
		  return(x * (x * (x * (x * (x * (x * (1.0 / 48.0) - 7.0 / 48.0) + 0.328125)
			- 35.0 / 288.0) - 91.0 / 256.0) - 7.0 / 768.0) + 7861.0 / 15360.0);
		if (x < 2.5)
		  return(x * (x * (x * (x * (x * (7.0 / 60.0 - x * (1.0 / 120.0)) - 0.65625)
			+ 133.0 / 72.0) - 2.5703125) + 1267.0 / 960.0) + 1379.0 / 7680.0);
		if (x < 3.5) {
		  x -= 3.5;
		  x *= x * x;
		  return(x * x * (1.0 / 720.0));
		}
		return(0.0);
} /* End of BsplnWght6 */

/************************************************************************/
static	double		BsplnWght7		(long				i,
									 double				x) {

		double				f;

		x = fabs(x - (double)i);
		if (x < 1.0) {
		  f = x * x;
		  return(f * (f * (f * (x * (1.0 / 144.0) - 1.0 / 36.0) + 1.0 / 9.0) - 1.0 / 3.0)
			+ 151.0 / 315.0);
		}
		if (x < 2.0)
		  return(x * (x * (x * (x * (x * (x * (0.05 - x * (1.0 / 240.0)) - 7.0 / 30.0) + 0.5)
			- 7.0 / 18.0) - 0.1) - 7.0 / 90.0) + 103.0 / 210.0);
		if (x < 3.0)
		  return(x * (x * (x * (x * (x * (x * (x * (1.0 / 720.0) - 1.0 / 36.0) + 7.0 / 30.0)
			- 19.0 / 18.0) + 49.0 / 18.0) - 23.0 / 6.0) + 217.0 / 90.0) - 139.0 / 630.0);
		if (x < 4.0) {
		  f = 4.0 - x;
		  x = f * f * f;
		  return(x * x * f * (1.0 / 5040.0));
		}
		return(0.0);
} /* End of BsplnWght7 */

/************************************************************************/
static	int			directBsplineMirror
									(float				*inPtr,
									 float				*outPtr,
									 long				nx,
									 long				ny,
									 long				nz,
									 int				degree) {

		double				*data;
		double				realPoles[3];
		double				gain;
		float				*p;
		long				np;
		long				x, y, z;

		p = allocLinF(nx);
		if (p == (float *)NULL) {
		  errorReport("ERROR - Not enough memory in directBsplineMirror");
		  return(ERROR);
		}
		for (z = 0L; (z < nz); z++)
		  for (y = 0L; (y < ny); y++) {
			getxF2F(inPtr, 0L, y, z, nx, ny, p, nx);
			putxF2F(outPtr, 0L, y, z, nx, ny, p, nx);
		  }
		freeLinF(p);
		switch (degree) {
		  case 0:
		  case 1:
			return(!ERROR);
		  case 2:
			np = 1L;
			realPoles[0] = sqrt(8.0) - 3.0;
			gain = 24.0 - sqrt(512.0);
			break;
		  case 3:
			np = 1L;
			realPoles[0] = sqrt(3.0) - 2.0;
			gain = 12.0 - sqrt(108.0);
			break;
		  case 4:
			np = 2L;
			realPoles[0] = sqrt(664.0 - sqrt(438976.0)) + sqrt(304.0) - 19.0;
			realPoles[1] = sqrt(664.0 + sqrt(438976.0)) - sqrt(304.0) - 19.0;
			gain = (1.0 - realPoles[0]) * (1.0 - realPoles[1]);
			gain *= gain;
			break;
		  case 5:
			np = 2L;
			realPoles[0] = 0.5 * (sqrt(270.0 - sqrt(70980.0)) + sqrt(105.0) - 13.0);
			realPoles[1] = 0.5 * (sqrt(270.0 + sqrt(70980.0)) - sqrt(105.0) - 13.0);
			gain = (1.0 - realPoles[0]) * (1.0 - realPoles[1]);
			gain *= gain;
			break;
		  case 6:
			np = 3L;
			realPoles[0] = -0.488294589303044755130118038883789062112279161239377608394;
			realPoles[1] = -0.081679271076237512597937765737059080653379610398148178525368;
			realPoles[2] = -0.00141415180832581775108724397655859252786416905534669851652709;
			gain = 2.598975999348577818390170516255374207847876853191217652822;
			break;
		  case 7:
			np = 3L;
			realPoles[0] = -0.5352804307964381655424037816816460718339231523426924148812;
			realPoles[1] = -0.122554615192326690515272264359357343605486549427295558490763;
			realPoles[2] = -0.0091486948096082769285930216516478534156925639545994482648003;
			gain = 3.0248282036441843886795463832305782146916878615537002580987;
			break;
		  default:
			errorReport("ERROR - Unknown degree");
			return(ERROR);
		}
		if (nx > 1L) {
		  data = allocLinD(nx);
		  if (data == (double *)NULL) {
			errorReport("ERROR - Unable to allocate x data");
			return(ERROR);
		  }
		  for (z = 0L; (z < nz); z++)
			for (y = 0L; (y < ny); y++) {
			  getxF2D(outPtr, 0L, y, z, nx, ny, data, nx);
			  iirConvolveMirror(data, data, nx, gain, realPoles, np);
			  putxD2F(outPtr, 0L, y, z, nx, ny, data, nx);
			}
		  freeLinD(data);
		}
		if (ny > 1L) {
		  data = allocLinD(ny);
		  if (data == (double *)NULL) {
			errorReport("ERROR - Unable to allocate y data");
			return(ERROR);
		  }
		  for (z = 0L; (z < nz); z++)
			for (x = 0L; (x < nx); x++) {
			  getyF2D(outPtr, x, 0L, z, nx, ny, data, ny);
			  iirConvolveMirror(data, data, ny, gain, realPoles, np);
			  putyD2F(outPtr, x, 0L, z, nx, ny, data, ny);
			}
		  freeLinD(data);
		}
		if (nz > 1L) {
		  data = allocLinD(nz);
		  if (data == (double *)NULL) {
			errorReport("ERROR - Unable to allocate z data");
			return(ERROR);
		  }
		  for (y = 0L; (y < ny); y++)
			for (x = 0L; (x < nx); x++) {
			  getzF2D(outPtr, x, y, 0L, nx, ny, data, nz);
			  iirConvolveMirror(data, data, nz, gain, realPoles, np);
			  putzD2F(outPtr, x, y, 0L, nx, ny, data, nz);
			}
		  freeLinD(data);
		}
		return(!ERROR);
} /* End of directBsplineMirror */

/************************************************************************/
static	long		fold			(long				i,
									 long				n) {

		ldiv_t				modOp;
		long				n2;

		i = labs(i);
		if (i < n)
		  return(i);
		if (n == 1L)
		  return(0L);
		n2 = (n << 1L) - 2L;
		modOp = ldiv(i, n2);
		return((modOp.rem < n) ? (modOp.rem) : (n2 - modOp.rem));
} /* End of fold */

/************************************************************************/
static	void		freeLinD		(double				*line) {

		free(line);
} /* End of freeLinD */

/************************************************************************/
static	void		freeLinF		(float				*line) {

		free(line);
} /* End of freeLinF */

/************************************************************************/
static	void		freeVolF		(float				*volume) {

		free(volume);
} /* End of freeVolF */

/************************************************************************/
static	int			gaussj			(double				a[][3],
									 double				b[][3]) {

		int					indxc[3], indxr[3];
		int					ipiv[3];
		double				big, dum, pivinv, swap;
		int					i, icol, irow, j, k, l, ll;
		int					n = 3, m = 3;

		for (j = 0; (j < n); j++)
		  ipiv[j] = 0;
		for (i = 0; (i < n); i++) {
		  big = 0.0;
		  for (j = 0; (j < n); j++)
			if (ipiv[j] != 1)
			  for (k = 0; (k < n); k++) {
				if (ipiv[k] == 0) {
				  if (fabs(a[j][k]) >= big) {
					big = fabs(a[j][k]);
					irow = j;
					icol = k;
				  }
				}
				else if (ipiv[k] > 1) {
				  errorReport("ERROR - Singular matrix in Gauss inversion");
				  return(ERROR);
				}
			  }
		  ++(ipiv[icol]);
		  if (irow != icol) {
			for (l = 0; (l < n); l++) {
			  swap = a[irow][l];
			  a[irow][l] = a[icol][l];
			  a[icol][l] = swap;
			}
			for (l = 0; (l < m); l++) {
			  swap = b[irow][l];
			  b[irow][l] = b[icol][l];
			  b[icol][l] = swap;
			}
		  }
		  indxr[i] = irow;
		  indxc[i] = icol;
		  if ((a[icol][icol] * a[icol][icol]) == 0.0) {
			errorReport("ERROR - Singular matrix in Gauss inversion");
			return(ERROR);
		  }
		  pivinv = 1.0 / a[icol][icol];
		  a[icol][icol] = 1.0;
		  for (l = 0; (l < n); l++)
			a[icol][l] *= pivinv;
		  for (l = 0; (l < m); l++)
			b[icol][l] *= pivinv;
		  for (ll = 0; (ll < n); ll++)
			if (ll != icol) {
			  dum = a[ll][icol];
			  a[ll][icol] = 0.0;
			  for (l = 0; (l < n); l++)
				a[ll][l] -= a[icol][l] * dum;
			  for (l = 0; (l < m); l++)
				b[ll][l] -= b[icol][l] * dum;
			}
		}
		for (l = n - 1; (l >= 0); l--) {
		  if (indxr[l] != indxc[l])
			for (k = 0; (k < n); k++) {
			  swap = a[k][indxr[l]];
			  a[k][indxr[l]] = a[k][indxc[l]];
			  a[k][indxc[l]] = swap;
			}
		}
		return(!ERROR);
} /* End of gaussj */

/************************************************************************/
static	void		getxF2D			(float				*volumeSrc,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny,
									 double				*rowDst,
									 long				n) {

		double				*end;

		volumeSrc += (ptrdiff_t)(x + nx * (y + ny * z));
		end = rowDst + (ptrdiff_t)n;
		while (rowDst < end)
		  *rowDst++ = (double)*volumeSrc++;
} /* End of getxF2D */

/************************************************************************/
static	void		getxF2F			(float				*volumeSrc,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny,
									 float				*rowDst,
									 long				n) {

		volumeSrc += (ptrdiff_t)(x + nx * (y + ny * z));
		rowDst = (float *)memcpy(rowDst, volumeSrc, (size_t)n * sizeof(float));
} /* End of getxF2F */

/************************************************************************/
static	void		getyF2D			(float				*volumeSrc,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny,
									 double				*columnDst,
									 long				n) {

		double				*end;

		volumeSrc += (ptrdiff_t)(x + nx * (y + ny * z));
		for (end = columnDst + (ptrdiff_t)n; (columnDst < end); volumeSrc += (ptrdiff_t)nx)
		  *columnDst++ = (double)*volumeSrc;
} /* End of getyF2D */

/************************************************************************/
static	void		getzF2D			(float				*volumeSrc,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny,
									 double				*pillarDst,
									 long				n) {

		double				*end;
		long				nxny = nx * ny;

		volumeSrc += (ptrdiff_t)(x + nx * y + nxny * z);
		for (end = pillarDst + (ptrdiff_t)n; (pillarDst < end); volumeSrc += (ptrdiff_t)nxny)
		  *pillarDst++ = (double)*volumeSrc;
} /* End of getzF2D */

/************************************************************************/
static	void		iirConvolveMirror
									(double				*inPtr,
									 double				*outPtr,
									 long				n,
									 double				gain,
									 double				*realPoles,
									 long				np) {

		double				*p, *q, *r;
		double				tolerance = log10((double)FLT_EPSILON);
		double				x0, pole;
		long				n2;
		long				i, j, k;

		p = inPtr;
		q = outPtr;
		r = outPtr + (ptrdiff_t)n;
		while (q < r)
		  *q++ = *p++ * gain;
		if (n == 1L)
		  return;
		n2 = 2L * (n - 1L);
		i = np;
		while (i-- > 0L) {
		  pole = *realPoles++;
		  j = (long)ceil(tolerance / log10(fabs(pole)));
		  k = j - n2 * (j / n2);
		  j -= k;
		  if (k < n) {
			p = outPtr + (ptrdiff_t)k;
			x0 = *p;
		  }
		  else {
			k = n2 - k;
			p = outPtr + (ptrdiff_t)k;
			x0 = *p;
			while (++p < r)
			  x0 = pole * x0 + *p;
			p--;
		  }
		  while (--p >= outPtr)
			x0 = pole * x0 + *p;
		  while (j > 0L) {
			p++;
			while (++p < r)
			  x0 = pole * x0 + *p;
			p--;
			while (--p >= outPtr)
			  x0 = pole * x0 + *p;
			j -= n2;
		  }
		  q = p++;
		  *p++ = x0;
		  x0 = *(q++ + (ptrdiff_t)n);
		  while (p < r)
			*p++ += *q++ * pole;
		  *q = (2.0 * *q - x0) / (1.0 - pole * pole);
		  p--;
		  while (--q >= outPtr)
			*q += *p-- * pole;
		}
} /* End of iirConvolveMirror */

/************************************************************************/
static	int			invertTrsf		(struct trsfStruct	*dirTrsf,
									 struct trsfStruct	*invTrsf) {

		double				dirSkew[3][3];
		int					i, j;

		for (j = 0; (j < 3); j++) {
		  for (i = 0; (i < 3); i++) {
			dirSkew[i][j] = dirTrsf->skew[i][j];
			invTrsf->skew[i][j] = 0.0;
		  }
		  invTrsf->skew[j][j] = 1.0;
		}
		if (gaussj(dirSkew, invTrsf->skew) == ERROR) {
		  errorReport("ERROR - Non reversible transformation");
		  return(ERROR);
		}
		for (j = 0; (j < 3); j++) {
		  invTrsf->dx[j] = 0.0;
		  for (i = 0; (i < 3); i++)
			invTrsf->dx[j] -= invTrsf->skew[j][i] * dirTrsf->dx[i];
		}
		return(!ERROR);
} /* End of invertTrsf */

/************************************************************************/
static	void		putxD2F			(float				*volumeDst,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny,
									 double				*rowSrc,
									 long				n) {

		double				*end;

		volumeDst += (ptrdiff_t)(x + nx * (y + ny * z));
		end = rowSrc + (ptrdiff_t)n;
		while (rowSrc < end)
		  *volumeDst++ = (float)*rowSrc++;
} /* End of putxD2F */

/************************************************************************/
static	void		putxF2F			(float				*volumeDst,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny,
									 float				*rowSrc,
									 long				n) {

		volumeDst += (ptrdiff_t)(x + nx * (y + ny * z));
		volumeDst = (float *)memcpy(volumeDst, rowSrc, (size_t)n * sizeof(float));
} /* End of putxF2F */

/************************************************************************/
static	void		putyD2F			(float				*volumeDst,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny,
									 double				*columnSrc,
									 long				n) {

		double				*end;

		volumeDst += (ptrdiff_t)(x + nx * (y + ny * z));
		for (end = columnSrc + (ptrdiff_t)n; (columnSrc < end); volumeDst += (ptrdiff_t)nx)
		  *volumeDst = (float)*columnSrc++;
} /* End of putyD2F */

/************************************************************************/
static	void		putzD2F			(float				*volumeDst,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny,
									 double				*pillarSrc,
									 long				n) {

		double				*end;
		long				nxny = nx * ny;

		volumeDst += (ptrdiff_t)(x + nx * y + nxny * z);
		for (end = pillarSrc + (ptrdiff_t)n; (pillarSrc < end); volumeDst += (ptrdiff_t)nxny)
		  *volumeDst = (float)*pillarSrc++;
} /* End of putzD2F */

/************************************************************************/
/* +===========================+										*/
/* | FUNCTION: affineTransform |										*/
/* +===========================+										*/
/*----------------------------------------------------------------------*/
/*																		*/
/* Objective:	Affine transformation of a volume						*/
/*				Resampling of a fitted spline model						*/
/*																		*/
/* Description:	- M. Unser, A. Aldroubi and M. Eden, "B-Spline Signal	*/
/*				Processing: Part I-Theory," IEEE Transactions on Signal	*/
/*				Processing, vol. 41, no. 2, pp. 821-832, February 1993.	*/
/*				- M. Unser, A. Aldroubi and M. Eden, "B-Spline Signal	*/
/*				Processing: Part II-Efficient Design and Applications,"	*/
/*				IEEE Transactions on Signal Processing, vol. 41, no. 2,	*/
/*				pp. 834-848, February 1993.								*/
/*																		*/
/* Conventions:	"transform" is a matrix of homogenous coordinates		*/
/*				"origin" is relative to the center of the				*/
/*					top-upper-left voxel								*/
/*				"inPtr" indexing scheme: first x, then y, last z		*/
/*				"outPtr" should be already allocated					*/
/*				"interpolationDegree" belongs to [0..7]					*/
/*					0 nearest neighbor									*/
/*					1 linear											*/
/*					2 quadratic											*/
/*					3 cubic												*/
/*					4 quartic											*/
/*					5 quintic											*/
/*					6 ?													*/
/*					7 ?													*/
/*																		*/
/* Comment:		Even in 2D, all entries of "transform" are considered	*/
/*				There is no antialiasing filter							*/
/*																		*/
/************************************************************************/
int					affineTransform	(double				transform[][4],
									 double				origin[],
									 float				*inPtr,
									 float				*outPtr,
									 long				nx,
									 long				ny,
									 long				nz,
									 int				interpolationDegree,
									 float				background) {

		double				*wtx, *wty, *wtz;
		float				*splinePtr;
		float				*pi, *pj;
		ptrdiff_t			*fdx, *fdy, *fdz;
		double				weightX[8], weightY[8], weightZ[8];
		ptrdiff_t			foldX[8], foldY[8], foldZ[8];
		struct trsfStruct	trsf, invTrsf;
		double				wx0, wx1, wy0, wy1;
		double				a11, a21, a31;
		double				xO, yO, zO;
		double				xz, yz, zz;
		double				xy, yy, zy;
		double				xIn, yIn, zIn;
		double				xOut, yOut, zOut;
		double				q, qi, qj;
		ptrdiff_t			fx0, fx1, fy0, fy1;
		long				half, width;
		long				nxy = nx * ny;
		long				x, y, z;
		long				i, j, k;

		if ((transform[3][0] != 0.0) || (transform[3][1] != 0.0) || (transform[3][2] != 0.0)
		  || (transform[3][3] != 1.0)){
		  errorReport("ERROR - Non-homogeneous transformation");
		  return(ERROR);
		}
		trsf.dx[0] = transform[0][3];
		trsf.dx[1] = transform[1][3];
		trsf.dx[2] = transform[2][3];
		trsf.skew[0][0] = transform[0][0];
		trsf.skew[0][1] = transform[0][1];
		trsf.skew[0][2] = transform[0][2];
		trsf.skew[1][0] = transform[1][0];
		trsf.skew[1][1] = transform[1][1];
		trsf.skew[1][2] = transform[1][2];
		trsf.skew[2][0] = transform[2][0];
		trsf.skew[2][1] = transform[2][1];
		trsf.skew[2][2] = transform[2][2];
		if (invertTrsf(&trsf, &invTrsf) == ERROR) {
		  errorReport("ERROR - Unable to compute the backward transformation");
		  return(ERROR);
		}
		xO = origin[0];
		yO = origin[1];
		zO = origin[2];
		a11 = invTrsf.skew[0][0];
		a21 = invTrsf.skew[1][0];
		a31 = invTrsf.skew[2][0];
		half = (long)interpolationDegree / 2L;
		width = (long)interpolationDegree + 1L;
		switch (interpolationDegree) {
		  case 0:
			for (zOut = -zO; (zOut < ((double)nz - zO)); zOut += 1.0) {
			  xz = invTrsf.skew[0][2] * zOut + invTrsf.dx[0] + xO;
			  yz = invTrsf.skew[1][2] * zOut + invTrsf.dx[1] + yO;
			  zz = invTrsf.skew[2][2] * zOut + invTrsf.dx[2] + zO;
			  for (yOut = -yO; (yOut < ((double)ny - yO)); yOut += 1.0) {
				xy = xz + invTrsf.skew[0][1] * yOut;
				yy = yz + invTrsf.skew[1][1] * yOut;
				zy = zz + invTrsf.skew[2][1] * yOut;
				for (xOut = -xO; (xOut < ((double)nx - xO)); xOut += 1.0) {
				  xIn = xy + a11 * xOut + 0.5;
				  yIn = yy + a21 * xOut + 0.5;
				  zIn = zy + a31 * xOut + 0.5;
				  x = (long)xIn;
				  if ((xIn < 0.0) && ((double)x != xIn))
					x--;
				  y = (long)yIn;
				  if ((yIn < 0.0) && ((double)y != yIn))
					y--;
				  z = (long)zIn;
				  if ((zIn < 0.0) && ((double)z != zIn))
					z--;
				  if ((0L <= x) && (x < nx) && (0L <= y) && (y < ny) && (0L <= z) && (z < nz))
					*outPtr++ = *(inPtr + (ptrdiff_t)(z * nxy + y * nx + x));
				  else
					*outPtr++ = background;
				}
			  }
			}
			break;
		  case 1:
			for (zOut = -zO; (zOut < ((double)nz - zO)); zOut += 1.0) {
			  xz = invTrsf.skew[0][2] * zOut + invTrsf.dx[0] + xO;
			  yz = invTrsf.skew[1][2] * zOut + invTrsf.dx[1] + yO;
			  zz = invTrsf.skew[2][2] * zOut + invTrsf.dx[2] + zO;
			  for (yOut = -yO; (yOut < ((double)ny - yO)); yOut += 1.0) {
				xy = xz + invTrsf.skew[0][1] * yOut;
				yy = yz + invTrsf.skew[1][1] * yOut;
				zy = zz + invTrsf.skew[2][1] * yOut;
				for (xOut = -xO; (xOut < ((double)nx - xO)); xOut += 1.0) {
				  xIn = xy + a11 * xOut + 0.5;
				  yIn = yy + a21 * xOut + 0.5;
				  zIn = zy + a31 * xOut + 0.5;
				  x = (long)xIn;
				  if ((xIn < 0.0) && ((double)x != xIn))
					x--;
				  y = (long)yIn;
				  if ((yIn < 0.0) && ((double)y != yIn))
					y--;
				  z = (long)zIn;
				  if ((zIn < 0.0) && ((double)z != zIn))
					z--;
				  if ((0L <= x) && (x < nx) && (0L <= y) && (y < ny) && (0L <= z) && (z < nz)) {
					yIn -= 0.5;
					y = (long)yIn;
					if ((yIn < 0.0) && ((double)y != yIn))
					  y--;
					wy0 = BsplnWght1(y, yIn);
					fy0 = (ptrdiff_t)(fold(y, ny) * nx);
					y++;
					wy1 = BsplnWght1(y, yIn);
					fy1 = (ptrdiff_t)(fold(y, ny) * nx);
					xIn -= 0.5;
					x = (long)xIn;
					if ((xIn < 0.0) && ((double)x != xIn))
					  x--;
					wx0 = BsplnWght1(x, xIn);
					fx0 = (ptrdiff_t)fold(x, nx);
					x++;
					wx1 = BsplnWght1(x, xIn);
					fx1 = (ptrdiff_t)fold(x, nx);
					if (nz == 1L) {
					  pj = inPtr;
					  pi = pj + fy0;
					  qi = wx0 * (double)*(pi + fx0);
					  qi += wx1 * (double)*(pi + fx1);
					  qj = wy0 * qi;
					  pi = pj + fy1;
					  qi = wx0 * (double)*(pi + fx0);
					  qi += wx1 * (double)*(pi + fx1);
					  q = qj + wy1 * qi;
					}
					else {
					  zIn -= 0.5;
					  z = (long)zIn;
					  if ((zIn < 0.0) && ((double)z != zIn))
						z--;
					  pj = inPtr + (ptrdiff_t)(fold(z, nz) * nxy);
					  pi = pj + fy0;
					  qi = wx0 * (double)*(pi + fx0);
					  qi += wx1 * (double)*(pi + fx1);
					  qj = wy0 * qi;
					  pi = pj + fy1;
					  qi = wx0 * (double)*(pi + fx0);
					  qi += wx1 * (double)*(pi + fx1);
					  qj += wy1 * qi;
					  q = BsplnWght1(z, zIn) * qj;
					  z++;
					  pj = inPtr + (ptrdiff_t)(fold(z, nz) * nxy);
					  pi = pj + fy0;
					  qi = wx0 * (double)*(pi + fx0);
					  qi += wx1 * (double)*(pi + fx1);
					  qj = wy0 * qi;
					  pi = pj + fy1;
					  qi = wx0 * (double)*(pi + fx0);
					  qi += wx1 * (double)*(pi + fx1);
					  qj += wy1 * qi;
					  q += BsplnWght1(z, zIn) * qj;
					}
				  }
				  else
					q = background;
				  *outPtr++ = (float)q;
				}
			  }
			}
			break;
		  case 2:
			splinePtr = allocVolF(nx, ny, nz);
			if (splinePtr == (float *)NULL) {
			  errorReport("ERROR - Not enough memory for B-spline coefficients");
			  return(ERROR);
			}
			if (directBsplineMirror(inPtr, splinePtr, nx, ny, nz, 2) == ERROR) {
			  freeVolF(splinePtr);
			  errorReport("ERROR - Unable to compute B-spline coefficients");
			  return(ERROR);
			}
			for (zOut = -zO; (zOut < ((double)nz - zO)); zOut += 1.0) {
			  xz = invTrsf.skew[0][2] * zOut + invTrsf.dx[0] + xO;
			  yz = invTrsf.skew[1][2] * zOut + invTrsf.dx[1] + yO;
			  zz = invTrsf.skew[2][2] * zOut + invTrsf.dx[2] + zO;
			  for (yOut = -yO; (yOut < ((double)ny - yO)); yOut += 1.0) {
				xy = xz + invTrsf.skew[0][1] * yOut;
				yy = yz + invTrsf.skew[1][1] * yOut;
				zy = zz + invTrsf.skew[2][1] * yOut;
				for (xOut = -xO; (xOut < ((double)nx - xO)); xOut += 1.0) {
				  xIn = xy + a11 * xOut + 0.5;
				  yIn = yy + a21 * xOut + 0.5;
				  zIn = zy + a31 * xOut + 0.5;
				  x = (long)xIn;
				  if ((xIn < 0.0) && ((double)x != xIn))
					x--;
				  y = (long)yIn;
				  if ((yIn < 0.0) && ((double)y != yIn))
					y--;
				  z = (long)zIn;
				  if ((zIn < 0.0) && ((double)z != zIn))
					z--;
				  if ((0L <= x) && (x < nx) && (0L <= y) && (y < ny) && (0L <= z) && (z < nz)) {
					if (nz == 1L) {
					  weightZ[0] = 1.0;
					  foldZ[0] = (ptrdiff_t)0;
					}
					else {
					  wtz = weightZ;
					  fdz = foldZ;
					  z -= half;
					  zIn -= 0.5;
					  for (k = z; (k < (z + width)); k++) {
						*wtz++ = BsplnWght2(k, zIn);
						*fdz++ = (ptrdiff_t)(fold(k, nz) * nxy);
					  }
					}
					wty = weightY;
					fdy = foldY;
					y -= half;
					yIn -= 0.5;
					for (j = y; (j < (y + width)); j++) {
					  *wty++ = BsplnWght2(j, yIn);
					  *fdy++ = (ptrdiff_t)(fold(j, ny) * nx);
					}
					wtx = weightX;
					fdx = foldX;
					x -= half;
					xIn -= 0.5;
					for (i = x; (i < (x + width)); i++) {
					  *wtx++ = BsplnWght2(i, xIn);
					  *fdx++ = (ptrdiff_t)fold(i, nx);
					}
					q = 0.0;
					wtz = weightZ;
					fdz = foldZ;
					for (k = (nz == 1L) ? (1L) : (width); (k-- > 0L);) {
					  wty = weightY;
					  fdy = foldY;
					  pj = splinePtr + *fdz++;
					  qj = 0.0;
					  for (j = width; (j-- > 0L);) {
						wtx = weightX;
						fdx = foldX;
						pi = pj + *fdy++;
						qi = 0.0;
						for (i = width; (i-- > 0L);)
						  qi += *wtx++ * (double)*(pi + *fdx++);
						qj += *wty++ * qi;
					  }
					  q += *wtz++ * qj;
					}
				  }
				  else
					q = background;
				  *outPtr++ = (float)q;
				}
			  }
			}
			freeVolF(splinePtr);
			break;
		  case 3:
			splinePtr = allocVolF(nx, ny, nz);
			if (splinePtr == (float *)NULL) {
			  errorReport("ERROR - Not enough memory for B-spline coefficients");
			  return(ERROR);
			}
			if (directBsplineMirror(inPtr, splinePtr, nx, ny, nz, 3) == ERROR) {
			  freeVolF(splinePtr);
			  errorReport("ERROR - Unable to compute B-spline coefficients");
			  return(ERROR);
			}
			for (zOut = -zO; (zOut < ((double)nz - zO)); zOut += 1.0) {
			  xz = invTrsf.skew[0][2] * zOut + invTrsf.dx[0] + xO;
			  yz = invTrsf.skew[1][2] * zOut + invTrsf.dx[1] + yO;
			  zz = invTrsf.skew[2][2] * zOut + invTrsf.dx[2] + zO;
			  for (yOut = -yO; (yOut < ((double)ny - yO)); yOut += 1.0) {
				xy = xz + invTrsf.skew[0][1] * yOut;
				yy = yz + invTrsf.skew[1][1] * yOut;
				zy = zz + invTrsf.skew[2][1] * yOut;
				for (xOut = -xO; (xOut < ((double)nx - xO)); xOut += 1.0) {
				  xIn = xy + a11 * xOut + 0.5;
				  yIn = yy + a21 * xOut + 0.5;
				  zIn = zy + a31 * xOut + 0.5;
				  x = (long)xIn;
				  if ((xIn < 0.0) && ((double)x != xIn))
					x--;
				  y = (long)yIn;
				  if ((yIn < 0.0) && ((double)y != yIn))
					y--;
				  z = (long)zIn;
				  if ((zIn < 0.0) && ((double)z != zIn))
					z--;
				  if ((0L <= x) && (x < nx) && (0L <= y) && (y < ny) && (0L <= z) && (z < nz)) {
					if (nz == 1L) {
					  weightZ[0] = 1.0;
					  foldZ[0] = (ptrdiff_t)0;
					}
					else {
					  wtz = weightZ;
					  fdz = foldZ;
					  zIn -= 0.5;
					  z = (long)zIn;
					  if ((zIn < 0.0) && ((double)z != zIn))
						z--;
					  z -= half;
					  for (k = z; (k < (z + width)); k++) {
						*wtz++ = BsplnWght3(k, zIn);
						*fdz++ = (ptrdiff_t)(fold(k, nz) * nxy);
					  }
					}
					wty = weightY;
					fdy = foldY;
					yIn -= 0.5;
					y = (long)yIn;
					if ((yIn < 0.0) && ((double)y != yIn))
					  y--;
					y -= half;
					for (j = y; (j < (y + width)); j++) {
					  *wty++ = BsplnWght3(j, yIn);
					  *fdy++ = (ptrdiff_t)(fold(j, ny) * nx);
					}
					wtx = weightX;
					fdx = foldX;
					xIn -= 0.5;
					x = (long)xIn;
					if ((xIn < 0.0) && ((double)x != xIn))
					  x--;
					x -= half;
					for (i = x; (i < (x + width)); i++) {
					  *wtx++ = BsplnWght3(i, xIn);
					  *fdx++ = (ptrdiff_t)fold(i, nx);
					}
					q = 0.0;
					wtz = weightZ;
					fdz = foldZ;
					for (k = (nz == 1L) ? (1L) : (width); (k-- > 0L);) {
					  wty = weightY;
					  fdy = foldY;
					  pj = splinePtr + *fdz++;
					  qj = 0.0;
					  for (j = width; (j-- > 0L);) {
						wtx = weightX;
						fdx = foldX;
						pi = pj + *fdy++;
						qi = 0.0;
						for (i = width; (i-- > 0L);)
						  qi += *wtx++ * (double)*(pi + *fdx++);
						qj += *wty++ * qi;
					  }
					  q += *wtz++ * qj;
					}
				  }
				  else
					q = background;
				  *outPtr++ = (float)q;
				}
			  }
			}
			freeVolF(splinePtr);
			break;
		  case 4:
			splinePtr = allocVolF(nx, ny, nz);
			if (splinePtr == (float *)NULL) {
			  errorReport("ERROR - Not enough memory for B-spline coefficients");
			  return(ERROR);
			}
			if (directBsplineMirror(inPtr, splinePtr, nx, ny, nz, 4) == ERROR) {
			  freeVolF(splinePtr);
			  errorReport("ERROR - Unable to compute B-spline coefficients");
			  return(ERROR);
			}
			for (zOut = -zO; (zOut < ((double)nz - zO)); zOut += 1.0) {
			  xz = invTrsf.skew[0][2] * zOut + invTrsf.dx[0] + xO;
			  yz = invTrsf.skew[1][2] * zOut + invTrsf.dx[1] + yO;
			  zz = invTrsf.skew[2][2] * zOut + invTrsf.dx[2] + zO;
			  for (yOut = -yO; (yOut < ((double)ny - yO)); yOut += 1.0) {
				xy = xz + invTrsf.skew[0][1] * yOut;
				yy = yz + invTrsf.skew[1][1] * yOut;
				zy = zz + invTrsf.skew[2][1] * yOut;
				for (xOut = -xO; (xOut < ((double)nx - xO)); xOut += 1.0) {
				  xIn = xy + a11 * xOut + 0.5;
				  yIn = yy + a21 * xOut + 0.5;
				  zIn = zy + a31 * xOut + 0.5;
				  x = (long)xIn;
				  if ((xIn < 0.0) && ((double)x != xIn))
					x--;
				  y = (long)yIn;
				  if ((yIn < 0.0) && ((double)y != yIn))
					y--;
				  z = (long)zIn;
				  if ((zIn < 0.0) && ((double)z != zIn))
					z--;
				  if ((0L <= x) && (x < nx) && (0L <= y) && (y < ny) && (0L <= z) && (z < nz)) {
					if (nz == 1L) {
					  weightZ[0] = 1.0;
					  foldZ[0] = (ptrdiff_t)0;
					}
					else {
					  wtz = weightZ;
					  fdz = foldZ;
					  z -= half;
					  zIn -= 0.5;
					  for (k = z; (k < (z + width)); k++) {
						*wtz++ = BsplnWght4(k, zIn);
						*fdz++ = (ptrdiff_t)(fold(k, nz) * nxy);
					  }
					}
					wty = weightY;
					fdy = foldY;
					y -= half;
					yIn -= 0.5;
					for (j = y; (j < (y + width)); j++) {
					  *wty++ = BsplnWght4(j, yIn);
					  *fdy++ = (ptrdiff_t)(fold(j, ny) * nx);
					}
					wtx = weightX;
					fdx = foldX;
					x -= half;
					xIn -= 0.5;
					for (i = x; (i < (x + width)); i++) {
					  *wtx++ = BsplnWght4(i, xIn);
					  *fdx++ = (ptrdiff_t)fold(i, nx);
					}
					q = 0.0;
					wtz = weightZ;
					fdz = foldZ;
					for (k = (nz == 1L) ? (1L) : (width); (k-- > 0L);) {
					  wty = weightY;
					  fdy = foldY;
					  pj = splinePtr + *fdz++;
					  qj = 0.0;
					  for (j = width; (j-- > 0L);) {
						wtx = weightX;
						fdx = foldX;
						pi = pj + *fdy++;
						qi = 0.0;
						for (i = width; (i-- > 0L);)
						  qi += *wtx++ * (double)*(pi + *fdx++);
						qj += *wty++ * qi;
					  }
					  q += *wtz++ * qj;
					}
				  }
				  else
					q = background;
				  *outPtr++ = (float)q;
				}
			  }
			}
			freeVolF(splinePtr);
			break;
		  case 5:
			splinePtr = allocVolF(nx, ny, nz);
			if (splinePtr == (float *)NULL) {
			  errorReport("ERROR - Not enough memory for B-spline coefficients");
			  return(ERROR);
			}
			if (directBsplineMirror(inPtr, splinePtr, nx, ny, nz, 5) == ERROR) {
			  freeVolF(splinePtr);
			  errorReport("ERROR - Unable to compute B-spline coefficients");
			  return(ERROR);
			}
			for (zOut = -zO; (zOut < ((double)nz - zO)); zOut += 1.0) {
			  xz = invTrsf.skew[0][2] * zOut + invTrsf.dx[0] + xO;
			  yz = invTrsf.skew[1][2] * zOut + invTrsf.dx[1] + yO;
			  zz = invTrsf.skew[2][2] * zOut + invTrsf.dx[2] + zO;
			  for (yOut = -yO; (yOut < ((double)ny - yO)); yOut += 1.0) {
				xy = xz + invTrsf.skew[0][1] * yOut;
				yy = yz + invTrsf.skew[1][1] * yOut;
				zy = zz + invTrsf.skew[2][1] * yOut;
				for (xOut = -xO; (xOut < ((double)nx - xO)); xOut += 1.0) {
				  xIn = xy + a11 * xOut + 0.5;
				  yIn = yy + a21 * xOut + 0.5;
				  zIn = zy + a31 * xOut + 0.5;
				  x = (long)xIn;
				  if ((xIn < 0.0) && ((double)x != xIn))
					x--;
				  y = (long)yIn;
				  if ((yIn < 0.0) && ((double)y != yIn))
					y--;
				  z = (long)zIn;
				  if ((zIn < 0.0) && ((double)z != zIn))
					z--;
				  if ((0L <= x) && (x < nx) && (0L <= y) && (y < ny) && (0L <= z) && (z < nz)) {
					if (nz == 1L) {
					  weightZ[0] = 1.0;
					  foldZ[0] = (ptrdiff_t)0;
					}
					else {
					  wtz = weightZ;
					  fdz = foldZ;
					  zIn -= 0.5;
					  z = (long)zIn;
					  if ((zIn < 0.0) && ((double)z != zIn))
						z--;
					  z -= half;
					  for (k = z; (k < (z + width)); k++) {
						*wtz++ = BsplnWght5(k, zIn);
						*fdz++ = (ptrdiff_t)(fold(k, nz) * nxy);
					  }
					}
					wty = weightY;
					fdy = foldY;
					yIn -= 0.5;
					y = (long)yIn;
					if ((yIn < 0.0) && ((double)y != yIn))
					  y--;
					y -= half;
					for (j = y; (j < (y + width)); j++) {
					  *wty++ = BsplnWght5(j, yIn);
					  *fdy++ = (ptrdiff_t)(fold(j, ny) * nx);
					}
					wtx = weightX;
					fdx = foldX;
					xIn -= 0.5;
					x = (long)xIn;
					if ((xIn < 0.0) && ((double)x != xIn))
					  x--;
					x -= half;
					for (i = x; (i < (x + width)); i++) {
					  *wtx++ = BsplnWght5(i, xIn);
					  *fdx++ = (ptrdiff_t)fold(i, nx);
					}
					q = 0.0;
					wtz = weightZ;
					fdz = foldZ;
					for (k = (nz == 1L) ? (1L) : (width); (k-- > 0L);) {
					  wty = weightY;
					  fdy = foldY;
					  pj = splinePtr + *fdz++;
					  qj = 0.0;
					  for (j = width; (j-- > 0L);) {
						wtx = weightX;
						fdx = foldX;
						pi = pj + *fdy++;
						qi = 0.0;
						for (i = width; (i-- > 0L);)
						  qi += *wtx++ * (double)*(pi + *fdx++);
						qj += *wty++ * qi;
					  }
					  q += *wtz++ * qj;
					}
				  }
				  else
					q = background;
				  *outPtr++ = (float)q;
				}
			  }
			}
			freeVolF(splinePtr);
			break;
		  case 6:
			splinePtr = allocVolF(nx, ny, nz);
			if (splinePtr == (float *)NULL) {
			  errorReport("ERROR - Not enough memory for B-spline coefficients");
			  return(ERROR);
			}
			if (directBsplineMirror(inPtr, splinePtr, nx, ny, nz, 6) == ERROR) {
			  freeVolF(splinePtr);
			  errorReport("ERROR - Unable to compute B-spline coefficients");
			  return(ERROR);
			}
			for (zOut = -zO; (zOut < ((double)nz - zO)); zOut += 1.0) {
			  xz = invTrsf.skew[0][2] * zOut + invTrsf.dx[0] + xO;
			  yz = invTrsf.skew[1][2] * zOut + invTrsf.dx[1] + yO;
			  zz = invTrsf.skew[2][2] * zOut + invTrsf.dx[2] + zO;
			  for (yOut = -yO; (yOut < ((double)ny - yO)); yOut += 1.0) {
				xy = xz + invTrsf.skew[0][1] * yOut;
				yy = yz + invTrsf.skew[1][1] * yOut;
				zy = zz + invTrsf.skew[2][1] * yOut;
				for (xOut = -xO; (xOut < ((double)nx - xO)); xOut += 1.0) {
				  xIn = xy + a11 * xOut + 0.5;
				  yIn = yy + a21 * xOut + 0.5;
				  zIn = zy + a31 * xOut + 0.5;
				  x = (long)xIn;
				  if ((xIn < 0.0) && ((double)x != xIn))
					x--;
				  y = (long)yIn;
				  if ((yIn < 0.0) && ((double)y != yIn))
					y--;
				  z = (long)zIn;
				  if ((zIn < 0.0) && ((double)z != zIn))
					z--;
				  if ((0L <= x) && (x < nx) && (0L <= y) && (y < ny) && (0L <= z) && (z < nz)) {
					if (nz == 1L) {
					  weightZ[0] = 1.0;
					  foldZ[0] = (ptrdiff_t)0;
					}
					else {
					  wtz = weightZ;
					  fdz = foldZ;
					  z -= half;
					  zIn -= 0.5;
					  for (k = z; (k < (z + width)); k++) {
						*wtz++ = BsplnWght6(k, zIn);
						*fdz++ = (ptrdiff_t)(fold(k, nz) * nxy);
					  }
					}
					wty = weightY;
					fdy = foldY;
					y -= half;
					yIn -= 0.5;
					for (j = y; (j < (y + width)); j++) {
					  *wty++ = BsplnWght6(j, yIn);
					  *fdy++ = (ptrdiff_t)(fold(j, ny) * nx);
					}
					wtx = weightX;
					fdx = foldX;
					x -= half;
					xIn -= 0.5;
					for (i = x; (i < (x + width)); i++) {
					  *wtx++ = BsplnWght6(i, xIn);
					  *fdx++ = (ptrdiff_t)fold(i, nx);
					}
					q = 0.0;
					wtz = weightZ;
					fdz = foldZ;
					for (k = (nz == 1L) ? (1L) : (width); (k-- > 0L);) {
					  wty = weightY;
					  fdy = foldY;
					  pj = splinePtr + *fdz++;
					  qj = 0.0;
					  for (j = width; (j-- > 0L);) {
						wtx = weightX;
						fdx = foldX;
						pi = pj + *fdy++;
						qi = 0.0;
						for (i = width; (i-- > 0L);)
						  qi += *wtx++ * (double)*(pi + *fdx++);
						qj += *wty++ * qi;
					  }
					  q += *wtz++ * qj;
					}
				  }
				  else
					q = background;
				  *outPtr++ = (float)q;
				}
			  }
			}
			freeVolF(splinePtr);
			break;
		  case 7:
			splinePtr = allocVolF(nx, ny, nz);
			if (splinePtr == (float *)NULL) {
			  errorReport("ERROR - Not enough memory for B-spline coefficients");
			  return(ERROR);
			}
			if (directBsplineMirror(inPtr, splinePtr, nx, ny, nz, 7) == ERROR) {
			  freeVolF(splinePtr);
			  errorReport("ERROR - Unable to compute B-spline coefficients");
			  return(ERROR);
			}
			for (zOut = -zO; (zOut < ((double)nz - zO)); zOut += 1.0) {
			  xz = invTrsf.skew[0][2] * zOut + invTrsf.dx[0] + xO;
			  yz = invTrsf.skew[1][2] * zOut + invTrsf.dx[1] + yO;
			  zz = invTrsf.skew[2][2] * zOut + invTrsf.dx[2] + zO;
			  for (yOut = -yO; (yOut < ((double)ny - yO)); yOut += 1.0) {
				xy = xz + invTrsf.skew[0][1] * yOut;
				yy = yz + invTrsf.skew[1][1] * yOut;
				zy = zz + invTrsf.skew[2][1] * yOut;
				for (xOut = -xO; (xOut < ((double)nx - xO)); xOut += 1.0) {
				  xIn = xy + a11 * xOut + 0.5;
				  yIn = yy + a21 * xOut + 0.5;
				  zIn = zy + a31 * xOut + 0.5;
				  x = (long)xIn;
				  if ((xIn < 0.0) && ((double)x != xIn))
					x--;
				  y = (long)yIn;
				  if ((yIn < 0.0) && ((double)y != yIn))
					y--;
				  z = (long)zIn;
				  if ((zIn < 0.0) && ((double)z != zIn))
					z--;
				  if ((0L <= x) && (x < nx) && (0L <= y) && (y < ny) && (0L <= z) && (z < nz)) {
					if (nz == 1L) {
					  weightZ[0] = 1.0;
					  foldZ[0] = (ptrdiff_t)0;
					}
					else {
					  wtz = weightZ;
					  fdz = foldZ;
					  zIn -= 0.5;
					  z = (long)zIn;
					  if ((zIn < 0.0) && ((double)z != zIn))
						z--;
					  z -= half;
					  for (k = z; (k < (z + width)); k++) {
						*wtz++ = BsplnWght7(k, zIn);
						*fdz++ = (ptrdiff_t)(fold(k, nz) * nxy);
					  }
					}
					wty = weightY;
					fdy = foldY;
					yIn -= 0.5;
					y = (long)yIn;
					if ((yIn < 0.0) && ((double)y != yIn))
					  y--;
					y -= half;
					for (j = y; (j < (y + width)); j++) {
					  *wty++ = BsplnWght7(j, yIn);
					  *fdy++ = (ptrdiff_t)(fold(j, ny) * nx);
					}
					wtx = weightX;
					fdx = foldX;
					xIn -= 0.5;
					x = (long)xIn;
					if ((xIn < 0.0) && ((double)x != xIn))
					  x--;
					x -= half;
					for (i = x; (i < (x + width)); i++) {
					  *wtx++ = BsplnWght7(i, xIn);
					  *fdx++ = (ptrdiff_t)fold(i, nx);
					}
					q = 0.0;
					wtz = weightZ;
					fdz = foldZ;
					for (k = (nz == 1L) ? (1L) : (width); (k-- > 0L);) {
					  wty = weightY;
					  fdy = foldY;
					  pj = splinePtr + *fdz++;
					  qj = 0.0;
					  for (j = width; (j-- > 0L);) {
						wtx = weightX;
						fdx = foldX;
						pi = pj + *fdy++;
						qi = 0.0;
						for (i = width; (i-- > 0L);)
						  qi += *wtx++ * (double)*(pi + *fdx++);
						qj += *wty++ * qi;
					  }
					  q += *wtz++ * qj;
					}
				  }
				  else
					q = background;
				  *outPtr++ = (float)q;
				}
			  }
			}
			freeVolF(splinePtr);
			break;
		  default:
			errorReport("ERROR - Unknown interpolation specification");
			return(ERROR);
		}
		return(!ERROR);
} /* End of affineTransform */
