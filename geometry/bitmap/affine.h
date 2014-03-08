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

/************************************************************************/
/* Public implementation follows										*/
/************************************************************************/

/*--- Defines ----------------------------------------------------------*/

#undef				FALSE
#define				FALSE			((int)(0 != 0))

#undef				TRUE
#define				TRUE			((int)(!FALSE))

#undef				ERROR
#define				ERROR			((int)TRUE)

/*--- Types ------------------------------------------------------------*/

struct	trsfStruct {
  double			dx[3];
  double			skew[3][3];
};

/*--- Functions --------------------------------------------------------*/

/**
 * Affine transformation software
Author: Philippe Thévenaz

This C routine is based on the following two papers:

    M. Unser, A. Aldroubi and M. Eden, "B-Spline Signal Processing: Part I--Theory," IEEE Transactions on Signal Processing, vol. 41, no. 2, pp. 821-832, February 1993.
    M. Unser, A. Aldroubi and M. Eden, "B-Spline Signal Processing: Part II--Efficient Design and Applications," IEEE Transactions on Signal Processing, vol. 41, no. 2, pp. 834-848, February 1993.

It implements the resampling of an image/volume under an affine transformation. The continuous model is based on splines of degree 0 (nearest neighbours), degree 1 (linear interpolation), degree 2 (quadratic), 3 (cubic), 4 (quartic), 5 (quintic), 6 and 7. By convention, the affine transformation is given by an homogenous matrix; the operation performed is output(A x) = input(x). In other words, a matrix given by
A = { {2,0,0,0}, {0,2,0,0}, {0,0,2,0}, {0,0,0,1} }
will result in a magnification by a factor 2 in linear dimensions. In case the desired operation would be output(x) = input(A x), it should be easy to modify the code (mainly: remove the call to invertTrsf() and assign invTRsf = trsf). The origin relative to which the transformation is performed is given with respect to the center of the top-upper-left voxel; the coordinate system is right-handed. Output values in need of extrapolation are set to the value background.

The input volume (the volume to transform) is given by inPtr, a pointer to an array of float values in raster order. More precisely, the values are ordered such that the x values are incremented first, then the y values, finally the z values. The size of the volume is given by nx, ny and nz, respectively. The output volume has necessarily the same size and follows the same organization. Its memory space cannot be shared with the input, and is supposed to be already allocated when the affineTransform() routine is called.

All routines are local, with the exception of the routine to call, named affineTransform(), and the routine errorReport(). The latter is not included in this distribution; its purpose is to display an error message given by a C-string. Else, the code is self-contained (provided a standard ANSI-C environment is available). It consists of only two files: affine.h and affine.c.
 */
extern	int			affineTransform	(double				transform[][4],
									 double				origin[],
									 float				*inPtr,
									 float				*outPtr,
									 long				nx,
									 long				ny,
									 long				nz,
									 int				interpolationDegree,
									 float				background);



