/***************************************************************************************************************:')

CoordSystem2Img.h

Coordinate system handling in images.

Fabrice Le Bars

Created: 2009-01-07

***************************************************************************************************************:)*/

#ifndef COORDSYSTEM2IMG_H
#define COORDSYSTEM2IMG_H

#include "CoordSystem.h"

#define KEEP_X_RATIO_COORDSYSTEM2IMG 0x00000001
#define KEEP_Y_RATIO_COORDSYSTEM2IMG 0x00000002
#define BEST_RATIO_COORDSYSTEM2IMG 0x00000004

/*
Structure that enables the use of a coordinate system in an image instead of 
using its lines and rows.
*/
struct _COORDSYSTEM2IMG
{
	COORDSYSTEM cs;
	UINT width;
	UINT height;
	double XJRatio;
	double YIRatio;
	double JXRatio;
	double IYRatio;
};
typedef struct _COORDSYSTEM2IMG COORDSYSTEM2IMG;

/*
Initialize a structure COORDSYSTEM2IMG.
See the description of the COORDSYSTEM2IMG structure for more information.

COORDSYSTEM2IMG* pCS2Img : (INOUT) Valid pointer to the structure.
COORDSYSTEM* pCS : (IN) Coordinate system to use.
UINT width : (IN) Width of the picture.
UINT height : (IN) Height of the picture.

Return : EXIT_SUCCESS or EXIT_FAILURE if there is an error.
*/
inline int InitCS2Img(COORDSYSTEM2IMG* pCS2Img, COORDSYSTEM* pCS, UINT width, UINT height)
{
	// x = xMin -> j = 0
	// x = xMax -> j = width-1
	// y = yMin -> i = height-1
	// y = yMax -> i = 0

	pCS2Img->cs = *pCS;
	pCS2Img->width = width;
	pCS2Img->height = height;

	pCS2Img->XJRatio = (pCS2Img->cs.xMax - pCS2Img->cs.xMin) / (double)(pCS2Img->width - 1);
	pCS2Img->YIRatio = (pCS2Img->cs.yMax - pCS2Img->cs.yMin) / (double)(pCS2Img->height - 1);

	pCS2Img->JXRatio = 1.0/pCS2Img->XJRatio;
	pCS2Img->IYRatio = 1.0/pCS2Img->YIRatio;

	return EXIT_SUCCESS;
}

/*
Initialize a structure COORDSYSTEM2IMG.
See the description of the COORDSYSTEM2IMG structure for more information.

COORDSYSTEM2IMG* pCS2Img : (INOUT) Valid pointer to the structure.
COORDSYSTEM* pCS : (IN) Coordinate system to use. //Can be modified to keep ratios.//
UINT width : (IN) Width of the picture.
UINT height : (IN) Height of the picture.
int flags : (IN) KEEP_X_RATIO_COORDSYSTEM2IMG to set the y ratio from x ratio, 
KEEP_Y_RATIO_COORDSYSTEM2IMG to set the x ratio from y ratio, BEST_RATIO_COORDSYSTEM2IMG 
to keep the best to avoid loosing any part of the image, 0 for free ratios.

Return : EXIT_SUCCESS or EXIT_FAILURE if there is an error.
*/
inline int InitCS2ImgEx(COORDSYSTEM2IMG* pCS2Img, COORDSYSTEM* pCS, UINT width, UINT height, int flags)
{
	// x = xMin -> j = 0
	// x = xMax -> j = width-1
	// y = yMin -> i = height-1
	// y = yMax -> i = 0

	double val = 0, temp = 0;

	pCS2Img->cs = *pCS;
	pCS2Img->width = width;
	pCS2Img->height = height;

	switch (flags)
	{
	case KEEP_X_RATIO_COORDSYSTEM2IMG:
		{
			pCS2Img->XJRatio = (pCS2Img->cs.xMax - pCS2Img->cs.xMin) / (double)(pCS2Img->width - 1);
			pCS2Img->YIRatio = pCS2Img->XJRatio;
			val = pCS2Img->YIRatio * (pCS2Img->height - 1);
			temp = (pCS2Img->cs.yMin + pCS2Img->cs.yMax);
			pCS2Img->cs.yMin = (temp-val)/2.0;
			pCS2Img->cs.yMax = (temp+val)/2.0;
			break;
		}
	case KEEP_Y_RATIO_COORDSYSTEM2IMG:
		{
			pCS2Img->YIRatio = (pCS2Img->cs.yMax - pCS2Img->cs.yMin) / (double)(pCS2Img->height - 1);
			pCS2Img->XJRatio = pCS2Img->YIRatio;
			val = pCS2Img->XJRatio * (pCS2Img->width - 1);
			temp = (pCS2Img->cs.xMin + pCS2Img->cs.xMax);
			pCS2Img->cs.xMin = (temp-val)/2.0;
			pCS2Img->cs.xMax = (temp+val)/2.0;
			break;
		}
	case BEST_RATIO_COORDSYSTEM2IMG:
		{
			pCS2Img->XJRatio = (pCS2Img->cs.xMax - pCS2Img->cs.xMin) / (double)(pCS2Img->width - 1);
			pCS2Img->YIRatio = (pCS2Img->cs.yMax - pCS2Img->cs.yMin) / (double)(pCS2Img->height - 1);

			if (pCS2Img->XJRatio < pCS2Img->YIRatio)
			{
				pCS2Img->XJRatio = pCS2Img->YIRatio;
				val = pCS2Img->XJRatio * (pCS2Img->width - 1);
				temp = (pCS2Img->cs.xMin + pCS2Img->cs.xMax);
				pCS2Img->cs.xMin = (temp-val)/2.0;
				pCS2Img->cs.xMax = (temp+val)/2.0;
			}
			else if (pCS2Img->XJRatio > pCS2Img->YIRatio)
			{
				pCS2Img->YIRatio = pCS2Img->XJRatio;
				val = pCS2Img->YIRatio * (pCS2Img->height - 1);
				temp = (pCS2Img->cs.yMin + pCS2Img->cs.yMax);
				pCS2Img->cs.yMin = (temp-val)/2.0;
				pCS2Img->cs.yMax = (temp+val)/2.0;
			}
			break;
		}
	default:
		{
			pCS2Img->XJRatio = (pCS2Img->cs.xMax - pCS2Img->cs.xMin) / (double)(pCS2Img->width - 1);
			pCS2Img->YIRatio = (pCS2Img->cs.yMax - pCS2Img->cs.yMin) / (double)(pCS2Img->height - 1);
			break;
		}
	}

	pCS2Img->JXRatio = 1.0/pCS2Img->XJRatio;
	pCS2Img->IYRatio = 1.0/pCS2Img->YIRatio;
	//*pCS = pCS2Img->cs;

	return EXIT_SUCCESS;
}

inline int XYCS2IJImg(COORDSYSTEM2IMG* pCS2Img, double x, double y, int* pI, int* pJ)
{
	*pI = (int)((pCS2Img->cs.yMax - y) * pCS2Img->IYRatio);
	*pJ = (int)((x - pCS2Img->cs.xMin) * pCS2Img->JXRatio);

	return EXIT_SUCCESS;
}

inline int IJImg2XYCS(COORDSYSTEM2IMG* pCS2Img, int i, int j, double* pX, double* pY)
{
	*pX = pCS2Img->XJRatio * j + pCS2Img->cs.xMin;
	*pY = pCS2Img->cs.yMax - pCS2Img->YIRatio * i;

	return EXIT_SUCCESS;
}

inline int XCS2JImg(COORDSYSTEM2IMG* pCS2Img, double x)
{
	return (int)((x - pCS2Img->cs.xMin) * pCS2Img->JXRatio);
}

inline int YCS2IImg(COORDSYSTEM2IMG* pCS2Img, double y)
{
	return (int)((pCS2Img->cs.yMax - y) * pCS2Img->IYRatio);
}

inline double IImg2YCS(COORDSYSTEM2IMG* pCS2Img, int i)
{
	return pCS2Img->cs.yMax - pCS2Img->YIRatio * i;
}

inline double JImg2XCS(COORDSYSTEM2IMG* pCS2Img, int j)
{
	return pCS2Img->XJRatio * j + pCS2Img->cs.xMin;
}

inline int GetCSPixelSize(COORDSYSTEM2IMG* pCS2Img, double* pSizeX, double* pSizeY)
{
	*pSizeX = pCS2Img->XJRatio;
	*pSizeY = pCS2Img->YIRatio;

	return EXIT_SUCCESS;
}

inline double GetCSPixelSizeX(COORDSYSTEM2IMG* pCS2Img)
{
	return pCS2Img->XJRatio;
}

inline double GetCSPixelSizeY(COORDSYSTEM2IMG* pCS2Img)
{
	return pCS2Img->YIRatio;
}

#endif // !COORDSYSTEM2IMG_H
