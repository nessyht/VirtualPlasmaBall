#pragma once
#include "GridXd.h"

typedef enum
{
	CIRCULAR,
	SPERICAL
} Geometry;

typedef enum
{
	ELECTRODE,
	GLASS,
	GAS,
	PLASMA,
	AIR,
	BORDER,
	FAKEELECTRODE,
	CONDUCTOR
} Material;

