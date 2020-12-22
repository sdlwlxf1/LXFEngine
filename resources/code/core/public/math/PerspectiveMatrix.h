// Copyright 1998-2019 Epic Games, Inc. All Rights Reserved.

#pragma once

#include "math/MathUtility.h"
#include "Plane.h"
#include "Matrix.h"

namespace LXF {

class PerspectiveMatrix
	: public Matrix
{
public:

#define Z_PRECISION	0.0f

	PerspectiveMatrix(float HalfFOVX, float HalfFOVY, float MultFOVX, float MultFOVY, float MinZ, float MaxZ);
	PerspectiveMatrix(float HalfFOV, float Width, float Height, float MinZ, float MaxZ);
	PerspectiveMatrix(float HalfFOV, float Width, float Height, float MinZ);
};

class FReversedZPerspectiveMatrix : public Matrix
{
public:
	FReversedZPerspectiveMatrix(float HalfFOVX, float HalfFOVY, float MultFOVX, float MultFOVY, float MinZ, float MaxZ);
	FReversedZPerspectiveMatrix(float HalfFOV, float Width, float Height, float MinZ, float MaxZ);
	FReversedZPerspectiveMatrix(float HalfFOV, float Width, float Height, float MinZ);
};

FORCEINLINE PerspectiveMatrix::PerspectiveMatrix(float HalfFOVX, float HalfFOVY, float MultFOVX, float MultFOVY, float MinZ, float MaxZ)
	: Matrix(
		Plane(MultFOVX / Math::Tan(HalfFOVX), 0.0f, 0.0f, 0.0f),
		Plane(0.0f, MultFOVY / Math::Tan(HalfFOVY), 0.0f, 0.0f),
		Plane(0.0f, 0.0f, ((MinZ == MaxZ) ? (1.0f - Z_PRECISION) : MaxZ / (MaxZ - MinZ)), 1.0f),
		Plane(0.0f, 0.0f, -MinZ * ((MinZ == MaxZ) ? (1.0f - Z_PRECISION) : MaxZ / (MaxZ - MinZ)), 0.0f)
	)
{ }


FORCEINLINE PerspectiveMatrix::PerspectiveMatrix(float HalfFOV, float Width, float Height, float MinZ, float MaxZ)
	: Matrix(
		Plane(1.0f / Math::Tan(HalfFOV), 0.0f, 0.0f, 0.0f),
		Plane(0.0f, Width / Math::Tan(HalfFOV) / Height, 0.0f, 0.0f),
		Plane(0.0f, 0.0f, ((MinZ == MaxZ) ? (1.0f - Z_PRECISION) : MaxZ / (MaxZ - MinZ)), 1.0f),
		Plane(0.0f, 0.0f, -MinZ * ((MinZ == MaxZ) ? (1.0f - Z_PRECISION) : MaxZ / (MaxZ - MinZ)), 0.0f)
	)
{ }


FORCEINLINE PerspectiveMatrix::PerspectiveMatrix(float HalfFOV, float Width, float Height, float MinZ)
	: Matrix(
		Plane(1.0f / Math::Tan(HalfFOV), 0.0f, 0.0f, 0.0f),
		Plane(0.0f, Width / Math::Tan(HalfFOV) / Height, 0.0f, 0.0f),
		Plane(0.0f, 0.0f, (1.0f - Z_PRECISION), 1.0f),
		Plane(0.0f, 0.0f, -MinZ * (1.0f - Z_PRECISION), 0.0f)
	)
{ }


FORCEINLINE FReversedZPerspectiveMatrix::FReversedZPerspectiveMatrix(float HalfFOVX, float HalfFOVY, float MultFOVX, float MultFOVY, float MinZ, float MaxZ)
	: Matrix(
		Plane(MultFOVX / Math::Tan(HalfFOVX), 0.0f, 0.0f, 0.0f),
		Plane(0.0f, MultFOVY / Math::Tan(HalfFOVY), 0.0f, 0.0f),
		Plane(0.0f, 0.0f, ((MinZ == MaxZ) ? 0.0f : MinZ / (MinZ - MaxZ)), 1.0f),
		Plane(0.0f, 0.0f, ((MinZ == MaxZ) ? MinZ : -MaxZ * MinZ / (MinZ - MaxZ)), 0.0f)
	)
{ }


FORCEINLINE FReversedZPerspectiveMatrix::FReversedZPerspectiveMatrix(float HalfFOV, float Width, float Height, float MinZ, float MaxZ)
	: Matrix(
		Plane(1.0f / Math::Tan(HalfFOV), 0.0f, 0.0f, 0.0f),
		Plane(0.0f, Width / Math::Tan(HalfFOV) / Height, 0.0f, 0.0f),
		Plane(0.0f, 0.0f, ((MinZ == MaxZ) ? 0.0f : MinZ / (MinZ - MaxZ)), 1.0f),
		Plane(0.0f, 0.0f, ((MinZ == MaxZ) ? MinZ : -MaxZ * MinZ / (MinZ - MaxZ)), 0.0f)
	)
{ }


FORCEINLINE FReversedZPerspectiveMatrix::FReversedZPerspectiveMatrix(float HalfFOV, float Width, float Height, float MinZ)
	: Matrix(
		Plane(1.0f / Math::Tan(HalfFOV), 0.0f, 0.0f, 0.0f),
		Plane(0.0f, Width / Math::Tan(HalfFOV) / Height, 0.0f, 0.0f),
		Plane(0.0f, 0.0f, 0.0f, 1.0f),
		Plane(0.0f, 0.0f, MinZ, 0.0f)
	)
{ }


};
