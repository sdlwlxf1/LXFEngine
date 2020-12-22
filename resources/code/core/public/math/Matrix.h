// this tiny game engine is rendered by CPU. I make it for learing the pipline and rasterization. This is the first step to realise my own game engine. so just do it!

#pragma once

#include "CoreTypes.h"
#include "MathUtility.h"
#include "Plane.h"
#include "MathSSE.h"

namespace LXF {

class Matrix
{
public:
	union {
		float M[4][4];
	};

	static const Matrix Identity;

	FORCEINLINE Matrix() { Memory::Memzero(this, sizeof(*this)); }
	FORCEINLINE Matrix(const Plane& InX, const Plane& InY, const Plane& InZ, const Plane& InW);
	FORCEINLINE Matrix(const Vector3f& InX, const Vector3f& InY, const Vector3f& InZ, const Vector3f& InW);
	FORCEINLINE void SetIdentity();

	FORCEINLINE void operator*=(const Matrix& Other);
	FORCEINLINE Matrix operator*(const Matrix& Other) const;
	FORCEINLINE Matrix operator+(const Matrix& Other) const;

	FORCEINLINE void operator+=(const Matrix& Other) { *this = *this + Other; }

	FORCEINLINE Matrix operator*(float Other) const;
	FORCEINLINE void operator*=(float Other) { *this = *this*Other; }

	inline bool operator==(const Matrix& Other) const;

	FORCEINLINE bool Equals(const Matrix& Other, float Tolerance = KINDA_SMALL_NUMBER) const;

	FORCEINLINE bool operator!=(const Matrix& Other) const { return !(*this == Other); }

	FORCEINLINE Matrix GetTransposed() const;

	inline float Determinant() const;
	inline Matrix Inverse() const;
	inline Matrix InverseFast() const;
	inline void ScaleTranslation(const Vector3f& InScale3D);
	inline Vector3f GetOrigin() const;
	inline Vector3f GetScaledAxis(EAxis::Type InAxis) const;
	inline void GetScaledAxes(Vector3f &X, Vector3f &Y, Vector3f &Z) const;
	inline Vector3f GetUnitAxis(EAxis::Type Axis) const;
	inline void GetUnitAxes(Vector3f &X, Vector3f &Y, Vector3f &Z) const;
	inline void SetAxis(int32 i, const Vector3f& Axis);
	inline void SetOrigin(const Vector3f& NewOrigin);
	inline void SetAxes(Vector3f* Axis0 = NULL, Vector3f* Axis1 = NULL, Vector3f* Axis2 = NULL, Vector3f* Origin = NULL);

	inline void RemoveScaling(float Tolerance = SMALL_NUMBER);
	inline Matrix GetMatrixWithoutScale(float Tolerance = SMALL_NUMBER) const;
	inline Vector3f ExtractScaling(float Tolerance = SMALL_NUMBER);
	inline Vector3f GetScaleVector(float Tolerance = SMALL_NUMBER) const;

	FORCEINLINE bool GetFrustumNearPlane(Plane& OutPlane) const;
	FORCEINLINE bool GetFrustumFarPlane(Plane& OutPlane) const;
	FORCEINLINE bool GetFrustumLeftPlane(Plane& OutPlane) const;
	FORCEINLINE bool GetFrustumRightPlane(Plane& OutPlane) const;
	FORCEINLINE bool GetFrustumTopPlane(Plane& OutPlane) const;
	FORCEINLINE bool GetFrustumBottomPlane(Plane& OutPlane) const;

	inline void Mirror(EAxis::Type MirrorAxis, EAxis::Type FlipAxis);

	FORCEINLINE Vector4f TransformVector4f(const Vector4f& V) const;
	FORCEINLINE Vector4f TransformPosition(const Vector3f &V) const;
	FORCEINLINE Vector3f InverseTransformPosition(const Vector3f &V) const;
	FORCEINLINE Vector4f TransformVector(const Vector3f& V) const;
	FORCEINLINE Vector3f InverseTransformVector(const Vector3f &V) const;
};

inline bool Matrix::operator==(const Matrix& Other) const
{
	for (int32 X = 0; X < 4; X++)
	{
		for (int32 Y = 0; Y < 4; Y++)
		{
			if (M[X][Y] != Other.M[X][Y])
			{
				return false;
			}
		}
	}

	return true;
}

FORCEINLINE Matrix::Matrix(const Plane& InX, const Plane& InY, const Plane& InZ, const Plane& InW)
{
	M[0][0] = InX.X; M[0][1] = InX.Y;  M[0][2] = InX.Z;  M[0][3] = InX.W;
	M[1][0] = InY.X; M[1][1] = InY.Y;  M[1][2] = InY.Z;  M[1][3] = InY.W;
	M[2][0] = InZ.X; M[2][1] = InZ.Y;  M[2][2] = InZ.Z;  M[2][3] = InZ.W;
	M[3][0] = InW.X; M[3][1] = InW.Y;  M[3][2] = InW.Z;  M[3][3] = InW.W;
}

FORCEINLINE Matrix::Matrix(const Vector3f& InX, const Vector3f& InY, const Vector3f& InZ, const Vector3f& InW)
{
	M[0][0] = InX.X; M[0][1] = InX.Y;  M[0][2] = InX.Z;  M[0][3] = 0.0f;
	M[1][0] = InY.X; M[1][1] = InY.Y;  M[1][2] = InY.Z;  M[1][3] = 0.0f;
	M[2][0] = InZ.X; M[2][1] = InZ.Y;  M[2][2] = InZ.Z;  M[2][3] = 0.0f;
	M[3][0] = InW.X; M[3][1] = InW.Y;  M[3][2] = InW.Z;  M[3][3] = 1.0f;
}

FORCEINLINE void Matrix::SetIdentity()
{
	M[0][0] = 1; M[0][1] = 0;  M[0][2] = 0;  M[0][3] = 0;
	M[1][0] = 0; M[1][1] = 1;  M[1][2] = 0;  M[1][3] = 0;
	M[2][0] = 0; M[2][1] = 0;  M[2][2] = 1;  M[2][3] = 0;
	M[3][0] = 0; M[3][1] = 0;  M[3][2] = 0;  M[3][3] = 1;
}

FORCEINLINE Matrix Matrix::operator*(float Other) const
{
	Matrix ResultMat;

	for (int32 X = 0; X < 4; X++)
	{
		for (int32 Y = 0; Y < 4; Y++)
		{
			ResultMat.M[X][Y] = M[X][Y] * Other;
		}
	}

	return ResultMat;
}

FORCEINLINE void Matrix::operator*=(const Matrix& Other)
{
	VectorMatrixMultiply(this, this, &Other);
}

FORCEINLINE Matrix Matrix::operator*(const Matrix& Other) const
{
	Matrix Result;
	VectorMatrixMultiply(&Result, this, &Other);
	return Result;
}

FORCEINLINE Matrix Matrix::operator+(const Matrix& Other) const
{
	Matrix ResultMat;

	for (int32 X = 0; X < 4; X++)
	{
		for (int32 Y = 0; Y < 4; Y++)
		{
			ResultMat.M[X][Y] = M[X][Y] + Other.M[X][Y];
		}
	}

	return ResultMat;
}

FORCEINLINE bool Matrix::Equals(const Matrix& Other, float Tolerance /*= KINDA_SMALL_NUMBER*/) const
{
	for (int32 X = 0; X < 4; X++)
	{
		for (int32 Y = 0; Y < 4; Y++)
		{
			if (Math::Abs(M[X][Y] - Other.M[X][Y]) > Tolerance)
			{
				return false;
			}
		}
	}

	return true;
}

FORCEINLINE Matrix Matrix::GetTransposed() const
{
	Matrix	Result;

	Result.M[0][0] = M[0][0];
	Result.M[0][1] = M[1][0];
	Result.M[0][2] = M[2][0];
	Result.M[0][3] = M[3][0];

	Result.M[1][0] = M[0][1];
	Result.M[1][1] = M[1][1];
	Result.M[1][2] = M[2][1];
	Result.M[1][3] = M[3][1];

	Result.M[2][0] = M[0][2];
	Result.M[2][1] = M[1][2];
	Result.M[2][2] = M[2][2];
	Result.M[2][3] = M[3][2];

	Result.M[3][0] = M[0][3];
	Result.M[3][1] = M[1][3];
	Result.M[3][2] = M[2][3];
	Result.M[3][3] = M[3][3];

	return Result;
}

inline float Matrix::Determinant() const
{
	return	M[0][0] * (
		M[1][1] * (M[2][2] * M[3][3] - M[2][3] * M[3][2]) -
		M[2][1] * (M[1][2] * M[3][3] - M[1][3] * M[3][2]) +
		M[3][1] * (M[1][2] * M[2][3] - M[1][3] * M[2][2])
		) -
		M[1][0] * (
			M[0][1] * (M[2][2] * M[3][3] - M[2][3] * M[3][2]) -
			M[2][1] * (M[0][2] * M[3][3] - M[0][3] * M[3][2]) +
			M[3][1] * (M[0][2] * M[2][3] - M[0][3] * M[2][2])
			) +
		M[2][0] * (
			M[0][1] * (M[1][2] * M[3][3] - M[1][3] * M[3][2]) -
			M[1][1] * (M[0][2] * M[3][3] - M[0][3] * M[3][2]) +
			M[3][1] * (M[0][2] * M[1][3] - M[0][3] * M[1][2])
			) -
		M[3][0] * (
			M[0][1] * (M[1][2] * M[2][3] - M[1][3] * M[2][2]) -
			M[1][1] * (M[0][2] * M[2][3] - M[0][3] * M[2][2]) +
			M[2][1] * (M[0][2] * M[1][3] - M[0][3] * M[1][2])
			);
}

inline Matrix Matrix::Inverse() const
{
	Matrix Result;

	// Check for zero scale matrix to invert
	if (GetScaledAxis(EAxis::X).IsNearlyZero(SMALL_NUMBER) &&
		GetScaledAxis(EAxis::Y).IsNearlyZero(SMALL_NUMBER) &&
		GetScaledAxis(EAxis::Z).IsNearlyZero(SMALL_NUMBER))
	{
		// just set to zero - avoids unsafe inverse of zero and duplicates what QNANs were resulting in before (scaling away all children)
		Result = Matrix::Identity;
	}
	else
	{
		const float	Det = Determinant();

		if (Det == 0.0f)
		{
			Result = Matrix::Identity;
		}
		else
		{
			VectorMatrixInverse(&Result, this);
		}
	}

	return Result;
}

inline Matrix Matrix::InverseFast() const
{
#if (!BUILD_SHIPPING)
	// Check for zero scale matrix to invert
	if (GetScaledAxis(EAxis::X).IsNearlyZero(SMALL_NUMBER) &&
		GetScaledAxis(EAxis::Y).IsNearlyZero(SMALL_NUMBER) &&
		GetScaledAxis(EAxis::Z).IsNearlyZero(SMALL_NUMBER))
	{
	}
	else
	{
		const float	Det = Determinant();

		if (Det == 0.0f || !Math::IsFinite(Det))
		{
		}
	}
#endif
	Matrix Result;
	VectorMatrixInverse(&Result, this);
	return Result;
}

inline void Matrix::ScaleTranslation(const Vector3f& InScale3D)
{
	M[3][0] *= InScale3D.X;
	M[3][1] *= InScale3D.Y;
	M[3][2] *= InScale3D.Z;
}

inline Vector3f Matrix::GetOrigin() const
{
	return Vector3f(M[3][0], M[3][1], M[3][2]);
}

inline Vector3f Matrix::GetScaledAxis(EAxis::Type InAxis) const
{
	switch (InAxis)
	{
	case EAxis::X:
		return Vector3f(M[0][0], M[0][1], M[0][2]);

	case EAxis::Y:
		return Vector3f(M[1][0], M[1][1], M[1][2]);

	case EAxis::Z:
		return Vector3f(M[2][0], M[2][1], M[2][2]);

	default:
		return Vector3f::ZeroVector;
	}
}

inline void Matrix::GetScaledAxes(Vector3f &X, Vector3f &Y, Vector3f &Z) const
{
	X.X = M[0][0]; X.Y = M[0][1]; X.Z = M[0][2];
	Y.X = M[1][0]; Y.Y = M[1][1]; Y.Z = M[1][2];
	Z.X = M[2][0]; Z.Y = M[2][1]; Z.Z = M[2][2];
}

inline Vector3f Matrix::GetUnitAxis(EAxis::Type InAxis) const
{
	return GetScaledAxis(InAxis).GetSafeNormal();
}

inline void Matrix::GetUnitAxes(Vector3f &X, Vector3f &Y, Vector3f &Z) const
{
	GetScaledAxes(X, Y, Z);
	X.Normalize();
	Y.Normalize();
	Z.Normalize();
}

inline void Matrix::SetAxis(int32 i, const Vector3f& Axis)
{
	checkSlow(i >= 0 && i <= 2);
	M[i][0] = Axis.X;
	M[i][1] = Axis.Y;
	M[i][2] = Axis.Z;
}

inline void Matrix::SetOrigin(const Vector3f& NewOrigin)
{
	M[3][0] = NewOrigin.X;
	M[3][1] = NewOrigin.Y;
	M[3][2] = NewOrigin.Z;
}

inline void Matrix::SetAxes(Vector3f* Axis0 /*= NULL*/, Vector3f* Axis1 /*= NULL*/, Vector3f* Axis2 /*= NULL*/, Vector3f* Origin /*= NULL*/)
{
	if (Axis0 != NULL)
	{
		M[0][0] = Axis0->X;
		M[0][1] = Axis0->Y;
		M[0][2] = Axis0->Z;
	}
	if (Axis1 != NULL)
	{
		M[1][0] = Axis1->X;
		M[1][1] = Axis1->Y;
		M[1][2] = Axis1->Z;
	}
	if (Axis2 != NULL)
	{
		M[2][0] = Axis2->X;
		M[2][1] = Axis2->Y;
		M[2][2] = Axis2->Z;
	}
	if (Origin != NULL)
	{
		M[3][0] = Origin->X;
		M[3][1] = Origin->Y;
		M[3][2] = Origin->Z;
	}
}

inline void Matrix::RemoveScaling(float Tolerance/*=SMALL_NUMBER*/)
{
	// For each row, find magnitude, and if its non-zero re-scale so its unit length.
	const float SquareSum0 = (M[0][0] * M[0][0]) + (M[0][1] * M[0][1]) + (M[0][2] * M[0][2]);
	const float SquareSum1 = (M[1][0] * M[1][0]) + (M[1][1] * M[1][1]) + (M[1][2] * M[1][2]);
	const float SquareSum2 = (M[2][0] * M[2][0]) + (M[2][1] * M[2][1]) + (M[2][2] * M[2][2]);
	const float Scale0 = Math::FloatSelect(SquareSum0 - Tolerance, Math::InvSqrt(SquareSum0), 1.0f);
	const float Scale1 = Math::FloatSelect(SquareSum1 - Tolerance, Math::InvSqrt(SquareSum1), 1.0f);
	const float Scale2 = Math::FloatSelect(SquareSum2 - Tolerance, Math::InvSqrt(SquareSum2), 1.0f);
	M[0][0] *= Scale0;
	M[0][1] *= Scale0;
	M[0][2] *= Scale0;
	M[1][0] *= Scale1;
	M[1][1] *= Scale1;
	M[1][2] *= Scale1;
	M[2][0] *= Scale2;
	M[2][1] *= Scale2;
	M[2][2] *= Scale2;
}

inline Matrix Matrix::GetMatrixWithoutScale(float Tolerance/*=SMALL_NUMBER*/) const
{
	Matrix Result = *this;
	Result.RemoveScaling(Tolerance);
	return Result;
}

inline Vector3f Matrix::ExtractScaling(float Tolerance/*=SMALL_NUMBER*/)
{
	Vector3f Scale3D(0, 0, 0);

	// For each row, find magnitude, and if its non-zero re-scale so its unit length.
	const float SquareSum0 = (M[0][0] * M[0][0]) + (M[0][1] * M[0][1]) + (M[0][2] * M[0][2]);
	const float SquareSum1 = (M[1][0] * M[1][0]) + (M[1][1] * M[1][1]) + (M[1][2] * M[1][2]);
	const float SquareSum2 = (M[2][0] * M[2][0]) + (M[2][1] * M[2][1]) + (M[2][2] * M[2][2]);

	if (SquareSum0 > Tolerance)
	{
		float Scale0 = Math::Sqrt(SquareSum0);
		Scale3D[0] = Scale0;
		float InvScale0 = 1.f / Scale0;
		M[0][0] *= InvScale0;
		M[0][1] *= InvScale0;
		M[0][2] *= InvScale0;
	}
	else
	{
		Scale3D[0] = 0;
	}

	if (SquareSum1 > Tolerance)
	{
		float Scale1 = Math::Sqrt(SquareSum1);
		Scale3D[1] = Scale1;
		float InvScale1 = 1.f / Scale1;
		M[1][0] *= InvScale1;
		M[1][1] *= InvScale1;
		M[1][2] *= InvScale1;
	}
	else
	{
		Scale3D[1] = 0;
	}

	if (SquareSum2 > Tolerance)
	{
		float Scale2 = Math::Sqrt(SquareSum2);
		Scale3D[2] = Scale2;
		float InvScale2 = 1.f / Scale2;
		M[2][0] *= InvScale2;
		M[2][1] *= InvScale2;
		M[2][2] *= InvScale2;
	}
	else
	{
		Scale3D[2] = 0;
	}

	return Scale3D;
}

/** return a 3D scale vector calculated from this matrix (where each component is the magnitude of a row vector). */
inline Vector3f Matrix::GetScaleVector(float Tolerance/*=SMALL_NUMBER*/) const
{
	Vector3f Scale3D(1, 1, 1);

	// For each row, find magnitude, and if its non-zero re-scale so its unit length.
	for (int32 i = 0; i < 3; i++)
	{
		const float SquareSum = (M[i][0] * M[i][0]) + (M[i][1] * M[i][1]) + (M[i][2] * M[i][2]);
		if (SquareSum > Tolerance)
		{
			Scale3D[i] = Math::Sqrt(SquareSum);
		}
		else
		{
			Scale3D[i] = 0.f;
		}
	}

	return Scale3D;
}

FORCEINLINE bool MakeFrustumPlane(float A, float B, float C, float D, Plane& OutPlane)
{
	const float	LengthSquared = A * A + B * B + C * C;
	if (LengthSquared > DELTA*DELTA)
	{
		const float	InvLength = Math::InvSqrt(LengthSquared);
		OutPlane = Plane(-A * InvLength, -B * InvLength, -C * InvLength, D * InvLength);
		return 1;
	}
	else
		return 0;
}

FORCEINLINE bool Matrix::GetFrustumNearPlane(Plane& OutPlane) const
{
	return MakeFrustumPlane(
		M[0][2],
		M[1][2],
		M[2][2],
		M[3][2],
		OutPlane
	);
}

FORCEINLINE bool Matrix::GetFrustumFarPlane(Plane& OutPlane) const
{
	return MakeFrustumPlane(
		M[0][3] - M[0][2],
		M[1][3] - M[1][2],
		M[2][3] - M[2][2],
		M[3][3] - M[3][2],
		OutPlane
	);
}

FORCEINLINE bool Matrix::GetFrustumLeftPlane(Plane& OutPlane) const
{
	return MakeFrustumPlane(
		M[0][3] + M[0][0],
		M[1][3] + M[1][0],
		M[2][3] + M[2][0],
		M[3][3] + M[3][0],
		OutPlane
	);
}

FORCEINLINE bool Matrix::GetFrustumRightPlane(Plane& OutPlane) const
{
	return MakeFrustumPlane(
		M[0][3] - M[0][0],
		M[1][3] - M[1][0],
		M[2][3] - M[2][0],
		M[3][3] - M[3][0],
		OutPlane
	);
}

FORCEINLINE bool Matrix::GetFrustumTopPlane(Plane& OutPlane) const
{
	return MakeFrustumPlane(
		M[0][3] - M[0][1],
		M[1][3] - M[1][1],
		M[2][3] - M[2][1],
		M[3][3] - M[3][1],
		OutPlane
	);
}

FORCEINLINE bool Matrix::GetFrustumBottomPlane(Plane& OutPlane) const
{
	return MakeFrustumPlane(
		M[0][3] + M[0][1],
		M[1][3] + M[1][1],
		M[2][3] + M[2][1],
		M[3][3] + M[3][1],
		OutPlane
	);
}


inline void Matrix::Mirror(EAxis::Type MirrorAxis, EAxis::Type FlipAxis)
{
	if (MirrorAxis == EAxis::X)
	{
		M[0][0] *= -1.f;
		M[1][0] *= -1.f;
		M[2][0] *= -1.f;

		M[3][0] *= -1.f;
	}
	else if (MirrorAxis == EAxis::Y)
	{
		M[0][1] *= -1.f;
		M[1][1] *= -1.f;
		M[2][1] *= -1.f;

		M[3][1] *= -1.f;
	}
	else if (MirrorAxis == EAxis::Z)
	{
		M[0][2] *= -1.f;
		M[1][2] *= -1.f;
		M[2][2] *= -1.f;

		M[3][2] *= -1.f;
	}

	if (FlipAxis == EAxis::X)
	{
		M[0][0] *= -1.f;
		M[0][1] *= -1.f;
		M[0][2] *= -1.f;
	}
	else if (FlipAxis == EAxis::Y)
	{
		M[1][0] *= -1.f;
		M[1][1] *= -1.f;
		M[1][2] *= -1.f;
	}
	else if (FlipAxis == EAxis::Z)
	{
		M[2][0] *= -1.f;
		M[2][1] *= -1.f;
		M[2][2] *= -1.f;
	}
}

// very high quality 4x4 matrix inverse
static FORCEINLINE void Inverse4x4(double* dst, const float* src)
{
	const double s0 = (double)(src[0]); const double s1 = (double)(src[1]); const double s2 = (double)(src[2]); const double s3 = (double)(src[3]);
	const double s4 = (double)(src[4]); const double s5 = (double)(src[5]); const double s6 = (double)(src[6]); const double s7 = (double)(src[7]);
	const double s8 = (double)(src[8]); const double s9 = (double)(src[9]); const double s10 = (double)(src[10]); const double s11 = (double)(src[11]);
	const double s12 = (double)(src[12]); const double s13 = (double)(src[13]); const double s14 = (double)(src[14]); const double s15 = (double)(src[15]);

	double inv[16];
	inv[0] = s5 * s10 * s15 - s5 * s11 * s14 - s9 * s6 * s15 + s9 * s7 * s14 + s13 * s6 * s11 - s13 * s7 * s10;
	inv[1] = -s1 * s10 * s15 + s1 * s11 * s14 + s9 * s2 * s15 - s9 * s3 * s14 - s13 * s2 * s11 + s13 * s3 * s10;
	inv[2] = s1 * s6  * s15 - s1 * s7  * s14 - s5 * s2 * s15 + s5 * s3 * s14 + s13 * s2 * s7 - s13 * s3 * s6;
	inv[3] = -s1 * s6  * s11 + s1 * s7  * s10 + s5 * s2 * s11 - s5 * s3 * s10 - s9 * s2 * s7 + s9 * s3 * s6;
	inv[4] = -s4 * s10 * s15 + s4 * s11 * s14 + s8 * s6 * s15 - s8 * s7 * s14 - s12 * s6 * s11 + s12 * s7 * s10;
	inv[5] = s0 * s10 * s15 - s0 * s11 * s14 - s8 * s2 * s15 + s8 * s3 * s14 + s12 * s2 * s11 - s12 * s3 * s10;
	inv[6] = -s0 * s6  * s15 + s0 * s7  * s14 + s4 * s2 * s15 - s4 * s3 * s14 - s12 * s2 * s7 + s12 * s3 * s6;
	inv[7] = s0 * s6  * s11 - s0 * s7  * s10 - s4 * s2 * s11 + s4 * s3 * s10 + s8 * s2 * s7 - s8 * s3 * s6;
	inv[8] = s4 * s9  * s15 - s4 * s11 * s13 - s8 * s5 * s15 + s8 * s7 * s13 + s12 * s5 * s11 - s12 * s7 * s9;
	inv[9] = -s0 * s9  * s15 + s0 * s11 * s13 + s8 * s1 * s15 - s8 * s3 * s13 - s12 * s1 * s11 + s12 * s3 * s9;
	inv[10] = s0 * s5  * s15 - s0 * s7  * s13 - s4 * s1 * s15 + s4 * s3 * s13 + s12 * s1 * s7 - s12 * s3 * s5;
	inv[11] = -s0 * s5  * s11 + s0 * s7  * s9 + s4 * s1 * s11 - s4 * s3 * s9 - s8 * s1 * s7 + s8 * s3 * s5;
	inv[12] = -s4 * s9  * s14 + s4 * s10 * s13 + s8 * s5 * s14 - s8 * s6 * s13 - s12 * s5 * s10 + s12 * s6 * s9;
	inv[13] = s0 * s9  * s14 - s0 * s10 * s13 - s8 * s1 * s14 + s8 * s2 * s13 + s12 * s1 * s10 - s12 * s2 * s9;
	inv[14] = -s0 * s5  * s14 + s0 * s6  * s13 + s4 * s1 * s14 - s4 * s2 * s13 - s12 * s1 * s6 + s12 * s2 * s5;
	inv[15] = s0 * s5  * s10 - s0 * s6  * s9 - s4 * s1 * s10 + s4 * s2 * s9 + s8 * s1 * s6 - s8 * s2 * s5;

	double det = s0 * inv[0] + s1 * inv[4] + s2 * inv[8] + s3 * inv[12];
	if (det != 0.0)
	{
		det = 1.0 / det;
	}
	for (int i = 0; i < 16; i++)
	{
		dst[i] = inv[i] * det;
	}
}

FORCEINLINE Vector4f Matrix::TransformVector4f(const Vector4f &P) const
{
	Vector4f Result;
	VectorRegister VecP = VectorLoadAligned(&P);
	VectorRegister VecR = VectorTransformVector(VecP, this);
	VectorStoreAligned(VecR, &Result);
	return Result;
}

FORCEINLINE Vector4f Matrix::TransformPosition(const Vector3f &V) const
{
	return TransformVector4f(Vector4f(V.X, V.Y, V.Z, 1.0f));
}

FORCEINLINE Vector3f Matrix::InverseTransformPosition(const Vector3f &V) const
{
	Matrix InvSelf = this->InverseFast();
	return InvSelf.TransformPosition(V);
}

FORCEINLINE Vector4f Matrix::TransformVector(const Vector3f& V) const
{
	return TransformVector4f(Vector4f(V.X, V.Y, V.Z, 0.0f));
}

FORCEINLINE Vector3f Matrix::InverseTransformVector(const Vector3f &V) const
{
	Matrix InvSelf = this->InverseFast();
	return InvSelf.TransformVector(V);
}

};
