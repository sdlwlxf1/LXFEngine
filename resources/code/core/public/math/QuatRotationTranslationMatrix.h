// this tiny game engine is rendered by CPU. I make it for learing the pipline and rasterization. This is the first step to realise my own game engine. so just do it!

#pragma once

#include "Matrix.h"

namespace LXF {

/** Rotation and translation matrix using quaternion rotation */
class QuatRotationTranslationMatrix
	: public Matrix
{
public:

	/** Constructor
	*
	* @param Q rotation
	* @param Origin translation to apply
	*/
	QuatRotationTranslationMatrix(const Quat& Q, const Vector3f& Origin);

	/** Matrix factory. Return an FMatrix so we don't have type conversion issues in expressions. */
	static Matrix Make(const Quat& Q, const Vector3f& Origin)
	{
		return QuatRotationTranslationMatrix(Q, Origin);
	}
};


/** Rotation matrix using quaternion rotation */
class QuatRotationMatrix
	: public QuatRotationTranslationMatrix
{
public:

	/** Constructor
	*
	* @param Q rotation
	*/
	QuatRotationMatrix(const Quat& Q)
		: QuatRotationTranslationMatrix(Q, Vector3f::ZeroVector)
	{
	}

	/** Matrix factory. Return an FMatrix so we don't have type conversion issues in expressions. */
	static Matrix Make(const Quat& Q)
	{
		return QuatRotationMatrix(Q);
}
};


FORCEINLINE QuatRotationTranslationMatrix::QuatRotationTranslationMatrix(const Quat& Q, const Vector3f& Origin)
{
#if !(BUILD_SHIPPING)
	// Make sure Quaternion is normalized
	check(Q.IsNormalized());
#endif
	const float x2 = Q.X + Q.X;  const float y2 = Q.Y + Q.Y;  const float z2 = Q.Z + Q.Z;
	const float xx = Q.X * x2;   const float xy = Q.X * y2;   const float xz = Q.X * z2;
	const float yy = Q.Y * y2;   const float yz = Q.Y * z2;   const float zz = Q.Z * z2;
	const float wx = Q.W * x2;   const float wy = Q.W * y2;   const float wz = Q.W * z2;

	M[0][0] = 1.0f - (yy + zz);	M[1][0] = xy - wz;				M[2][0] = xz + wy;			M[3][0] = Origin.X;
	M[0][1] = xy + wz;			M[1][1] = 1.0f - (xx + zz);		M[2][1] = yz - wx;			M[3][1] = Origin.Y;
	M[0][2] = xz - wy;			M[1][2] = yz + wx;				M[2][2] = 1.0f - (xx + yy);	M[3][2] = Origin.Z;
	M[0][3] = 0.0f;				M[1][3] = 0.0f;					M[2][3] = 0.0f;				M[3][3] = 1.0f;
}
};
