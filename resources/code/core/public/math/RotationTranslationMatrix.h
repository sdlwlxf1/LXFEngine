// this tiny game engine is rendered by CPU. I make it for learing the pipline and rasterization. This is the first step to realise my own game engine. so just do it!

#pragma once

#include "Matrix.h"
#include "Rotator.h"

namespace LXF {

/** Combined rotation and translation matrix */
class RotationTranslationMatrix
	: public Matrix
{
public:

	RotationTranslationMatrix(const Rotator& Rot, const Vector3f& Origin);

	static Matrix Make(const Rotator& Rot, const Vector3f& Origin)
	{
		return RotationTranslationMatrix(Rot, Origin);
	}
};


FORCEINLINE RotationTranslationMatrix::RotationTranslationMatrix(const Rotator& Rot, const Vector3f& Origin)
{
#if PLATFORM_ENABLE_VECTORINTRINSICS

	const VectorRegister Angles = MakeVectorRegister(Rot.Pitch, Rot.Yaw, Rot.Roll, 0.0f);
	const VectorRegister HalfAngles = VectorMultiply(Angles, GlobalVectorConstants::DEG_TO_RAD);

	union { VectorRegister v; float f[4]; } SinAngles, CosAngles;
	VectorSinCos(&SinAngles.v, &CosAngles.v, &HalfAngles);

	const float	SP = SinAngles.f[0];
	const float	SY = SinAngles.f[1];
	const float	SR = SinAngles.f[2];
	const float	CP = CosAngles.f[0];
	const float	CY = CosAngles.f[1];
	const float	CR = CosAngles.f[2];

#else

	float SP, SY, SR;
	float CP, CY, CR;
	FMath::SinCos(&SP, &CP, FMath::DegreesToRadians(Rot.Pitch));
	FMath::SinCos(&SY, &CY, FMath::DegreesToRadians(Rot.Yaw));
	FMath::SinCos(&SR, &CR, FMath::DegreesToRadians(Rot.Roll));

#endif // PLATFORM_ENABLE_VECTORINTRINSICS

	M[0][0] = CP * CY;
	M[0][1] = CP * SY;
	M[0][2] = SP;
	M[0][3] = 0.f;

	M[1][0] = SR * SP * CY - CR * SY;
	M[1][1] = SR * SP * SY + CR * CY;
	M[1][2] = -SR * CP;
	M[1][3] = 0.f;

	M[2][0] = -(CR * SP * CY + SR * SY);
	M[2][1] = CY * SR - CR * SP * SY;
	M[2][2] = CR * CP;
	M[2][3] = 0.f;

	M[3][0] = Origin.X;
	M[3][1] = Origin.Y;
	M[3][2] = Origin.Z;
	M[3][3] = 1.f;
}

};
