// this tiny game engine is rendered by CPU. I make it for learing the pipline and rasterization. This is the first step to realise my own game engine. so just do it!
#include "math/MathUtility.h"
#include "math/Vector.h"
#include "math/Vector2.h"
#include "math/Vector4.h"
#include "math/Plane.h"
#include "math/Rotator.h"
#include "math/Quat.h"
#include "math/Matrix.h"
#include "math/RotationMatrix.h"
#include "math/QuatRotationTranslationMatrix.h"

using namespace LXF;
const Vector3f Vector3f::ZeroVector(0.0f, 0.0f, 0.0f);
const Vector3f Vector3f::OneVector(1.0f, 1.0f, 1.0f);
const Vector3f Vector3f::UpVector(0.0f, 0.0f, 1.0f);
const Vector3f Vector3f::DownVector(0.0f, 0.0f, -1.0f);
const Vector3f Vector3f::ForwardVector(1.0f, 0.0f, 0.0f);
const Vector3f Vector3f::BackwardVector(-1.0f, 0.0f, 0.0f);
const Vector3f Vector3f::RightVector(0.0f, 1.0f, 0.0f);
const Vector3f Vector3f::LeftVector(0.0f, -1.0f, 0.0f);

const Vector3i Vector3i::ZeroVector(0, 0, 0);
const Vector3i Vector3i::OneVector(1, 1, 1);
const Vector3i Vector3i::UpVector(0, 0, 1);
const Vector3i Vector3i::DownVector(0, 0, -1);
const Vector3i Vector3i::ForwardVector(1, 0, 0);
const Vector3i Vector3i::BackwardVector(-1, 0, 0);
const Vector3i Vector3i::RightVector(0, 1, 0);
const Vector3i Vector3i::LeftVector(0, -1, 0);

const Vector2f Vector2f::ZeroVector(0.0f, 0.0f);
const Vector2f Vector2f::OneVector(1.0f, 1.0f);
const Vector2f Vector2f::ForwardVector(1.0f, 0.0f);
const Vector2f Vector2f::BackwardVector(-1.0f, 0.0f);
const Vector2f Vector2f::RightVector(0.0f, 1.0f);
const Vector2f Vector2f::LeftVector(0.0f, -1.0f);

const Vector2i Vector2i::ZeroVector(0, 0);
const Vector2i Vector2i::OneVector(1, 1);
const Vector2i Vector2i::ForwardVector(1, 0);
const Vector2i Vector2i::BackwardVector(-1, 0);
const Vector2i Vector2i::RightVector(0, 1);
const Vector2i Vector2i::LeftVector(0, -1);

const Vector4f Vector4f::ZeroVector(0.0f, 0.0f, 0.0f);
const Vector4f Vector4f::OneVector(1.0f, 1.0f, 1.0f);
const Vector4f Vector4f::UpVector(0.0f, 0.0f, 1.0f);
const Vector4f Vector4f::DownVector(0.0f, 0.0f, -1.0f);
const Vector4f Vector4f::ForwardVector(1.0f, 0.0f, 0.0f);
const Vector4f Vector4f::BackwardVector(-1.0f, 0.0f, 0.0f);
const Vector4f Vector4f::RightVector(0.0f, 1.0f, 0.0f);
const Vector4f Vector4f::LeftVector(0.0f, -1.0f, 0.0f);

const Vector4i Vector4i::ZeroVector(0, 0, 0);
const Vector4i Vector4i::OneVector(1, 1, 1);
const Vector4i Vector4i::UpVector(0, 0, 1);
const Vector4i Vector4i::DownVector(0, 0, -1);
const Vector4i Vector4i::ForwardVector(1, 0, 0);
const Vector4i Vector4i::BackwardVector(-1, 0, 0);
const Vector4i Vector4i::RightVector(0, 1, 0);
const Vector4i Vector4i::LeftVector(0, -1, 0);

const Rotator Rotator::ZeroRotator(0.f, 0.f, 0.f);

const Matrix Matrix::Identity(Plane(1, 0, 0, 0), Plane(0, 1, 0, 0), Plane(0, 0, 1, 0), Plane(0, 0, 0, 1));

const Quat Quat::Identity(0, 0, 0, 1);

template <typename T>
Vector2<T>::Vector2(const Vector4<T>& V) : X(V.X) , Y(V.Y) { CheckNaN(); }

template <typename T>
Vector2<T>::Vector2(const Vector<T>& V) : X(V.X) , Y(V.Y) { CheckNaN(); }

template <typename T>
Vector<T>::Vector(const Vector2<T> V, float InZ) : X(V.X), Y(V.Y), Z(InZ) { CheckNaN(); }

template <typename T>
Vector<T>::Vector(const Vector4<T>& V) : X(V.X), Y(V.Y), Z(V.Z) { CheckNaN(); }

template <typename T>
Vector4<T>::Vector4(const Vector<T>& InVector, T InW) : X(InVector.X), Y(InVector.Y), Z(InVector.Z), W(InW) { CheckNaN(); }


Matrix Quat::operator*(const Matrix& M) const
{
	Matrix Result;
	Quat VT, VR;
	Quat Inv = Inverse();
	for (int32 I = 0; I < 4; ++I)
	{
		Quat VQ(M.M[I][0], M.M[I][1], M.M[I][2], M.M[I][3]);
		VectorQuaternionMultiply(&VT, this, &VQ);
		VectorQuaternionMultiply(&VR, &VT, &Inv);
		Result.M[I][0] = VR.X;
		Result.M[I][1] = VR.Y;
		Result.M[I][2] = VR.Z;
		Result.M[I][3] = VR.W;
	}

	return Result;
}


Vector3f Rotator::Vector() const
{
	const float PitchNoWinding = Math::Fmod(Pitch, 360.f);
	const float YawNoWinding = Math::Fmod(Yaw, 360.f);

	float CP, SP, CY, SY;
	Math::SinCos(&SP, &CP, Math::DegreesToRadians(PitchNoWinding));
	Math::SinCos(&SY, &CY, Math::DegreesToRadians(YawNoWinding));
	Vector3f V = Vector3f(CP*CY, CP*SY, SP);

#if ENABLE_NAN_DIAGNOSTIC
	if (V.ContainsNaN())
	{
		V = Vector3f::ForwardVector;
	}
#endif

	return V;
}

Quat Rotator::Quaternion() const
{
	//SCOPE_CYCLE_COUNTER(STAT_MathConvertRotatorToQuat);

	CheckNaN();

#if PLATFORM_ENABLE_VECTORINTRINSICS
	const VectorRegister Angles = MakeVectorRegister(Pitch, Yaw, Roll, 0.0f);
	const VectorRegister AnglesNoWinding = VectorMod(Angles, GlobalVectorConstants::Float360);
	const VectorRegister HalfAngles = VectorMultiply(AnglesNoWinding, GlobalVectorConstants::DEG_TO_RAD_HALF);

	VectorRegister SinAngles, CosAngles;
	VectorSinCos(&SinAngles, &CosAngles, &HalfAngles);

	// Vectorized conversion, measured 20% faster than using scalar version after VectorSinCos.
	// Indices within VectorRegister (for shuffles): P=0, Y=1, R=2
	const VectorRegister SR = VectorReplicate(SinAngles, 2);
	const VectorRegister CR = VectorReplicate(CosAngles, 2);

	const VectorRegister SY_SY_CY_CY_Temp = VectorShuffle(SinAngles, CosAngles, 1, 1, 1, 1);

	const VectorRegister SP_SP_CP_CP = VectorShuffle(SinAngles, CosAngles, 0, 0, 0, 0);
	const VectorRegister SY_CY_SY_CY = VectorShuffle(SY_SY_CY_CY_Temp, SY_SY_CY_CY_Temp, 0, 2, 0, 2);

	const VectorRegister CP_CP_SP_SP = VectorShuffle(CosAngles, SinAngles, 0, 0, 0, 0);
	const VectorRegister CY_SY_CY_SY = VectorShuffle(SY_SY_CY_CY_Temp, SY_SY_CY_CY_Temp, 2, 0, 2, 0);

	const uint32 Neg = uint32(1 << 31);
	const uint32 Pos = uint32(0);
	const VectorRegister SignBitsLeft = MakeVectorRegister(Pos, Neg, Pos, Pos);
	const VectorRegister SignBitsRight = MakeVectorRegister(Neg, Neg, Neg, Pos);
	const VectorRegister LeftTerm = VectorBitwiseXor(SignBitsLeft, VectorMultiply(CR, VectorMultiply(SP_SP_CP_CP, SY_CY_SY_CY)));
	const VectorRegister RightTerm = VectorBitwiseXor(SignBitsRight, VectorMultiply(SR, VectorMultiply(CP_CP_SP_SP, CY_SY_CY_SY)));

	Quat RotationQuat;
	const VectorRegister Result = VectorAdd(LeftTerm, RightTerm);
	VectorStoreAligned(Result, &RotationQuat);
#else
	const float DEG_TO_RAD = PI / (180.f);
	const float RADS_DIVIDED_BY_2 = DEG_TO_RAD / 2.f;
	float SP, SY, SR;
	float CP, CY, CR;

	const float PitchNoWinding = Math::Fmod(Pitch, 360.0f);
	const float YawNoWinding = Math::Fmod(Yaw, 360.0f);
	const float RollNoWinding = Math::Fmod(Roll, 360.0f);

	Math::SinCos(&SP, &CP, PitchNoWinding * RADS_DIVIDED_BY_2);
	Math::SinCos(&SY, &CY, YawNoWinding * RADS_DIVIDED_BY_2);
	Math::SinCos(&SR, &CR, RollNoWinding * RADS_DIVIDED_BY_2);

	Quat RotationQuat;
	RotationQuat.X = CR * SP*SY - SR * CP*CY;
	RotationQuat.Y = -CR * SP*CY - SR * CP*SY;
	RotationQuat.Z = CR * CP*SY - SR * SP*CY;
	RotationQuat.W = CR * CP*CY + SR * SP*SY;
#endif // PLATFORM_ENABLE_VECTORINTRINSICS

#if ENABLE_NAN_DIAGNOSTIC
	// Very large inputs can cause NaN's. Want to catch this here
	if (RotationQuat.ContainsNaN())
	{
		RotationQuat = Quat::Identity;
	}
#endif

	return RotationQuat;
}

class Rotator Quat::Rotator() const
{
	CheckNaN();
	const float SingularityTest = Z * X - W * Y;
	const float YawY = 2.f*(W*Z + X * Y);
	const float YawX = (1.f - 2.f*(Math::Square(Y) + Math::Square(Z)));

	// reference 
	// http://en.wikipedia.org/wiki/Conversion_between_quaternions_and_Euler_angles
	// http://www.euclideanspace.com/maths/geometry/rotations/conversions/quaternionToEuler/

	// this value was found from experience, the above websites recommend different values
	// but that isn't the case for us, so I went through different testing, and finally found the case 
	// where both of world lives happily. 
	const float SINGULARITY_THRESHOLD = 0.4999995f;
	const float RAD_TO_DEG = (180.f) / PI;
	class Rotator RotatorFromQuat;

	if (SingularityTest < -SINGULARITY_THRESHOLD)
	{
		RotatorFromQuat.Pitch = -90.f;
		RotatorFromQuat.Yaw = Math::Atan2(YawY, YawX) * RAD_TO_DEG;
		RotatorFromQuat.Roll = Rotator::NormalizeAxis(-RotatorFromQuat.Yaw - (2.f * Math::Atan2(X, W) * RAD_TO_DEG));
	}
	else if (SingularityTest > SINGULARITY_THRESHOLD)
	{
		RotatorFromQuat.Pitch = 90.f;
		RotatorFromQuat.Yaw = Math::Atan2(YawY, YawX) * RAD_TO_DEG;
		RotatorFromQuat.Roll = Rotator::NormalizeAxis(RotatorFromQuat.Yaw - (2.f * Math::Atan2(X, W) * RAD_TO_DEG));
	}
	else
	{
		RotatorFromQuat.Pitch = Math::FastAsin(2.f*(SingularityTest)) * RAD_TO_DEG;
		RotatorFromQuat.Yaw = Math::Atan2(YawY, YawX) * RAD_TO_DEG;
		RotatorFromQuat.Roll = Math::Atan2(-2.f*(W*X + Y * Z), (1.f - 2.f*(Math::Square(X) + Math::Square(Y)))) * RAD_TO_DEG;
	}

#if ENABLE_NAN_DIAGNOSTIC
	if (RotatorFromQuat.ContainsNaN())
	{
		logOrEnsureNanError(TEXT("Quat::Rotator(): Rotator result %s contains NaN! Quat = %s, YawY = %.9f, YawX = %.9f"), *RotatorFromQuat.ToString(), *this->ToString(), YawY, YawX);
		RotatorFromQuat = Rotator::ZeroRotator;
	}
#endif

	return RotatorFromQuat;
}

Quat::Quat(const Matrix& M)
{
	// If Matrix is NULL, return Identity quaternion. If any of them is 0, you won't be able to construct rotation
	// if you have two plane at least, we can reconstruct the frame using cross product, but that's a bit expensive op to do here
	// for now, if you convert to matrix from 0 scale and convert back, you'll lose rotation. Don't do that. 
	if (M.GetScaledAxis(EAxis::X).IsNearlyZero() || M.GetScaledAxis(EAxis::Y).IsNearlyZero() || M.GetScaledAxis(EAxis::Z).IsNearlyZero())
	{
		*this = Quat::Identity;
		return;
	}

#if !(BUILD_SHIPPING)
	// Make sure the Rotation part of the Matrix is unit length.
	// Changed to this (same as RemoveScaling) from RotDeterminant as using two different ways of checking unit length matrix caused inconsistency. 
	if (!((Math::Abs(1.f - M.GetScaledAxis(EAxis::X).SizeSquared()) <= KINDA_SMALL_NUMBER) && (Math::Abs(1.f - M.GetScaledAxis(EAxis::Y).SizeSquared()) <= KINDA_SMALL_NUMBER) && (Math::Abs(1.f - M.GetScaledAxis(EAxis::Z).SizeSquared()) <= KINDA_SMALL_NUMBER)))
	{
		*this = Quat::Identity;
		return;
	}
#endif

	//const MeReal *const t = (MeReal *) tm;
	float	s;

	// Check diagonal (trace)
	const float tr = M.M[0][0] + M.M[1][1] + M.M[2][2];

	if (tr > 0.0f)
	{
		float InvS = Math::InvSqrt(tr + 1.f);
		this->W = 0.5f * (1.f / InvS);
		s = 0.5f * InvS;

		this->X = (M.M[1][2] - M.M[2][1]) * s;
		this->Y = (M.M[2][0] - M.M[0][2]) * s;
		this->Z = (M.M[0][1] - M.M[1][0]) * s;
	}
	else
	{
		// diagonal is negative
		int32 i = 0;

		if (M.M[1][1] > M.M[0][0])
			i = 1;

		if (M.M[2][2] > M.M[i][i])
			i = 2;

		static const int32 nxt[3] = { 1, 2, 0 };
		const int32 j = nxt[i];
		const int32 k = nxt[j];

		s = M.M[i][i] - M.M[j][j] - M.M[k][k] + 1.0f;

		float InvS = Math::InvSqrt(s);

		float qt[4];
		qt[i] = 0.5f * (1.f / InvS);

		s = 0.5f * InvS;

		qt[3] = (M.M[j][k] - M.M[k][j]) * s;
		qt[j] = (M.M[i][j] + M.M[j][i]) * s;
		qt[k] = (M.M[i][k] + M.M[k][i]) * s;

		this->X = qt[0];
		this->Y = qt[1];
		this->Z = qt[2];
		this->W = qt[3];

		CheckNaN();
	}
}

Quat::Quat(const class Rotator& R)
{
	*this = R.Quaternion();
	CheckNaN();
}

Quat::Quat(Vector3f Axis, float AngleRad)
{
	const float half_a = 0.5f * AngleRad;
	float s, c;
	Math::SinCos(&s, &c, half_a);

	X = s * Axis.X;
	Y = s * Axis.Y;
	Z = s * Axis.Z;
	W = c;

	CheckNaN();
}


Vector3f Rotator::RotateVector(const Vector3f& V) const
{
	return RotationMatrix(*this).TransformVector(V);
}

Vector3f Rotator::UnRotateVector(const Vector3f& V) const
{
	return RotationMatrix(*this).GetTransposed().TransformVector(V);
}

Matrix RotationMatrix::Make(Quat const& Rot)
{
	return QuatRotationTranslationMatrix(Rot, Vector3f::ZeroVector);
}

Matrix RotationMatrix::MakeFromX(Vector3f const& XAxis)
{
	Vector3f const NewX = XAxis.GetSafeNormal();

	// try to use up if possible
	Vector3f const UpVector = (Math::Abs(NewX.Z) < (1.f - KINDA_SMALL_NUMBER)) ? Vector3f(0, 0, 1.f) : Vector3f(1.f, 0, 0);

	const Vector3f NewY = (UpVector ^ NewX).GetSafeNormal();
	const Vector3f NewZ = NewX ^ NewY;

	return Matrix(NewX, NewY, NewZ, Vector3f::ZeroVector);
}

Matrix RotationMatrix::MakeFromY(Vector3f const& YAxis)
{
	Vector3f const NewY = YAxis.GetSafeNormal();

	// try to use up if possible
	Vector3f const UpVector = (Math::Abs(NewY.Z) < (1.f - KINDA_SMALL_NUMBER)) ? Vector3f(0, 0, 1.f) : Vector3f(1.f, 0, 0);

	const Vector3f NewZ = (UpVector ^ NewY).GetSafeNormal();
	const Vector3f NewX = NewY ^ NewZ;

	return Matrix(NewX, NewY, NewZ, Vector3f::ZeroVector);
}

Matrix RotationMatrix::MakeFromZ(Vector3f const& ZAxis)
{
	Vector3f const NewZ = ZAxis.GetSafeNormal();

	// try to use up if possible
	Vector3f const UpVector = (Math::Abs(NewZ.Z) < (1.f - KINDA_SMALL_NUMBER)) ? Vector3f(0, 0, 1.f) : Vector3f(1.f, 0, 0);

	const Vector3f NewX = (UpVector ^ NewZ).GetSafeNormal();
	const Vector3f NewY = NewZ ^ NewX;

	return Matrix(NewX, NewY, NewZ, Vector3f::ZeroVector);
}

Matrix RotationMatrix::MakeFromXY(Vector3f const& XAxis, Vector3f const& YAxis)
{
	Vector3f NewX = XAxis.GetSafeNormal();
	Vector3f Norm = YAxis.GetSafeNormal();

	// if they're almost same, we need to find arbitrary vector
	if (Math::IsNearlyEqual(Math::Abs(NewX | Norm), 1.f))
	{
		// make sure we don't ever pick the same as NewX
		Norm = (Math::Abs(NewX.Z) < (1.f - KINDA_SMALL_NUMBER)) ? Vector3f(0, 0, 1.f) : Vector3f(1.f, 0, 0);
	}

	const Vector3f NewZ = (NewX ^ Norm).GetSafeNormal();
	const Vector3f NewY = NewZ ^ NewX;

	return Matrix(NewX, NewY, NewZ, Vector3f::ZeroVector);
}

Matrix RotationMatrix::MakeFromXZ(Vector3f const& XAxis, Vector3f const& ZAxis)
{
	Vector3f const NewX = XAxis.GetSafeNormal();
	Vector3f Norm = ZAxis.GetSafeNormal();

	// if they're almost same, we need to find arbitrary vector
	if (Math::IsNearlyEqual(Math::Abs(NewX | Norm), 1.f))
	{
		// make sure we don't ever pick the same as NewX
		Norm = (Math::Abs(NewX.Z) < (1.f - KINDA_SMALL_NUMBER)) ? Vector3f(0, 0, 1.f) : Vector3f(1.f, 0, 0);
	}

	const Vector3f NewY = (Norm ^ NewX).GetSafeNormal();
	const Vector3f NewZ = NewX ^ NewY;

	return Matrix(NewX, NewY, NewZ, Vector3f::ZeroVector);
}

Matrix RotationMatrix::MakeFromYX(Vector3f const& YAxis, Vector3f const& XAxis)
{
	Vector3f const NewY = YAxis.GetSafeNormal();
	Vector3f Norm = XAxis.GetSafeNormal();

	// if they're almost same, we need to find arbitrary vector
	if (Math::IsNearlyEqual(Math::Abs(NewY | Norm), 1.f))
	{
		// make sure we don't ever pick the same as NewX
		Norm = (Math::Abs(NewY.Z) < (1.f - KINDA_SMALL_NUMBER)) ? Vector3f(0, 0, 1.f) : Vector3f(1.f, 0, 0);
	}

	const Vector3f NewZ = (Norm ^ NewY).GetSafeNormal();
	const Vector3f NewX = NewY ^ NewZ;

	return Matrix(NewX, NewY, NewZ, Vector3f::ZeroVector);
}

Matrix RotationMatrix::MakeFromYZ(Vector3f const& YAxis, Vector3f const& ZAxis)
{
	Vector3f const NewY = YAxis.GetSafeNormal();
	Vector3f Norm = ZAxis.GetSafeNormal();

	// if they're almost same, we need to find arbitrary vector
	if (Math::IsNearlyEqual(Math::Abs(NewY | Norm), 1.f))
	{
		// make sure we don't ever pick the same as NewX
		Norm = (Math::Abs(NewY.Z) < (1.f - KINDA_SMALL_NUMBER)) ? Vector3f(0, 0, 1.f) : Vector3f(1.f, 0, 0);
	}

	const Vector3f NewX = (NewY ^ Norm).GetSafeNormal();
	const Vector3f NewZ = NewX ^ NewY;

	return Matrix(NewX, NewY, NewZ, Vector3f::ZeroVector);
}

Matrix RotationMatrix::MakeFromZX(Vector3f const& ZAxis, Vector3f const& XAxis)
{
	Vector3f const NewZ = ZAxis.GetSafeNormal();
	Vector3f Norm = XAxis.GetSafeNormal();

	// if they're almost same, we need to find arbitrary vector
	if (Math::IsNearlyEqual(Math::Abs(NewZ | Norm), 1.f))
	{
		// make sure we don't ever pick the same as NewX
		Norm = (Math::Abs(NewZ.Z) < (1.f - KINDA_SMALL_NUMBER)) ? Vector3f(0, 0, 1.f) : Vector3f(1.f, 0, 0);
	}

	const Vector3f NewY = (NewZ ^ Norm).GetSafeNormal();
	const Vector3f NewX = NewY ^ NewZ;

	return Matrix(NewX, NewY, NewZ, Vector3f::ZeroVector);
}

Matrix RotationMatrix::MakeFromZY(Vector3f const& ZAxis, Vector3f const& YAxis)
{
	Vector3f const NewZ = ZAxis.GetSafeNormal();
	Vector3f Norm = YAxis.GetSafeNormal();

	// if they're almost same, we need to find arbitrary vector
	if (Math::IsNearlyEqual(Math::Abs(NewZ | Norm), 1.f))
	{
		// make sure we don't ever pick the same as NewX
		Norm = (Math::Abs(NewZ.Z) < (1.f - KINDA_SMALL_NUMBER)) ? Vector3f(0, 0, 1.f) : Vector3f(1.f, 0, 0);
	}

	const Vector3f NewX = (Norm ^ NewZ).GetSafeNormal();
	const Vector3f NewY = NewZ ^ NewX;

	return Matrix(NewX, NewY, NewZ, Vector3f::ZeroVector);
}

float Math::SRand()
{
	return 0;
}

Vector3f Math::RayPlaneIntersection(const Vector3f& RayOrigin, const Vector3f& RayDirection, const Plane& Plane)
{
	const Vector3f PlaneNormal = Vector3f(Plane.X, Plane.Y, Plane.Z);
	const Vector3f PlaneOrigin = PlaneNormal * Plane.W;

	const float Distance = Vector3f::DotProduct((PlaneOrigin - RayOrigin), PlaneNormal) / Vector3f::DotProduct(RayDirection, PlaneNormal);
	return RayOrigin + RayDirection * Distance;
}

Vector3f Math::LinePlaneIntersection(const Vector3f &Point1, const Vector3f &Point2, const Plane &Plane)
{
	return
		Point1
		+ (Point2 - Point1)
		*	((Plane.W - (Point1 | Plane)) / ((Point2 - Point1) | Plane));
}

Vector3f Math::LinePlaneIntersection(const Vector3f& Point1, const Vector3f& Point2, const Vector3f& PlaneOrigin, const Vector3f& PlaneNormal)
{
	return
		Point1
		+ (Point2 - Point1)
		*	(((PlaneOrigin - Point1) | PlaneNormal) / ((Point2 - Point1) | PlaneNormal));
}

float Math::Atan2(float Y, float X)
{
	//return atan2f(Y,X);
	// atan2f occasionally returns NaN with perfectly valid input (possibly due to a compiler or library bug).
	// We are replacing it with a minimax approximation with a max relative error of 7.15255737e-007 compared to the C library function.
	// On PC this has been measured to be 2x faster than the std C version.

	const float absX = Math::Abs(X);
	const float absY = Math::Abs(Y);
	const bool yAbsBigger = (absY > absX);
	float t0 = yAbsBigger ? absY : absX; // Max(absY, absX)
	float t1 = yAbsBigger ? absX : absY; // Min(absX, absY)

	if (t0 == 0.f)
		return 0.f;

	float t3 = t1 / t0;
	float t4 = t3 * t3;

	static const float c[7] = {
		+7.2128853633444123e-03f,
		-3.5059680836411644e-02f,
		+8.1675882859940430e-02f,
		-1.3374657325451267e-01f,
		+1.9856563505717162e-01f,
		-3.3324998579202170e-01f,
		+1.0f
	};

	t0 = c[0];
	t0 = t0 * t4 + c[1];
	t0 = t0 * t4 + c[2];
	t0 = t0 * t4 + c[3];
	t0 = t0 * t4 + c[4];
	t0 = t0 * t4 + c[5];
	t0 = t0 * t4 + c[6];
	t3 = t0 * t3;

	t3 = yAbsBigger ? (0.5f * PI) - t3 : t3;
	t3 = (X < 0.0f) ? PI - t3 : t3;
	t3 = (Y < 0.0f) ? -t3 : t3;

	return t3;
}
