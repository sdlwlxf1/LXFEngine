// this tiny game engine is rendered by CPU. I make it for learing the pipline and rasterization. floathis is the first step to realise my own game engine. so just do it!

#pragma once

#include "CoreTypes.h"
#include "MathUtility.h"
#include "VectorRegister.h"

namespace LXF {

class Quat
{
public:
	float X, Y, Z, W;

	static const Quat Identity;

#if ENABLE_NAN_DIAGNOSTIC
	void CheckNaN() const {
		if (ContainsNaN())
		{
			*const_cast<Quat*>(this) = Identity;
		}
	}
#else
	void CheckNaN() const {}
#endif	

	FORCEINLINE Quat() : X(0), Y(0), Z(0), W(1) { };
	FORCEINLINE Quat(float InF) : X(InF), Y(InF), Z(InF) { CheckNaN(); }
	FORCEINLINE Quat(float InX, float InY, float InZ, float InW) : X(InX), Y(InY), Z(InZ), W(InW) { CheckNaN(); }
	explicit Quat(const Matrix& M);
	explicit Quat(const class Rotator& R);

	Quat(Vector3f Axis, float AngleRad);

	FORCEINLINE Quat operator+(const Quat& Q) const { return Quat(X + Q.X, Y + Q.Y, Z + Q.Z, W + Q.W); }
	FORCEINLINE Quat operator-(const Quat& Q) const { return Quat(X - Q.X, Y - Q.Y, Z - Q.Z, W - Q.W); }
	FORCEINLINE Quat operator*(float Scale) const { return Quat(X * Scale, Y * Scale, Z * Scale, W * Scale); }
	FORCEINLINE Quat operator/(float Scale) const { const float RScale = 1.0f / Scale; return Quat(X * RScale, Y * RScale, Z * RScale, W * RScale); }
	FORCEINLINE bool operator==(const Quat& Q) const { return X == Q.X && Y == Q.Y && Z == Q.Z && W == Q.W; }
	FORCEINLINE bool operator!=(const Quat& Q) const { return X != Q.X || Y != Q.Y || Z != Q.Z || W != Q.W; }
	FORCEINLINE Quat operator-() const { return Quat(-X, -Y, -Z, -W); }

	FORCEINLINE Quat operator+=(const Quat& Q) { X += Q.X; Y += Q.Y; Z += Q.Z; W += Q.W; CheckNaN(); return *this; }
	FORCEINLINE Quat operator-=(const Quat& Q) { X -= Q.X; Y -= Q.Y; Z -= Q.Z; W -= Q.W; CheckNaN(); return *this; }
	FORCEINLINE Quat operator*=(float Scale) { X *= Scale; Y *= Scale; Z *= Scale; W *= Scale; CheckNaN(); return *this; }
	FORCEINLINE Quat operator/=(float Scale) { const float RScale = 1.0f / Scale; X *= RScale; Y *= RScale; Z *= RScale; W *= RScale; CheckNaN(); return *this; }
	FORCEINLINE Quat operator/=(const Quat& Q) { X /= Q.X; Y /= Q.Y; Z /= Q.Z; W /= Q.W; CheckNaN(); return *this; }
	FORCEINLINE Quat operator*(const Quat& Q) const;
	FORCEINLINE Quat operator*=(const Quat& Q);

	FORCEINLINE Vector3f operator*(const Vector3f& V) const;
	FORCEINLINE Matrix operator*(const Matrix& M) const;
	FORCEINLINE Quat Inverse() const;

	FORCEINLINE Quat GetNormalized() const { Quat Rot = *this; Rot.Normalize(); return Rot; }

	FORCEINLINE bool Equals(const Quat& Q, float Tolerance = KINDA_SMALL_NUMBER) const;

	FORCEINLINE bool ContainsNaN() const { return (!Math::IsFinite(X) || !Math::IsFinite(Y) || !Math::IsFinite(Z) || !Math::IsFinite(W)); }
	FORCEINLINE void Normalize(float Tolerance = SMALL_NUMBER);

	FORCEINLINE Quat GetNormalized(float Tolerance = SMALL_NUMBER) const;
	bool IsNormalized() const;
	FORCEINLINE float Size() const;
	FORCEINLINE float SizeSquared() const;
	FORCEINLINE float GetAngle() const;

	FORCEINLINE Vector3f GetAxisX() const;
	FORCEINLINE Vector3f GetAxisY() const;
	FORCEINLINE Vector3f GetAxisZ() const;
	FORCEINLINE Vector3f GetForwardVector() const;
	FORCEINLINE Vector3f GetRightVector() const;
	FORCEINLINE Vector3f GetUpVector() const;
	FORCEINLINE Vector3f Vector() const;
	class Rotator Rotator() const;

	FORCEINLINE Vector3f RotateVector(Vector3f V) const;
	FORCEINLINE Vector3f UnrotateVector(Vector3f V) const;

}; 


FORCEINLINE Quat Quat::operator*=(const Quat& Q)
{
	VectorRegister A = VectorLoadAligned(this);
	VectorRegister B = VectorLoadAligned(&Q);
	VectorRegister Result;
	VectorQuaternionMultiply(&Result, &A, &B);
	VectorStoreAligned(Result, this);

	CheckNaN();

	return *this;
}

FORCEINLINE Vector3f Quat::operator*(const Vector3f& V) const
{
	return RotateVector(V);
}

FORCEINLINE Quat Quat::operator*(const Quat& Q) const
{
	Quat Result;
	VectorQuaternionMultiply(&Result, this, &Q);

	Result.CheckNaN();

	return Result;
}

FORCEINLINE Quat Quat::Inverse() const
{
	return Quat(-X, -Y, -Z, W);
}

FORCEINLINE void Quat::Normalize(float Tolerance /*= SMALL_NUMBER*/)
{
#if PLATFORM_ENABLE_VECTORINTRINSICS
	const VectorRegister Vector = VectorLoadAligned(this);

	const VectorRegister SquareSum = VectorDot4(Vector, Vector);
	const VectorRegister NonZeroMask = VectorCompareGE(SquareSum, VectorLoadFloat1(&Tolerance));
	const VectorRegister InvLength = VectorReciprocalSqrtAccurate(SquareSum);
	const VectorRegister NormalizedVector = VectorMultiply(InvLength, Vector);
	VectorRegister Result = VectorSelect(NonZeroMask, NormalizedVector, GlobalVectorConstants::Float0001);

	VectorStoreAligned(Result, this);
#else
	const float SquareSum = X * X + Y * Y + Z * Z + W * W;

	if (SquareSum >= Tolerance)
	{
		const float Scale = Math::InvSqrt(SquareSum);

		X *= Scale;
		Y *= Scale;
		Z *= Scale;
		W *= Scale;
	}
	else
	{
		*this = Quat::Identity;
	}
#endif // PLATFORM_ENABLE_VECTORINTRINSICS
}

FORCEINLINE Quat Quat::GetNormalized(float Tolerance /*= SMALL_NUMBER*/) const
{
	Quat Result(*this);
	Result.Normalize(Tolerance);
	return Result;
}

FORCEINLINE bool Quat::Equals(const Quat& Q, float Tolerance /*= KINDA_SMALL_NUMBER*/) const
{
#if PLATFORM_ENABLE_VECTORINTRINSICS
	const VectorRegister ToleranceV = VectorLoadFloat1(&Tolerance);
	const VectorRegister A = VectorLoadAligned(this);
	const VectorRegister B = VectorLoadAligned(&Q);

	const VectorRegister RotationSub = VectorAbs(VectorSubtract(A, B));
	const VectorRegister RotationAdd = VectorAbs(VectorAdd(A, B));
	return !VectorAnyGreaterThan(RotationSub, ToleranceV) || !VectorAnyGreaterThan(RotationAdd, ToleranceV);
#else
	return (FMath::Abs(X - Q.X) <= Tolerance && FMath::Abs(Y - Q.Y) <= Tolerance && FMath::Abs(Z - Q.Z) <= Tolerance && FMath::Abs(W - Q.W) <= Tolerance)
		|| (FMath::Abs(X + Q.X) <= Tolerance && FMath::Abs(Y + Q.Y) <= Tolerance && FMath::Abs(Z + Q.Z) <= Tolerance && FMath::Abs(W + Q.W) <= Tolerance);
#endif // PLATFORM_ENABLE_VECTORINTRINSICS
}

FORCEINLINE bool Quat::IsNormalized() const
{
	return (Math::Abs(1.f - SizeSquared()) < THRESH_QUAT_NORMALIZED);
}

FORCEINLINE float Quat::Size() const
{
	return Math::Sqrt(X * X + Y * Y + Z * Z + W * W);
}

FORCEINLINE float Quat::SizeSquared() const
{
	return (X * X + Y * Y + Z * Z + W * W);
}

FORCEINLINE float Quat::GetAngle() const
{
	return 2.f * Math::Acos(W);
}

FORCEINLINE Vector3f Quat::GetAxisX() const
{
	return RotateVector(Vector3f(1.f, 0.f, 0.f));
}

FORCEINLINE Vector3f Quat::GetAxisY() const
{
	return RotateVector(Vector3f(0.f, 1.f, 0.f));
}

FORCEINLINE Vector3f Quat::GetAxisZ() const
{
	return RotateVector(Vector3f(0.f, 0.f, 1.f));
}

FORCEINLINE Vector3f Quat::GetForwardVector() const
{
	return GetAxisX();
}

FORCEINLINE Vector3f Quat::GetRightVector() const
{
	return GetAxisY();
}

FORCEINLINE Vector3f Quat::GetUpVector() const
{
	return GetAxisZ();
}

FORCEINLINE Vector3f Quat::Vector() const
{
	return GetAxisX();
}

FORCEINLINE Vector3f Quat::RotateVector(Vector3f V) const
{
	// http://people.csail.mit.edu/bkph/articles/Quaternions.pdf
	// V' = V + 2w(Q x V) + (2Q x (Q x V))
	// refactor:
	// V' = V + w(2(Q x V)) + (Q x (2(Q x V)))
	// T = 2(Q x V);
	// V' = V + w*(T) + (Q x T)

	const Vector3f Q(X, Y, Z);
	const Vector3f T = 2.f * Vector3f::CrossProduct(Q, V);
	const Vector3f Result = V + (W * T) + Vector3f::CrossProduct(Q, T);
	return Result;
}

FORCEINLINE Vector3f Quat::UnrotateVector(Vector3f V) const
{
	const Vector3f Q(-X, -Y, -Z); // Inverse
	const Vector3f T = 2.f * Vector3f::CrossProduct(Q, V);
	const Vector3f Result = V + (W * T) + Vector3f::CrossProduct(Q, T);
	return Result;
}

};
