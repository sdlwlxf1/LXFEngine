// this tiny game engine is rendered by CPU. I make it for learing the pipline and rasterization. This is the first step to realise my own game engine. so just do it!

#pragma once

#include "platform/Platform.h"
#include "MathUtility.h"

namespace LXF {

template <typename T>
class Vector4
{
public:
	T X, Y, Z, W;

	static const Vector4 ZeroVector;
	static const Vector4 OneVector;
	static const Vector4 UpVector;
	static const Vector4 DownVector;
	static const Vector4 ForwardVector;
	static const Vector4 BackwardVector;
	static const Vector4 RightVector;
	static const Vector4 LeftVector;

#if ENABLE_NAN_DIAGNOSTIC
	void CheckNaN() const
	{
		if (ContainsNaN())
		{
			*const_cast<Vector4*>(this) = ZeroVector;
		}
	}
#else
	void CheckNaN() const {}
#endif

	Vector4(T InF) : X(InF), Y(InF), Z(InF), W(InF) { CheckNaN(); }
	Vector4(T InX = 0, T InY = 0, T InZ = 0, T InW = 1) : X(InX), Y(InY), Z(InZ), W(InW) { CheckNaN(); }
	Vector4(const Vector<T>& InVector, T InW = 1);

	Vector4 operator^(const Vector4& V) const { return Vector4( Y * V.Z - Z * V.Y, Z * V.X - X * V.Z, X * V.Y - Y * V.X, 0.0f ); }
	static Vector4 CrossProduct(const Vector4& A, const Vector4& B) { return A ^ B; }

	FORCEINLINE Vector4 operator+(const Vector4& V) const { return Vector4(X + V.X, Y + V.Y, Z + V.Z, W + V.W); }
	FORCEINLINE Vector4 operator-(const Vector4& V) const { return Vector4(X - V.X, Y - V.Y, Z - V.Z, W - V.W); }
	FORCEINLINE Vector4 operator*(T Scale) const { return Vector4(X * Scale, Y * Scale, Z * Scale, W * Scale); }
	FORCEINLINE Vector4 operator/(T Scale) const { const T RScale = 1.0f / Scale; return Vector4(X * RScale, Y * RScale, Z * RScale, W * RScale); }
	FORCEINLINE bool operator==(const Vector4& V) const { return X == V.X && Y == V.Y && Z == V.Z && W == V.W; }
	FORCEINLINE bool operator!=(const Vector4& V) const { return X != V.X || Y != V.Y || Z != V.Z || W != V.W; }
	FORCEINLINE Vector4 operator-() const { return Vector4(-X, -Y, -Z, -W); }

	FORCEINLINE Vector4 operator+=(const Vector4& V) { X += V.X; Y += V.Y; Z += V.Z; W += V.W; CheckNaN(); return *this; }
	FORCEINLINE Vector4 operator-=(const Vector4& V) { X -= V.X; Y -= V.Y; Z -= V.Z; W -= V.W; CheckNaN(); return *this; }
	FORCEINLINE Vector4 operator*=(T Scale) { X *= Scale; Y *= Scale; Z *= Scale; W *= Scale; CheckNaN(); return *this; }
	FORCEINLINE Vector4 operator/=(T Scale) { const T RScale = 1.0f / Scale; X *= RScale; Y *= RScale; Z *= RScale; W *= RScale; CheckNaN(); return *this; }
	FORCEINLINE Vector4 operator*=(const Vector4& V) { X *= V.X; Y *= V.Y; Z *= V.Z; W *= V.W; CheckNaN(); return *this; }
	FORCEINLINE Vector4 operator/=(const Vector4& V) { X /= V.X; Y /= V.Y; Z /= V.Z; W /= V.W; CheckNaN(); return *this; }

	bool Equals(const Vector4& V, float Tolerance) const { return Math::Abs(X - V.X) <= Tolerance && Math::Abs(Y - V.Y) <= Tolerance && Math::Abs(Z - V.Z) <= Tolerance; }

	T& operator[](int32 Index) { check(Index >= 0 && Index < 4); return (&X)[Index]; }

	float Size() const { return Math::Sqrt(X * X + Y * Y + Z * Z); };
	float SizeSquared() const { return X * X + Y * Y + Z * X; };
	float Size2D() const { return Math::Sqrt(X * X + Y * Y); };
	float SizeSquared2D() const { return X * X + Y * Y; };

	bool IsNearlyZero(float Tolerance=KINDA_SMALL_NUMBER) const { return Math::Abs(X) <= Tolerance && Math::Abs(Y) <= Tolerance && Math::Abs(Z) <= Tolerance; }

	bool Normalize(float Tolerance=SMALL_NUMBER);
	Vector4 GetSafeNormal(float Tolerance=SMALL_NUMBER) const;
	bool ContainsNaN() const { return (!Math::IsFinite(X) || !Math::IsFinite(Y) || !Math::IsFinite(Z)); }

}; 

template <typename T>
FORCEINLINE Vector4<T> operator*(T Scale, const Vector4<T>& V) { return V.operator*(Scale); }

template <typename T>
Vector4<T> Vector4<T>::GetSafeNormal(float Tolerance/*=SMALL_NUMBER*/) const
{
	const float SquareSum = X * X + Y * Y + Z * Z;
	if (SquareSum == 1.0f)
	{
		return *this;
	}
	else if (SquareSum < Tolerance)
	{
		return Vector4::ZeroVector;
	}
	const float Scale = Math::InvSqrt(SquareSum);
	return Vector4(X * Scale, Y * Scale, Z * Scale);
}

template <typename T>
bool Vector4<T>::Normalize(float Tolerance/*=SMALL_NUMBER*/)
{
	const float SquareSum = X * X + Y * Y + Z * Z;
	if (SquareSum > Tolerance)
	{
		const float Scale = Math::InvSqrt(SquareSum);
		X *= Scale; Y *= Scale; Z *= Scale;
		return true;
	}
	return false;
}

};
