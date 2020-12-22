// this tiny game engine is rendered by CPU. I make it for learing the pipline and rasterization. This is the first step to realise my own game engine. so just do it!

#pragma once

#include "platform/Platform.h"
#include "MathUtility.h"

namespace LXF {

template <typename T>
class Vector
{
public:
	T X, Y, Z;

	static const Vector ZeroVector;
	static const Vector OneVector;
	static const Vector UpVector;
	static const Vector DownVector;
	static const Vector ForwardVector;
	static const Vector BackwardVector;
	static const Vector RightVector;
	static const Vector LeftVector;

#if ENABLE_NAN_DIAGNOSTIC
	void CheckNaN() const
	{
		if (ContainsNaN())
		{
			*const_cast<Vector*>(this) = ZeroVector;
		}
	}
#else
	void CheckNaN() const {}
#endif

	Vector() : X(0), Y(0), Z(0) { };

	Vector(T InF) : X(InF), Y(InF), Z(InF) { CheckNaN(); }
	Vector(T InX, T InY, T InZ) : X(InX), Y(InY), Z(InZ) { CheckNaN(); }
	explicit Vector(const Vector2<T> V, float InZ);
	Vector(const Vector4<T>& V);


	Vector operator^(const Vector& V) const { return Vector(Y * V.Z - Z * V.Y, Z * V.X - X * V.Z, X * V.Y - Y * V.X); }
	static Vector CrossProduct(const Vector& A, const Vector& B) { return A ^ B; }

	T operator|(const Vector& V) const { return X * V.X + Y * V.Y + Z * V.Z; }
	static T DotProduct(const Vector& A, const Vector& B) { return A | B; }

	FORCEINLINE Vector operator+(const Vector& V) const { return Vector(X + V.X, Y + V.Y, Z + V.Z); }
	FORCEINLINE Vector operator-(const Vector& V) const { return Vector(X - V.X, Y - V.Y, Z - V.Z); }
	FORCEINLINE Vector operator*(T Scale) const { return Vector(X * Scale, Y * Scale, Z * Scale); }
	FORCEINLINE Vector operator/(T Scale) const { const T RScale = 1.0f / Scale; return Vector(X * RScale, Y * RScale, Z * RScale); }
	FORCEINLINE bool operator==(const Vector& V) const { return X == V.X && Y == V.Y && Z == V.Z; }
	FORCEINLINE bool operator!=(const Vector& V) const { return X != V.X || Y != V.Y || Z != V.Z; }
	FORCEINLINE Vector operator-() const { return Vector(-X, -Y, -Z); }

	FORCEINLINE Vector operator+=(const Vector& V) { X += V.X; Y += V.Y; Z += V.Z; CheckNaN(); return *this; }
	FORCEINLINE Vector operator-=(const Vector& V) { X -= V.X; Y -= V.Y; Z -= V.Z; CheckNaN(); return *this; }
	FORCEINLINE Vector operator*=(T Scale) { X *= Scale; Y *= Scale; Z *= Scale; CheckNaN(); return *this; }
	FORCEINLINE Vector operator/=(T Scale) { const T RScale = 1.0f / Scale; X *= RScale; Y *= RScale; Z *= RScale; CheckNaN(); return *this; }
	FORCEINLINE Vector operator*=(const Vector& V) { X *= V.X; Y *= V.Y; Z *= V.Z; CheckNaN(); return *this; }
	FORCEINLINE Vector operator/=(const Vector& V) { X /= V.X; Y /= V.Y; Z /= V.Z; CheckNaN(); return *this; }

	bool Equals(const Vector& V, float Tolerance) const { return Math::Abs(X - V.X) <= Tolerance && Math::Abs(Y - V.Y) <= Tolerance && Math::Abs(Z - V.Z) <= Tolerance; }

	T& operator[](int32 Index) { check(Index >= 0 && Index < 3); return (&X)[Index]; }

	float Size() const { return Math::Sqrt(X * X + Y * Y + Z * Z); };
	float SizeSquared() const { return X * X + Y * Y + Z * Z; };
	float Size2D() const { return Math::Sqrt(X * X + Y * Y); };
	float SizeSquared2D() const { return X * X + Y * Y; };

	bool IsNearlyZero(float Tolerance=KINDA_SMALL_NUMBER) const { return Math::Abs(X) <= Tolerance && Math::Abs(Y) <= Tolerance && Math::Abs(Z) <= Tolerance; }

	bool Normalize(float Tolerance=SMALL_NUMBER);
	Vector GetSafeNormal(float Tolerance=SMALL_NUMBER) const;
	bool ContainsNaN() const { return (!Math::IsFinite(X) || !Math::IsFinite(Y) || !Math::IsFinite(Z)); }

}; 

template <typename T>
FORCEINLINE Vector<T> operator*(T Scale, const Vector<T>& V) { return V.operator*(Scale); }

template <typename T>
Vector<T> Vector<T>::GetSafeNormal(float Tolerance/*=SMALL_NUMBER*/) const
{
	const float SquareSum = X * X + Y * Y + Z * Z;
	if (SquareSum == 1.0f)
	{
		return *this;
	}
	else if (SquareSum < Tolerance)
	{
		return Vector::ZeroVector;
	}
	const float Scale = Math::InvSqrt(SquareSum);
	return Vector(X * Scale, Y * Scale, Z * Scale);
}

template <typename T>
bool Vector<T>::Normalize(float Tolerance/*=SMALL_NUMBER*/)
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
