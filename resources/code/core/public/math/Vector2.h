// this tiny game engine is rendered by CPU. I make it for learing the pipline and rasterization. This is the first step to realise my own game engine. so just do it!

#pragma once

#include "platform/Platform.h"
#include "MathUtility.h"

namespace LXF {

template <typename T>
class Vector2
{
public:
	T X, Y;

	static const Vector2 ZeroVector;
	static const Vector2 OneVector;
	static const Vector2 ForwardVector;
	static const Vector2 BackwardVector;
	static const Vector2 RightVector;
	static const Vector2 LeftVector;

#if ENABLE_NAN_DIAGNOSTIC
	void CheckNaN() const
	{
		if (ContainsNaN())
		{
			*const_cast<Vector2*>(this) = ZeroVector;
		}
	}
#else
	void CheckNaN() const {}
#endif

	Vector2() : X(0), Y(0) { };

	Vector2(T InF) : X(InF), Y(InF) { CheckNaN(); }
	Vector2(T InX, T InY) : X(InX), Y(InY) { CheckNaN(); }
	explicit Vector2(const Vector<T>& V);
	explicit Vector2(const Vector4<T>& V);

	T operator^(const Vector2& V) const { return X * V.Y - Y * V.X; }
	static T CrossProduct(const Vector2& A, const Vector2& B) { return A ^ B; }

	T operator|(const Vector2& V) const { return X * V.X + Y * V.Y; }
	static T DotProduct(const Vector2& A, const Vector2& B) { return A | B; }

	FORCEINLINE Vector2 operator+(const Vector2& V) const { return Vector2(X + V.X, Y + V.Y); }
	FORCEINLINE Vector2 operator-(const Vector2& V) const { return Vector2(X - V.X, Y - V.Y); }
	FORCEINLINE Vector2 operator*(T Scale) const { return Vector2(X * Scale, Y * Scale); }
	FORCEINLINE Vector2 operator/(T Scale) const { const T RScale = 1.0f / Scale; return Vector2(X * RScale, Y * RScale); }
	FORCEINLINE bool operator==(const Vector2& V) const { return X == V.X && Y == V.Y; }
	FORCEINLINE bool operator!=(const Vector2& V) const { return X != V.X || Y != V.Y; }
	FORCEINLINE Vector2 operator-() const { return Vector2(-X, -Y); }

	FORCEINLINE Vector2 operator+=(const Vector2& V) { X += V.X; Y += V.Y; CheckNaN(); return *this; }
	FORCEINLINE Vector2 operator-=(const Vector2& V) { X -= V.X; Y -= V.Y; CheckNaN(); return *this; }
	FORCEINLINE Vector2 operator*=(T Scale) { X *= Scale; Y *= Scale; CheckNaN(); return *this; }
	FORCEINLINE Vector2 operator/=(T Scale) { const T RScale = 1.0f / Scale; X *= RScale; Y *= RScale; CheckNaN(); return *this; }
	FORCEINLINE Vector2 operator*=(const Vector2& V) { X *= V.X; Y *= V.Y; CheckNaN(); return *this; }
	FORCEINLINE Vector2 operator/=(const Vector2& V) { X /= V.X; Y /= V.Y; CheckNaN(); return *this; }

	bool Equals(const Vector2& V, float Tolerance) const { return Math::Abs(X - V.X) <= Tolerance && Math::Abs(Y - V.Y) <= Tolerance; }

	T& operator[](int32 Index) { check(Index >= 0 && Index < 3); return (&X)[Index]; }

	float Size() const { return Math::Sqrt(X * X + Y * Y ); };
	float SizeSquared() const { return X * X + Y * Y; };

	bool IsNearlyZero(float Tolerance=KINDA_SMALL_NUMBER) const { return Math::Abs(X) <= Tolerance && Math::Abs(Y) <= Tolerance; }

	bool Normalize(float Tolerance=SMALL_NUMBER);
	Vector2 GetSafeNormal(float Tolerance=SMALL_NUMBER) const;
	bool ContainsNaN() const { return (!Math::IsFinite(X) || !Math::IsFinite(Y)); }

}; 

template <typename T>
FORCEINLINE Vector2<T> operator*(T Scale, const Vector2<T>& V) { return V.operator*(Scale); }

template <typename T>
Vector2<T> Vector2<T>::GetSafeNormal(float Tolerance/*=SMALL_NUMBER*/) const
{
	const float SquareSum = X * X + Y * Y;
	if (SquareSum == 1.0f)
	{
		return *this;
	}
	else if (SquareSum < Tolerance)
	{
		return Vector2::ZeroVector;
	}
	const float Scale = Math::InvSqrt(SquareSum);
	return Vector2(X * Scale, Y * Scale);
}

template <typename T>
bool Vector2<T>::Normalize(float Tolerance/*=SMALL_NUMBER*/)
{
	const float SquareSum = X * X + Y * Y;
	if (SquareSum > Tolerance)
	{
		const float Scale = Math::InvSqrt(SquareSum);
		X *= Scale; Y *= Scale;
		return true;
	}
	return false;
}

};
