// this tiny game engine is rendered by CPU. I make it for learing the pipline and rasterization. This is the first step to realise my own game engine. so just do it!

#pragma once
#include "CoreTypes.h"
#include <math.h>
#include <stdlib.h>

namespace LXF {

#ifndef ENABLE_NAN_DIAGNOSTIC
#if BUILD_DEBUG
#define ENABLE_NAN_DIAGNOSTIC 1
#else
#define ENABLE_NAN_DIAGNOSTIC 0
#endif
#endif

#undef PI 
#define PI (3.1415926535897932f)
#define INV_PI			(0.31830988618f)
#define HALF_PI			(1.57079632679f)
#define SMALL_NUMBER (1.e-8f)
#define KINDA_SMALL_NUMBER (1.e-4f)
#define BIG_NUMBER (3.4e+38f)
#define FLOAT_NON_FRACTIONAL (8388608.f)

#define MAX_FLT 3.402823466e+38F

#define DELTA			(0.00001f)

#define THRESH_POINT_ON_PLANE			(0.10f)		/* Thickness of plane for front/back/inside test */
#define THRESH_POINT_ON_SIDE			(0.20f)		/* Thickness of polygon side's side-plane for point-inside/outside/on side test */
#define THRESH_POINTS_ARE_SAME			(0.00002f)	/* Two points are same if within this distance */
#define THRESH_POINTS_ARE_NEAR			(0.015f)	/* Two points are near if within this distance and can be combined if imprecise math is ok */
#define THRESH_NORMALS_ARE_SAME			(0.00002f)	/* Two normal points are same if within this distance */
#define THRESH_UVS_ARE_SAME			    (0.0009765625f)/* Two UV are same if within this threshold (1.0f/1024f) */
	/* Making this too large results in incorrect CSG classification and disaster */
#define THRESH_VECTORS_ARE_NEAR			(0.0004f)	/* Two vectors are near if within this distance and can be combined if imprecise math is ok */
													/* Making this too large results in lighting problems due to inaccurate texture coordinates */
#define THRESH_SPLIT_POLY_WITH_PLANE	(0.25f)		/* A plane splits a polygon in half */
#define THRESH_SPLIT_POLY_PRECISELY		(0.01f)		/* A plane exactly splits a polygon */
#define THRESH_ZERO_NORM_SQUARED		(0.0001f)	/* Size of a unit normal that is considered "zero", squared */
#define THRESH_NORMALS_ARE_PARALLEL		(0.999845f)	/* Two unit vectors are parallel if abs(A dot B) is greater than or equal to this. This is roughly cosine(1.0 degrees). */
#define THRESH_NORMALS_ARE_ORTHOGONAL	(0.017455f)	/* Two unit vectors are orthogonal (perpendicular) if abs(A dot B) is less than or equal this. This is roughly cosine(89.0 degrees). */

#define THRESH_VECTOR_NORMALIZED		(0.01f)		/** Allowed error for a normalized vector (against squared magnitude) */
#define THRESH_QUAT_NORMALIZED			(0.01f)		/** Allowed error for a normalized quaternion (against squared magnitude) */

// Forward declarations.
template <typename T>
class Vector;
template <typename T>
class Vector2;
template <typename T>
class Vector4;
class Plane;
class Rotator;
class Matrix;
class Quat;


typedef Vector<float> Vector3f;
typedef Vector<int> Vector3i;
typedef Vector2<float> Vector2f;
typedef Vector2<int> Vector2i;
typedef Vector4<float> Vector4f;
typedef Vector4<int> Vector4i;

namespace EAxis
{
	enum Type
	{
		None,
		X,
		Y,
		Z,
	};
}

class Math {
public:

	static float Sqrt(float Value) { return sqrtf(Value); }
	static float InvSqrt(float Value) { return 1.0f/sqrtf(Value); }
	static float Pow(float A, float B) { return powf(A,B); }
	static float Sin(float Value) { return sinf(Value); }
	static float Asin(float Value) { return asinf((Value < -1.f) ? -1.f : ((Value<1.f) ? Value : 1.f)); }
	static float Cos(float Value) { return cosf(Value); }
	static float Acos(float Value) { return cosf((Value < -1.f) ? -1.f : ((Value<1.f) ? Value : 1.f)); }
	static float Tan(float Value) { return tanf(Value); }
	static float Atan(float Value) { return atanf(Value); }
	static float Atan2(float Y, float X);

#define FASTASIN_HALF_PI (1.5707963050f)
	static FORCEINLINE float FastAsin(float Value)
	{
		// Clamp input to [-1,1].
		bool nonnegative = (Value >= 0.0f);
		float x = Math::Abs(Value);
		float omx = 1.0f - x;
		if (omx < 0.0f)
		{
			omx = 0.0f;
		}
		float root = Math::Sqrt(omx);
		// 7-degree minimax approximation
		float result = ((((((-0.0012624911f * x + 0.0066700901f) * x - 0.0170881256f) * x + 0.0308918810f) * x - 0.0501743046f) * x + 0.0889789874f) * x - 0.2145988016f) * x + FASTASIN_HALF_PI;
		result *= root;  // acos(|x|)
		// acos(x) = pi - acos(-x) when x < 0, asin(x) = pi/2 - acos(x)
		return (nonnegative ? FASTASIN_HALF_PI - result : result - FASTASIN_HALF_PI);
	}
#undef FASTASIN_HALF_PI


	static bool IsNaN(float A) { return ((*(uint32*)&A) & 0x7FFFFFFFU) > 0x7F800000U; }
	static bool IsNaN(double A) { return ((*(uint64*)&A) & 0x7FFFFFFFFFFFFFFFULL) > 0x7FF0000000000000ULL; }
	static bool IsFinite(float A) { return ((*(uint32*)&A) & 0x7FFFFFFFU) != 0x7F800000U; }
	static bool IsFinite(double A) { return ((*(uint64*)&A) & 0x7FFFFFFFFFFFFFFFULL) != 0x7FF0000000000000ULL; }

	// Returns e^Value
	static FORCEINLINE float Exp(float Value) { return expf(Value); }
	// Returns 2^Value
	static FORCEINLINE float Exp2(float Value) { return powf(2.f, Value); /*exp2f(Value);*/ }
	static FORCEINLINE float Loge(float Value) { return logf(Value); }
	static FORCEINLINE float LogX(float Base, float Value) { return Loge(Value) / Loge(Base); }
	// 1.0 / Loge(2) = 1.4426950f
	static FORCEINLINE float Log2(float Value) { return Loge(Value) * 1.4426950f; }

	static FORCEINLINE bool IsNearlyEqual(float A, float B, float ErrorTolerance = SMALL_NUMBER)
	{
		return Abs<float>(A - B) <= ErrorTolerance;
	}

	static FORCEINLINE bool IsNearlyEqual(double A, double B, double ErrorTolerance = SMALL_NUMBER)
	{
		return Abs<double>(A - B) <= ErrorTolerance;
	}

	static FORCEINLINE bool IsNearlyZero(float Value, float ErrorTolerance = SMALL_NUMBER)
	{
		return Abs<float>(Value) <= ErrorTolerance;
	}

	static int32 Rand() { return rand(); }
	static void RandInit(int32 Seed) { srand(Seed); }
	static float FRand() { return Rand() / (float)RAND_MAX; }
	static float SRand();

	template <class T>
	static T Clamp(const T X, const T Min, const T Max) { return X < Min ? Min : X < Max ? X : Max; }

	static FORCEINLINE float FloatSelect(float Comparand, float ValueGEZero, float ValueLTZero)
	{
		return Comparand >= 0.f ? ValueGEZero : ValueLTZero;
	}
	static FORCEINLINE double FloatSelect(double Comparand, double ValueGEZero, double ValueLTZero)
	{
		return Comparand >= 0.f ? ValueGEZero : ValueLTZero;
	}

	template<class T>
	static T Abs(const T A) { return (A >= (T)0) ? A : -A; }

	template<>
	static float Abs(const float A) { return fabsf(A); }

	static float TruncToFloat(float F) { return truncf(F); }

	static float Fmod(float X, float Y)
	{
		const float AbsY = fabsf(Y);
		if (AbsY <= 1.e-8f)
		{
			return 0.f;
		}
		const float Div = (X / Y);
		const float Quotient = fabsf(Div) < FLOAT_NON_FRACTIONAL ? TruncToFloat(Div) : Div;
		float IntPortion = Y * Quotient;

		if (fabsf(IntPortion) > fabsf(X))
		{
			IntPortion = X;
		}

		const float Result = X - IntPortion;
		return Math::Clamp(Result, -AbsY, AbsY);
	}

	/** Performs a linear interpolation between two values, Alpha ranges from 0-1 */
	template< class T, class U >
	static T Lerp(const T& A, const T& B, const U& Alpha)
	{
		return (T)(A + Alpha * (B - A));
	}

	/** Performs a linear interpolation between two values, Alpha ranges from 0-1. Handles full numeric range of T */
	template< class T >
	static T LerpStable(const T& A, const T& B, double Alpha)
	{
		return (T)((A * (1.0 - Alpha)) + (B * Alpha));
	}

	/** Performs a linear interpolation between two values, Alpha ranges from 0-1. Handles full numeric range of T */
	template< class T >
	static T LerpStable(const T& A, const T& B, float Alpha)
	{
		return (T)((A * (1.0f - Alpha)) + (B * Alpha));
	}

	/** Performs a 2D linear interpolation between four values values, FracX, FracY ranges from 0-1 */
	template< class T, class U >
	static T BiLerp(const T& P00, const T& P10, const T& P01, const T& P11, const U& FracX, const U& FracY)
	{
		return Lerp(
			Lerp(P00, P10, FracX),
			Lerp(P01, P11, FracX),
			FracY
		);
	}

	template<class T>
	static auto RadiansToDegrees(T const& RadVal) -> decltype(RadVal * (180.f / PI))
	{
		return RadVal * (180.f / PI);
	}

	template<class T>
	static auto DegreesToRadians(T const& DegVal) -> decltype(DegVal * (PI / 180.f))
	{
		return DegVal * (PI / 180.f);
	}

	static void SinCos(float* ScalarSin, float* ScalarCos, float Value)
	{
		// Map Value to y in [-pi,pi], x = 2*pi*quotient + remainder.
		float quotient = (INV_PI*0.5f)*Value;
		if (Value >= 0.0f)
		{
			quotient = (float)((int)(quotient + 0.5f));
		}
		else
		{
			quotient = (float)((int)(quotient - 0.5f));
		}
		float y = Value - (2.0f*PI)*quotient;

		// Map y to [-pi/2,pi/2] with sin(y) = sin(Value).
		float sign;
		if (y > HALF_PI)
		{
			y = PI - y;
			sign = -1.0f;
		}
		else if (y < -HALF_PI)
		{
			y = -PI - y;
			sign = -1.0f;
		}
		else
		{
			sign = +1.0f;
		}

		float y2 = y * y;

		// 11-degree minimax approximation
		*ScalarSin = (((((-2.3889859e-08f * y2 + 2.7525562e-06f) * y2 - 0.00019840874f) * y2 + 0.0083333310f) * y2 - 0.16666667f) * y2 + 1.0f) * y;

		// 10-degree minimax approximation
		float p = ((((-2.6051615e-07f * y2 + 2.4760495e-05f) * y2 - 0.0013888378f) * y2 + 0.041666638f) * y2 - 0.5f) * y2 + 1.0f;
		*ScalarCos = sign * p;
	}

	/** Returns 1, 0, or -1 depending on relation of T to 0 */
	template< class T >
	static FORCEINLINE T Sign(const T A)
	{
		return (A > (T)0) ? (T)1 : ((A < (T)0) ? (T)-1 : (T)0);
	}

	/** Returns higher value in a generic way */
	template< class T >
	static FORCEINLINE T Max(const T A, const T B)
	{
		return (A >= B) ? A : B;
	}

	/** Returns lower value in a generic way */
	template< class T >
	static FORCEINLINE T Min(const T A, const T B)
	{
		return (A <= B) ? A : B;
	}

	/** Returns highest of 3 values */
	template< class T >
	static FORCEINLINE T Max3(const T A, const T B, const T C)
	{
		return Max(Max(A, B), C);
	}

	/** Returns lowest of 3 values */
	template< class T >
	static FORCEINLINE T Min3(const T A, const T B, const T C)
	{
		return Min(Min(A, B), C);
	}

	/** Multiples value by itself */
	template< class T >
	static FORCEINLINE T Square(const T A)
	{
		return A * A;
	}

	static Vector3f RayPlaneIntersection(const Vector3f& RayOrigin, const Vector3f& RayDirection, const Plane& Plane);
	static Vector3f LinePlaneIntersection(const Vector3f& Point1, const Vector3f& Point2, const Vector3f& PlaneOrigin, const Vector3f& PlaneNormal);
	static Vector3f LinePlaneIntersection(const Vector3f& Point1, const Vector3f& Point2, const Plane& Plane);

};

};
