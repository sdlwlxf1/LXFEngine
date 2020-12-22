// this tiny game engine is rendered by CPU. I make it for learing the pipline and rasterization. floathis is the first step to realise my own game engine. so just do it!

#pragma once

#include "platform/Platform.h"
#include "Vector.h"
#include "MathUtility.h"

namespace LXF {

class Rotator
{
public:
	float Pitch, Yaw, Roll;

	static const Rotator ZeroRotator;

#if ENABLE_NAN_DIAGNOSTIC
	void CheckNaN() const {
		if (ContainsNaN())
		{
			*const_cast<Rotator*>(this) = ZeroRotator;
		}
	}
#else
	void CheckNaN() const {}
#endif	

	FORCEINLINE Rotator() : Pitch(0), Yaw(0), Roll(0) { };
	FORCEINLINE Rotator(float InF) : Pitch(InF), Yaw(InF), Roll(InF) { CheckNaN(); }
	FORCEINLINE Rotator(float InPitch, float InYaw, float InRoll) : Pitch(InPitch), Yaw(InYaw), Roll(InRoll) { CheckNaN(); }

	FORCEINLINE Rotator operator+(const Rotator& R) const { return Rotator(Pitch + R.Pitch, Yaw + R.Yaw, Roll + R.Roll); }
	FORCEINLINE Rotator operator-(const Rotator& R) const { return Rotator(Pitch - R.Pitch, Yaw - R.Yaw, Roll - R.Roll); }
	FORCEINLINE Rotator operator*(float Scale) const { return Rotator(Pitch * Scale, Yaw * Scale, Roll * Scale); }
	FORCEINLINE Rotator operator/(float Scale) const { const float RScale = 1.0f / Scale; return Rotator(Pitch * RScale, Yaw * RScale, Roll * RScale); }
	FORCEINLINE bool operator==(const Rotator& R) const { return Pitch == R.Pitch && Yaw == R.Yaw && Roll == R.Roll; }
	FORCEINLINE bool operator!=(const Rotator& R) const { return Pitch != R.Pitch || Yaw != R.Yaw || Roll != R.Roll; }
	FORCEINLINE Rotator operator-() const { return Rotator(-Pitch, -Yaw, -Roll); }

	FORCEINLINE Rotator operator+=(const Rotator& R) { Pitch += R.Pitch; Yaw += R.Yaw; Roll += R.Roll; CheckNaN(); return *this; }
	FORCEINLINE Rotator operator-=(const Rotator& R) { Pitch -= R.Pitch; Yaw -= R.Yaw; Roll -= R.Roll; CheckNaN(); return *this; }
	FORCEINLINE Rotator operator*=(float Scale) { Pitch *= Scale; Yaw *= Scale; Roll *= Scale; CheckNaN(); return *this; }
	FORCEINLINE Rotator operator/=(float Scale) { const float RScale = 1.0f / Scale; Pitch *= RScale; Yaw *= RScale; Roll *= RScale; CheckNaN(); return *this; }
	FORCEINLINE Rotator operator*=(const Rotator& R) { Pitch *= R.Pitch; Yaw *= R.Yaw; Roll *= R.Roll; CheckNaN(); return *this; }
	FORCEINLINE Rotator operator/=(const Rotator& R) { Pitch /= R.Pitch; Yaw /= R.Yaw; Roll /= R.Roll; CheckNaN(); return *this; }

	FORCEINLINE float& operator[](int32 Index) { check(Index >= 0 && Index < 3); return (&Pitch)[Index]; }

	FORCEINLINE static float ClampAxis(float Angle);
	FORCEINLINE static float NormalizeAxis(float Angle);

	FORCEINLINE Rotator Clamp() const { Rotator(ClampAxis(Pitch), ClampAxis(Yaw), ClampAxis(Roll)); }

	FORCEINLINE Rotator GetNormalized() const { Rotator Rot = *this; Rot.Normalize(); return Rot; }

	FORCEINLINE bool IsNearlyZero(float Tolerance) const {
		return (Math::Abs(NormalizeAxis(Pitch)) <= Tolerance)
			&& (Math::Abs(NormalizeAxis(Yaw)) <= Tolerance)
			&& (Math::Abs(NormalizeAxis(Roll)) <= Tolerance);
	}

	FORCEINLINE bool IsZero() const { return (ClampAxis(Pitch) == 0.f) && (ClampAxis(Yaw) == 0.f) && (ClampAxis(Roll) == 0.f); }

	FORCEINLINE bool Equals(const Rotator& R, float Tolerance) const
	{
		return (Math::Abs(NormalizeAxis(Pitch - R.Pitch)) <= Tolerance)
			&& (Math::Abs(NormalizeAxis(Yaw - R.Yaw)) <= Tolerance)
			&& (Math::Abs(NormalizeAxis(Roll - R.Roll)) <= Tolerance);
	}

	FORCEINLINE bool ContainsNaN() const { return (!Math::IsFinite(Pitch) || !Math::IsFinite(Yaw) || !Math::IsFinite(Roll)); }

	FORCEINLINE void Normalize() { Pitch = NormalizeAxis(Pitch); Yaw = NormalizeAxis(Yaw); Roll = NormalizeAxis(Roll); }

	Vector3f Vector() const;
	Quat Quaternion() const;

	FORCEINLINE Vector3f RotateVector(const Vector3f& V) const;
	FORCEINLINE Vector3f UnRotateVector(const Vector3f& V) const;

}; 

FORCEINLINE float Rotator::ClampAxis(float Angle)
{
	Angle = Math::Fmod(Angle, 360.f);
	if (Angle < 0.f)
	{
		Angle += 360.f;
	}
	return Angle;
}

FORCEINLINE float Rotator::NormalizeAxis(float Angle)
{
	Angle = ClampAxis(Angle);
	if (Angle > 180.f)
	{
		Angle -= 360.f;
	}
	return Angle;
}

};
