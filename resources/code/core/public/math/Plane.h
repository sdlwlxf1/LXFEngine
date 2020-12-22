// this tiny game engine is rendered by CPU. I make it for learing the pipline and rasterization. This is the first step to realise my own game engine. so just do it!

#pragma once

#include "platform/Platform.h"
#include "MathUtility.h"
#include "Vector.h"

namespace LXF {

class Plane : public Vector3f
{
public:
	float W;

	Plane() : W(0) { }

	Plane(float InX, float InY, float InZ, float InW) : Vector3f(InX, InY, InZ), W(InW) { }
	Plane(Vector3f InNormal, float InW) : Vector3f(InNormal), W(InW) {}
	Plane(Vector3f A, Vector3f B, Vector3f C) : Vector3f(((B - A) ^ (C - A)).GetSafeNormal()) { W = A | (Vector3f)(*this); }


	Plane operator+(const Plane& V) const { return Plane(X + V.X, Y + V.Y, Z + V.Z, W + V.W); }
	Plane operator-(const Plane& V) const { return Plane(X - V.X, Y - V.Y, Z - V.Z, W - V.W); }
	Plane operator*(float Scale) const { return Plane(X * Scale, Y * Scale, Z * Scale, W * Scale); }
	Plane operator/(float Scale) const { const float RScale = 1.0f / Scale; return Plane(X * RScale, Y * RScale, Z * RScale, W * RScale); }
	bool operator==(const Plane& V) const { return X == V.X && Y == V.Y && Z == V.Z && W == V.W; }
	bool operator!=(const Plane& V) const { return X != V.X || Y != V.Y || Z != V.Z || W != V.W; }
	Plane operator-() const { return Plane(-X, -Y, -Z, -W); }

	Plane operator+=(const Plane& V) { X += V.X; Y += V.Y; Z += V.Z; W += V.W; return *this; }
	Plane operator-=(const Plane& V) { X -= V.X; Y -= V.Y; Z -= V.Z; W -= V.W; return *this; }
	Plane operator*=(float Scale) { X *= Scale; Y *= Scale; Z *= Scale; W *= Scale; return *this; }
	Plane operator/=(float Scale) { const float RScale = 1.0f / Scale; X *= RScale; Y *= RScale; Z *= RScale; W *= RScale; return *this; }
	Plane operator*=(const Plane& V) { X *= V.X; Y *= V.Y; Z *= V.Z; W *= V.W; return *this; }
	Plane operator/=(const Plane& V) { X /= V.X; Y /= V.Y; Z /= V.Z; W /= V.W; return *this; }

	bool Equals(const Plane& V, float Tolerance) const { return Math::Abs(X - V.X) <= Tolerance && Math::Abs(Y - V.Y) <= Tolerance && Math::Abs(Z - V.Z) <= Tolerance && Math::Abs(W - V.W) <= Tolerance; }
}; 

};
