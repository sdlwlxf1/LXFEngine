// this tiny game engine is rendered by CPU. I make it for learing the pipline and rasterization. This is the first step to realise my own game engine. so just do it!

#pragma once

#include "RotationTranslationMatrix.h"

namespace LXF {

/** Rotation matrix no translation */
class RotationMatrix
	: public RotationTranslationMatrix
{
public:

	/**
	 * Constructor.
	 *
	 * @param Rot rotation
	 */
	RotationMatrix(const Rotator& Rot)
		: RotationTranslationMatrix(Rot, Vector3f::ZeroVector)
	{ }

	/** Matrix factory. Return an Matrix so we don't have type conversion issues in expressions. */
	static Matrix Make(Rotator const& Rot)
	{
		return RotationMatrix(Rot);
	}

	/** Matrix factory. Return an Matrix so we don't have type conversion issues in expressions. */
	static Matrix Make(Quat const& Rot);

	/** Builds a rotation matrix given only a XAxis. Y and Z are unspecified but will be orthonormal. XAxis need not be normalized. */
	static Matrix MakeFromX(Vector3f const& XAxis);

	/** Builds a rotation matrix given only a YAxis. X and Z are unspecified but will be orthonormal. YAxis need not be normalized. */
	static Matrix MakeFromY(Vector3f const& YAxis);

	/** Builds a rotation matrix given only a ZAxis. X and Y are unspecified but will be orthonormal. ZAxis need not be normalized. */
	static Matrix MakeFromZ(Vector3f const& ZAxis);

	/** Builds a matrix with given X and Y axes. X will remain fixed, Y may be changed minimally to enforce orthogonality. Z will be computed. Inputs need not be normalized. */
	static Matrix MakeFromXY(Vector3f const& XAxis, Vector3f const& YAxis);

	/** Builds a matrix with given X and Z axes. X will remain fixed, Z may be changed minimally to enforce orthogonality. Y will be computed. Inputs need not be normalized. */
	static Matrix MakeFromXZ(Vector3f const& XAxis, Vector3f const& ZAxis);

	/** Builds a matrix with given Y and X axes. Y will remain fixed, X may be changed minimally to enforce orthogonality. Z will be computed. Inputs need not be normalized. */
	static Matrix MakeFromYX(Vector3f const& YAxis, Vector3f const& XAxis);

	/** Builds a matrix with given Y and Z axes. Y will remain fixed, Z may be changed minimally to enforce orthogonality. X will be computed. Inputs need not be normalized. */
	static Matrix MakeFromYZ(Vector3f const& YAxis, Vector3f const& ZAxis);

	/** Builds a matrix with given Z and X axes. Z will remain fixed, X may be changed minimally to enforce orthogonality. Y will be computed. Inputs need not be normalized. */
	static Matrix MakeFromZX(Vector3f const& ZAxis, Vector3f const& XAxis);

	/** Builds a matrix with given Z and Y axes. Z will remain fixed, Y may be changed minimally to enforce orthogonality. X will be computed. Inputs need not be normalized. */
	static Matrix MakeFromZY(Vector3f const& ZAxis, Vector3f const& YAxis);
};
};
