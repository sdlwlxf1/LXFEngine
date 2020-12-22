// Copyright 1998-2019 Epic Games, Inc. All Rights Reserved.

#pragma once

#include "CoreTypes.h"
#include "Vector.h"
#include "VectorRegister.h"
#include "Rotator.h"
#include "Vector4.h"
#include "Matrix.h"
#include "Quat.h"
#include "ScalarRegister.h"

namespace LXF {

#if ENABLE_VECTORIZED_TRANSFORM

/**
 * Transform composed of Scale, Rotation (as a quaternion), and Translation.
 *
 * Transforms can be used to convert from one space to another, for example by transforming
 * positions and directions from local space to world space.
 *
 * Transformation of position vectors is applied in the order:  Scale -> Rotate -> Translate.
 * Transformation of direction vectors is applied in the order: Scale -> Rotate.
 *
 * Order matters when composing transforms: C = A * B will yield a transform C that logically
 * first applies A then B to any subsequent transformation. Note that this is the opposite order of quaternion (Quat) multiplication.
 *
 * Example: LocalToWorld = (DeltaRotation * LocalToWorld) will change rotation in local space by DeltaRotation.
 * Example: LocalToWorld = (LocalToWorld * DeltaRotation) will change rotation in world space by DeltaRotation.
 */

class Transform
{

protected:
	/** Rotation of this transformation, as a quaternion */
	VectorRegister	Rotation;
	/** Translation of this transformation, as a vector */
	VectorRegister	Translation;
	/** 3D scale (always applied in local space) as a vector */
	VectorRegister Scale3D;	
public:
	/**
	 * The identity transformation (Rotation = Quat::Identity, Translation = Vector3f::ZeroVector, Scale3D = (1,1,1))
	 */
	static const Transform Identity;

#if ENABLE_NAN_DIAGNOSTIC
	FORCEINLINE void CheckNaN_Scale3D() const
	{
		if (VectorContainsNaNOrInfinite(Scale3D))
		{
			const_cast<Transform*>(this)->Scale3D = VectorSet_W0( VectorOne() );
		}
	}

	FORCEINLINE void CheckNaN_Translate() const
	{
		if (VectorContainsNaNOrInfinite(Translation))
		{
			const_cast<Transform*>(this)->Translation = VectorZero();
		}
	}

	FORCEINLINE void CheckNaN_Rotate() const
	{
		if (VectorContainsNaNOrInfinite(Rotation))
		{
			const_cast<Transform*>(this)->Rotation = VectorSet_W1( VectorZero() );
		}
	}

	FORCEINLINE void CheckNaN_All() const
	{
		CheckNaN_Scale3D();
		CheckNaN_Rotate();
		CheckNaN_Translate();
	}

	FORCEINLINE void Check_IsValid() const
	{
		CheckNaN_All();
		if (!IsValid())
		{
		}
		
	}
#else
	FORCEINLINE void CheckNaN_Translate() const {}
	FORCEINLINE void CheckNaN_Rotate() const {}
	FORCEINLINE void CheckNaN_Scale3D() const {}
	FORCEINLINE void CheckNaN_All() const {}
	FORCEINLINE void Check_IsValid() const {}
#endif

	/**
	 * Constructor with initialization to the identity transform.
	 */
	FORCEINLINE Transform()
	{
		// Rotation = {0,0,0,1)
		Rotation = VectorSet_W1( VectorZero() );
		// Translation = {0,0,0,0)
		Translation = VectorZero();
		// Scale3D = {1,1,1,0);
		Scale3D = VectorSet_W0( VectorOne() );
	}

	/**
	 * Constructor with an initial translation
	 *
	 * @param InTranslation The value to use for the translation component
	 */
	FORCEINLINE explicit Transform(const Vector3f& InTranslation) 
	{
		// Rotation = {0,0,0,1) quaternion identity
		Rotation =  VectorSet_W1( VectorZero() );
		//Translation = InTranslation;
		Translation = MakeVectorRegister(InTranslation.X, InTranslation.Y, InTranslation.Z, 0.0f );
		// Scale3D = {1,1,1,0);
		Scale3D = VectorSet_W0( VectorOne() );

		CheckNaN_All();
	}

	/**
	 * Constructor with an initial rotation
	 *
	 * @param InRotation The value to use for rotation component
	 */
	FORCEINLINE explicit Transform(const Quat& InRotation) 
	{
		// Rotation = InRotation
		Rotation =  VectorLoadAligned( &InRotation.X );
		// Translation = {0,0,0,0)
		Translation = VectorZero();
		// Scale3D = {1,1,1,0);
		Scale3D = VectorSet_W0( VectorOne() );

		CheckNaN_All();
	}

	/**
	 * Constructor with an initial rotation
	 *
	 * @param InRotation The value to use for rotation component  (after being converted to a quaternion)
	 */
	FORCEINLINE explicit Transform(const Rotator& InRotation) 
	{
		Quat InQuatRotation = InRotation.Quaternion();
		// Rotation = InRotation
		Rotation =  VectorLoadAligned( &InQuatRotation.X );
		// Translation = {0,0,0,0)
		Translation = VectorZero();
		// Scale3D = {1,1,1,0);
		Scale3D = VectorSet_W0( VectorOne() );

		CheckNaN_All();
	}

	/**
	 * Constructor with all components initialized
	 *
	 * @param InRotation The value to use for rotation component
	 * @param InTranslation The value to use for the translation component
	 * @param InScale3D The value to use for the scale component
	 */
	FORCEINLINE Transform(const Quat& InRotation, const Vector3f& InTranslation, const Vector3f& InScale3D = Vector3f::OneVector)
	{
		// Rotation = InRotation
		Rotation =  VectorLoadAligned( &InRotation.X );
		// Translation = InTranslation
		Translation = MakeVectorRegister(InTranslation.X, InTranslation.Y, InTranslation.Z, 0.0f );
		// Scale3D = InScale3D
		Scale3D = MakeVectorRegister(InScale3D.X, InScale3D.Y, InScale3D.Z, 0.0f );

		CheckNaN_All();
	}

	/**
	 * Constructor with all components initialized as VectorRegisters
	 *
	 * @param InRotation The value to use for rotation component
	 * @param InTranslation The value to use for the translation component
	 * @param InScale3D The value to use for the scale component
	 */
	FORCEINLINE Transform(const VectorRegister& InRotation, const VectorRegister& InTranslation, const VectorRegister& InScale3D) 
		: Rotation(InRotation),
		Translation(InTranslation),
		Scale3D(InScale3D)
	{
		CheckNaN_All();
	}

	/**
	 * Constructor with all components initialized, taking a Rotator as the rotation component
	 *
	 * @param InRotation The value to use for rotation component (after being converted to a quaternion)
	 * @param InTranslation The value to use for the translation component
	 * @param InScale3D The value to use for the scale component
	 */
	FORCEINLINE Transform(const Rotator& InRotation, const Vector3f& InTranslation, const Vector3f& InScale3D = Vector3f::OneVector)
	{
		Quat InQuatRotation = InRotation.Quaternion();
		// Rotation = InRotation
		Rotation =  VectorLoadAligned( &InQuatRotation.X );
		// Translation = InTranslation
		Translation = MakeVectorRegister(InTranslation.X, InTranslation.Y, InTranslation.Z, 0.0f );
		// Scale3D = InScale3D
		Scale3D = MakeVectorRegister(InScale3D.X, InScale3D.Y, InScale3D.Z, 0.0f );

		CheckNaN_All();
	}

	/**
	 * Constructor with leaving uninitialized memory
	 */
	FORCEINLINE explicit Transform(ENoInit) 
	{
		// Note: This can be used to track down initialization issues with bone transform arrays; but it will
		// cause issues with transient fields such as RootMotionDelta that get initialized to 0 by default
#if ENABLE_NAN_DIAGNOSTIC
		float qnan = Math::Log2(-5.3f);
		check(Math::IsNaN(qnan));
		Translation = MakeVectorRegister(qnan, qnan, qnan, qnan);
		Rotation = MakeVectorRegister(qnan, qnan, qnan, qnan);
		Scale3D = MakeVectorRegister(qnan, qnan, qnan, qnan);
#endif
	}

	/**
	 * Constructor for converting a Matrix (including scale) into a Transform.
	 */
	FORCEINLINE explicit Transform(const Matrix& InMatrix)
	{
		SetFromMatrix(InMatrix);
		CheckNaN_All();
	}

	/** Constructor that takes basis axes and translation */
	FORCEINLINE Transform(const Vector3f& InX, const Vector3f& InY, const Vector3f& InZ, const Vector3f& InTranslation)
	{
		SetFromMatrix(Matrix(InX, InY, InZ, InTranslation));
		CheckNaN_All();
	}

	/**
	* Copy another Transform into this one
	*/
	FORCEINLINE Transform& operator=(const Transform& Other)
	{
		this->Rotation = Other.Rotation;
		this->Translation = Other.Translation;
		this->Scale3D = Other.Scale3D;

		return *this;
	}


	FORCEINLINE Matrix ToMatrixWithScale() const
	{

		Matrix OutMatrix;
		VectorRegister DiagonalsXYZ;
		VectorRegister Adds;
		VectorRegister Subtracts;

		ToMatrixInternal( DiagonalsXYZ, Adds, Subtracts );
		const VectorRegister DiagonalsXYZ_W0 = VectorSet_W0(DiagonalsXYZ);

		// OutMatrix.M[0][0] = (1.0f - (yy2 + zz2)) * Scale.X;    // Diagonal.X
		// OutMatrix.M[0][1] = (xy2 + wz2) * Scale.X;             // Adds.X
		// OutMatrix.M[0][2] = (xz2 - wy2) * Scale.X;             // Subtracts.Z
		// OutMatrix.M[0][3] = 0.0f;                              // DiagonalsXYZ_W0.W
		const VectorRegister AddX_DC_DiagX_DC = VectorShuffle(Adds, DiagonalsXYZ_W0, 0, 0, 0, 0);
		const VectorRegister SubZ_DC_DiagW_DC = VectorShuffle(Subtracts, DiagonalsXYZ_W0, 2, 0, 3, 0);
		const VectorRegister Row0 = VectorShuffle(AddX_DC_DiagX_DC, SubZ_DC_DiagW_DC, 2, 0, 0, 2);

		// OutMatrix.M[1][0] = (xy2 - wz2) * Scale.Y;             // Subtracts.X
		// OutMatrix.M[1][1] = (1.0f - (xx2 + zz2)) * Scale.Y;    // Diagonal.Y
		// OutMatrix.M[1][2] = (yz2 + wx2) * Scale.Y;             // Adds.Y
		// OutMatrix.M[1][3] = 0.0f;                            // DiagonalsXYZ_W0.W
		const VectorRegister SubX_DC_DiagY_DC = VectorShuffle(Subtracts, DiagonalsXYZ_W0, 0, 0, 1, 0);
		const VectorRegister AddY_DC_DiagW_DC = VectorShuffle(Adds, DiagonalsXYZ_W0, 1, 0, 3, 0);
		const VectorRegister Row1 = VectorShuffle(SubX_DC_DiagY_DC, AddY_DC_DiagW_DC, 0, 2, 0, 2);

		// OutMatrix.M[2][0] = (xz2 + wy2) * Scale.Z;             // Adds.Z
		// OutMatrix.M[2][1] = (yz2 - wx2) * Scale.Z;             // Subtracts.Y
		// OutMatrix.M[2][2] = (1.0f - (xx2 + yy2)) * Scale.Z;    // Diagonals.Z
		// OutMatrix.M[2][3] = 0.0f;                              // DiagonalsXYZ_W0.W
		const VectorRegister AddZ_DC_SubY_DC = VectorShuffle(Adds, Subtracts, 2, 0, 1, 0);
		const VectorRegister Row2 = VectorShuffle(AddZ_DC_SubY_DC, DiagonalsXYZ_W0, 0, 2, 2, 3);

		VectorStoreAligned(Row0, &(OutMatrix.M[0][0]));
		VectorStoreAligned(Row1, &(OutMatrix.M[1][0]));
		VectorStoreAligned(Row2, &(OutMatrix.M[2][0]));

		// OutMatrix.M[3][0] = Translation.X;
		// OutMatrix.M[3][1] = Translation.Y;
		// OutMatrix.M[3][2] = Translation.Z;
		// OutMatrix.M[3][3] = 1.0f;
		const VectorRegister Row3 = VectorSet_W1(Translation);
		VectorStoreAligned(Row3, &(OutMatrix.M[3][0]));

		return OutMatrix;
	}

	/**
	* Convert this Transform to matrix with scaling and compute the inverse of that.
	*/
	FORCEINLINE Matrix ToInverseMatrixWithScale() const
	{
		// todo: optimize
		return ToMatrixWithScale().Inverse();
	}

	/**
	* Convert this Transform to inverse.
	*/
	FORCEINLINE Transform Inverse() const
	{
		// Replacement of Inverse of Matrix
		if (VectorAnyGreaterThan(VectorAbs(Scale3D), GlobalVectorConstants::SmallNumber))
		{
			return InverseFast();
		}
		else
		{
			return Transform::Identity;
		}
	}

	/**
	* Convert this Transform to a transformation matrix, ignoring its scaling
	*/
	FORCEINLINE Matrix ToMatrixNoScale() const
	{
		Matrix OutMatrix;
		VectorRegister DiagonalsXYZ;
		VectorRegister Adds;
		VectorRegister Subtracts;

		ToMatrixInternalNoScale( DiagonalsXYZ, Adds, Subtracts );
		const VectorRegister DiagonalsXYZ_W0 = VectorSet_W0(DiagonalsXYZ);

		// OutMatrix.M[0][0] = (1.0f - (yy2 + zz2));			// Diagonal.X
		// OutMatrix.M[0][1] = (xy2 + wz2);						// Adds.X
		// OutMatrix.M[0][2] = (xz2 - wy2);						// Subtracts.Z
		// OutMatrix.M[0][3] = 0.0f;                            // DiagonalsXYZ_W0.W
		const VectorRegister AddX_DC_DiagX_DC = VectorShuffle(Adds, DiagonalsXYZ_W0, 0, 0, 0, 0);
		const VectorRegister SubZ_DC_DiagW_DC = VectorShuffle(Subtracts, DiagonalsXYZ_W0, 2, 0, 3, 0);
		const VectorRegister Row0 = VectorShuffle(AddX_DC_DiagX_DC, SubZ_DC_DiagW_DC, 2, 0, 0, 2);

		// OutMatrix.M[1][0] = (xy2 - wz2);			            // Subtracts.X
		// OutMatrix.M[1][1] = (1.0f - (xx2 + zz2));		    // Diagonal.Y
		// OutMatrix.M[1][2] = (yz2 + wx2);						// Adds.Y
		// OutMatrix.M[1][3] = 0.0f;                            // DiagonalsXYZ_W0.W
		const VectorRegister SubX_DC_DiagY_DC = VectorShuffle(Subtracts, DiagonalsXYZ_W0, 0, 0, 1, 0);
		const VectorRegister AddY_DC_DiagW_DC = VectorShuffle(Adds, DiagonalsXYZ_W0, 1, 0, 3, 0);
		const VectorRegister Row1 = VectorShuffle(SubX_DC_DiagY_DC, AddY_DC_DiagW_DC, 0, 2, 0, 2);

		// OutMatrix.M[2][0] = (xz2 + wy2);						// Adds.Z
		// OutMatrix.M[2][1] = (yz2 - wx2);						// Subtracts.Y
		// OutMatrix.M[2][2] = (1.0f - (xx2 + yy2));		    // Diagonals.Z
		// OutMatrix.M[2][3] = 0.0f;                            // DiagonalsXYZ_W0.W
		const VectorRegister AddZ_DC_SubY_DC = VectorShuffle(Adds, Subtracts, 2, 0, 1, 0);
		const VectorRegister Row2 = VectorShuffle(AddZ_DC_SubY_DC, DiagonalsXYZ_W0, 0, 2, 2, 3);

		VectorStoreAligned(Row0, &(OutMatrix.M[0][0]));
		VectorStoreAligned(Row1, &(OutMatrix.M[1][0]));
		VectorStoreAligned(Row2, &(OutMatrix.M[2][0]));

		// OutMatrix.M[3][0] = Translation.X;
		// OutMatrix.M[3][1] = Translation.Y;
		// OutMatrix.M[3][2] = Translation.Z;
		// OutMatrix.M[3][3] = 1.0f;
		const VectorRegister Row3 = VectorSet_W1(Translation);
		VectorStoreAligned(Row3, &(OutMatrix.M[3][0]));

		return OutMatrix;
	}

	/** Set this transform to the weighted blend of the supplied two transforms. */
	FORCEINLINE void Blend(const Transform& Atom1, const Transform& Atom2, float Alpha)
	{
#if !(BUILD_SHIPPING)
		// Check that all bone atoms coming from animation are normalized
		check( Atom1.IsRotationNormalized() );
		check( Atom2.IsRotationNormalized() );
#endif

		if( Math::Abs(Alpha) <= ZERO_ANIMWEIGHT_THRESH )
		{
			// if blend is all the way for child1, then just copy its bone atoms
			(*this) = Atom1;
		}
		else if( Math::Abs(Alpha - 1.0f) <= ZERO_ANIMWEIGHT_THRESH )
		{
			// if blend is all the way for child2, then just copy its bone atoms
			(*this) = Atom2;
		}
		else
		{
			// Simple linear interpolation for translation and scale.			
			ScalarRegister BlendWeight = ScalarRegister(Alpha);

			Translation = Math::Lerp(Atom1.Translation, Atom2.Translation, BlendWeight.Value);			
			Scale3D = Math::Lerp(Atom1.Scale3D, Atom2.Scale3D, BlendWeight.Value);			

			VectorRegister VRotation = VectorLerpQuat(Atom1.Rotation, Atom2.Rotation, BlendWeight.Value);

			// ..and renormalize
			Rotation = VectorNormalizeQuaternion(VRotation);

			CheckNaN_All(); // MR
		}
	}

	/** Set this Transform to the weighted blend of it and the supplied Transform. */
	FORCEINLINE void BlendWith(const Transform& OtherAtom, float Alpha)
	{
#if !(BUILD_SHIPPING)
		// Check that all bone atoms coming from animation are normalized
		check( IsRotationNormalized() );
		check( OtherAtom.IsRotationNormalized() );
#endif

		if( Alpha > ZERO_ANIMWEIGHT_THRESH )
		{
			if( Alpha >= 1.f - ZERO_ANIMWEIGHT_THRESH )
			{
				// if blend is all the way for child2, then just copy its bone atoms
				(*this) = OtherAtom;
			}
			else 
			{
				// Simple linear interpolation for translation and scale.				
				ScalarRegister BlendWeight = ScalarRegister(Alpha);
				Translation = Math::Lerp(Translation, OtherAtom.Translation, BlendWeight.Value);
				
				Scale3D = Math::Lerp(Scale3D, OtherAtom.Scale3D, BlendWeight.Value);
				
				VectorRegister VRotation = VectorLerpQuat(Rotation, OtherAtom.Rotation, BlendWeight.Value);
				
				// ..and renormalize
				Rotation = VectorNormalizeQuaternion(VRotation);

				CheckNaN_All(); 
			}
		}
	}


	/**
	 * Quaternion addition is wrong here. This is just a special case for linear interpolation.
	 * Use only within blends!!
	 * Rotation part is NOT normalized!!
	 */
	FORCEINLINE Transform operator+(const Transform& Atom) const
	{
		return Transform( VectorAdd(Rotation, Atom.Rotation), VectorAdd(Translation, Atom.Translation ), VectorAdd( Scale3D, Atom.Scale3D ) );
	}

	FORCEINLINE Transform& operator+=(const Transform& Atom)
	{ 
		Translation = VectorAdd(Translation, Atom.Translation);
		Rotation = VectorAdd(Rotation, Atom.Rotation);
		Scale3D = VectorAdd(Scale3D, Atom.Scale3D);

		return *this;
	}

	FORCEINLINE Transform operator*(const ScalarRegister& Mult) const
	{		
		return Transform( VectorMultiply(Rotation, Mult), VectorMultiply(Translation, Mult), VectorMultiply(Scale3D, Mult) );
	}

	FORCEINLINE Transform& operator*=(const ScalarRegister& Mult)
	{			
		Translation= VectorMultiply(Translation, Mult);
		Rotation = VectorMultiply(Rotation, Mult);
		Scale3D = VectorMultiply(Scale3D, Mult);

		return *this;
	}

	FORCEINLINE Transform		operator*(const Transform& Other) const;
	FORCEINLINE void			operator*=(const Transform& Other);
	FORCEINLINE Transform		operator*(const Quat& Other) const;
	FORCEINLINE void			operator*=(const Quat& Other);

	FORCEINLINE static bool AnyHasNegativeScale(const Vector3f& InScale3D, const Vector3f& InOtherScale3D);
	FORCEINLINE void ScaleTranslation(const Vector3f& InScale3D);
	FORCEINLINE void ScaleTranslation(const float& Scale);
	FORCEINLINE void RemoveScaling(float Tolerance=SMALL_NUMBER);
	FORCEINLINE float GetMaximumAxisScale() const;
	FORCEINLINE float GetMinimumAxisScale() const;
	// Inverse does not work well with VQS format(in particular non-uniform), so removing it, but made two below functions to be used instead. 

	/*******************************************************************************************
	 * The below 2 functions are the ones to get delta transform and return Transform format that can be concatenated
	 * Inverse itself can't concatenate with VQS format(since VQS always transform from S->Q->T, where inverse happens from T(-1)->Q(-1)->S(-1))
	 * So these 2 provides ways to fix this
	 * GetRelativeTransform returns this*Other(-1) and parameter is Other(not Other(-1))
	 * GetRelativeTransformReverse returns this(-1)*Other, and parameter is Other. 
	 *******************************************************************************************/
	Transform GetRelativeTransform(const Transform& Other) const;
	Transform GetRelativeTransformReverse(const Transform& Other) const;
	/**
	 * Set current transform and the relative to ParentTransform.
	 * Equates to This = This->GetRelativeTransform(Parent), but saves the intermediate Transform storage and copy.
	 */
	void		SetToRelativeTransform(const Transform& ParentTransform);

	FORCEINLINE Vector4f	TransformVector4(const Vector4f& V) const;
	FORCEINLINE Vector4f	TransformVector4NoScale(const Vector4f& V) const;
	FORCEINLINE Vector3f	TransformPosition(const Vector3f& V) const;
	FORCEINLINE Vector3f	TransformPositionNoScale(const Vector3f& V) const;


	/** Inverts the transform and then transforms V - correctly handles scaling in this transform. */
	FORCEINLINE Vector3f		InverseTransformPosition(const Vector3f &V) const;
	FORCEINLINE Vector3f		InverseTransformPositionNoScale(const Vector3f &V) const;
	FORCEINLINE Vector3f		TransformVector(const Vector3f& V) const;
	FORCEINLINE Vector3f		TransformVectorNoScale(const Vector3f& V) const;

	/** 
	 *	Transform a direction vector by the inverse of this matrix - will not take into account translation part.
	 *	If you want to transform a surface normal (or plane) and correctly account for non-uniform scaling you should use TransformByUsingAdjointT with adjoint of matrix inverse.
	 */
	FORCEINLINE Vector3f InverseTransformVector(const Vector3f &V) const;
	FORCEINLINE Vector3f InverseTransformVectorNoScale(const Vector3f &V) const;

	/**
	* Transform a rotation.
	* For example if this is a LocalToWorld transform, TransformRotation(Q) would transform Q from local to world space.
	*/
	FORCEINLINE Quat TransformRotation(const Quat& Q) const;

	/**
	* Inverse transform a rotation.
	* For example if this is a LocalToWorld transform, InverseTransformRotation(Q) would transform Q from world to local space.
	*/
	FORCEINLINE Quat InverseTransformRotation(const Quat& Q) const;

	FORCEINLINE Transform	GetScaled(float Scale) const;
	FORCEINLINE Transform	GetScaled(Vector3f Scale) const;
	FORCEINLINE Vector3f		GetScaledAxis(EAxis::Type InAxis) const;
	FORCEINLINE Vector3f		GetUnitAxis(EAxis::Type InAxis) const;
	FORCEINLINE void		Mirror(EAxis::Type MirrorAxis, EAxis::Type FlipAxis);
	FORCEINLINE static Vector3f	GetSafeScaleReciprocal(const Vector3f& InScale, float Tolerance=SMALL_NUMBER);


	FORCEINLINE Vector3f GetLocation() const
	{
		return GetTranslation();
	}

	FORCEINLINE Rotator Rotator() const
	{
		Quat OutRotation;
		VectorStoreAligned(Rotation, &OutRotation);
		return OutRotation.Rotator();
	}

	/** Calculate the  */
	FORCEINLINE float GetDeterminant() const
	{
		//#todo - vectorized version of this
		Vector4f OutScale3D;
		VectorStoreAligned(Scale3D, &OutScale3D);
		return OutScale3D.X * OutScale3D.Y * OutScale3D.Z;
	}

	/** Set the translation of this transformation */
	FORCEINLINE void SetLocation(const Vector3f& Origin)
	{		
		Translation = VectorLoadFloat3_W0(&Origin);
		CheckNaN_Translate();
	}

	/**
	 * Checks the components for NaN's
	 * @return Returns true if any component (rotation, translation, or scale) is a NAN
	 */
	bool ContainsNaN() const
	{
		if (VectorContainsNaNOrInfinite(Rotation))
		{
			return true;
		}
		if (VectorContainsNaNOrInfinite(Translation))
		{
			return true;
		}

		if (VectorContainsNaNOrInfinite(Scale3D))
		{
			return true;
		}
		return false;
	}

	FORCEINLINE bool IsValid() const
	{
		if ( ContainsNaN() )
		{
			return false;
		}

		if ( !IsRotationNormalized() )
		{
			return false;
		}

		return true;
	}

private:

	FORCEINLINE static bool Private_AnyHasNegativeScale(const VectorRegister& InScale3D, const  VectorRegister& InOtherScale3D)
	{
		return !!VectorAnyLesserThan(VectorMin(InScale3D, InOtherScale3D), GlobalVectorConstants::FloatZero);
	}

	FORCEINLINE bool Private_RotationEquals( const VectorRegister& InRotation, const ScalarRegister& Tolerance = ScalarRegister(GlobalVectorConstants::KindaSmallNumber)) const
	{			
		// !( (Math::Abs(X-Q.X) > Tolerance) || (Math::Abs(Y-Q.Y) > Tolerance) || (Math::Abs(Z-Q.Z) > Tolerance) || (Math::Abs(W-Q.W) > Tolerance) )
		const VectorRegister RotationSub = VectorAbs(VectorSubtract(Rotation, InRotation));		
		// !( (Math::Abs(X+Q.X) > Tolerance) || (Math::Abs(Y+Q.Y) > Tolerance) || (Math::Abs(Z+Q.Z) > Tolerance) || (Math::Abs(W+Q.W) > Tolerance) )
		const VectorRegister RotationAdd = VectorAbs(VectorAdd(Rotation, InRotation));
		return !VectorAnyGreaterThan(RotationSub, Tolerance.Value) || !VectorAnyGreaterThan(RotationAdd, Tolerance.Value);
	}

	FORCEINLINE bool Private_TranslationEquals( const VectorRegister& InTranslation, const ScalarRegister& Tolerance = ScalarRegister(GlobalVectorConstants::KindaSmallNumber)) const
	{			
		// !( (Math::Abs(X-V.X) > Tolerance) || (Math::Abs(Y-V.Y) > Tolerance) || (Math::Abs(Z-V.Z) > Tolerance) )
		const VectorRegister TranslationDiff = VectorAbs(VectorSubtract(Translation, InTranslation));		
		return !VectorAnyGreaterThan(TranslationDiff, Tolerance.Value);
	}

	FORCEINLINE bool Private_Scale3DEquals( const VectorRegister& InScale3D, const ScalarRegister& Tolerance = ScalarRegister(GlobalVectorConstants::KindaSmallNumber)) const
	{
		// !( (Math::Abs(X-V.X) > Tolerance) || (Math::Abs(Y-V.Y) > Tolerance) || (Math::Abs(Z-V.Z) > Tolerance) )
		const VectorRegister ScaleDiff = VectorAbs(VectorSubtract(Scale3D, InScale3D));
		return !VectorAnyGreaterThan(ScaleDiff, Tolerance.Value);
	}

public:

	// Test if A's rotation equals B's rotation, within a tolerance. Preferred over "A.GetRotation().Equals(B.GetRotation())" because it is faster on some platforms.
	FORCEINLINE static bool AreRotationsEqual(const Transform& A, const Transform& B, float Tolerance=KINDA_SMALL_NUMBER)
	{
		return A.Private_RotationEquals(B.Rotation, ScalarRegister(Tolerance));
	}

	// Test if A's translation equals B's translation, within a tolerance. Preferred over "A.GetTranslation().Equals(B.GetTranslation())" because it avoids VectorRegister->Vector3f conversion.
	FORCEINLINE static bool AreTranslationsEqual(const Transform& A, const Transform& B, float Tolerance=KINDA_SMALL_NUMBER)
	{
		return A.Private_TranslationEquals(B.Translation, ScalarRegister(Tolerance));
	}

	// Test if A's scale equals B's scale, within a tolerance. Preferred over "A.GetScale3D().Equals(B.GetScale3D())" because it avoids VectorRegister->Vector3f conversion.
	FORCEINLINE static bool AreScale3DsEqual(const Transform& A, const Transform& B, float Tolerance=KINDA_SMALL_NUMBER)
	{
		return A.Private_Scale3DEquals(B.Scale3D, ScalarRegister(Tolerance));
	}


	// Test if this Transform's rotation equals another's rotation, within a tolerance. Preferred over "GetRotation().Equals(Other.GetRotation())" because it is faster on some platforms.
	FORCEINLINE bool RotationEquals(const Transform& Other, float Tolerance=KINDA_SMALL_NUMBER) const
	{
		return AreRotationsEqual(*this, Other, Tolerance);
	}

	// Test if this Transform's translation equals another's translation, within a tolerance. Preferred over "GetTranslation().Equals(Other.GetTranslation())" because it avoids VectorRegister->Vector3f conversion.
	FORCEINLINE bool TranslationEquals(const Transform& Other, float Tolerance=KINDA_SMALL_NUMBER) const
	{
		return AreTranslationsEqual(*this, Other, Tolerance);
	}

	// Test if this Transform's scale equals another's scale, within a tolerance. Preferred over "GetScale3D().Equals(Other.GetScale3D())" because it avoids VectorRegister->Vector3f conversion.
	FORCEINLINE bool Scale3DEquals(const Transform& Other, float Tolerance=KINDA_SMALL_NUMBER) const
	{
		return AreScale3DsEqual(*this, Other, Tolerance);
	}


	// Test if all components of the transforms are equal, within a tolerance.
	FORCEINLINE bool Equals(const Transform& Other, float Tolerance=KINDA_SMALL_NUMBER) const
	{
		const ScalarRegister ToleranceRegister(Tolerance);
		return Private_TranslationEquals(Other.Translation, ToleranceRegister) && Private_RotationEquals(Other.Rotation, ToleranceRegister) && Private_Scale3DEquals(Other.Scale3D, ToleranceRegister);
	}

	// Test if rotation and translation components of the transforms are equal, within a tolerance.
	FORCEINLINE bool EqualsNoScale(const Transform& Other, float Tolerance=KINDA_SMALL_NUMBER) const
	{
		const ScalarRegister ToleranceRegister(Tolerance);
		return Private_TranslationEquals(Other.Translation, ToleranceRegister) && Private_RotationEquals(Other.Rotation, ToleranceRegister);
	}

	FORCEINLINE static void Multiply(Transform* OutTransform, const Transform* A, const Transform* B);
	/**
	 * Sets the components
	 * @param InRotation The new value for the Rotation component
	 * @param InTranslation The new value for the Translation component
	 * @param InScale3D The new value for the Scale3D component
	 */
	FORCEINLINE void SetComponents(const Quat& InRotation, const Vector3f& InTranslation, const Vector3f& InScale3D) 
	{
		Rotation = VectorLoadAligned(&InRotation);
		Translation = VectorLoadFloat3_W0(&InTranslation);
		Scale3D = VectorLoadFloat3_W0(&InScale3D);

		CheckNaN_All();
	}

	/**
	 * Sets the components to the identity transform:
	 *   Rotation = (0,0,0,1)
	 *   Translation = (0,0,0)
	 *   Scale3D = (1,1,1)
	 */
	FORCEINLINE void SetIdentity()
	{
		// Rotation = {0,0,0,1)
		Rotation = VectorSet_W1( VectorZero() );
		// Translation = {0,0,0,0)
		Translation = VectorZero();
		// Scale3D = {1,1,1,0);
		Scale3D = VectorSet_W0 ( VectorOne() );
	}

	/**
	 * Scales the Scale3D component by a new factor
	 * @param Scale3DMultiplier The value to multiply Scale3D with
	 */
	FORCEINLINE void MultiplyScale3D(const Vector3f& Scale3DMultiplier)
	{
		Scale3D = VectorMultiply(Scale3D, VectorLoadFloat3_W0(&Scale3DMultiplier));
		CheckNaN_Scale3D();
	}

	/**
	 * Sets the translation component
	 * @param NewTranslation The new value for the translation component
	 */
	FORCEINLINE void SetTranslation(const Vector3f& NewTranslation)
	{
		Translation = VectorLoadFloat3_W0(&NewTranslation);
		CheckNaN_Translate();
	}

	/** Copy translation from another Transform. */
	FORCEINLINE void CopyTranslation(const Transform& Other)
	{
		Translation = Other.Translation;
	}

	/**
	 * Concatenates another rotation to this transformation 
	 * @param DeltaRotation The rotation to concatenate in the following fashion: Rotation = Rotation * DeltaRotation
	 */
	FORCEINLINE void ConcatenateRotation(const Quat& DeltaRotation)
	{		
		Rotation = VectorQuaternionMultiply2(Rotation, VectorLoadAligned(&DeltaRotation));		
		CheckNaN_Rotate();
	}

	/**
	 * Adjusts the translation component of this transformation 
	 * @param DeltaTranslation The translation to add in the following fashion: Translation += DeltaTranslation
	 */
	FORCEINLINE void AddToTranslation(const Vector3f& DeltaTranslation)
	{		
		Translation = VectorAdd(Translation, VectorLoadFloat3_W0(&DeltaTranslation));
		CheckNaN_Translate();
	}

	/**
	 * Add the translations from two Transforms and return the result.
	 * @return A.Translation + B.Translation
	 */
	FORCEINLINE static Vector3f AddTranslations(const Transform& A, const Transform& B)
	{
		Vector3f Result;
		VectorStoreFloat3(VectorAdd(A.Translation, B.Translation), &Result);
		return Result;
	}

	/**
	 * Subtract translations from two Transforms and return the difference.
	 * @return A.Translation - B.Translation.
	 */
	FORCEINLINE static Vector3f SubtractTranslations(const Transform& A, const Transform& B)
	{
		Vector3f Result;
		VectorStoreFloat3(VectorSubtract(A.Translation, B.Translation), &Result);
		return Result;
	}

	/**
	 * Sets the rotation component
	 * @param NewRotation The new value for the rotation component
	 */
	FORCEINLINE void SetRotation(const Quat& NewRotation)
	{
		Rotation = VectorLoadAligned(&NewRotation);
		CheckNaN_Rotate();
	}

	/** Copy rotation from another Transform. */
	FORCEINLINE void CopyRotation(const Transform& Other)
	{
		Rotation = Other.Rotation;
	}

	/**
	 * Sets the Scale3D component
	 * @param NewScale3D The new value for the Scale3D component
	 */
	FORCEINLINE void SetScale3D(const Vector3f& NewScale3D)
	{
		Scale3D = VectorLoadFloat3_W0(&NewScale3D);
		CheckNaN_Scale3D();
	}

	/** Copy scale from another Transform. */
	FORCEINLINE void CopyScale3D(const Transform& Other)
	{
		Scale3D = Other.Scale3D;
	}

	/**
	 * Sets both the translation and Scale3D components at the same time
	 * @param NewTranslation The new value for the translation component
	 * @param NewScale3D The new value for the Scale3D component
	 */
	FORCEINLINE void SetTranslationAndScale3D(const Vector3f& NewTranslation, const Vector3f& NewScale3D)
	{
		Translation = VectorLoadFloat3_W0(&NewTranslation);
		Scale3D = VectorLoadFloat3_W0(&NewScale3D);

		CheckNaN_Translate();
		CheckNaN_Scale3D();
	}
	
	/** @note : Added template type function for Accumulate
	  * The template type isn't much useful yet, but it is with the plan to move forward
	  * to unify blending features with just type of additive or full pose
	  * Eventually it would be nice to just call blend and it all works depending on full pose
	  * or additive, but right now that is a lot more refactoring
	  * For now this types only defines the different functionality of accumulate
	  */

	/**
	* Accumulates another transform with this one
	*
	* Rotation is accumulated multiplicatively (Rotation = SourceAtom.Rotation * Rotation)
	* Translation is accumulated additively (Translation += SourceAtom.Translation)
	* Scale3D is accumulated multiplicatively (Scale3D *= SourceAtom.Scale3D)
	*
	* @param SourceAtom The other transform to accumulate into this one
	*/
	FORCEINLINE void Accumulate(const Transform& SourceAtom)
	{
		const VectorRegister BlendedRotation = SourceAtom.Rotation;
		const VectorRegister RotationW = VectorReplicate(BlendedRotation, 3);

		// if( Square(SourceAtom.Rotation.W) < 1.f - DELTA * DELTA )
		if (VectorAnyGreaterThan(GlobalVectorConstants::RotationSignificantThreshold, VectorMultiply(RotationW, RotationW)))
		{
			// Rotation = SourceAtom.Rotation * Rotation;
			Rotation = VectorQuaternionMultiply2(BlendedRotation, Rotation);
		}

		// Translation += SourceAtom.Translation;
		// Scale *= SourceAtom.Scale;
		Translation = VectorAdd(Translation, SourceAtom.Translation);
		Scale3D = VectorMultiply(Scale3D, SourceAtom.Scale3D);

		CheckNaN_All();

		checkSlow(IsRotationNormalized());
	}

	/**
	* Accumulates another transform with this one, with a blending weight
	*
	* Let SourceAtom = Atom * BlendWeight
	* Rotation is accumulated multiplicatively (Rotation = SourceAtom.Rotation * Rotation).
	* Translation is accumulated additively (Translation += SourceAtom.Translation)
	* Scale3D is accumulated multiplicatively (Scale3D *= SourceAtom.Scale3D)
	*
	* Note: Rotation will not be normalized! Will have to be done manually.
	*
	* @param Atom The other transform to accumulate into this one
	* @param BlendWeight The weight to multiply Atom by before it is accumulated.
	*/
	FORCEINLINE void Accumulate(const Transform& Atom, const ScalarRegister& BlendWeight)
	{
		// SourceAtom = Atom * BlendWeight;
		const VectorRegister BlendedRotation = VectorMultiply(Atom.Rotation, BlendWeight.Value);
		const VectorRegister BlendedTranslation = VectorMultiply(Atom.Translation, BlendWeight.Value);
		const VectorRegister BlendedScale = VectorMultiply(Atom.Scale3D, BlendWeight.Value);

		const VectorRegister RotationW = VectorReplicate(BlendedRotation, 3);

		// Add ref pose relative animation to base animation, only if rotation is significant.
		// if( Square(SourceAtom.Rotation.W) < 1.f - DELTA * DELTA )
		if (VectorAnyGreaterThan(GlobalVectorConstants::RotationSignificantThreshold, VectorMultiply(RotationW, RotationW)))
		{
			// Rotation = SourceAtom.Rotation * Rotation;
			Rotation = VectorQuaternionMultiply2(BlendedRotation, Rotation);
		}

		// Translation += SourceAtom.Translation;
		// Scale *= SourceAtom.Scale;
		Translation = VectorAdd(Translation, BlendedTranslation);
		Scale3D = VectorMultiply(Scale3D, BlendedScale);

		CheckNaN_All();
	}
	/**
	 * Accumulates another transform with this one, with an optional blending weight
	 *
	 * Rotation is accumulated additively, in the shortest direction (Rotation = Rotation +/- DeltaAtom.Rotation * Weight)
	 * Translation is accumulated additively (Translation += DeltaAtom.Translation * Weight)
	 * Scale3D is accumulated additively (Scale3D += DeltaAtom.Scale * Weight)
	 *
	 * @param DeltaAtom The other transform to accumulate into this one
	 * @param Weight The weight to multiply DeltaAtom by before it is accumulated.
	 */
	FORCEINLINE void AccumulateWithShortestRotation(const Transform& DeltaAtom, const ScalarRegister& BlendWeight)
	{
		const VectorRegister BlendedRotation = VectorMultiply(DeltaAtom.Rotation, BlendWeight.Value);

		Rotation = VectorAccumulateQuaternionShortestPath(Rotation, BlendedRotation);

		Translation = VectorMultiplyAdd(DeltaAtom.Translation, BlendWeight, Translation);
		Scale3D = VectorMultiplyAdd(DeltaAtom.Scale3D, BlendWeight, Scale3D);

		CheckNaN_All();
	}

	/** Accumulates another transform with this one, with a blending weight
	*
	* Let SourceAtom = Atom * BlendWeight
	* Rotation is accumulated multiplicatively (Rotation = SourceAtom.Rotation * Rotation).
	* Translation is accumulated additively (Translation += SourceAtom.Translation)
	* Scale3D is accumulated assuming incoming scale is additive scale (Scale3D *= (1 + SourceAtom.Scale3D))
	*
	* When we create additive, we create additive scale based on [TargetScale/SourceScale -1]
	* because that way when you apply weight of 0.3, you don't shrink. We only saves the % of grow/shrink
	* when we apply that back to it, we add back the 1, so that it goes back to it.
	* This solves issue where you blend two additives with 0.3, you don't come back to 0.6 scale, but 1 scale at the end
	* because [1 + [1-1]*0.3 + [1-1]*0.3] becomes 1, so you don't shrink by applying additive scale
	*
	* Note: Rotation will not be normalized! Will have to be done manually.
	*
	* @param Atom The other transform to accumulate into this one
	* @param BlendWeight The weight to multiply Atom by before it is accumulated.
	*/
	FORCEINLINE void AccumulateWithAdditiveScale(const Transform& Atom, const ScalarRegister& BlendWeight)
	{
		const VectorRegister DefaultScale = MakeVectorRegister(1.f, 1.f, 1.f, 0.f);

		// SourceAtom = Atom * BlendWeight;
		const VectorRegister BlendedRotation = VectorMultiply(Atom.Rotation, BlendWeight.Value);
		const VectorRegister BlendedScale = VectorMultiply(Atom.Scale3D, BlendWeight.Value);
		const VectorRegister BlendedTranslation = VectorMultiply(Atom.Translation, BlendWeight.Value);

		const VectorRegister RotationW = VectorReplicate(BlendedRotation, 3);

		// Add ref pose relative animation to base animation, only if rotation is significant.
		// if( Square(SourceAtom.Rotation.W) < 1.f - DELTA * DELTA )
		if (VectorAnyGreaterThan(GlobalVectorConstants::RotationSignificantThreshold, VectorMultiply(RotationW, RotationW)))
		{
			// Rotation = SourceAtom.Rotation * Rotation;
			Rotation = VectorQuaternionMultiply2(BlendedRotation, Rotation);
		}

		// Translation += SourceAtom.Translation;
		// Scale *= SourceAtom.Scale;
		Translation = VectorAdd(Translation, BlendedTranslation);
		Scale3D = VectorMultiply(Scale3D, VectorAdd(DefaultScale, BlendedScale));

		CheckNaN_All();
	}

	/**
	 * Set the translation and Scale3D components of this transform to a linearly interpolated combination of two other transforms
	 *
	 * Translation = Math::Lerp(SourceAtom1.Translation, SourceAtom2.Translation, Alpha)
	 * Scale3D = Math::Lerp(SourceAtom1.Scale3D, SourceAtom2.Scale3D, Alpha)
	 *
	 * @param SourceAtom1 The starting point source atom (used 100% if Alpha is 0)
	 * @param SourceAtom2 The ending point source atom (used 100% if Alpha is 1)
	 * @param Alpha The blending weight between SourceAtom1 and SourceAtom2
	 */
	FORCEINLINE void LerpTranslationScale3D(const Transform& SourceAtom1, const Transform& SourceAtom2, const ScalarRegister& Alpha)
	{
		Translation	= Math::Lerp(SourceAtom1.Translation, SourceAtom2.Translation, Alpha.Value);
		Scale3D = Math::Lerp(SourceAtom1.Scale3D, SourceAtom2.Scale3D, Alpha.Value);

		CheckNaN_Translate();
		CheckNaN_Scale3D();
	}

	/**
	 * Normalize the rotation component of this transformation
	 */
	FORCEINLINE void NormalizeRotation()
	{
		Rotation = VectorNormalizeQuaternion(Rotation);
		CheckNaN_Rotate();
	}

	/**
	 * Checks whether the rotation component is normalized or not
	 *
	 * @return true if the rotation component is normalized, and false otherwise.
	 */
	FORCEINLINE bool IsRotationNormalized() const
	{		
		const VectorRegister TestValue = VectorAbs(VectorSubtract(VectorOne(), VectorDot4(Rotation, Rotation)));
		return !VectorAnyGreaterThan(TestValue, GlobalVectorConstants::ThreshQuatNormalized);
	}

	/**
	 * Blends the Identity transform with a weighted source transform and accumulates that into a destination transform
	 *
	 * SourceAtom = Blend(Identity, SourceAtom, BlendWeight)
	 * FinalAtom.Rotation = SourceAtom.Rotation * FinalAtom.Rotation
	 * FinalAtom.Translation += SourceAtom.Translation
	 * FinalAtom.Scale3D *= SourceAtom.Scale3D
	 *
	 * @param FinalAtom [in/out] The atom to accumulate the blended source atom into
	 * @param SourceAtom The target transformation (used when BlendWeight = 1)
	 * @param Alpha The blend weight between Identity and SourceAtom
	 */
	FORCEINLINE static void BlendFromIdentityAndAccumulate(Transform& FinalAtom, Transform& SourceAtom, const ScalarRegister& BlendWeight)
	{
		const VectorRegister Const0001 = GlobalVectorConstants::Float0001;
		const VectorRegister ConstNegative0001 = VectorSubtract(VectorZero(), Const0001);
		const VectorRegister VOneMinusAlpha = VectorSubtract(VectorOne(), BlendWeight.Value);
		const VectorRegister DefaultScale = MakeVectorRegister(1.f, 1.f, 1.f, 0.f);

		// Blend rotation
		//     To ensure the 'shortest route', we make sure the dot product between the both rotations is positive.
		//     const float Bias = (|A.B| >= 0 ? 1 : -1)
		//     BlendedAtom.Rotation = (B * Alpha) + (A * (Bias * (1.f - Alpha)));
		//     BlendedAtom.Rotation.QuaternionNormalize();
		//  Note: A = (0,0,0,1), which simplifies things a lot; only care about sign of B.W now, instead of doing a dot product
		const VectorRegister RotationB = SourceAtom.Rotation;

		const VectorRegister QuatRotationDirMask = VectorCompareGE(RotationB, VectorZero());
		const VectorRegister BiasTimesA = VectorSelect(QuatRotationDirMask, Const0001, ConstNegative0001);
		const VectorRegister RotateBTimesWeight = VectorMultiply(RotationB, BlendWeight.Value);
		const VectorRegister UnnormalizedRotation = VectorMultiplyAdd(BiasTimesA, VOneMinusAlpha, RotateBTimesWeight);

		// Normalize blended rotation ( result = (Q.Q >= 1e-8) ? (Q / |Q|) : (0,0,0,1) )
		const VectorRegister BlendedRotation = VectorNormalizeSafe(UnnormalizedRotation, Const0001);

		// FinalAtom.Rotation = BlendedAtom.Rotation * FinalAtom.Rotation;
		FinalAtom.Rotation = VectorQuaternionMultiply2(BlendedRotation, FinalAtom.Rotation);

		// Blend translation and scale
		//    BlendedAtom.Translation = Lerp(Zero, SourceAtom.Translation, Alpha);
		//    BlendedAtom.Scale = Lerp(0, SourceAtom.Scale, Alpha);
		const VectorRegister BlendedTranslation	= Math::Lerp(VectorZero(), SourceAtom.Translation, BlendWeight.Value);
		const VectorRegister BlendedScale3D	= Math::Lerp(VectorZero(), SourceAtom.Scale3D, BlendWeight.Value);

		// Apply translation and scale to final atom
		//     FinalAtom.Translation += BlendedAtom.Translation
		//     FinalAtom.Scale *= BlendedAtom.Scale
		FinalAtom.Translation = VectorAdd( FinalAtom.Translation, BlendedTranslation );
		FinalAtom.Scale3D = VectorMultiply( FinalAtom.Scale3D, VectorAdd(DefaultScale, BlendedScale3D));
		checkSlow( FinalAtom.IsRotationNormalized() );
	}


	/**
	 * Returns the rotation component
	 *
	 * @return The rotation component
	 */
	FORCEINLINE Quat GetRotation() const
	{
		CheckNaN_Rotate();
		Quat OutRotation;
		VectorStoreAligned(Rotation, &OutRotation);
		return OutRotation;
	}

	/**
	 * Returns the translation component
	 *
	 * @return The translation component
	 */
	FORCEINLINE Vector3f GetTranslation() const
	{
		CheckNaN_Translate();
		Vector3f OutTranslation;
		VectorStoreFloat3(Translation, &OutTranslation);
		return OutTranslation;
	}

	/**
	 * Returns the Scale3D component
	 *
	 * @return The Scale3D component
	 */
	FORCEINLINE Vector3f GetScale3D() const
	{
		CheckNaN_Scale3D();
		Vector3f OutScale3D;
		VectorStoreFloat3(Scale3D, &OutScale3D);
		return OutScale3D;
	}

	/**
	 * Sets the Rotation and Scale3D of this transformation from another transform
	 *
	 * @param SrcBA The transform to copy rotation and Scale3D from
	 */
	FORCEINLINE void CopyRotationPart(const Transform& SrcBA)
	{
		Rotation = SrcBA.Rotation;
		Scale3D = SrcBA.Scale3D;

		CheckNaN_Rotate();
		CheckNaN_Scale3D();
	}

	/**
	 * Sets the Translation and Scale3D of this transformation from another transform
	 *
	 * @param SrcBA The transform to copy translation and Scale3D from
	 */
	FORCEINLINE void CopyTranslationAndScale3D(const Transform& SrcBA)
	{
		Translation = SrcBA.Translation;
		Scale3D = SrcBA.Scale3D;

		CheckNaN_Translate();
		CheckNaN_Scale3D();
	}

	void SetFromMatrix(const Matrix& InMatrix)
	{
		Matrix M = InMatrix;

		// Get the 3D scale from the matrix
		Vector3f InScale = M.ExtractScaling();
		Scale3D = VectorLoadFloat3_W0(&InScale);

		// If there is negative scaling going on, we handle that here
		if(InMatrix.Determinant() < 0.f)
		{
			// Assume it is along X and modify transform accordingly. 
			// It doesn't actually matter which axis we choose, the 'appearance' will be the same			
			Scale3D = VectorMultiply(Scale3D, GlobalVectorConstants::FloatMinus1_111 );			
			M.SetAxis(0, -M.GetScaledAxis( EAxis::X ));
		}

		Quat InRotation = Quat(M);
		Rotation = VectorLoadAligned(&InRotation);
		Vector3f InTranslation = InMatrix.GetOrigin();
		Translation = VectorLoadFloat3_W0(&InTranslation);

		// Normalize rotation
		Rotation = VectorNormalizeQuaternion(Rotation);		
	}

private:
	FORCEINLINE void ToMatrixInternal( VectorRegister& OutDiagonals, VectorRegister& OutAdds, VectorRegister& OutSubtracts ) const
	{
#if !(BUILD_SHIPPING)
		// Make sure Rotation is normalized when we turn it into a matrix.
		check( IsRotationNormalized() );
#endif		

		const VectorRegister RotationX2Y2Z2 = VectorAdd(Rotation, Rotation);	// x2, y2, z2
		const VectorRegister RotationXX2YY2ZZ2 = VectorMultiply(RotationX2Y2Z2, Rotation);	// xx2, yy2, zz2		

		// The diagonal terms of the rotation matrix are:
		//   (1 - (yy2 + zz2)) * scale
		//   (1 - (xx2 + zz2)) * scale
		//   (1 - (xx2 + yy2)) * scale
		const VectorRegister yy2_xx2_xx2 = VectorSwizzle(RotationXX2YY2ZZ2, 1, 0, 0, 0);
		const VectorRegister zz2_zz2_yy2 = VectorSwizzle(RotationXX2YY2ZZ2, 2, 2, 1, 0);
		const VectorRegister DiagonalSum = VectorAdd(yy2_xx2_xx2, zz2_zz2_yy2);
		const VectorRegister Diagonals = VectorSubtract(VectorOne(), DiagonalSum);
		OutDiagonals = VectorMultiply(Diagonals, Scale3D);

		// Grouping the non-diagonal elements in the rotation block by operations:
		//    ((x*y2,y*z2,x*z2) + (w*z2,w*x2,w*y2)) * scale.xyz and
		//    ((x*y2,y*z2,x*z2) - (w*z2,w*x2,w*y2)) * scale.yxz
		// Rearranging so the LHS and RHS are in the same order as for +
		//    ((x*y2,y*z2,x*z2) - (w*z2,w*x2,w*y2)) * scale.yxz

		// RotBase = x*y2, y*z2, x*z2
		// RotOffset = w*z2, w*x2, w*y2
		const VectorRegister x_y_x = VectorSwizzle(Rotation, 0, 1, 0, 0);
		const VectorRegister y2_z2_z2 = VectorSwizzle(RotationX2Y2Z2, 1, 2, 2, 0);
		const VectorRegister RotBase = VectorMultiply(x_y_x, y2_z2_z2);

		const VectorRegister w_w_w = VectorReplicate(Rotation, 3);
		const VectorRegister z2_x2_y2 = VectorSwizzle(RotationX2Y2Z2, 2, 0, 1, 0);
		const VectorRegister RotOffset = VectorMultiply(w_w_w, z2_x2_y2);

		// Adds = (RotBase + RotOffset)*Scale3D :  (x*y2 + w*z2) * Scale3D.X , (y*z2 + w*x2) * Scale3D.Y, (x*z2 + w*y2) * Scale3D.Z
		// Subtracts = (RotBase - RotOffset)*Scale3DYZX :  (x*y2 - w*z2) * Scale3D.Y , (y*z2 - w*x2) * Scale3D.Z, (x*z2 - w*y2) * Scale3D.X
		const VectorRegister Adds = VectorAdd(RotBase, RotOffset);
		OutAdds = VectorMultiply(Adds, Scale3D);
		const VectorRegister Scale3DYZXW = VectorSwizzle( Scale3D, 1, 2, 0, 3);
		const VectorRegister Subtracts = VectorSubtract(RotBase, RotOffset);
		OutSubtracts = VectorMultiply(Subtracts , Scale3DYZXW);
	}

	FORCEINLINE void ToMatrixInternalNoScale( VectorRegister& OutDiagonals, VectorRegister& OutAdds, VectorRegister& OutSubtracts ) const
	{
#if !(UE_BUILD_SHIPPING || UE_BUILD_TEST) && WITH_EDITORONLY_DATA
		// Make sure Rotation is normalized when we turn it into a matrix.
		ensure( IsRotationNormalized() );
#endif		
		const VectorRegister RotationX2Y2Z2 = VectorAdd(Rotation, Rotation);	// x2, y2, z2
		const VectorRegister RotationXX2YY2ZZ2 = VectorMultiply(RotationX2Y2Z2, Rotation);	// xx2, yy2, zz2		

		// The diagonal terms of the rotation matrix are:
		//   (1 - (yy2 + zz2))
		//   (1 - (xx2 + zz2))
		//   (1 - (xx2 + yy2))
		const VectorRegister yy2_xx2_xx2 = VectorSwizzle(RotationXX2YY2ZZ2, 1, 0, 0, 0);
		const VectorRegister zz2_zz2_yy2 = VectorSwizzle(RotationXX2YY2ZZ2, 2, 2, 1, 0);
		const VectorRegister DiagonalSum = VectorAdd(yy2_xx2_xx2, zz2_zz2_yy2);
		OutDiagonals = VectorSubtract(VectorOne(), DiagonalSum);

		// Grouping the non-diagonal elements in the rotation block by operations:
		//    ((x*y2,y*z2,x*z2) + (w*z2,w*x2,w*y2)) and
		//    ((x*y2,y*z2,x*z2) - (w*z2,w*x2,w*y2))
		// Rearranging so the LHS and RHS are in the same order as for +
		//    ((x*y2,y*z2,x*z2) - (w*z2,w*x2,w*y2))

		// RotBase = x*y2, y*z2, x*z2
		// RotOffset = w*z2, w*x2, w*y2
		const VectorRegister x_y_x = VectorSwizzle(Rotation, 0, 1, 0, 0);
		const VectorRegister y2_z2_z2 = VectorSwizzle(RotationX2Y2Z2, 1, 2, 2, 0);
		const VectorRegister RotBase = VectorMultiply(x_y_x, y2_z2_z2);

		const VectorRegister w_w_w = VectorReplicate(Rotation, 3);
		const VectorRegister z2_x2_y2 = VectorSwizzle(RotationX2Y2Z2, 2, 0, 1, 0);
		const VectorRegister RotOffset = VectorMultiply(w_w_w, z2_x2_y2);

		// Adds = (RotBase + RotOffset):  (x*y2 + w*z2) , (y*z2 + w*x2), (x*z2 + w*y2)
		// Subtracts = (RotBase - RotOffset) :  (x*y2 - w*z2) , (y*z2 - w*x2), (x*z2 - w*y2)
		OutAdds = VectorAdd(RotBase, RotOffset);		
		OutSubtracts = VectorSubtract(RotBase, RotOffset);
	}

	/** 
	 * mathematically if you have 0 scale, it should be infinite, 
	 * however, in practice if you have 0 scale, and relative transform doesn't make much sense 
	 * anymore because you should be instead of showing gigantic infinite mesh
	 * also returning BIG_NUMBER causes sequential NaN issues by multiplying 
	 * so we hardcode as 0
	 */
	static FORCEINLINE VectorRegister		GetSafeScaleReciprocal(const VectorRegister& InScale, const ScalarRegister& Tolerance = ScalarRegister(GlobalVectorConstants::SmallNumber))
	{		
		// SafeReciprocalScale.X = (InScale.X == 0) ? 0.f : 1/InScale.X; // same for YZW
		VectorRegister SafeReciprocalScale;

		/// VectorRegister( 1.0f / InScale.x, 1.0f / InScale.y, 1.0f / InScale.z, 1.0f / InScale.w )
		const VectorRegister ReciprocalScale = VectorReciprocalAccurate(InScale);
		
		//VectorRegister( Vec1.x == Vec2.x ? 0xFFFFFFFF : 0, same for yzw )
		const VectorRegister ScaleZeroMask = VectorCompareGE(Tolerance.Value, VectorAbs(InScale));

		//const VectorRegister ScaleZeroMask = VectorCompareEQ(InScale, VectorZero());

		// VectorRegister( for each bit i: Mask[i] ? Vec1[i] : Vec2[i] )
		SafeReciprocalScale = VectorSelect(ScaleZeroMask, VectorZero(), ReciprocalScale);

		return SafeReciprocalScale;
	}

	/** Returns Inverse Transform of this Transform **/
	FORCEINLINE Transform InverseFast() const
	{
		// Inverse QST (A) = QST (~A)
		// Since A*~A = Identity, 
		// A(P) = Q(A)*S(A)*P*-Q(A) + T(A)
		// ~A(A(P)) = Q(~A)*S(~A)*(Q(A)*S(A)*P*-Q(A) + T(A))*-Q(~A) + T(~A) = Identity
		// Q(~A)*Q(A)*S(~A)*S(A)*P*-Q(A)*-Q(~A) + Q(~A)*S(~A)*T(A)*-Q(~A) + T(~A) = Identity
		// [Q(~A)*Q(A)]*[S(~A)*S(A)]*P*-[Q(~A)*Q(A)] + [Q(~A)*S(~A)*T(A)*-Q(~A) + T(~A)] = I

		// Identity Q = (0, 0, 0, 1) = Q(~A)*Q(A)
		// Identity Scale = 1 = S(~A)*S(A)
		// Identity Translation = (0, 0, 0) = [Q(~A)*S(~A)*T(A)*-Q(~A) + T(~A)]

		//	Q(~A) = Q(~A)
		//	S(~A) = 1.f/S(A)
		//	T(~A) = - (Q(~A)*S(~A)*T(A)*Q(A))	
		checkSlow(IsRotationNormalized());
		checkSlow(VectorAnyGreaterThan(VectorAbs(Scale3D), GlobalVectorConstants::SmallNumber));

		// Invert the scale
		const VectorRegister InvScale = VectorSet_W0(GetSafeScaleReciprocal(VectorSet_W1(Scale3D), ScalarRegister(GlobalVectorConstants::SmallNumber)));

		// Invert the rotation
		const VectorRegister InvRotation = VectorQuaternionInverse(Rotation);

		// Invert the translation
		const VectorRegister ScaledTranslation = VectorMultiply(InvScale, Translation);
		const VectorRegister t2 = VectorQuaternionRotateVector(InvRotation, ScaledTranslation);
		const VectorRegister InvTranslation = VectorSet_W0(VectorNegate(t2));

		return Transform(InvRotation, InvTranslation, InvScale);
	}

	/**
	* Create a new transform: OutTransform = A * B using the matrix while keeping the scale that's given by A and B
	* Please note that this operation is a lot more expensive than normal Multiply
	*
	* Order matters when composing transforms : A * B will yield a transform that logically first applies A then B to any subsequent transformation.
	*
	* @param  OutTransform pointer to transform that will store the result of A * B.
	* @param  A Transform A.
	* @param  B Transform B.
	*/
	FORCEINLINE static void MultiplyUsingMatrixWithScale(Transform* OutTransform, const Transform* A, const Transform* B);
	/**
	* Create a new transform from multiplications of given to matrices (AMatrix*BMatrix) using desired scale
	* This is used by MultiplyUsingMatrixWithScale and GetRelativeTransformUsingMatrixWithScale
	* This is only used to handle negative scale
	*
	* @param	AMatrix first Matrix of operation
	* @param	BMatrix second Matrix of operation
	* @param	DesiredScale - there is no check on if the magnitude is correct here. It assumes that is correct. 
	* @param	OutTransform the constructed transform 
	*/
	FORCEINLINE static void ConstructTransformFromMatrixWithDesiredScale(const Matrix& AMatrix, const Matrix& BMatrix, const VectorRegister& DesiredScale, Transform& OutTransform);
	/**
	* Create a new transform: OutTransform = Base * Relative(-1) using the matrix while keeping the scale that's given by Base and Relative
	* Please note that this operation is a lot more expensive than normal GetRelativeTrnasform
	*
	* @param  OutTransform pointer to transform that will store the result of Base * Relative(-1).
	* @param  BAse Transform Base.
	* @param  Relative Transform Relative.
	*/
	static void GetRelativeTransformUsingMatrixWithScale(Transform* OutTransform, const Transform* Base, const Transform* Relative);
};

FORCEINLINE bool Transform::AnyHasNegativeScale(const Vector3f& InScale3D, const Vector3f& InOtherScale3D)
{
	VectorRegister VectorInScale3D = VectorLoadFloat3_W0(&InScale3D);
	VectorRegister VectorInOtherScale3D = VectorLoadFloat3_W0(&InOtherScale3D);

	return Private_AnyHasNegativeScale(VectorInScale3D, VectorInOtherScale3D);
}

/** Scale the translation part of the Transform by the supplied vector. */
FORCEINLINE void Transform::ScaleTranslation(const Vector3f& InScale3D)
{
	VectorRegister VectorInScale3D = VectorLoadFloat3_W0(&InScale3D);
	Translation = VectorMultiply( Translation, VectorInScale3D );
	CheckNaN_Translate();
}

/** Scale the translation part of the Transform by the supplied value. */
FORCEINLINE void Transform::ScaleTranslation(const float& InScale)
{
	ScaleTranslation( Vector3f(InScale) );
}

// this function is from matrix, and all it does is to normalize rotation portion
FORCEINLINE void Transform::RemoveScaling(float Tolerance/*=SMALL_NUMBER*/)
{
	Scale3D = VectorSet_W0( VectorOne() );
	NormalizeRotation();	

	CheckNaN_Rotate();
	CheckNaN_Scale3D();
}

FORCEINLINE void Transform::MultiplyUsingMatrixWithScale(Transform* OutTransform, const Transform* A, const Transform* B)
{
	ConstructTransformFromMatrixWithDesiredScale(A->ToMatrixWithScale(), B->ToMatrixWithScale(), VectorMultiply(A->Scale3D, B->Scale3D), *OutTransform);
}

FORCEINLINE void Transform::ConstructTransformFromMatrixWithDesiredScale(const Matrix& AMatrix, const Matrix& BMatrix, const VectorRegister& DesiredScale, Transform& OutTransform)
{
	// the goal of using M is to get the correct orientation
	// but for translation, we still need scale
	Matrix M = AMatrix * BMatrix;
	M.RemoveScaling();

	// apply negative scale back to axes
	Vector3f SignedScale;
	VectorStoreFloat3(VectorSign(DesiredScale), &SignedScale);

	M.SetAxis(0, SignedScale.X * M.GetScaledAxis(EAxis::X));
	M.SetAxis(1, SignedScale.Y * M.GetScaledAxis(EAxis::Y));
	M.SetAxis(2, SignedScale.Z * M.GetScaledAxis(EAxis::Z));

	// @note: if you have negative with 0 scale, this will return rotation that is identity
	// since matrix loses that axes
	Quat Rotation = Quat(M);
	Rotation.Normalize();

	// set values back to output
	OutTransform.Scale3D = DesiredScale;
	OutTransform.Rotation = VectorLoadAligned(&Rotation);

	// technically I could calculate this using Transform but then it does more quat multiplication 
	// instead of using Scale in matrix multiplication
	// it's a question of between RemoveScaling vs using Transform to move translation
	Vector3f Translation = M.GetOrigin();
	OutTransform.Translation = VectorLoadFloat3_W0(&Translation);
}

/** Returns Multiplied Transform of 2 Transforms **/
FORCEINLINE void Transform::Multiply(Transform* OutTransform, const Transform* A, const Transform* B)
{
	A->CheckNaN_All();
	B->CheckNaN_All();

	checkSlow(A->IsRotationNormalized());
	checkSlow(B->IsRotationNormalized());

	//	When Q = quaternion, S = single scalar scale, and T = translation
	//	QST(A) = Q(A), S(A), T(A), and QST(B) = Q(B), S(B), T(B)

	//	QST (AxB) 

	// QST(A) = Q(A)*S(A)*P*-Q(A) + T(A)
	// QST(AxB) = Q(B)*S(B)*QST(A)*-Q(B) + T(B)
	// QST(AxB) = Q(B)*S(B)*[Q(A)*S(A)*P*-Q(A) + T(A)]*-Q(B) + T(B)
	// QST(AxB) = Q(B)*S(B)*Q(A)*S(A)*P*-Q(A)*-Q(B) + Q(B)*S(B)*T(A)*-Q(B) + T(B)
	// QST(AxB) = [Q(B)*Q(A)]*[S(B)*S(A)]*P*-[Q(B)*Q(A)] + Q(B)*S(B)*T(A)*-Q(B) + T(B)

	//	Q(AxB) = Q(B)*Q(A)
	//	S(AxB) = S(A)*S(B)
	//	T(AxB) = Q(B)*S(B)*T(A)*-Q(B) + T(B)
	checkSlow(VectorGetComponent(A->Scale3D, 3) == 0.f);
	checkSlow(VectorGetComponent(B->Scale3D, 3) == 0.f);

	if (Private_AnyHasNegativeScale(A->Scale3D, B->Scale3D))
	{
		// @note, if you have 0 scale with negative, you're going to lose rotation as it can't convert back to quat
		MultiplyUsingMatrixWithScale(OutTransform, A, B);
	}
	else
	{
		const VectorRegister QuatA = A->Rotation;
		const VectorRegister QuatB = B->Rotation;
		const VectorRegister TranslateA = A->Translation;
		const VectorRegister TranslateB = B->Translation;
		const VectorRegister ScaleA = A->Scale3D;
		const VectorRegister ScaleB = B->Scale3D;

		// RotationResult = B.Rotation * A.Rotation
		OutTransform->Rotation = VectorQuaternionMultiply2(QuatB, QuatA);

		// TranslateResult = B.Rotate(B.Scale * A.Translation) + B.Translate
		const VectorRegister ScaledTransA = VectorMultiply(TranslateA, ScaleB);
		const VectorRegister RotatedTranslate = VectorQuaternionRotateVector(QuatB, ScaledTransA);
		OutTransform->Translation = VectorAdd(RotatedTranslate, TranslateB);

		// ScaleResult = Scale.B * Scale.A
		OutTransform->Scale3D = VectorMultiply(ScaleA, ScaleB);;
	}
}
/** 
 * Apply Scale to this transform
 */
FORCEINLINE Transform Transform::GetScaled(float InScale) const
{
	Transform A(*this);
	
	VectorRegister VScale = VectorLoadFloat1(&InScale);
	A.Scale3D = VectorMultiply( A.Scale3D, VScale);

	A.CheckNaN_Scale3D();

	return A;
}

/** 
 * Apply Scale to this transform
 */
FORCEINLINE Transform Transform::GetScaled(Vector3f InScale) const
{
	Transform A(*this);

	VectorRegister VScale = VectorLoadFloat3_W0(&InScale);
	A.Scale3D = VectorMultiply( A.Scale3D, VScale);

	A.CheckNaN_Scale3D();

	return A;
}

FORCEINLINE Vector4f Transform::TransformVector4NoScale(const Vector4f& V) const
{
	CheckNaN_All();

	// if not, this won't work
	checkSlow (V.W == 0.f || V.W == 1.f);

	const VectorRegister InputVector = VectorLoadAligned(&V);

	//Transform using QST is following
	//QST(P) = Q.Rotate(S*P) + T where Q = quaternion, S = 1.0f, T = translation

	//RotatedVec = Q.Rotate(V.X, V.Y, V.Z, 0.f)
	const VectorRegister InputVectorW0 = VectorSet_W0(InputVector);	
	const VectorRegister RotatedVec = VectorQuaternionRotateVector(Rotation, InputVectorW0);

	// NewVect.XYZ += Translation * W
	// NewVect.W += 1 * W
	const VectorRegister WWWW = VectorReplicate(InputVector, 3);
	const VectorRegister TranslatedVec = VectorMultiplyAdd(Translation, WWWW, RotatedVec);

	Vector4f NewVectOutput;
	VectorStoreAligned(TranslatedVec, &NewVectOutput);
	return NewVectOutput;
}

FORCEINLINE Vector4f Transform::TransformVector4(const Vector4f& V) const
{
	CheckNaN_All();

	// if not, this won't work
	checkSlow (V.W == 0.f || V.W == 1.f);

	const VectorRegister InputVector = VectorLoadAligned(&V);

	//Transform using QST is following
	//QST(P) = Q.Rotate(S*P) + T where Q = quaternion, S = scale, T = translation

	//RotatedVec = Q.Rotate(Scale*V.X, Scale*V.Y, Scale*V.Z, 0.f)
	const VectorRegister InputVectorW0 = VectorSet_W0(InputVector);
	const VectorRegister ScaledVec = VectorMultiply(Scale3D, InputVectorW0);
	const VectorRegister RotatedVec = VectorQuaternionRotateVector(Rotation, ScaledVec);

	// NewVect.XYZ += Translation * W
	// NewVect.W += 1 * W
	const VectorRegister WWWW = VectorReplicate(InputVector, 3);
	const VectorRegister TranslatedVec = VectorMultiplyAdd(Translation, WWWW, RotatedVec);

	Vector4f NewVectOutput;
	VectorStoreAligned(TranslatedVec, &NewVectOutput);
	return NewVectOutput;
}


FORCEINLINE Vector3f Transform::TransformPosition(const Vector3f& V) const
{
	CheckNaN_All();

	const VectorRegister InputVectorW0 = VectorLoadFloat3_W0(&V);

	//Transform using QST is following
	//QST(P) = Q.Rotate(S*P) + T where Q = quaternion, S = scale, T = translation
	
	//RotatedVec = Q.Rotate(Scale*V.X, Scale*V.Y, Scale*V.Z, 0.f)
	const VectorRegister ScaledVec = VectorMultiply(Scale3D, InputVectorW0);
	const VectorRegister RotatedVec = VectorQuaternionRotateVector(Rotation, ScaledVec);

	const VectorRegister TranslatedVec = VectorAdd(RotatedVec, Translation);

	Vector3f Result;
	VectorStoreFloat3(TranslatedVec, &Result);
	return Result;
}

FORCEINLINE Vector3f Transform::TransformPositionNoScale(const Vector3f& V) const
{
	CheckNaN_All();

	const VectorRegister InputVectorW0 = VectorLoadFloat3_W0(&V);

	//Transform using QST is following
	//QST(P) = Q.Rotate(S*P) + T where Q = quaternion, S = 1.0f, T = translation

	//RotatedVec = Q.Rotate(V.X, V.Y, V.Z, 0.f)
	const VectorRegister RotatedVec = VectorQuaternionRotateVector(Rotation, InputVectorW0);

	const VectorRegister TranslatedVec = VectorAdd(RotatedVec, Translation);

	Vector3f Result;
	VectorStoreFloat3(TranslatedVec, &Result);
	return Result;
}

FORCEINLINE Vector3f Transform::TransformVector(const Vector3f& V) const
{
	CheckNaN_All();

	const VectorRegister InputVectorW0 = VectorLoadFloat3_W0(&V);

	//RotatedVec = Q.Rotate(Scale*V.X, Scale*V.Y, Scale*V.Z, 0.f)
	const VectorRegister ScaledVec = VectorMultiply(Scale3D, InputVectorW0);
	const VectorRegister RotatedVec = VectorQuaternionRotateVector(Rotation, ScaledVec);

	Vector3f Result;
	VectorStoreFloat3(RotatedVec, &Result);
	return Result;
}

FORCEINLINE Vector3f Transform::TransformVectorNoScale(const Vector3f& V) const
{
	CheckNaN_All();

	const VectorRegister InputVectorW0 = VectorLoadFloat3_W0(&V);

	//RotatedVec = Q.Rotate(V.X, V.Y, V.Z, 0.f)
	const VectorRegister RotatedVec = VectorQuaternionRotateVector(Rotation, InputVectorW0);

	Vector3f Result;
	VectorStoreFloat3(RotatedVec, &Result);
	return Result;
}

// do backward operation when inverse, translation -> rotation -> scale
FORCEINLINE Vector3f Transform::InverseTransformPosition(const Vector3f &V) const
{
	CheckNaN_All();

	const VectorRegister InputVector = VectorLoadFloat3_W0(&V);

	// (V-Translation)
	const VectorRegister TranslatedVec = VectorSet_W0(VectorSubtract(InputVector, Translation));

	// ( Rotation.Inverse() * (V-Translation) )
	const VectorRegister VR = VectorQuaternionInverseRotateVector(Rotation, TranslatedVec);

	// GetSafeScaleReciprocal(Scale3D);
	const VectorRegister SafeReciprocal = GetSafeScaleReciprocal(Scale3D);	

	// ( Rotation.Inverse() * (V-Translation) ) * GetSafeScaleReciprocal(Scale3D);
	const VectorRegister VResult = VectorMultiply(VR, SafeReciprocal);

	Vector3f Result;
	VectorStoreFloat3(VResult, &Result);
	return Result;
}

// do backward operation when inverse, translation -> rotation
FORCEINLINE Vector3f Transform::InverseTransformPositionNoScale(const Vector3f &V) const
{
	CheckNaN_All();

	const VectorRegister InputVector = VectorLoadFloat3_W0(&V);

	// (V-Translation)
	const VectorRegister TranslatedVec = VectorSet_W0(VectorSubtract(InputVector, Translation));

	// ( Rotation.Inverse() * (V-Translation) )
	const VectorRegister VResult = VectorQuaternionInverseRotateVector(Rotation, TranslatedVec);

	Vector3f Result;
	VectorStoreFloat3(VResult, &Result);
	return Result;
}


// do backward operation when inverse, translation -> rotation -> scale
FORCEINLINE Vector3f Transform::InverseTransformVector(const Vector3f &V) const
{
	CheckNaN_All();

	const VectorRegister InputVector = VectorLoadFloat3_W0(&V);

	// ( Rotation.Inverse() * V ) aka. Vector3f Quat::operator*( const Vector3f& V ) const
	const VectorRegister VR = VectorQuaternionInverseRotateVector(Rotation, InputVector);

	// GetSafeScaleReciprocal(Scale3D);
	const VectorRegister SafeReciprocal = GetSafeScaleReciprocal(Scale3D);

	// ( Rotation.Inverse() * V) * GetSafeScaleReciprocal(Scale3D);
	const VectorRegister VResult = VectorMultiply(VR, SafeReciprocal);

	Vector3f Result;
	VectorStoreFloat3(VResult, &Result);
	return Result;
}

// do backward operation when inverse, translation -> rotation
FORCEINLINE Vector3f Transform::InverseTransformVectorNoScale(const Vector3f &V) const
{
	CheckNaN_All();

	VectorRegister InputVector = VectorLoadFloat3_W0(&V);

	// ( Rotation.Inverse() * V )
	VectorRegister VResult = VectorQuaternionInverseRotateVector(Rotation, InputVector);

	Vector3f Result;
	VectorStoreFloat3(VResult, &Result);
	return Result;
}

FORCEINLINE Quat Transform::TransformRotation(const Quat& Q) const
{
	return GetRotation() * Q;
}

FORCEINLINE Quat Transform::InverseTransformRotation(const Quat& Q) const
{
	return GetRotation().Inverse() * Q;
}

FORCEINLINE Transform Transform::operator*(const Transform& Other) const
{
	Transform Output;
	Multiply(&Output, this, &Other);
	return Output;
}

FORCEINLINE void Transform::operator*=(const Transform& Other)
{
	Multiply(this, this, &Other);
}

FORCEINLINE Transform Transform::operator*(const Quat& Other) const
{
	Transform Output, OtherTransform(Other, Vector3f::ZeroVector, Vector3f::OneVector);
	Multiply(&Output, this, &OtherTransform);
	return Output;
}

FORCEINLINE void Transform::operator*=(const Quat& Other)
{
	Transform OtherTransform(Other, Vector3f::ZeroVector, Vector3f::OneVector);
	Multiply(this, this, &OtherTransform);
}

// x = 0, y = 1, z = 2
FORCEINLINE Vector3f Transform::GetScaledAxis( EAxis::Type InAxis ) const
{
	if ( InAxis == EAxis::X )
	{
		return TransformVector(Vector3f(1.f, 0.f, 0.f));
	}
	else if ( InAxis == EAxis::Y )
	{
		return TransformVector(Vector3f(0.f, 1.f, 0.f));
	}

	return TransformVector(Vector3f(0.f, 0.f, 1.f));
}

// x = 0, y = 1, z = 2
FORCEINLINE Vector3f Transform::GetUnitAxis( EAxis::Type InAxis ) const
{
	if ( InAxis == EAxis::X )
	{
		return TransformVectorNoScale(Vector3f(1.f, 0.f, 0.f));
	}
	else if ( InAxis == EAxis::Y )
	{
		return TransformVectorNoScale(Vector3f(0.f, 1.f, 0.f));
	}

	return TransformVectorNoScale(Vector3f(0.f, 0.f, 1.f));
}

FORCEINLINE void Transform::Mirror(EAxis::Type MirrorAxis, EAxis::Type FlipAxis)
{
	// We do convert to Matrix for mirroring. 
	Matrix M = ToMatrixWithScale();
	M.Mirror(MirrorAxis, FlipAxis);
	SetFromMatrix(M);
}

/** same version of Matrix::GetMaximumAxisScale function **/
/** @return the maximum magnitude of any row of the matrix. */
FORCEINLINE float Transform::GetMaximumAxisScale() const
{
	CheckNaN_Scale3D();

	float Scale3DAbsMax;
	// Scale3DAbsXYZ1 = { Abs(X), Abs(Y)), Abs(Z), 0 }
	const VectorRegister Scale3DAbsXYZ0 =  VectorAbs(Scale3D);
	// Scale3DAbsYZX1 = { Abs(Y),Abs(Z)),Abs(X), 0 }
	const VectorRegister Scale3DAbsYZX0 = VectorSwizzle(Scale3DAbsXYZ0, 1,2,0,3);
	// Scale3DAbsZXY1 = { Abs(Z),Abs(X)),Abs(Y), 0 }
	const VectorRegister Scale3DAbsZXY0 = VectorSwizzle(Scale3DAbsXYZ0, 2,0,1,3);
	// t0 = { Max(Abs(X), Abs(Y)),  Max(Abs(Y), Abs(Z)), Max(Abs(Z), Abs(X)), 0 }
	const VectorRegister t0 = VectorMax(Scale3DAbsXYZ0, Scale3DAbsYZX0);
	// t1 = { Max(Abs(X), Abs(Y), Abs(Z)), Max(Abs(Y), Abs(Z), Abs(X)), Max(Abs(Z), Abs(X), Abs(Y)), 0 }
	const VectorRegister t2 = VectorMax(t0, Scale3DAbsZXY0);
	// Scale3DAbsMax = Max(Abs(X), Abs(Y), Abs(Z));
	VectorStoreFloat1(t2, &Scale3DAbsMax);

	return Scale3DAbsMax;
}

/** @return the minimum magnitude of all components of the 3D scale. */
FORCEINLINE float Transform::GetMinimumAxisScale() const
{
	CheckNaN_Scale3D();

	float Scale3DAbsMin;
	// Scale3DAbsXYZ1 = { Abs(X), Abs(Y)), Abs(Z), 0 }
	const VectorRegister Scale3DAbsXYZ0 =  VectorAbs(Scale3D);
	// Scale3DAbsYZX1 = { Abs(Y),Abs(Z)),Abs(X), 0 }
	const VectorRegister Scale3DAbsYZX0 = VectorSwizzle(Scale3DAbsXYZ0, 1,2,0,3);
	// Scale3DAbsZXY1 = { Abs(Z),Abs(X)),Abs(Y), 0 }
	const VectorRegister Scale3DAbsZXY0 = VectorSwizzle(Scale3DAbsXYZ0, 2,0,1,3);
	// t0 = { Min(Abs(X), Abs(Y)),  Min(Abs(Y), Abs(Z)), Min(Abs(Z), Abs(X)), 0 }
	const VectorRegister t0 = VectorMin(Scale3DAbsXYZ0, Scale3DAbsYZX0);
	// t1 = { Min(Abs(X), Abs(Y), Abs(Z)), Min(Abs(Y), Abs(Z), Abs(X)), Min(Abs(Z), Abs(X), Abs(Y)), 0 }
	const VectorRegister t2 = VectorMin(t0, Scale3DAbsZXY0);
	// Scale3DAbsMax = Min(Abs(X), Abs(Y), Abs(Z));
	VectorStoreFloat1(t2, &Scale3DAbsMin);

	return Scale3DAbsMin;
}

/** 
 * mathematically if you have 0 scale, it should be infinite, 
 * however, in practice if you have 0 scale, and relative transform doesn't make much sense 
 * anymore because you should be instead of showing gigantic infinite mesh
 * also returning BIG_NUMBER causes sequential NaN issues by multiplying 
 * so we hardcode as 0
 */
FORCEINLINE Vector3f Transform::GetSafeScaleReciprocal(const Vector3f& InScale, float Tolerance)
{
	Vector3f SafeReciprocalScale;
	if (Math::Abs(InScale.X) <= Tolerance)
	{
		SafeReciprocalScale.X = 0.f;
	}
	else
	{
		SafeReciprocalScale.X = 1/InScale.X;
	}

	if (Math::Abs(InScale.Y) <= Tolerance)
	{
		SafeReciprocalScale.Y = 0.f;
	}
	else
	{
		SafeReciprocalScale.Y = 1/InScale.Y;
	}

	if (Math::Abs(InScale.Z) <= Tolerance)
	{
		SafeReciprocalScale.Z = 0.f;
	}
	else
	{
		SafeReciprocalScale.Z = 1/InScale.Z;
	}

	return SafeReciprocalScale;
}

#endif
};
