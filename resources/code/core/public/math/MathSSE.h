// this tiny game engine is rendered by CPU. I make it for learing the pipline and rasterization. This is the first step to realise my own game engine. so just do it!

#pragma once

#include "platform/Platform.h"
#include <math.h>

// We require SSE2
#include <emmintrin.h>

namespace LXF {

typedef __m128	VectorRegister;
typedef __m128i VectorRegisterInt;
typedef __m128d VectorRegisterDouble;


/**
 * @param A0	Selects which element (0-3) from 'A' into 1st slot in the result
 * @param A1	Selects which element (0-3) from 'A' into 2nd slot in the result
 * @param B2	Selects which element (0-3) from 'B' into 3rd slot in the result
 * @param B3	Selects which element (0-3) from 'B' into 4th slot in the result
 */
#define SHUFFLEMASK(A0,A1,B2,B3) ( (A0) | ((A1)<<2) | ((B2)<<4) | ((B3)<<6) )


FORCEINLINE VectorRegister MakeVectorRegister(uint32 X, uint32 Y, uint32 Z, uint32 W)
{
	union { VectorRegister v; VectorRegisterInt i; } Tmp;
	Tmp.i = _mm_setr_epi32(X, Y, Z, W);
	return Tmp.v;
}

FORCEINLINE VectorRegister MakeVectorRegister(float X, float Y, float Z, float W)
{
	return _mm_setr_ps(X, Y, Z, W);
}

FORCEINLINE VectorRegisterInt MakeVectorRegisterInt(int32 X, int32 Y, int32 Z, int32 W)
{
	return _mm_setr_epi32(X, Y, Z, W);
}

FORCEINLINE VectorRegister VectorLoad(const void* Ptr)
{
	return _mm_loadu_ps((float*)(Ptr));
}

#include "VectorConstants.h"

FORCEINLINE VectorRegister VectorZero(void)
{
	return _mm_setzero_ps();
}

FORCEINLINE VectorRegister VectorOne(void)
{
	return (GlobalVectorConstants::FloatOne);
}

FORCEINLINE float VectorGetComponent(VectorRegister Vec, uint32 ComponentIndex)
{
	return (((float*)&(Vec))[ComponentIndex]);
}

#define VectorLoadFloat3( Ptr )			MakeVectorRegister( ((const float*)(Ptr))[0], ((const float*)(Ptr))[1], ((const float*)(Ptr))[2], 0.0f )
#define VectorLoadFloat3_W0( Ptr )		MakeVectorRegister( ((const float*)(Ptr))[0], ((const float*)(Ptr))[1], ((const float*)(Ptr))[2], 0.0f )
#define VectorLoadFloat3_W1( Ptr )		MakeVectorRegister( ((const float*)(Ptr))[0], ((const float*)(Ptr))[1], ((const float*)(Ptr))[2], 1.0f )
#define VectorLoadAligned( Ptr )		_mm_load_ps( (const float*)(Ptr) )
#define VectorLoadFloat1( Ptr )			_mm_load1_ps( (const float*)(Ptr) )
#define VectorLoadFloat2( Ptr )			_mm_castpd_ps(_mm_load1_pd((const double*)(Ptr)))
#define VectorSetFloat3( X, Y, Z )		MakeVectorRegister( X, Y, Z, 0.0f )
#define VectorSetFloat1( F )	_mm_set1_ps( F )

#define VectorSet_W0( Vec )		_mm_and_ps( Vec, GlobalVectorConstants::XYZMask )
FORCEINLINE VectorRegister VectorSet_W1(const VectorRegister& Vector)
{
	// Temp = (Vector[2]. Vector[3], 1.0f, 1.0f)
	VectorRegister Temp = _mm_movehl_ps(VectorOne(), Vector);

	// Return (Vector[0], Vector[1], Vector[2], 1.0f)
	return _mm_shuffle_ps(Vector, Temp, SHUFFLEMASK(0, 1, 0, 3));
}

FORCEINLINE VectorRegister VectorSet(float X, float Y, float Z, float W)
{
	return MakeVectorRegister(X, Y, Z, W);
}
#define VectorStoreAligned( Vec, Ptr )	_mm_store_ps( (float*)(Ptr), Vec )
#define VectorStoreAlignedStreamed( Vec, Ptr )	_mm_stream_ps( (float*)(Ptr), Vec )
FORCEINLINE void VectorStore(const VectorRegister& Vec, void* Ptr)
{
	_mm_storeu_ps((float*)(Ptr), Vec);
}
FORCEINLINE void VectorStoreFloat3(const VectorRegister& Vec, void* Ptr)
{
	union { VectorRegister v; float f[4]; } Tmp;
	Tmp.v = Vec;
	float* FloatPtr = (float*)(Ptr);
	FloatPtr[0] = Tmp.f[0];
	FloatPtr[1] = Tmp.f[1];
	FloatPtr[2] = Tmp.f[2];
}
#define VectorStoreFloat1( Vec, Ptr )	_mm_store_ss((float*)(Ptr), Vec)
#define VectorStoreFloat1( Vec, Ptr )	_mm_store_ss((float*)(Ptr), Vec)
#define VectorReplicate( Vec, ElementIndex )	_mm_shuffle_ps( Vec, Vec, SHUFFLEMASK(ElementIndex,ElementIndex,ElementIndex,ElementIndex) )
#define VectorAbs( Vec )				_mm_and_ps(Vec, GlobalVectorConstants::SignMask)
#define VectorNegate( Vec )				_mm_sub_ps(_mm_setzero_ps(),Vec)

FORCEINLINE VectorRegister VectorAdd(const VectorRegister& Vec1, const VectorRegister& Vec2)
{
	return _mm_add_ps(Vec1, Vec2);
}

FORCEINLINE VectorRegister VectorSubtract(const VectorRegister& Vec1, const VectorRegister& Vec2)
{
	return _mm_sub_ps(Vec1, Vec2);
}

FORCEINLINE VectorRegister VectorMultiply(const VectorRegister& Vec1, const VectorRegister& Vec2)
{
	return _mm_mul_ps(Vec1, Vec2);
}

#define VectorMin( Vec1, Vec2 )			_mm_min_ps( Vec1, Vec2 )
#define VectorMax( Vec1, Vec2 )			_mm_max_ps( Vec1, Vec2 )
#define VectorSwizzle( Vec, X, Y, Z, W )	_mm_shuffle_ps( Vec, Vec, SHUFFLEMASK(X,Y,Z,W) )
#define VectorShuffle( Vec1, Vec2, X, Y, Z, W )	_mm_shuffle_ps( Vec1, Vec2, SHUFFLEMASK(X,Y,Z,W) )
#define VectorMultiplyAdd( Vec1, Vec2, Vec3 )	_mm_add_ps( _mm_mul_ps(Vec1, Vec2), Vec3 )

FORCEINLINE VectorRegister VectorDot3(const VectorRegister& Vec1, const VectorRegister& Vec2)
{
	VectorRegister Temp = VectorMultiply(Vec1, Vec2);
	return VectorAdd(VectorReplicate(Temp, 0), VectorAdd(VectorReplicate(Temp, 1), VectorReplicate(Temp, 2)));
}

FORCEINLINE VectorRegister VectorDot4(const VectorRegister& Vec1, const VectorRegister& Vec2)
{
	VectorRegister Temp1, Temp2;
	Temp1 = VectorMultiply(Vec1, Vec2);
	Temp2 = _mm_shuffle_ps(Temp1, Temp1, SHUFFLEMASK(2, 3, 0, 1));	// (Z,W,X,Y).
	Temp1 = VectorAdd(Temp1, Temp2);								// (X*X + Z*Z, Y*Y + W*W, Z*Z + X*X, W*W + Y*Y)
	Temp2 = _mm_shuffle_ps(Temp1, Temp1, SHUFFLEMASK(1, 2, 3, 0));	// Rotate left 4 bytes (Y,Z,W,X).
	return VectorAdd(Temp1, Temp2);								// (X*X + Z*Z + Y*Y + W*W, Y*Y + W*W + Z*Z + X*X, Z*Z + X*X + W*W + Y*Y, W*W + Y*Y + X*X + Z*Z)
}

#define VectorCompareEQ( Vec1, Vec2 )			_mm_cmpeq_ps( Vec1, Vec2 )
#define VectorCompareNE( Vec1, Vec2 )			_mm_cmpneq_ps( Vec1, Vec2 )
#define VectorCompareGT( Vec1, Vec2 )			_mm_cmpgt_ps( Vec1, Vec2 )
#define VectorCompareGE( Vec1, Vec2 )			_mm_cmpge_ps( Vec1, Vec2 )
#define VectorCompareLT( Vec1, Vec2 )			_mm_cmplt_ps( Vec1, Vec2 )
#define VectorCompareLE( Vec1, Vec2 )			_mm_cmple_ps( Vec1, Vec2 )
FORCEINLINE VectorRegister VectorSelect(const VectorRegister& Mask, const VectorRegister& Vec1, const VectorRegister& Vec2)
{
	return _mm_xor_ps(Vec2, _mm_and_ps(Mask, _mm_xor_ps(Vec1, Vec2)));
}
#define VectorBitwiseOr(Vec1, Vec2)	_mm_or_ps(Vec1, Vec2)
#define VectorBitwiseAnd(Vec1, Vec2) _mm_and_ps(Vec1, Vec2)
#define VectorBitwiseXor(Vec1, Vec2) _mm_xor_ps(Vec1, Vec2)

#define VectorMask_LT( Vec1, Vec2 )			_mm_cmplt_ps(Vec1, Vec2)
#define VectorMask_LE( Vec1, Vec2 )			_mm_cmple_ps(Vec1, Vec2)
#define VectorMask_GT( Vec1, Vec2 )			_mm_cmpgt_ps(Vec1, Vec2)
#define VectorMask_GE( Vec1, Vec2 )			_mm_cmpge_ps(Vec1, Vec2)
#define VectorMask_EQ( Vec1, Vec2 )			_mm_cmpeq_ps(Vec1, Vec2)
#define VectorMask_NE( Vec1, Vec2 )			_mm_cmpneq_ps(Vec1, Vec2)

#define VectorMaskBits( VecMask )			_mm_movemask_ps( VecMask )


FORCEINLINE VectorRegister VectorCross(const VectorRegister& Vec1, const VectorRegister& Vec2)
{
	VectorRegister A_YZXW = _mm_shuffle_ps(Vec1, Vec1, SHUFFLEMASK(1, 2, 0, 3));
	VectorRegister B_ZXYW = _mm_shuffle_ps(Vec2, Vec2, SHUFFLEMASK(2, 0, 1, 3));
	VectorRegister A_ZXYW = _mm_shuffle_ps(Vec1, Vec1, SHUFFLEMASK(2, 0, 1, 3));
	VectorRegister B_YZXW = _mm_shuffle_ps(Vec2, Vec2, SHUFFLEMASK(1, 2, 0, 3));
	return VectorSubtract(VectorMultiply(A_YZXW, B_ZXYW), VectorMultiply(A_ZXYW, B_YZXW));
}

FORCEINLINE VectorRegister VectorPow(const VectorRegister& Base, const VectorRegister& Exponent)
{
	//@TODO: Optimize
	union { VectorRegister v; float f[4]; } B, E;
	B.v = Base;
	E.v = Exponent;
	return _mm_setr_ps(powf(B.f[0], E.f[0]), powf(B.f[1], E.f[1]), powf(B.f[2], E.f[2]), powf(B.f[3], E.f[3]));
}

#define VectorReciprocalSqrt(Vec)		_mm_rsqrt_ps( Vec )
#define VectorReciprocal(Vec)			_mm_rcp_ps(Vec)
FORCEINLINE VectorRegister VectorReciprocalLen(const VectorRegister& Vector)
{
	VectorRegister RecipLen = VectorDot4(Vector, Vector);
	return VectorReciprocalSqrt(RecipLen);
}
FORCEINLINE VectorRegister VectorReciprocalSqrtAccurate(const VectorRegister& Vec)
{
	// Perform two passes of Newton-Raphson iteration on the hardware estimate
	//    v^-0.5 = x
	// => x^2 = v^-1
	// => 1/(x^2) = v
	// => F(x) = x^-2 - v
	//    F'(x) = -2x^-3

	//    x1 = x0 - F(x0)/F'(x0)
	// => x1 = x0 + 0.5 * (x0^-2 - Vec) * x0^3
	// => x1 = x0 + 0.5 * (x0 - Vec * x0^3)
	// => x1 = x0 + x0 * (0.5 - 0.5 * Vec * x0^2)

	const VectorRegister OneHalf = GlobalVectorConstants::FloatOneHalf;
	const VectorRegister VecDivBy2 = VectorMultiply(Vec, OneHalf);

	// Initial estimate
	const VectorRegister x0 = VectorReciprocalSqrt(Vec);

	// First iteration
	VectorRegister x1 = VectorMultiply(x0, x0);
	x1 = VectorSubtract(OneHalf, VectorMultiply(VecDivBy2, x1));
	x1 = VectorMultiplyAdd(x0, x1, x0);

	// Second iteration
	VectorRegister x2 = VectorMultiply(x1, x1);
	x2 = VectorSubtract(OneHalf, VectorMultiply(VecDivBy2, x2));
	x2 = VectorMultiplyAdd(x1, x2, x1);

	return x2;
}
FORCEINLINE VectorRegister VectorReciprocalAccurate(const VectorRegister& Vec)
{
	// Perform two passes of Newton-Raphson iteration on the hardware estimate
	//   x1 = x0 - f(x0) / f'(x0)
	//
	//    1 / Vec = x
	// => x * Vec = 1 
	// => F(x) = x * Vec - 1
	//    F'(x) = Vec
	// => x1 = x0 - (x0 * Vec - 1) / Vec
	//
	// Since 1/Vec is what we're trying to solve, use an estimate for it, x0
	// => x1 = x0 - (x0 * Vec - 1) * x0 = 2 * x0 - Vec * x0^2 

	// Initial estimate
	const VectorRegister x0 = VectorReciprocal(Vec);

	// First iteration
	const VectorRegister x0Squared = VectorMultiply(x0, x0);
	const VectorRegister x0Times2 = VectorAdd(x0, x0);
	const VectorRegister x1 = VectorSubtract(x0Times2, VectorMultiply(Vec, x0Squared));

	// Second iteration
	const VectorRegister x1Squared = VectorMultiply(x1, x1);
	const VectorRegister x1Times2 = VectorAdd(x1, x1);
	const VectorRegister x2 = VectorSubtract(x1Times2, VectorMultiply(Vec, x1Squared));

	return x2;
}

FORCEINLINE VectorRegister VectorNormalize(const VectorRegister& Vector)
{
	return VectorMultiply(Vector, VectorReciprocalLen(Vector));
}

FORCEINLINE VectorRegister VectorTransformVector(const VectorRegister&  VecP, const void* MatrixM)
{
	const VectorRegister *M = (const VectorRegister *)MatrixM;
	VectorRegister VTempX, VTempY, VTempZ, VTempW;

	VTempX = VectorReplicate(VecP, 0);
	VTempY = VectorReplicate(VecP, 1);
	VTempZ = VectorReplicate(VecP, 2);
	VTempW = VectorReplicate(VecP, 3);
	VTempX = VectorMultiply(VTempX, M[0]);
	VTempY = VectorMultiply(VTempY, M[1]);
	VTempZ = VectorMultiply(VTempZ, M[2]);
	VTempW = VectorMultiply(VTempW, M[3]);
	VTempX = VectorAdd(VTempX, VTempY);
	VTempZ = VectorAdd(VTempZ, VTempW);
	VTempX = VectorAdd(VTempX, VTempZ);

	return VTempX;
}

FORCEINLINE void VectorMatrixMultiply(void *Result, const void* Matrix1, const void* Matrix2)
{
	const VectorRegister *A = (const VectorRegister *)Matrix1;
	const VectorRegister *B = (const VectorRegister *)Matrix2;
	VectorRegister *R = (VectorRegister *)Result;
	VectorRegister Temp, R0, R1, R2, R3;

	// First row of result (Matrix1[0] * Matrix2).
	Temp = VectorMultiply(VectorReplicate(A[0], 0), B[0]);
	Temp = VectorMultiplyAdd(VectorReplicate(A[0], 1), B[1], Temp);
	Temp = VectorMultiplyAdd(VectorReplicate(A[0], 2), B[2], Temp);
	R0 = VectorMultiplyAdd(VectorReplicate(A[0], 3), B[3], Temp);

	// Second row of result (Matrix1[1] * Matrix2).
	Temp = VectorMultiply(VectorReplicate(A[1], 0), B[0]);
	Temp = VectorMultiplyAdd(VectorReplicate(A[1], 1), B[1], Temp);
	Temp = VectorMultiplyAdd(VectorReplicate(A[1], 2), B[2], Temp);
	R1 = VectorMultiplyAdd(VectorReplicate(A[1], 3), B[3], Temp);

	// Third row of result (Matrix1[2] * Matrix2).
	Temp = VectorMultiply(VectorReplicate(A[2], 0), B[0]);
	Temp = VectorMultiplyAdd(VectorReplicate(A[2], 1), B[1], Temp);
	Temp = VectorMultiplyAdd(VectorReplicate(A[2], 2), B[2], Temp);
	R2 = VectorMultiplyAdd(VectorReplicate(A[2], 3), B[3], Temp);

	// Fourth row of result (Matrix1[3] * Matrix2).
	Temp = VectorMultiply(VectorReplicate(A[3], 0), B[0]);
	Temp = VectorMultiplyAdd(VectorReplicate(A[3], 1), B[1], Temp);
	Temp = VectorMultiplyAdd(VectorReplicate(A[3], 2), B[2], Temp);
	R3 = VectorMultiplyAdd(VectorReplicate(A[3], 3), B[3], Temp);

	// Store result
	R[0] = R0;
	R[1] = R1;
	R[2] = R2;
	R[3] = R3;
}

FORCEINLINE void VectorMatrixInverse(void* DstMatrix, const void* SrcMatrix)
{
	typedef float Float4x4[4][4];
	const Float4x4& M = *((const Float4x4*)SrcMatrix);
	Float4x4 Result;
	float Det[4];
	Float4x4 Tmp;

	Tmp[0][0] = M[2][2] * M[3][3] - M[2][3] * M[3][2];
	Tmp[0][1] = M[1][2] * M[3][3] - M[1][3] * M[3][2];
	Tmp[0][2] = M[1][2] * M[2][3] - M[1][3] * M[2][2];

	Tmp[1][0] = M[2][2] * M[3][3] - M[2][3] * M[3][2];
	Tmp[1][1] = M[0][2] * M[3][3] - M[0][3] * M[3][2];
	Tmp[1][2] = M[0][2] * M[2][3] - M[0][3] * M[2][2];

	Tmp[2][0] = M[1][2] * M[3][3] - M[1][3] * M[3][2];
	Tmp[2][1] = M[0][2] * M[3][3] - M[0][3] * M[3][2];
	Tmp[2][2] = M[0][2] * M[1][3] - M[0][3] * M[1][2];

	Tmp[3][0] = M[1][2] * M[2][3] - M[1][3] * M[2][2];
	Tmp[3][1] = M[0][2] * M[2][3] - M[0][3] * M[2][2];
	Tmp[3][2] = M[0][2] * M[1][3] - M[0][3] * M[1][2];

	Det[0] = M[1][1] * Tmp[0][0] - M[2][1] * Tmp[0][1] + M[3][1] * Tmp[0][2];
	Det[1] = M[0][1] * Tmp[1][0] - M[2][1] * Tmp[1][1] + M[3][1] * Tmp[1][2];
	Det[2] = M[0][1] * Tmp[2][0] - M[1][1] * Tmp[2][1] + M[3][1] * Tmp[2][2];
	Det[3] = M[0][1] * Tmp[3][0] - M[1][1] * Tmp[3][1] + M[2][1] * Tmp[3][2];

	float Determinant = M[0][0] * Det[0] - M[1][0] * Det[1] + M[2][0] * Det[2] - M[3][0] * Det[3];
	const float	RDet = 1.0f / Determinant;

	Result[0][0] = RDet * Det[0];
	Result[0][1] = -RDet * Det[1];
	Result[0][2] = RDet * Det[2];
	Result[0][3] = -RDet * Det[3];
	Result[1][0] = -RDet * (M[1][0] * Tmp[0][0] - M[2][0] * Tmp[0][1] + M[3][0] * Tmp[0][2]);
	Result[1][1] = RDet * (M[0][0] * Tmp[1][0] - M[2][0] * Tmp[1][1] + M[3][0] * Tmp[1][2]);
	Result[1][2] = -RDet * (M[0][0] * Tmp[2][0] - M[1][0] * Tmp[2][1] + M[3][0] * Tmp[2][2]);
	Result[1][3] = RDet * (M[0][0] * Tmp[3][0] - M[1][0] * Tmp[3][1] + M[2][0] * Tmp[3][2]);
	Result[2][0] = RDet * (
		M[1][0] * (M[2][1] * M[3][3] - M[2][3] * M[3][1]) -
		M[2][0] * (M[1][1] * M[3][3] - M[1][3] * M[3][1]) +
		M[3][0] * (M[1][1] * M[2][3] - M[1][3] * M[2][1])
		);
	Result[2][1] = -RDet * (
		M[0][0] * (M[2][1] * M[3][3] - M[2][3] * M[3][1]) -
		M[2][0] * (M[0][1] * M[3][3] - M[0][3] * M[3][1]) +
		M[3][0] * (M[0][1] * M[2][3] - M[0][3] * M[2][1])
		);
	Result[2][2] = RDet * (
		M[0][0] * (M[1][1] * M[3][3] - M[1][3] * M[3][1]) -
		M[1][0] * (M[0][1] * M[3][3] - M[0][3] * M[3][1]) +
		M[3][0] * (M[0][1] * M[1][3] - M[0][3] * M[1][1])
		);
	Result[2][3] = -RDet * (
		M[0][0] * (M[1][1] * M[2][3] - M[1][3] * M[2][1]) -
		M[1][0] * (M[0][1] * M[2][3] - M[0][3] * M[2][1]) +
		M[2][0] * (M[0][1] * M[1][3] - M[0][3] * M[1][1])
		);
	Result[3][0] = -RDet * (
		M[1][0] * (M[2][1] * M[3][2] - M[2][2] * M[3][1]) -
		M[2][0] * (M[1][1] * M[3][2] - M[1][2] * M[3][1]) +
		M[3][0] * (M[1][1] * M[2][2] - M[1][2] * M[2][1])
		);
	Result[3][1] = RDet * (
		M[0][0] * (M[2][1] * M[3][2] - M[2][2] * M[3][1]) -
		M[2][0] * (M[0][1] * M[3][2] - M[0][2] * M[3][1]) +
		M[3][0] * (M[0][1] * M[2][2] - M[0][2] * M[2][1])
		);
	Result[3][2] = -RDet * (
		M[0][0] * (M[1][1] * M[3][2] - M[1][2] * M[3][1]) -
		M[1][0] * (M[0][1] * M[3][2] - M[0][2] * M[3][1]) +
		M[3][0] * (M[0][1] * M[1][2] - M[0][2] * M[1][1])
		);
	Result[3][3] = RDet * (
		M[0][0] * (M[1][1] * M[2][2] - M[1][2] * M[2][1]) -
		M[1][0] * (M[0][1] * M[2][2] - M[0][2] * M[2][1]) +
		M[2][0] * (M[0][1] * M[1][2] - M[0][2] * M[1][1])
		);

	memcpy(DstMatrix, &Result, 16 * sizeof(float));
}

FORCEINLINE VectorRegister VectorQuaternionMultiply2(const VectorRegister& Quat1, const VectorRegister& Quat2)
{
	VectorRegister Result = VectorMultiply(VectorReplicate(Quat1, 3), Quat2);
	Result = VectorMultiplyAdd(VectorMultiply(VectorReplicate(Quat1, 0), VectorSwizzle(Quat2, 3, 2, 1, 0)), GlobalVectorConstants::QMULTI_SIGN_MASK0, Result);
	Result = VectorMultiplyAdd(VectorMultiply(VectorReplicate(Quat1, 1), VectorSwizzle(Quat2, 2, 3, 0, 1)), GlobalVectorConstants::QMULTI_SIGN_MASK1, Result);
	Result = VectorMultiplyAdd(VectorMultiply(VectorReplicate(Quat1, 2), VectorSwizzle(Quat2, 1, 0, 3, 2)), GlobalVectorConstants::QMULTI_SIGN_MASK2, Result);

	return Result;
}

FORCEINLINE void VectorQuaternionMultiply(void* RESTRICT Result, const void* RESTRICT Quat1, const void* RESTRICT Quat2)
{
	*((VectorRegister *)Result) = VectorQuaternionMultiply2(*((const VectorRegister *)Quat1), *((const VectorRegister *)Quat2));
}

FORCEINLINE VectorRegister VectorDivide(const VectorRegister& Vec1, const VectorRegister& Vec2)
{
	return _mm_div_ps(Vec1, Vec2);
}

// Returns true if the vector contains a component that is either NAN or +/-infinite.
FORCEINLINE bool VectorContainsNaNOrInfinite(const VectorRegister& Vec)
{
	// https://en.wikipedia.org/wiki/IEEE_754-1985
	// Infinity is represented with all exponent bits set, with the correct sign bit.
	// NaN is represented with all exponent bits set, plus at least one fraction/significand bit set.
	// This means finite values will not have all exponent bits set, so check against those bits.

	union { float F; uint32 U; } InfUnion;
	InfUnion.U = 0x7F800000;
	const float Inf = InfUnion.F;
	const VectorRegister FloatInfinity = MakeVectorRegister(Inf, Inf, Inf, Inf);

	// Mask off Exponent
	VectorRegister ExpTest = VectorBitwiseAnd(Vec, FloatInfinity);
	// Compare to full exponent. If any are full exponent (not finite), the signs copied to the mask are non-zero, otherwise it's zero and finite.
	bool IsFinite = VectorMaskBits(VectorCompareEQ(ExpTest, FloatInfinity)) == 0;
	return !IsFinite;
}

FORCEINLINE VectorRegister VectorTruncate(const VectorRegister& X)
{
	return _mm_cvtepi32_ps(_mm_cvttps_epi32(X));
}

FORCEINLINE VectorRegister VectorFractional(const VectorRegister& X)
{
	return VectorSubtract(X, VectorTruncate(X));
}

FORCEINLINE VectorRegister VectorCeil(const VectorRegister& X)
{
	VectorRegister Trunc = VectorTruncate(X);
	VectorRegister PosMask = VectorCompareGE(X, GlobalVectorConstants::FloatZero);
	VectorRegister Add = VectorSelect(PosMask, GlobalVectorConstants::FloatOne, (GlobalVectorConstants::FloatZero));
	return VectorAdd(Trunc, Add);
}

FORCEINLINE VectorRegister VectorFloor(const VectorRegister& X)
{
	VectorRegister Trunc = VectorTruncate(X);
	VectorRegister PosMask = VectorCompareGE(X, (GlobalVectorConstants::FloatZero));
	VectorRegister Sub = VectorSelect(PosMask, (GlobalVectorConstants::FloatZero), (GlobalVectorConstants::FloatOne));
	return VectorSubtract(Trunc, Sub);
}

FORCEINLINE VectorRegister VectorMod(const VectorRegister& X, const VectorRegister& Y)
{
	VectorRegister Div = VectorDivide(X, Y);
	// Floats where abs(f) >= 2^23 have no fractional portion, and larger values would overflow VectorTruncate.
	VectorRegister NoFractionMask = VectorCompareGE(VectorAbs(Div), GlobalVectorConstants::FloatNonFractional);
	VectorRegister Temp = VectorSelect(NoFractionMask, Div, VectorTruncate(Div));
	VectorRegister Result = VectorSubtract(X, VectorMultiply(Y, Temp));
	// Clamp to [-AbsY, AbsY] because of possible failures for very large numbers (>1e10) due to precision loss.
	VectorRegister AbsY = VectorAbs(Y);
	return VectorMax(VectorNegate(AbsY), VectorMin(Result, AbsY));
}

FORCEINLINE VectorRegister VectorSign(const VectorRegister& X)
{
	VectorRegister Mask = VectorCompareGE(X, (GlobalVectorConstants::FloatZero));
	return VectorSelect(Mask, (GlobalVectorConstants::FloatOne), (GlobalVectorConstants::FloatMinusOne));
}

FORCEINLINE VectorRegister VectorStep(const VectorRegister& X)
{
	VectorRegister Mask = VectorCompareGE(X, (GlobalVectorConstants::FloatZero));
	return VectorSelect(Mask, (GlobalVectorConstants::FloatOne), (GlobalVectorConstants::FloatZero));
}

//TODO: Vectorize
FORCEINLINE VectorRegister VectorExp(const VectorRegister& X)
{
	return MakeVectorRegister(Math::Exp(VectorGetComponent(X, 0)), Math::Exp(VectorGetComponent(X, 1)), Math::Exp(VectorGetComponent(X, 2)), Math::Exp(VectorGetComponent(X, 3)));
}

//TODO: Vectorize
FORCEINLINE VectorRegister VectorExp2(const VectorRegister& X)
{
	return MakeVectorRegister(Math::Exp2(VectorGetComponent(X, 0)), Math::Exp2(VectorGetComponent(X, 1)), Math::Exp2(VectorGetComponent(X, 2)), Math::Exp2(VectorGetComponent(X, 3)));
}

//TODO: Vectorize
FORCEINLINE VectorRegister VectorLog(const VectorRegister& X)
{
	return MakeVectorRegister(Math::Loge(VectorGetComponent(X, 0)), Math::Loge(VectorGetComponent(X, 1)), Math::Loge(VectorGetComponent(X, 2)), Math::Loge(VectorGetComponent(X, 3)));
}

//TODO: Vectorize
FORCEINLINE VectorRegister VectorLog2(const VectorRegister& X)
{
	return MakeVectorRegister(Math::Log2(VectorGetComponent(X, 0)), Math::Log2(VectorGetComponent(X, 1)), Math::Log2(VectorGetComponent(X, 2)), Math::Log2(VectorGetComponent(X, 3)));
}

#define VectorMergeVecXYZ_VecW( VecXYZ, VecW )	VectorSelect( GlobalVectorConstants::XYZMask, VecXYZ, VecW )
#define VectorLoadByte4( Ptr )			_mm_cvtepi32_ps(_mm_unpacklo_epi16(_mm_unpacklo_epi8(_mm_cvtsi32_si128(*(int32*)Ptr), _mm_setzero_si128()), _mm_setzero_si128()))
FORCEINLINE VectorRegister VectorLoadSignedByte4(const void* Ptr)
{
	auto Temp = _mm_unpacklo_epi16(_mm_unpacklo_epi8(_mm_cvtsi32_si128(*(int32*)Ptr), _mm_setzero_si128()), _mm_setzero_si128());
	auto Mask = _mm_cmpgt_epi32(Temp, _mm_set1_epi32(127));
	auto Comp = _mm_and_si128(Mask, _mm_set1_epi32(~127));
	return _mm_cvtepi32_ps(_mm_or_si128(Comp, Temp));
}
FORCEINLINE VectorRegister VectorLoadByte4Reverse(void* Ptr)
{
	VectorRegister Temp = VectorLoadByte4(Ptr);
	return _mm_shuffle_ps(Temp, Temp, SHUFFLEMASK(3, 2, 1, 0));
}
FORCEINLINE void VectorStoreByte4(const VectorRegister& Vec, void* Ptr)
{
	// Looks complex but is really quite straightforward:
	// Convert 4x floats to 4x 32-bit ints, then pack into 4x 16-bit ints, then into 4x 8-bit unsigned ints, then store as a 32-bit value
	*(int32*)Ptr = _mm_cvtsi128_si32(_mm_packus_epi16(_mm_packs_epi32(_mm_cvttps_epi32(Vec), _mm_setzero_si128()), _mm_setzero_si128()));
}
FORCEINLINE void VectorStoreSignedByte4(const VectorRegister& Vec, void* Ptr)
{
	// Looks complex but is really quite straightforward:
	// Convert 4x floats to 4x 32-bit ints, then pack into 4x 16-bit ints, then into 4x 8-bit unsigned ints, then store as a 32-bit value
	*(int32*)Ptr = _mm_cvtsi128_si32(_mm_packs_epi16(_mm_packs_epi32(_mm_cvttps_epi32(Vec), _mm_setzero_si128()), _mm_setzero_si128()));
}
FORCEINLINE VectorRegister VectorLoadURGB10A2N(void* Ptr)
{
	VectorRegister Tmp;

	Tmp = _mm_and_ps(_mm_load_ps1((const float *)Ptr), MakeVectorRegister(0x3FFu, 0x3FFu << 10, 0x3FFu << 20, 0x3u << 30));
	Tmp = _mm_xor_ps(Tmp, VectorSet(0, 0, 0, 0x80000000));
	Tmp = _mm_cvtepi32_ps(*(const VectorRegisterInt*)&Tmp);
	Tmp = _mm_add_ps(Tmp, VectorSet(0, 0, 0, 32768.0f*65536.0f));
	Tmp = _mm_mul_ps(Tmp, VectorSet(1.0f / 1023.0f, 1.0f / (1023.0f*1024.0f), 1.0f / (1023.0f*1024.0f*1024.0f), 1.0f / (3.0f*1024.0f*1024.0f*1024.0f)));

	return Tmp;
}
FORCEINLINE void VectorStoreURGB10A2N(const VectorRegister& Vec, void* Ptr)
{
	VectorRegister Tmp;
	Tmp = _mm_max_ps(Vec, MakeVectorRegister(0.0f, 0.0f, 0.0f, 0.0f));
	Tmp = _mm_min_ps(Tmp, MakeVectorRegister(1.0f, 1.0f, 1.0f, 1.0f));
	Tmp = _mm_mul_ps(Tmp, MakeVectorRegister(1023.0f, 1023.0f*1024.0f*0.5f, 1023.0f*1024.0f*1024.0f, 3.0f*1024.0f*1024.0f*1024.0f*0.5f));

	VectorRegisterInt TmpI;
	TmpI = _mm_cvttps_epi32(Tmp);
	TmpI = _mm_and_si128(TmpI, MakeVectorRegisterInt(0x3FFu, 0x3FFu << (10 - 1), 0x3FFu << 20, 0x3u << (30 - 1)));

	VectorRegisterInt TmpI2;
	TmpI2 = _mm_shuffle_epi32(TmpI, _MM_SHUFFLE(3, 2, 3, 2));
	TmpI = _mm_or_si128(TmpI, TmpI2);

	TmpI2 = _mm_shuffle_epi32(TmpI, _MM_SHUFFLE(1, 1, 1, 1));
	TmpI2 = _mm_add_epi32(TmpI2, TmpI2);
	TmpI = _mm_or_si128(TmpI, TmpI2);

	_mm_store_ss((float *)Ptr, *(const VectorRegister*)&TmpI);
}
#define VectorLoadURGBA16N( Ptr ) _mm_cvtepi32_ps(_mm_unpacklo_epi16(_mm_loadl_epi64((const __m128i*)Ptr), _mm_setzero_si128()))
FORCEINLINE VectorRegister VectorLoadSRGBA16N(const void* Ptr)
{
	auto Temp = _mm_unpacklo_epi16(_mm_loadl_epi64((const __m128i*)Ptr), _mm_setzero_si128());
	auto Mask = _mm_cmpgt_epi32(Temp, _mm_set1_epi32(32767));
	auto Comp = _mm_and_si128(Mask, _mm_set1_epi32(~32767));
	return _mm_cvtepi32_ps(_mm_or_si128(Comp, Temp));
}
FORCEINLINE void VectorStoreURGBA16N(const VectorRegister& Vec, void* Ptr)
{

	VectorRegister Tmp;
	Tmp = _mm_max_ps(Vec, MakeVectorRegister(0.0f, 0.0f, 0.0f, 0.0f));
	Tmp = _mm_min_ps(Tmp, MakeVectorRegister(1.0f, 1.0f, 1.0f, 1.0f));
	Tmp = _mm_mul_ps(Tmp, MakeVectorRegister(65535.0f, 65535.0f, 65535.0f, 65535.0f));

	VectorRegisterInt TmpI = _mm_cvtps_epi32(Tmp);

	uint16* Out = (uint16*)Ptr;
	Out[0] = static_cast<int16>(_mm_extract_epi16(TmpI, 0));
	Out[1] = static_cast<int16>(_mm_extract_epi16(TmpI, 2));
	Out[2] = static_cast<int16>(_mm_extract_epi16(TmpI, 4));
	Out[3] = static_cast<int16>(_mm_extract_epi16(TmpI, 6));
}
#define VectorAnyGreaterThan( Vec1, Vec2 )		_mm_movemask_ps( _mm_cmpgt_ps(Vec1, Vec2) )
#define VectorResetFloatRegisters()
#define VectorGetControlRegister()		_mm_getcsr()
#define	VectorSetControlRegister(ControlStatus) _mm_setcsr( ControlStatus )
#define VECTOR_ROUND_TOWARD_ZERO		_MM_ROUND_TOWARD_ZERO

/**
 * Using "static const float ..." or "static const VectorRegister ..." in functions creates the branch and code to construct those constants.
 * Doing this in FORCEINLINE not only means you introduce a branch per static, but you bloat the FORCEINLINEd code immensely.
 * Defining these constants at the global scope causes them to be created at startup, and avoids the cost at the function level.
 * Doing it at the function level is okay for anything that is a simple "const float", but usage of "sqrt()" here forces actual function calls.
 */
namespace VectorSinConstantsSSE
{
	static const float p = 0.225f;
	static const float a = (16 * sqrtf(p));
	static const float b = ((1 - p) / sqrtf(p));
	static const VectorRegister A = MakeVectorRegister(a, a, a, a);
	static const VectorRegister B = MakeVectorRegister(b, b, b, b);
}

FORCEINLINE VectorRegister VectorSin(const VectorRegister& X)
{
	//Sine approximation using a squared parabola restrained to f(0) = 0, f(PI) = 0, f(PI/2) = 1.
	//based on a good discussion here http://forum.devmaster.net/t/fast-and-accurate-sine-cosine/9648
	//After approx 2.5 million tests comparing to sin(): 
	//Average error of 0.000128
	//Max error of 0.001091

	VectorRegister y = VectorMultiply(X, GlobalVectorConstants::OneOverTwoPi);
	y = VectorSubtract(y, VectorFloor(VectorAdd(y, GlobalVectorConstants::FloatOneHalf)));
	y = VectorMultiply(VectorSinConstantsSSE::A, VectorMultiply(y, VectorSubtract(GlobalVectorConstants::FloatOneHalf, VectorAbs(y))));
	return VectorMultiply(y, VectorAdd(VectorSinConstantsSSE::B, VectorAbs(y)));
}

FORCEINLINE VectorRegister VectorCos(const VectorRegister& X)
{
	return VectorSin(VectorAdd(X, GlobalVectorConstants::PiByTwo));
}


/**
* Computes the sine and cosine of each component of a Vector.
*
* @param VSinAngles	VectorRegister Pointer to where the Sin result should be stored
* @param VCosAngles	VectorRegister Pointer to where the Cos result should be stored
* @param VAngles VectorRegister Pointer to the input angles
*/
FORCEINLINE void VectorSinCos(VectorRegister* RESTRICT VSinAngles, VectorRegister* RESTRICT VCosAngles, const VectorRegister* RESTRICT VAngles)
{
	// Map to [-pi, pi]
	// X = A - 2pi * round(A/2pi)
	// Note the round(), not truncate(). In this case round() can round halfway cases using round-to-nearest-even OR round-to-nearest.

	// Quotient = round(A/2pi)
	VectorRegister Quotient = VectorMultiply(*VAngles, GlobalVectorConstants::OneOverTwoPi);
	Quotient = _mm_cvtepi32_ps(_mm_cvtps_epi32(Quotient)); // round to nearest even is the default rounding mode but that's fine here.
	// X = A - 2pi * Quotient
	VectorRegister X = VectorSubtract(*VAngles, VectorMultiply(GlobalVectorConstants::TwoPi, Quotient));

	// Map in [-pi/2,pi/2]
	VectorRegister sign = VectorBitwiseAnd(X, GlobalVectorConstants::SignBit);
	VectorRegister c = VectorBitwiseOr(GlobalVectorConstants::Pi, sign);  // pi when x >= 0, -pi when x < 0
	VectorRegister absx = VectorAbs(X);
	VectorRegister rflx = VectorSubtract(c, X);
	VectorRegister comp = VectorCompareGT(absx, GlobalVectorConstants::PiByTwo);
	X = VectorSelect(comp, rflx, X);
	sign = VectorSelect(comp, GlobalVectorConstants::FloatMinusOne, GlobalVectorConstants::FloatOne);

	const VectorRegister XSquared = VectorMultiply(X, X);

	// 11-degree minimax approximation
	//*ScalarSin = (((((-2.3889859e-08f * y2 + 2.7525562e-06f) * y2 - 0.00019840874f) * y2 + 0.0083333310f) * y2 - 0.16666667f) * y2 + 1.0f) * y;
	const VectorRegister SinCoeff0 = MakeVectorRegister(1.0f, -0.16666667f, 0.0083333310f, -0.00019840874f);
	const VectorRegister SinCoeff1 = MakeVectorRegister(2.7525562e-06f, -2.3889859e-08f, /*unused*/ 0.f, /*unused*/ 0.f);

	VectorRegister S;
	S = VectorReplicate(SinCoeff1, 1);
	S = VectorMultiplyAdd(XSquared, S, VectorReplicate(SinCoeff1, 0));
	S = VectorMultiplyAdd(XSquared, S, VectorReplicate(SinCoeff0, 3));
	S = VectorMultiplyAdd(XSquared, S, VectorReplicate(SinCoeff0, 2));
	S = VectorMultiplyAdd(XSquared, S, VectorReplicate(SinCoeff0, 1));
	S = VectorMultiplyAdd(XSquared, S, VectorReplicate(SinCoeff0, 0));
	*VSinAngles = VectorMultiply(S, X);

	// 10-degree minimax approximation
	//*ScalarCos = sign * (((((-2.6051615e-07f * y2 + 2.4760495e-05f) * y2 - 0.0013888378f) * y2 + 0.041666638f) * y2 - 0.5f) * y2 + 1.0f);
	const VectorRegister CosCoeff0 = MakeVectorRegister(1.0f, -0.5f, 0.041666638f, -0.0013888378f);
	const VectorRegister CosCoeff1 = MakeVectorRegister(2.4760495e-05f, -2.6051615e-07f, /*unused*/ 0.f, /*unused*/ 0.f);

	VectorRegister C;
	C = VectorReplicate(CosCoeff1, 1);
	C = VectorMultiplyAdd(XSquared, C, VectorReplicate(CosCoeff1, 0));
	C = VectorMultiplyAdd(XSquared, C, VectorReplicate(CosCoeff0, 3));
	C = VectorMultiplyAdd(XSquared, C, VectorReplicate(CosCoeff0, 2));
	C = VectorMultiplyAdd(XSquared, C, VectorReplicate(CosCoeff0, 1));
	C = VectorMultiplyAdd(XSquared, C, VectorReplicate(CosCoeff0, 0));
	*VCosAngles = VectorMultiply(C, sign);
}

FORCEINLINE VectorRegister VectorTan(const VectorRegister& X)
{
	return MakeVectorRegister(Math::Tan(VectorGetComponent(X, 0)), Math::Tan(VectorGetComponent(X, 1)), Math::Tan(VectorGetComponent(X, 2)), Math::Tan(VectorGetComponent(X, 3)));
}

FORCEINLINE VectorRegister VectorASin(const VectorRegister& X)
{
	return MakeVectorRegister(Math::Asin(VectorGetComponent(X, 0)), Math::Asin(VectorGetComponent(X, 1)), Math::Asin(VectorGetComponent(X, 2)), Math::Asin(VectorGetComponent(X, 3)));
}

FORCEINLINE VectorRegister VectorACos(const VectorRegister& X)
{
	return MakeVectorRegister(Math::Acos(VectorGetComponent(X, 0)), Math::Acos(VectorGetComponent(X, 1)), Math::Acos(VectorGetComponent(X, 2)), Math::Acos(VectorGetComponent(X, 3)));
}

FORCEINLINE VectorRegister VectorATan(const VectorRegister& X)
{
	return MakeVectorRegister(Math::Atan(VectorGetComponent(X, 0)), Math::Atan(VectorGetComponent(X, 1)), Math::Atan(VectorGetComponent(X, 2)), Math::Atan(VectorGetComponent(X, 3)));
}

FORCEINLINE VectorRegister VectorATan2(const VectorRegister& X, const VectorRegister& Y)
{
	return MakeVectorRegister(Math::Atan2(VectorGetComponent(X, 0), VectorGetComponent(Y, 0)),
		Math::Atan2(VectorGetComponent(X, 1), VectorGetComponent(Y, 1)),
		Math::Atan2(VectorGetComponent(X, 2), VectorGetComponent(Y, 2)),
		Math::Atan2(VectorGetComponent(X, 3), VectorGetComponent(Y, 3)));
}

// To be continued...


//////////////////////////////////////////////////////////////////////////
//Integer ops

//Bitwise
/** = a & b */
#define VectorIntAnd(A, B)		_mm_and_si128(A, B)
/** = a | b */
#define VectorIntOr(A, B)		_mm_or_si128(A, B)
/** = a ^ b */
#define VectorIntXor(A, B)		_mm_xor_si128(A, B)
/** = (~a) & b */
#define VectorIntAndNot(A, B)	_mm_andnot_si128(A, B)
/** = ~a */
#define VectorIntNot(A)	_mm_xor_si128(A, GlobalVectorConstants::IntAllMask)

//Comparison
#define VectorIntCompareEQ(A, B)	_mm_cmpeq_epi32(A,B)
#define VectorIntCompareNEQ(A, B)	VectorIntNot(_mm_cmpeq_epi32(A,B))
#define VectorIntCompareGT(A, B)	_mm_cmpgt_epi32(A,B)
#define VectorIntCompareLT(A, B)	_mm_cmplt_epi32(A,B)
#define VectorIntCompareGE(A, B)	VectorIntNot(VectorIntCompareLT(A,B))
#define VectorIntCompareLE(A, B)	VectorIntNot(VectorIntCompareGT(A,B))


FORCEINLINE VectorRegisterInt VectorIntSelect(const VectorRegisterInt& Mask, const VectorRegisterInt& Vec1, const VectorRegisterInt& Vec2)
{
	return _mm_xor_si128(Vec2, _mm_and_si128(Mask, _mm_xor_si128(Vec1, Vec2)));
}

//Arithmetic
#define VectorIntAdd(A, B)	_mm_add_epi32(A, B)
#define VectorIntSubtract(A, B)	_mm_sub_epi32(A, B)

FORCEINLINE VectorRegisterInt VectorIntMultiply(const VectorRegisterInt& A, const VectorRegisterInt& B)
{
	//SSE2 doesn't have a multiply op for 4 32bit ints. Ugh.
	__m128i Temp0 = _mm_mul_epu32(A, B);
	__m128i Temp1 = _mm_mul_epu32(_mm_srli_si128(A, 4), _mm_srli_si128(B, 4));
	return _mm_unpacklo_epi32(_mm_shuffle_epi32(Temp0, _MM_SHUFFLE(0, 0, 2, 0)), _mm_shuffle_epi32(Temp1, _MM_SHUFFLE(0, 0, 2, 0)));
}

#define VectorIntNegate(A) VectorIntSubtract( GlobalVectorConstants::IntZero, A)

FORCEINLINE VectorRegisterInt VectorIntMin(const VectorRegisterInt& A, const VectorRegisterInt& B)
{
	VectorRegisterInt Mask = VectorIntCompareLT(A, B);
	return VectorIntSelect(Mask, A, B);
}

FORCEINLINE VectorRegisterInt VectorIntMax(const VectorRegisterInt& A, const VectorRegisterInt& B)
{
	VectorRegisterInt Mask = VectorIntCompareGT(A, B);
	return VectorIntSelect(Mask, A, B);
}

FORCEINLINE VectorRegisterInt VectorIntAbs(const VectorRegisterInt& A)
{
	VectorRegisterInt Mask = VectorIntCompareGE(A, GlobalVectorConstants::IntZero);
	return VectorIntSelect(Mask, A, VectorIntNegate(A));
}

#define VectorIntSign(A) VectorIntSelect( VectorIntCompareGE(A, GlobalVectorConstants::IntZero), GlobalVectorConstants::IntOne, GlobalVectorConstants::IntMinusOne )

#define VectorIntToFloat(A) _mm_cvtepi32_ps(A)
#define VectorFloatToInt(A) _mm_cvttps_epi32(A)

//Loads and stores

/**
* Stores a vector to memory (aligned or unaligned).
*
* @param Vec	Vector to store
* @param Ptr	Memory pointer
*/
#define VectorIntStore( Vec, Ptr )			_mm_storeu_si128( (VectorRegisterInt*)(Ptr), Vec )

/**
* Loads 4 int32s from unaligned memory.
*
* @param Ptr	Unaligned memory pointer to the 4 int32s
* @return		VectorRegisterInt(Ptr[0], Ptr[1], Ptr[2], Ptr[3])
*/
#define VectorIntLoad( Ptr )				_mm_loadu_si128( (VectorRegisterInt*)(Ptr) )

/**
* Stores a vector to memory (aligned).
*
* @param Vec	Vector to store
* @param Ptr	Aligned Memory pointer
*/
#define VectorIntStoreAligned( Vec, Ptr )			_mm_store_si128( (VectorRegisterInt*)(Ptr), Vec )

/**
* Loads 4 int32s from aligned memory.
*
* @param Ptr	Aligned memory pointer to the 4 int32s
* @return		VectorRegisterInt(Ptr[0], Ptr[1], Ptr[2], Ptr[3])
*/
#define VectorIntLoadAligned( Ptr )				_mm_load_si128( (VectorRegisterInt*)(Ptr) )

/**
* Loads 1 int32 from unaligned memory into all components of a vector register.
*
* @param Ptr	Unaligned memory pointer to the 4 int32s
* @return		VectorRegisterInt(*Ptr, *Ptr, *Ptr, *Ptr)
*/
#define VectorIntLoad1( Ptr )	_mm_shuffle_epi32(_mm_loadu_si128((VectorRegisterInt*)Ptr),_MM_SHUFFLE(0,0,0,0))

};
