// this tiny game engine is rendered by CPU. I make it for learing the pipline and rasterization. This is the first step to realise my own game engine. so just do it!

#pragma once

#include "platform/Platform.h"
#include "MathUtility.h"
#define SIMD_ALIGNMENT (16)
#include "MathSSE.h"

#include "VectorCommon.h"

namespace LXF {

/** Vector that represents (1/255,1/255,1/255,1/255) */
extern const VectorRegister VECTOR_INV_255;

/**
* Below this weight threshold, animations won't be blended in.
*/
#define ZERO_ANIMWEIGHT_THRESH (0.00001f)

namespace GlobalVectorConstants
{
	static const VectorRegister AnimWeightThreshold = MakeVectorRegister(ZERO_ANIMWEIGHT_THRESH, ZERO_ANIMWEIGHT_THRESH, ZERO_ANIMWEIGHT_THRESH, ZERO_ANIMWEIGHT_THRESH);
	static const VectorRegister RotationSignificantThreshold = MakeVectorRegister(1.0f - DELTA * DELTA, 1.0f - DELTA * DELTA, 1.0f - DELTA * DELTA, 1.0f - DELTA * DELTA);
}

};
