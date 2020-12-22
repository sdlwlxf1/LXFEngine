// this tiny game engine is rendered by CPU. I make it for learing the pipline and rasterization. This is the first step to realise my own game engine. so just do it!

#pragma once
#include <string.h>

namespace LXF {

#define FORCEINLINE __forceinline
#define PLATFORM_ENABLE_VECTORINTRINSICS					1

template<typename T32BITS, typename T64BITS, int PointerSize>
struct SelectIntPointerType
{
	// nothing here are is it an error if the partial specializations fail
};

template<typename T32BITS, typename T64BITS>
struct SelectIntPointerType<T32BITS, T64BITS, 8>
{
	typedef T64BITS TIntPointer; // select the 64 bit type
};

template<typename T32BITS, typename T64BITS>
struct SelectIntPointerType<T32BITS, T64BITS, 4>
{
	typedef T32BITS TIntPointer; // select the 32 bit type
};


typedef signed int int32;
typedef unsigned int uint32;
typedef signed long long int64;
typedef unsigned long long uint64;
typedef signed char int8;
typedef unsigned char uint8;
typedef unsigned short int	uint16;		// 16-bit unsigned.
typedef signed short int	int16;		// 16-bit signed.

typedef SelectIntPointerType<uint32, uint64, sizeof(void*)>::TIntPointer UPTRINT;	// unsigned int the same size as a pointer
typedef SelectIntPointerType<int32, int64, sizeof(void*)>::TIntPointer PTRINT;		// signed int the same size as a pointer
typedef UPTRINT SIZE_T;																// unsigned int the same size as a pointer
typedef PTRINT SSIZE_T;																// signed int the same size as a pointer

#ifndef RESTRICT
#define RESTRICT __restrict						/* no alias hint */
#endif

#if DO_GUARD_SLOW
#define checkSlow(expr)					check(expr)
#define checkfSlow(expr, format, ...)	checkf(expr, format, ##__VA_ARGS__)
#define verifySlow(expr)				check(expr)
#else
#define checkSlow(expr)					{  }
#define checkfSlow(expr, format, ...)	{  }
#define verifySlow(expr)				{  }
#endif

#define ensure(           InExpression                ) ()
#define ensureMsgf(       InExpression, InFormat, ... ) ()
#define ensureAlways(     InExpression                ) ()
#define ensureAlwaysMsgf( InExpression, InFormat, ... ) ()

#ifndef checkCode
#define checkCode( Code )		do { Code; } while ( false );
#endif
#ifndef verify
#define verify(expr)		{}	
#endif
#ifndef check
#define check(expr)			{}	
#endif


class Memory
{
public:
	static void* Memzero(void* Dest, SIZE_T Count) { return memset(Dest, 0, Count); }
	static void* Memset(void* Dest, uint8 Char, SIZE_T Count) { return memset(Dest, Char, Count); }
	static void* Memcpy(void* Dest, const void* Src, SIZE_T Count) { return memcpy(Dest, Src, Count); }
	static void* Memmove(void* Dest, const void* Src, SIZE_T Count) { return memmove(Dest, Src, Count); }
	static int32 Memcmp(const void* Buf1, const void* Buf2, SIZE_T Count) { return memcmp(Buf1, Buf2, Count); }

};
};
