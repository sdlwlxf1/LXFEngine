// Copyright 1998-2019 Epic Games, Inc. All Rights Reserved.
#pragma once

namespace LXF {

#if USING_CODE_ANALYSIS
#if !defined( CA_IN ) || !defined( CA_OUT ) || !defined( CA_READ_ONLY ) || !defined( CA_WRITE_ONLY ) || !defined( CA_VALID_POINTER ) || !defined( CA_CHECK_RETVAL ) || !defined( CA_NO_RETURN ) || !defined( CA_SUPPRESS ) || !defined( CA_ASSUME )
#error Code analysis macros are not configured correctly for this platform
#endif
#else
	// Just to be safe, define all of the code analysis macros to empty macros
#define CA_IN 
#define CA_OUT
#define CA_READ_ONLY
#define CA_WRITE_ONLY
#define CA_VALID_POINTER
#define CA_CHECK_RETVAL
#define CA_NO_RETURN
#define CA_SUPPRESS( WarningNumber )
#define CA_ASSUME( Expr )
#define CA_CONSTANT_IF(Condition) if (Condition)
#endif

#ifndef USING_THREAD_SANITISER
#define USING_THREAD_SANITISER 0
#endif

#if USING_THREAD_SANITISER
#if !defined( TSAN_SAFE ) || !defined( TSAN_BEFORE ) || !defined( TSAN_AFTER ) || !defined( TSAN_ATOMIC )
#error Thread Sanitiser macros are not configured correctly for this platform
#endif
#else
	// Define TSAN macros to empty when not enabled
#define TSAN_SAFE
#define TSAN_BEFORE(Addr)
#define TSAN_AFTER(Addr)
#define TSAN_ATOMIC(Type) Type
#endif

enum { INDEX_NONE = -1 };
enum { UNICODE_BOM = 0xfeff };

enum EForceInit
{
	ForceInit,
	ForceInitToZero
};
enum ENoInit { NoInit };
enum EInPlace { InPlace };

// Handle type to stably track users on a specific platform
typedef int32 FPlatformUserId;
const FPlatformUserId PLATFORMUSERID_NONE = INDEX_NONE;

// Push and pop macro definitions
#ifdef __clang__
#define PUSH_MACRO(name) _Pragma(PREPROCESSOR_TO_STRING(push_macro(PREPROCESSOR_TO_STRING(name))))
#define POP_MACRO(name) _Pragma(PREPROCESSOR_TO_STRING(pop_macro(PREPROCESSOR_TO_STRING(name))))
#else
#define PUSH_MACRO(name) __pragma(push_macro(PREPROCESSOR_TO_STRING(name)))
#define POP_MACRO(name) __pragma(pop_macro(PREPROCESSOR_TO_STRING(name)))
#endif


#ifdef __COUNTER__
	// Created a variable with a unique name
#define ANONYMOUS_VARIABLE( Name ) PREPROCESSOR_JOIN(Name, __COUNTER__)
#else
	// Created a variable with a unique name.
	// Less reliable than the __COUNTER__ version.
#define ANONYMOUS_VARIABLE( Name ) PREPROCESSOR_JOIN(Name, __LINE__)
#endif

};
