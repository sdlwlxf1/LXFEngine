// Copyright 1998-2019 Epic Games, Inc. All Rights Reserved.

#pragma once

namespace LXF {

#ifndef BUILD_DEBUG
#define BUILD_DEBUG				0
#endif
#ifndef BUILD_DEVELOPMENT
#define BUILD_DEVELOPMENT		0
#endif
#ifndef BUILD_TEST
#define BUILD_TEST				0
#endif
#ifndef BUILD_SHIPPING
#define BUILD_SHIPPING			0
#endif
#ifndef GAME
#define GAME						0
#endif
#ifndef EDITOR
#define EDITOR					0
#endif
#ifndef BUILD_SHIPPING_WITH_EDITOR
#define BUILD_SHIPPING_WITH_EDITOR 0
#endif
#ifndef BUILD_DOCS
#define BUILD_DOCS				0
#endif

/*--------------------------------------------------------------------------------
	Basic options that by default depend on the build configuration and platform

	DO_GUARD_SLOW									If true, then checkSlow, checkfSlow and verifySlow are compiled into the executable.
	DO_CHECK										If true, then checkCode, checkf, verify, check, checkNoEntry, checkNoReentry, checkNoRecursion, verifyf, checkf, ensure, ensureAlways, ensureMsgf and ensureAlwaysMsgf are compiled into the executables
	STATS											If true, then the stats system is compiled into the executable.
	ALLOW_DEBUG_FILES								If true, then debug files like screen shots and profiles can be saved from the executable.
	NO_LOGGING										If true, then no logs or text output will be produced

--------------------------------------------------------------------------------*/

#if BUILD_DEBUG
#ifndef DO_GUARD_SLOW
#define DO_GUARD_SLOW									1
#endif
#ifndef DO_CHECK
#define DO_CHECK										1
#endif
#ifndef STATS
#define STATS											((WITH_UNREAL_DEVELOPER_TOOLS || !WITH_EDITORONLY_DATA || USE_STATS_WITHOUT_ENGINE || USE_MALLOC_PROFILER || FORCE_USE_STATS) && !ENABLE_STATNAMEDEVENTS)
#endif
#ifndef ALLOW_DEBUG_FILES
#define ALLOW_DEBUG_FILES								1
#endif
#ifndef ALLOW_CONSOLE
#define ALLOW_CONSOLE									1
#endif
#ifndef NO_LOGGING
#define NO_LOGGING										0
#endif
#elif BUILD_DEVELOPMENT
#ifndef DO_GUARD_SLOW
#define DO_GUARD_SLOW									0
#endif
#ifndef DO_CHECK
#define DO_CHECK										1
#endif
#ifndef STATS
#define STATS											((WITH_UNREAL_DEVELOPER_TOOLS || !WITH_EDITORONLY_DATA || USE_STATS_WITHOUT_ENGINE || USE_MALLOC_PROFILER || FORCE_USE_STATS) && !ENABLE_STATNAMEDEVENTS)
#endif
#ifndef ALLOW_DEBUG_FILES
#define ALLOW_DEBUG_FILES								1
#endif
#ifndef ALLOW_CONSOLE
#define ALLOW_CONSOLE									1
#endif
#ifndef NO_LOGGING
#define NO_LOGGING										0
#endif
#elif BUILD_TEST
#ifndef DO_GUARD_SLOW
#define DO_GUARD_SLOW									0
#endif
#ifndef DO_CHECK
#define DO_CHECK										USE_CHECKS_IN_SHIPPING
#endif
#ifndef STATS
#define STATS											((USE_MALLOC_PROFILER || FORCE_USE_STATS) && !ENABLE_STATNAMEDEVENTS)
#endif
#ifndef ALLOW_DEBUG_FILES
#define ALLOW_DEBUG_FILES								1
#endif
#ifndef ALLOW_CONSOLE
#define ALLOW_CONSOLE									1
#endif
#ifndef NO_LOGGING
#define NO_LOGGING										!USE_LOGGING_IN_SHIPPING
#endif
#elif BUILD_SHIPPING
#if WITH_EDITOR
#ifndef DO_GUARD_SLOW
#define DO_GUARD_SLOW								0
#endif
#ifndef DO_CHECK
#define DO_CHECK									1
#endif
#ifndef STATS
#define STATS										1
#endif
#ifndef ALLOW_DEBUG_FILES
#define ALLOW_DEBUG_FILES							1
#endif
#ifndef ALLOW_CONSOLE
#define ALLOW_CONSOLE								0
#endif
#ifndef NO_LOGGING
#define NO_LOGGING									0
#endif
#else
#ifndef DO_GUARD_SLOW
#define DO_GUARD_SLOW								0
#endif
#ifndef DO_CHECK
#define DO_CHECK									USE_CHECKS_IN_SHIPPING
#endif
#ifndef STATS
#define STATS										(FORCE_USE_STATS && !ENABLE_STATNAMEDEVENTS)
#endif
#ifndef ALLOW_DEBUG_FILES
#define ALLOW_DEBUG_FILES							0
#endif
#ifndef ALLOW_CONSOLE
#define ALLOW_CONSOLE								ALLOW_CONSOLE_IN_SHIPPING
#endif
#ifndef NO_LOGGING
#define NO_LOGGING									!USE_LOGGING_IN_SHIPPING
#endif
#endif
#else
#endif
};
