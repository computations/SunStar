//debug.h
//Ben Bettisworth
//Macros to assist in debugging. When using this, and the DEBUG preprocessor 
//flag is defined, then verbosity is increased extremely. Its useful to put
//this here, so that I can't accidently add functionality into a debug section.
//Secondly, when DEBUG_IF_FLAG is defined to be zero, and an optmizer is turned
//on, then that code gets pruned.
//
//I really don't think this will ever be pretty

#pragma once

#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <execinfo.h>
#include <time.h>

const clock_t CLOCK_START = clock();

#ifdef DEBUG
#define DEBUG_IF_FLAG 1
#else
#define DEBUG_IF_FLAG 0
#endif

#ifdef EMIT_DEBUG
#define EMIT_DEBUG_FLAG 1
#else
#define EMIT_DEBUG_FLAG 0
#endif

#define debug_print(fmt, ...) {if(DEBUG_IF_FLAG && EMIT_DEBUG_FLAG) fprintf(stderr, "[%f] %s:%d:%s(): " fmt "\n",\
        ((double)clock() - CLOCK_START)/CLOCKS_PER_SEC, __FILE__,\
        __LINE__, __func__, __VA_ARGS__); }

#define debug_string(x) { debug_print("%s", x)}

#define debug_matrix(desc,x,s) { if(DEBUG_IF_FLAG && EMIT_DEBUG_FLAG){\
    debug_string(desc);\
    for(size_t i=0;i<s;++i){fprintf(stderr, "[%f]\t", ((double)clock() - CLOCK_START)/CLOCKS_PER_SEC); for(size_t j=0;j<s;++j){fprintf(stderr, "%.1f\t", x[i*s+j]);} fprintf(stderr, "\n");}}}

#define debug_print_map(desc, map){ if(DEBUG_IF_FLAG && EMIT_DEBUG_FLAG){\
    debug_string(desc);for(auto kv:map){fprintf(stderr, "(%s : %lu) ", kv.first.c_str(), kv.second);} fprintf\
}}

#define print_trace(){ if(DEBUG_IF_FLAG) {\
        void* callstack[128];\
        int frames = backtrace(callstack, 128);\
        char** bt_symbols = backtrace_symbols(callstack, frames);\
        fprintf(stderr, "BACKTRACE AT %s:%d:%s():\n", __FILE__, __LINE__, __func__);\
        for(int i=0;i<frames;++i){\
            fprintf(stderr, "%s\n", bt_symbols[i]);\
        }\
    }\
}

#define assert_string(cond, comment) if(DEBUG_IF_FLAG){ {\
    if(!(cond)){\
        fprintf(stderr, "assertion \"%s\" failed: file: %s, line: %d, comment: %s\n",#cond, __FILE__, __LINE__, comment);\
        abort();\
    }\
} }
