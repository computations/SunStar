#pragma once

#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <execinfo.h>

#ifdef DEBUG
#define DEBUG 1
#else
#define DEBUG 0
#endif

#define debug_print(fmt, ...) {if(DEBUG) fprintf(stderr, "%s:%d:%s(): " fmt "\n", __FILE__,\
        __LINE__, __func__, __VA_ARGS__); }
#define debug_string(x) { debug_print("%s", x)}

#define print_trace(){ if(DEBUG) {\
        void* callstack[128];\
        int frames = backtrace(callstack, 128);\
        char** bt_symbols = backtrace_symbols(callstack, frames);\
        fprintf(stderr, "BACKTRACE AT %s:%d:%s():\n", __FILE__, __LINE__, __func__);\
        for(int i=0;i<frames;++i){\
            fprintf(stderr, "%s\n", bt_symbols[i]);\
        }\
    }\
}

#define assert_string(cond, comment) {\
    if(!(cond)){\
        fprintf(stderr, "assertion \"%s\" failed: file: %s, line: %d, comment: %s\n",#cond, __FILE__, __LINE__, comment);\
        abort();\
    }\
}
