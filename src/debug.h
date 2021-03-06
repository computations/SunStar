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
#include <sys/ioctl.h>
#include <unistd.h>

const clock_t CLOCK_START = clock();
extern bool __PROGRESS_BAR_FLAG__;

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

#define print_clock {if(DEBUG_IF_FLAG && EMIT_DEBUG_FLAG){fprintf(stderr, "[%f] ", \
        ((double)clock() - CLOCK_START)/CLOCKS_PER_SEC);}};

#define debug_print(fmt, ...) {if(DEBUG_IF_FLAG && EMIT_DEBUG_FLAG){print_clock; fprintf(stderr, "%s:%d:%s(): " fmt "\n",\
         __FILE__, __LINE__, __func__, __VA_ARGS__); }}

#define debug_string(x) { debug_print("%s", x)}

#define debug_vector(desc,x){ if(DEBUG_IF_FLAG && EMIT_DEBUG_FLAG){\
    debug_string(desc);\
    print_clock; for(size_t i=0;i<x.size();++i){fprintf(stderr,"%.1f\t",x[i]);}fprintf(stderr,"\n");}}

#define debug_matrix(desc,x,s) { if(DEBUG_IF_FLAG && EMIT_DEBUG_FLAG){\
    debug_string(desc);\
    for(size_t i=0;i<s;++i){fprintf(stderr, "[%f]\t", ((double)clock() - CLOCK_START)/CLOCKS_PER_SEC); for(size_t j=0;j<s;++j){fprintf(stderr, "%.1f\t", x[i*s+j]);} fprintf(stderr, "\n");}}}

#define debug_d2vector_t(desc,m) { if(DEBUG_IF_FLAG && EMIT_DEBUG_FLAG){\
    debug_string(desc);\
    for(size_t i=0;i<m.size();++i){print_clock; for(size_t j=0;j<m.size();++j){fprintf(stderr, "%.1f\t", m[i][j]);} fprintf(stderr, "\n");}}}

#define debug_print_map(desc, map){ if(DEBUG_IF_FLAG && EMIT_DEBUG_FLAG){\
    print_clock; fprintf(stderr, desc);for(auto kv:map){fprintf(stderr, "(%s : %lu) ", kv.first.c_str(), kv.second);} fprintf(stderr, "\n");\
}}

#define print_trace(){ if(DEBUG_IF_FLAG) {\
        void* callstack[128];\
        int frames = backtrace(callstack, 128);\
        char** bt_symbols = backtrace_symbols(callstack, frames);\
        print_clock;\
        fprintf(stderr, "BACKTRACE AT %s:%d:%s():\n", __FILE__, __LINE__, __func__);\
        for(int i=0;i<frames;++i){\
            print_clock;\
            fprintf(stderr, "%s\n", bt_symbols[i]);\
        }\
    }\
}

#define assert_string(cond, comment) if(DEBUG_IF_FLAG&&EMIT_DEBUG_FLAG){ {\
    if(!(cond)){\
        print_clock;\
        fprintf(stderr, "assertion \"%s\" failed: file: %s, line: %d, comment: %s\n",#cond, __FILE__, __LINE__, comment);\
        abort();\
    }\
} }

#define turn_on_progress(){ __PROGRESS_BAR_FLAG__ = true; }

#define turn_off_progress(){ __PROGRESS_BAR_FLAG__ = false; } 

#define toggle_progress(){ __PROGRESS_BAR_FLAG__ = !__PROGRESS_BAR_FLAG__; }


#define print_progress(done, total) { if(__PROGRESS_BAR_FLAG__){ struct winsize w; ioctl(STDOUT_FILENO, TIOCGWINSZ, &w);\
    print_progress_cols(done, total, w.ws_col); }}

#define print_progress_cols(done, total, cols) { if(__PROGRESS_BAR_FLAG__){\
    size_t adj_cols = cols - 5; if(done == 0) adj_cols--;\
    size_t total_padding=0; size_t tmp = total; while(tmp!=0){ tmp/=10; adj_cols--; total_padding++;}\
    size_t done_padding=total_padding;tmp = done; while(tmp!=0){ tmp/=10; adj_cols--; done_padding--;} \
    if(done==0) done_padding=total_padding-1;\
    adj_cols-=done_padding;\
    size_t done_cols = (done*adj_cols/total); size_t left_cols = adj_cols - done_cols; printf("\r");\
    fprintf(stdout,"[");\
    fprintf(stdout,"\e[32m");for(size_t i = 0; i<done_cols; ++i){ fprintf(stdout,"="); }\
    fprintf(stdout,"\e[31m");for(size_t i = 0; i<left_cols; ++i){ fprintf(stdout,"-"); }\
    fprintf(stdout, "\e[0m][");\
    for(size_t _i_ =0; _i_<done_padding; ++_i_){ fprintf(stdout," "); }\
    fprintf(stdout,"\e[33m%lu\e[0m/\e[33m%lu\e[0m]", done, total); fprintf(stdout, "\e[0m");fflush(stdout);\
}}

#define finish_progress(){if(__PROGRESS_BAR_FLAG__){fprintf(stdout, "\n"); fflush(stdout);}}
