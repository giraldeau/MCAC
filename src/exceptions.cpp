/*
 * Copyright (c) 2009-2017, Farooq Mela
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.

 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
 * ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */
//#include <execinfo.h> // for backtrace
//#include <dlfcn.h>    // for dladdr
//#include <cxxabi.h>   // for __cxa_demangle
//
#include <cstdio>
#include <cstdlib>
#include <string>
#include <sstream>
#include "exceptions.hpp"
//
//
//// This function produces a stack backtrace with demangled function & method names.
//std::string Backtrace(int skip)
//{
//    void *callstack[128];
//    const int nMaxFrames = sizeof(callstack) / sizeof(callstack[0]);
//    char buf[1024];
//    int nFrames = backtrace(callstack, nMaxFrames);
//    char **symbols = backtrace_symbols(callstack, nFrames);
//
//    std::ostringstream trace_buf;
//    for (int i = skip; i < nFrames; i++) {
//        printf("%s\n", symbols[i]);
//
//        Dl_info info;
//        if (dladdr(callstack[i], &info) && info.dli_sname) {
//            char *demangled = NULL;
//            int status = -1;
//            if (info.dli_sname[0] == '_')
//                demangled = abi::__cxa_demangle(info.dli_sname, NULL, 0, &status);
//            snprintf(buf, sizeof(buf), "%-3d %*p %s + %zd\n",
//                     i, int(2 + sizeof(void*) * 2), callstack[i],
//                     status == 0 ? demangled :
//                     info.dli_sname == 0 ? symbols[i] : info.dli_sname,
//                     (char *)callstack[i] - (char *)info.dli_saddr);
//            free(demangled);
//        } else {
//            snprintf(buf, sizeof(buf), "%-3d %*p %s\n",
//                     i, int(2 + sizeof(void*) * 2), callstack[i], symbols[i]);
//        }
//        trace_buf << buf;
//    }
//    free(symbols);
//    if (nFrames == nMaxFrames)
//        trace_buf << "[truncated]\n";
//    return trace_buf.str();
//}


#define UNW_LOCAL_ONLY
#include <cxxabi.h>
#include <libunwind.h>
#include <cstdio>
#include <cstdlib>

std::string Backtrace(int skip){
    std::ostringstream trace_buf;
    char buf[1024];

    unw_cursor_t cursor;
    unw_context_t context;

    // Initialize cursor to current frame for local unwinding.
    unw_getcontext(&context);
    unw_init_local(&cursor, &context);

    // Unwind frames one by one, going up the frame stack.
    while (unw_step(&cursor) > 0) {
        if (skip-- > 0) {
            continue;
        }
        unw_word_t offset, pc;
        unw_get_reg(&cursor, UNW_REG_IP, &pc);
        if (pc == 0) {
            break;
        }
        std::sprintf(buf, "0x%lx:", pc);
        trace_buf << std::string(buf);

        char sym[256];
        if (unw_get_proc_name(&cursor, sym, sizeof(sym), &offset) == 0) {
            char* nameptr = sym;
            int status;
            char* demangled = abi::__cxa_demangle(sym, nullptr, nullptr, &status);
            if (status == 0) {
                nameptr = demangled;
            }
            std::sprintf(buf, " (%s+0x%lx)\n", nameptr, offset);
            trace_buf << std::string(buf);
            std::free(demangled);
        } else {
//            std::printf(" -- error: unable to obtain symbol name for this frame\n");
        }
    }
    return trace_buf.str();

}
