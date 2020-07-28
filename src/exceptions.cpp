/*
 * MCAC
 * Copyright (C) 2020 CORIA
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
#include <cstdio>
#include <cstdlib>
#include <string>
#include "exceptions.hpp"


#define UNW_LOCAL_ONLY
#include <cxxabi.h>
#include <libunwind.h>


std::string Backtrace(int skip) {
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
            char *nameptr = sym;
            int status;
            char *demangled = abi::__cxa_demangle(sym, nullptr, nullptr, &status);
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
