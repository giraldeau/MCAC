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
#ifndef INCLUDE_EXCEPTIONS_HPP
#define INCLUDE_EXCEPTIONS_HPP
#include <exception>
#include <string>
#include <sstream>
#include <iostream>
#include <cstring>
#include "constants.hpp"

std::string Backtrace(int skip = 2);

namespace mcac {
class BaseException : public std::exception {
public:
    std::string name{"Unkown Error"};
    std::string detail{""};
    std::string message{""};
    ErrorCodes code{ErrorCodes::UNKNOWN_ERROR};
    void write_message() {
        message = Backtrace() + "\n" + name;
        if (!detail.empty()) {
            message += ": " + detail;
        }
    }
    BaseException() {
        write_message();
    }
    explicit BaseException(const std::string &new_detail) : BaseException() {
        detail = new_detail;
        write_message();
    }
    [[nodiscard]] const char *what() const noexcept override {
        return message.c_str();
    }
};

class IOError : public BaseException {
public:
    IOError() {
        name = "I/O Error";
        code = ErrorCodes::IO_ERROR;
        write_message();
    }
    explicit IOError(const std::string &new_detail) : IOError() {
        detail = new_detail;
        write_message();
    }
};

class VerletError : public BaseException {
public:
    VerletError() {
        name = "Verlet Error";
        code = ErrorCodes::VERLET_ERROR;
        write_message();
    }
    explicit VerletError(const std::string &new_detail) : VerletError() {
        detail = new_detail;
        write_message();
    }
};

class InputError : public BaseException {
public:
    InputError() {
        name = "Input Error";
        code = ErrorCodes::INPUT_ERROR;
        write_message();
    }
    explicit InputError(const std::string &new_detail) : InputError() {
        detail = new_detail;
        write_message();
    }
};

class AbandonError : public BaseException {
public:
    AbandonError() {
        name = "Abandon";
        code = ErrorCodes::ABANDON_ERROR;
        write_message();
    }
    explicit AbandonError(const std::string &new_detail) : AbandonError() {
        detail = new_detail;
        write_message();
    }
};

class TooDenseError : public BaseException {
public:
    TooDenseError() {
        name = "Too dense error";
        code = ErrorCodes::TOO_DENSE_ERROR;
        std::stringstream standard_error;
        standard_error << "Impossible to add new aggregate without overlap." << std::endl;
        detail = standard_error.str();
        write_message();
    }
    explicit TooDenseError(const std::string &new_detail) : TooDenseError() {
        detail = new_detail;
        write_message();
    }
};

class VolSurfError : public BaseException {
public:
    VolSurfError() {
        name = "Volume/Surface error";
        code = ErrorCodes::VOL_SURF_ERROR;
        std::stringstream standard_error;
        standard_error << "An Aggregate has a negative volume or surface." << std::endl;
        detail = standard_error.str();
        write_message();
    }
    explicit VolSurfError(const std::string &new_detail) : VolSurfError() {
        detail = new_detail;
        write_message();
    }
};

class MergeError : public BaseException {
public:
    MergeError() {
        name = "Merge error";
        code = ErrorCodes::MERGE_ERROR;
        write_message();
    }
    explicit MergeError(const std::string &new_detail) : MergeError() {
        detail = new_detail;
        write_message();
    }
};
}  //namespace mcac
#endif //INCLUDE_EXCEPTIONS_HPP
