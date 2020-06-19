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
#include "constants.hpp"


namespace mcac {
class BaseException : public std::exception {
public:
    std::string name{"Unkown Error"};
    std::string detail{""};
    std::string message{""};
    ErrorCodes code{ErrorCodes::UNKNOWN_ERROR};
    void write_message() {
        if (detail.empty()) {
            message = name;
        } else {
            message = name + ": " + detail;
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
        char *cstr = new char[message.length() + 1];
        strcpy(cstr, message.c_str());
        return cstr;
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
}  //namespace mcac
#endif //INCLUDE_EXCEPTIONS_HPP
