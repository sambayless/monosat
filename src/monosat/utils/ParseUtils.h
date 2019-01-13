/************************************************************************************[ParseUtils.h]
 Copyright (c) 2003-2006, Niklas Een, Niklas Sorensson
 Copyright (c) 2007-2010, Niklas Sorensson

 Permission is hereby granted, free of charge, to any person obtaining a copy of this software and
 associated documentation files (the "Software"), to deal in the Software without restriction,
 including without limitation the rights to use, copy, modify, merge, publish, distribute,
 sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
 furnished to do so, subject to the following conditions:

 The above copyright notice and this permission notice shall be included in all copies or
 substantial portions of the Software.

 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT
 NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
 DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT
 OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 **************************************************************************************************/

#ifndef Minisat_ParseUtils_h
#define Minisat_ParseUtils_h

#include <cstdlib>
#include <cstdio>
#include <string>

#include <zlib.h>
#include "monosat/mtl/Vec.h"
#include <stdexcept>
#include <cstdarg>

namespace Monosat {

//-------------------------------------------------------------------------------------------------
// A simple buffered character stream class:
class parse_error : public std::runtime_error {
public:
    explicit parse_error(const std::string& arg) : std::runtime_error(arg){}
};


//Supporting function for throwing parse errors
inline void parse_errorf(const char* fmt, ...){
    va_list args;
    va_start(args, fmt);
    char buf[1000];
    vsnprintf(buf, sizeof buf, fmt, args);
    va_end(args);
    throw parse_error(buf);
}

static const int buffer_size = 1048576;

class StreamBuffer {
    gzFile in;
    //unsigned char buf[buffer_size];
    vec<unsigned char> buf;//storing the buffer in a vector, to ensure it is heap allocated instead of stack allocated,
    //to avoid overflowing the default JVM stack size (which is very small on some platforms).
    int pos;
    int size;

    void assureLookahead(){
        if(pos >= size){
            pos = 0;
            size = gzread(in, &buf[0], buf.size());
        }
    }

public:
    explicit StreamBuffer(gzFile i) :
            in(i), pos(0), size(0){
        buf.growTo(buffer_size);
        assureLookahead();
    }

    int operator*() const{
        return (pos >= size) ? EOF : buf[pos];
    }

    void operator++(){
        pos++;
        assureLookahead();
    }

    void operator+=(int n){
        assert(n >= 0);
        for(int i = 0; i < n; i++)
            this->operator++();
    }

    int position() const{
        return pos;
    }
};

//-------------------------------------------------------------------------------------------------
// End-of-file detection functions for StreamBuffer and char*:

static inline bool isEof(StreamBuffer& in){
    return *in == EOF;
}

static inline bool isEof(const char* in){
    return *in == '\0';
}

static inline bool isNumber(const char in){
    return ((in >= '0' && in <= '9') || in == '-');
}

static inline bool isWhitespace(const char in){
    return (in >= 9 && in <= 13) || in == 32;
}

//-------------------------------------------------------------------------------------------------
// Generic parse functions parametrized over the input-stream type.

template<class B>
static void skipWhitespace(B& in){
    while(!isEof(in) && ((*in >= 9 && *in <= 13) || *in == 32))
        ++in;
}

template<class B>
static void skipWhitespaceNoNewLines(B& in){
    while(!isEof(in) && ((*in == 9 || *in == 11 || *in == 12) || *in == 32))
        ++in;
}

template<class B>
static void skipLine(B& in){
    for(;;){
        if(isEof(in))
            return;
        if(*in == '\n'){
            ++in;
            return;
        }
        ++in;
    }
}

template<class B>
static int parseInt(B& in){
    int val = 0;
    bool neg = false;
    skipWhitespace(in);
    if(*in == '-')
        neg = true, ++in;
    else if(*in == '+')
        ++in;
    if(*in < '0' || *in > '9')
        parse_errorf("PARSE ERROR! Unexpected char while parsing int: %c\n", *in);
    while(*in >= '0' && *in <= '9')
        val = val * 10 + (*in - '0'), ++in;
    return neg ? -val : val;
}

template<class B>
static int64_t parseLong(B& in){
    int64_t val = 0;
    bool neg = false;
    skipWhitespace(in);
    if(*in == '-')
        neg = true, ++in;
    else if(*in == '+')
        ++in;
    if(*in < '0' || *in > '9')
        parse_errorf("PARSE ERROR! Unexpected char while parsing long: %c\n", *in);
    while(*in >= '0' && *in <= '9')
        val = val * 10 + (*in - '0'), ++in;
    return neg ? -val : val;
}

template<class B>
static double parseDouble(B& in, vec<char>& tmp){
    int val = 0;
    bool neg = false;
    skipWhitespace(in);
    tmp.clear();
    while(*in == '+' || *in == '-' || *in == '.' || *in == 'E' || *in == 'e' || (*in >= '0' && *in <= '9')){
        tmp.push(*in);
        ++in;
    }
    tmp.push(0);

    return strtod(&tmp[0], nullptr);
}

// String matching: in case of a match the input iterator will be advanced the corresponding
// number of characters.
template<class B>
static bool match(B& in, const char* str){
    int i;
    for(i = 0; str[i] != '\0'; i++)
        if(in[i] != str[i])
            return false;

    in += i;

    return true;
}

// String matching: consumes characters eagerly, but does not require random access iterator.
template<class B>
static bool eagerMatch(B& in, const char* str){
    for(; *str != '\0'; ++str, ++in)
        if(*str != *in)
            return false;
    return true;
}

// String matching: in case of a match the input iterator will be advanced the corresponding
// number of characters. The matching much end the input string.
template<class B>
static bool match_end(B& in, const char* str){
    int i;
    for(i = 0; str[i] != '\0'; i++)
        if(in[i] != str[i])
            return false;
    if(in[i] != '\0'){
        return false;
    }
    in += i;

    return true;
}
//=================================================================================================
}

#endif
