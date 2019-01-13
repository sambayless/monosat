/***************************************************************************************[Global.cc]
Copyright (c) 2005-2010, Niklas Een, Niklas Sorensson

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

#include <cassert>
#include <cstring>
#include <cstdio>
#include "Global.h"

namespace Monosat {
namespace PB {
char* vnsprintf(const char* format, va_list args){
    static FILE* dummy = fopen("/dev/null", "wb");
    unsigned chars_written;
    char* ret;
    va_list args_copy;

#ifdef __va_copy
    __va_copy (args_copy, args);
#else
    args_copy = args;
#endif
    chars_written = vfprintf(dummy, format, args);
    ret = xmalloc<char>(chars_written + 1);
    ret[chars_written] = (char) 255;
    args = args_copy;
    vsprintf(ret, format, args);
    assert(ret[chars_written] == 0);
    return ret;
}


char* nsprintf(const char* format, ...){
    va_list args;
    va_start(args, format);
    char* ret = vnsprintf(format, args);
    va_end(args);
    return ret;
}

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

macro bool hasChar(cchar* text, int chr){
    while(*text != 0) if(*text++ == chr) return true;
    return false;
}

void splitString(cchar* text, cchar* seps, vec<char*>& out){
    while(hasChar(seps, *text)) text++;
    if(*text == 0) return;
    cchar* start = text;
    for(;;){
        if(*text == 0 || hasChar(seps, *text)){
            out.push(xstrndup(start, text - start));
            while(hasChar(seps, *text)) text++;
            if(*text == 0) return;
            start = text;
        }else
            text++;
    }
}


//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
bool opt_satlive = true;
bool opt_ansi = true;

void reportf(const char* format, ...){
    static bool col0 = true;
    static bool bold = false;
    va_list args;
    va_start(args, format);
    char* text = vnsprintf(format, args);
    va_end(args);

    for(char* p = text; *p != 0; p++){
        if(col0 && opt_satlive)
            putchar('c'), putchar(' ');

        if(*p == '\b'){
            bold = !bold;
            if(opt_ansi)
                putchar(27), putchar('['), putchar(bold ? '1' : '0'), putchar('m');
            col0 = false;
        }else{
            putchar(*p);
            col0 = (*p == '\n' || *p == '\r');
        }
    }
    xfree(text);
    fflush(stdout);
}
}
}