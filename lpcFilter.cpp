#include "lpcFilter.h"

float lpcFilter(float inputSample, float *a, float *memory, int num_coeff)
{
    float output;

    output = 0;
    for (int i=0; i <num_coeff; i++)
    {
        output += a[i] * memory[i];
    }
    output = inputSample - output;

    for (int i=num_coeff-1; i>0; i--)
    {
        memory[i] = memory[i-1];
    }

    memory[0] = output;

    return output;
}

//--------------------- License ------------------------------------------------

// Copyright (c) 2016 Finn Bayer, Christoph Eike, Uwe Simmer

// Permission is hereby granted, free of charge, to any person obtaining 
// a copy of this software and associated documentation files 
// (the "Software"), to deal in the Software without restriction, 
// including without limitation the rights to use, copy, modify, merge, 
// publish, distribute, sublicense, and/or sell copies of the Software, 
// and to permit persons to whom the Software is furnished to do so, 
// subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included 
// in all copies or substantial portions of the Software.

// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, 
// EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF 
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. 
// IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY 
// CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, 
// TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE 
// SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

//------------------------------------------------------------------------------