#include <float.h>
#include "calcAutoCoeff.h"

void calcAutoCoeff(float *acf, int num_acf, float *signal, int num_signal)
{
    for (int k = 0; k < num_acf; k++)
    {
        acf[k] = 0;
    }

    for (int k = 0; k <num_acf; k++)
    {
        for (int i = 0; i < num_signal - k; i++)
        {
            acf[k] += (signal[i + k] * signal[i]);
        }
    }

    if (acf[0] < FLT_MIN)
    {
        acf[0] = 1;

        for (int k = 1; k < num_acf; k++)
        {
            acf[k] = 0;
        }
        return;
    }

    float acf0 = acf[0];
    for (int i = 0; i <num_acf; i++)
    {
        acf[i] = acf[i] / acf0;
    }
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