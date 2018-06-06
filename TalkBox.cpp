#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "TalkBox.h"
#include "calcAutoCoeff.h"
#include "durbin.h"
#include "lpcFilter.h"

const float k_max = 0.99f;

#define M_PI    3.14159265358979323846

/* a tunable high-pass filter based on a first order allpass */

inline double highpass(double in, double coeff, double *mem)
{
    double out;

    in *= 0.5;

    out = coeff * (in - mem[1]);
    out += mem[0];

    mem[0] = in;                    // non-recursive state
    mem[1] = out;                   // recursive state

    return (in - out);
}

TalkBox::TalkBox(double fs, int num_coeffs, int block_length)
{
    this->fs = fs;
    this->num_coeffs = num_coeffs;
    this->block_length = block_length;
    memory_rms_size = 4;

    input_buffer0 = new float[block_length];
    input_buffer1 = new float[block_length];
    memory_rms = new float[memory_rms_size];

    for (int i = 0; i<num_acf; i++)
        acf[i] = new float[num_coeffs + 1]();

    acf_smooth = new float[num_coeffs + 1];
    a_temp = new float[num_coeffs];
    a_lpc = new float[num_coeffs];
    memory_lpc = new float[num_coeffs];

    // parameter for smoothing of acf
    float acf_tau = 0.03f;  // in seconds
    acf_alpha = 1 - (block_length / ( acf_tau * float(fs)));

    // gate off
    gate_level = 0;

    // high pass design
    double ftan = tan(M_PI * 20000. / fs);
    high_pass_coeff = (int32_t)((ftan - 1) / (ftan + 1) * 0x7FFFFFFF);

    // set states to null
    resetStates();

    sample_buffer = input_buffer0;
    block_buffer  = input_buffer1;

    acf_index = 0;
}

TalkBox::~TalkBox(void)
{
    delete[] input_buffer0;
    delete[] input_buffer1;
    delete[] memory_rms;
    for (int i = 0; i<num_acf; i++)
        delete[] acf[i];
    delete[] acf_smooth;
    delete[] a_temp;
    delete[] a_lpc;
    delete[] memory_lpc;
}

void TalkBox::process(float samples[])
{
   float temp32;

    // synthesizer signal
    temp32 = samples[0];

    // input * gain
    temp32 *= error_gain;

    // input * voice_rms
    temp32 *= voice_rms;

    // all-pole filter
    std::unique_lock<std::mutex> locker(a_coeff_mutex, std::defer_lock);
    locker.lock();

    samples[0] = lpcFilter(temp32, a_lpc, memory_lpc, num_coeffs);

    locker.unlock();

    // voice signal
    sample_buffer[buffer_position++] = samples[1];

    if (buffer_position >= block_length)
    {
        buffer_position = 0;

        if (block_ready == true)
            printf("timing error\n");

        // swap buffer
        float *tmp_ptr = block_buffer;
        block_buffer = sample_buffer;
        sample_buffer = tmp_ptr;

        block_ready = true;
    }
}

void TalkBox::calculateLPCcoefficients(void)
{
    float temp32;
    float abs_voice;
    float error_power;

    // new input block available?
    if (block_ready == false)
        return;

    abs_voice = 0;
    for (int i=0; i<block_length; i++)
    {
        // voice
        temp32 = block_buffer[i];

        abs_voice += fabsf(temp32);

        // high pass
        block_buffer[i] = (float) highpass(temp32, high_pass_coeff, memory_hp);
    }

    abs_voice /= block_length;

    // RMS (FIR)
    for (int i = memory_rms_size - 1; i > 0; i--)
    {
        memory_rms[i] =  memory_rms[i - 1];
    }
    memory_rms[0] = abs_voice;

    voice_rms = 0;
    for (int i = 0; i < memory_rms_size; i++)
    {
        voice_rms += memory_rms[i];
    }
    voice_rms /= memory_rms_size;

    // +12 dB
    voice_rms *= 4;

    if (voice_rms > 1)
        voice_rms = 1;

    if (voice_rms < gate_level)     // gate
    {
        voice_rms = 0;
    }

    calcAutoCoeff(acf[acf_index], num_coeffs+1, block_buffer, block_length);

    // averaging of acfs
    for (int i = 0; i < num_coeffs + 1; i++)
        acf[acf_index][i] = (acf[0][i] + acf[1][i] + acf[2][i] + acf[3][i]) / num_acf;

    for (int i = 0; i < num_coeffs + 1; i++)
    {
        acf_smooth[i] = acf_alpha * acf_smooth[i] + (1.0f - acf_alpha) * acf[acf_index][i];
    }

    if (voice_rms)
    {
        error_power = durbin(acf_smooth, a_temp, num_coeffs, k_max);

        error_gain = sqrtf(error_power);

        std::unique_lock<std::mutex> locker(a_coeff_mutex, std::defer_lock);
        locker.lock();

        for (int i = 0; i < num_coeffs; i++)
            a_lpc[i] = a_temp[i];

        locker.unlock();
    }
    else
    {
        error_gain = 0;
    }

    acf_index++;
    if (acf_index >= num_acf)
        acf_index = 0;

    block_ready = false;
}

void TalkBox::resetStates(void)
{
    voice_rms = 0;
    error_gain = 0;
    buffer_position = 0;
    block_ready = false;

    memory_hp[0] = memory_hp[1] = 0;

    for (int i=0; i<memory_rms_size; i++)
        memory_rms[i] = 0;

    for (int i=0; i<num_coeffs + 1; i++)
        acf_smooth[i] = 0;

    for (int i=0; i<num_coeffs; i++)
        memory_lpc[i] = 0;
}

void TalkBox::checkLpcFilter(float *memory, int num_coeff)
{
    int nan = 0;
    int infinite = 0;

    for (int i=0; i <num_coeff; i++)
    {
        if (_isnan(memory[i]))
            nan = 1;
    }

    for (int i=0; i <num_coeff; i++)
    {
        if (!_finite(memory[i]))
            infinite = 1;
    }

    if (nan || infinite)
    {
        printf("filter is unstable\n");

        for (int i=0; i <num_coeff; i++)
        {
            memory[i] = 0;
        }
    }
}

void TalkBox::setSmoothingTime(float acf_tau)
{
    if (acf_tau > 0)
        acf_alpha = 1 - (block_length / ( acf_tau * float(fs)));
    else
        acf_alpha = 0;

    if (acf_alpha < 0)
        acf_alpha = 0;
}

void TalkBox::setGateLevel(float level)
{
    gate_level = level;
}

void TalkBox::setPreemphasis(float fcuttoff)
{
    double ftan = tan(M_PI * fcuttoff / fs);
    high_pass_coeff = (ftan - 1) / (ftan + 1);
}

int TalkBox::getNumCoeffs(void)
{
    return num_coeffs;
}

void TalkBox::getCoefficients(float all_pole_coefficients[])
{
    for (int i=0; i<num_coeffs; i++)
        all_pole_coefficients[i] = a_lpc[i];
}

float TalkBox::getPreemphasis(void)
{
    return float(high_pass_coeff);
}

float TalkBox::getErrorGain(void)
{
    return ( error_gain );
}

float TalkBox::getVoiceGain(void)
{
    return ( voice_rms );
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
