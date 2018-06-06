#ifndef _TALK_BOX
#define _TALK_BOX

#include <mutex>

const int num_acf = 4;

class TalkBox
{
protected:
    double fs;
    int num_coeffs;
    int block_length;
    int memory_rms_size;
    float voice_rms;
    float error_gain;
    int buffer_position;
    float *input_buffer0;
    float *input_buffer1;
    float *sample_buffer;
    float *block_buffer;
    bool block_ready;
    double high_pass_coeff;
    double memory_hp[2];
    float *memory_rms;
    float acf_alpha;
    float gate_level;
    int16_t acf_index;
    float *acf[num_acf];
    float *acf_smooth;
    float *a_temp;
    float *a_lpc;
    float *memory_lpc;
    std::mutex a_coeff_mutex;

public:
    TalkBox(double fs, int num_coeffs, int block_length);
    ~TalkBox(void);
    void process(float samples[]);
    void calculateLPCcoefficients(void);
    void resetStates(void);
    void checkLpcFilter(float *memory, int num_coeff);
    void setSmoothingTime(float acf_tau);
    void setGateLevel(float level);
    void setPreemphasis(float fcuttoff);
    int  getNumCoeffs(void);
    void getCoefficients(float all_pole_coefficients[]);
    float getPreemphasis(void);
    float getErrorGain(void);
    float getVoiceGain(void);
};

#endif  // _TALK_BOX

