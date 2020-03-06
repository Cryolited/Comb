#ifndef ANALYSISBANK_H
#define ANALYSISBANK_H

#include <iostream>
#include <vector>
#include <cmath>
#include <complex>
#include <algorithm>
#include <fftw3.h>


using namespace std;

struct params
{
    uint16_t WIN_H_RADIX = 18;
    uint16_t FB_OVERLAP_RATIO = 2; // перекрытие
    uint16_t NFFT = 128; // колво фильтров
    uint16_t WIN_OVERLAP_RATIO = 8; // длина фильтров
};

class AnalysisBank
{
public:


    class signal {
    public:
        vector<int16_t> si;
        vector<int16_t> sq;
        void resizeS(uint16_t length) // с нулями
        {
           si.resize(length);
           sq.resize(length);
        }
    };

    AnalysisBank();
    //void setFs(uint32_t Fs);
    //float getFs();
    void openSignal(uint32_t size);
    void readSignal(vector<int16_t>& inVecQ, vector<int16_t>& inVecI);
    void readSignal(int16_t inVecQ[], int16_t inVecI[], uint32_t size);
    void saveSignal();
    void setParams(params pFun);
    void createNpr();
    void filterAnalyze();
    void saveAnalyzeFB();
    vector<vector<complex<double>>>& getAnalyzeFB();
    void saveCreateFB();

    signal sig;
    vector<int32_t> h_fb_win_fxp;
    vector<vector<complex<double>> > filt;
    params p;
private:

    //uint16_t WIN_H_RADIX = 18;
   // uint16_t FB_OVERLAP_RATIO = 2; // перекрытие
   // uint16_t NFFT = 128; // колво фильтров
   // uint16_t WIN_OVERLAP_RATIO = 8; // длина фильтров
    //uint32_t Fs = 512e6;
    //double t_us = 0.2 * 1e-6;
    //uint16_t first_period_part = p.NFFT ;
    //uint16_t second_period_part = p.NFFT * WIN_OVERLAP_RATIO;
    vector<complex<double>> pulse_sig_phase_n;
    int32_t maxSummLog;
    vector<complex<double>> filtered;
    vector<complex<double>> sigOut;
    uint16_t coeff_radix;

    void npr_coeff(int16_t N,int16_t L);
    void npr_synthesis();
    void non_maximally_decimated_fb();
    void maximally_decimated_fb(int16_t ovRat);
    int32_t fbWinMaxGain();
    void fft( fftw_complex* in, fftw_complex* out, uint32_t fft_size, bool backward );
    complex<double> fpga_round( complex<double> din, int32_t shval );
    void circshift( complex<double>* x,unsigned sizex, unsigned pos );

};

#endif // ANALYSISBANK_H
