
#include "analysisbank.h"

class Generator {
    public:

    struct params
    {
        uint16_t WIN_H_RADIX = 18;
        uint16_t FB_OVERLAP_RATIO = 2; // перекрытие
        uint16_t NFFT = 128; // колво фильтров
        uint16_t WIN_OVERLAP_RATIO = 8; // длина фильтров
    };
    class signal {
    public:
        vector<float> t;
        vector<int16_t> i;
        vector<int16_t> q;
        vector<int16_t> si;
        vector<int16_t> sq;
        void resize(uint16_t length) // без нулей
        {
           t.resize(length);
           i.resize(length);
           q.resize(length);
        }
        void resizeS(uint16_t length) // с нулями
        {
           si.resize(length);
           sq.resize(length);
        }
    };

    void genSignal(float LFM) // генерирование сигнала
    {
        params p;
        first_period_part = p.NFFT ;
        second_period_part = p.NFFT * p.WIN_OVERLAP_RATIO;
        uint16_t DATA_AMPL           = 30000;
        uint16_t THRESHOLD           = DATA_AMPL/2;
        double time_sec            = t_us;
        uint16_t f0                  = 0;
        uint16_t len = Fs*time_sec ;
        uint16_t len2 = ((len+ first_period_part + second_period_part)/p.NFFT)*p.NFFT ;
        sig.resize(len);
        sig.resizeS(len2);
        float F0 = -LFM/2;
        uint32_t l =0;
        float b = LFM/time_sec;
        for (int n = 0; n < len ; ++n) {
            sig.t[n] = (float)n/Fs;
            if ( LFM == 0)
            {
                sig.q[n] = round(DATA_AMPL*sin(2.0*M_PI*1e6*sig.t[n])) +rand()*1000 ;
                sig.i[n] = round(DATA_AMPL*cos(2.0*M_PI*1e6*sig.t[n])) +rand()*1000;
            }else {
                sig.q[n] = round(DATA_AMPL*sin(2.0*M_PI*(F0*sig.t[n]+b/2.0*sig.t[n]*sig.t[n]))) ;
                sig.i[n] = round(DATA_AMPL*cos(2.0*M_PI*(F0*sig.t[n]+b/2.0*sig.t[n]*sig.t[n]))) ;
            }

        }
        for (int n = 0; n < len2; ++n) {
            if (n < first_period_part)
            {
                sig.sq[n] = 0;
                sig.si[n] = 0;
            }
            else
            if (n < len + first_period_part )
            {
                sig.sq[n] = sig.q[l] ;
                sig.si[n] = sig.i[l] ;
                l++;
            }
            else
            {
                sig.sq[n] = 0;
                sig.si[n] = 0;
            }
        }

    }
    signal sig;
    double t_us = 0.2 * 1e-6;
    uint16_t first_period_part;
    uint16_t second_period_part;
    uint32_t Fs = 512e6;
};
