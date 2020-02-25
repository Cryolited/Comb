#include <iostream>
#include <vector>
#include <cmath>
#include <complex> // mb not need
#include <fftw3.h>

#define WIN_H_RADIX        18
#define FB_OVERLAP_RATIO  2
#define NFFT  128
#define WIN_OVERLAP_RATIO  8


using namespace std;

class AnalysisBank {
    public:

        void setFs(uint32_t Fs)
        {
            this->Fs = Fs;
        }

        float getFs()
        {
            return Fs;
        }

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

        void genSignal(float LFM)
        {
            uint16_t DATA_AMPL           = 30000;
            uint16_t THRESHOLD           = DATA_AMPL/2;
            double time_sec            = t_us;
            uint16_t f0                  = 0;
            uint16_t len = Fs*time_sec ;
            uint16_t len2 = ((len+ first_period_part + second_period_part)/NFFT)*NFFT ;
            sig.resize(len);
            sig.resizeS(len2);
            float F0 = -LFM/2;
            uint32_t l =0;
            float b = LFM/time_sec;
            for (int n = 0; n < len ; ++n) {
                sig.t[n] = (float)n/Fs;
                if ( LFM == 0)
                {
                    sig.q[n] = round(DATA_AMPL*sin(2.0*M_PI*1e6*sig.t[n])) ;
                    sig.i[n] = round(DATA_AMPL*cos(2.0*M_PI*1e6*sig.t[n])) ;
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

        void createNpr()
        {
            npr_coeff(NFFT,2*WIN_OVERLAP_RATIO);
            get_signal();
        }



        vector<int32_t> h_fb_win_fxp;
        signal sig;
    private:
        uint32_t Fs = 512e6;


        uint32_t BAND              = 500e6;
        uint32_t FW                = 4e6;
        double t_us = 0.2 * 1e-6;
        double period_us = 10 * 1e-6;
        uint16_t first_period_part = NFFT ;
        uint16_t second_period_part = NFFT * WIN_OVERLAP_RATIO;
        vector<complex<double>> pulse_sig_phase_n;
        int32_t maxSummLog;
        vector<complex<double>> filtered;



        uint16_t coeff_radix;

        void npr_coeff(int16_t N,int16_t L)
        {
            //if (L == 16)
              double K=5.856;
            int16_t M = N /2;            
            vector<complex<double> > A(1024, 0);
            vector<complex<double> > B(1024, 0);
            //fftw_complex B[1024];
            h_fb_win_fxp.resize(L*M);
            double F,x;
             for (int n = 0; n < L*M ; ++n) {
                 F = (double)n / L/M;
                 x = K*(2*M*F-0.5); // rrerf
                 if (n < L*M/2)
                 A[n]= sqrt(0.5*erfc(x)) / 1023.9911405565; // 1024 вес
                 else
                 A[n]= A[L*M-n].real(); // Для симметрии мб индексация
             }
            fft((fftw_complex*) &A[0],(fftw_complex*) &B[0],1024, true);//(fftw_complex*)
            double max_coeff_val = B[0].real();
            coeff_radix = log2(pow(2,WIN_H_RADIX-1)/max_coeff_val);
            for (int n = 0; n < L*M ; ++n) { // fftshift
                  A[n] = B[(512+n)%1024].real() ;
                 // cout << B[(512+n)%1024].real()<< endl;
                  h_fb_win_fxp[n] = round( A[n].real() * pow(2,coeff_radix) );
            }
            N = 1024;

        }

        void get_signal()
        {
            int32_t maxSumm = fbWinMaxGain();
            //fb_analysis_win_max_gain_bit = ceil(max(log2(sum(abs(buffer(h_fb_win_fxp,NFFT)),2))));
            maxSummLog = log2(maxSumm) + 0.5; // ceil
            int16_t round_fft = coeff_radix-maxSummLog ;
            non_maximally_decimated_fb();
        }

        void non_maximally_decimated_fb()
        {
            //Y_sum = zeros(NFFT,overlapped_ratio*ceil(length(pulse_sig_round)/NFFT));
            int16_t size = FB_OVERLAP_RATIO * ceil(sig.si.size() / NFFT);
            vector<vector<double> > a(NFFT, vector<double>(size, 0));
            filtered.resize(sig.si.size() * FB_OVERLAP_RATIO );
            for(int n = 0; n < FB_OVERLAP_RATIO ; ++n)
            {
                pulse_sig_phase_n.resize(sig.si.size());
                for (int k = NFFT/FB_OVERLAP_RATIO*n; k < sig.si.size() ; ++k)
                {
                    //pulse_sig_phase_n[k] = sig.si[k];
                    pulse_sig_phase_n[k] = complex<double>(sig.si[k],sig.sq[k]);
                }
                maximally_decimated_fb();
            }

        }

        void maximally_decimated_fb()
        {
            vector<complex<double>> f(sig.si.size(),0);
            vector<complex<double>> f2(sig.si.size(),0);
            int longSize = 0;
            for(int n = 0; n < NFFT ; ++n)
            {                
                vector<double> F(sig.si.size(),0);

                longSize = sig.si.size()/NFFT ;
                for( int k = 0; k < longSize; ++k ) // filter 9
                {
                    complex<double> fiq(0);
                    for( int m = 0; m < WIN_OVERLAP_RATIO; ++m ) // 8
                    {
                        if( k - m > 0 )
                        {                            
                            int indH = 128*m +127 -n;
                            int indS =  (k - m)*128 + n ;
                            int32_t h_fir = h_fb_win_fxp[indH]; //m*128 + n
                            complex<double> sigPh = pulse_sig_phase_n[indS];
                            sigPh *= h_fir;
                            fiq +=  sigPh ;

                        }
                    }
                    f[k*128 + n] = fpga_round(fiq, maxSummLog);
                }
            }            
            for( int k = 0; k < longSize; ++k ) // filter 9
            {
                fft( (fftw_complex*) &f[128*k],(fftw_complex*) &f2[128*k],128, false);
                for( int n = 0; n < 128; ++n )
                {
                    f2[128*k + n] = fpga_round(f2[128*k + n], coeff_radix-maxSummLog);
                    filtered[]
                }
            }
            int t = 1;

            //Y_sum(:,k+(n-1)*overlapped_ratio) = Y_phase_n(:,n);
        }


        int32_t fbWinMaxGain()
        {
            int32_t summ = 1,maxSumm = 0; // =1 учитывает первый эл-т
            for (int n = 0; n < NFFT ; ++n) { // 128
                for (int m = 0; m < WIN_OVERLAP_RATIO ; ++m) { //8
                    summ = summ + abs(h_fb_win_fxp[n + m*128]);
                }
                if (maxSumm < summ)
                    maxSumm = summ;
                summ = 1;
            }
            return maxSumm;
        }

        void fft( fftw_complex* in, fftw_complex* out, uint32_t fft_size, bool backward )
        {
            if (backward)
            {
                fftw_plan plan = fftw_plan_dft_1d( fft_size, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);
                fftw_execute(plan);
                fftw_destroy_plan(plan);
            }
                else
            {
                fftw_plan plan = fftw_plan_dft_1d( fft_size, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
                fftw_execute(plan);
                fftw_destroy_plan(plan);
            }
            fftw_cleanup();
        }

        complex<double> fpga_round( complex<double> din, int32_t shval )
        {
           return complex<double>(floor(din.real()/ pow(2,shval) + 0.5), floor(din.imag()/ pow(2,shval) + 0.5)); // floor ??
        }

};


int main()
{
    float LFM_dev_hz = 50*1e6;
    //uint32_t LFM_dev_hz = 0;
    AnalysisBank comb;
    comb.genSignal(LFM_dev_hz);
    comb.createNpr();

    cout << "Hello World!" << endl;
    return 0;
}
