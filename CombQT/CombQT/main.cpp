#include <iostream>
#include <vector>
#include <cmath>
#include <complex>
#include <algorithm>
#include <fftw3.h>
#include "time.h"
/*
#ifdef _DEBUG
#include <crtdbg.h>
#define _CRTDBG_MAP_ALLOC
#endif
*/
#define WIN_H_RADIX        18
#define FB_OVERLAP_RATIO  2 // перекрытие
#define NFFT  128 // колво фильтров
#define WIN_OVERLAP_RATIO  8 // длина фильтров


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

        void genSignal(float LFM) // генерирование сигнала
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
                    sig.q[n] = round(DATA_AMPL*sin(2.0*M_PI*1e6*sig.t[n]));// +rand()*10000 ;
                    sig.i[n] = round(DATA_AMPL*cos(2.0*M_PI*1e6*sig.t[n]));// +rand()*10000;
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
           /* FILE * fp = fopen("/home/anatoly/hub/Comb/I", "wb"); // Запись в файл
            FILE * fp2 = fopen("/home/anatoly/hub/Comb/Q", "wb");
            fwrite(&sig.si[0], sizeof(int16_t), sig.si.size(), fp) ;
            fwrite(&sig.sq[0], sizeof(int16_t), sig.sq.size(), fp2) ;
            fclose(fp);
            fclose(fp2);*/
        }

        void createNpr() //фильтрация сигнала
        {
            npr_coeff(NFFT,2*WIN_OVERLAP_RATIO); // создание коэф.
            get_signal(); // сама фильтрация
        }



        vector<int32_t> h_fb_win_fxp;
        signal sig;
    private:
        uint32_t Fs = 512e6;
        //uint32_t BAND              = 500e6;
        //uint32_t FW                = 4e6;
        double t_us = 0.2 * 1e-6;
        //double period_us = 10 * 1e-6;
        uint16_t first_period_part = NFFT ;
        uint16_t second_period_part = NFFT * WIN_OVERLAP_RATIO;
        vector<complex<double>> pulse_sig_phase_n;
        int32_t maxSummLog;
        vector<complex<double>> filtered;
        vector<complex<double>> sigOut;



        uint16_t coeff_radix;

        void npr_coeff(int16_t N,int16_t L) // генерирование коэф. фильтра
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
            int32_t maxSumm = fbWinMaxGain(); // максимальное значение
            //fb_analysis_win_max_gain_bit = ceil(max(log2(sum(abs(buffer(h_fb_win_fxp,NFFT)),2))));
            maxSummLog = log2(maxSumm) + 0.5; // ceil
            int16_t round_fft = coeff_radix-maxSummLog ;
            non_maximally_decimated_fb(); // создание ан. гребенки фильтров

            /*FILE * fp = fopen("/home/anatoly/hub/Comb/Anal", "wb");
            fwrite(&filtered[0], sizeof(complex<double>), filtered.size(), fp) ;
            fclose(fp);*/

            //npr_synthesis(); //Синтезирующая гребенка
        }

        void npr_synthesis()  // Синтезирующая гребенка
        {
            uint16_t sizeFiltered = filtered.size();
            uint16_t rows = sizeFiltered / NFFT;
            vector<complex<double> > yfft(sizeFiltered , 0);
            for( int k = 0; k < rows ; ++k ) // for each row = 18 Фурье
            {
                fft( (fftw_complex*) &filtered[NFFT*k],(fftw_complex*) &yfft[NFFT*k],NFFT, true);
            }
            for(int n = 0; n < NFFT ; ++n)  // Фильтрация (свертка)
            {
                int16_t longSize = sig.si.size()/NFFT ;
                for( int k = 0; k < longSize; ++k ) // filter 9
                {
                    complex<double> fiq1(0);
                    complex<double> fiq2(0);
                    for( int m = 0; m < WIN_OVERLAP_RATIO; ++m ) // 8
                    {
                        if( k - m >= 0 )
                        {
                            //int indH = 128*m +127 -n;
                            int indH = 128*m + n;
                            int indS1 =  (k - m)*256 + n ;
                            int indS2 =  (k - m)*256 + (64 + n)%128 +128; //192 is cyclic shift
                            int32_t h_fir = h_fb_win_fxp[indH];
                            complex<double> sigPh1 = fpga_round(yfft[indS1],1);
                            complex<double> sigPh2 = fpga_round(yfft[indS2],1);
                            if (k*NFFT + n + longSize*NFFT > 2237) // дебаг
                                int fl=1;
                            sigPh1 *= h_fir;
                            sigPh2 *= h_fir;
                            fiq1 +=  sigPh1 ;
                            fiq2 +=  sigPh2 ;

                        }
                    }
                    filtered[k*NFFT + n] = fpga_round(fiq1, 16);
                    filtered[k*NFFT + n + longSize*NFFT] = fpga_round(fiq2, 16); //For debug only

                }
            }
            // combine filter results
            sigOut.resize(sig.si.size());
            int ind=0;
            for(int n = 0; n < sig.si.size() ; ++n) // Сложение двух перекрытий
            {
                if (n < NFFT/2)
                    sigOut[n] = filtered[n] ;
                else
                    sigOut[n] = filtered[n] + filtered[sizeFiltered/2 + ind++];
            }
            int flag = 1;
            output();


        }

        void output() // вывод в файл
        {
            FILE * fp = fopen("/home/anatoly/hub/Comb/Signal", "wb");
            fwrite(&sigOut[0], sizeof(complex<double>), sigOut.size(), fp) ;
            fclose(fp);
        }

        void non_maximally_decimated_fb() // Анализирующая гребенка
        {
            //Y_sum = zeros(NFFT,overlapped_ratio*ceil(length(pulse_sig_round)/NFFT));
            int16_t size = FB_OVERLAP_RATIO * ceil(sig.si.size() / NFFT);
            vector<vector<double> > a(NFFT, vector<double>(size, 0));
            filtered.resize(sig.si.size() * FB_OVERLAP_RATIO );
            //int *x= new int ();
            for(int n = 0; n < FB_OVERLAP_RATIO ; ++n)
            {
                pulse_sig_phase_n.clear();
                pulse_sig_phase_n.resize(sig.si.size());
                int32_t ind = 0;
                for (int k = NFFT/FB_OVERLAP_RATIO*n; k < sig.si.size() ; ++k)
                {
                    //pulse_sig_phase_n[k] = sig.si[k];
                    pulse_sig_phase_n[ind++] = complex<double>(sig.si[k],sig.sq[k]);
                }
                maximally_decimated_fb(n);
            }

        }

        void maximally_decimated_fb(int16_t ovRat)
        {
            vector<complex<double>> f(sig.si.size(),0);
            vector<complex<double>> f2(sig.si.size(),0);
            int longSize = 0;
            for(int n = 0; n < NFFT ; ++n) // Фильтрация
            {                
                vector<double> F(sig.si.size(),0);

                longSize = sig.si.size()/NFFT ;
                for( int k = 0; k < longSize; ++k ) // filter 9
                {
                    complex<double> fiq(0);
                    for( int m = 0; m < WIN_OVERLAP_RATIO; ++m ) // 8
                    {
                        if( k - m >= 0 )
                        {                            
                            int indH = 128*m +127 -n;

                            int indS =  (k - m)*128 + n ;
                            int32_t h_fir = h_fb_win_fxp[indH]; //m*128 + n
                            complex<double> sigPh = pulse_sig_phase_n[indS];
                            sigPh *= h_fir;
                            fiq +=  sigPh ;

                        }
                    }
                    f[k*NFFT + n] = fpga_round(fiq, maxSummLog); // округление

                }
            }
            //circshift(&f[0],sig.si.size(),131);
            //filtered.clear();
            //filtered.resize(sig.si.size() * FB_OVERLAP_RATIO );
            for( int k = 0; k < longSize; ++k ) // Фильтрация(свертка) filter 9
            {
                if (ovRat>0)
                rotate(&f[NFFT*k],&f[NFFT*(k+1)-NFFT/FB_OVERLAP_RATIO*ovRat],&f[NFFT*(k+1)]);

                fft( (fftw_complex*) &f[NFFT*k],(fftw_complex*) &f2[NFFT*k],NFFT, false);
                for( int n = 0; n < NFFT; ++n )
                {
                    //f2[128*k + n] = fpga_round(f2[128*k + n], coeff_radix-maxSummLog);
                    int ind1 = NFFT*(FB_OVERLAP_RATIO*k+ovRat) + n;
                    int ind2 = NFFT*k + n;
                    filtered[ind1] = fpga_round(f2[ind2], coeff_radix-maxSummLog);

                }
            }
            int flag = 1;

        }


        int32_t fbWinMaxGain() // максимум усиления
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

        void fft( fftw_complex* in, fftw_complex* out, uint32_t fft_size, bool backward ) //Фурье
        {
            if (backward)
            {
                fftw_plan plan = fftw_plan_dft_1d( fft_size, in, out, FFTW_BACKWARD,  FFTW_ESTIMATE);
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

        complex<double> fpga_round( complex<double> din, int32_t shval ) // округление
        {
           return complex<double>(floor(din.real()/ pow(2,shval) + 0.5), floor(din.imag()/ pow(2,shval) + 0.5)); // floor ??
        }

        void circshift( complex<double>* x,unsigned sizex, unsigned pos ) // собственный сдвиг
        {
            complex<double> temp;
            unsigned iter = ( sizex & 0x1 ) ? sizex/2 + 1 : sizex/2;
            //unsigned pos = iter;
            for( unsigned i = 0; i < iter; ++i )
            {
                temp = x[(pos + i)%sizex];
                x[(pos + i)%sizex] = x[i];
                x[i] = temp;
            }
        }

};


int main()
{

    float LFM_dev_hz = 50*1e6;
     LFM_dev_hz = 0;
    AnalysisBank comb;
    comb.genSignal(LFM_dev_hz);
    //clock_t t = clock();

    comb.createNpr();

    //t=clock()-t;
    //t=float(t);//t время в секундах )
    //cout << "Time:  "<< t << endl;
    return 0;
}
