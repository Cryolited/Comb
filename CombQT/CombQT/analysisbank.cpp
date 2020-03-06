#include "analysisbank.h"

AnalysisBank::AnalysisBank()
{

}

void AnalysisBank::setParams(params pFun) // Установить параметры
{
    p = pFun;
}

void AnalysisBank::openSignal(uint32_t size) //1152 открыть записанный сигнал
{
    sig.si.resize(size);
    sig.sq.resize(size);

    FILE * fp = fopen("/home/anatoly/hub/Comb/I", "rb"); // Чтение из файла
    FILE * fp2 = fopen("/home/anatoly/hub/Comb/Q", "rb");
    fread(&sig.si[0], sizeof(int16_t), size, fp) ;
    fread(&sig.sq[0], sizeof(int16_t), size, fp2) ;
    fclose(fp);
    fclose(fp2);
}

void AnalysisBank::readSignal(vector<int16_t>& inVecQ, vector<int16_t>& inVecI) // считать вектор входного сигнала
{
    sig.si.resize(inVecI.size());
    sig.sq.resize(inVecQ.size());

    copy ( inVecQ.begin(), inVecQ.end(), sig.sq.begin() );
    copy ( inVecI.begin(), inVecI.end(), sig.si.begin() );
}

void AnalysisBank::readSignal(int16_t inVecQ[], int16_t inVecI[], uint32_t size) // считать массив входного сигнала
{
    sig.si.resize(size);
    sig.sq.resize(size);

    copy ( inVecQ, inVecQ+size, sig.sq.begin() );
    copy ( inVecI, inVecI+size, sig.si.begin() );
}


void AnalysisBank::saveSignal() // записать в файл сгенерированный сигнал
{
     FILE * fp = fopen("/home/anatoly/hub/Comb/I", "wb"); // Запись в файл
     FILE * fp2 = fopen("/home/anatoly/hub/Comb/Q", "wb");
     fwrite(&sig.si[0], sizeof(int16_t), sig.si.size(), fp) ;
     fwrite(&sig.sq[0], sizeof(int16_t), sig.sq.size(), fp2) ;
     fclose(fp);
     fclose(fp2);
}

void AnalysisBank::createNpr() //фильтрация сигнала
{
    npr_coeff(p.NFFT,2*p.WIN_OVERLAP_RATIO); // создание коэф.
}


void AnalysisBank::saveAnalyzeFB() // записать в файл спек-му Анализирующей Гребенки
{
    FILE * fp = fopen("/home/anatoly/hub/Comb/Anal", "wb");
    for (int n = 0; n < filt.size() ; ++n) {
     fwrite(&filt[n][0], sizeof(complex<double>), filt[0].size(), fp) ;
    }
    fclose(fp);
}

vector<vector<complex<double>>>& AnalysisBank::getAnalyzeFB() // получить спек-му Анализирующей Гребенки
{
    return filt;
}

void AnalysisBank::npr_coeff(int16_t N,int16_t L) // генерирование коэф. фильтра
{
    //if (L == 16)
      double K=5.856;
    int16_t M = N /2;
    N = p.NFFT * p.WIN_OVERLAP_RATIO;
    vector<complex<double> > A(N, 0);
    vector<complex<double> > B(N, 0);
    //fftw_complex B[1024];
    h_fb_win_fxp.resize(L*M);
    double F,x;
     for (int n = 0; n < L*M ; ++n) {
         F = double(n) / L/M;
         x = K*(2*M*F-0.5); // rrerf
         if (n < L*M/2)
         A[n]= sqrt(0.5*erfc(x)) /1024;// 1023.9911405565; // 1024 вес
         else
         A[n]= A[L*M-n].real(); // Для симметрии (мб индексация?)
     }
    fft((fftw_complex*) &A[0],(fftw_complex*) &B[0],N, true);//(fftw_complex*)
    double max_coeff_val = B[0].real();
    coeff_radix = log2(pow(2,p.WIN_H_RADIX-1)/max_coeff_val);
    for (int n = 0; n < L*M ; ++n) { // fftshift
          A[n] = B[(N/2+n)%N].real() ;
         // cout << B[(512+n)%1024].real()<< endl;
          h_fb_win_fxp[n] = round( A[n].real() * pow(2,coeff_radix) );
    }


}

void AnalysisBank::filterAnalyze()
{
    int32_t maxSumm = fbWinMaxGain(); // максимальное значение
    maxSummLog = log2(maxSumm) + 0.5; // ceil
    //int16_t round_fft = coeff_radix-maxSummLog ;
    non_maximally_decimated_fb(); // фильтрация ан. гребенкой фильтров
    //npr_synthesis(); //Синтезирующая гребенка
    //saveAnalyzeFB();
}

void AnalysisBank::npr_synthesis()  // Синтезирующая гребенка
{
    uint16_t sizeFiltered = filtered.size();
    uint16_t rows = sizeFiltered / p.NFFT;
    vector<complex<double> > yfft(sizeFiltered , 0);
    for( int k = 0; k < rows ; ++k ) // for each row = 18 Фурье
    {
        fft( (fftw_complex*) &filtered[p.NFFT*k],(fftw_complex*) &yfft[p.NFFT*k],p.NFFT, true);
    }
    for(int n = 0; n < p.NFFT ; ++n)  // Фильтрация (свертка)
    {
        int16_t longSize = sig.si.size()/p.NFFT ;
        for( int k = 0; k < longSize; ++k ) // filter 9
        {
            complex<double> fiq1(0);
            complex<double> fiq2(0);
            for( int m = 0; m < p.WIN_OVERLAP_RATIO; ++m ) // 8
            {
                if( k - m >= 0 )
                {
                    //int indH = 128*m +127 -n;
                    int indH = p.NFFT*m + n;
                    int indS1 =  (k - m)*p.NFFT*2 + n ;
                    int indS2 =  (k - m)*p.NFFT*2 + (p.NFFT/2 + n)%p.NFFT +p.NFFT; //192 is cyclic shift
                    int32_t h_fir = h_fb_win_fxp[indH];
                    complex<double> sigPh1 = fpga_round(yfft[indS1],1);
                    complex<double> sigPh2 = fpga_round(yfft[indS2],1);
                    //if (k*NFFT + n + longSize*NFFT > 2237) // дебаг
                      //  int fl=1;
                    sigPh1 *= h_fir;
                    sigPh2 *= h_fir;
                    fiq1 +=  sigPh1 ;
                    fiq2 +=  sigPh2 ;

                }
            }
            filtered[k*p.NFFT + n] = fpga_round(fiq1, 16);
            filtered[k*p.NFFT + n + longSize*p.NFFT] = fpga_round(fiq2, 16); //For debug only

        }
    }
    // combine filter results
    sigOut.resize(sig.si.size());
    int ind=0;
    for(int n = 0; n < sig.si.size() ; ++n) // Сложение двух перекрытий
    {
        if (n < p.NFFT/2)
            sigOut[n] = filtered[n] ;
        else
            sigOut[n] = filtered[n] + filtered[sizeFiltered/2 + ind++];
    }
    int flag = 1;
    saveCreateFB();


}

void AnalysisBank::saveCreateFB() // вывод в файл сигнала после Синетезирующей Гребенки
{
    FILE * fp = fopen("/home/anatoly/hub/Comb/Signal", "wb");
    fwrite(&sigOut[0], sizeof(complex<double>), sigOut.size(), fp) ;
    fclose(fp);
}

void AnalysisBank::non_maximally_decimated_fb() // Анализирующая гребенка
{
    //Y_sum = zeros(NFFT,overlapped_ratio*ceil(length(pulse_sig_round)/NFFT));
    uint16_t size = p.FB_OVERLAP_RATIO * ceil(sig.si.size() / p.NFFT);
    //vector<vector<double> > a(NFFT, vector<double>(size, 0));
    filt.assign(p.NFFT, vector<complex<double>>(size));
    //filtered.resize(sig.si.size() * FB_OVERLAP_RATIO );
    //int *x= new int (); // ПРОВЕРКА УТЕЧЕК
    for(int16_t n = 0; n < p.FB_OVERLAP_RATIO ; ++n)
    {
        //pulse_sig_phase_n.clear();
        //pulse_sig_phase_n.resize(sig.si.size());
        //int32_t ind = 0;
        /*for (int k = NFFT/FB_OVERLAP_RATIO*n; k < sig.si.size() ; ++k)
        {
            //pulse_sig_phase_n[k] = sig.si[k];
            pulse_sig_phase_n[ind++] = complex<double>(sig.si[k],sig.sq[k]);
        }*/
        maximally_decimated_fb(n);
    }

}

void AnalysisBank::maximally_decimated_fb(int16_t ovRat)
{
    vector<complex<double>> f(sig.si.size(),0);
    vector<complex<double>> f2(sig.si.size(),0);
    int longSize = 0;
    vector<double> F(sig.si.size(),0);
    int32_t h_fir = 0;
    uint16_t indH,indS;
    complex<double> sigPh;
    for(int n = 0; n < p.NFFT ; ++n) // Фильтрация(свертка)
    {
        longSize = sig.si.size()/p.NFFT ;
        for( int k = 0; k < longSize; ++k ) // filter 9
        {
            complex<double> fiq(0);
            for( int m = 0; m < p.WIN_OVERLAP_RATIO; ++m ) // 8
            {
                if( k - m >= 0 )
                {
                    indH = p.NFFT*m +p.NFFT -1 -n;
                    indS = ( (k - m)*p.NFFT + n + p.NFFT/p.FB_OVERLAP_RATIO*ovRat )%sig.si.size();
                    //indS =  (k - m)*128 + n ;
                    h_fir = h_fb_win_fxp[indH]; //m*128 + n
                    sigPh = complex<double>(sig.si[indS],sig.sq[indS]);
                    sigPh *= h_fir;
                    fiq +=  sigPh ;

                }
            }
            f[k*p.NFFT + n] = fpga_round(fiq, maxSummLog); // округление

        }
    }
    //circshift(&f[0],sig.si.size(),131);
    //filtered.clear();
    //filtered.resize(sig.si.size() * FB_OVERLAP_RATIO );
    uint16_t ind1, ind2;
    for( int k = 0; k < longSize; ++k ) // Фурье
    {
        if (ovRat>0)
        rotate(&f[p.NFFT*k],&f[p.NFFT*(k+1)-p.NFFT/p.FB_OVERLAP_RATIO*ovRat],&f[p.NFFT*(k+1)]);

        fft( (fftw_complex*) &f[p.NFFT*k], (fftw_complex*) &f2[p.NFFT*k], p.NFFT, false);
        for( int n = 0; n < p.NFFT; ++n )
        {
            //f2[128*k + n] = fpga_round(f2[128*k + n], coeff_radix-maxSummLog);
            //ind1 = NFFT*(FB_OVERLAP_RATIO*k+ovRat) + n;
            ind2 = p.NFFT*k + n;
            //filtered[ind1] = fpga_round(f2[ind2], coeff_radix-maxSummLog);
            filt[n][p.FB_OVERLAP_RATIO*k+ovRat] = fpga_round(f2[ind2], coeff_radix-maxSummLog);

        }
    }
    int flag = 1; // for debug only

}


int32_t AnalysisBank::fbWinMaxGain() // максимум усиления
{
    int32_t summ = 1,maxSumm = 0; // =1 учитывает первый эл-т
    for (int n = 0; n < p.NFFT ; ++n) { // 128
        for (int m = 0; m < p.WIN_OVERLAP_RATIO ; ++m) { //8
            summ = summ + abs(h_fb_win_fxp[n + m*p.NFFT]);
        }
        if (maxSumm < summ)
            maxSumm = summ;
        summ = 1;
    }
    return maxSumm;
}

void AnalysisBank::fft( fftw_complex* in, fftw_complex* out, uint32_t fft_size, bool backward ) //Фурье
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

complex<double> AnalysisBank::fpga_round( complex<double> din, int32_t shval ) // округление
{
   return complex<double>(floor(din.real()/ pow(2,shval) + 0.5), floor(din.imag()/ pow(2,shval) + 0.5)); // floor ??
}

void AnalysisBank::circshift( complex<double>* x,unsigned sizex, unsigned pos ) // собственный сдвиг
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




