#include <iostream>
#include <complex>
#include <vector>

#include <fftw3.h>

using namespace std;
/////////////////////////////////////////////////////////////////////////////
int main()
{
    int k =1;
   // создаем одномерную выборку, все значения которой равны 1
   vector<complex<double> > data(64, 1.);
   vector<complex<double> > data2(64, 0);

   // создаем план для библиотеки fftw
   fftw_plan plan=fftw_plan_dft_1d(data.size(), (fftw_complex*) &data[0], (fftw_complex*) &data2[0], FFTW_BACKWARD, FFTW_ESTIMATE);

   // преобразование Фурье
   fftw_execute(plan);
   fftw_destroy_plan(plan);

   // выводим в файл результат преобразования Фурье (должна получиться Дельта-функция)

   for(size_t i=0; i<data.size(); ++i)
   {
      cout<<data[i]<<endl;
   }
   return 0;
}
