#include "analysisbank.cpp"
#include "generator.cpp"
#include "time.h"


int main()
{

    float LFM_dev_hz = 50*1e6;
     LFM_dev_hz = 0;
    AnalysisBank comb;
    Generator gen;
    gen.genSignal(LFM_dev_hz); // сгенерировать сигнал
    comb.readSignal(gen.sig.sq, gen.sig.si); // передать в гребенку
    comb.saveSignal(); // Сохранить сигнал в файл
    //comb.openSignal(1152);

    clock_t t = clock();
   for (int n=0;n < 1000;++n)
        comb.createNpr();    
    t=clock()-t;
    double tSec;
    tSec=(double(t)) / CLOCKS_PER_SEC;//t время в секундах )
    cout << "Time:  "<< tSec << endl;

    vector<vector<complex<double>>> F = comb.getAnalyzeFB(); // получить результат с анал. гребенки
    comb.saveAnalyzeFB(); // сохранить анал. гребенку
    return 0;
}
