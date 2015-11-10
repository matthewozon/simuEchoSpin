#include <iostream>
//#include <C_thread_bench.h>

#include "C_matrix.h"

#include "c_el_hahn.h"

//using namespace std;

int main()
{
    C_el_hahn spin;
    spin.Bz0=1.0;
    spin.Bxy=0.02;
    spin.gamma=1000;
    spin.T1=0.5; spin.T2=2.0; spin.tw=Pi/(2.0*spin.gamma*spin.Bxy); spin.TE=30.0*spin.tw;
    double dw=1.0*spin.gamma*spin.Bz0/100.0;
    spin.muz0 = 1.0;//1.0*(1.0+(dw/(spin.gamma*spin.Bz0))); the shift is computed in the methods of C_el_hahn
    spin.DW=spin.gamma*spin.Bz0/20.0;

    spin.D = 0.0*0.0000005;
    spin.G = 1.0;//0.001;

    //spin.monochromatic(0.0, 0.0, 1.0, 0.0*dw);
    //spin.biChromatic(0.0, 0.0, 1.0, dw);
    //spin.multiChrome(0.0, 0.0, 1.0, 200);
    //spin.monochromaticDiffusion(0.0, 0.0, 1.0, 0.0*dw);
    spin.multiChromeDiffusion(0.0, 0.0, 1.0, 200);
    return 0;
}



//ok pour echo
//C_el_hahn spin;
//spin.Bz0=1.0;
//spin.Bxy=0.02;
//spin.gamma=1000;
//spin.T1=0.5; spin.T2=2.0; spin.tw=Pi/(2.0*spin.gamma*spin.Bxy); spin.TE=30.0*spin.tw;
//double dw=1.0*spin.gamma*spin.Bz0/100.0;
//spin.muz0 = 1.0;//1.0*(1.0+(dw/(spin.gamma*spin.Bz0))); the shift is computed in the methods of C_el_hahn
//spin.DW=spin.gamma*spin.Bz0/10.0;

////spin.monochromatic(0.0, 0.0, 1.0, 0.0*dw, true);
////spin.biChromatic(0.0, 0.0, 1.0, dw, true);
//spin.multiChrome(0.0, 0.0, 1.0, 1000, true);
