#ifndef C_EL_HAHN_H
#define C_EL_HAHN_H

#include "C_sigRMN.h"


class C_el_hahn : public C_sigRMN
{
public:
    C_el_hahn();

    //simulation of the monochrmatic momentum during a two pi/2-pulse sequence with a constant phase shift induced by the field inhomogeneities
    void monochromatic(double mx0, double my0, double mz0, double dw, bool rotFrame=true);
    void biChromatic(double mx0, double my0, double mz0, double dw, bool rotFrame=true);
    void multiChrome(double mx0, double my0, double mz0, int N, bool rotFrame=true);

    void monochromaticDiffusion(double mx0, double my0, double mz0, double dw, bool rotFrame=true);
    void multiChromeDiffusion(double mx0, double my0, double mz0, int N, bool rotFrame=true);

};

#endif // C_EL_HAHN_H
