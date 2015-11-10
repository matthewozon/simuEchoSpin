#include "C_sigRMN.h"

C_sigRMN::C_sigRMN()
{
}

void C_sigRMN::frameToFix(double t00)//TODO fix the continuity pb due to frame change
{
    double tmpX, tmpY;
    for(int i=0 ; i<MU.getNbRow() ; i++)
    {
        //std::cout << T(i,0)+t00 << std::endl;
        tmpX = cos(omega*(T(i,0)+t00))*MU(i,0)-sin(omega*(T(i,0)+t00))*MU(i,1);
        tmpY = sin(omega*(T(i,0)+t00))*MU(i,0)+cos(omega*(T(i,0)+t00))*MU(i,1);
        MU(i,0) = tmpX;
        MU(i,1) = tmpY;
    }
    return;
}

void C_sigRMN::frameToRotating(double t00)
{
    double tmpX, tmpY;
    for(int i=0 ; i<MU.getNbRow() ; i++)
    {
        tmpX = cos(omega*(T(i,0)+t00))*MU(i,0)+sin(omega*(T(i,0)+t00))*MU(i,1);
        tmpY = -sin(omega*(T(i,0)+t00))*MU(i,0)+cos(omega*(T(i,0)+t00))*MU(i,1);
        MU(i,0) = tmpX;
        MU(i,1) = tmpY;
    }
    return;
}

double C_sigRMN::generateRMNsignal(double mx0, double my0, double mz0, double dt, double dw)
{
    //init the arrays
    MU.resize(sigLength,3);
    T.resize(sigLength,1);
    MU(0,0)=mx0; MU(0,1)=my0; MU(0,2)=mz0;
    double t=0.0;
    T(0,0) = 0.0;

    //the constant of the problem
    double A, B, phi;
    //double Bxy = sqrt(Bx*Bx+By*By);
    double omega1 = gamma*Bxy;
    double wd = sqrt(SQR(omega1) + SQR(dw));

    if(ABS(omega1)<SMALL_NUM_F)
    {
        //the amplitude of the perturbating field is to small to account for any effect of RMN
        if(dw<SMALL_NUM_F)
        {
            //all components of the momentum are set to the initial value, it is a steady state
            for(int i=1 ; i<sigLength ; i++)
            {
                t += dt;
                T(i,0)=t;
                MU(i,0) = mx0;
                MU(i,1) = my0;
                MU(i,2) = mz0;
            }
            return t;
        }
        else
        {
            //the rotating frame is not well synchronized with the Larmor frequency imposedd by the magnetic field
            A=sqrt(SQR(mx0)+SQR(my0));
            phi = atan2(my0,mx0);
            for(int i=1 ; i<sigLength ; i++)
            {
                t += dt;
                T(i,0)=t;
                MU(i,0) = A*cos(dw*t+phi);
                MU(i,1) = A*sin(dw*t+phi);
                MU(i,2) = mz0;
            }
            return t;
        }
    }

    //sample the signal for the general case

    A=sqrt(SQR((dw/wd)*mx0 - (omega1/wd)*mz0)+SQR(my0));
    phi=atan2(my0,(dw/wd)*mx0 - (omega1/wd)*mz0);
    B=(mx0+(dw/omega1)*mz0)/(1+SQR(dw/omega1));

    //std::cout << "A " << A << " B " << B << " phi " << phi << std::endl;

    for(int i=1 ; i<sigLength ; i++)
    {
        t += dt;
        T(i,0)=t;
        MU(i,0) = (dw/wd)*A*cos(wd*t+phi) + B;
        MU(i,1) = A*sin(wd*t+phi);
        MU(i,2) = -(omega1/wd)*A*cos(wd*t+phi) + (dw/omega1)*B;
    }
    return t;
}

double C_sigRMN::generateFIDsignal(double mx0, double my0, double mz0, double dt, double dw)
{
    //init the arrays
    MU.resize(sigLength,3);
    T.resize(sigLength,1);
    MU(0,0)=mx0; MU(0,1)=my0; MU(0,2)=mz0;
    double t=0.0;
    T(0,0) = 0.0;

    //sample the signal for the general case
    double A = sqrt(SQR(mx0)+SQR(my0));
    double phi = atan2(my0,mx0);

    for(int i=1 ; i<sigLength ; i++)
    {
        t += dt;
        T(i,0)=t;
        MU(i,0) = A*exp(-t/T2)*cos(dw*t+phi);
        MU(i,1) = A*exp(-t/T2)*sin(dw*t+phi);
        MU(i,2) = muz0 + (mz0-muz0)*exp(-t/T1);
    }
    return t;
}

double C_sigRMN::generateFIDsignalWithDiffusion(double mx0, double my0, double mz0, double dt, double dw)
{
    //init the arrays
    MU.resize(sigLength,3);
    T.resize(sigLength,1);
    MU(0,0)=mx0; MU(0,1)=my0; MU(0,2)=mz0;
    double t=0.0;
    T(0,0) = 0.0;

    //sample the signal for the general case
    double A = sqrt(SQR(mx0)+SQR(my0));
    double phi = atan2(my0,mx0);

    for(int i=1 ; i<sigLength ; i++)
    {
        t += dt;
        T(i,0)=t;
        MU(i,0) = A*exp(-((t/T2) +(SQR(gamma*G)*D*t*t*t)) )*cos(dw*t+phi);
        MU(i,1) = A*exp(-((t/T2) +(SQR(gamma*G)*D*t*t*t)) )*sin(dw*t+phi);
        MU(i,2) = muz0 + (mz0-muz0)*exp(-t/T1);
    }
    return t;
}
