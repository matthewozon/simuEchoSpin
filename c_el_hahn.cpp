#include "c_el_hahn.h"

C_el_hahn::C_el_hahn()
{
}

void C_el_hahn::monochromatic(double mx0, double my0, double mz0, double dw, bool rotFrame)
{
    //the perturbating magnetic field pulsation is gamma*Bz0, but the actual constant field is Bz0+(dw/gamma)
    double dt;
    double tmpBz0 = Bz0;
    double tmpBxy = Bxy;
    //gamma=1000.0;//for the water the actual value is about 42.58 MHz/T
    //bulshit value for water in a magnetic field of 1T at temperature 10^-30K, but it is way easier if the purpose is to visualize the simulated signal mu0is equivalent to (hbar*gamma)^2*B0/(4kT)
    //mux0=0.0; muy0=0.0; muz0=1.0;
    //T1=0.5; T2=0.01;
    double tmpX=mx0, tmpY=my0, tmpZ=mz0;
    double tempMuz0=muz0;
    muz0 = tempMuz0*(1.0+(dw/(gamma*tmpBz0)));

    /******************************************************/
    /** initialization */
    //apply the constant magnetic field during t0
    Bz0=tmpBz0+(dw/gamma);
    Bxy=0.0; omega=0;
    double t0=MAX(T1,T2);
    sigLength=200;
    dt=t0/((double) sigLength);
    double t00=0.0;

    //suppose that the magnetic momentum at rest is null
    t00+=generateFIDsignal(tmpX, tmpY, tmpZ, dt, dw);
    tmpX=MU(MU.endL,0); tmpY=MU(MU.endL,1); tmpZ=MU(MU.endL,2);
    if(!rotFrame) frameToFix(t00);
    MU.save("initSig.txt");
    T.save("initSigTime.txt");
    std::cout << t00 << std::endl;


    /******************************************************/
    /** first pulse */
    //play one pi/2-pulse
    Bz0=tmpBz0+(dw/gamma);
    Bxy=tmpBxy; omega=gamma*Bz0;
    sigLength=4000;
    //t90 = Pi/(2.0*gamma*Bxy);
    dt=tw/((double) sigLength);

    //suppose that the magnetic momentum at rest is null
    t00+=generateRMNsignal(tmpX, tmpY, tmpZ, dt, dw);
    tmpX=MU(MU.endL,0); tmpY=MU(MU.endL,1); tmpZ=MU(MU.endL,2);
    if(!rotFrame) frameToFix(t00);
    MU.save("p90.txt");
    T.save("p90Time.txt");
    std::cout << t00 << std::endl;


    /******************************************************/
    /** free induction decay */
    //stop the radiofrequence pulse until the next pulse
    Bz0=tmpBz0+(dw/gamma);
    Bxy=0.0; omega=0;
    sigLength=4000;
    dt=(0.5*TE - tw)/((double) sigLength);
    //suppose that the magnetic momentum at rest is null
    t00+=generateFIDsignal(tmpX, tmpY, tmpZ, dt, dw);
    tmpX=MU(MU.endL,0); tmpY=MU(MU.endL,1); tmpZ=MU(MU.endL,2);
    if(!rotFrame) frameToFix(t00);
    MU.save("fid.txt");
    T.save("fidTime.txt");



    /******************************************************/
    /** second pulse */
    //play the pi-pulse
    Bz0=tmpBz0+(dw/gamma);
    Bxy=tmpBxy; omega=gamma*Bz0;
    sigLength=4000;
    dt=2.0*tw/((double) sigLength);

    //suppose that the magnetic momentum at rest is null
    t00+=generateRMNsignal(tmpX, tmpY, tmpZ, dt, dw);
    tmpX=MU(MU.endL,0); tmpY=MU(MU.endL,1); tmpZ=MU(MU.endL,2);
    if(!rotFrame) frameToFix(t00);
    MU.save("p180.txt");
    T.save("p180Time.txt");
    std::cout << t00 << std::endl;


    /******************************************************/
    /** spin echo */
    //stop the radiofrequence pulse until the next pulse
    Bz0=tmpBz0+(dw/gamma);
    Bxy=0.0; omega=0;
    sigLength=2000;
    dt=(TE-tw)/((double) sigLength);

    //suppose that the magnetic momentum at rest is null
    t00+=generateFIDsignal(tmpX, tmpY, tmpZ, dt, dw);
    if(!rotFrame) frameToFix(t00);
    MU.save("spinEcho.txt");
    T.save("spinEchoTime.txt");

    //restore the parameters
    Bz0=tmpBz0;
    Bxy=tmpBxy; omega=gamma*Bz0;
    muz0=tempMuz0;
    return;
}

void C_el_hahn::biChromatic(double mx0, double my0, double mz0, double dw, bool rotFrame)
{
    C_matrix<double> MUsteady1;
    C_matrix<double> MUfisrtPulse1;
    C_matrix<double> MUfirstFID1;
    C_matrix<double> MUsecondPulse1;
    C_matrix<double> MUsecondFID1;
    //the perturbating magnetic field pulsation is gamma*Bz0, but the actual constant field is Bz0+(dw/gamma)
    double dt;
    double tmpBz0 = Bz0;
    double tmpBxy = Bxy;
    double tempMuz0=muz0;
    muz0 = tempMuz0*(1.0+(dw/(gamma*tmpBz0)));
    std::cout << 1.0+(dw/(gamma*tmpBz0)) << " " << tempMuz0 << " " << muz0 << std::endl;
    //gamma=1000.0;//for the water the actual value is about 42.58 MHz/T
    //bulshit value for water in a magnetic field of 1T at temperature 10^-30K, but it is way easier if the purpose is to visualize the simulated signal mu0is equivalent to (hbar*gamma)^2*B0/(4kT)
    //mux0=0.0; muy0=0.0; muz0=1.0;
    //T1=0.5; T2=0.01;
    double tmpX=mx0, tmpY=my0, tmpZ=mz0;

    /*****************************************************************/
    /** We first consider the simulation with a slightly larger field*/
    /******************************************************/
    /** initialization */
    //apply the constant magnetic field during t0
    Bz0=tmpBz0+(dw/gamma);
    Bxy=0.0; omega=0;
    double t0=MAX(T1,T2);
    sigLength=200;
    dt=t0/((double) sigLength);
    double t00=0.0;

    //suppose that the magnetic momentum at rest is null
    t00+=generateFIDsignal(tmpX, tmpY, tmpZ, dt, dw);
    tmpX=MU(MU.endL,0); tmpY=MU(MU.endL,1); tmpZ=MU(MU.endL,2);
    if(!rotFrame) frameToFix(t00);
    MUsteady1=MU;
    //MU.save("initSig.txt");
    //T.save("initSigTime.txt");
    std::cout << t00 << std::endl;


    /******************************************************/
    /** first pulse */
    //play one pi/2-pulse
    Bz0=tmpBz0+(dw/gamma);
    Bxy=tmpBxy; omega=gamma*Bz0;
    sigLength=4000;
    //t90 = Pi/(2.0*gamma*Bxy);
    dt=tw/((double) sigLength);

    //suppose that the magnetic momentum at rest is null
    t00+=generateRMNsignal(tmpX, tmpY, tmpZ, dt, dw);
    tmpX=MU(MU.endL,0); tmpY=MU(MU.endL,1); tmpZ=MU(MU.endL,2);
    if(!rotFrame) frameToFix(t00);
    MUfisrtPulse1=MU;
    //MU.save("p90.txt");
    //T.save("p90Time.txt");
    std::cout << t00 << std::endl;


    /******************************************************/
    /** free induction decay */
    //stop the radiofrequence pulse until the next pulse
    Bz0=tmpBz0+(dw/gamma);
    Bxy=0.0; omega=0;
    sigLength=4000;
    dt=(0.5*TE - tw)/((double) sigLength);
    //suppose that the magnetic momentum at rest is null
    t00+=generateFIDsignal(tmpX, tmpY, tmpZ, dt, dw);
    tmpX=MU(MU.endL,0); tmpY=MU(MU.endL,1); tmpZ=MU(MU.endL,2);
    if(!rotFrame) frameToFix(t00);
    MUfirstFID1=MU;
    //MU.save("fid.txt");
    //T.save("fidTime.txt");



    /******************************************************/
    /** second pulse */
    //play the pi-pulse
    Bz0=tmpBz0+(dw/gamma);
    Bxy=tmpBxy; omega=gamma*Bz0;
    sigLength=4000;
    dt=2.0*tw/((double) sigLength);

    //suppose that the magnetic momentum at rest is null
    t00+=generateRMNsignal(tmpX, tmpY, tmpZ, dt, dw);
    tmpX=MU(MU.endL,0); tmpY=MU(MU.endL,1); tmpZ=MU(MU.endL,2);
    if(!rotFrame) frameToFix(t00);
    MUsecondPulse1=MU;
    //MU.save("p180.txt");
    //T.save("p180Time.txt");
    std::cout << t00 << std::endl;


    /******************************************************/
    /** spin echo */
    //stop the radiofrequence pulse until the next pulse
    Bz0=tmpBz0+(dw/gamma);
    Bxy=0.0; omega=0;
    sigLength=2000;
    dt=(TE-tw)/((double) sigLength);

    //suppose that the magnetic momentum at rest is null
    t00+=generateFIDsignal(tmpX, tmpY, tmpZ, dt, dw);
    if(!rotFrame) frameToFix(t00);
    MUsecondFID1=MU;
    //MU.save("spinEcho.txt");
    //T.save("spinEchoTime.txt");

    //restore the parameters
    Bz0=tmpBz0;
    Bxy=tmpBxy; omega=gamma*Bz0;
    muz0=tempMuz0;




    /*****************************************************************/
    /** Then we consider the simulation with a slightly smaller field*/
    /******************************************************/
    /** initialization */
    tmpX=mx0; tmpY=my0; tmpZ=mz0;
    dw=-dw;
    muz0 = tempMuz0*(1.0+(dw/(gamma*tmpBz0)));
    std::cout << 1.0+(dw/(gamma*tmpBz0)) << " " << tempMuz0 << " " << muz0 << std::endl;
    //apply the constant magnetic field during t0
    Bz0=tmpBz0+(dw/gamma);
    Bxy=0.0; omega=0;
    //double t0=MAX(T1,T2);
    sigLength=200;
    dt=t0/((double) sigLength);
    t00=0.0;

    //suppose that the magnetic momentum at rest is null
    t00+=generateFIDsignal(tmpX, tmpY, tmpZ, dt, dw);
    tmpX=MU(MU.endL,0); tmpY=MU(MU.endL,1); tmpZ=MU(MU.endL,2);
    if(!rotFrame) frameToFix(t00);
    MU=(MU+MUsteady1)*0.5;
    MU.save("initSig.txt");
    T.save("initSigTime.txt");
    std::cout << t00 << std::endl;



    /******************************************************/
    /** first pulse */
    //play one pi/2-pulse
    Bz0=tmpBz0+(dw/gamma);
    Bxy=tmpBxy; omega=gamma*Bz0;
    sigLength=4000;
    //t90 = Pi/(2.0*gamma*Bxy);
    dt=tw/((double) sigLength);

    //suppose that the magnetic momentum at rest is null
    t00+=generateRMNsignal(tmpX, tmpY, tmpZ, dt, dw);
    tmpX=MU(MU.endL,0); tmpY=MU(MU.endL,1); tmpZ=MU(MU.endL,2);
    if(!rotFrame) frameToFix(t00);
    std::cout << MU.endL << " " << MUfisrtPulse1.endL << std::endl;
    MU=(MU+MUfisrtPulse1)*0.5;
    MU.save("p90.txt");
    T.save("p90Time.txt");
    std::cout << t00 << std::endl;




    /******************************************************/
    /** free induction decay */
    //stop the radiofrequence pulse until the next pulse
    Bz0=tmpBz0+(dw/gamma);
    Bxy=0.0; omega=0;
    sigLength=4000;
    dt=(0.5*TE - tw)/((double) sigLength);
    //suppose that the magnetic momentum at rest is null
    t00+=generateFIDsignal(tmpX, tmpY, tmpZ, dt, dw);
    tmpX=MU(MU.endL,0); tmpY=MU(MU.endL,1); tmpZ=MU(MU.endL,2);
    if(!rotFrame) frameToFix(t00);
    std::cout << MU.endL << " " << MUfirstFID1.endL << std::endl;
    MU=(MU+MUfirstFID1)*0.5;
    MU.save("fid.txt");
    T.save("fidTime.txt");




    /******************************************************/
    /** second pulse */
    //play the pi-pulse
    Bz0=tmpBz0+(dw/gamma);
    Bxy=tmpBxy; omega=gamma*Bz0;
    sigLength=4000;
    dt=2.0*tw/((double) sigLength);

    //suppose that the magnetic momentum at rest is null
    t00+=generateRMNsignal(tmpX, tmpY, tmpZ, dt, dw);
    tmpX=MU(MU.endL,0); tmpY=MU(MU.endL,1); tmpZ=MU(MU.endL,2);
    if(!rotFrame) frameToFix(t00);
    std::cout << MU.endL << " " << MUsecondPulse1.endL << std::endl;
    MU=(MU+MUsecondPulse1)*0.5;
    MU.save("p180.txt");
    T.save("p180Time.txt");
    std::cout << t00 << std::endl;



    /******************************************************/
    /** spin echo */
    //stop the radiofrequence pulse until the next pulse
    Bz0=tmpBz0+(dw/gamma);
    Bxy=0.0; omega=0;
    sigLength=2000;
    dt=(TE-tw)/((double) sigLength);

    //suppose that the magnetic momentum at rest is null
    t00+=generateFIDsignal(tmpX, tmpY, tmpZ, dt, dw);
    if(!rotFrame) frameToFix(t00);
    std::cout << MU.endL << " " << MUsecondFID1.endL << std::endl;
    MU=(MU+MUsecondFID1)*0.5;
    MU.save("spinEcho.txt");
    T.save("spinEchoTime.txt");

    //restore the parameters
    Bz0=tmpBz0;
    Bxy=tmpBxy; omega=gamma*Bz0;
    muz0=tempMuz0;


    return;
}

void C_el_hahn::multiChrome(double mx0, double my0, double mz0, int N, bool rotFrame)
{
    C_matrix<double> MUsteady1;
    C_matrix<double> MUfisrtPulse1;
    C_matrix<double> MUfirstFID1;
    C_matrix<double> MUsecondPulse1;
    C_matrix<double> MUsecondFID1;
    //the perturbating magnetic field pulsation is gamma*Bz0, but the actual constant field is Bz0+(dw/gamma)
    double dt;
    double tmpBz0 = Bz0;
    double tmpBxy = Bxy;
    double tempMuz0=muz0;
    double dw=0.0;
    double tmpX=mx0, tmpY=my0, tmpZ=mz0;

    muz0 = tempMuz0*(1.0+(dw/(gamma*tmpBz0)));
    //std::cout << 1.0+(dw/(gamma*tmpBz0)) << " " << tempMuz0 << " " << muz0 << std::endl;

    /*****************************************************************/
    /** We first consider the simulation with a slightly larger field*/
    /******************************************************/
    /** initialization */
    //apply the constant magnetic field during t0
    Bz0=tmpBz0+(dw/gamma);
    Bxy=0.0; omega=0;
    double t0=MAX(T1,T2);
    sigLength=200;
    dt=t0/((double) sigLength);
    double t00=0.0;

    //suppose that the magnetic momentum at rest is null
    t00+=generateFIDsignal(tmpX, tmpY, tmpZ, dt, dw);
    tmpX=MU(MU.endL,0); tmpY=MU(MU.endL,1); tmpZ=MU(MU.endL,2);
    if(!rotFrame) frameToFix(t00);
    MUsteady1=MU;

    T.save("initSigTime.txt");


    /******************************************************/
    /** first pulse */
    //play one pi/2-pulse
    Bz0=tmpBz0+(dw/gamma);
    Bxy=tmpBxy; omega=gamma*Bz0;
    sigLength=4000;
    dt=tw/((double) sigLength);

    //suppose that the magnetic momentum at rest is null
    t00+=generateRMNsignal(tmpX, tmpY, tmpZ, dt, dw);
    tmpX=MU(MU.endL,0); tmpY=MU(MU.endL,1); tmpZ=MU(MU.endL,2);
    if(!rotFrame) frameToFix(t00);
    MUfisrtPulse1=MU;

    T.save("p90Time.txt");


    /******************************************************/
    /** free induction decay */
    //stop the radiofrequence pulse until the next pulse
    Bz0=tmpBz0+(dw/gamma);
    Bxy=0.0; omega=0;
    sigLength=4000;
    dt=(0.5*TE - tw)/((double) sigLength);

    //suppose that the magnetic momentum at rest is null
    t00+=generateFIDsignal(tmpX, tmpY, tmpZ, dt, dw);
    tmpX=MU(MU.endL,0); tmpY=MU(MU.endL,1); tmpZ=MU(MU.endL,2);
    if(!rotFrame) frameToFix(t00);
    MUfirstFID1=MU;

    T.save("fidTime.txt");



    /******************************************************/
    /** second pulse */
    //play the pi-pulse
    Bz0=tmpBz0+(dw/gamma);
    Bxy=tmpBxy; omega=gamma*Bz0;
    sigLength=4000;
    dt=tw/((double) sigLength);

    //suppose that the magnetic momentum at rest is null
    t00+=generateRMNsignal(tmpX, tmpY, tmpZ, dt, dw);
    tmpX=MU(MU.endL,0); tmpY=MU(MU.endL,1); tmpZ=MU(MU.endL,2);
    if(!rotFrame) frameToFix(t00);
    MUsecondPulse1=MU;

    T.save("p180Time.txt");


    /******************************************************/
    /** spin echo */
    //stop the radiofrequence pulse until the next pulse
    Bz0=tmpBz0+(dw/gamma);
    Bxy=0.0; omega=0;
    sigLength=4000;
    dt=1.5*(TE-tw)/((double) sigLength);

    //suppose that the magnetic momentum at rest is null
    t00+=generateFIDsignal(tmpX, tmpY, tmpZ, dt, dw);
    if(!rotFrame) frameToFix(t00);
    MUsecondFID1=MU;

    T.save("spinEchoTime.txt");

    //restore parameters
    Bz0=tmpBz0;
    Bxy=tmpBxy; omega=gamma*Bz0;
    muz0=tempMuz0;
    tmpX=mx0; tmpY=my0; tmpZ=mz0;



    for(int n=0 ; n<N ; n++)
    {
        t00=0.0;
        dw=-DW+2.0*DW*((double) n)/((double) (N-1));//
        //dw=2.0*DW*((((double) rand())/((double) RAND_MAX))-0.5);
        muz0 = tempMuz0*(1.0+(dw/(gamma*tmpBz0)));

        /*****************************************************************/
        /** We first consider the simulation with a slightly larger field*/
        /******************************************************/
        /** initialization */
        //apply the constant magnetic field during t0
        Bz0=tmpBz0+(dw/gamma);
        Bxy=0.0; omega=0;
        sigLength=200;
        dt=t0/((double) sigLength);

        //suppose that the magnetic momentum at rest is null
        t00+=generateFIDsignal(tmpX, tmpY, tmpZ, dt, dw);
        tmpX=MU(MU.endL,0); tmpY=MU(MU.endL,1); tmpZ=MU(MU.endL,2);
        if(!rotFrame) frameToFix(t00);
        MUsteady1=(MU+MUsteady1);
        //std::cout << t00 << std::endl;


        /******************************************************/
        /** first pulse */
        //play one pi/2-pulse
        Bz0=tmpBz0+(dw/gamma);
        Bxy=tmpBxy; omega=gamma*Bz0;
        sigLength=4000;
        dt=tw/((double) sigLength);

        //suppose that the magnetic momentum at rest is null
        t00+=generateRMNsignal(tmpX, tmpY, tmpZ, dt, dw);
        tmpX=MU(MU.endL,0); tmpY=MU(MU.endL,1); tmpZ=MU(MU.endL,2);
        if(!rotFrame) frameToFix(t00);
        MUfisrtPulse1=(MU+MUfisrtPulse1);
        //std::cout << t00 << std::endl;


        /******************************************************/
        /** free induction decay */
        //stop the radiofrequence pulse until the next pulse
        Bz0=tmpBz0+(dw/gamma);
        Bxy=0.0; omega=0;
        sigLength=4000;
        dt=(0.5*TE - tw)/((double) sigLength);

        //suppose that the magnetic momentum at rest is null
        t00+=generateFIDsignal(tmpX, tmpY, tmpZ, dt, dw);
        tmpX=MU(MU.endL,0); tmpY=MU(MU.endL,1); tmpZ=MU(MU.endL,2);
        if(!rotFrame) frameToFix(t00);
        MUfirstFID1=(MU+MUfirstFID1);



        /******************************************************/
        /** second pulse */
        //play the pi-pulse
        Bz0=tmpBz0+(dw/gamma);
        Bxy=tmpBxy; omega=gamma*Bz0;
        sigLength=4000;
        dt=tw/((double) sigLength);//2.0*

        //suppose that the magnetic momentum at rest is null
        t00+=generateRMNsignal(tmpX, tmpY, tmpZ, dt, dw);
        tmpX=MU(MU.endL,0); tmpY=MU(MU.endL,1); tmpZ=MU(MU.endL,2);
        if(!rotFrame) frameToFix(t00);
        MUsecondPulse1=(MU+MUsecondPulse1);
        //std::cout << t00 << std::endl;


        /******************************************************/
        /** spin echo */
        //stop the radiofrequence pulse until the next pulse
        Bz0=tmpBz0+(dw/gamma);
        Bxy=0.0; omega=0;
        sigLength=4000;
        dt=1.5*(TE-tw)/((double) sigLength);

        //suppose that the magnetic momentum at rest is null
        t00+=generateFIDsignal(tmpX, tmpY, tmpZ, dt, dw);
        if(!rotFrame) frameToFix(t00);
        MUsecondFID1=(MU+MUsecondFID1);

        //restore the parameters
        Bz0=tmpBz0;
        Bxy=tmpBxy; omega=gamma*Bz0;
        muz0=tempMuz0;
        tmpX=mx0; tmpY=my0; tmpZ=mz0;
    }

    MUsteady1=MUsteady1*(1.0/((double) (N+1)));
    MUfisrtPulse1=MUfisrtPulse1*(1.0/((double) (N+1)));
    MUfirstFID1=MUfirstFID1*(1.0/((double) (N+1)));
    MUsecondPulse1=MUsecondPulse1*(1.0/((double) (N+1)));
    MUsecondFID1=MUsecondFID1*(1.0/((double) (N+1)));
    MUsteady1.save("initSig.txt");
    MUfisrtPulse1.save("p90.txt");
    MUfirstFID1.save("fid.txt");
    MUsecondPulse1.save("p180.txt");
    MUsecondFID1.save("spinEcho.txt");

    return;
}

void C_el_hahn::monochromaticDiffusion(double mx0, double my0, double mz0, double dw, bool rotFrame)
{
    //the perturbating magnetic field pulsation is gamma*Bz0, but the actual constant field is Bz0+(dw/gamma)
    double dt;
    double tmpBz0 = Bz0;
    double tmpBxy = Bxy;
    double tmpX=mx0, tmpY=my0, tmpZ=mz0;
    double tempMuz0=muz0;
    muz0 = tempMuz0*(1.0+(dw/(gamma*tmpBz0)));

    /******************************************************/
    /** initialization */
    //apply the constant magnetic field during t0
    Bz0=tmpBz0+(dw/gamma);
    Bxy=0.0; omega=0;
    double t0=MAX(T1,T2);
    sigLength=200;
    dt=t0/((double) sigLength);
    double t00=0.0;

    //suppose that the magnetic momentum at rest is null
    t00+=generateFIDsignalWithDiffusion(tmpX, tmpY, tmpZ, dt, dw);
    tmpX=MU(MU.endL,0); tmpY=MU(MU.endL,1); tmpZ=MU(MU.endL,2);
    if(!rotFrame) frameToFix(t00);
    MU.save("initSig.txt");
    T.save("initSigTime.txt");
    std::cout << t00 << std::endl;


    /******************************************************/
    /** first pulse */
    //play one pi/2-pulse
    Bz0=tmpBz0+(dw/gamma);
    Bxy=tmpBxy; omega=gamma*Bz0;
    sigLength=4000;
    //t90 = Pi/(2.0*gamma*Bxy);
    dt=tw/((double) sigLength);

    //suppose that the magnetic momentum at rest is null
    t00+=generateRMNsignal(tmpX, tmpY, tmpZ, dt, dw);
    tmpX=MU(MU.endL,0); tmpY=MU(MU.endL,1); tmpZ=MU(MU.endL,2);
    if(!rotFrame) frameToFix(t00);
    MU.save("p90.txt");
    T.save("p90Time.txt");
    std::cout << t00 << std::endl;


    /******************************************************/
    /** free induction decay */
    //stop the radiofrequence pulse until the next pulse
    Bz0=tmpBz0+(dw/gamma);
    Bxy=0.0; omega=0;
    sigLength=4000;
    dt=(0.5*TE - tw)/((double) sigLength);
    //suppose that the magnetic momentum at rest is null
    t00+=generateFIDsignalWithDiffusion(tmpX, tmpY, tmpZ, dt, dw);
    tmpX=MU(MU.endL,0); tmpY=MU(MU.endL,1); tmpZ=MU(MU.endL,2);
    if(!rotFrame) frameToFix(t00);
    MU.save("fid.txt");
    T.save("fidTime.txt");



    /******************************************************/
    /** second pulse */
    //play the pi-pulse
    Bz0=tmpBz0+(dw/gamma);
    Bxy=tmpBxy; omega=gamma*Bz0;
    sigLength=4000;
    dt=2.0*tw/((double) sigLength);

    //suppose that the magnetic momentum at rest is null
    t00+=generateRMNsignal(tmpX, tmpY, tmpZ, dt, dw);
    tmpX=MU(MU.endL,0); tmpY=MU(MU.endL,1); tmpZ=MU(MU.endL,2);
    if(!rotFrame) frameToFix(t00);
    MU.save("p180.txt");
    T.save("p180Time.txt");
    std::cout << t00 << std::endl;


    /******************************************************/
    /** spin echo */
    //stop the radiofrequence pulse until the next pulse
    Bz0=tmpBz0+(dw/gamma);
    Bxy=0.0; omega=0;
    sigLength=2000;
    dt=(TE-tw)/((double) sigLength);

    //suppose that the magnetic momentum at rest is null
    t00+=generateFIDsignalWithDiffusion(tmpX, tmpY, tmpZ, dt, dw);
    if(!rotFrame) frameToFix(t00);
    MU.save("spinEcho.txt");
    T.save("spinEchoTime.txt");

    //restore the parameters
    Bz0=tmpBz0;
    Bxy=tmpBxy; omega=gamma*Bz0;
    muz0=tempMuz0;
    return;
    return;
}

void C_el_hahn::multiChromeDiffusion(double mx0, double my0, double mz0, int N, bool rotFrame)
{
    C_matrix<double> MUsteady1;
    C_matrix<double> MUfisrtPulse1;
    C_matrix<double> MUfirstFID1;
    C_matrix<double> MUsecondPulse1;
    C_matrix<double> MUsecondFID1;
    //the perturbating magnetic field pulsation is gamma*Bz0, but the actual constant field is Bz0+(dw/gamma)
    double dt;
    double tmpBz0 = Bz0;
    double tmpBxy = Bxy;
    double tempMuz0=muz0;
    double dw=0.0;
    double tmpX=mx0, tmpY=my0, tmpZ=mz0;

    muz0 = tempMuz0*(1.0+(dw/(gamma*tmpBz0)));
    //std::cout << 1.0+(dw/(gamma*tmpBz0)) << " " << tempMuz0 << " " << muz0 << std::endl;

    /*****************************************************************/
    /** We first consider the simulation with a slightly larger field*/
    /******************************************************/
    /** initialization */
    //apply the constant magnetic field during t0
    Bz0=tmpBz0+(dw/gamma);
    Bxy=0.0; omega=0;
    double t0=MAX(T1,T2);
    sigLength=200;
    dt=t0/((double) sigLength);
    double t00=0.0;

    //suppose that the magnetic momentum at rest is null
    t00+=generateFIDsignalWithDiffusion(tmpX, tmpY, tmpZ, dt, dw);
    tmpX=MU(MU.endL,0); tmpY=MU(MU.endL,1); tmpZ=MU(MU.endL,2);
    if(!rotFrame) frameToFix(t00);
    MUsteady1=MU;

    T.save("initSigTime.txt");


    /******************************************************/
    /** first pulse */
    //play one pi/2-pulse
    Bz0=tmpBz0+(dw/gamma);
    Bxy=tmpBxy; omega=gamma*Bz0;
    sigLength=4000;
    dt=tw/((double) sigLength);

    //suppose that the magnetic momentum at rest is null
    t00+=generateRMNsignal(tmpX, tmpY, tmpZ, dt, dw);
    tmpX=MU(MU.endL,0); tmpY=MU(MU.endL,1); tmpZ=MU(MU.endL,2);
    if(!rotFrame) frameToFix(t00);
    MUfisrtPulse1=MU;

    T.save("p90Time.txt");


    /******************************************************/
    /** free induction decay */
    //stop the radiofrequence pulse until the next pulse
    Bz0=tmpBz0+(dw/gamma);
    Bxy=0.0; omega=0;
    sigLength=4000;
    dt=(0.5*TE - tw)/((double) sigLength);

    //suppose that the magnetic momentum at rest is null
    t00+=generateFIDsignalWithDiffusion(tmpX, tmpY, tmpZ, dt, dw);
    tmpX=MU(MU.endL,0); tmpY=MU(MU.endL,1); tmpZ=MU(MU.endL,2);
    if(!rotFrame) frameToFix(t00);
    MUfirstFID1=MU;

    T.save("fidTime.txt");



    /******************************************************/
    /** second pulse */
    //play the pi-pulse
    Bz0=tmpBz0+(dw/gamma);
    Bxy=tmpBxy; omega=gamma*Bz0;
    sigLength=4000;
    dt=tw/((double) sigLength);

    //suppose that the magnetic momentum at rest is null
    t00+=generateRMNsignal(tmpX, tmpY, tmpZ, dt, dw);
    tmpX=MU(MU.endL,0); tmpY=MU(MU.endL,1); tmpZ=MU(MU.endL,2);
    if(!rotFrame) frameToFix(t00);
    MUsecondPulse1=MU;

    T.save("p180Time.txt");


    /******************************************************/
    /** spin echo */
    //stop the radiofrequence pulse until the next pulse
    Bz0=tmpBz0+(dw/gamma);
    Bxy=0.0; omega=0;
    sigLength=4000;
    dt=1.5*(TE-tw)/((double) sigLength);

    //suppose that the magnetic momentum at rest is null
    t00+=generateFIDsignalWithDiffusion(tmpX, tmpY, tmpZ, dt, dw);
    if(!rotFrame) frameToFix(t00);
    MUsecondFID1=MU;

    T.save("spinEchoTime.txt");

    //restore parameters
    Bz0=tmpBz0;
    Bxy=tmpBxy; omega=gamma*Bz0;
    muz0=tempMuz0;
    tmpX=mx0; tmpY=my0; tmpZ=mz0;



    for(int n=0 ; n<N ; n++)
    {
        t00=0.0;
        dw=-DW+2.0*DW*((double) n)/((double) (N-1));//
        //dw=2.0*DW*((((double) rand())/((double) RAND_MAX))-0.5);
        muz0 = tempMuz0*(1.0+(dw/(gamma*tmpBz0)));

        /*****************************************************************/
        /** We first consider the simulation with a slightly larger field*/
        /******************************************************/
        /** initialization */
        //apply the constant magnetic field during t0
        Bz0=tmpBz0+(dw/gamma);
        Bxy=0.0; omega=0;
        sigLength=200;
        dt=t0/((double) sigLength);

        //suppose that the magnetic momentum at rest is null
        t00+=generateFIDsignalWithDiffusion(tmpX, tmpY, tmpZ, dt, dw);
        tmpX=MU(MU.endL,0); tmpY=MU(MU.endL,1); tmpZ=MU(MU.endL,2);
        if(!rotFrame) frameToFix(t00);
        MUsteady1=(MU+MUsteady1);
        //std::cout << t00 << std::endl;


        /******************************************************/
        /** first pulse */
        //play one pi/2-pulse
        Bz0=tmpBz0+(dw/gamma);
        Bxy=tmpBxy; omega=gamma*Bz0;
        sigLength=4000;
        dt=tw/((double) sigLength);

        //suppose that the magnetic momentum at rest is null
        t00+=generateRMNsignal(tmpX, tmpY, tmpZ, dt, dw);
        tmpX=MU(MU.endL,0); tmpY=MU(MU.endL,1); tmpZ=MU(MU.endL,2);
        if(!rotFrame) frameToFix(t00);
        MUfisrtPulse1=(MU+MUfisrtPulse1);
        //std::cout << t00 << std::endl;


        /******************************************************/
        /** free induction decay */
        //stop the radiofrequence pulse until the next pulse
        Bz0=tmpBz0+(dw/gamma);
        Bxy=0.0; omega=0;
        sigLength=4000;
        dt=(0.5*TE - tw)/((double) sigLength);

        //suppose that the magnetic momentum at rest is null
        t00+=generateFIDsignalWithDiffusion(tmpX, tmpY, tmpZ, dt, dw);
        tmpX=MU(MU.endL,0); tmpY=MU(MU.endL,1); tmpZ=MU(MU.endL,2);
        if(!rotFrame) frameToFix(t00);
        MUfirstFID1=(MU+MUfirstFID1);



        /******************************************************/
        /** second pulse */
        //play the pi-pulse
        Bz0=tmpBz0+(dw/gamma);
        Bxy=tmpBxy; omega=gamma*Bz0;
        sigLength=4000;
        dt=tw/((double) sigLength);//2.0*

        //suppose that the magnetic momentum at rest is null
        t00+=generateRMNsignal(tmpX, tmpY, tmpZ, dt, dw);
        tmpX=MU(MU.endL,0); tmpY=MU(MU.endL,1); tmpZ=MU(MU.endL,2);
        if(!rotFrame) frameToFix(t00);
        MUsecondPulse1=(MU+MUsecondPulse1);
        //std::cout << t00 << std::endl;


        /******************************************************/
        /** spin echo */
        //stop the radiofrequence pulse until the next pulse
        Bz0=tmpBz0+(dw/gamma);
        Bxy=0.0; omega=0;
        sigLength=4000;
        dt=1.5*(TE-tw)/((double) sigLength);

        //suppose that the magnetic momentum at rest is null
        t00+=generateFIDsignalWithDiffusion(tmpX, tmpY, tmpZ, dt, dw);
        if(!rotFrame) frameToFix(t00);
        MUsecondFID1=(MU+MUsecondFID1);

        //restore the parameters
        Bz0=tmpBz0;
        Bxy=tmpBxy; omega=gamma*Bz0;
        muz0=tempMuz0;
        tmpX=mx0; tmpY=my0; tmpZ=mz0;
    }

    MUsteady1=MUsteady1*(1.0/((double) (N+1)));
    MUfisrtPulse1=MUfisrtPulse1*(1.0/((double) (N+1)));
    MUfirstFID1=MUfirstFID1*(1.0/((double) (N+1)));
    MUsecondPulse1=MUsecondPulse1*(1.0/((double) (N+1)));
    MUsecondFID1=MUsecondFID1*(1.0/((double) (N+1)));
    MUsteady1.save("initSig.txt");
    MUfisrtPulse1.save("p90.txt");
    MUfirstFID1.save("fid.txt");
    MUsecondPulse1.save("p180.txt");
    MUsecondFID1.save("spinEcho.txt");
    return;
}
