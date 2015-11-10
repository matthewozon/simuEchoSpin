#ifndef C_SIGRMN_H
#define C_SIGRMN_H
#include "C_matrix.h"

class C_sigRMN
{
public:
    C_sigRMN();

    //the physical parameters
    double muz0;//initial state of the magnetization or external resultant
    double T1, T2;//caracteristic time response of the tissue
    double Bz0;//magnet field: the constant componant
    double Bxy, omega;//amplitude, pulsation of the RF wave (it can be either in the rotating or lab frame)
    /** **************************************/
    /** the total magnetic field applied is  */
    /** Bx = Bx*cos(omega*t)                 */
    /** By = By*sin(omega*t)                 */
    /** Bz = Bz0                             */
    /** rotating frame R=omega z             */
    /** **************************************/
    double gamma;//gyromagnetic monentum of water
    double tw;//duration of the rf pulse (90 degrees) it is assumed that tw<< T1, T2, TE
    double TE;//the echo time. TE>> 1/DW
    double DW;//the bandwidth of the Larmor frequencies. DW<<gamma*sqrt(Bx^2+By^2)=gamma*Bxy
    double D;//the diffusion coeficient
    double G;//the gradient of the local magnetic field

    //the simulated signal
    C_matrix<double> MU;//the magnetic momentum
    C_matrix<double> T;//for each signal generation, the time start at 0.0
    int sigLength;//length of the signal

    void frameToFix(double t00);//compute the change of reference frame: from the rotating system to the lab frame
    void frameToRotating(double t00);


    //generation of the signals in the rotating frame
    /** monochromatic resonnance without damping or diffusion */
    /** *******************************************************/
    /** solve the equation system is the rotating frame:      */
    /** d(mux)/dt = -dw muy                                   */
    /** d(muy)/dt =  dw mux - omega1 muz                      */
    /** d(muz)/dt =  omega1 muy                               */
    /** where dw=omega-omega0, omega1=gamma*Bxy, omega0 the   */
    /**                                 mean Larmor frequency */
    /** *******************************************************/
    double generateRMNsignal(double mx0, double my0, double mz0, double dt, double dw /** pulsation shift due to variation in the magnetic field Bz0*/);

    /** monochromatic free induction decay with damping but no diffusion*/
    /** *******************************************************/
    /** solve the equation system is the rotating frame:      */
    /** d(mux)/dt = -mux/T2 -dw muy                           */
    /** d(muy)/dt =  dw mux -muy/T2                           */
    /** d(muz)/dt =  -(muz-muz0)/T1                           */
    /** where dw=omega-omega0, omega1=gamma*Bxy, omega0 the   */
    /**                                 mean Larmor frequency */
    /** *******************************************************/
    double generateFIDsignal(double mx0, double my0, double mz0, double dt, double dw /** pulsation shift due to variation in the magnetic field Bz0*/);


    /** monochromatic free induction decay with damping but no diffusion*/
    /** *******************************************************/
    /** solve the equation system is the rotating frame:      */
    /** d(mux)/dt = -mux/T2 -(dw+delta(t)) muy                           */
    /** d(muy)/dt =  (dw+delta(t)) mux -muy/T2                           */
    /** d(muz)/dt =  -(muz-muz0)/T1                           */
    /** where dw=omega-omega0, omega1=gamma*Bxy, omega0 the   */
    /**                                 mean Larmor frequency */
    /** *******************************************************/
    double generateFIDsignalWithDiffusion(double mx0, double my0, double mz0, double dt, double dw /** pulsation shift due to variation in the magnetic field Bz0*/);



};

#endif // C_SIGRMN_H
