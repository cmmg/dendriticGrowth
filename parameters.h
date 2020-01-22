//problem geometry, mesh control
#define DIMS 2
#define problemWidth 500.0
#define refinementFactor 8

//phase field properties
#define InterfaceEnergyParameter {1.0e-2, 1.0e-2, 1.0e-2} //{Kx, Ky, Kz}
#define dFdC  400*c[q]*(c[q]-1.0)*(c[q]-0.5) //derivative of the free energy
#define Mobility 1.0

//time step controls
#define TimeStep 0.04
#define TotalTime 35000*TimeStep
#define PSTEPS 1000 

//output controls
#define outputFileName "solution"

//solidifcation kobayashi parameters 

#define D 1.0

#define tau0 1.0
#define W0 1.0
#define mm 4.0
#define em 0.05
#define theta0 0.125 
#define lam  D*tau0/W0/W0/0.6267

//#define L 3333.33
//#define tau0 1.0/L   // (1/L) where L = 3333.33
//#define KK 2.0
//#define amplitude 0.01
