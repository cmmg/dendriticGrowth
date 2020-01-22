//problem geometry, mesh control
#define DIMS 2
#define problemWidth 9.0
#define refinementFactor 8

//phase field properties
#define InterfaceEnergyParameter {1.0e-2, 1.0e-2, 1.0e-2} //{Kx, Ky, Kz}
#define dFdC  400*c[q]*(c[q]-1.0)*(c[q]-0.5) //derivative of the free energy
#define Mobility 1.0

//time step controls
#define TimeStep 0.0002
#define TotalTime 10000*TimeStep

//output controls
#define outputFileName "solution"

//solidifcation kobayashi parameters 
#define ebar 0.01
#define delta 0.040
#define mm 6.0
#define theta0 0 //3.1416/2.0
#define L 3333.33
#define tau0 1.0/L   // (1/L) where L = 3333.33
#define KK 2.0
#define amplitude 0.01
