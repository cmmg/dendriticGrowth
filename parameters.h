//problem geometry, mesh control
#define DIMS 2
#define problem_Height 615
#define problem_Width 615
#define refinementFactor 6

//phase field properties
#define InterfaceEnergyParameter {1.0e-2, 1.0e-2, 1.0e-2} //{Kx, Ky, Kz}
#define dFdC  400*c[q]*(c[q]-1.0)*(c[q]-0.5) //derivative of the free energy
#define Mobility 1.0

//time step controls
#define TimeStep 5.0
#define TotalTime 15000*TimeStep

//output controls
#define outputFileName "solution"


//soldification zhu
#define aone 0.8839
#define atwo  0.6267
#define zhi 10.5
#define lam aone*zhi  // 0.8839*54 = 47.73

//#define D_tilda 30.0
//#define thermalby1_k 28451.0

#define thermal 14794
#define VV 0.515
#define ep 0.01
#define mm 4.0
#define theta0 0.0
#define ke 0.48

#define D 5.82

#define DS 1.0
#define DL 3000.0

//randomness parameters
#define OMEGA1 0.3
#define OMEGA2 0 //0.1 

