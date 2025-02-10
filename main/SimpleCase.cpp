#include "TimeEvolution.h"
#include <cstring>

/*
TimeEvolutionF:            runs usual RK4 time evolution
PointConvergenceEvolution: runs pointwise convergence evolution. Poisson convergence is calculated here as well
PointConvergence:          runs pointwise convergence at the end of the simulation
PartialPointConvergence:   runs pointwise convergence for each variable, at the end of the simulation
*/

//For multithreading when writing
//Define a mutex to synchronize file access
std::mutex fileMutex;

double gaussian(double x) {
    double sigma = 0.2;
    double mu2 = 0.25;
    return 0.01*x*x*exp(-(x*x - mu2*mu2) * (x*x - mu2*mu2) / (4. * sigma * sigma * sigma * sigma));
}

double gaussianA3r(double x) {
    double sigma = 0.2;
    double mu2 = 0.25;
    return 0.1*x*x*exp(-(x*x - mu2*mu2) * (x*x - mu2*mu2) / (4. * sigma * sigma * sigma * sigma));
}

double gaussianCphi(double x) {
    double sigma = 0.2;
    double mu2 = 0.4;
    return 1.*x*x*exp(-(x*x - mu2*mu2) * (x*x - mu2*mu2) / (4. * sigma * sigma * sigma * sigma));
}

int main(int argc, char* argv[]){

    std::string Choice = "te"; //te, pc, l2. This chooses what type of run you want. Check the Makefile if you need

    // Check if there are command line arguments -> This is what we will use to run the different cases
    if(argc <= 1) std::cout << "No command line arguments provided, using default value for Time Evolution" << std::endl;
    else          Choice = std::string(argv[1]);

    int    variables(16)  ; //Number of variables; Take care with this number
    int    sizex    (3006); //9000, 3000, 1000
    int    sizet    (108*StepInstantWrite); //Ensure that this is always an integer multiple of StepInstantWrite
    double dx       (rscri/sizex);
    double dt       (0.00001)     ; //0.000018, 0.00006, 0.00002

    std::cout << "CFL condition: " << dt/dx << std::endl;

    //Right Hand Sides
    void (*RHS[])(double**, double**, double*, int, double, int*, double*) = {ElectricField, PsiElectric, Beta, Alpha, Phi, A, cphi, dphi, cPi, dPi, trK, Gammarr, Chi, Arr, Lambdar, Theta};

    std::cout << "RHS Set" << std::endl;

    ////////////////////////////////////////
    ///////////INITIAL CONDITIONS///////////
    ////////////////////////////////////////

    double*  x     = new double [sizex];
    double** data  = new double*[variables];
    double** ddata = new double*[variables];

    for(int j = 0; j < sizex; j++)          x[j] = dx/2. + j*dx;
    for(int j = 0; j < variables; j++){  data[j] = new double[sizex]();
                                        ddata[j] = new double[sizex]();}

    ////////////////////////////////////////
    /////////////RETRIEVING CHI/////////////
    ////////////////////////////////////////

    // Read data from file which has the output of the Mathematica notebook
    std::string chi_file;
    
    chi_file = "/home/machado/Desktop/Universidade/Tese/thecode/output/chiA3r_values.txt";

    std::ifstream inputFile(chi_file);
    if (!inputFile) {
        std::cerr << "Failed to open file." << std::endl;
        return 1;
    }

    //Array for the conformal factor. The values are retrieved from 
    //mathematica notebook RNinidata.nb, because there is no analytical
    //solution

    int RawSize = 20000;

    double* RawChi = new double[RawSize];
    double* RawX   = new double[RawSize];

    //Retrieving from the file which is the output of Mathematica
    double xMATH, yMATH;
    int ctr = 0;
    while (inputFile >> xMATH >> yMATH) {
        RawChi[ctr] = yMATH;
        RawX  [ctr] = xMATH;
        ctr++;
    }

    inputFile.close();

    SpaceDerivator sder(Ghost_Left, Ghost_Right);

    double* chi = sder.FourthOrderInterpolator(RawChi, RawSize, RawX, 0, x, sizex);

    double* EI       = data[Eindex        ];
    double* PsiI     = data[Psiindex      ];
    double* BetaI    = data[Betaindex     ];
    double* AlphaI   = data[Alphaindex    ];
    double* PhiI     = data[Phiindex      ];
    double* AI       = data[Aindex        ];
    double* cphiI    = data[cphiindex     ];
    double* dphiI    = data[dphiindex     ];
    double* cPiI     = data[cPiindex      ];
    double* dPiI     = data[dPiindex      ];
    double* trKI     = data[trKindex      ];
    double* GammarrI = data[Gammarrindex  ];
    double* ChiI     = data[Chiindex      ];
    double* ArrI     = data[Arrindex      ];
    double* LambdarI = data[Lambdarindex  ];
    double* ThetaI   = data[Thetaindex    ];

    //Extra just for initial conditions
    double* RHSArI   = new double[sizex];

    //Parity conditions: vectors are odd, scalars are even
                               //E, Psi, Beta, Alpha, Phi, A, cphi, dphi, cPi, dPi, trK, Gammarr, Chi, Arr, Lambdar, Theta
    int* p = new int[variables]{-1,   1,   -1,     1,   1,-1,    1,    1,   1,   1,   1,       1,   1,  -1,      -1,     1};

    double* Omegaf  = Omega(data, ddata, x, sizex, 0., p);
    double* dOmega  = sder.FirstDerSpaceCenteredDiff4(Omegaf, sizex, x, 1.);

    for(int i = 0; i < sizex; i++){
        PsiI     [i] = EM*(0.0);
        BetaI    [i] = Kcmc*x[i]/3.;
        AlphaI   [i] = sqrt(Omegaf[i]*Omegaf[i]+(Kcmc*Kcmc*x[i]*x[i])/9.);
        AI       [i] = EM*(gaussianA3r(x[i]));
        PhiI     [i] = EM*0.0;
        cphiI    [i] = SF*(gaussianCphi(x[i]));
        dphiI    [i] = SF*0.0;
        trKI     [i] = 0.;
        GammarrI [i] = 1.;
        ChiI     [i] = chi[i];
        //ChiI     [i] = 1.;
        ArrI     [i] = 0.;
        LambdarI [i] = 0.;
        ThetaI   [i] = 0.;
        RHSArI   [i] = 0.;
    }

    double* dcphi    = sder.FirstDerSpaceCenteredDiff4(cphiI   , sizex, x, p[cphiindex   ]);
    double* ddphi    = sder.FirstDerSpaceCenteredDiff4(dphiI   , sizex, x, p[dphiindex   ]);
    double* dChi     = sder.FirstDerSpaceCenteredDiff4(ChiI    , sizex, x, p[Chiindex    ]);
    double* dAlpha   = sder.FirstDerSpaceCenteredDiff4(AlphaI  , sizex, x, p[Alphaindex  ]);
    double* dGammarr = sder.FirstDerSpaceCenteredDiff4(GammarrI, sizex, x, p[Gammarrindex]);
    double* dArr     = sder.FirstDerSpaceCenteredDiff4(ArrI    , sizex, x, p[Arrindex    ]);
    double* dTheta   = sder.FirstDerSpaceCenteredDiff4(ThetaI  , sizex, x, p[Thetaindex  ]);
    double* dtrK     = sder.FirstDerSpaceCenteredDiff4(trKI    , sizex, x, p[trKindex    ]);
    double* dA       = sder.FirstDerSpaceCenteredDiff4(AI      , sizex, x, p[Aindex      ]);

    for(int i = 0; i < sizex; i++){
        cPiI     [i] = SF*((-2.*Kcmc*AI[i]*AlphaI[i]*ChiI[i]*Omegaf[i]*Omegaf[i]*dA[i]- Kcmc*AI[i]*AI[i]*ChiI[i]*ChiI[i]*Omegaf[i]*Omegaf[i]*dAlpha[i] -Kcmc*AI[i]*AI[i]*AlphaI[i]*Omegaf[i]*dChi[i] +Kcmc*AI[i]*AI[i]*AlphaI[i]*ChiI[i]*Omegaf[i]*dChi[i] + 3.*Kcmc*AI[i]*AI[i]*AlphaI[i]*Omegaf[i]*Omegaf[i]*dChi[i] - 12.*M_PI*cphiI[i]*BetaI[i]*pow(ChiI[i],3)*Omegaf[i]*dcphi[i]*dOmega[i] - 12.*M_PI*cphiI[i]*cphiI[i]*BetaI[i]*pow(ChiI[i],3)*dOmega[i]*dOmega[i])/(12.*M_PI*AlphaI[i]*pow(ChiI[i],3)*Omegaf[i]*(Omegaf[i]*dcphi[i] + cphiI[i]*dOmega[i])));
        dPiI     [i] = SF*0.0;
        EI       [i] = EM*0.0;
    }

    std::cout << "Initial Conditions Set" << std::endl;
    
    ////////////////////////////////////
    ///////////TIME EVOLUTION///////////
    ////////////////////////////////////

    TimeEvolution te;

    if(Choice == "te"){

        System sys_ev(data, ddata, dt, x, variables, sizet, sizex, p, StepInstantWrite);

        std::cout << "Time Evolution System Set" << std::endl;

        te.RungeKutta4(sys_ev, RHS);

        std::cout << "Time Evolution Complete, use Evolution Vid or Evolution Vid Binary python file" << std::endl;

        sys_ev.WriteBin("output/output3");     //Writing all of data
        //sys_ev.Write("output/output.txt");     //Writing all of data
        sys_ev.WriteConstraint("output/Constraints"); //Writing Poisson constraint data
    }

    /////////////////////////////////////////////////////
    ///////////POINT-WISE CONVERGENCE EVOLUTION//////////
    /////////////////////////////////////////////////////

    if(Choice == "pc"){

        System sys_point_conv_ev(data, ddata, dt, x, variables, sizet, sizex, p, StepInstantWrite);

        std::cout << "Pointwise Convergence Evolution System Set" << std::endl;

        te.PointConvergenceEv(sys_point_conv_ev, 4., 3., true, RHS, "output/pointconvergence.txt");

        std::cout << "Pointwise Convergence Evolution Complete, use PointConv_Video.py to see the output video or PointConv_Partial.py to see pointwise convergence evolution for one of the variables (you have to choose in the python file). Use Constraint_Conv_Vid to see the convergence of the constraints used." << std::endl;
    }

    ///////////////////////////////////////////
    ////////////////L2 CONVERGENCE/////////////
    ///////////////////////////////////////////

    if(Choice == "l2total"){

        System sys_point_conv(data, ddata, dt, x, variables, sizet, sizex, p, StepInstantWrite);

        std::cout << "L2 Convergence System Set" << std::endl;

        te.L2Convergence(sys_point_conv, 4., 3., RHS, "output/convergence.txt");

        std::cout << "L2 Convergence Calculated, use L2Conv_Plot.py to see the output" << std::endl;
    }

    ///////////////////////////////////////////////////
    ///////////////PARTIAL L2 CONVERGENCE//////////////
    ///////////////////////////////////////////////////

    if(Choice == "l2"){

        System sys_partial_point(data, ddata, dt, x, variables, sizet, sizex, p, StepInstantWrite);

        std::cout << "Partial L2 Convergence System Set" << std::endl;

        //Don't forget to change to L2PartialConvergence only, take out the 2
        te.L2PartialConvergence(sys_partial_point, 4., 3., RHS, "output/partialconvergence.txt");

        std::cout << "Partial L2 Convergence Calculated, use L2_Individual_Plot.py to see the output" << std::endl;
    }

    for (int i = 0; i < variables; i++) delete[] data [i];
    delete[] data    ;
    delete[] p       ;
    delete[] x       ;
    delete[] Omegaf  ;
    delete[] dOmega  ;
    delete[] dcphi   ;
    delete[] ddphi   ;
    delete[] dA      ;
    delete[] dGammarr;
    delete[] dChi    ;
    delete[] dAlpha  ;
    delete[] dArr    ;
    delete[] dTheta  ;
    delete[] dtrK    ;
    delete[] RHSArI  ;

    return 0;
}
