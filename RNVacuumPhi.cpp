#include "TimeEvolution.h"
#include <matplotlibcpp.h>

//For Aconf

void Yaconf(quad* data, quad x, int sizex, quad& result, quad* aconf, BackgroundParameters<quad> BP){
    #ifdef DEBUG
        printf("[%s]\n", __PRETTY_FUNCTION__);
    #endif

    quad Y = data[0];

    quad Ccmc = BP.Ccmc;
    quad R0   = BP.R0;
    quad Q    = BP.QBH;
    quad mass = BP.M;
 
    result = sqrt(9*pow(Ccmc,2) + 6*Ccmc*Kcmc*pow(R0,3) + pow(R0,2)*(9*pow(Q,2) - 18*mass*R0 + 9*pow(R0,2) + pow(Kcmc,2)*pow(R0,4)) - 18*(3*pow(Ccmc,2) + Ccmc*Kcmc*pow(R0,3) + pow(R0,2)*(2*pow(Q,2) - 3*mass*R0 + pow(R0,2)))*Y + 9*(15*pow(Ccmc,2) + 2*Ccmc*Kcmc*pow(R0,3) + pow(R0,2)*(6*pow(Q,2) - 6*mass*R0 + pow(R0,2)))*pow(Y,2) - 6*(30*pow(Ccmc,2) + 6*pow(Q,2)*pow(R0,2) + Ccmc*Kcmc*pow(R0,3) - 3*mass*pow(R0,3))*pow(Y,3) + 9*(15*pow(Ccmc,2) + pow(Q,2)*pow(R0,2))*pow(Y,4) - 54*pow(Ccmc,2)*pow(Y,5) + 9*pow(Ccmc,2)*pow(Y,6))/(3.*x*pow(R0,2));
}

int main(int argc, char* argv[]){

    std::string Choice = "te"; //te, pc, l2. This chooses what type of run you want. Check the Makefile if you need

    // Check if there are command line arguments -> This is what we will use to run the different cases
    if(argc <= 1) std::cout << "No command line arguments provided, using default value for Time Evolution" << std::endl;
    else          Choice = std::string(argv[1]);

    int  variables(16);
    int  sizex    (200); //207 for single value
    int  sizet    (540*StepInstantWrite);
    quad dx       (rscri/sizex);
    quad dt       (0.001);

    /*For convergence tests
    int  sizex    (1800); //207 for single value
    int  sizet    (540*StepInstantWrite);
    quad dx       (rscri/sizex);
    quad dt       (0.00001);
    */

    std::cout << "CFL condition: " << dt/dx << std::endl;

    SpaceDerivator<quad> sder(Ghost_Left, Ghost_Right);
    
    quad M  (1.0);
    quad QBH(0.8);
    quad Ccmc =  sder.CcmcShootingMethod((quad) 0.01, (quad) 25.0, M, QBH);
    quad R0   =  sder.R0ShootingMethod  ((quad) 1.0  , M, QBH, Ccmc);

    sder.SetBackgroundParameters(M, QBH, Ccmc, R0);

    std::cout << "Background Parameters Set (Mass, Charge, Ccmc, R0)" << std::endl;

    using precision = long double;
    
    void (*RHS[])(precision**, precision**, precision*, int, precision, int*, precision*, BackgroundParameters<precision>) = {ElectricField, PsiElectric, Beta, Alpha, Phi, A, cphi, dphi, cPi, dPi, trK, Gammarr, Chi, Arr, Lambdar, Theta};

    //Parity conditions
    int* p = new int[variables]{pE, pPsi, pBeta, pAlpha, pPhi, pA, pcphi, pdphi, pcPi, pdPi, ptrK, pGammarr,  pChi,  pArr, pLambdar, pTheta};

    std::cout << "RHS Set" << std::endl;

    ////////////////////////////////////////
    ///////////INITIAL CONDITIONS///////////
    ////////////////////////////////////////

    quad* x = sder.Grid(sizex, true);

    quad valatscri    (1.0); //Value of the variable at scri (both aconf and chi should have this value)
    quad in_guess1  (200.0); //First guess for aconf at x=dx
    quad in_guess2 (-0.002); //Second guess for aconf at x=dx

    InitialData         <quad> id(sizex, variables, p, x);
    BackgroundParameters<quad> BP = createBackgroundParameters(M, QBH, Ccmc, R0);
    
    quad** data = id.RNVacuumWithPhi(valatscri, in_guess1, in_guess2, Yaconf, BP);

    std::cout << "Initial Conditions Set" << std::endl;

    ////////////////////////////////////
    ///////////TIME EVOLUTION///////////
    ////////////////////////////////////

    //We have to take care in changing from quad to precision or double
    //to make the simulation less expensive
    TimeEvolution<precision> te;

    if(Choice == "te"){
        SpaceDerivator<precision> sderld(Ghost_Left, Ghost_Right);

        //Converting to precision
        precision*  xld    = sderld.Grid(sizex, true);
        precision** datald = ConvertData<quad, precision>(data, variables, sizex);

        BackgroundParameters<precision> BPld = ConvertData<quad, precision>(BP);
        
        System<precision> sys_ev(datald, static_cast<precision>(dt), xld, variables, sizet, sizex, p, StepInstantWrite, BPld);

        std::cout << "Time Evolution System Set" << std::endl;

        te.RungeKutta4(sys_ev, RHS);

        std::cout << "Time Evolution Complete, use Evolution Vid or Evolution Vid Binary python file" << std::endl;

        sys_ev.WriteBin("output/output3");     //Writing all of data
        sys_ev.WriteConstraint("output/Constraints"); //Writing Poisson constraint data

        for(int i = 0; i < variables; i++) delete[] datald[i];
        delete[] datald;
        delete[] xld;
        CleanBP(BPld);
        CleanBP(BP);
    }
    
    /////////////////////////////////////////////////////
    ///////////POINT-WISE CONVERGENCE EVOLUTION//////////
    /////////////////////////////////////////////////////

    if(Choice == "pc"){
        //Resolution decrease factor
        int f = 3;

        quad* xMR = sder.Grid(sizex/f    , true);
        quad* xLR = sder.Grid(sizex/(f*f), true);

        InitialData<quad> idMR(sizex/f    , variables, p, xMR);
        InitialData<quad> idLR(sizex/(f*f), variables, p, xLR);

        quad** dataMR = idMR.RNVacuumWithPhi(valatscri, in_guess1, in_guess2, Yaconf, BP);
        quad** dataLR = idLR.RNVacuumWithPhi(valatscri, in_guess1, in_guess2, Yaconf, BP);

        //Converting to precision
        SpaceDerivator<precision> sderld(Ghost_Left, Ghost_Right);

        precision*  xHRld = sderld.Grid(sizex      , true);
        precision*  xMRld = sderld.Grid(sizex/f    , true);
        precision*  xLRld = sderld.Grid(sizex/(f*f), true);

        precision** dataHRld = ConvertData<quad, precision>(data  , variables,       sizex);
        precision** dataMRld = ConvertData<quad, precision>(dataMR, variables,     sizex/f);
        precision** dataLRld = ConvertData<quad, precision>(dataLR, variables, sizex/(f*f));

        BackgroundParameters<precision> BPld = ConvertData<quad, precision>(BP);

        System<precision> sysHR(dataHRld, static_cast<precision>(dt)    , xHRld, variables, sizet      , sizex      , p, StepInstantWrite      , BPld);
        System<precision> sysMR(dataMRld, static_cast<precision>(dt*f)  , xMRld, variables, sizet/f    , sizex/f    , p, StepInstantWrite/f    , BPld);
        System<precision> sysLR(dataLRld, static_cast<precision>(dt*f*f), xLRld, variables, sizet/(f*f), sizex/(f*f), p, StepInstantWrite/(f*f), BPld);

        std::cout << "Pointwise Convergence Evolution System Set" << std::endl;

        te.PointConvergenceEv(sysHR, sysMR, sysLR, 4., 3., true, RHS);

        std::cout << "Pointwise Convergence Evolution Complete, use PointConv_Video.py to see the output video or PointConv_Partial.py to see pointwise convergence evolution for one of the variables (you have to choose in the python file). Use Constraint_Conv_Vid to see the convergence of the constraints used." << std::endl;

        for(int i = 0; i < variables; i++) { delete[] dataMR[i]; delete[] dataLR[i]; delete[] dataHRld[i]; delete[] dataMRld[i]; delete[] dataLRld[i]; }
        delete[] dataMR;
        delete[] dataLR;
        delete[] dataHRld;
        delete[] dataMRld;
        delete[] dataLRld;
        delete[] xMR;
        delete[] xLR;
        delete[] xHRld;
        delete[] xMRld;
        delete[] xLRld;
        CleanBP(BPld);
        CleanBP(BP);
    }

    for(int i = 0; i < variables; i++) { delete[] data[i]; }    
    delete[] data;
    delete[] p   ;
    delete[] x   ;
    
    return 0;
}