#include "TimeEvolution.h"
#include <matplotlibcpp.h>

/*
TimeEvolutionF:            runs usual RK4 time evolution
PointConvergenceEvolution: runs pointwise convergence evolution. Poisson convergence is calculated here as well
PointConvergence:          runs pointwise convergence at the end of the simulation
PartialPointConvergence:   runs pointwise convergence for each variable, at the end of the simulation
*/

//For multithreading when writing
//Define a mutex to synchronize file access
std::mutex fileMutex;

quad cphigaussian(quad x) {
    quad sigma = 0.2;
    quad mu2   = 0.5;
    quad Acphi = 1.0;
    
    return Acphi*x*x*exp(-(x*x - mu2*mu2) * (x*x - mu2*mu2) / (4. * sigma * sigma * sigma * sigma));
}

quad cphigaussian_derivative(quad x) {
    quad sigma = 0.2;
    quad mu2   = 0.5;
    quad Acphi = 1.0;
    
    quad exponent = -(x*x - mu2*mu2) * (x*x - mu2*mu2) / (4. * sigma * sigma * sigma * sigma);
    quad exp_term = exp(exponent);

    quad term1 = 2 * Acphi * x * exp_term; // g'(x) * h(x)
    quad term2 = -Acphi * x * x * (x*x - mu2*mu2) * x / (sigma * sigma * sigma * sigma) * exp_term; // g(x) * h'(x)

    return term1 + term2;
}

//These functions come from initial conditions for chi with A3r and cphi working
void Y(quad* data, quad x, int sizex, quad& result, quad* aconf, BackgroundParameters<quad> BP){
    #ifdef DEBUG
        printf("[%s]\n", __PRETTY_FUNCTION__);
    #endif

    result = data[1];
}

void U(quad* data, quad x, int sizex, quad& result, quad* aconf, BackgroundParameters<quad> BP){
    #ifdef DEBUG
        printf("[%s]\n", __PRETTY_FUNCTION__);
    #endif

    //Giving names
    quad Y  = data[0];
    quad U  = data[1];

    quad dcphi = cphigaussian_derivative(x);
    quad cphi  = cphigaussian(x);
    
    result = (-3*Y)/pow(-1 + pow(x,2),2) - (M_PI*pow(x,2)*pow(cphi,2)*Y)/(9.*pow(-1 + pow(x,2),2)) + (2*M_PI*pow(x,4)*pow(cphi,2)*Y)/(9.*pow(-1 + pow(x,2),2)) - (M_PI*pow(x,6)*pow(cphi,2)*Y)/(9.*pow(-1 + pow(x,2),2)) + (3*pow(Y,5))/pow(-1 + pow(x,2),2) + (M_PI*x*cphi*Y*dcphi)/(9.*pow(-1 + pow(x,2),2)) - (M_PI*pow(x,3)*cphi*Y*dcphi)/(3.*pow(-1 + pow(x,2),2)) + (M_PI*pow(x,5)*cphi*Y*dcphi)/(3.*pow(-1 + pow(x,2),2)) - (M_PI*pow(x,7)*cphi*Y*dcphi)/(9.*pow(-1 + pow(x,2),2)) - (M_PI*Y*pow(dcphi,2))/(36.*pow(-1 + pow(x,2),2)) + (M_PI*pow(x,2)*Y*pow(dcphi,2))/(9.*pow(-1 + pow(x,2),2)) - (M_PI*pow(x,4)*Y*pow(dcphi,2))/(6.*pow(-1 + pow(x,2),2)) + (M_PI*pow(x,6)*Y*pow(dcphi,2))/(9.*pow(-1 + pow(x,2),2)) - (M_PI*pow(x,8)*Y*pow(dcphi,2))/(36.*pow(-1 + pow(x,2),2)) - (2*U)/(x*pow(-1 + pow(x,2),2)) + (2*x*U)/pow(-1 + pow(x,2),2);
}

void U0(quad* data, quad x, int sizex, quad& result, BackgroundParameters<quad> BP){
    #ifdef DEBUG
        printf("[%s]\n", __PRETTY_FUNCTION__);
    #endif

    //Giving names
    quad Y  = data[0];
    quad U  = data[1];

    quad dcphi = cphigaussian_derivative(x);
    quad cphi  = cphigaussian(x);
    
    result = ((-3*Y)/pow(-1 + pow(x,2),2) - (M_PI*pow(x,2)*pow(cphi,2)*Y)/(9.*pow(-1 + pow(x,2),2)) + (2*M_PI*pow(x,4)*pow(cphi,2)*Y)/(9.*pow(-1 + pow(x,2),2)) - (M_PI*pow(x,6)*pow(cphi,2)*Y)/(9.*pow(-1 + pow(x,2),2)) + (3*pow(Y,5))/pow(-1 + pow(x,2),2) + (M_PI*x*cphi*Y*dcphi)/(9.*pow(-1 + pow(x,2),2)) - (M_PI*pow(x,3)*cphi*Y*dcphi)/(3.*pow(-1 + pow(x,2),2)) + (M_PI*pow(x,5)*cphi*Y*dcphi)/(3.*pow(-1 + pow(x,2),2)) - (M_PI*pow(x,7)*cphi*Y*dcphi)/(9.*pow(-1 + pow(x,2),2)) - (M_PI*Y*pow(dcphi,2))/(36.*pow(-1 + pow(x,2),2)) + (M_PI*pow(x,2)*Y*pow(dcphi,2))/(9.*pow(-1 + pow(x,2),2)) - (M_PI*pow(x,4)*Y*pow(dcphi,2))/(6.*pow(-1 + pow(x,2),2)) + (M_PI*pow(x,6)*Y*pow(dcphi,2))/(9.*pow(-1 + pow(x,2),2)) - (M_PI*pow(x,8)*Y*pow(dcphi,2))/(36.*pow(-1 + pow(x,2),2)) + (2*x*U)/pow(-1 + pow(x,2),2))/3.0;
}

int main(int argc, char* argv[]){

    std::string Choice = "te"; //te, pc, l2. This chooses what type of run you want. Check the Makefile if you need

    // Check if there are command line arguments -> This is what we will use to run the different cases
    if(argc <= 1) std::cout << "No command line arguments provided, using default value for Time Evolution" << std::endl;
    else          Choice = std::string(argv[1]);

    int  variables(16) ; //Number of variables; Take care with this number
    int  sizex    (1800); //9000, 3000, 1000
    int  sizet    (18*StepInstantWrite); //Ensure that this is always an integer multiple of StepInstantWrite
    quad dx       (rscri/sizex);
    quad dt       (0.00001)    ; //0.000018, 0.00006, 0.00002

    std::cout << "CFL condition: " << dt/dx << std::endl;

    using precision = long double;

    //Right Hand Sides
    void (*RHS[])(precision**, precision**, precision*, int, precision, int*, precision*, BackgroundParameters<precision>) = {ElectricField, PsiElectric, Beta, Alpha, Phi, A, cphi, dphi, cPi, dPi, trK, Gammarr, Chi, Arr, Lambdar, Theta};

    //Parity conditions: vectors are odd, scalars are even
                               //E, Psi, Beta, Alpha, Phi, A, cphi, dphi, cPi, dPi, trK, Gammarr, Chi, Arr, Lambdar, Theta
    int* p = new int[variables]{pE, pPsi, pBeta, pAlpha, pPhi, pA, pcphi, pdphi, pcPi, pdPi, ptrK, pGammarr,  pChi,  pArr, pLambdar, pTheta};

    std::cout << "RHS Set" << std::endl;

    ////////////////////////////////////////
    ///////////INITIAL CONDITIONS///////////
    ////////////////////////////////////////
    SpaceDerivator<quad> sder(Ghost_Left, Ghost_Right);

    quad* x = sder.Grid(sizex, true);

    quad valatscri = 1.0  ;  //Value that chi should have at scri
    quad in_guess1 = 0.999;  //First  guess for the value of chi at x=dx/2
    quad in_guess2 = 1.015;  //Second guess for the value of chi at x=dx/2

    //Initial Conditions

    InitialData<quad> id(sizex, variables, p, x);

    quad** data = id.RealScalar(valatscri, in_guess1, in_guess2, Y, U, U0, cphigaussian);

    std::cout << "Initial Conditions Set" << std::endl;
    
    ////////////////////////////////////
    ///////////TIME EVOLUTION///////////
    ////////////////////////////////////

    TimeEvolution<precision> te;

    if(Choice == "te"){

        SpaceDerivator<precision> sderld(Ghost_Left, Ghost_Right);

        //Converting to precision
        precision*  xld    = sderld.Grid(static_cast<precision>(sizex), true);
        precision** datald = ConvertData<quad, precision>(data, variables, sizex);

        System<precision> sys_ev(datald, static_cast<precision>(dt), xld, variables, static_cast<precision>(sizet), static_cast<precision>(sizex), p, StepInstantWrite);

        std::cout << "Time Evolution System Set" << std::endl;

        te.RungeKutta4(sys_ev, RHS);

        std::cout << "Time Evolution Complete, use Evolution Vid or Evolution Vid Binary python file" << std::endl;

        sys_ev.WriteBin("output/output3");     //Writing all of data
        sys_ev.WriteConstraint("output/Constraints"); //Writing Poisson constraint data
    }

    /////////////////////////////////////////////////////
    ///////////POINT-WISE CONVERGENCE EVOLUTION//////////
    /////////////////////////////////////////////////////

    if(Choice == "pc"){
        //Resolution decrease factor
        int f = 3;

        SpaceDerivator<precision> sderld(Ghost_Left, Ghost_Right);

        quad* xMR = sder.Grid(sizex/f    , true);
        quad* xLR = sder.Grid(sizex/(f*f), true);

        InitialData<quad> idMR(sizex/f    , variables, p, xMR);
        InitialData<quad> idLR(sizex/(f*f), variables, p, xLR);

        quad** dataMR = idMR.RealScalar(valatscri, in_guess1, in_guess2, Y, U, U0, cphigaussian);
        quad** dataLR = idLR.RealScalar(valatscri, in_guess1, in_guess2, Y, U, U0, cphigaussian);

        precision*  xHRld = sderld.Grid(static_cast<precision>(sizex)      , true);
        precision*  xMRld = sderld.Grid(static_cast<precision>(sizex/f)    , true);
        precision*  xLRld = sderld.Grid(static_cast<precision>(sizex/(f*f)), true);

        precision** dataHRld = ConvertData<quad, precision>(data  , variables,       sizex);
        precision** dataMRld = ConvertData<quad, precision>(dataMR, variables,     sizex/f);
        precision** dataLRld = ConvertData<quad, precision>(dataLR, variables, sizex/(f*f));

        System<precision> sysHR(dataHRld, static_cast<precision>(dt)    , xHRld, variables, static_cast<precision>(sizet)      , static_cast<precision>(sizex)      , p, StepInstantWrite);
        System<precision> sysMR(dataMRld, static_cast<precision>(dt*f)  , xMRld, variables, static_cast<precision>(sizet/f)    , static_cast<precision>(sizex/f)    , p, StepInstantWrite/f);
        System<precision> sysLR(dataLRld, static_cast<precision>(dt*f*f), xLRld, variables, static_cast<precision>(sizet/(f*f)), static_cast<precision>(sizex/(f*f)), p, StepInstantWrite/(f*f));

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
    }

    ///////////////////////////////////////////
    ////////////////L2 CONVERGENCE/////////////
    ///////////////////////////////////////////

    if(Choice == "l2total"){
        //Resolution decrease factor
        int f = 3;

        SpaceDerivator<precision> sderld(Ghost_Left, Ghost_Right);

        quad* xMR = sder.Grid(sizex/f    , true);
        quad* xLR = sder.Grid(sizex/(f*f), true);

        InitialData<quad> idMR(sizex/f    , variables, p, xMR);
        InitialData<quad> idLR(sizex/(f*f), variables, p, xLR);

        quad** dataMR = idMR.RealScalar(valatscri, in_guess1, in_guess2, Y, U, U0, cphigaussian);
        quad** dataLR = idLR.RealScalar(valatscri, in_guess1, in_guess2, Y, U, U0, cphigaussian);

        precision*  xHRld = sderld.Grid(static_cast<precision>(sizex)      , true);
        precision*  xMRld = sderld.Grid(static_cast<precision>(sizex/f)    , true);
        precision*  xLRld = sderld.Grid(static_cast<precision>(sizex/(f*f)), true);

        precision** dataHRld = ConvertData<quad, precision>(data  , variables,       sizex);
        precision** dataMRld = ConvertData<quad, precision>(dataMR, variables,     sizex/f);
        precision** dataLRld = ConvertData<quad, precision>(dataLR, variables, sizex/(f*f));

        System<precision> sysHR(dataHRld, static_cast<precision>(dt)    , xHRld, variables, static_cast<precision>(sizet)      , static_cast<precision>(sizex)      , p, StepInstantWrite);
        System<precision> sysMR(dataMRld, static_cast<precision>(dt*f)  , xMRld, variables, static_cast<precision>(sizet/f)    , static_cast<precision>(sizex/f)    , p, StepInstantWrite/f);
        System<precision> sysLR(dataLRld, static_cast<precision>(dt*f*f), xLRld, variables, static_cast<precision>(sizet/(f*f)), static_cast<precision>(sizex/(f*f)), p, StepInstantWrite/(f*f));

        std::cout << "L2 Convergence System Set" << std::endl;

        te.L2Convergence(sysHR, sysMR, sysLR, 4., 3., true, RHS);

        std::cout << "L2 Convergence Calculated, use L2Conv_Plot.py to see the output" << std::endl;

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
    }

    ///////////////////////////////////////////////////
    ///////////////PARTIAL L2 CONVERGENCE//////////////
    ///////////////////////////////////////////////////

    if(Choice == "l2"){
        //Resolution decrease factor
        int f = 3;

        SpaceDerivator<precision> sderld(Ghost_Left, Ghost_Right);

        quad* xMR = sder.Grid(sizex/f    , true);
        quad* xLR = sder.Grid(sizex/(f*f), true);

        InitialData<quad> idMR(sizex/f    , variables, p, xMR);
        InitialData<quad> idLR(sizex/(f*f), variables, p, xLR);

        quad** dataMR = idMR.RealScalar(valatscri, in_guess1, in_guess2, Y, U, U0, cphigaussian);
        quad** dataLR = idLR.RealScalar(valatscri, in_guess1, in_guess2, Y, U, U0, cphigaussian);

        precision*  xHRld = sderld.Grid(static_cast<precision>(sizex)      , true);
        precision*  xMRld = sderld.Grid(static_cast<precision>(sizex/f)    , true);
        precision*  xLRld = sderld.Grid(static_cast<precision>(sizex/(f*f)), true);

        precision** dataHRld = ConvertData<quad, precision>(data  , variables,       sizex);
        precision** dataMRld = ConvertData<quad, precision>(dataMR, variables,     sizex/f);
        precision** dataLRld = ConvertData<quad, precision>(dataLR, variables, sizex/(f*f));

        System<precision> sysHR(dataHRld, static_cast<precision>(dt)    , xHRld, variables, static_cast<precision>(sizet)      , static_cast<precision>(sizex)      , p, StepInstantWrite);
        System<precision> sysMR(dataMRld, static_cast<precision>(dt*f)  , xMRld, variables, static_cast<precision>(sizet/f)    , static_cast<precision>(sizex/f)    , p, StepInstantWrite/f);
        System<precision> sysLR(dataLRld, static_cast<precision>(dt*f*f), xLRld, variables, static_cast<precision>(sizet/(f*f)), static_cast<precision>(sizex/(f*f)), p, StepInstantWrite/(f*f));

        std::cout << "Partial L2 Convergence System Set" << std::endl;

        //Don't forget to change to L2PartialConvergence only, take out the 2
        te.L2PartialConvergence(sysHR, sysMR, sysLR, 4., 3., true, RHS);

        std::cout << "Partial L2 Convergence Calculated, use L2_Individual_Plot.py to see the output" << std::endl;
        
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
    }

    for (int i = 0; i < variables; i++){ delete[] data [i];}
    delete[] data;
    delete[] p   ;
    delete[] x   ;

    return 0;
}
