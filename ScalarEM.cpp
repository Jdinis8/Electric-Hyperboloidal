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

double gaussianA3r(double x) {
    double sigma = 0.2;
    double mu2 = 0.25;
    double A = 0.1;
    return ASw*A*x*x*exp(-(x*x - mu2*mu2) * (x*x - mu2*mu2) / (4. * sigma * sigma * sigma * sigma));
}

double gaussianCphi(double x) {
    double sigma = 0.2;
    double mu2 = 0.4;
    double A = 1.;
    return cphiSw*A*x*x*exp(-(x*x - mu2*mu2) * (x*x - mu2*mu2) / (4. * sigma * sigma * sigma * sigma));
}

double dgaussianA3r(double x) {
    double sigma = 0.2;
    double mu2 = 0.25;
    double A = 0.1;

    double exp_term = exp(-(x*x - mu2*mu2) * (x*x - mu2*mu2) / (4. * sigma * sigma * sigma * sigma));
    double common_term = ASw * A * x * exp_term * (x*x - mu2*mu2) / (sigma * sigma * sigma * sigma);

    return 2 * ASw * A * x * exp_term - common_term * x;
}

double dgaussianCphi(double x) {
    double sigma = 0.2;
    double mu2 = 0.4;
    double A = 1.0;

    double exp_term = exp(-(x*x - mu2*mu2) * (x*x - mu2*mu2) / (4. * sigma * sigma * sigma * sigma));
    double common_term = cphiSw * A * x * exp_term * (x*x - mu2*mu2) / (sigma * sigma * sigma * sigma);

    return 2 * cphiSw * A * x * exp_term - common_term * x;
}

//These functions come from initial conditions for chi with A3r and cphi working
void Y(double* data, double x, int sizex, double& result, double* aconf, BackgroundParameters<double> BP){
    #ifdef DEBUG
        printf("[%s]\n", __PRETTY_FUNCTION__);
    #endif

    result = data[1];
}

void U(double* data, double x, int sizex, double& result, double* aconf, BackgroundParameters<double> BP){
    #ifdef DEBUG
        printf("[%s]\n", __PRETTY_FUNCTION__);
    #endif

    //Giving names
    double Y  = data[0];
    double U  = data[1];

    result = (3.*x*Y)/(pow(-1. + pow(x,2),2)*(-1.*x + 1.*pow(x,3))) - (3.*pow(x,3)*Y)/(pow(-1. + pow(x,2),2)*(-1.*x + 1.*pow(x,3))) + (0.0001170985474858286*pow(exp(1),99.99999999999999*pow(x,2) - 312.4999999999999*pow(x,4))*pow(x,3)*Y)/(pow(-1. + pow(x,2),2)*(-1.*x + 1.*pow(x,3))) + (0.01089016491618206*pow(exp(1),99.99999999999999*pow(x,2) - 312.4999999999999*pow(x,4))*pow(x,5)*Y)/(pow(-1. + pow(x,2),2)*(-1.*x + 1.*pow(x,3))) + (0.1515255204466622*pow(exp(1),99.99999999999999*pow(x,2) - 312.4999999999999*pow(x,4))*pow(x,7)*Y)/(pow(-1. + pow(x,2),2)*(-1.*x + 1.*pow(x,3))) - (4.522931396640128*pow(exp(1),99.99999999999999*pow(x,2) - 312.4999999999999*pow(x,4))*pow(x,9)*Y)/(pow(-1. + pow(x,2),2)*(-1.*x + 1.*pow(x,3))) + (31.449420370000812*pow(exp(1),99.99999999999999*pow(x,2) - 312.4999999999999*pow(x,4))*pow(x,11)*Y)/(pow(-1. + pow(x,2),2)*(-1.*x + 1.*pow(x,3))) - (95.42187914262627*pow(exp(1),99.99999999999999*pow(x,2) - 312.4999999999999*pow(x,4))*pow(x,13)*Y)/(pow(-1. + pow(x,2),2)*(-1.*x + 1.*pow(x,3))) + (151.72897917291877*pow(exp(1),99.99999999999999*pow(x,2) - 312.4999999999999*pow(x,4))*pow(x,15)*Y)/(pow(-1. + pow(x,2),2)*(-1.*x + 1.*pow(x,3))) - (132.79707150814744*pow(exp(1),99.99999999999999*pow(x,2) - 312.4999999999999*pow(x,4))*pow(x,17)*Y)/(pow(-1. + pow(x,2),2)*(-1.*x + 1.*pow(x,3))) + (60.83635474849686*pow(exp(1),99.99999999999999*pow(x,2) - 312.4999999999999*pow(x,4))*pow(x,19)*Y)/(pow(-1. + pow(x,2),2)*(-1.*x + 1.*pow(x,3))) - (11.435405027912942*pow(exp(1),99.99999999999999*pow(x,2) - 312.4999999999999*pow(x,4))*pow(x,21)*Y)/(pow(-1. + pow(x,2),2)*(-1.*x + 1.*pow(x,3))) - (3.*x*pow(Y,5))/(pow(-1. + pow(x,2),2)*(-1.*x + 1.*pow(x,3))) + (3.*pow(x,3)*pow(Y,5))/(pow(-1. + pow(x,2),2)*(-1.*x + 1.*pow(x,3))) + (5.52747592215012e-8*pow(exp(1),139.06249999999997*pow(x,2) - 624.9999999999998*pow(x,4))*pow(x,9)*pow(Y,9))/(pow(-1. + pow(x,2),2)*(-1.*x + 1.*pow(x,3))) - (2.76373796107506e-7*pow(exp(1),139.06249999999997*pow(x,2) - 624.9999999999998*pow(x,4))*pow(x,11)*pow(Y,9))/(pow(-1. + pow(x,2),2)*(-1.*x + 1.*pow(x,3))) + (5.52747592215012e-7*pow(exp(1),139.06249999999997*pow(x,2) - 624.9999999999998*pow(x,4))*pow(x,13)*pow(Y,9))/(pow(-1. + pow(x,2),2)*(-1.*x + 1.*pow(x,3))) - (5.52747592215012e-7*pow(exp(1),139.06249999999997*pow(x,2) - 624.9999999999998*pow(x,4))*pow(x,15)*pow(Y,9))/(pow(-1. + pow(x,2),2)*(-1.*x + 1.*pow(x,3))) + (2.76373796107506e-7*pow(exp(1),139.06249999999997*pow(x,2) - 624.9999999999998*pow(x,4))*pow(x,17)*pow(Y,9))/(pow(-1. + pow(x,2),2)*(-1.*x + 1.*pow(x,3))) - (5.52747592215012e-8*pow(exp(1),139.06249999999997*pow(x,2) - 624.9999999999998*pow(x,4))*pow(x,19)*pow(Y,9))/(pow(-1. + pow(x,2),2)*(-1.*x + 1.*pow(x,3))) + (2.0000000000000004*U)/(pow(-1. + pow(x,2),2)*(-1.*x + 1.*pow(x,3))) - (4.000000000000001*pow(x,2)*U)/(pow(-1. + pow(x,2),2)*(-1.*x + 1.*pow(x,3))) + (2.0000000000000004*pow(x,4)*U)/(pow(-1. + pow(x,2),2)*(-1.*x + 1.*pow(x,3)));
}

void U0(double* data, double x, int sizex, double& result, BackgroundParameters<double> BP){
    #ifdef DEBUG
        printf("[%s]\n", __PRETTY_FUNCTION__);
    #endif

    //Giving names
    double Y  = data[0];
    double U  = data[1];

    result = ((3.*x*Y)/(pow(-1. + pow(x,2),2)*(-1.*x + 1.*pow(x,3))) - (3.*pow(x,3)*Y)/(pow(-1. + pow(x,2),2)*(-1.*x + 1.*pow(x,3))) + (0.0001170985474858286*pow(exp(1),99.99999999999999*pow(x,2) - 312.4999999999999*pow(x,4))*pow(x,3)*Y)/(pow(-1. + pow(x,2),2)*(-1.*x + 1.*pow(x,3))) + (0.01089016491618206*pow(exp(1),99.99999999999999*pow(x,2) - 312.4999999999999*pow(x,4))*pow(x,5)*Y)/(pow(-1. + pow(x,2),2)*(-1.*x + 1.*pow(x,3))) + (0.1515255204466622*pow(exp(1),99.99999999999999*pow(x,2) - 312.4999999999999*pow(x,4))*pow(x,7)*Y)/(pow(-1. + pow(x,2),2)*(-1.*x + 1.*pow(x,3))) - (4.522931396640128*pow(exp(1),99.99999999999999*pow(x,2) - 312.4999999999999*pow(x,4))*pow(x,9)*Y)/(pow(-1. + pow(x,2),2)*(-1.*x + 1.*pow(x,3))) + (31.449420370000812*pow(exp(1),99.99999999999999*pow(x,2) - 312.4999999999999*pow(x,4))*pow(x,11)*Y)/(pow(-1. + pow(x,2),2)*(-1.*x + 1.*pow(x,3))) - (95.42187914262627*pow(exp(1),99.99999999999999*pow(x,2) - 312.4999999999999*pow(x,4))*pow(x,13)*Y)/(pow(-1. + pow(x,2),2)*(-1.*x + 1.*pow(x,3))) + (151.72897917291877*pow(exp(1),99.99999999999999*pow(x,2) - 312.4999999999999*pow(x,4))*pow(x,15)*Y)/(pow(-1. + pow(x,2),2)*(-1.*x + 1.*pow(x,3))) - (132.79707150814744*pow(exp(1),99.99999999999999*pow(x,2) - 312.4999999999999*pow(x,4))*pow(x,17)*Y)/(pow(-1. + pow(x,2),2)*(-1.*x + 1.*pow(x,3))) + (60.83635474849686*pow(exp(1),99.99999999999999*pow(x,2) - 312.4999999999999*pow(x,4))*pow(x,19)*Y)/(pow(-1. + pow(x,2),2)*(-1.*x + 1.*pow(x,3))) - (11.435405027912942*pow(exp(1),99.99999999999999*pow(x,2) - 312.4999999999999*pow(x,4))*pow(x,21)*Y)/(pow(-1. + pow(x,2),2)*(-1.*x + 1.*pow(x,3))) - (3.*x*pow(Y,5))/(pow(-1. + pow(x,2),2)*(-1.*x + 1.*pow(x,3))) + (3.*pow(x,3)*pow(Y,5))/(pow(-1. + pow(x,2),2)*(-1.*x + 1.*pow(x,3))) + (5.52747592215012e-8*pow(exp(1),139.06249999999997*pow(x,2) - 624.9999999999998*pow(x,4))*pow(x,9)*pow(Y,9))/(pow(-1. + pow(x,2),2)*(-1.*x + 1.*pow(x,3))) - (2.76373796107506e-7*pow(exp(1),139.06249999999997*pow(x,2) - 624.9999999999998*pow(x,4))*pow(x,11)*pow(Y,9))/(pow(-1. + pow(x,2),2)*(-1.*x + 1.*pow(x,3))) + (5.52747592215012e-7*pow(exp(1),139.06249999999997*pow(x,2) - 624.9999999999998*pow(x,4))*pow(x,13)*pow(Y,9))/(pow(-1. + pow(x,2),2)*(-1.*x + 1.*pow(x,3))) - (5.52747592215012e-7*pow(exp(1),139.06249999999997*pow(x,2) - 624.9999999999998*pow(x,4))*pow(x,15)*pow(Y,9))/(pow(-1. + pow(x,2),2)*(-1.*x + 1.*pow(x,3))) + (2.76373796107506e-7*pow(exp(1),139.06249999999997*pow(x,2) - 624.9999999999998*pow(x,4))*pow(x,17)*pow(Y,9))/(pow(-1. + pow(x,2),2)*(-1.*x + 1.*pow(x,3))) - (5.52747592215012e-8*pow(exp(1),139.06249999999997*pow(x,2) - 624.9999999999998*pow(x,4))*pow(x,19)*pow(Y,9))/(pow(-1. + pow(x,2),2)*(-1.*x + 1.*pow(x,3))))/3.0;
}

int main(int argc, char* argv[]){

    std::string Choice = "te"; //te, pc, l2. This chooses what type of run you want. Check the Makefile if you need

    // Check if there are command line arguments -> This is what we will use to run the different cases
    if(argc <= 1) std::cout << "No command line arguments provided, using default value for Time Evolution" << std::endl;
    else          Choice = std::string(argv[1]);

    int    variables(16) ; //Number of variables; Take care with this number
    int    sizex    (900); //9000, 3000, 1000
    int    sizet    (45*StepInstantWrite); //Ensure that this is always an integer multiple of StepInstantWrite
    double dx       (rscri/sizex);
    double dt       (0.00001)    ; //0.000018, 0.00006, 0.00002

    std::cout << "CFL condition: " << dt/dx << std::endl;

    //Right Hand Sides
    void (*RHS[])(double**, double**, double*, int, double, int*, double*, BackgroundParameters<double>) = {ElectricField, PsiElectric, Beta, Alpha, Phi, A, cphi, dphi, cPi, dPi, trK, Gammarr, Chi, Arr, Lambdar, Theta};

    //Parity conditions: vectors are odd, scalars are even
                               //E, Psi, Beta, Alpha, Phi, A, cphi, dphi, cPi, dPi, trK, Gammarr, Chi, Arr, Lambdar, Theta
    int* p = new int[variables]{pE, pPsi, pBeta, pAlpha, pPhi, pA, pcphi, pdphi, pcPi, pdPi, ptrK, pGammarr,  pChi,  pArr, pLambdar, pTheta};

    std::cout << "RHS Set" << std::endl;

    ////////////////////////////////////////
    ///////////INITIAL CONDITIONS///////////
    ////////////////////////////////////////
    SpaceDerivator<double> sder(Ghost_Left, Ghost_Right);

    double* x = sder.Grid(sizex, true);

    double valatscri  = 1.0  ;  //Value that chi should have at scri
    double in_guess1  = 0.995;  //First  guess for the value of chi at x=dx/2
    double in_guess2  = 1.015;  //Second guess for the value of chi at x=dx/2

    //Initial Conditions

    InitialData<double> id(sizex, variables, p, x);

    double** data = id.RealScalar_A3r(valatscri, in_guess1, in_guess2, Y, U, U0, gaussianCphi, gaussianA3r);

    std::cout << "Initial Conditions Set" << std::endl;
    
    ////////////////////////////////////
    ///////////TIME EVOLUTION///////////
    ////////////////////////////////////

    TimeEvolution<double> te;

    if(Choice == "te"){

        System<double> sys_ev(data, dt, x, variables, sizet, sizex, p, StepInstantWrite);

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

        double* xMR = sder.Grid(sizex/f    , true);
        double* xLR = sder.Grid(sizex/(f*f), true);

        InitialData<double> idMR(sizex/f    , variables, p, xMR);
        InitialData<double> idLR(sizex/(f*f), variables, p, xLR);

        double** dataMR = idMR.RealScalar_A3r(valatscri, in_guess1, in_guess2, Y, U, U0, gaussianCphi, gaussianA3r);
        double** dataLR = idLR.RealScalar_A3r(valatscri, in_guess1, in_guess2, Y, U, U0, gaussianCphi, gaussianA3r);

        System<double> sysHR(data  , dt    , x  , variables, sizet      , sizex      , p, StepInstantWrite);
        System<double> sysMR(dataMR, dt*f  , xMR, variables, sizet/f    , sizex/f    , p, StepInstantWrite/f);
        System<double> sysLR(dataLR, dt*f*f, xLR, variables, sizet/(f*f), sizex/(f*f), p, StepInstantWrite/(f*f));

        std::cout << "Pointwise Convergence Evolution System Set" << std::endl;

        te.PointConvergenceEv(sysHR, sysMR, sysLR, 4., 3., true, RHS);

        std::cout << "Pointwise Convergence Evolution Complete, use PointConv_Video.py to see the output video or PointConv_Partial.py to see pointwise convergence evolution for one of the variables (you have to choose in the python file). Use Constraint_Conv_Vid to see the convergence of the constraints used." << std::endl;
    }

    ///////////////////////////////////////////
    ////////////////L2 CONVERGENCE/////////////
    ///////////////////////////////////////////

    if(Choice == "l2total"){
        //Resolution decrease factor
        int f = 3;

        double* xMR = sder.Grid(sizex/f    , true);
        double* xLR = sder.Grid(sizex/(f*f), true);

        InitialData<double> idMR(sizex/f    , variables, p, xMR);
        InitialData<double> idLR(sizex/(f*f), variables, p, xLR);

        double** dataMR = idMR.RealScalar_A3r(valatscri, in_guess1, in_guess2, Y, U, U0, gaussianCphi, gaussianA3r);
        double** dataLR = idLR.RealScalar_A3r(valatscri, in_guess1, in_guess2, Y, U, U0, gaussianCphi, gaussianA3r);

        System<double> sysHR(data  , dt    , x  , variables, sizet      , sizex      , p, StepInstantWrite);
        System<double> sysMR(dataMR, dt*f  , xMR, variables, sizet/f    , sizex/f    , p, StepInstantWrite/f);
        System<double> sysLR(dataLR, dt*f*f, xLR, variables, sizet/(f*f), sizex/(f*f), p, StepInstantWrite/(f*f));

        std::cout << "L2 Convergence System Set" << std::endl;

        te.L2Convergence(sysHR, sysMR, sysLR, 4., 3., true, RHS);

        std::cout << "L2 Convergence Calculated, use L2Conv_Plot.py to see the output" << std::endl;
    }

    ///////////////////////////////////////////////////
    ///////////////PARTIAL L2 CONVERGENCE//////////////
    ///////////////////////////////////////////////////

    if(Choice == "l2"){
        //Resolution decrease factor
        int f = 3;

        double* xMR = sder.Grid(sizex/f    , true);
        double* xLR = sder.Grid(sizex/(f*f), true);

        InitialData<double> idMR(sizex/f    , variables, p, xMR);
        InitialData<double> idLR(sizex/(f*f), variables, p, xLR);

        double** dataMR = idMR.RealScalar_A3r(valatscri, in_guess1, in_guess2, Y, U, U0, gaussianCphi, gaussianA3r);
        double** dataLR = idLR.RealScalar_A3r(valatscri, in_guess1, in_guess2, Y, U, U0, gaussianCphi, gaussianA3r);

        System<double> sysHR(data  , dt    , x  , variables, sizet      , sizex      , p, StepInstantWrite);
        System<double> sysMR(dataMR, dt*f  , xMR, variables, sizet/f    , sizex/f    , p, StepInstantWrite/f);
        System<double> sysLR(dataLR, dt*f*f, xLR, variables, sizet/(f*f), sizex/(f*f), p, StepInstantWrite/(f*f));

        std::cout << "Partial L2 Convergence System Set" << std::endl;

        //Don't forget to change to L2PartialConvergence only, take out the 2
        te.L2PartialConvergence(sysHR, sysMR, sysLR, 4., 3., true, RHS);

        std::cout << "Partial L2 Convergence Calculated, use L2_Individual_Plot.py to see the output" << std::endl;
    }

    for (int i = 0; i < variables; i++){ delete[] data [i];}
    delete[] data    ;
    delete[] p       ;
    delete[] x       ;

    return 0;
}