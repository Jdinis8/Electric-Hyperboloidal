#include "TimeEvolution.h"
#include <matplotlibcpp.h>

////////////////////////////////////////
////////Evolves A3r and Phi/////////////
////////////////////////////////////////

/*
TimeEvolutionF:            runs usual RK4 time evolution
PointConvergenceEvolution: runs pointwise convergence evolution. Poisson convergence is calculated here as well
PointConvergence:          runs pointwise convergence at the end of the simulation
PartialPointConvergence:   runs pointwise convergence for each variable, at the end of the simulation
*/

//For multithreading when writing
//Define a mutex to synchronize file access
std::mutex fileMutex;

double gaussianPhi(double x) {
    double sigma = 0.2;
    double mu2 = 0.25;
    double A = 0.5;
    return ASw*A*x*x*exp(-(x*x - mu2*mu2) * (x*x - mu2*mu2) / (4. * sigma * sigma * sigma * sigma));
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
    void (*RHS[])(double**, double**, double*, int, double, int*, double*) = {ElectricField, PsiElectric, Beta, Alpha, Phi, A, cphi, dphi, cPi, dPi, trK, Gammarr, Chi, Arr, Lambdar, Theta};

    //Parity conditions: vectors are odd, scalars are even
                               //E, Psi, Beta, Alpha, Phi, A, cphi, dphi, cPi, dPi, trK, Gammarr, Chi, Arr, Lambdar, Theta
    int* p = new int[variables]{-1,   1,   -1,     1,   1,-1,    1,    1,   1,   1,   1,       1,   1,  -1,      -1,     1};

    std::cout << "RHS Set" << std::endl;

    ////////////////////////////////////////
    ///////////INITIAL CONDITIONS///////////
    ////////////////////////////////////////
    SpaceDerivator<double> sder(Ghost_Left, Ghost_Right);

    double* x = sder.Grid(sizex, true);

    //Initial Conditions

    InitialData<double> id(sizex, variables, p, x);

    double** data = id.A3r_Phi(gaussianPhi);

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

        double** dataMR = idMR.A3r_Phi(gaussianPhi);
        double** dataLR = idLR.A3r_Phi(gaussianPhi);

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

        double** dataMR = idMR.A3r_Phi(gaussianPhi);
        double** dataLR = idLR.A3r_Phi(gaussianPhi);

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

        double** dataMR = idMR.A3r_Phi(gaussianPhi);
        double** dataLR = idLR.A3r_Phi(gaussianPhi);

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