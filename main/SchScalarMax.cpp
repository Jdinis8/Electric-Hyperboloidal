#include "TimeEvolution.h"
#include <matplotlibcpp.h>

quad cphigaussian(quad x) {
    quad sigma = 0.2;
    quad mu2   = 0.4;
    quad Acphi = 1.0;
    
    return Acphi*x*x*exp(-(x*x - mu2*mu2) * (x*x - mu2*mu2) / (4. * sigma * sigma * sigma * sigma));
}

quad cphigaussian_derivative(quad x) {
    quad sigma = 0.2;
    quad mu2   = 0.4;
    quad Acphi = 1.0;
    
    quad exponent = -(x*x - mu2*mu2) * (x*x - mu2*mu2) / (4. * sigma * sigma * sigma * sigma);
    quad exp_term = exp(exponent);

    quad term1 = 2 * Acphi * x * exp_term; // g'(x) * h(x)
    quad term2 = -Acphi * x * x * (x*x - mu2*mu2) * x / (sigma * sigma * sigma * sigma) * exp_term; // g(x) * h'(x)

    return term1 + term2;
}

quad gaussianA3r(quad x) {
    quad sigma = 0.2;
    quad mu2 = 0.25;
    quad A = 0.1;
    return ASw*A*x*x*exp(-(x*x - mu2*mu2) * (x*x - mu2*mu2) / (4. * sigma * sigma * sigma * sigma));
}

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

//For Chi
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
   
    quad Ccmc = BP.Ccmc;
    quad M    = BP.M;
    quad R0   = BP.R0;

    quad dx = rscri/sizex;
    
    // Find the index corresponding to the given x
    int index = static_cast<int>(round(2.0 * x / dx - 1.0));

    quad aconfi = aconf[index];
    
    quad dcphi = cphigaussian_derivative(x);
    quad cphi  = cphigaussian(x);
    quad A3    = gaussianA3r(x);
    
    if(x < 1e-6) result = (- (pow(Kcmc,2)*M_PI*pow(x,2)*pow(cphi,2)*Y)/9. - (pow(Kcmc,2)*M_PI*pow(q,2)*pow(A3,2)*pow(cphi,2)*Y)/36. + (pow(Kcmc,2)*M_PI*pow(q,2)*pow(x,2)*pow(A3,2)*pow(cphi,2)*Y)/18. - (pow(Kcmc,2)*M_PI*pow(q,2)*pow(x,4)*pow(A3,2)*pow(cphi,2)*Y)/36. + (pow(Kcmc,2)*M_PI*x*cphi*Y*dcphi)/9. - (pow(Kcmc,2)*M_PI*pow(x,3)*cphi*Y*dcphi)/9. - (pow(Kcmc,2)*M_PI*Y*pow(dcphi,2))/36. + (pow(Kcmc,2)*M_PI*pow(x,2)*Y*pow(dcphi,2))/18. - (pow(Kcmc,2)*M_PI*pow(x,4)*Y*pow(dcphi,2))/36.)/2.0;
    else result = (-3*pow(Ccmc,2)*pow(aconfi,4))/(4.*pow(x,6)*pow(Y,7)) - (pow(Kcmc,2)*Y)/(12.*pow(aconfi,2)) + (3*pow(Ccmc,2)*pow(aconfi,4)*Y)/(4.*pow(x,6)) - (pow(Kcmc,2)*M_PI*pow(x,2)*pow(cphi,2)*Y)/9. - (pow(Kcmc,2)*M_PI*pow(q,2)*pow(A3,2)*pow(cphi,2)*Y)/36. + (pow(Kcmc,2)*M_PI*pow(q,2)*pow(x,2)*pow(A3,2)*pow(cphi,2)*Y)/18. - (pow(Kcmc,2)*M_PI*pow(q,2)*pow(x,4)*pow(A3,2)*pow(cphi,2)*Y)/36. + (pow(Kcmc,2)*pow(Y,5))/(12.*pow(aconfi,2)) + (pow(Kcmc,2)*M_PI*x*cphi*Y*dcphi)/9. - (pow(Kcmc,2)*M_PI*pow(x,3)*cphi*Y*dcphi)/9. - (pow(Kcmc,2)*M_PI*Y*pow(dcphi,2))/36. + (pow(Kcmc,2)*M_PI*pow(x,2)*Y*pow(dcphi,2))/18. - (pow(Kcmc,2)*M_PI*pow(x,4)*Y*pow(dcphi,2))/36. - U/x - (sqrt(pow(Kcmc,2)*pow(x,6) + 9*pow(x,4)*pow(aconfi,2) + 6*(Ccmc*Kcmc - 3*M)*pow(x,3)*pow(aconfi,3) + 9*pow(Ccmc,2)*pow(aconfi,6))*U)/(3.*pow(x,3)*aconfi);
}

void U0(quad* data, quad x, int sizex, quad& result, BackgroundParameters<quad> BP){
    #ifdef DEBUG
        printf("[%s]\n", __PRETTY_FUNCTION__);
    #endif

    //Giving names
    quad Y    = data[0];

    quad Ccmc = BP.Ccmc;
    quad R0   = BP.R0  ;
    quad M    = BP.M   ;

    quad dx    = rscri/sizex;
    quad dcphi = cphigaussian_derivative(dx/2.0);
    quad cphi  = cphigaussian(dx/2.0);
    quad A3    = gaussianA3r(dx/2.0);

    result = (- (pow(Kcmc,2)*M_PI*pow(x,2)*pow(cphi,2)*Y)/9. - (pow(Kcmc,2)*M_PI*pow(q,2)*pow(A3,2)*pow(cphi,2)*Y)/36. + (pow(Kcmc,2)*M_PI*pow(q,2)*pow(x,2)*pow(A3,2)*pow(cphi,2)*Y)/18. - (pow(Kcmc,2)*M_PI*pow(q,2)*pow(x,4)*pow(A3,2)*pow(cphi,2)*Y)/36. + (pow(Kcmc,2)*M_PI*x*cphi*Y*dcphi)/9. - (pow(Kcmc,2)*M_PI*pow(x,3)*cphi*Y*dcphi)/9. - (pow(Kcmc,2)*M_PI*Y*pow(dcphi,2))/36. + (pow(Kcmc,2)*M_PI*pow(x,2)*Y*pow(dcphi,2))/18. - (pow(Kcmc,2)*M_PI*pow(x,4)*Y*pow(dcphi,2))/36.)/2.0;
}

int main(int argc, char* argv[]){

    std::string Choice = "te"; //te, pc, l2. This chooses what type of run you want. Check the Makefile if you need

    // Check if there are command line arguments -> This is what we will use to run the different cases
    if(argc <= 1) std::cout << "No command line arguments provided, using default value for Time Evolution" << std::endl;
    else          Choice = std::string(argv[1]);

    int  variables(16);
    int  sizex    (1080);
    int  sizet    (9*StepInstantWrite);
    quad dx       (rscri/sizex);
    quad dt       (0.00001);

    std::cout << "CFL condition: " << dt/dx << std::endl;

    SpaceDerivator<quad> sder(Ghost_Left, Ghost_Right);
    
    quad M    =  1.0;
    quad QBH  =  0.0;    
    quad Ccmc =  sder.CcmcShootingMethod((quad) 0.01, (quad) 25.0, M, QBH);
    quad R0   =  sder.R0ShootingMethod  ((quad) 1.0  , M, QBH, Ccmc);

    sder.SetBackgroundParameters(M, QBH, Ccmc, R0);

    std::cout << "Background Parameters Set (Mass, Charge, Ccmc, R0)" << std::endl;

    void (*RHS[])(long double**, long double**, long double*, int, long double, int*, long double*) = {ElectricField, PsiElectric, Beta, Alpha, Phi, A, cphi, dphi, cPi, dPi, trK, Gammarr, Chi, Arr, Lambdar, Theta};

    //Parity conditions
    int* p = new int[variables]{pE, pPsi, pBeta, pAlpha, pPhi, pA, pcphi, pdphi, pcPi, pdPi, ptrK, pGammarr,  pChi,  pArr, pLambdar, pTheta};

    std::cout << "RHS Set" << std::endl;

    ////////////////////////////////////////
    ///////////INITIAL CONDITIONS///////////
    ////////////////////////////////////////

    quad* x = sder.Grid(sizex, true);

    quad valatscri  =    1.0; //Value of the variable at scri (both aconf and chi should have this value)
    quad in_guess1  =  200.0; //First guess for aconf at x=dx
    quad in_guess2  = -0.002; //Second guess for aconf at x=dx
    quad chi_guess1 =    1.0; //First guess for chi at x=dx/2
    quad chi_guess2 =   1.01; //Second guess for chi at x=dx/2

    InitialData         <quad> id(sizex, variables, p, x);
    BackgroundParameters<quad> BP = createBackgroundParameters(M, QBH, Ccmc, R0);

    quad** data = id.SchRealScalarA3r(valatscri, in_guess1, in_guess2, Yaconf, valatscri, chi_guess1, chi_guess2, Y, U, U0, cphigaussian, cphigaussian_derivative, gaussianA3r, BP);

    std::cout << "Initial Conditions Set" << std::endl;

    ////////////////////////////////////
    ///////////TIME EVOLUTION///////////
    ////////////////////////////////////

    //We have to take care in changing from quad to long double or double
    //to make the simulation less expensive
    TimeEvolution<long double> te;

    if(Choice == "te"){

        SpaceDerivator<long double> sderld(Ghost_Left, Ghost_Right);

        //Converting to long double
        long double*  xld    = sderld.Grid(static_cast<long double>(sizex), true);
        long double** datald = new long double*[variables];

        for(int i = 0; i < variables; i++) datald[i] = new long double[sizex]();
        for(int i = 0; i < variables; i++) for(int j = 0; j < sizex; j++) datald[i][j] = static_cast<long double>(data[i][j]);

        BackgroundParameters<long double> BPld = ConvertData<quad, long double>(BP);

        System<long double> sys_ev(datald, static_cast<long double>(dt), xld, variables, static_cast<long double>(sizet), static_cast<long double>(sizex), p, StepInstantWrite, BPld);

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

        quad** dataMR = idMR.SchRealScalarA3r(valatscri, in_guess1, in_guess2, Yaconf, valatscri, chi_guess1, chi_guess2, Y, U, U0, cphigaussian, cphigaussian_derivative, gaussianA3r, BP);
        quad** dataLR = idLR.SchRealScalarA3r(valatscri, in_guess1, in_guess2, Yaconf, valatscri, chi_guess1, chi_guess2, Y, U, U0, cphigaussian, cphigaussian_derivative, gaussianA3r,BP);

        //Converting to long double
        SpaceDerivator<long double> sderld(Ghost_Left, Ghost_Right);

        long double*  xHRld = sderld.Grid(static_cast<long double>(sizex)      , true);
        long double*  xMRld = sderld.Grid(static_cast<long double>(sizex/f)    , true);
        long double*  xLRld = sderld.Grid(static_cast<long double>(sizex/(f*f)), true);

        long double** dataHRld = ConvertData<quad, long double>(data  , variables,       sizex);
        long double** dataMRld = ConvertData<quad, long double>(dataMR, variables,     sizex/f);
        long double** dataLRld = ConvertData<quad, long double>(dataLR, variables, sizex/(f*f));

        BackgroundParameters<long double> BPld = ConvertData<quad, long double>(BP);

        System<long double> sysHR(dataHRld, static_cast<long double>(dt)    , xHRld, variables, static_cast<long double>(sizet)      , static_cast<long double>(sizex)      , p, StepInstantWrite      , BPld);
        System<long double> sysMR(dataMRld, static_cast<long double>(dt*f)  , xMRld, variables, static_cast<long double>(sizet/f)    , static_cast<long double>(sizex/f)    , p, StepInstantWrite/f    , BPld);
        System<long double> sysLR(dataLRld, static_cast<long double>(dt*f*f), xLRld, variables, static_cast<long double>(sizet/(f*f)), static_cast<long double>(sizex/(f*f)), p, StepInstantWrite/(f*f), BPld);

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
