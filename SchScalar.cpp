#include "TimeEvolution.h"
#include <matplotlibcpp.h>

quad cphigaussian(quad x) {
    quad sigma = 0.1;
    quad mu2   = 0.4;
    quad Acphi = 0.0001;
    
    return Acphi*exp(-(x*x - mu2*mu2) * (x*x - mu2*mu2) / (4. * sigma * sigma * sigma * sigma));
}

quad cphigaussian_derivative(quad x) {
    quad sigma = 0.1;
    quad mu2   = 0.4;
    quad Acphi = 0.0001;

    quad exponent = -(x*x - mu2*mu2) * (x*x - mu2*mu2) / (4. * sigma * sigma * sigma * sigma);
    quad exp_term = exp(exponent);

    quad term1 = -Acphi * (x*x - mu2*mu2) * x / (sigma * sigma * sigma * sigma) * exp_term;

    return term1;
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
    quad M    = BP.M   ;
    
    quad dx = rscri/sizex;
    
    // Find the index corresponding to the given x
    int index = static_cast<int>(round(2.0 * x / dx - 1.0));

    quad aconfi = aconf[index];
    
    quad dcphi = cphigaussian_derivative(x);
    quad cphi  = cphigaussian(x);

    result = (-3*pow(Ccmc,2)*pow(aconfi,4))/(4.*pow(x,6)*pow(Y,7)) - (pow(Kcmc,2)*Y)/(12.*pow(aconfi,2)) + (3*pow(Ccmc,2)*pow(aconfi,4)*Y)/(4.*pow(x,6)) - (pow(Kcmc,2)*M_PI*pow(x,2)*pow(cphi,2)*Y)/9. + (pow(Kcmc,2)*pow(Y,5))/(12.*pow(aconfi,2)) + (pow(Kcmc,2)*M_PI*x*cphi*Y*dcphi)/9. - (pow(Kcmc,2)*M_PI*pow(x,3)*cphi*Y*dcphi)/9. - (pow(Kcmc,2)*M_PI*Y*pow(dcphi,2))/36. + (pow(Kcmc,2)*M_PI*pow(x,2)*Y*pow(dcphi,2))/18. - (pow(Kcmc,2)*M_PI*pow(x,4)*Y*pow(dcphi,2))/36. - U/x - (sqrt(pow(Kcmc,2)*pow(x,6) + 9*pow(x,4)*pow(aconfi,2) + 6*(Ccmc*Kcmc - 3*M)*pow(x,3)*pow(aconfi,3) + 9*pow(Ccmc,2)*pow(aconfi,6))*U)/(3.*pow(x,3)*aconfi);
}

void U0(quad* data, quad x, int sizex, quad& result, BackgroundParameters<quad> BP){
    #ifdef DEBUG
        printf("[%s]\n", __PRETTY_FUNCTION__);
    #endif

    //Giving names
    quad Y    = data[0];

    quad cphi  = cphigaussian(0.0);
    quad dcphi = cphigaussian_derivative(0.0);

    result =  (- (pow(Kcmc,2)*M_PI*pow(x,2)*pow(cphi,2)*Y)/9. + (pow(Kcmc,2)*M_PI*x*cphi*Y*dcphi)/9. - (pow(Kcmc,2)*M_PI*pow(x,3)*cphi*Y*dcphi)/9. - (pow(Kcmc,2)*M_PI*Y*pow(dcphi,2))/36. + (pow(Kcmc,2)*M_PI*pow(x,2)*Y*pow(dcphi,2))/18. - (pow(Kcmc,2)*M_PI*pow(x,4)*Y*pow(dcphi,2))/36.)/2.0;
}

int main(int argc, char* argv[]){

    std::string Choice = "te"; //te, pc, l2. This chooses what type of run you want. Check the Makefile if you need

    // Check if there are command line arguments -> This is what we will use to run the different cases
    if(argc <= 1) std::cout << "No command line arguments provided, using default value for Time Evolution" << std::endl;
    else          Choice = std::string(argv[1]);

    int  variables(16);
    int  sizex    (540);
    int  sizet    (2000*StepInstantWrite);
    quad dx       (rscri/sizex);
    quad dt       (0.0002);

    std::cout << "CFL condition: " << dt/dx << std::endl;

    using precision = double;

    SpaceDerivator<quad> sder(Ghost_Left, Ghost_Right);
    
    quad M    =  1.0;
    quad QBH  =  0.0;    
    quad Ccmc =  sder.CcmcShootingMethod((quad) 0.01, (quad) 25.0, M, QBH);
    quad R0   =  sder.R0ShootingMethod  ((quad) 1.0  , M, QBH, Ccmc);

    sder.SetBackgroundParameters(M, QBH, Ccmc, R0);

    std::cout << "Background Parameters Set (Mass, Charge, Ccmc, R0)" << std::endl;

    void (*RHS[])(precision**, precision**, precision*, int, precision, int*, precision*, BackgroundParameters<precision>) = {ElectricField, PsiElectric, Beta, Alpha, Phi, A, cphi, dphi, cPi, dPi, trK, Gammarr, Chi, Arr, Lambdar, Theta};

    //Parity conditions
    int* p = new int[variables]{pE, pPsi, pBeta, pAlpha, pPhi, pA, pcphi, pdphi, pcPi, pdPi, ptrK, pGammarr,  pChi,  pArr, pLambdar, pTheta};

    std::cout << "RHS Set" << std::endl;

    ////////////////////////////////////////
    ///////////INITIAL CONDITIONS///////////
    ////////////////////////////////////////

    quad* x = sder.Grid(sizex, true);

    quad valatscri(1.0); //Value of the variable at scri (both aconf and chi should have this value)
    quad in_guess1(250.0); //First guess for aconf at x=dx
    quad in_guess2(-0.002); //Second guess for aconf at x=dx
    quad chi_guess1(1.0); //First guess for chi at x=dx/2
    quad chi_guess2(1.01); //Second guess for chi at x=dx/2

    InitialData         <quad> id(sizex, variables, p, x);
    BackgroundParameters<quad> BP = createBackgroundParameters(M, QBH, Ccmc, R0);

    quad** data = id.SchRealScalar(valatscri, in_guess1, in_guess2, Yaconf, valatscri, chi_guess1, chi_guess2, Y, U, U0, cphigaussian, cphigaussian_derivative, BP);

    std::cout << "Initial Conditions Set" << std::endl;

    ////////////////////////////////////
    ///////////TIME EVOLUTION///////////
    ////////////////////////////////////

    //We have to take care in changing from quad to long double or double
    //to make the simulation less expensive
    TimeEvolution<precision> te;

    if(Choice == "te"){

        SpaceDerivator<precision> sderld(Ghost_Left, Ghost_Right);

        //Converting to precision
        precision*  xld    = sderld.Grid(static_cast<precision>(sizex), true);
        precision** datald = ConvertData<quad, precision>(data, variables, sizex);

        BackgroundParameters<precision> BPld = ConvertData<quad, precision>(BP);

        System<precision> sys_ev(datald, static_cast<precision>(dt), xld, variables, static_cast<precision>(sizet), static_cast<precision>(sizex), p, StepInstantWrite, BPld);

        std::cout << "Time Evolution System Set" << std::endl;

        te.RungeKutta4(sys_ev, RHS);

        std::cout << "Time Evolution Complete, use Evolution Vid or Evolution Vid Binary python file" << std::endl;

        sys_ev.WriteBin("output/output3");     //Writing all of data
        sys_ev.WriteConstraint("output/Constraints"); //Writing Poisson constraint data
        sys_ev.WriteProperties("output/Properties"); //Writing Misner-Sharp Mass and Apparent Horizon Radius

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

        quad** dataMR = idMR.SchRealScalar(valatscri, in_guess1, in_guess2, Yaconf, valatscri, chi_guess1, chi_guess2, Y, U, U0, cphigaussian, cphigaussian_derivative, BP);
        quad** dataLR = idLR.SchRealScalar(valatscri, in_guess1, in_guess2, Yaconf, valatscri, chi_guess1, chi_guess2, Y, U, U0, cphigaussian, cphigaussian_derivative, BP);

        //Converting to precision
        SpaceDerivator<precision> sderld(Ghost_Left, Ghost_Right);

        precision*  xHRld = sderld.Grid(static_cast<precision>(sizex)      , true);
        precision*  xMRld = sderld.Grid(static_cast<precision>(sizex/f)    , true);
        precision*  xLRld = sderld.Grid(static_cast<precision>(sizex/(f*f)), true);

        precision** dataHRld = ConvertData<quad, precision>(data  , variables,       sizex);
        precision** dataMRld = ConvertData<quad, precision>(dataMR, variables,     sizex/f);
        precision** dataLRld = ConvertData<quad, precision>(dataLR, variables, sizex/(f*f));

        BackgroundParameters<precision> BPld = ConvertData<quad, precision>(BP);

        System<precision> sysHR(dataHRld, static_cast<precision>(dt)    , xHRld, variables, static_cast<precision>(sizet)      , static_cast<precision>(sizex)      , p, StepInstantWrite, BPld);
        System<precision> sysMR(dataMRld, static_cast<precision>(dt*f)  , xMRld, variables, static_cast<precision>(sizet/f)    , static_cast<precision>(sizex/f)    , p, StepInstantWrite/f, BPld);
        System<precision> sysLR(dataLRld, static_cast<precision>(dt*f*f), xLRld, variables, static_cast<precision>(sizet/(f*f)), static_cast<precision>(sizex/(f*f)), p, StepInstantWrite/(f*f), BPld);

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