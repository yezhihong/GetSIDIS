//This is a class for sidis generator or cross section calculation
//Tested with gcc v4.4.7
//LHAPDF is required. Developed with v6.1
//ROOT is required. Developed with v5.34.21
//Potential error with other versions not tested
//Last update date 6 Dec. 2017 version 2.2

#ifndef _LSIDIS_H_
#define _LSIDIS_H_

#include <iostream>
#include <fstream>
#include <cstring>
#include <cmath>

#include "LHAPDF/LHAPDF.h"

#include "TLorentzVector.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "TMatrixD.h"
#include "TH1D.h"

//using namespace LHAPDF;

class Lsidis{
    protected:/*{{{*/
        static constexpr double PI = 3.14159265358979323846264;//math Pi
        static constexpr double GeVfm = 1.0 / 0.1973269718;//GeV times fm in natural unit
        static constexpr double GeV2nb = 1.0 / 389379.338;//GeV^2 times nb in natural unit
        static constexpr double alpha_em = 1.0 / 137.035999074;//EM fine-structure constant
        static constexpr double eu2 = 4.0 / 9.0;//up type quark charge square
        static constexpr double ed2 = 1.0 / 9.0;//down type quark charge square
        static constexpr double Me = 0.000510998928;//electron mass in GeV
        static constexpr double Mmu = 0.1056583715;//muon mass in GeV
        static constexpr double Mp = 0.938272081;//proton mass in GeV
        static constexpr double Mn = 0.939565413;//neutron mass in GeV
        static constexpr double MLambda = 1.115683;//Lambda mass in GeV
        static constexpr double MSigmap = 1.18937;//Sigma+ mass in GeV
        static constexpr double MSigmam = 1.197449;//Sigma- mass in GeV
        static constexpr double Mpion = 0.13957018;//charged pion mass in GeV
        static constexpr double Mpi0 = 0.1349766;//neutral pion mass in GeV
        static constexpr double Mkaon = 0.493677;//charged kaon mass in GeV
        static constexpr double MK0 = 0.497614;//neutral kaon mass in GeV
        static constexpr double Meta = 0.77526;//rho mass in GeV
        static constexpr int NTERMS=96;//for NLO calculation ZYE 05/03/2022}}}

    private:/*{{{*/
        bool st_nucleus;
        bool st_pdf;
        bool st_ff;
        bool st_hadron;
        bool st_l;
        bool st_P;
        bool st_lp;
        bool st_Ph;
        bool st_gibbs;
        double Np;//number of protons in the nucleus
        double Nn;//number of neutrons in the nucleus
        double Polp;//effective polarization of the proton in the nucleus
        double Poln;//effective polarization of the neutron in the nucleus
        double Mh;//final hadron mass
        double MXminp;//minimum of the invariant mass of the final X from proton
        double MXminn;//minimum of the invariant mass of the final X from neutron
        char hname[10];//final hadron name
        int lid;//incoming lepton id
        int lic;//incoming lepton charge
        int hid;//final hadron id
        int hic;//final hadron charge
        LHAPDF::PDF * pdfs;//pointer to collinear unpolarized PDFs
        LHAPDF::PDF * ffs;//pointer to collinear unpolarized FFs
        double XI[NTERMS];//for NLO calculation, ZYE 05/03/2022
        double WI[NTERMS];//for NLO calculation, ZYE 05/03/2022
        double XX[NTERMS+1];//for NLO calculation, ZYE 05/03/2022
        TLorentzVector Pl;//incoming lepton 
        TLorentzVector Plp;//outgoing lepton
        TLorentzVector Pq;//virtual boson
        TLorentzVector PP;//incoming nucleus
        TLorentzVector PPh;//outgoing hadron
        TLorentzVector PX;//outgoing undetected particles
        double Slepton;//incoming lepton polarization rate
        double SNL;//incoming nucleus longitudinal polarization rate
        double SNT;//incoming nucleus transverse polarization rate
        double x;//Bjorken scaling variable x
        double xn;//Nachtmann scaling variable xn
        double y;//lepton energy transferred fraction scalar variable
        double z;//fragmentation scalar variable
        double Q2;//transferred momentum square, also used as factorization scale
        double Pt;//transverse momentum of final hadron in Trento convention
        double phih;//azimuthal angle of final hadron in Trento convention
        double phiS;//azimuthal angle of transverse polarization in Trento convention
        double W;//invariant mass of final state except scattered lepton
        double Wp;//invariant mass of final state except scattered lepton and the hadron
        double gamma;//
        double epsilon;//longi-trans-photon flux ratio
        double Rfactor;//a criterism for current fragmentation judgement
        double jacobian;//Jacobian from simulation space to cross section defined space
        bool physics_control;//
        double f1[9];//PDFs of u, d, s, ubar, dbar, sbar, ccbar,bbbar, gluon, added last three by ZYE 05/03/2022
        double D1[9];//FFs of u, d, s, ubar, dbar, sbar, ccbar,bbbar, gluon, added last three by ZYE 05/03/2022
        double g1[6];//Polarized PDF of u, d, s, ubar, dbar, sbar
        double h1[6];//Transversity PDF of u, d, s, ubar, dbar, sbar
        double TMDpars[2];//Model parameters of gaussian type TMDs
        double Xmin[6];//simulation variables lower limits
        double Xmax[6];//simulation variables upper limits
        double Xseed[6];//start point for MCMC
        double volume;//simulation variables volume
        double sigmatotal;//total cross section within the simulation range [Xmin, Xmax]
        TH1D * Xhisto[6];//for Gibbs sampler}}}
    
    public:/*{{{*/
        Lsidis();
        Lsidis(const TLorentzVector l, const TLorentzVector P);
        static double DEG(const double rad);//convert from rad to deg
        static double RAD(const double deg);//convert form deg to rad  
        int ChangeTMDpars(const double kt2, const double pt2);//change the TMD model parameters TMDpars
        int SetNucleus(const double np, const double nn, const double pp, const double pn);//Set incoming nucleus
        int CheckNucleus();
        int SetHadron(const char * hadron);//Set detected hadron
        int GetHadronID();//Get pid of the detected hadron
        int GetHadronCharge();//Get charge of the detected hadron
        int CheckHadron();
        int SetPDFset(const char * set);//Set unpolarized collinear PDF set
        int SetFFset(const char * set);//Set unpolarized collinear FF set
        int CheckPDFset();
        int CheckFFset();
        int SetInitialState(const TLorentzVector l, const TLorentzVector P);//Set incoming lepton and nucleon 4-momenta
        int SetLeptonBeam(const TLorentzVector l);//Set incoming lepton 4-momentum
        int SetLeptonBeam(const double Ebeam, const double m);//Set incoming lepton along z-direction
        int SetIonBeam(const TLorentzVector P);//Set incoming nucleon 4-momentum
        int SetIonBeam(const double Ebeam, const double M);//Set incoming nucleon along negative z-direction
        int SetTarget(const double np, const double nn, const double pp, const double pn);//Set fixed target
        int SetTarget();//Set fixed target
        int SetFinalState(const TLorentzVector lp, const TLorentzVector Ph);//Set final lepton and detected hadron 4-momenta
        int SetFinalLepton(const TLorentzVector lp);//Set final scattered lepton 4-momenta
        int SetFinalHadron(const TLorentzVector Ph);//Set final detected hadron 4-momenta
        int CalculateVariables();//Calculate Lorentz scalar variables
        int CalculateFinalState();//Calculate scattered lepton and detected hadron from x, y, z, Pt, phih, phiS
        int CalculateRfactor(const double kT2, const double MiT2, const double MfT2);//Calculate the Rfactor for current fragmentation criteria
        int SetVariables(const double x0, const double y0, const double z0, const double Pt0, const double phih0, const double phiS0);//Set variables x, y, z, Pt2, phih, phiS
        double GetVariable(const char * var);//Get particular variable of current event
        TLorentzVector GetLorentzVector(const char * part);//Get 4-momentum of a particle of current event
        int GetPDFs();//Get PDFs at current x, Q2
        int GetFFs();//Get FFs at current z, Q2
        double FUUT();//Calculate structure function F_UU,T at current x, z, Q2, Pt
        double dsigma(const int mode);//Differential cross section in dsigma/dx dy dz dPt2 dphih dphiS
        int SetRange(const double * mins, const double * maxs);//Uniform sampling range
        double GenerateEvent(const int mode, const int method);//Generate a sidis event
        double GetEventWeight(const int mode, const int method);//Get event weight
        int CalculateSigmaTotal(const int mode, const int method, const int precision);//Calculate total cross section within the simulation range [Xmin, Xmax]
        int GibbsStarter(const int mode, const int method);
        double GibbsSampler(const int mode, const int method);//Generate an event with Gibbs sampler
        int WATE96();
        double FUUT_NLO(bool bNLO, bool bNeutron);//calculate the NLO structure function F_UU,T
        int CheckIsPhy();
        int Test();
};/*}}}*/

Lsidis::Lsidis(){//constructor{{{
    st_nucleus = false;
    st_hadron = false;
    st_pdf = false;
    st_ff = false;
    st_l = false;
    st_P = false;
    st_lp = false;
    st_Ph = false;
    st_gibbs = false;
    physics_control = false;
    TMDpars[0] = 0.54;//kt2
    TMDpars[1] = 0.13;//pt2
    Slepton = 0;
    SNL = 0;
    SNT = 0;
    //for NLO calculation, ZYE 05/03/2022
    WATE96();
}/*}}}*/

Lsidis::Lsidis(const TLorentzVector l, const TLorentzVector P){//constructor with initial state{{{
    Lsidis();
    SetInitialState(l, P);
}/*}}}*/

double Lsidis::DEG(const double rad){//convert rad to deg{{{
    return rad / PI * 180.0;
}/*}}}*/

double Lsidis::RAD(const double deg){//convert deg to rad{{{
    return deg / 180.0 * PI;
}/*}}}*/

int Lsidis::ChangeTMDpars(const double kt2, const double pt2){//change TMD model parameters{{{
    TMDpars[0] = kt2;
    TMDpars[1] = pt2;
    return 0;
}/*}}}*/

int Lsidis::SetNucleus(const double np, const double nn, const double pp = 0, const double pn = 0){//set number of protons and neutrons in the incoming nucleus{{{
    Np = np;
    Nn = nn;
    Polp = pp;
    Poln = pn;
    st_nucleus = true;
    return 0;
}/*}}}*/

int Lsidis::CheckNucleus(){//check the status of nucleus{{{
    if (st_nucleus){
        printf("Number of protons:  %.1f\n", Np);
        printf("Number of neutrons: %.1f\n", Nn);
        return 0;
    }
    else {
        printf("No nucleus initialized!\n");
        return -1;
    }
}/*}}}*/

int Lsidis::CheckIsPhy(){/*{{{*/
   if(physics_control) return 1;
   else return 0;
}/*}}}*/

int Lsidis::SetHadron(const char * hadron){//set final hadron{{{
    if (strcmp(hadron, "pi+") == 0){
        hic = 1; hid = 211;
        Mh = Mpion;
        MXminp = Mp;
        MXminn = Mp + Mpion;
    }
    else if (strcmp(hadron, "pi-") == 0){
        hic = -1; hid = -211;
        Mh = Mpion;
        MXminp = Mp + Mpion;
        MXminn = Mp;
    }
    else if (strcmp(hadron, "pi0") == 0){
        hic = 0; hid = 111;
        Mh = Mpi0;
        MXminp = Mp;
        MXminn = Mp;
    }
    else if (strcmp(hadron, "K+") == 0){
        hic = 1; hid = 321;
        Mh = Mkaon;
        MXminp = MLambda;
        MXminn = MSigmam;
    }
    else if (strcmp(hadron, "K-") == 0){
        hic = -1; hid = -321;
        Mh = Mkaon;
        MXminp = Mp + Mkaon;
        MXminn = Mp + Mkaon;
    }
    else if (strcmp(hadron, "K0") == 0){
        hic = 0; hid = 310;//Ks
        Mh = MK0;
        MXminp = MSigmap;
        MXminn = MLambda;
    }
    else if (strcmp(hadron, "p") == 0){
        hic = 1; hid = 2212;
        Mh = Mp;
        MXminp = Mpi0;
        MXminn = Mpion;
    }
    else {
        perror("Invalid hadron option in SetHadron!");
        return -1;
    }
    strcpy(hname, hadron);
    st_hadron = true;
    return 0;
}/*}}}*/

int Lsidis::GetHadronID(){//get pid of the detected hadron{{{
    return hid;
}/*}}}*/

int Lsidis::GetHadronCharge(){//get charge of the detected hadron{{{
    return hic;
}/*}}}*/

int Lsidis::CheckHadron(){//check the status of final hadron{{{
    if (st_hadron){
        printf("final hadron: %s\n", hname);
        printf("charge: %d,\t pid: %d\n", hic, hid);
        return 0;
    }
    else {
        perror("No hadron initialized!");
        return -1;
    }
}/*}}}*/

int Lsidis::SetPDFset(const char * set){//set collinear unpolarized pdf set{{{
    pdfs = LHAPDF::mkPDF(set, 0);
    st_pdf = true;
    return 0;
}/*}}}*/

int Lsidis::SetFFset(const char * set = "DSSFFlo"){//set collinear unpolarized ff set{{{
    if (!st_hadron){
        perror("No hadron initialized!");
        return -1;
    }
    if (hid >= 0)
        ffs = LHAPDF::mkPDF(set, hid);
    else
        ffs = LHAPDF::mkPDF(set, -hid+1000);
    st_ff = true;
    return 0;
}   /*}}}*/

int Lsidis::CheckPDFset(){//check the status of collinear unpolarized pdf set{{{
    if (st_pdf){
        pdfs->print(std::cout, 1);
        printf("x range:  %.4E ~ %.4E\n", pdfs->xMin(), pdfs->xMax());
        printf("Q2 range: %.4E ~ %.4E GeV^2\n", pdfs->q2Min(), pdfs->q2Max());
        return 0;
    }
    else {
        printf("No pdf set initialized!\n");
        return -1;
    }
}/*}}}*/

int Lsidis::SetInitialState(const TLorentzVector l, const TLorentzVector P){//set incoming lepton and nucleon 4-momentums{{{
    Pl = l;
    PP = P;
    st_l = true;
    st_P = true;
    return 0;
}/*}}}*/

int Lsidis::SetLeptonBeam(const TLorentzVector l){//set incoming lepton 4-momentum{{{
    Pl = l;
    st_l = true;
    return 0;
}/*}}}*/

int Lsidis::SetLeptonBeam(const double Ebeam, const double m = Me){//set incoming lepton with beam energy Ebeam along z-direction{{{
    Pl.SetXYZM(0, 0, sqrt(Ebeam * Ebeam - m * m), m);
    st_l = true;
    return 0;
}/*}}}*/

int Lsidis::SetIonBeam(const TLorentzVector P){//set incoming nucleon 4-momentum{{{
    PP = P;
    st_P = true;
    return 0;
}/*}}}*/

int Lsidis::SetIonBeam(const double Ebeam, const double M = Mp){//set incoming nucleon with beam energy Ebeam along negative z-direction{{{
    PP.SetXYZM(0, 0, -sqrt(Ebeam * Ebeam - M * M), Ebeam);
    st_P = true;
    return 0;
}/*}}}*/

int Lsidis::SetTarget(const double np, const double nn, const double pp, const double pn){//set fixed target{{{
    SetNucleus(np, nn, pp, pn);
    PP.SetXYZM(0, 0, 0, Mp);
    return 0;
}/*}}}*/

int Lsidis::SetTarget(){//set fixed target{{{
    if (st_nucleus){
        PP.SetXYZM(0, 0, 0, Mp);
        return 0;
    }
    else {
        perror("No nucleus initialized!");
        return -1;
    }
}/*}}}*/

int Lsidis::SetFinalState(const TLorentzVector lp, const TLorentzVector Ph){//set final state{{{
    SetFinalLepton(lp);
    SetFinalHadron(Ph);
    return 0;
}/*}}}*/

int Lsidis::SetFinalLepton(const TLorentzVector lp){//set scattered lepton{{{
    Plp = lp;
    st_lp = true;
    return 0;
}/*}}}*/

int Lsidis::SetFinalHadron(const TLorentzVector Ph){//set final hadron{{{
    PPh = Ph;
    st_Ph = true;
    return 0;
}/*}}}*/

int Lsidis::CalculateVariables(){//calculate lorentz scalar variables{{{
    if (st_l && st_P && st_lp && st_Ph && st_nucleus && st_hadron){
        Pq = Pl - Plp;//calculate virtual photon 4-momentum
        PX = Pq + PP - PPh;//calculate final X 4-momentum
        Q2 = - (Pq * Pq);
        if (Q2 < Mpion * Mpion){//Q2 below a cut
            physics_control = false;
            return 0;
        }
        double W2 = Mp * Mp + 2.0 * PP * Pq - Q2;
        if (W2 < pow(Mp + Mpi0, 2)){//below the lowest threshold
            physics_control = false;
            return 0;
        }
        W = sqrt(W2);
        double Wp2 = PX * PX;
        if (Wp2 < Mpi0 * Mpi0){//below the lowest threshold
            physics_control = false;
            return 0;
        }
        Wp = sqrt(Wp2);
        x = Q2 / (PP * Pq) / 2.0;
        if (x > 1.0 || x < 1e-6){//x below a cut
            physics_control = false;
            return 0;
        }
        y = (PP * Pq) / (PP * Pl);
        if (y > 1.0 || y < 0.0){
            physics_control = false;
            return 0;
        }
        z = (PP * PPh) / (PP * Pq);
        if (z > 1.0 || z < 0.05){//z below a cut
            physics_control = false;
            return 0;
        }
        gamma = 2.0 * x * Mp / sqrt(Q2);
        epsilon = (1.0 - y - 0.25 * gamma * gamma * y * y) / (1.0 - y + 0.5 * y * y + 0.25 * gamma * gamma * y * y);
        xn = 2.0 * x / (1.0 + sqrt(1.0 + gamma * gamma));
        double Pt2 = - (PPh * PPh) + 2.0 * (PP * PPh) * (Pq * PPh) / (PP * Pq) / (1.0 + gamma * gamma) - gamma * gamma / (1.0 + gamma * gamma) * (pow(Pq * PPh, 2) / Q2 - pow(PP * PPh, 2) / (Mp * Mp));
        if (Pt2 < 0) Pt2 = 0.0;
        Pt = sqrt(Pt2);
        double lt2 = - (Pl * Pl) + 2.0 * (PP * Pl) * (Pq * Pl) / (PP * Pq) / (1.0 + gamma * gamma) - gamma * gamma / (1.0 + gamma * gamma) * (pow(Pq * Pl, 2) / Q2 - pow(PP * Pl, 2) / (Mp * Mp));
        double ch = - 1.0 / sqrt(lt2 * Pt * Pt) * ( (Pl * PPh) - ((Pq * Pl) * (PP * PPh) + (PP * Pl) * (Pq * PPh)) / (1.0 + gamma * gamma) / (PP * Pq) + gamma * gamma / (1.0 + gamma * gamma) * ( (Pq * Pl) * (Pq * PPh) / Q2 - (PP * Pl) * (PP * PPh) / (Mp * Mp)));
        TMatrixD eg5(4,4);
        eg5(0,0) = Pl.E(); eg5(0,1) = -Pl.X(); eg5(0,2) = -Pl.Y(); eg5(0,3) = -Pl.Z();
        eg5(1,0) = PPh.E(); eg5(1,1) = -PPh.X(); eg5(1,2) = -PPh.Y(); eg5(1,3) = -PPh.Z();
        eg5(2,0) = PP.E(); eg5(2,1) = -PP.X(); eg5(2,2) = -PP.Y(); eg5(2,3) = -PP.Z();
        eg5(3,0) = Pq.E(); eg5(3,1) = -Pq.X(); eg5(3,2) = -Pq.Y(); eg5(3,3) = -Pq.Z();
        double sh = -1.0 / sqrt(lt2 * Pt * Pt) * eg5.Determinant() / ((PP * Pq) * sqrt(1.0 + gamma * gamma));
        if (sh > 0) phih = acos(ch);
        else phih = -acos(ch);
        physics_control = true;
        GetPDFs();
        GetFFs();
        return 0;
    }
    else {
        physics_control = false;
        perror("Initial or final state missing!");
        return -1;
    }
}/*}}}*/

int Lsidis::CalculateFinalState(){//Calculate scattered electron and detected hadron from x, y, z, Pt, phih, phiS{{{
    if (st_l && st_P && st_nucleus && st_hadron){
        if (x > 1.0 || x < 1e-6){//x below a cut
            physics_control = false;
            return 0;
        }
        if (y > 1.0 || y < 0.0){//y
            physics_control = false;
            return 0;
        }
        if (z > 1.0 || z < 0.05){//z below a cut
            physics_control = false;
            return 0;
        }
        if (Pt < 0.0){//unphysical Pt
            physics_control = false;
            return 0;
        }
        Q2 = 2.0 * x * y * (PP * Pl);
        if (Q2 < Mpion * Mpion){//Q2 below a cut
            physics_control = false;
            return 0;
        }
        TLorentzVector Pl_2 = Pl;
        Pl_2.Boost(-PP.BoostVector());
        TLorentzVector Pl_1 = Pl_2;
        Pl_1.RotateZ(-Pl_2.Phi());
        Pl_1.RotateY(-Pl_2.Theta());
        double Elp_1 = Pl_1.E() - Q2 / (2.0 * x * Mp);
        double theta_lp_1 = acos(1.0 - Q2 / (2.0 * Pl_1.E() * Elp_1));
        TLorentzVector Plp_1(Elp_1 * sin(theta_lp_1), 0, Elp_1 * cos(theta_lp_1), Elp_1);
        TLorentzVector Pq_1 = Pl_1 - Plp_1;
        double shl = -cos(Pq_1.Theta()) * sin(phiS) / sqrt(1.0 - pow(sin(Pq_1.Theta()) * sin (phiS), 2));
        double chl = cos(phiS) / sqrt(1.0 - pow(sin(Pq_1.Theta()) * sin (phiS), 2));
        double phil;
        if (shl > 0) phil = acos(chl);
        else phil = -acos(chl);
        Plp_1.RotateZ(phil);
        Pq_1.RotateZ(phil);
        Plp = Plp_1;
        Plp.RotateY(Pl_2.Theta());
        Plp.RotateZ(Pl_2.Phi());
        Plp.Boost(PP.BoostVector());
        Pq = Pl - Plp;
        if (z < sqrt(Mh * Mh + Pt * Pt) / Pq_1.E()){//below hadron threshold
            physics_control = false;
            return 0;
        }
        TLorentzVector PPh_0;
        TLorentzVector Plp_0 = Plp_1;
        Plp_0.RotateZ(-Pq_1.Phi());
        Plp_0.RotateY(-Pq_1.Theta());
        PPh_0.SetXYZT(Pt * cos(phih + Plp_0.Phi()), Pt * sin(phih + Plp_0.Phi()), sqrt(pow(z * Pq_1.E(), 2) - Mh * Mh - Pt * Pt), z * Pq_1.E());
        TLorentzVector PPh_1 = PPh_0;
        PPh_1.RotateY(Pq_1.Theta());
        PPh_1.RotateZ(Pq_1.Phi());
        PPh = PPh_1;
        PPh.RotateY(Pl_2.Theta());
        PPh.RotateZ(Pl_2.Phi());
        PPh.Boost(PP.BoostVector());
        //CalculateVariables();
        double W2 = Mp * Mp + 2.0 * (PP * Pq) - Q2;
        if (W2 < pow(Mp + Mpi0, 2)){//below the lowest threshold
            physics_control = false;
            return 0;
        }
        W = sqrt(W2);
        PX = Pq + PP - PPh;
        double Wp2 = PX * PX;
        if (Wp2 < Mpi0 * Mpi0){//below the lowest threshold
            physics_control = false;
            return 0;
        }
        Wp = sqrt(Wp2);
        gamma = 2.0 * x * Mp / sqrt(Q2);
        epsilon = (1.0 - y - 0.25 * gamma * gamma * y * y) / (1.0 - y + 0.5 * y * y + 0.25 * gamma * gamma * y * y);
        xn = 2.0 * x / (1.0 + sqrt(1.0 + gamma * gamma));
        physics_control = true;
        GetPDFs();
        GetFFs();
        return 0;
    }
    else {
        physics_control = false;
        perror("Initial state missing for final state calculation!");
        return -1;
    }
}/*}}}*/

int Lsidis::CalculateRfactor(const double kT2 = 0.5, const double MiT2 = 0.5, const double MfT2 = 0.5){/*{{{*/
    if (physics_control){
        double yi = 0.5 * log(Q2 / MiT2);
        double yf = -0.5 * log(Q2 / MfT2);
        double yh = log( (sqrt(Q2) * z * (Q2 - xn * xn * Mp * Mp)) / ( 2.0 * xn * xn * Mp * Mp * sqrt(Mh * Mh + Pt * Pt)) - sqrt(Q2) / (xn * Mp) * sqrt(pow(z * (Q2 - xn * xn * Mp * Mp), 2) / (4.0 * xn * xn  * Mp * Mp * (Mh * Mh + Pt * Pt)) - 1.0));
        double Rf = 0.5 * sqrt(Pt * Pt + Mh * Mh) * sqrt(MfT2) * (exp(yf - yh) + exp(yh - yf)) - sqrt(kT2) * Pt;
        double Ri = 0.5 * sqrt(Pt * Pt + Mh * Mh) * sqrt(MiT2) * (exp(yi - yh) - exp(yh - yi)) - sqrt(kT2) * Pt;
        Rfactor = std::abs(Rf / Ri);
        return 0;
    }
    else {
        Rfactor = 1. / 0.;//NaN
        return -1;
    }
}/*}}}*/

double Lsidis::GetVariable(const char * var){//get current variable{{{
    if (!physics_control){
        std::cout << "Warning: unphysical event!" << std::endl;
        return -1;
    }
    if (strcmp(var, "x") == 0) return x;
    else if (strcmp(var, "y") == 0) return y;
    else if (strcmp(var, "z") == 0) return z;
    else if (strcmp(var, "Q2") == 0) return Q2;
    else if (strcmp(var, "Pt") == 0) return Pt;
    else if (strcmp(var, "phih") == 0) return phih;
    else if (strcmp(var, "phiS") == 0) return phiS;
    else if (strcmp(var, "W") == 0) return W;
    else if (strcmp(var, "Wp") == 0) return Wp;
    else if (strcmp(var, "xn") == 0) return xn;
    else if (strcmp(var, "gamma") == 0) return gamma;
    else if (strcmp(var, "epsilon") == 0) return epsilon;
    else if (strcmp(var, "Rfactor") == 0) return Rfactor;
    else if (strcmp(var, "sigmatotal") == 0) return sigmatotal;
    else {
        perror("No such variable!");
        return -1;
    }
}/*}}}*/

int Lsidis::SetVariables(const double x0, const double y0, const double z0, const double Pt0, const double phih0, const double phiS0){//Set variables x, y, z, Pt2, phih, phiS{{{
    x = x0;
    y = y0;
    z = z0;
    Pt = Pt0;
    phih = phih0;
    phiS = phiS0;
    st_lp = true;
    st_Ph = true;
    return 0;
}/*}}}*/

TLorentzVector Lsidis::GetLorentzVector(const char * part){//Get 4-momentum of a particle in current event{{{
    if (!physics_control){
        std::cout << "Warning: unphysical event!" << std::endl;
        TLorentzVector tt(0,0,0,0);
        return tt;
    }
    if (strcmp(part, "P") == 0) return PP;
    else if (strcmp(part, "l") == 0) return Pl;
    else if (strcmp(part, "lp") == 0) return Plp;
    else if (strcmp(part, "Ph") == 0) return PPh;
    else {
        perror("No such variable!");
        TLorentzVector tt(0,0,0,0);
        return tt;
    }
}/*}}}*/

int Lsidis::GetPDFs(){//get PDFs of proton at current x, Q2{{{
    if (!st_pdf){
        perror("No PDF set initialized!");
        return -1;
    }
    if (physics_control){
        f1[0] = pdfs->xfxQ2(2, x, Q2) / x;
        f1[1] = pdfs->xfxQ2(1, x, Q2) / x;
        f1[2] = pdfs->xfxQ2(3, x, Q2) / x;
        f1[3] = pdfs->xfxQ2(-2, x, Q2) / x;
        f1[4] = pdfs->xfxQ2(-1, x, Q2) / x;
        f1[5] = pdfs->xfxQ2(-3, x, Q2) / x;
        f1[6] =(pdfs->xfxQ2(4, x, Q2)+pdfs->xfxQ2(-4, x, Q2)) / x;
        f1[7] =(pdfs->xfxQ2(5, x, Q2)+pdfs->xfxQ2(-5, x, Q2)) / x;
        f1[8] = pdfs->xfxQ2(0, x, Q2) / x;
    }
    else {
        for (int i = 0; i < 9; i++)
            f1[i] = 0;
    }
    return 0;
}/*}}}*/

int Lsidis::GetFFs(){//get FFs to hadron at current z, Q2{{{
    if (!st_ff){
        perror("No FF set initialized!");
        return -1;
    }
    if (physics_control){
        D1[0] = ffs->xfxQ2(2, z, Q2) / z;//u
        D1[1] = ffs->xfxQ2(1, z, Q2) / z;//d
        D1[2] = ffs->xfxQ2(3, z, Q2) / z;//s
        D1[3] = ffs->xfxQ2(-2, z, Q2) / z;//ubar
        D1[4] = ffs->xfxQ2(-1, z, Q2) / z;//dbar
        D1[5] = ffs->xfxQ2(-3, z, Q2) / z;//sbar
    }
    else {
        for (int i = 0; i < 6; i++)
            D1[i] = 0;
    }
    return 0;
}/*}}}*/

int Lsidis::WATE96(){/*{{{*/
    /*
       IMPLICIT DOUBLE PRECISION (A-H, O-Z) ! G.W. 15/02/2007
       CC**********************************************************************
       CC*****							                                   *****
       CC***** THE X[i]	AND W[i] ARE THE DIRECT	OUTPUT FROM A PROGRAM  *****
       CC***** USING NAG ROUTINE D01BCF	TO CALCULATE THE	           *****
       CC***** GAUSS-LEGENDRE WEIGHTS FOR 96 POINT INTEGRATION.	       *****
       CC***** THEY AGREE TO TYPICALLY 14 DECIMAL PLACES WITH THE         *****
       CC***** TABLE IN	ABRAMOWITZ & STEGUN, PAGE 919.		           *****
       CC*****							                                   *****
       CC***** ---->   PETER HARRIMAN, APRIL 3RD 1990.	        	       *****
       CC*****							                                   *****
       CC**********************************************************************
       */
    double X[48],W[48];
    /*{{{*/
    X[ 1-1]=   0.01627674484960183561;
    X[ 2-1]=   0.04881298513604856015;
    X[ 3-1]=   0.08129749546442434360;
    X[ 4-1]=   0.11369585011066471632;
    X[ 5-1]=   0.14597371465489567682;
    X[ 6-1]=   0.17809688236761733390;
    X[ 7-1]=   0.21003131046056591064;
    X[ 8-1]=   0.24174315616383866556;
    X[ 9-1]=   0.27319881259104774468;
    X[10-1]=   0.30436494435449495954;
    X[11-1]=   0.33520852289262397655;
    X[12-1]=   0.36569686147231213885;
    X[13-1]=   0.39579764982890709712;
    X[14-1]=   0.42547898840729897474;
    X[15-1]=   0.45470942216774136446;
    X[16-1]=   0.48345797392059470382;
    X[17-1]=   0.51169417715466604391;
    X[18-1]=   0.53938810832435567233;
    X[19-1]=   0.56651041856139533470;
    X[20-1]=   0.59303236477757022282;
    X[21-1]=   0.61892584012546672523;
    X[22-1]=   0.64416340378496526886;
    X[23-1]=   0.66871831004391424358;
    X[24-1]=   0.69256453664216964528;
    X[25-1]=   0.71567681234896561582;
    X[26-1]=   0.73803064374439816819;
    X[27-1]=   0.75960234117664555964;
    X[28-1]=   0.78036904386743123629;
    X[29-1]=   0.80030874413913884180;
    X[30-1]=   0.81940031073792957139;
    X[31-1]=   0.83762351122818502758;
    X[32-1]=   0.85495903343459936363;
    X[33-1]=   0.87138850590929436968;
    X[34-1]=   0.88689451740241818933;
    X[35-1]=   0.90146063531585023110;
    X[36-1]=   0.91507142312089592706;
    X[37-1]=   0.92771245672230655266;
    X[38-1]=   0.93937033975275308073;
    X[39-1]=   0.95003271778443564022;
    X[40-1]=   0.95968829144874048809;
    X[41-1]=   0.96832682846326217918;
    X[42-1]=   0.97593917458513455843;
    X[43-1]=   0.98251726356301274934;
    X[44-1]=   0.98805412632962202890;
    X[45-1]=   0.99254390032376081654;
    X[46-1]=   0.99598184298720747465;
    X[47-1]=   0.99836437586317963722;
    X[48-1]=   0.99968950388322870559;
    W[ 1-1]=   0.03255061449236316962;
    W[ 2-1]=   0.03251611871386883307;
    W[ 3-1]=   0.03244716371406427668;
    W[ 4-1]=   0.03234382256857594104;
    W[ 5-1]=   0.03220620479403026124;
    W[ 6-1]=   0.03203445623199267876;
    W[ 7-1]=   0.03182875889441101874;
    W[ 8-1]=   0.03158933077072719007;
    W[ 9-1]=   0.03131642559686137819;
    W[10-1]=   0.03101033258631386231;
    W[11-1]=   0.03067137612366917839;
    W[12-1]=   0.03029991542082762553;
    W[13-1]=   0.02989634413632842385;
    W[14-1]=   0.02946108995816795100;
    W[15-1]=   0.02899461415055528410;
    W[16-1]=   0.02849741106508543861;
    W[17-1]=   0.02797000761684838950;
    W[18-1]=   0.02741296272602931385;
    W[19-1]=   0.02682686672559184485;
    W[20-1]=   0.02621234073567250055;
    W[21-1]=   0.02557003600534944960;
    W[22-1]=   0.02490063322248370695;
    W[23-1]=   0.02420484179236479915;
    W[24-1]=   0.02348339908592633665;
    W[25-1]=   0.02273706965832950717;
    W[26-1]=   0.02196664443874448477;
    W[27-1]=   0.02117293989219144572;
    W[28-1]=   0.02035679715433347898;
    W[29-1]=   0.01951908114014518992;
    W[30-1]=   0.01866067962741165898;
    W[31-1]=   0.01778250231604547316;
    W[32-1]=   0.01688547986424539715;
    W[33-1]=   0.01597056290256253144;
    W[34-1]=   0.01503872102699521608;
    W[35-1]=   0.01409094177231515264;
    W[36-1]=   0.01312822956696188190;
    W[37-1]=   0.01215160467108866759;
    W[38-1]=   0.01116210209983888144;
    W[39-1]=   0.01016077053500880978;
    W[40-1]=   0.00914867123078384552;
    W[41-1]=   0.00812687692569928101;
    W[42-1]=   0.00709647079115442616;
    W[43-1]=   0.00605854550423662775;
    W[44-1]=   0.00501420274292825661;
    W[45-1]=   0.00396455433844564804;
    W[46-1]=   0.00291073181793626202;
    W[47-1]=   0.00185396078894924657;
    W[48-1]=   0.00079679206555731759;
    /*}}}*/

    //double XI[96],WI[96], XX[97];
    for(int i=0; i<48;i++){
        XI[i]=-X[47-i];
        WI[i]=W[47-i];
        XI[i+48]=X[i];
        WI[i+48]=W[i];
    }

    for (int i=0;i<96;i++){
        XX[i]=0.5*(XI[i]+1.);
    }

    XX[96]=1.0;
    double EXPON=2.0;//! G.W. 04/07/2007
    //EXPON=1.50;
    //     EXPON=3.0;
    double YI=0.0;
    for(int i=0;i<96;i++){
        YI=2.*pow((0.5*(1.+XI[i])),EXPON)-1.;
        WI[i]=WI[i]/(1.+XI[i])*(1.+YI)*EXPON;
        XI[i]=YI;
        XX[i]=0.5*(1.+YI);
    }
    return 0;
}
/*}}}*/

double Lsidis::FUUT_NLO(bool bNLO=false, bool bNeutron=false){//calculate the NLO structure function F_UU,T {{{

    /*Initialization{{{*/
    const double CF=4.0/3.0;
    double D1_tmp[9],f1_tmp[9], pdf_save[9][96];
    double fF1x_sidis=0.;
    double fFLx_sidis=0.;

    // nucleon mass and inelasticity, QUESTION: y_inel=y, right? ZYE 05/03/2022
    //double mnuc=(ZDIS*938.27208816/1000.+(ADIS-ZDIS)*939.56542052/
    //        1000.)/(ADIS*1.);
    //double y_inel=Q2/x/(2.*mnuc*ebeam);

    if(bNeutron){//in GetPFs(), PDFs were divided by x already, get it back here
        f1_tmp[0] = x*f1[1] * eu2 ;//u
        f1_tmp[1] = x*f1[0] * ed2 ;//d
        f1_tmp[2] = x*f1[2] * ed2 ;//s
        f1_tmp[3] = x*f1[4] * eu2 ;//ubar
        f1_tmp[4] = x*f1[3] * ed2 ;//dbar
        f1_tmp[5] = x*f1[5] * ed2 ;//sbar
        f1_tmp[6] = x*f1[6] * eu2 ;//c+cbar
        f1_tmp[7] = x*f1[7] * ed2 ;//b+bbar
        f1_tmp[8] = x*f1[8] ;//g
    }
    else{
        f1_tmp[0] = x*f1[0] * eu2 ;//u
        f1_tmp[1] = x*f1[1] * ed2 ;//d
        f1_tmp[2] = x*f1[2] * ed2 ;//s
        f1_tmp[3] = x*f1[3] * eu2 ;//ubar
        f1_tmp[4] = x*f1[4] * ed2 ;//dbar
        f1_tmp[5] = x*f1[5] * ed2 ;//sbar
        f1_tmp[6] = x*f1[6] * eu2 ;//c+cbar
        f1_tmp[7] = x*f1[7] * ed2 ;//b+bbar
        f1_tmp[8] = x*f1[8] ;//g
    }
    f1_tmp[6] = 0.0;
    f1_tmp[7] = 0.0;

    //call alpha_s from LHAPDF
    double AL= pdfs->alphasQ2(Q2);
    double AL1x=log(1.-x);
    double AL1z=log(1.-z);
    /*}}}*/

    /* start computing the LO cross-sections{{{*/
    double fqqbar=(f1_tmp[0]*D1[0]+f1_tmp[1]*D1[1]+f1_tmp[2]*D1[2]+f1_tmp[3]*D1[3]+f1_tmp[4]*D1[4]
            +f1_tmp[5]*D1[5]+f1_tmp[6]*D1[6]+f1_tmp[7]*D1[7])*z;
    fF1x_sidis=fqqbar;
    /*}}}*/

    // start computing the NLO cross-sections
    if(bNLO){/*{{{*/
        fF1x_sidis=fF1x_sidis*(1. + AL*CF*(-16.+2.*((7.-6.*z-z*z)/4.+ (-3.+2.*z+z*z)*AL1z/2.+AL1z*AL1z)
                +2.*((7.-6.*x-x*x)/4.+(-3.+2.*x+x*x)*AL1x/2.+AL1x*AL1x)+4.*(AL1x*AL1z)));

        for (int i=1;i<NTERMS;i++){/*{{{*/
            double ytemp=0.50*(1.0-x)*XI[i]+0.50*(1.0+x);
            double xy=x/ytemp;

            if(bNeutron){
                pdf_save[0][i]= pdfs->xfxQ2(1,  xy, Q2) * eu2;//u
                pdf_save[1][i]= pdfs->xfxQ2(2,  xy, Q2) * ed2;//d
                pdf_save[2][i]= pdfs->xfxQ2(3,  xy, Q2) * ed2;//s
                pdf_save[3][i]= pdfs->xfxQ2(-1, xy, Q2) * eu2;//ubar
                pdf_save[4][i]= pdfs->xfxQ2(-2, xy, Q2) * ed2;//dbar
                pdf_save[5][i]= pdfs->xfxQ2(-3, xy, Q2) * ed2;//sbar
                pdf_save[6][i]= pdfs->xfxQ2(-4, xy, Q2) * eu2;//c+cbar
                pdf_save[7][i]= pdfs->xfxQ2(-5, xy, Q2) * ed2;//b+bbar
                pdf_save[8][i]= pdfs->xfxQ2(0,  xy, Q2);//g
            }
            else{
                pdf_save[0][i]= pdfs->xfxQ2(2,  xy, Q2) * eu2;//u
                pdf_save[1][i]= pdfs->xfxQ2(1,  xy, Q2) * ed2;//d
                pdf_save[2][i]= pdfs->xfxQ2(3,  xy, Q2) * ed2;//s
                pdf_save[3][i]= pdfs->xfxQ2(-2, xy, Q2) * eu2;//ubar
                pdf_save[4][i]= pdfs->xfxQ2(-1, xy, Q2) * ed2;//dbar
                pdf_save[5][i]= pdfs->xfxQ2(-3, xy, Q2) * ed2;//sbar
                pdf_save[6][i]= pdfs->xfxQ2(-4, xy, Q2) * eu2;//c+cbar
                pdf_save[7][i]= pdfs->xfxQ2(-5, xy, Q2) * ed2;//b+bbar
                pdf_save[8][i]= pdfs->xfxQ2(0,  xy, Q2);//g
            }
            pdf_save[6][i] = 0.0;
            pdf_save[7][i] = 0.0;

        }//for(int i=0;i<NTERNS;i++)}}}


        for(int i=1;i<NTERMS;i++){                 /*{{{*/
            double yz=0.50*(1.0-z)*XI[i]+0.50*(1.0+z);
            double zy=z/yz;
            double AL1=log(1.0-yz);

            D1_tmp[0] = ffs->xfxQ2(2,  zy, Q2);//u
            D1_tmp[1] = ffs->xfxQ2(1,  zy, Q2);//d
            D1_tmp[2] = ffs->xfxQ2(3,  zy, Q2);//s
            D1_tmp[3] = ffs->xfxQ2(-2, zy, Q2);//ubar
            D1_tmp[4] = ffs->xfxQ2(-1, zy, Q2);//dbar
            D1_tmp[5] = ffs->xfxQ2(-3, zy, Q2);//sbar
            D1_tmp[6] = ffs->xfxQ2(-4, zy, Q2);//c+cbar
            D1_tmp[7] = ffs->xfxQ2(-5, zy, Q2);//b+bbar
            D1_tmp[8] = ffs->xfxQ2(0,  zy, Q2);//g

            double fqqbarx=f1_tmp[0]*D1_tmp[0]+f1_tmp[1]*D1_tmp[1]+f1_tmp[2]*D1_tmp[2]+f1_tmp[3]*D1_tmp[3]
                          +f1_tmp[4]*D1_tmp[4]+f1_tmp[5]*D1_tmp[5]+f1_tmp[6]*D1_tmp[6]+f1_tmp[7]*D1_tmp[7];
            double fgqx=(f1_tmp[0]+f1_tmp[1]+f1_tmp[2]+f1_tmp[3]+f1_tmp[4]+f1_tmp[5]+f1_tmp[6]+f1_tmp[7])*D1_tmp[8];

            double C12=CF*2.*((1.+yz*yz)*log(yz)/(1.-yz)+(1.-yz));
            double C13=CF*2.*AL1/(1.-yz)*(1.+yz*yz);
            double C1c=-2.*CF*AL1x*(1.+yz);
            double C1d=CF*4.*AL1x/(1.-yz);
            double C12gq=CF*2.*(yz+(1.+(1.-yz)*(1.-yz))/yz*(log(yz)+AL1));
            double C13gq=CF*2.*AL1x*(1.+(1.-yz)*(1.-yz))/yz;

            fF1x_sidis=fF1x_sidis+.50*(1.0-z)*WI[i]*AL*
                (C12*fqqbarx+C13*(fqqbarx-fqqbar)+C1c*fqqbarx+C1d*(fqqbarx-fqqbar)+(C12gq+C13gq)*fgqx);

            for(int j=1;j<NTERMS;j++){/*{{{*/

                double yx=0.50*(1.0-x)*XI[j]+0.50*(1.0+x);

                double fqqbarx=pdf_save[0][j]*D1_tmp[0]+pdf_save[1][j]*D1_tmp[1]+pdf_save[2][j]*D1_tmp[2]
                              +pdf_save[3][j]*D1_tmp[3]+pdf_save[4][j]*D1_tmp[4]+pdf_save[5][j]*D1_tmp[5]
                              +pdf_save[6][j]*D1_tmp[6]+pdf_save[7][j]*D1_tmp[7];
                double fgqx=(pdf_save[0][j]+pdf_save[1][j]+pdf_save[2][j]+pdf_save[3][j]+pdf_save[4][j]
                            +pdf_save[5][j]+pdf_save[6][j]+pdf_save[7][j])*D1_tmp[8];
                double fqgx=pdf_save[8][j]*(4.*(D1_tmp[0]+D1_tmp[1]+2.*D1_tmp[6])+D1_tmp[2]+D1_tmp[3]+D1_tmp[4]+D1_tmp[5]+2.*D1_tmp[7])/9.;

                double CLqq=CF*8.*yx*yz;
                double CLgq=CF*8.*yx*(1.-yz);
                double CLqg=8.*yx*(1.-yx);

                fFLx_sidis=fFLx_sidis+0.5*(1.0-x)*0.5*(1.0-z)*WI[i]*WI[j]*AL*(CLqq*fqqbarx+CLgq*fgqx+CLqg*fqgx);

                double fqqbarxauxx=fqqbarx-(pdf_save[0][j]*D1[0]+pdf_save[1][j]*D1[1]+pdf_save[2][j]*D1[2]
                        +pdf_save[3][j]*D1[3]+pdf_save[4][j]*D1[4]+pdf_save[5][j]*D1[5]
                        +pdf_save[6][j]*D1[6]+pdf_save[7][j]*D1[7])*z;
                double fqqbarxauxz=fqqbarx-(f1_tmp[0]*D1_tmp[0]+f1_tmp[1]*D1_tmp[1]+f1_tmp[2]*D1_tmp[2]
                        +f1_tmp[3]*D1_tmp[3]+f1_tmp[4]*D1_tmp[4]+f1_tmp[5]*D1_tmp[5]
                        +f1_tmp[6]*D1_tmp[6]+f1_tmp[7]*D1_tmp[7]);
                double fgqxauxz=fgqx-(f1_tmp[0]+f1_tmp[1]+f1_tmp[2]+f1_tmp[3]+f1_tmp[4]+f1_tmp[5]+f1_tmp[6]+f1_tmp[7])*D1_tmp[8];
                double fqgxauxx=fqgx-pdf_save[8][j]*(4.*(D1[0]+D1[1]+2.*D1[6])+D1[2]+D1[3]+D1[4]+D1[5]+2.*D1[7])*z/9.;

                double C1a=CF*4.*(1.+yx*yz);
                double C1b=CF*(-2.)*(1.+yx)/(1.-yz);
                double C1c=CF*(-2.)*(1.+yz)/(1.-yx);
                double C1d=CF*4./(1.-yx)/(1.-yz);
                double C12gq=CF*2.*(2.*(1.+yx-yx*yz)-(1.+yx)/yz);
                double C13gq=CF*2.*(1.+pow((1.-yz),2))/yz/(1.-yx);
                double C12qg=(yx*yx+pow((1.-yx),2))*(1./yz-2.);
                double C13qg=(yx*yx+pow((1.-yx),2))/(1-yz);

                fF1x_sidis=fF1x_sidis+0.5*(1.0-x)*0.5*(1.0-z)*WI[i]*WI[j]*AL*
                    (C1a*fqqbarx+C1b*fqqbarxauxx+C1c*fqqbarxauxz+C1d*(-fqqbarx+fqqbarxauxx+fqqbarxauxz+fqqbar)
                     +C12gq*fgqx+C13gq*fgqxauxz+C12qg*fqgx+C13qg*fqgxauxx);
            }//for(int j=0;j<NTERNS;j++)}}}
        }//for(int i=0;i<NTERNS;i++)}}}
    }//if(bNLO)/*}}}*/

    return (fF1x_sidis/z+2.*(1.-y)/(1+pow((1-y),2))*fFLx_sidis/z)/x;
}/*}}}*/

double Lsidis::FUUT(){//calculate the structure function F_UU,T at current x, z, Q2, Pt{{{
    if(!(st_nucleus && st_hadron)){
        perror("Missing initialization of nucleus and hadron!");
        return -1;
    }
    double FFp = 0.0, FFn = 0.0;
    if (physics_control){
        double PhT2 = z * z * TMDpars[0] + TMDpars[1];
        if (Np > 0 && Wp > MXminp){
            FFp = Np * x * (eu2 * (f1[0] * D1[0] + f1[3] * D1[3]) + ed2 * (f1[1] * D1[1] + f1[2] * D1[2] + f1[4] * D1[4] + f1[5] * D1[5])) * exp(-Pt * Pt / PhT2) / (PI * PhT2);
            double FF_LO = Np * x * FUUT_NLO(false, false);
            double FF_NLO = Np * x * FUUT_NLO(true, false);
            cout<<Form("+++ FFp = %f,  FF_LO=%f, FF_NLO=%f", FFp, FF_LO, FF_NLO)<<endl;
        }
        if (Nn > 0 && Wp > MXminn){
            FFn = Nn * x * (eu2 * (f1[1] * D1[0] + f1[4] * D1[3]) + ed2 * (f1[0] * D1[1] + f1[2] * D1[2] + f1[3] * D1[4] + f1[5] * D1[5])) * exp(-Pt * Pt / PhT2) / (PI * PhT2);
            double FF_LO = Nn * x * FUUT_NLO(false, true);
            double FF_NLO = Nn * x * FUUT_NLO(true, true);
            cout<<Form("--- FFn = %f,  FF_LO=%f, FF_NLO=%f", FFp, FF_LO, FF_NLO)<<endl;        }
    }
    return FFp + FFn;
}/*}}}*/

double Lsidis::dsigma(const int mode = 0){//calculate the differential cross section at current x, y, z, Pt2, phih, phiS{{{
    double ds = 1e-16;
    if (physics_control){
        if (mode == 0){//No azimuthal modulations
            ds = alpha_em * alpha_em * y / (2.0 * x * Q2 * (1.0 - epsilon)) * (1.0 + gamma * gamma / (2.0 * x)) * FUUT();
        }
        else {
            perror("No such mode option!");
            return -1;
        }
    }
    return ds; //in GeV^-4
}/*}}}*/

int Lsidis::SetRange(const double * mins, const double * maxs){//Set uniform sampling range{{{
    volume = 1.0;
    for (int i = 0; i < 6; i++){
        Xmin[i] = mins[i];
        Xmax[i] = maxs[i];
        volume = volume * (maxs[i] - mins[i]);
    }
    return 0;
}/*}}}*/

double Lsidis::GenerateEvent(const int mode = 0, const int method = 0){//Generate an event and return the cross section weight{{{
    double var[6];
    double weight = 0;
    for (int i = 0; i < 6; i++){
        var[i] = gRandom->Uniform(Xmin[i], Xmax[i]);
    }
    if (method == 0){//generate in x, y, z, Pt, phih, phiS
        jacobian = 2.0 * var[3];
        SetVariables(var[0], var[1], var[2], var[3], var[4], var[5]);
        CalculateFinalState();
        weight = dsigma(mode);
    }
    else if (method == 1){//generate in x, Q2, z, Pt, phih, phiS
        double y0 = var[1] / (2.0 * var[0] * (PP * Pl));
        jacobian = 2.0 * var[3] * y0 / var[1];
        SetVariables(var[0], y0, var[2], var[3], var[4], var[5]);
        CalculateFinalState();
        weight = dsigma(mode);
    }
    else {
        perror("No such method choice!");
        return -1;
    }
    return weight * volume * jacobian;
}/*}}}*/

double Lsidis::GetEventWeight(const int mode = 0, const int method = 0){//Generate an event and return the cross section weight{{{
    double weight = 0;
    if (method == 0){//generate in x, y, z, Pt, phih, phiS
        jacobian = 2.0 * Pt;
        weight = dsigma(mode);
    }
    else if (method == 1){//generate in x, Q2, z, Pt, phih, phiS
        double y0 = Q2 / (2.0 * x * (PP * Pl));
        jacobian = 2.0 * Pt * y0 / Q2;
        weight = dsigma(mode);
    }
    else {
        perror("No such method choice!");
        return -1;
    }
    return weight * volume * jacobian;
}/*}}}*/

int Lsidis::CalculateSigmaTotal(const int mode = 0, const int method = 0, const int precision = 3){//Calculate the total cross section within the simulation range defined by [Xmin, Xmax] in unit of GeV^-2{{{
    double Nsim = pow(10, precision * 2);
    double sum = 0.0;
    for (Long64_t i = 0; i < Nsim; i++){
        sum += GenerateEvent(mode, method);
    }
    sigmatotal = sum / Nsim;
    return 0;
}/*}}}*/

int Lsidis::GibbsStarter(const int mode = 0, const int method = 0){//prepare for Gibbs sampler{{{
    CalculateSigmaTotal(mode, method);//using default precision
    double var[6];
    double weight = 0;
    double judge = 0;
    const int itermax = 10;
    int iter = 0;
    while (judge < 1.0e-9 && iter < itermax){
        iter++;
        for (int ib = 0; ib < 1000; ib++){//search initial Xseed
            for (int i = 0; i < 6; i++){
                var[i] = gRandom->Uniform(Xmin[i], Xmax[i]);
            }
            if (method == 0){//generate in x, y, z, Pt, phih, phiS
                jacobian = 2.0 * var[3];
                SetVariables(var[0], var[1], var[2], var[3], var[4], var[5]);
                CalculateFinalState();
                weight = dsigma(mode) * jacobian;
            }
            else if (method == 1){//generate in x, Q2, z, Pt, phih, phiS
                double y0 = var[1] / (2.0 * var[0] * (PP * Pl));
                jacobian = 2.0 * var[3] * y0 / var[1];
                SetVariables(var[0], y0, var[2], var[3], var[4], var[5]);
                CalculateFinalState();
                weight = dsigma(mode) * jacobian;
            }
            else {
                perror("No such method choice!");
                return -1;
            }
            if (weight > judge){
                judge = weight;
                for (int ii = 0; ii < 6; ii++)
                    Xseed[ii] = var[ii];
            }
        }
    }
    if (judge == 0.0){
        perror("Gibbs starter failed!");
        return -1;
    }
    if (iter >= itermax){
        printf("Warning: max iteration reached by Gibbs starter!  %.3E\n", judge);
    }
    st_gibbs = true;
    for (int ih = 0; ih < 6; ih++){
        Xhisto[ih] = new TH1D(Form("X%d", ih), "", 1000, Xmin[ih], Xmax[ih]);
    }
    return 0;
}/*}}}*/

double Lsidis::GibbsSampler(const int mode = 0, const int method = 0){//Generate an event with Gibbs sampler{{{
    if (!st_gibbs){
        GibbsStarter(mode, method);
    }
    if (!st_gibbs) return -1;
    int dim = 6;
    if (mode == 0){
        dim = 4;
        Xseed[4] = gRandom->Uniform(-M_PI, M_PI);
        Xseed[5] = gRandom->Uniform(-M_PI, M_PI);
        for (int ih = 0; ih < dim; ih++){
            for (int ib = 1; ib <= Xhisto[ih]->GetNbinsX(); ib++){//set histogram
                Xseed[ih] = Xhisto[ih]->GetBinCenter(ib);
                if (method == 0){//generate in x, y, z, Pt, phih, phiS
                    jacobian = 2.0 * Xseed[3];
                    SetVariables(Xseed[0], Xseed[1], Xseed[2], Xseed[3], Xseed[4], Xseed[5]);
                }
                else if (method == 1){//generate in x, Q2, z, Pt, phih, phiS
                    double y0 = Xseed[1] / (2.0 * Xseed[0] * (PP * Pl));
                    jacobian = 2.0 * Xseed[3] * y0 / Xseed[1];
                    SetVariables(Xseed[0], y0, Xseed[2], Xseed[3], Xseed[4], Xseed[5]);
                }
                CalculateFinalState();
                Xhisto[ih]->SetBinContent(ib, dsigma(mode) * jacobian);
            }
            do {//Sample and double check
                if (Xhisto[ih]->Integral(1,-1) > 0){
                    Xseed[ih] = Xhisto[ih]->GetRandom();
                }
                else {
                    Xseed[ih] = gRandom->Uniform(Xmin[ih], Xmax[ih]);
                    //printf("Warning: GibbsSampler to zero point!\n");
                }
                if (method == 0){//generate in x, y, z, Pt, phih, phiS
                    jacobian = 2.0 * Xseed[3];
                    SetVariables(Xseed[0], Xseed[1], Xseed[2], Xseed[3], Xseed[4], Xseed[5]);
                }
                else if (method == 1){//generate in x, Q2, z, Pt, phih, phiS
                    double y0 = Xseed[1] / (2.0 * Xseed[0] * (PP * Pl));
                    jacobian = 2.0 * Xseed[3] * y0 / Xseed[1];
                    SetVariables(Xseed[0], y0, Xseed[2], Xseed[3], Xseed[4], Xseed[5]);
                }
                CalculateFinalState();
            } while(dsigma(mode) < 1.e-20);
        }
    }
    return sigmatotal;
}/*}}}*/

int Lsidis::Test(){//Test code{{{
    std::cout << x << "  " << y << "  " << z << "  " << Pt << "  " << phih << "  " << phiS << std::endl;
    return 0;
}/*}}}*/

#endif

/* update log from v1.1 to 2.0
   implement Gibbs sampler for sidis event generator

   update log from v2.0 to 2.1
   fix pid bug
   */


// Local Variables:
// mode: c++
// End:
