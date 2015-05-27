//-------------------------------------------------------------------------------------------
//
// IsomerizationReaction.cpp
//
// Author: Struan Robertson
// Date:   30/Dec/2007
//
// This file contains the implementation of the IsomerizationReaction class.
//
//-------------------------------------------------------------------------------------------
#include <limits>
#include "IsomerizationReaction.h"
#include "gWellProperties.h"

using namespace Constants ;
using namespace std;
using namespace mesmer;

namespace mesmer
{
  //
  // Read the Molecular data from input stream.
  //
  bool IsomerizationReaction::InitializeReaction(PersistPtr ppReac)
  {
    m_ppPersist = ppReac;

    // Read reactant details.
    PersistPtr ppReactantList = ppReac->XmlMoveTo("reactantList");
    if(!ppReactantList)
      ppReactantList=ppReac; //Be forgiving; we can get by without a reactantList element

    PersistPtr ppReactant1  = ppReactantList->XmlMoveTo("reactant");
    Molecule* pMol1 = GetMolRef(ppReactant1);
    if(!pMol1){
      cerr << "Cannot get reactant definition for Isomerization reaction " << getName() << ".";
      return false;
    }

    // Save reactant as Molecule.
    Molecule* pColMol = pMol1;
    if(pColMol){
      m_rct1 = pColMol;
    } else {
      cerr << "Isomer reactant must be a colliding molecule";
      return false;
    }

    //Read product details.
    PersistPtr ppProductList = ppReac->XmlMoveTo("productList");
    if(!ppProductList)
      ppProductList=ppReac; //Be forgiving; we can get by without a productList element
    PersistPtr ppProduct1 = ppProductList->XmlMoveTo("product");
    pMol1 = GetMolRef(ppProduct1);
    if (!pMol1) {
      cerr << "Cannot get product definition for Isomerization reaction " << getName() << ".";
      return false;
    }

    // Save product as Molecule.
    pColMol = pMol1;
    if(pColMol){
      m_pdt1 = pColMol;
    } else {
      cerr << "Isomer product must be a colliding molecule";
      return false;
    }

    // Read heat of reaction and rate parameters.
    return ReadRateCoeffParameters(ppReac);

  }

  // Is reaction equilibrating and therefore contributes
  // to the calculation of equilibrium fractions.
  bool IsomerizationReaction::isEquilibratingReaction(double &Keq, Molecule **rct, Molecule **pdt) {

    Keq = calcEquilibriumConstant() ;

    *rct = m_rct1 ;
    *pdt = m_pdt1 ;

    return true ;
  }

  //
  // Calculate reaction equilibrium constant.
  //
  double IsomerizationReaction::calcEquilibriumConstant() {

    double Keq(0.0) ;

    // Get Canonical partition functions.
    double Qrct1 = m_rct1->getDOS().rovibronicGrnCanPrtnFn() ;
    double Qpdt1 = m_pdt1->getDOS().rovibronicGrnCanPrtnFn() ;

    double beta = getEnv().beta ;

    double HeatOfReaction = getHeatOfReaction();
    Keq = (Qpdt1 / Qrct1)*exp(-beta * HeatOfReaction) ;

    return Keq ;
  }

  //
  // Add isomer reaction terms to reaction matrix.
  //
  void IsomerizationReaction::AddReactionTerms(qdMatrix         *CollOptr,
    molMapType       &isomermap,
    const double    rMeanOmega)
  {
    // Get densities of states for detailed balance.
    vector<double> rctDOS;
    vector<double> pdtDOS;
    m_rct1->getDOS().getGrainDensityOfStates(rctDOS) ;
    m_pdt1->getDOS().getGrainDensityOfStates(pdtDOS) ;

    // Locate isomers in system matrix.
    const int rctLocation = isomermap[m_rct1] ;
    const int pdtLocation = isomermap[m_pdt1] ;

    // Need to know the number of grouped grains in both wells.
    const int rShiftedGrains(m_rct1->getColl().reservoirShift());
    const int pShiftedGrains(m_pdt1->getColl().reservoirShift());

    const int colloptrsize = m_pdt1->getColl().get_colloptrsize() + pShiftedGrains ;

    const int forwardThreshE = get_EffGrnFwdThreshold();
    const int reverseThreshE = get_EffGrnRvsThreshold();
    const int fluxStartIdx   = get_fluxFirstNonZeroIdx();

    for ( int i=fluxStartIdx, j = reverseThreshE, k = forwardThreshE; j < colloptrsize; ++i, ++j, ++k) {
      int ii(rctLocation + k - rShiftedGrains) ;
      int jj(pdtLocation + j - pShiftedGrains) ;
      qd_real Flux(m_GrainFlux[i]), qMeanOmega(rMeanOmega), rDos(rctDOS[k]), pDos(pdtDOS[j]) ;
      (*CollOptr)[ii][ii] -= qMeanOmega * Flux / rDos;               // Forward loss reaction.
      (*CollOptr)[jj][jj] -= qMeanOmega * Flux / pDos ;              // Backward loss reaction from detailed balance.
      (*CollOptr)[ii][jj]  = qMeanOmega * Flux / sqrt(rDos * pDos) ; // Reactive gain.
      (*CollOptr)[jj][ii]  = (*CollOptr)[ii][jj] ;                   // Reactive gain.
    }

  }

  //
  // Add isomer reaction terms to contracted basis reaction matrix.
  //
  void IsomerizationReaction::AddContractedBasisReactionTerms(qdMatrix *CollOptr, molMapType &isomermap)
  {
    // Get densities of states for detailed balance.
    vector<double> rctDOS;
    vector<double> pdtDOS;
    m_rct1->getDOS().getGrainDensityOfStates(rctDOS) ;
    m_pdt1->getDOS().getGrainDensityOfStates(pdtDOS) ;

    const int rctColloptrsize = m_rct1->getColl().get_colloptrsize();
    const int pdtColloptrsize = m_pdt1->getColl().get_colloptrsize();

    const int forwardThreshE = get_EffGrnFwdThreshold();
    const int reverseThreshE = get_EffGrnRvsThreshold();
    const int fluxStartIdx   = get_fluxFirstNonZeroIdx();

    vector<double> fwdMicroRateCoef(rctColloptrsize, 0.0) ;
    vector<double> RvsMicroRateCoef(pdtColloptrsize, 0.0) ;
    vector<double> CrsMicroRateCoef(pdtColloptrsize, 0.0) ;
    for ( int i=fluxStartIdx, j = reverseThreshE, k=0; j < rctColloptrsize; ++i, ++j, ++k) {
      int ll = k + forwardThreshE;
      int mm = k + reverseThreshE;
      fwdMicroRateCoef[ll] = m_GrainFlux[i] / rctDOS[ll] ;                    // Forward loss reaction.
      RvsMicroRateCoef[mm] = m_GrainFlux[i] / pdtDOS[mm] ;                    // Backward loss reaction. 
      CrsMicroRateCoef[mm] = m_GrainFlux[i] / sqrt(rctDOS[ll] * pdtDOS[mm]) ; // Reactive gain from detailed balance.
    }

    // Locate isomers in system matrix.
    const int rctLocation = isomermap[m_rct1] ;
    const int pdtLocation = isomermap[m_pdt1] ;

    // Calculate the elements of the reactant block.

    const int rctBasisSize = static_cast<int>(m_rct1->getColl().get_nbasis()) ;

    for (int i=0, ii(rctLocation), egvI(rctColloptrsize-1) ; i < rctBasisSize ; i++, ii++, --egvI) {
      (*CollOptr)[ii][ii] -= m_rct1->getColl().matrixElement(egvI, egvI, fwdMicroRateCoef) ;
      for (int j=i+1, jj(rctLocation + j), egvJ(rctColloptrsize-j-1) ; j < rctBasisSize ; j++, jj++, --egvJ) {
        qd_real tmp = m_rct1->getColl().matrixElement(egvI, egvJ, fwdMicroRateCoef) ;
        (*CollOptr)[ii][jj] -= tmp ;
        (*CollOptr)[jj][ii] -= tmp ;
      }
    }

    // Calculate the elements of the product block.

    const int pdtBasisSize = static_cast<int>(m_pdt1->getColl().get_nbasis());

    for (int i=0, ii(pdtLocation), egvI(pdtColloptrsize-1) ; i < pdtBasisSize ; i++, ii++, --egvI) {
      (*CollOptr)[ii][ii] -= m_pdt1->getColl().matrixElement(egvI, egvI, RvsMicroRateCoef) ;
      for (int j=i+1, jj(pdtLocation + j), egvJ(pdtColloptrsize-j-1)  ; j < pdtBasisSize ; j++, jj++, --egvJ) {
        qd_real tmp = m_pdt1->getColl().matrixElement(egvI, egvJ, RvsMicroRateCoef) ;
        (*CollOptr)[ii][jj] -= tmp ;
        (*CollOptr)[jj][ii] -= tmp ;
      }
    }

    // Calculate the elements of the cross block.

    vector<double> rctBasisVector(rctColloptrsize, 0.0) ;
    vector<double> pdtBasisVector(pdtColloptrsize, 0.0) ;
    for (int i=0, rctEgv(rctColloptrsize-1) ; i < rctBasisSize ; i++, --rctEgv) {
      int ii(rctLocation + i) ;
      m_rct1->getColl().eigenVector(rctEgv, rctBasisVector) ;
      for (int j=0, pdtEgv(pdtColloptrsize-1)  ; j < pdtBasisSize ; j++, --pdtEgv) {
        int jj(pdtLocation + j) ;
        qd_real tmp(0.0) ;

        if (i==0 && j==0 ) {

          // Special case for equilibrium eigenvectors which obey a detailed balance relation.
          // SHR, 8/Mar/2009: are there other relations like this I wonder.

          qd_real elmti = (*CollOptr)[ii][ii] ;
          qd_real elmtj = (*CollOptr)[jj][jj] ;
          tmp = sqrt(elmti*elmtj) ;

        } else {

          // General case.

          m_pdt1->getColl().eigenVector(pdtEgv, pdtBasisVector) ;
          double sum = 0.0 ;
          for (int k(rctColloptrsize-reverseThreshE), m(rctColloptrsize-1), n(pdtColloptrsize-1); k >=0 ; --m, --n, --k) {
            // tmp += rctBasisVector[m]*rctBasisVector[m]*fwdMicroRateCoef[m];
            // tmp += pdtBasisVector[n]*pdtBasisVector[n]*RvsMicroRateCoef[n];
            // tmp += rctBasisVector[m]*pdtBasisVector[n]*RvsMicroRateCoef[n]*sqrt(pdtDOS[n]/rctDOS[m]);
            sum += rctBasisVector[m]*pdtBasisVector[n]*CrsMicroRateCoef[n];
          }

          tmp = qd_real(sum) ;
        }
        (*CollOptr)[ii][jj] += tmp ;
        (*CollOptr)[jj][ii] += tmp ;
      }
    }

  }

  //
  // Calculate grained forward and reverse k(E)s from transition state flux
  //
  void IsomerizationReaction::calcGrainRateCoeffs(){
    vector<double> rctGrainDOS;
    vector<double> pdtGrainDOS;
    m_rct1->getDOS().getGrainDensityOfStates(rctGrainDOS) ;
    m_pdt1->getDOS().getGrainDensityOfStates(pdtGrainDOS) ;

    calcEffGrnThresholds();
    const int forwardTE = get_EffGrnFwdThreshold();
    int reverseTE = get_EffGrnRvsThreshold();
    calcFluxFirstNonZeroIdx();
    const int fluxStartIdx = get_fluxFirstNonZeroIdx();

    const int MaximumGrain = (getEnv().MaxGrn-fluxStartIdx);
    m_GrainKfmc.clear();
    m_GrainKfmc.resize(MaximumGrain , 0.0);
    m_GrainKbmc.clear();
    m_GrainKbmc.resize(MaximumGrain , 0.0);

    for (int i = reverseTE, j = fluxStartIdx; i < MaximumGrain; ++i, ++j){
      m_GrainKbmc[i] = m_GrainFlux[j] / pdtGrainDOS[i];
    }
    for (int i = forwardTE, j = fluxStartIdx; i < MaximumGrain; ++i, ++j){
      m_GrainKfmc[i] = m_GrainFlux[j] / rctGrainDOS[i];
    }

    // the code that follows is for printing of the f & r k(E)s
    if (getFlags().kfEGrainsEnabled){
      ctest << "\nk_f(e) grains for " << getName() << ":\n{\n";
      for (int i = 0; i < MaximumGrain; ++i){
        ctest << m_GrainKfmc[i] << endl;
      }
      ctest << "}\n";
    }
    if (getFlags().kbEGrainsEnabled){
      ctest << "\nk_b(e) grains for " << getName() << ":\n{\n";
      for (int i = 0; i < MaximumGrain; ++i){
        ctest << m_GrainKbmc[i] << endl;
      }
      ctest << "}\n";
    }
    if (getFlags().grainTSsosEnabled){
      ctest << "\nN(e) for TS of " << getName() << " (referenced to " << (this->get_reactant())->getName() << " energy):\n{\n";
      for (int i = 0; i < MaximumGrain; ++i){
        ctest << m_GrainKfmc[i]*rctGrainDOS[i]/SpeedOfLight_in_cm << endl;
      }
      ctest << "}\n";
    }
    if (getFlags().testRateConstantEnabled)
      HighPresRateCoeffs(NULL);
  }

  // Calculate high pressure rate coefficients at current T.
  void IsomerizationReaction::HighPresRateCoeffs(vector<double> *pCoeffs) {

    double k_forward(0.0), k_backward(0.0);
    vector<double> rctGrainDOS, rctGrainEne, pdtGrainDOS, pdtGrainEne ;
    m_rct1->getDOS().getGrainDensityOfStates(rctGrainDOS);
    m_pdt1->getDOS().getGrainDensityOfStates(pdtGrainDOS);
    m_rct1->getDOS().getGrainEnergies(rctGrainEne);
    m_pdt1->getDOS().getGrainEnergies(pdtGrainEne);
    const int MaximumGrain = (getEnv().MaxGrn-get_fluxFirstNonZeroIdx());
    const double beta = getEnv().beta;
    const double temperature = 1. / (boltzmann_RCpK * beta);

    const int forwardTE = get_EffGrnFwdThreshold();
    const int reverseTE = get_EffGrnRvsThreshold();
    calcFluxFirstNonZeroIdx();
    const int fluxStartIdx = get_fluxFirstNonZeroIdx();
    for ( int i=fluxStartIdx, j = reverseTE, k = forwardTE; max(j,k) < MaximumGrain ; ++i, ++j, ++k) {
      k_forward  += m_GrainFlux[i] * exp(-beta * rctGrainEne[k]);
      k_backward += m_GrainFlux[i] * exp(-beta * pdtGrainEne[j]);
    }

    const double rctprtfn = canonicalPartitionFunction(rctGrainDOS, rctGrainEne, beta);
    const double pdtprtfn = canonicalPartitionFunction(pdtGrainDOS, pdtGrainEne, beta);
    k_forward /= rctprtfn;
    k_backward /= pdtprtfn;

    if (pCoeffs) {
      const double Keq = calcEquilibriumConstant() ;
      pCoeffs->push_back(k_forward) ;
      pCoeffs->push_back(k_backward) ;
      pCoeffs->push_back(Keq) ;
    } else {
      ctest << endl << "Canonical first order forward rate constant of isomerization reaction " 
				<< getName() << " = " << k_forward << " s-1 (" << temperature << " K)" << endl;
      ctest << "Canonical first order backward rate constant of isomerization reaction " 
				<< getName() << " = " << k_backward << " s-1 (" << temperature << " K)" << endl;
    }
  }

  void IsomerizationReaction::calcEffGrnThresholds(void){  // see the comments in
    double thresh = get_ThresholdEnergy();    // calcEffGrnThresholds under AssociationReaction.cpp
    double RxnHeat = getHeatOfReaction();
    if (thresh < RxnHeat && m_pMicroRateCalculator->getID() == "MesmerILT"){
      cerr << "E_infinity should be equal to or greater than the heat of reaction in ILT.";
      exit(1);
    }
    int TS_en = this->get_fluxGrnZPE();
    int pdt_en = m_pdt1->getColl().get_grnZPE();
    int rct_en = m_rct1->getColl().get_grnZPE();
    int GrainedRxnHeat = pdt_en - rct_en;
    if(thresh<0.0){
      set_EffGrnFwdThreshold(0);
      set_EffGrnRvsThreshold(-GrainedRxnHeat);
    }
    else if(thresh>0.0 && thresh<RxnHeat){
      set_EffGrnFwdThreshold(GrainedRxnHeat);
      set_EffGrnRvsThreshold(0);
    }
    else{
      set_EffGrnFwdThreshold(TS_en-rct_en);
      set_EffGrnRvsThreshold(TS_en-pdt_en);
    }
  }

}//namespace
