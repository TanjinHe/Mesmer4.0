//-------------------------------------------------------------------------------------------
//
// IrreversibleUnimolecularReaction.cpp
//
// Author: Struan Robertson
// Date:   30/Dec/2007
//
// This file contains the implementation of the IrreversibleUnimolecularReaction class.
//
//-------------------------------------------------------------------------------------------
#include <limits>
#include "IrreversibleUnimolecularReaction.h"
#include "gWellProperties.h"
#include "gStructure.h"

using namespace Constants ;
using namespace std;
using namespace mesmer;

namespace mesmer
{
  //
  // Read the Molecular data from input stream.
  //
  bool IrreversibleUnimolecularReaction::InitializeReaction(PersistPtr ppReac)
  {
    m_ppPersist = ppReac;

    // Read reactant details.
    PersistPtr ppReactantList = ppReac->XmlMoveTo("reactantList");
    if(!ppReactantList)
      ppReactantList=ppReac; //Be forgiving; we can get by without a reactantList element

    PersistPtr ppReactant1  = ppReactantList->XmlMoveTo("reactant");
    Molecule* pMol1 = GetMolRef(ppReactant1);
    if(!pMol1){
      cerr << "Cannot find reactant molecule definition for irreversible reaction " << getName() << ".";
      return false;
    }

    // Save reactant as Molecule.

    Molecule* pColMol = pMol1;
    if(pColMol){
      m_rct1 = pColMol;
    } else {
      cerr << "Reactant of a irreversible reaction must be a colliding molecule";
      return false;
    }

    // Read product details. The detail of products may be absent or, may be needed
    // to calculate the microcanonical rates. If there are products, save them as
    // type Molecule.
    PersistPtr ppProductList = ppReac->XmlMoveTo("productList");
    if(!ppProductList)
      ppProductList=ppReac; //Be forgiving; we can get by without a productList element

    PersistPtr ppProduct1 = ppProductList->XmlMoveTo("product");
    if (ppProduct1) {
      pMol1 = GetMolRef(ppProduct1);
      if (pMol1){
        m_pdt1 = pMol1;
      }
      else {
        cerr << "Irreversible reaction" << getName() << " has no product defined." << endl;
        return false;
      }

      Molecule* pMol2 = NULL ;
      PersistPtr ppProduct2  = ppProduct1->XmlMoveTo("product");
      if (ppProduct2) {
        pMol2 = GetMolRef(ppProduct2);
        if (pMol2){
          m_pdt2 = pMol2;
          Molecule* pMol3 = NULL ;
          PersistPtr ppProduct3  = ppProduct2->XmlMoveTo("product");
          if (ppProduct3) {
            pMol3 = GetMolRef(ppProduct3);
            if (pMol3)
              m_pdt3 = pMol3;
          }
        }
        else {
          cinfo << "Irreversible reaction " << getName() << " has only one product defined." << endl;
        }
      }
    }

    // Read heat of reaction and rate parameters.
    return ReadRateCoeffParameters(ppReac) ;

  }

  //
  // Calculate reaction equilibrium constant.
  //
  double IrreversibleUnimolecularReaction::calcEquilibriumConstant() {

    // equilibrium constant:
    double Keq(0.0) ;
    const double beta = getEnv().beta ;

    // partition function for each products
    double Qpdts = m_pdt1->getDOS().rovibronicGrnCanPrtnFn() ;

    // rovibronic partition function for products multiplied by translation contribution
    if (m_pdt2){
      Qpdts *= m_pdt2->getDOS().rovibronicGrnCanPrtnFn();
      Qpdts *= translationalContribution(m_pdt1->getStruc().getMass(), m_pdt2->getStruc().getMass(), beta);
    }

    // rovibronic partition function for reactant
    const double Qrct1 = m_rct1->getDOS().rovibronicGrnCanPrtnFn() ;

    Keq = Qpdts / Qrct1;

    // Heat of reaction: use heat of reaction to calculate the zpe weighing of different wells
    const double HeatOfReaction = getHeatOfReaction() ;
    Keq *= exp(-beta * HeatOfReaction) ;

    return Keq ;
  }

  //
  // Add dissociation reaction terms to collision matrix.
  //
  void IrreversibleUnimolecularReaction::AddReactionTerms(qdMatrix *CollOptr, molMapType &isomermap, const double rMeanOmega) {
    // Get densities of states for detailed balance.
    vector<double> rctDOS;
    m_rct1->getDOS().getGrainDensityOfStates(rctDOS) ;

    // Locate reactant in system matrix.
    const int rctLocation = isomermap[m_rct1] ;
    const int rShiftedGrains(m_rct1->getColl().reservoirShift());

    const int colloptrsize   = m_rct1->getColl().get_colloptrsize();
    const int forwardThreshE = get_EffGrnFwdThreshold();
    const int fluxStartIdx   = get_fluxFirstNonZeroIdx();

    m_MtxGrnKf.clear();
    m_MtxGrnKf.resize(colloptrsize , 0.0);

    for ( int i=fluxStartIdx, j = forwardThreshE; j < colloptrsize + rShiftedGrains; ++i, ++j) {
      int ii(rctLocation + j - rShiftedGrains) ;
      qd_real rtcf = qd_real(m_GrainFlux[i]) / qd_real(rctDOS[j]) ;
      (*CollOptr)[ii][ii] -= qd_real(rMeanOmega) * rtcf ;  // Forward loss reaction.
      m_MtxGrnKf[j - rShiftedGrains] = to_double(rtcf) ;
    }
  }

  //
  // Add contracted basis set reaction terms to the reaction matrix.
  //
  void IrreversibleUnimolecularReaction::AddContractedBasisReactionTerms(qdMatrix *CollOptr, molMapType &isomermap) {

    // Get densities of states for detailed balance.
    vector<double> rctDOS;
    m_rct1->getDOS().getGrainDensityOfStates(rctDOS) ;

    const int rctColloptrsize = m_rct1->getColl().get_colloptrsize();
    const int forwardThreshE  = get_EffGrnFwdThreshold();
    const int fluxStartIdx    = get_fluxFirstNonZeroIdx();

    vector<double> fwdMicroRateCoef(rctColloptrsize, 0.0) ;
    for ( int i=fluxStartIdx, j = forwardThreshE ; j < rctColloptrsize; ++i, ++j) {
      fwdMicroRateCoef[j] = m_GrainFlux[i] / rctDOS[j] ;    // Forward loss reaction.
    }

    // Locate isomers in system matrix.
    const int rctLocation = isomermap[m_rct1] ;

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
  }

  //
  // Calculate grained forward k(E)s from transition state flux
  //
  void IrreversibleUnimolecularReaction::calcGrainRateCoeffs(){
    vector<double> rctGrainDOS;
    m_rct1->getDOS().getGrainDensityOfStates(rctGrainDOS) ;

    calcEffGrnThresholds();
    const int forwardTE = get_EffGrnFwdThreshold();
    calcFluxFirstNonZeroIdx();
    const int fluxStartIdx = get_fluxFirstNonZeroIdx();

    const int MaximumGrain = (getEnv().MaxGrn-fluxStartIdx);
    m_GrainKfmc.clear();
    m_GrainKfmc.resize(MaximumGrain , 0.0);

    // calculate forward k(E)s from flux
    for (int i = forwardTE, j = fluxStartIdx; i < MaximumGrain; ++i, ++j){
      m_GrainKfmc[i] = m_GrainFlux[j] / rctGrainDOS[i];
    }

    // the code that follows is for printing the forward k(E)s
    if (getFlags().kfEGrainsEnabled){
      ctest << "\nk_f(e) grains for " << getName() << ":\n{\n";
      for (int i = 0; i < MaximumGrain; ++i){
        ctest << m_GrainKfmc[i] << endl;
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

  void IrreversibleUnimolecularReaction::calcFluxFirstNonZeroIdx(void) {
    double thresh = get_ThresholdEnergy();
    if(thresh<0.0){m_GrnFluxFirstNonZeroIdx = int(-thresh/getEnv().GrainSize);}
    else{m_GrnFluxFirstNonZeroIdx = 0;}
  }

  // Calculate high pressure rate coefficients at current T.
  void IrreversibleUnimolecularReaction::HighPresRateCoeffs(vector<double> *pCoeffs) {

    vector<double> rctGrainDOS, rctGrainEne;
    m_rct1->getDOS().getGrainDensityOfStates(rctGrainDOS);
    m_rct1->getDOS().getGrainEnergies(rctGrainEne);
    const size_t MaximumGrain = (getEnv().MaxGrn-get_fluxFirstNonZeroIdx());
    const double beta = getEnv().beta;

    double kf(0.0);
    for(size_t i(0); i < MaximumGrain; ++i)
      kf += m_GrainKfmc[i] * exp( log(rctGrainDOS[i]) - beta * rctGrainEne[i]);

    const double rctprtfn = canonicalPartitionFunction(rctGrainDOS, rctGrainEne, beta);
    kf /= rctprtfn;

    if (pCoeffs) {
      pCoeffs->push_back(kf) ;
      const double Keq = calcEquilibriumConstant() ;
      if(!IsNan(Keq)) {
        pCoeffs->push_back(kf/Keq) ;
        pCoeffs->push_back(Keq) ;
      }
    } else {
      const double temperature = 1./(boltzmann_RCpK * beta);
      ctest << endl << "Canonical pseudo first order forward rate constant of irreversible reaction "
				<< getName() << " = " << kf << " s-1 (" << temperature << " K)" << endl;
    }
  }

  void IrreversibleUnimolecularReaction::calcEffGrnThresholds(void){       
    int TS_en = get_fluxGrnZPE();
    int rct_en = m_rct1->getColl().get_grnZPE();

    // Handle the case of the reverse threshold energy being negative only if the 
    // reaction can use product properties. (Probably set because ofreverse ILT).

    double RxnHeat(0.0);
    int pdtsGrnZPE = 0;
    double threshold   = get_ThresholdEnergy();// see the comments in
    if(UsesProductProperties()) {
      RxnHeat          = getHeatOfReaction();  // calcEffGrnThresholds under AssociationReaction.cpp
      pdtsGrnZPE       = get_pdtsGrnZPE();
    }
    int rctGrnZPE      = m_rct1->getColl().get_grnZPE();
    int GrainedRxnHeat = pdtsGrnZPE - rctGrnZPE;

    if(threshold<0.0){
      set_EffGrnFwdThreshold(0);
    }
    else if(UsesProductProperties() && threshold>0.0 && threshold<RxnHeat){// if the reverse threshold energy is negative
      set_EffGrnFwdThreshold( GrainedRxnHeat);  // forward grained flux threshold energy = heat of reaction
    }
    else{
      set_EffGrnFwdThreshold(TS_en-rct_en);
    }
  }

  const int IrreversibleUnimolecularReaction::get_pdtsGrnZPE(){
    double zpe = m_pdt1->getDOS().get_zpe() - getEnv().EMin;
    if (m_pdt2) zpe += m_pdt2->getDOS().get_zpe();
    double grnZpe = zpe / getEnv().GrainSize ; //convert to grain
    if (grnZpe < 0.0)
      cinfo << "Grain zero point energy is negative in " << getName() << ".";

    return int(grnZpe);
  }

}//namespace
