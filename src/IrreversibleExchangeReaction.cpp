//-------------------------------------------------------------------------------------------
//
// IrreversibleExchangeReaction.cpp
//
// Author: Struan Robertson
// Date:   30/Dec/2007
//
// This file contains the implementation of the IrreversibleExchangeReaction class.
//
//-------------------------------------------------------------------------------------------
#include <limits>
#include "IrreversibleExchangeReaction.h"
#include "gStructure.h"

using namespace Constants;
using namespace std;
using namespace mesmer;

namespace mesmer
{
	// Read the Molecular data from input stream.
	bool IrreversibleExchangeReaction::InitializeReaction(PersistPtr ppReac)
	{
		m_ppPersist = ppReac;

		PersistPtr ppReactantList = ppReac->XmlMoveTo("reactantList");
		if (!ppReactantList)
			ppReactantList = ppReac; //Be forgiving; we can get by without a reactantList element

		PersistPtr ppReactant1 = ppReactantList->XmlMoveTo("reactant");      // Read reactant details.
		Molecule* pMol1 = GetMolRef(ppReactant1);
		if (!pMol1){
			cerr << "Cannot find 1st reactant molecule definition for association reaction " << getName() << ".";
			return false;
		}
		PersistPtr ppReactant2 = ppReactant1->XmlMoveTo("reactant");
		Molecule* pMol2 = GetMolRef(ppReactant2);
		if (!pMol2)
		{
			cerr << "Cannot find 2nd reactant molecule definition for association reaction " << getName() << ".";
			return false;
		}

		// if deficientReactantLocation=true, then pMol1 (the first rct
		// in the XML input) is the deficient reactant (m_rct1)

		Molecule* tmp_rct1 = pMol1;
		Molecule* tmp_rct2 = pMol2;

		if (deficientReactantLocation){
			m_rct1 = tmp_rct1;
			m_rct2 = tmp_rct2;
		}
		else {
			m_rct1 = tmp_rct2;
			m_rct2 = tmp_rct1;
		}

		if (!m_rct1){
			cerr << "the deficient reactant in the association reaction is undefined" << endl;
			return false;
		}
		if (!m_rct2){
			cerr << "the excess reactant in the association reaction is undefined" << endl;
			return false;
		}

		PersistPtr ppProductList = ppReac->XmlMoveTo("productList");
		if (!ppProductList)
			ppProductList = ppReac; //Be forgiving; we can get by without a productList element

		PersistPtr ppProduct1 = ppProductList->XmlMoveTo("product");     // Read product details. Save them as type Molecule
		if (ppProduct1) {
			pMol1 = GetMolRef(ppProduct1);
			if (pMol1) {
				m_pdt1 = pMol1;
			}
			else {
				cerr << "Exchange reaction" << getName() << " has no products defined.";
			}

			PersistPtr ppProduct2 = ppProduct1->XmlMoveTo("product");
			if (ppProduct2) {
				pMol2 = GetMolRef(ppProduct2);
				if (pMol2) {
					m_pdt2 = pMol2;
				}
				else {
					cerr << "Exchange reaction " << getName() << " has only one product defined.";
				}
			}
		}

		return ReadRateCoeffParameters(ppReac);       // Read heat of reaction and rate parameters.

	}

	//
	// Get Grain canonical partition function for rotational, vibrational, and electronic contributions.
	//
	double IrreversibleExchangeReaction::rctsRovibronicGrnCanPrtnFn() {
		vector<double> rctGrainDOS;
		vector<double> rctGrainEne;
		calcRctsGrainDensityOfStates(rctGrainDOS, rctGrainEne);

		// Calculate the rovibronic partition function based on the grain DOS
		return canonicalPartitionFunction(rctGrainDOS, rctGrainEne, getEnv().beta);
	}

	//
	// Get Grain canonical partition function for rotational, vibrational, and electronic contributions.
	//
	double IrreversibleExchangeReaction::pdtsRovibronicGrnCanPrtnFn() {
		vector<double> pdtGrainDOS;
		vector<double> pdtGrainEne;
		calcPdtsGrainDensityOfStates(pdtGrainDOS, pdtGrainEne);

		// Calculate the rovibronic partition function based on the grain DOS
		return canonicalPartitionFunction(pdtGrainDOS, pdtGrainEne, getEnv().beta);
	}

	double IrreversibleExchangeReaction::calcEquilibriumConstant() {   // Calculate reaction equilibrium constant.
		// equilibrium constant:
		double Keq(0.0);
		const double beta = getEnv().beta;

		// rovibronic partition function for products/reactants
		double Qrcts = rctsRovibronicGrnCanPrtnFn();
		double Qpdts = pdtsRovibronicGrnCanPrtnFn();

		// rovibronic partition function for reactants/products multiplied by translation contribution
		Qrcts *= translationalContribution(m_rct1->getStruc().getMass(), m_rct2->getStruc().getMass(), beta);
		Qpdts *= translationalContribution(m_pdt1->getStruc().getMass(), m_pdt2->getStruc().getMass(), beta);

		Keq = Qpdts / Qrcts;

		// Heat of reaction: use heat of reaction to calculate the zpe weighing of different wells
		const double HeatOfReaction = getHeatOfReaction();
		Keq *= exp(-beta * HeatOfReaction);

		return Keq;
	}

	// Add exchange reaction terms to collision matrix.
	void IrreversibleExchangeReaction::AddReactionTerms(qdMatrix *CollOptr, molMapType &isomermap, const double rMeanOmega)
	{
		const int jj = (*m_sourceMap)[get_pseudoIsomer()];
		(*CollOptr)[jj][jj] -= qd_real(rMeanOmega) * qd_real(m_MtxGrnKf[0]);
	}

	bool IrreversibleExchangeReaction::calcRctsGrainDensityOfStates(std::vector<double>& grainDOS, std::vector<double>& grainEne)    // Calculate rovibrational reactant DOS
	{
		std::vector<double> rctsCellDOS;
		getRctsCellDensityOfStates(rctsCellDOS);

		const int MaximumCell = getEnv().MaxCell;
		const int cellOffset = get_pseudoIsomer()->getDOS().get_cellOffset();
		std::vector<double> rctsCellEne;
		getCellEnergies(MaximumCell, getEnv().CellSize, rctsCellEne);

		string catName = m_rct1->getName() + " + " + m_rct2->getName();

		if (getFlags().cyclePrintCellDOS){
			ctest << endl << "Cell rovibronic density of states of " << catName << endl << "{" << endl;
			for (int i = 0; i < MaximumCell; ++i){
				formatFloat(ctest, rctsCellEne[i], 6, 15);
				formatFloat(ctest, rctsCellDOS[i], 6, 15);
				ctest << endl;
			}
			ctest << "}" << endl;
			getFlags().cyclePrintCellDOS = false;
		}

		calcGrainAverages(getEnv().MaxGrn, getEnv().cellPerGrain(), cellOffset, rctsCellDOS, rctsCellEne, grainDOS, grainEne);

		if (getFlags().cyclePrintGrainDOS){
			ctest << endl << "Grain rovibronic density of states of " << catName << endl << "{" << endl;
			for (size_t i(0); i < getEnv().MaxGrn; ++i){
				formatFloat(ctest, grainEne[i], 6, 15);
				formatFloat(ctest, grainDOS[i], 6, 15);
				ctest << endl;
			}
			ctest << "}" << endl;
			getFlags().cyclePrintGrainDOS = false;
		}

		return true;
	}

	//
	// Calculate the rovibrational density of states of products.
	//
	bool IrreversibleExchangeReaction::calcPdtsGrainDensityOfStates(std::vector<double>& grainDOS, std::vector<double>& grainEne)
	{
		std::vector<double> pdtsCellDOS;
		getPdtsCellDensityOfStates(pdtsCellDOS);

		const int MaximumCell = getEnv().MaxCell;
		const int cellOffset = m_pdt1->getDOS().get_cellOffset(); // ** temporary statement to get cellOffset from one of the molecules.
		std::vector<double> pdtsCellEne;
		getCellEnergies(MaximumCell, getEnv().CellSize, pdtsCellEne);

		const string catName = m_pdt1->getName() + " + " + m_pdt2->getName();

		if (getFlags().cyclePrintCellDOS){
			ctest << endl << "Cell rovibronic density of states of " << catName << endl << "{" << endl;
			for (int i = 0; i < MaximumCell; ++i){
				formatFloat(ctest, pdtsCellEne[i], 6, 15);
				formatFloat(ctest, pdtsCellDOS[i], 6, 15);
				ctest << endl;
			}
			ctest << "}" << endl;
			getFlags().cyclePrintCellDOS = false;
		}

		calcGrainAverages(getEnv().MaxGrn, getEnv().cellPerGrain(), cellOffset, pdtsCellDOS, pdtsCellEne, grainDOS, grainEne);

		if (getFlags().cyclePrintGrainDOS){
			ctest << endl << "Grain rovibronic density of states of " << catName << endl << "{" << endl;
			for (size_t i(0); i < getEnv().MaxGrn; ++i){
				formatFloat(ctest, grainEne[i], 6, 15);
				formatFloat(ctest, grainDOS[i], 6, 15);
				ctest << endl;
			}
			ctest << "}" << endl;
			getFlags().cyclePrintGrainDOS = false;
		}

		return true;
	}

	//
	// Get reactants cell density of states.
	//
	void IrreversibleExchangeReaction::getRctsCellDensityOfStates(vector<double> &cellDOS) {
		countDimerCellDOS(m_rct1->getDOS(), m_rct2->getDOS(), cellDOS);
	}

	//
	// Get products cell density of states.
	//
	void IrreversibleExchangeReaction::getPdtsCellDensityOfStates(vector<double> &cellDOS) {
		countDimerCellDOS(m_pdt1->getDOS(), m_pdt2->getDOS(), cellDOS);
	}

	// Calculate grained forward and reverse k(E)s from transition state flux
	void IrreversibleExchangeReaction::calcGrainRateCoeffs(){

		vector<double> rctGrainDOS;
		vector<double> rctGrainEne;
		calcRctsGrainDensityOfStates(rctGrainDOS, rctGrainEne);

		calcEffGrnThresholds();
		const int forwardTE = get_EffGrnFwdThreshold();
		calcFluxFirstNonZeroIdx();
		const int fluxStartIdx = get_fluxFirstNonZeroIdx();

		const int MaximumGrain = (getEnv().MaxGrn - fluxStartIdx);
		m_GrainKfmc.clear();
		m_GrainKfmc.resize(MaximumGrain, 0.0);

		for (int i = forwardTE, j = fluxStartIdx; i < MaximumGrain; ++i, ++j){
			if (rctGrainDOS[i] == 0.0){ m_GrainKfmc[i] = 0.0; }
			else{ m_GrainKfmc[i] = m_GrainFlux[j] / rctGrainDOS[i]; }
		}

		if (getFlags().kfEGrainsEnabled){               // printing of the forward k(E)s
			ctest << "\nk_f(e) grains for " << getName() << ":\n{\n";
			for (int i = 0; i < MaximumGrain; ++i){
				ctest << m_GrainKfmc[i] << endl;
			}
			ctest << "}\n";
		}
		if (getFlags().grainTSsosEnabled){
			ctest << "\nN(e) for TS of " << getName() << " (referenced to " << (this->get_pseudoIsomer())->getName() << " energy):\n{\n";
			for (int i = 0; i < MaximumGrain; ++i){
				ctest << m_GrainKfmc[i] * rctGrainDOS[i] / SpeedOfLight_in_cm << endl;
			}
			ctest << "}\n";
		}

		// Always execute this routine for irreversible exchange reactions
		// in order to derive the irreversible collision matrix element.
		HighPresRateCoeffs(NULL);
	}

	// Calculate high pressure rate coefficients at current T.
	void IrreversibleExchangeReaction::HighPresRateCoeffs(vector<double> *pCoeffs) {

		const int MaximumGrain = (getEnv().MaxGrn - get_fluxFirstNonZeroIdx());
		const double beta = getEnv().beta;
		vector<double> rctsGrainDOS;
		vector<double> rctsGrainEne;
		double k_forward(0.0);

		calcRctsGrainDensityOfStates(rctsGrainDOS, rctsGrainEne);
		for (int i(0); i < MaximumGrain; ++i){
			k_forward += m_GrainKfmc[i] * exp(log(rctsGrainDOS[i]) - beta * rctsGrainEne[i]);
		}
		const double prtfn = canonicalPartitionFunction(rctsGrainDOS, rctsGrainEne, beta);
		const double trans = translationalContribution(m_rct1->getStruc().getMass(), m_rct2->getStruc().getMass(), beta);
		k_forward /= prtfn;
		k_forward /= trans;

		// Save high pressure rate coefficient for use in the construction of the collision operator.
		m_MtxGrnKf.clear();
		m_MtxGrnKf.push_back(k_forward*m_ERConc);

		if (pCoeffs) {
			pCoeffs->push_back(k_forward);
			const double Keq = calcEquilibriumConstant();
      if(!IsNan(Keq)) {
				pCoeffs->push_back(k_forward / Keq);
				pCoeffs->push_back(Keq);
      }
		}
		else if (getFlags().testRateConstantEnabled) {
			const double temperature = 1. / (boltzmann_RCpK * beta);
			ctest << endl << "Canonical pseudo first order rate constant of irreversible reaction "
				<< getName() << " = " << k_forward*m_ERConc << " s-1 (" << temperature << " K)" << endl;
			ctest << "Canonical bimolecular rate constant of irreversible reaction "
				<< getName() << " = " << k_forward << " cm^3/mol/s (" << temperature << " K)" << endl;
		}
	}

	const int IrreversibleExchangeReaction::get_rctsGrnZPE(){
		double grnZpe = (m_rct1->getDOS().get_zpe() + m_rct2->getDOS().get_zpe() - getEnv().EMin) / getEnv().GrainSize; //convert to grains
		if (grnZpe < 0.0)
			cerr << "Grain zero point energy has to be a non-negative value.";

		return int(grnZpe);
	}

	void IrreversibleExchangeReaction::calcEffGrnThresholds(void){
		int TS_en = get_fluxGrnZPE();// see the comments in calcEffGrnThresholds under AssociationReaction.cpp  
		int rct_en = get_rctsGrnZPE();

		double thresh = get_ThresholdEnergy();
		if (thresh < 0.0){
			set_EffGrnFwdThreshold(0 / getEnv().GrainSize);
		}
		else{
			set_EffGrnFwdThreshold((TS_en - rct_en) / getEnv().GrainSize);
		}
	}

	void IrreversibleExchangeReaction::calcFluxFirstNonZeroIdx(void) {
		double thresh = get_ThresholdEnergy();
		if (thresh < 0.0){ m_GrnFluxFirstNonZeroIdx = int(-thresh / getEnv().GrainSize); }
		else{ m_GrnFluxFirstNonZeroIdx = 0; }
	}

}//namespace
