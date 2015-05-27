//-------------------------------------------------------------------------------------------
//
// System.cpp
//
// Author: Struan Robertson
// Date:   11/Feb/2003
//
// This file contains the implementation of the System class.
//
//-------------------------------------------------------------------------------------------
#include "System.h"
#include <fstream>
#include "ParseForPlugin.h"

using namespace std ;
using namespace Constants ;

namespace mesmer
{
  //Global varaiable and function declared in System.h
  std::string libfile;

  PersistPtr GetFromLibrary(const std::string molName, PersistPtr ppMolList)
  {
    //Search the library of molecules, copy to the main XML file and return a pointer to the copy
    // (or library version if ot copied).
    PersistPtr ppNewMol;
    if(molName.empty())
      return ppNewMol;
    static PersistPtr ppLib; //initiallized by system to NULL
    if(!ppLib)
    {
      ppLib = XMLPersist::XmlLoad(libfile,"");
      if(libfile.empty() || !ppLib)
      {
        cwarn << "Could not find Library file to search it for missing molecule(s)."<<endl;
        return false;
      }
    }
    PersistPtr ppLibMolList   = ppLib->XmlMoveTo("moleculeList");
    if(!ppLibMolList)
      ppLibMolList = ppLib; //Can do without <moleculeList>

    PersistPtr ppMol = ppLibMolList->XmlMoveTo("molecule");
    string tmolName(molName);
    const char* libId = NULL;
    while(ppMol)
    {
      if(tmolName==ppMol->XmlReadValue("id", false))
      {
        //ignore library molecules with attribute active="false"
        const char* active = ppMol->XmlReadValue("active", optional);
        if(!active || strcmp(active, "false"))
        {
          //Check whether this match is an alias, e.g.
          // <molecule id="aliasName" ref="libraryName"/> 
          const char* txt = ppMol->XmlReadValue("ref",optional);
          if(txt)
          {
            libId = txt;
            tmolName = libId; //continue looking for real molecule
          }
          else //not an alias
          {            
            if(ppMolList) //no copy if not specified
            {
              //Delete a molecule of same name in datafile, if present
              PersistPtr ppOldMol = ppMolList;
              while(ppOldMol = ppOldMol->XmlMoveTo("molecule"))
              {
                if(molName == ppOldMol->XmlReadValue("id", false))
                  break;
              }

              //Copy a matching molecule from library to the main XML file
              //Replace old version if present
              ppMol = ppMolList->XmlCopy(ppMol, ppOldMol);

              cinfo << molName << " copied from " << libfile;
              //Write its provenance, under metadatList if present
              PersistPtr pp = ppMol->XmlMoveTo("metadataList");
              pp = pp ? pp : ppMol;
              pp->XmlWriteMetadata("copiedFrom", libfile);
              if(libId)//originally an alias 
              {
                ppMol->XmlWriteAttribute("id",molName);
                ppMol->XmlWriteAttribute("libId",libId);
                cinfo << " Original library id = " << libId;
              }
              cinfo << endl;
            }
            return ppMol;
          }
        }
      }
      ppMol = ppMol->XmlMoveTo("molecule");
    }
    cinfo << "Could not find " << molName << " in " << libfile << endl;
    return ppMol; // empty, no suitable molecule found
  }


  System::System(const std::string& libraryfilename) : m_pMoleculeManager(0), 
    m_pReactionManager(0), 
    m_pTitle(NULL),
    m_pDescription(NULL)
  {
    libfile = libraryfilename;
  }

  System::~System() {
    delete m_pReactionManager;
    delete m_pMoleculeManager;
  }

  // Initialize the System object.
  bool System::initialize(void) {
    m_pMoleculeManager = new MoleculeManager() ;
    m_pReactionManager = new ReactionManager(m_pMoleculeManager) ;
    m_pConditionsManager = new ConditionsManager(this) ;
    return m_collisionOperator.initialize(m_pMoleculeManager, m_pReactionManager);
  }

  //
  // Parse an input data file.
  //
  bool System::parse(PersistPtr ppIOPtr)
  {
    initializeConversionMaps();
    m_ppIOPtr = ppIOPtr;

    m_pTitle       = ppIOPtr->XmlReadValue("title", false);
    if(!m_pTitle)
      m_pTitle = ppIOPtr->XmlReadValue("me:title", false);
    m_pDescription = ppIOPtr->XmlReadValue("description", false);
    if(!m_pDescription)
      m_pDescription = ppIOPtr->XmlReadValue("me:description", false);

    //Molecule List
    PersistPtr ppMolList = ppIOPtr->XmlMoveTo("moleculeList");
    m_pMoleculeManager->set_PersistPtr(ppMolList ? ppMolList : ppIOPtr);

    // What are we going to do?

    PersistPtr ppControl = ppIOPtr;//first time at top, then gets next block
    while(ppControl = ppControl->XmlMoveTo("me:control"))
    {
      //Ignore a <control> that has an attribute  active="false","no" or "0"
      const char* txt = ppControl->XmlReadValue("active",optional);
      if(txt && !ppControl->XmlReadBoolean("active"))
        continue;

      m_CalcMethod = ParseForPlugin<CalcMethod>(this, "me:calcMethod", ppControl);
      if(!m_CalcMethod)
        return false;
      m_CalcMethodsForEachControl.push_back(m_CalcMethod); //for each <control> block

      static PersistPtr ppConditions;
      // Parse molecules, reactions, parameters and conditions only
      // in conjunction with the first control block.
      if(m_CalcMethodsForEachControl.size()==1)
      {
        if (m_CalcMethod->DoesOwnParsing()) {
          // Thermodynamic table, UnitTests,etc. do their own file parsing.
          m_FlagsForEachControl.push_back(m_Flags);
          return true ;
        }

        // Model Parameters 
        PersistPtr ppParams;
        // Add this section if it does not exist, to contain defaults
        while(!(ppParams = ppIOPtr->XmlMoveTo("me:modelParameters")))
          ppIOPtr->XmlWriteElement("me:modelParameters");

        // The grain size and grain number are linked via the total maximum energy,
        // so only on of then is independent. Look for grain size first and if that
        // fails look for the Max. number of grains.

        m_Env.GrainSize = ppParams->XmlReadInteger("me:grainSize", optional);
        if (m_Env.GrainSize==0) {      
          m_Env.MaxGrn = ppParams->XmlReadInteger("me:numberOfGrains",optional);
          if (IsNan(m_Env.MaxGrn))
            m_Env.MaxGrn=0;
        }

        m_Env.CellSize = ppParams->XmlReadDouble("me:cellSize", optional);
        if (IsNan(m_Env.CellSize)) {
          m_Env.CellSize = 1.0 ; // Default cell size in cm-1.
        }

        m_Env.MaximumTemperature = ppParams->XmlReadDouble("me:maxTemperature",optional);
        if(IsNan(m_Env.MaximumTemperature))
          m_Env.MaximumTemperature = 0.0;
        m_Env.EAboveHill         = ppParams->XmlReadDouble("me:energyAboveTheTopHill");
        m_Env.useBasisSetMethod  = ppParams->XmlReadBoolean("me:runBasisSetMethodroutines");
        if (m_Env.useBasisSetMethod) {
          PersistPtr ppBasisSet = ppParams->XmlMoveTo("me:runBasisSetMethodroutines");
          if(ppBasisSet) {
            m_Env.nBasisSet = ppBasisSet->XmlReadInteger("me:numberBasisFunctions");
          } else {
            cerr << "Basis set method requested but number of basis functions unspecified.";
            return false;
          }
        }
        PersistPtr pp = ppParams->XmlMoveTo("me:automaticallySetMaxEne");
        if (pp) {
          m_Flags.autoSetMaxEne = true;
          m_Flags.popThreshold = ppParams->XmlReadDouble("me:automaticallySetMaxEne");
        }
        cinfo.flush();

        //Reaction Conditions
        ppConditions = ppIOPtr->XmlMoveTo("me:conditions");
        if(!ppConditions)
        {
          cerr << "No conditions specified";
          return false;
        }

        m_ConditionsForEachControl.push_back(m_pConditionsManager);
        if(!m_pConditionsManager->ParseBathGas(ppConditions))
          return false;;

        //Reaction List
        PersistPtr ppReacList = ppIOPtr->XmlMoveTo("reactionList");
        if(!ppReacList)
        {
          cerr << "No reactions have been specified";
          return false;
        }
        if(!m_pReactionManager->addreactions(ppReacList, m_Env, m_Flags))
          return false;

        if(!m_pConditionsManager->ParseConditions())
          return false;
      }
      else // second and subsequent <control> blocks
      {
        // If there is another <conditions> block parse it and save a pointer to it.
        PersistPtr pp = ppConditions->XmlMoveTo("me:conditions");
        if(pp)
        {
          m_pConditionsManager = new ConditionsManager(this);
          ppConditions = pp;
          if(!m_pConditionsManager->ParseBathGas(pp))
            return false;
          if(!m_pConditionsManager->ParseConditions())
            return false;
          m_ConditionsForEachControl.push_back(m_pConditionsManager);
        }
        else //no additional <control> block found
          m_ConditionsForEachControl.push_back(NULL);
      }

      if(ppControl)
      {
        m_Flags.testDOSEnabled                 = ppControl->XmlReadBoolean("me:testDOS");
        m_Flags.testRateConstantEnabled        = ppControl->XmlReadBoolean("me:testRateConstants");
        m_Flags.microRateEnabled               = ppControl->XmlReadBoolean("me:testMicroRates");
        m_Flags.grainBoltzmannEnabled          = ppControl->XmlReadBoolean("me:printGrainBoltzmann");
        m_Flags.grainDOSEnabled                = ppControl->XmlReadBoolean("me:printGrainDOS");
        m_Flags.grainTSsosEnabled              = ppControl->XmlReadBoolean("me:printTSsos");
        m_Flags.cellDOSEnabled                 = ppControl->XmlReadBoolean("me:printCellDOS");
        m_Flags.reactionOCSEnabled             = ppControl->XmlReadBoolean("me:printReactionOperatorColumnSums");
        m_Flags.kfEGrainsEnabled               = ppControl->XmlReadBoolean("me:printGrainkfE");
        m_Flags.kbEGrainsEnabled               = ppControl->XmlReadBoolean("me:printGrainkbE");
        // Both Tunnelling and Tunneling will work
        m_Flags.TunnellingCoeffEnabled         = ppControl->XmlReadBoolean("me:printTunnellingCoefficients");
        if (!m_Flags.TunnellingCoeffEnabled)
          m_Flags.TunnellingCoeffEnabled       = ppControl->XmlReadBoolean("me:printTunnelingCoefficients");
        m_Flags.CrossingCoeffEnabled           = ppControl->XmlReadBoolean("me:printCrossingCoefficients");
        m_Flags.cellFluxEnabled                = ppControl->XmlReadBoolean("me:printCellTransitionStateFlux");
        m_Flags.grainFluxEnabled               = ppControl->XmlReadBoolean("me:printGrainTransitionStateFlux");
        m_Flags.rateCoefficientsOnly           = ppControl->XmlReadBoolean("me:calculateRateCoefficientsOnly");
        m_Flags.useTheSameCellNumber           = ppControl->XmlReadBoolean("me:useTheSameCellNumberForAllConditions");
        m_Flags.grainedProfileEnabled          = ppControl->XmlReadBoolean("me:printGrainedSpeciesProfile");
        m_Flags.speciesProfileEnabled          = ppControl->XmlReadBoolean("me:printSpeciesProfile");
        m_Flags.printPhenomenologicalEvolution = ppControl->XmlReadBoolean("me:printPhenomenologicalEvolution");
        m_Flags.printSinkFluxes                = ppControl->XmlReadBoolean("me:printSinkFluxes");
        m_Flags.InitialDistEnabled             = ppControl->XmlReadBoolean("me:printInitialDistribution");
        m_Flags.viewEvents                     = ppControl->XmlReadBoolean("me:printEventsTimeStamps");
        m_Flags.allowSmallerDEDown             = ppControl->XmlReadBoolean("me:allowSmallerDeltaEDown");
        m_Flags.print_TabbedMatrices           = ppControl->XmlReadBoolean("me:printTabbedMatrices");
        m_Flags.useDOSweightedDT               = ppControl->XmlReadBoolean("me:useDOSweighedDownWardTransition");
        m_Flags.bForceMacroDetailedBalance     = ppControl->XmlReadBoolean("me:ForceMacroDetailedBalance");

        // System configuration information
        if (ppControl->XmlReadBoolean("me:runPlatformDependentPrecisionCheck")) configuration();


        //if (m_Flags.grainedProfileEnabled && (m_Flags.speciesProfileEnabled)){
        //  cinfo << "Turn off grained species profile to prevent disk flooding." << endl;
        //  m_Flags.grainedProfileEnabled = false;
        //}

        const char* txtPCOP = ppControl->XmlReadValue("me:printCollisionOperatorLevel",false);
        if(txtPCOP) {
          istringstream ss(txtPCOP);
          ss >> m_Flags.showCollisionOperator;
        }

        const char* txt = ppControl->XmlReadValue("me:eigenvalues", optional);
        if(txt && strcmp(txt, "all")==0)
          m_Flags.printEigenValuesNum = -1;
        else
          //no longer uses defaults.xml
          m_Flags.printEigenValuesNum = ppControl->XmlReadInteger("me:eigenvalues",optional);

        // Should eigenvectors be printed?
        if (m_Flags.printEigenValuesNum!=0) {
          m_Flags.printEigenVectors = ppControl->XmlReadBoolean("me:printEigenVectors");
        }

        const char* txtROS = ppControl->XmlReadValue("me:printReactionOperatorSize",optional);
        if(txtROS) {
          istringstream sROS(txtROS);
          sROS >> m_Flags.printReactionOperatorNum;
        }

        const char* txtSTI = ppControl->XmlReadValue("me:shortestTimeOfInterest", optional);
        if (txtSTI){
          istringstream ss(txtSTI);
          ss >> m_Flags.shortestTimeOfInterest;
        }

        const char* txtMET = ppControl->XmlReadValue("me:MaximumEvolutionTime", optional);
        if (txtMET){
          istringstream ss(txtMET);
          ss >> m_Flags.maxEvolutionTime;
        }

        PersistPtr pp = ppControl->XmlMoveTo("me:printGrainProfileAtTime");
        if (pp)
          if(!m_collisionOperator.parseDataForGrainProfileAtTime(pp))
            return false;

        //Copy flags to vector and make a fresh m_Flags
        m_FlagsForEachControl.push_back(m_Flags);
        m_Flags = MesmerFlags();
      } 
    }
    if(!Rdouble::SetUpLinkedVars())
      return false;

    return true;
  }

  //
  // Main calculation method.
  //
  void System::executeCalculation()
  {
    unsigned nConditionBlocks = m_ConditionsForEachControl.size()
      - count(m_ConditionsForEachControl.begin(), m_ConditionsForEachControl.end(), (ConditionsManager*)NULL);
    for (unsigned i = 0; i<m_CalcMethodsForEachControl.size(); ++i)
    {
      m_CalcMethod = m_CalcMethodsForEachControl[i];
      assert(m_CalcMethod);
      m_Flags = m_FlagsForEachControl[i];
      cinfo << "Execute calcMethod " << m_CalcMethod->getID();

      if (nConditionBlocks>1)
      {
        if (m_ConditionsForEachControl[i]) //update if non-NULL
          m_pConditionsManager = m_ConditionsForEachControl[i];

        const char a[4][7] = { "first", "second", "third", "fourth" };
        cinfo << " using " << a[min(i, nConditionBlocks)] << " conditions block ";
      }
      cinfo << endl;
      m_CalcMethod->DoCalculation(this);

      //Calls Finish() for each Reaction. Usually does nothing except in AssociationReaction
      m_pReactionManager->finish();
    }
    //Write the energy convention as an attribute on <moleculeList>
    m_pMoleculeManager->WriteEnergyConvention();
  }

  //
  // Begin calculation.
  // over all PT values, constant parameters
  bool System::calculate(double& chiSquare, vector<double> &residuals, bool writeReport)
  {
    // Controls the print-out of grain/cell DOS in each cycle (This is only for source term)
    if (m_Flags.cellDOSEnabled) m_Flags.cyclePrintCellDOS = true;
    if (m_Flags.grainDOSEnabled) m_Flags.cyclePrintGrainDOS = true;

    TimeCount events; unsigned int timeElapsed =0;

    //
    // Reset microcanonical rate re-calculation flag as parameters, such
    // as reaction threshold may have been altered between invocations of
    // this method.
    //
    for (size_t i(0) ; i < m_pReactionManager->size() ; ++i) {
      (*m_pReactionManager)[i]->resetCalcFlag();
    }

    // Find the highest temperature
    m_Env.MaximumTemperature = m_pConditionsManager->getMaxTemperature();


    //---------------------------------------------
    // looping over temperatures and concentrations
    unsigned int calPoint(0);
    //XML output
    //Considered putting this output under each PT pair in me:conditions.
    //But doesn't work with a range of Ps or Ts. So has to have its own section.

    //There will usually be an <analysis> section for every calculate()
    //When fitting set m_Flags.overwriteXmlAnalysis true
    string comment = m_Flags.overwriteXmlAnalysis ?
      "Only selected calculations shown here" : "All calculations shown";
    PersistPtr ppAnalysis = m_ppIOPtr->XmlWriteMainElement("me:analysis", comment, m_Flags.overwriteXmlAnalysis);
    m_ppAnalysis = ppAnalysis;
    if(Rdouble::withRange().size()!=0)
    {
      PersistPtr ppParams = ppAnalysis->XmlWriteElement("me:parameters");
      for(size_t i=0;i!=Rdouble::withRange().size();++i)
      {
        string varname = Rdouble::withRange()[i]->get_varname();
        //The varnames contain ':' or '(' or')' which is incompatible with XML. replace by '-'.
        replace(varname.begin(), varname.end(), ':', '-');
        replace(varname.begin(), varname.end(), '(', '-');
        replace(varname.begin(), varname.end(), ')', '-');
        //Tolerate spaces in varnames
        replace(varname.begin(), varname.end(), ' ', '_');
        stringstream ss;
        ss << *Rdouble::withRange()[i];
        ppParams->XmlWriteAttribute(varname, ss.str());
      }
    }

    stringstream rateCoeffTable ;

    rateCoeffTable << endl ;
    rateCoeffTable << "    Temperature  Concentration    Exp. Coeff.    Cal. Coeff." << endl ;
    rateCoeffTable << endl ;

    chiSquare = 0.0; // reset the value to zero
    residuals.clear() ;

    for (calPoint = 0; calPoint < m_pConditionsManager->getNumPTPoints(); ++calPoint) {

      m_Env.beta = 1.0 / (boltzmann_RCpK * m_pConditionsManager->PTPointTemp(calPoint)) ; 
      m_Env.conc = m_pConditionsManager->PTPointConc(calPoint); // unit of conc: particles per cubic centimeter
      m_Env.bathGasName = m_pConditionsManager->PTPointBathGas(calPoint);
      map<Reaction*,double> excessConcs = m_pConditionsManager->PTPointExcessConcs(calPoint);
      //Set excess Concentrations for all reactions
      for(map<Reaction*,double>::iterator it=excessConcs.begin();it!=excessConcs.end();++it)
        it->first->set_concExcessReactant(it->second);

      if (writeReport) {cinfo << "PT Grid " << calPoint << endl;}
      Precision precision = m_pConditionsManager->PTPointPrecision(calPoint);
      m_collisionOperator.setPrecision(precision) ;

      ctest << "PT Grid " << calPoint << " Condition: conc = " << m_Env.conc << ", temp = " 
        << m_pConditionsManager->PTPointTemp(calPoint);

      switch (precision) {
      case DOUBLE_DOUBLE:
        ctest << ", diagonalization precision: double-double\n{\n"; 
        break;
      case QUAD_DOUBLE: 
        ctest << ", diagonalization precision: quad-double\n{\n"; 
        break;
      default: 
        ctest << ", diagonalization precision: double\n{\n";
      }

      // Build collison matrix for system.
      if (!m_collisionOperator.BuildReactionOperator(m_Env, m_Flags, writeReport))
        throw (std::runtime_error("Failed building system collison operator.")); 

      if (writeReport) {string thisEvent = "Build Collison Operator" ;
      events.setTimeStamp(thisEvent, timeElapsed) ;
      cinfo << thisEvent;
      if(timeElapsed>0)
        cinfo << " -- Time elapsed: " << timeElapsed << " seconds.";
      cinfo << endl;
      }

      if (!m_Flags.rateCoefficientsOnly){

        if (!m_collisionOperator.calculateEquilibriumFractions())
          throw (std::runtime_error("Failed calculating equilibrium fractions.")); 

        // Diagonalise the collision operator.
        m_collisionOperator.diagReactionOperator(m_Flags, m_Env, ppAnalysis) ;

        if (writeReport) {
          string thisEvent = "Diagonalize the Reaction Operator" ;
          events.setTimeStamp(thisEvent, timeElapsed) ;
          cinfo << thisEvent;
          if(timeElapsed>0)
            cinfo << " -- Time elapsed: " << timeElapsed << " seconds.";
          cinfo << endl;
        }

        if (!m_Env.useBasisSetMethod) {

          if (!m_Flags.bIsSystemSecondOrder) {
            PersistPtr ppPopList;
            if(m_Flags.speciesProfileEnabled)
            {
              ppPopList  = ppAnalysis->XmlWriteElement("me:populationList");
              ppPopList->XmlWriteAttribute("T", toString(m_pConditionsManager->PTPointTemp(calPoint)));
              ppPopList->XmlWriteAttribute("conc", toString(m_Env.conc));
            }

            PersistPtr ppAvEList;
            if(m_Flags.grainedProfileEnabled)
            {
              ppAvEList = ppAnalysis->XmlWriteElement("me:avEnergyList");
              ppAvEList->XmlWriteAttribute("T", toString(m_pConditionsManager->PTPointTemp(calPoint)));
              ppAvEList->XmlWriteAttribute("conc", toString(m_Env.conc));
              ppAvEList->XmlWriteAttribute("energyUnits", "kJ/mol");
            }

            // Calculate time-dependent properties.
            m_collisionOperator.timeEvolution(m_Flags, ppAnalysis, ppPopList, ppAvEList);
            if(m_collisionOperator.hasGrainProfileData())
            {
              PersistPtr ppGrainList  = ppAnalysis->XmlWriteElement("me:grainPopulationList");
              ppGrainList->XmlWriteAttribute("T", toString(m_pConditionsManager->PTPointTemp(calPoint)));
              ppGrainList->XmlWriteAttribute("conc", toString(m_Env.conc));
              m_collisionOperator.printGrainProfileAtTime(ppGrainList);
            }
          } else {
			cinfo << "At present it is not possible to print profiles for systems with a second order term." << endl ;
		  }


          // Calculate rate coefficients. 
          PersistPtr ppList = ppAnalysis->XmlWriteElement("me:rateList");
          ppList->XmlWriteAttribute("T", toString(m_pConditionsManager->PTPointTemp(calPoint)));
          ppList->XmlWriteAttribute("conc", toString(m_Env.conc));
          ppList->XmlWriteAttribute("bathGas", m_Env.bathGasName);
          ppList->XmlWriteAttribute("me:units", "s-1");
          qdMatrix mesmerRates(1);
          qdMatrix lossRates(1);
          m_collisionOperator.BartisWidomPhenomenologicalRates(mesmerRates, lossRates, m_Flags, ppList);

          rateCoeffTable << formatFloat(m_pConditionsManager->PTPointTemp(calPoint), 6, 15) ;
          rateCoeffTable << formatFloat(m_pConditionsManager->PTPointConc(calPoint), 6, 15) ;

          // For these conditions calculate the contribution to the chi^2 merit function
          // for any of the experimentally observable data types. 

          chiSquare += calcChiSqRateCoefficients(mesmerRates, calPoint, rateCoeffTable, residuals);

          chiSquare += calcChiSqYields(calPoint, rateCoeffTable, residuals);

          chiSquare += calcChiSqEigenvalues(calPoint, rateCoeffTable, residuals);

          ctest << "}\n";

        } else {

          qdMatrix mesmerRates(1);
          m_collisionOperator.BartisWidomBasisSetRates(mesmerRates, m_Flags);

        }
      }

    } // End of temperature and concentration loop. 

    rateCoeffTable << endl ;

    if (writeReport) cinfo << rateCoeffTable.str() ;

    string thisEvent = "Finish Calculation of P-T case";
    events.setTimeStamp(thisEvent, timeElapsed);
    cinfo << thisEvent;
    if(timeElapsed>0)
      cinfo << " -- Time elapsed: " << timeElapsed << " seconds.";
    cinfo << endl;
    if (writeReport) cwarn << calPoint << " temperature/concentration-pressure points calculated." << endl;

    if (m_Flags.viewEvents) cinfo << events;

    return true;
  }

  //
  // Calculate rate coefficients etc. for a specific condition.
  //
  bool System::calculate(size_t nConds, vector<double> &Quantities, bool writeReport )
  {

    //
    // Reset microcanonical rate re-calculation flag as parameters, such
    // as reaction threshold may have been altered between invocations of
    // this method.
    //
    for (size_t i(0) ; i < m_pReactionManager->size() ; ++i) {
      (*m_pReactionManager)[i]->resetCalcFlag();
    }

    m_Flags.printEigenValuesNum = 0 ;

    double temp = m_pConditionsManager->PTPointTemp(nConds) ; 
    m_Env.beta = 1.0 / (boltzmann_RCpK * temp) ;

    // unit of conc: particles per cubic centimetre
    m_Env.conc = m_pConditionsManager->PTPointConc(nConds);

    m_Env.bathGasName = m_pConditionsManager->PTPointBathGas(nConds);

    Precision precision = m_pConditionsManager->PTPointPrecision(nConds);
    m_collisionOperator.setPrecision(precision) ;

    // Build collison matrix for system.
    if (!m_collisionOperator.BuildReactionOperator(m_Env, m_Flags))
      throw (std::runtime_error("Failed building system collison operator.")); 

    if (!m_collisionOperator.calculateEquilibriumFractions())
      throw (std::runtime_error("Failed calculating equilibrium fractions.")); 

    // Diagonalise the collision operator.
    m_collisionOperator.diagReactionOperator(m_Flags, m_Env) ;

    // Calculate rate coefficients.
    qdMatrix  mesmerRates(1) ;
    qdMatrix lossRates(1);
    m_collisionOperator.BartisWidomPhenomenologicalRates(mesmerRates, lossRates, m_Flags);

    // For this conditions calculate the experimentally observable data types. 

    stringstream rateCoeffTable ;
    vector<double> residuals ;
    calcChiSqRateCoefficients(mesmerRates, nConds, rateCoeffTable, residuals);

    calcChiSqYields(nConds, rateCoeffTable, residuals);

    calcChiSqEigenvalues(nConds, rateCoeffTable, residuals);

    double data(0.0) ;
    while (rateCoeffTable >> data) {
      Quantities.push_back(data) ;
    }

    return true ;
  }

  //
  // Calculate rate coefficients for conditions other than those 
  // supplied directly by user e.g. for analytical representation.
  //
  bool System::calculate(double &Temperature, double &Concentration, Precision precision,
    map<string, double> &phenRates, double &MaxT, const string& bathGas)
  {
    assert(!bathGas.empty());
    m_Env.bathGasName = bathGas;
    m_Flags.printEigenValuesNum = 0 ;

    m_Env.beta = 1.0 / (boltzmann_RCpK * Temperature) ;
    m_Env.MaximumTemperature = MaxT;
    m_Env.conc = Concentration; // Unit of conc: particles per cubic centimeter.

    m_collisionOperator.setPrecision(precision) ;

    // Build collison matrix for system.
    if (!m_collisionOperator.BuildReactionOperator(m_Env, m_Flags))
      throw (std::runtime_error("Failed building system collison operator.")); 

    if (!m_collisionOperator.calculateEquilibriumFractions())
      throw (std::runtime_error("Failed calculating equilibrium fractions.")); 

    // Diagonalise the collision operator.
    m_collisionOperator.diagReactionOperator(m_Flags, m_Env) ;

    // Calculate rate coefficients.
    qdMatrix BWrates(1);
    qdMatrix lossRates(1);
    m_collisionOperator.BartisWidomPhenomenologicalRates(BWrates, lossRates, m_Flags);
    m_collisionOperator.get_phenRates(phenRates);

    return true ;
  }

  double System::calcChiSqRateCoefficients(const qdMatrix& mesmerRates, const unsigned calPoint, stringstream &rateCoeffTable, vector<double> &residuals){

    double chiSquare(0.0) ;

    vector<conditionSet> expRates;
    m_pConditionsManager->get_experimentalRates(calPoint, expRates);
    for (size_t i(0); i < expRates.size(); ++i){

      string ref1, ref2, refReaction; 
      double expRate(0.0), expErr(0.0); 
      expRates[i].get_conditionSet(ref1, ref2, refReaction, expRate, expErr);

      // Get the position of ref1 and ref2 inside m_SpeciesSequence vector
      int seqMatrixLoc1(-1), seqMatrixLoc2(-1);
      seqMatrixLoc1 = m_collisionOperator.getSpeciesSequenceIndex(ref1);
      seqMatrixLoc2 = m_collisionOperator.getSpeciesSequenceIndex(ref2);

      if(seqMatrixLoc1<0 || seqMatrixLoc2<0)
        throw(std::runtime_error("Failed to locate species in rate coefficient matrix.")) ;

      // 
      // In the following it is assumed that experimental rate coefficients will always 
      // be quoted as a absolute values. Since the diagonal values of the BW matrix are
      // negative, their absolute value is required for comparision with experimental
      // values hence the fabs invocation.
      //
      double rateCoeff = fabs(to_double(mesmerRates[seqMatrixLoc2][seqMatrixLoc1])) ;

      // Is a bimolecular rate coefficient required?

      Reaction *reaction = m_pReactionManager->find(refReaction) ;
      ReactionType reactionType = UNDEFINED_REACTION ;
      if (reaction)
        reactionType = reaction->getReactionType() ;
      if (reactionType == ASSOCIATION || reactionType == PSEUDOISOMERIZATION || reactionType == SECONDORDERASSOCIATION) {
        double concExcessReactant = reaction->get_concExcessReactant() ;

        // Test concentration and reaction sense.

        if (concExcessReactant > 0.0 && (reaction->get_reactant()->getName() == ref1)) {
          rateCoeff /= concExcessReactant ;
        }

        // Second order reactions need to account for a geometric factor

        if (reactionType == SECONDORDERASSOCIATION) {
          rateCoeff /= 4.0 ;
        }

      } else {
        // No reference reaction. Assume reaction is unimolecular.
      }

      double diff = (expRate - rateCoeff)/expErr ;
      residuals.push_back(diff) ;
      chiSquare +=  diff * diff ;

      rateCoeffTable << formatFloat(expRate, 6, 15) << formatFloat(rateCoeff, 6, 15) << endl ;
      AddCalcValToXml(calPoint, i, rateCoeff);

    }
    return chiSquare;
  }

  double System::calcChiSqYields(const unsigned calPoint, stringstream &rateCoeffTable, vector<double> &residuals) {

    double chiSquare(0.0) ;
    vector<conditionSet> expYields;
    m_pConditionsManager->get_experimentalYields(calPoint, expYields);
    if (expYields.size() == 0)
      return chiSquare ;

    //
    // Calculate yields for these conditions. Assume all yields are measured at the same time.
    //
    YieldMap yieldMap ;
    try {
      string product, yieldTime, ref; 
      double expYield(0.0), expErr(0.0); 
      expYields[0].get_conditionSet(product, yieldTime, ref, expYield, expErr);

      double time ;
      stringstream s1(yieldTime); s1 >> time ;

      m_collisionOperator.calculateYields(yieldMap, time) ;
    } catch (std::runtime_error& e) {
      cerr << "Error: during calculation of Chi^2 for yields:" << endl ;
      cerr << e.what() << endl ;
      return 0.0 ;
    }

    for (size_t i(0); i < expYields.size(); ++i){

      string ref1, ref2, refReaction; 
      double expYield(0.0), expErr(0.0); 
      expYields[i].get_conditionSet(ref1, ref2, refReaction, expYield, expErr);

      double yield(0.0) ; 
      YieldMap::const_iterator yielditr = yieldMap.begin();
      for (; yielditr != yieldMap.end() ; ++yielditr) {

        Reaction* sinkReaction = yielditr->first ;
        vector<Molecule*> pdts;
        sinkReaction->get_products(pdts);
        for (size_t i(0) ; i < pdts.size() ; i++) {
          if (ref1 == pdts[i]->getName())
            yield += yielditr->second ;
        }

      }

      double diff = (expYield - yield)/expErr ;
      residuals.push_back(diff) ;
      chiSquare += (diff * diff);

      rateCoeffTable << formatFloat(expYield, 6, 15) << formatFloat(yield, 6, 15) << endl ;
      AddCalcValToXml(calPoint, i, yield);
    }

    return chiSquare;
  }

  double System::calcChiSqEigenvalues(const unsigned calPoint, stringstream &rateCoeffTable, vector<double> &residuals) {

    double chiSquare(0.0) ;

    vector<conditionSet> expEigenvalues;
    m_pConditionsManager->get_experimentalEigenvalues(calPoint, expEigenvalues);
    for (size_t i(0); i < expEigenvalues.size(); ++i){

      string ref1, ref2, strIdEigenvalue; 
      double expEigenvalue(0.0), expErr(0.0); 
      expEigenvalues[i].get_conditionSet(ref1, ref2, strIdEigenvalue, expEigenvalue, expErr);

      size_t idEigenvalue ;
      stringstream s1(strIdEigenvalue); s1 >> idEigenvalue ;

      double eigenvalue(0.0) ;
      try {
        eigenvalue = m_collisionOperator.getEigenvalue(idEigenvalue) ;
      } catch (std::runtime_error& e) {
        cerr << e.what() << endl ;
        continue ;
      }

      double diff = (expEigenvalue - eigenvalue)/expErr ;
      residuals.push_back(diff) ;
      chiSquare += (diff * diff); 

      rateCoeffTable << formatFloat(expEigenvalue, 6, 15) << formatFloat(eigenvalue, 6, 15) << endl ;
      AddCalcValToXml(calPoint, i, eigenvalue);
    }

    return chiSquare;
  }

  void System::AddCalcValToXml(const unsigned calPoint, size_t i, double val) const
  {
    // Add extra attribute(s) containing calculated value and timestamp to <me:experimentalRate> (or similar element).
    PersistPtr pp = m_pConditionsManager->get_experimentalDataPtr(calPoint, i);
    TimeCount events;
    string timeString;
    pp->XmlWriteAttribute("calculated", events.setTimeStamp(timeString));
    stringstream ss;
    ss << val;
    pp->XmlWriteAttribute("calcVal", ss.str());
    pp->XmlReadDouble("calcVal",false);
  }

  void System::WriteMetadata(const string& infilename)
  {
    PersistPtr ppList = m_ppIOPtr->XmlWriteMainElement("metadataList", "");
    ppList->XmlWriteAttribute("xmlns:dc", "http://purl.org/dc/elements/1.1/");
    if(m_pTitle)
      ppList->XmlWriteValueElement("dc:title", m_pTitle);
    if(m_pDescription)
      ppList->XmlWriteValueElement("dc:description", m_pDescription);
    ppList->XmlWriteValueElement("dc:source", infilename);

    ppList->XmlWriteValueElement("dc:creator","Mesmer v"+string(MESMER_VERSION));
    TimeCount events;
    string timeString = events.setTimeStamp("");
    cinfo << "Write metadata " << timeString << endl;
    ppList->XmlWriteValueElement("dc:date", timeString);

    const char* author = getenv("MESMER_AUTHOR");
    if(!author)
      author = getenv("USERNAME"); //Windows
    if(!author)
      author = getenv("LOGNAME"); //Unix?
    if(!author)
      author = "unknown";
    ppList->XmlWriteValueElement("dc:contributor", author);
  }

  void System::getAllBathGases(set<string>& bathGases)
  {
    for(size_t i(0) ; i!=m_ConditionsForEachControl.size() ; ++i)
      if(m_ConditionsForEachControl[i])
        m_ConditionsForEachControl[i]->getAllBathGases(bathGases) ;

    bathGases.insert(getMoleculeManager()->get_BathGasName());
  }


  void System::configuration(void){
    clog << "\nPrinting system precision configuration:" << endl;
    clog << "Size of float = " << sizeof(float) << endl;
    clog << "Size of double = " << sizeof(double) << endl;
    clog << "Size of long double = " << sizeof(long double) << endl;
    clog << "Size of double-double = " << sizeof(dd_real) << endl;
    clog << "Size of quad-double = " << sizeof(qd_real) << endl;

    clog << "\nEpsilon is the difference between 1 and the smallest value greater than 1 that is representable for the data type." << endl;
    clog << "float epsilon == " << numeric_limits<float>::epsilon() << endl;
    clog << "double epsilon == " << numeric_limits<double>::epsilon() << endl;
    clog << "long double epsilon == " << numeric_limits<long double>::epsilon() << endl;
    clog << "dd_real epsilon == " << numeric_limits<dd_real>::epsilon() << endl;
    clog << "qd_real epsilon == " << numeric_limits<qd_real>::epsilon() << endl;

    clog << "\nfloat max == " << numeric_limits<float>::max() << endl;
    clog << "double max == " << numeric_limits<double>::max() << endl;
    clog << "long double max == " << numeric_limits<long double>::max() << endl;
    clog << "dd_real max == " << numeric_limits<dd_real>::max() << endl;
    clog << "qd_real max == " << numeric_limits<qd_real>::max() << endl;

    clog << "\nfloat min == " << numeric_limits<float>::min() << endl;
    clog << "double min == " << numeric_limits<double>::min() << endl;
    clog << "long double min == " << numeric_limits<long double>::min() << endl;
    clog << "dd_real min == " << numeric_limits<dd_real>::min() << endl;
    clog << "qd_real min == " << numeric_limits<qd_real>::min() << endl;
  }

  /* 
  System::parse() loops to read all <control> blocks.
  CalcMethod is pushed to a vector<CalcMethod> calcMethods
  Make a MesmerFlags, populate it, push to vector<MesmerFlags> flags

  System::executeCalculation()
  for loop over size of calcMethods
  m_CalcMethod set to current member
  copy from flags vector to m_Flags
  */
} //namespace

