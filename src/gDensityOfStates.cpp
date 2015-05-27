#include <cmath>
#include <iomanip>
#include <numeric>
#include <string>
#include "Molecule.h"
#include "System.h"
#include "ParseForPlugin.h"
#include "gStructure.h"
#include "gDensityOfStates.h"

using namespace std;
using namespace Constants;
using namespace OpenBabel;
namespace mesmer
{

  //-------------------------------------------------------------------------------------------------
  // Cell density of states related properties
  //-------------------------------------------------------------------------------------------------
  string MolecularComponent::m_energyConvention;

  //
  // Constructor, destructor and initialization
  //
  gDensityOfStates::~gDensityOfStates()
  {
    //Delete the density of state calculators because they are cloned instances
    for (unsigned i = 0; i < m_DOSCalculators.size(); ++i)
      delete m_DOSCalculators[i];

    if (m_Hessian)
      delete m_Hessian;
    if (m_Modes)
      delete m_Modes;
  }

  gDensityOfStates::gDensityOfStates(Molecule* pMol)
    :m_RotCstA(0.0),
    m_RotCstB(0.0),
    m_RotCstC(0.0),
    m_Sym(1.0),
    m_ZPE(NaN,"cm-1"),
    m_scaleFactor(1.0),
    m_SpinMultiplicity(1),
    m_RC_chk(-1),
    m_Sym_chk(-1),
    m_ZPE_chk(-1),
    m_scaleFactor_chk(-1),
    m_SpinMultiplicity_chk(-1),
    m_eleExc(),
    m_VibFreq(),
    m_Hessian(NULL),
    m_Modes(NULL),
    m_nModes(0),
    m_HessianUnits(),
    m_MaximumCell(0),
    m_cellDOS(),
    m_grainEne(),
    m_grainDOS() {
    m_host = pMol;
  }

  bool gDensityOfStates::initialization() {

    // Define basic molecular parameters. 

    Molecule* pMol = m_host;

    ErrorContext c(pMol->getName());

    PersistPtr pp = pMol->get_PersistentPointer();

    PersistPtr ppPropList = pp->XmlMoveTo("propertyList");
    if (!ppPropList)
      ppPropList = pp; // A propertyList element is not essential.

    // Vibrational frequencies. Test for Hessain first and, if absent,
    // try to read frequencies.

    bool hasVibFreq(true);
    const char *txt;
    if ( (m_Hessian = ReadPropertyMatrix<double>("me:hessian", ppPropList)) ){
      PersistPtr pMtrx = ppPropList->XmlMoveToProperty("me:hessian");
      txt = pMtrx->XmlReadValue("units", false);
      m_HessianUnits = (txt) ? string(txt) : "kJ/mol/amu/Ang^2";
      FrqsFromHessian();
    }
    else if ((txt = ppPropList->XmlReadProperty("me:vibFreqs", optional))) {
      istringstream idata(txt);
      double x;
      while (idata >> x)
        m_VibFreq.push_back(x);
    }
    else {
      hasVibFreq = false;
      if (!pMol->getStruc().IsAtom())
        cinfo << "Cannot find me:vibFreqs or me:hessian. Assuming an atom or a sink molecule." << endl;
    }

    m_scaleFactor = ppPropList->XmlReadPropertyDouble("me:frequenciesScaleFactor");
    m_scaleFactor_chk = 0;

    // Rotational constants.

    std::vector<double> rCnst(3, 0.0);
    bool hasRotConst(false);
    txt = ppPropList->XmlReadProperty("me:rotConsts", optional);
    if (txt) {
      //data from <me:rotConsts>
      istringstream idata(txt);
      idata >> rCnst[0] >> rCnst[1] >> rCnst[2];
      txt = ppPropList->XmlReadPropertyAttribute("me:rotConsts", "units", optional);
      if (string(txt) == "amuA^2") {
        // The supplied data are moments of inertia with units amuAng^2.
        rCnst[0] = conMntInt2RotCnt / rCnst[0];
        rCnst[1] = conMntInt2RotCnt / rCnst[1];
        rCnst[2] = conMntInt2RotCnt / rCnst[2];
      }
      hasRotConst = true;
      m_RC_chk = 0;
    }
    else {
      // Attempt to calculate rotational constants from atomic coordinates.
      gStructure& gs = pMol->getStruc();
      if (gs.ReadStructure()) {
        rCnst = gs.CalcRotConsts();
        cinfo << "Rotational constants were calculated from atom coordinates: "
          << rCnst[2] << ' ' << rCnst[1] << ' ' << rCnst[0] << " cm-1" << endl;
        hasRotConst = true;
      }
    }
    if (hasRotConst) {
      // Check rotational constants are valid.
      if (rCnst[0] < 0.0 || rCnst[1] < 0.0 || rCnst[2] < 0.0) {
        cerr << "A negative rotational constant found " << endl;
        throw (std::runtime_error("Fatal error"));
      }
      // Make sure the rotational constants are in ascending order.
      std::sort(rCnst.begin(), rCnst.end());
      m_RotCstA = rCnst[2];
      m_RotCstB = rCnst[1];
      m_RotCstC = rCnst[0];
      m_RC_chk = 0;
    }
    else if (!pMol->getStruc().IsAtom()){
      cinfo << "No rotational constants from <me:rotConsts> or structure. "
        "Assuming an atom or a sink molecule." << endl;
    }
    else
      m_RC_chk = 0;

    m_Sym = ppPropList->XmlReadPropertyDouble("me:symmetryNumber");
    m_Sym_chk = 0;

    if (hasVibFreq != hasRotConst){
      cerr << "Improper setting on vibrational frequencies or rotational constants." << endl;
    }

    // Electronic excitations.

    txt = ppPropList->XmlReadProperty("me:electronicExcitation", optional);
    if (txt) {
      istringstream idata(txt);
      m_eleExc.clear();
      double _iele;
      while (idata >> _iele) m_eleExc.push_back(_iele);
    }

    // Spin multiplicity.

    m_SpinMultiplicity = ppPropList->XmlReadPropertyInteger("me:spinMultiplicity", optional);
    if (m_SpinMultiplicity == 0)
      m_SpinMultiplicity = pp->XmlReadInteger("spinMultiplicity");
    m_SpinMultiplicity_chk = 0;

    if (!ReadDOSMethods())
      throw(std::runtime_error(""));

    // Check whether the molecule has the correct number of degrees of freedom.
    if (!m_host->checkDegOfFreedom()){
      string errorMsg = "Incorrect number of degrees of freedom compared with atom count for " + m_host->getName();
      throw (std::runtime_error(errorMsg));
    }

    return ReadZeroPointEnergy(ppPropList);
  }

  bool gDensityOfStates::ReadZeroPointEnergy(PersistPtr ppPropList)
  {
    /* 
    The convention attribute on a <me:ZPE> element of a molecule describes
    what zero energy means. Any string can be used for a convention but the
    built-in ones are:
    arbitrary       - energy zero is chosen by the user, often the ZPE of the
                      lowest species. This is the default if there is no
                      convention attribute of me:ZPE.
    computational   - energy zero is separated nuclei and electrons from computational 
                      chemistry programs (but see below for ZPE correction).
    thermodynamic   - energy zero is the heat of formation at 0K
                      (with a 0K reference temperature).

    All molecules taking part in the reactions must have the same convention,
    which is set either by a convention attribute on <moleculeList> or 
    by the first molecule in the data file that has an energy
    specified. (Bath gas molecules and sink molecules may not have a
    specified energy.) It is stored (as a static variable) in
    MolecularComponent::m_energyConvention.
    
    When reading a datafile the energy of a molecule is looked for 
    first in <me:ZPE>, and if not found, successively in <me:Hf0>
    and <me:Hf298>. If the convention has not been set in <moleculeList>,
    it is set by the convention on <me:ZPE>, if there is one, and to arbitrary
    otherwise. If there is a <me:Hf0> 0r <me:Hf298> rather than <me:ZPE> the
    convention is set to thermodynamic. Once the convention is set, the value
    of any <me:ZPE> with a convention attribute and any <me:Hf0> or <me:Hf298>
    is converted by Mesmer to the previously set energy convention. A convention
    thermodynamic298K is used in this process but is for internal use only.

    When the energy for a molecule comes from <me:Hf0> or <me:Hf298>, the
    output XML file will have a new <me:ZPE> element with the energy in the set
    convention (possibly converted) and an attribute giving the original element.
    This means that if this file is subsequently used as an input, all the
    energies are read from the <me:ZPE> elements, with no conversion being necessary.

    The safest practice is to add a convention attribute to <moleculeList>.

    But if the energies of all the molecules are of the same type there is no
    need for any convention attributes. The commonest case is when they are
    all provided in <me:ZPE> elements, when the convention is arbitrary. 
    If they are all in either <me:Hf0> or <me:Hf298> the convention is
    thermodynamic.

    Only if there are both computational and thermodynamic energies is
    an attribute, convention="computational", necessary on each <me:ZPE>.

    The conversion to and from computational is done via the
    properties of the dissociated atoms of the molecule and may not be
    as accurate as alternative calculation methods via molecules with
    similar bonding.

    me:HfAT0 is no longer supported.
    */

    // Look successively for alternative molecular energy inputs
    return ReadEnergy("me:ZPE", "computational") ||
      ReadEnergy("me:Hf0", "thermodynamic") ||
      ReadEnergy("me:Hf298", "thermodynamic298K");
      ReadEnergy("DHf(0K)", "thermodynamic") || /*OpenBabel styles*/
      ReadEnergy("DHf(298.15K)", "thermodynamic298K");
   }

  double gDensityOfStates::ConvertEnergyConvention(
    const std::string& fromConvention, const std::string& toConvention, double fromValue, const string& units)
  {
    // Converts between "computational", "thermodynamic"  and "thermodynamic298K" conventions
    // The last parameter specifies the units of the input and returned values.
    if (IsNan(fromValue))
      return fromValue;
    if (fromConvention == toConvention)
      return fromValue;

    ErrorContext c(getHost()->getName());
    double H, S, G;
    thermodynamicsFunctions(298, 1.0 / kJPerMol_in_RC, H, S, G); //kJ/mol
    double atomZPE, atomHf0, atomHf298, atomdH298, stddH298;
    if (!(getHost()->getStruc()).GetConstituentAtomThermo(atomZPE, atomHf0, atomHf298, atomdH298,stddH298))
    {
      cerr << "Missing thermo data for atoms" << endl;
      return false;
    }

      /*     X(298K)  =>    ElemsStdStates(298K)      Hf(298K)
      H298(X)  |                  |    H298 for elems
             X(0K)    =>    ElemsStdStates(0K)        Hf(0K)
      compE    |                  |    Hf0 for atoms
         nuclei,electrons =>  Atoms(0K)               compEatoms

       Hf(298K)  +  H298els  -  Hf(0K)  -  H298 = 0     
       compE -compEatoms - Hf(0K) + Hfatoms(0K) = 0   */

    //Use Hf0 as intermediate
    double val = ConvertEnergy(units, "kJ/mol", fromValue);
    if (fromConvention == "thermodynamic298K")
      val += stddH298 - H; //all in kJ/mol
    else if (fromConvention == "computational")
      val += atomHf0 - atomZPE ;
    else if (fromConvention != "thermodynamic")
      val = NaN;

    if (toConvention == "thermodynamic298K")
      val += H - stddH298; //all in kJ/mol
    else if (toConvention == "computational")
      val += atomZPE - atomHf0;
    else if (toConvention != "thermodynamic")
      val = NaN;

    if (IsNan(val))
      cerr << "energy cannot be converted from "
           << fromConvention << " to " << toConvention << endl;
    else
      cinfo <<"energy converted from "
            << fromConvention << " to " << toConvention << endl;

    return ConvertEnergy("kJ/mol", units, val);
  }

  bool gDensityOfStates::ReadEnergy(std::string elName, std::string nativeConvention)
  {
    string unitsInput;
    PersistPtr ppMol = getHost()->get_PersistentPointer();
    PersistPtr ppPropList = ppMol->XmlMoveTo("propertyList");
    if (!ppPropList) //do without <propertyList>
      ppPropList = ppMol;
    double tempzpe = ppMol->XmlReadPropertyDouble(elName, optional);
    if (IsNan(tempzpe))
      return false; //element elName not found

    //Read units; default from defaults.xml, probably kJ/mol.
    unitsInput = ppPropList->XmlReadPropertyAttribute(elName, "units");

    if (elName != "me:ZPE")
    {
      //If not set already, set the energy convention for all molecules
      if (m_energyConvention.empty())
      {
        m_energyConvention = "thermodynamic";
        cinfo << "The energy convention is " << m_energyConvention
              << ", derived from " << elName << endl;
      }
    }
    else //<me:ZPE>
    {
      const char* txt = ppPropList->XmlReadPropertyAttribute(elName, "convention", optional);
      //If not set already, set the energy convention for all molecules
      if (m_energyConvention.empty())
        m_energyConvention = txt ? txt : "arbitrary";

      double zpCorrection = 0.0; //cm-1
      if (m_energyConvention == "computational" || (txt && strcmp(txt, "computational")==0))
      {
        // Mesmer assumes that the ZPE is the true zero point energy, unless
        // an attribute zeroPointVibEnergyAdded="false" 
        // (or me:zeroPointVibEnergyAdded="false", an old version which also works with the schema)
        // is present. This indicates that the value provided is 
        // that the value provided is a computed energy at the bottom of the well
        // and it is corrected here by adding 0.5*Sum(vib freqs).
        string zpAttName("zeroPointVibEnergyAdded");
        txt = ppPropList->XmlReadPropertyAttribute(elName, zpAttName, optional);
        if (!txt)
        {
          zpAttName = "me:" + zpAttName; //possible older version
          txt = ppPropList->XmlReadPropertyAttribute(elName, zpAttName, optional);
        }
        if (txt && strcmp(txt, "false") == 0)
        {
          if (m_VibFreq.size() > 0)
          {
            zpCorrection = accumulate(m_VibFreq.begin(), m_VibFreq.end(), 0.0);
            zpCorrection *= 0.5;
          }
          tempzpe += ConvertFromWavenumbers(unitsInput, zpCorrection);
          //Write back a corrected value and change attribute to zeroPointVibEnergyAdded="true"
          PersistPtr ppScalar = ppMol->XmlMoveToProperty(elName);
          ppScalar->XmlWrite(toString(tempzpe)); //in original units
          ppScalar->XmlWriteAttribute(zpAttName, "true");
          cinfo << "Zero point correction made by adding " << zpCorrection << " cm-1" << endl;
        }
      }
      bool rangeSet;
      PersistPtr ppProp = ppMol->XmlMoveToProperty(elName);
      ReadRdoubleRange(m_host->getName() + elName.substr(2), ppProp, m_ZPE, rangeSet,
        getConvertedEnergy(unitsInput, 1.0), zpCorrection);
    }

    if (m_energyConvention != "arbitrary")
      tempzpe = ConvertEnergyConvention(nativeConvention, m_energyConvention, tempzpe, unitsInput);

    m_ZPE = getConvertedEnergy(unitsInput, tempzpe);

    if (elName != "me:ZPE")
    {
      // Write the value converted from Hf0 or Hf298 back to a me:ZPE element in the
      // XML file, mainly so it can be displayed in the Firefox energy level diagram
      // and to allow the energy to be a range variable. 
      // <me:Hf0 and <me:Hf298> cannot have range attributes.
      stringstream ss;
      ss << ConvertFromWavenumbers(unitsInput, m_ZPE) ;
      PersistPtr ppScalar = ppPropList->XmlWriteProperty("me:ZPE", ss.str(), unitsInput);
      ppScalar->XmlWriteAttribute("source", elName);
      cinfo << "New me:ZPE element written with data from " << elName << endl;
    }
    return true;
  }

  double gDensityOfStates::get_Hf298Thermo() //kJ/mol
  {
    if (IsNan(m_ZPE))
      return NaN;
    cinfo << "enthalpy of formation at 298K:";
    return ConvertEnergyConvention(m_energyConvention, "thermodynamic298K",
                                   ConvertFromWavenumbers("kJ/mol", m_ZPE),"kJ/mol");
  }

  //
  // Get the number of degrees of freedom for this species.
  //
  unsigned int gDensityOfStates::getNoOfDegOfFreeedom() {
    unsigned int nDOF(0);
    for (vector<DensityOfStatesCalculator*>::size_type j = 0; j < m_DOSCalculators.size(); ++j) {
      nDOF += m_DOSCalculators[j]->NoDegOfFreedom(this);
    }
    return nDOF;
  }

  //
  // Get cell density of states.
  //
  bool gDensityOfStates::getCellDensityOfStates(vector<double> &cellDOS, int startingCell, bool bcalc) {
    // If density of states have not already been calcualted then do so.
    if (bcalc && !calcDensityOfStates())
    {
      cerr << "Failed calculating DOS" << endl;
      return false;
    }
    if (startingCell == 0)
      cellDOS = m_cellDOS;
    else{
      int MaximumCell = m_host->getEnv().MaxCell;
      for (int i(startingCell); i < MaximumCell; ++i){
        cellDOS.push_back(m_cellDOS[i]);
      }
    }
    return true;
  }

  bool gDensityOfStates::ReadDOSMethods() {
    // There must be only one <me:DOSCMethod> specifying a method which
    // includes the rotations. If not present the default is used.
    // There can be multiple additional methods:
    // <me:ExtraDOSCMethod name="..."> optional parameters </me:ExtraDOSCMethod>

    ErrorContext c(getHost()->getName());
    PersistPtr pp = getHost()->get_PersistentPointer();

    //Get the method which includes rotations, or use the default
    const char* name = pp->XmlReadValue("me:DOSCMethod");
    if (!*name) // Must be alt form e.g. <me:DOSCMethod name="QMRotors"/>
      name = pp->XmlMoveTo("me:DOSCMethod")->XmlReadValue("name");

    m_DOSCalculators.push_back(DensityOfStatesCalculator::Find(string(name)));
    if (!m_DOSCalculators[0] || !m_DOSCalculators[0]->ReadParameters(this, pp))
      return false; //error message already output
    if (!m_DOSCalculators[0]->includesRotations())
    {
      cerr << "The calculator specified in <me:DOSCMethod>"
        " should be one that includes the rotations."
        " Use  <me:ExtraDOSCMethod name=\""
        << m_DOSCalculators[0]->getID() << "\"> </me:ExtraDOSCMethod> instead."
        << endl;
      return false;
    }

    // Beyer-Swinehart object added by default at m_DOSCalculators[1]
    m_DOSCalculators.push_back(DensityOfStatesCalculator::Find("BeyerSwinehart"));
    if (!m_DOSCalculators[1])
    {
      cerr << "Beyer-Swinehart algorithm failed to initialize correctly" << endl;
      return false;
    }

    //Read any additional methods
    PersistPtr pp2 = pp;
    while (pp2 = pp2->XmlMoveTo("me:ExtraDOSCMethod"))
    {
      string dosMethod;
      const char* name;
      if ((name = pp2->XmlRead())
        || (name = pp2->XmlReadValue("name", optional))
        || (name = pp2->XmlReadValue("xsi:type", optional)))
      {
        if (name[2] == ':')
          name += 3; //Remove prefix "me:"
        dosMethod = name;
      }

      DensityOfStatesCalculator* pDOSCalculator = DensityOfStatesCalculator::Find(dosMethod);
      if (!pDOSCalculator || !pDOSCalculator->ReadParameters(this, pp2))
        return false;
      m_DOSCalculators.push_back(pDOSCalculator);
    }

    //Check there is only one <me:DOSCMethod>
    PersistPtr pp1 = pp->XmlMoveTo("me:DOSCMethod"); //to the element just found
    if (pp1->XmlMoveTo("me:DOSCMethod"))
      cerr << "Too many <me:DOSCMethod> elements on this molecule. "
      << "Only the first is used. Additional methods should be under <me:ExtraDOSCMethod>." << endl;
    return true;
  }

  bool gDensityOfStates::RemoveDOSCalculator(const string& id)
  {
    vector<DensityOfStatesCalculator*>::iterator iter;
    for (iter = m_DOSCalculators.begin(); iter != m_DOSCalculators.end(); ++iter)
    {
      if (id == (*iter)->getID())
      {
        delete *iter; //because plugin was a new instance made with Clone()
        m_DOSCalculators.erase(iter);
        return true;
      }
    }
    return false;
  }
  bool gDensityOfStates::AddDOSCalculator(const string& id)
  {
    m_DOSCalculators.push_back(DensityOfStatesCalculator::Find(id));
    return m_DOSCalculators.back();
  }

  DensityOfStatesCalculator* gDensityOfStates::GetDOSCalculator(const string& id)
  {
    vector<DensityOfStatesCalculator*>::iterator iter;
    for (iter = m_DOSCalculators.begin(); iter != m_DOSCalculators.end(); ++iter)
    {
      if (id == (*iter)->getID())
        return *iter;
    }
    return NULL;
  }

  RotationalTop gDensityOfStates::test_rotConsts()
  {
    std::vector<double> mmtsInt;
    return get_rotConsts(mmtsInt);
  }

  RotationalTop gDensityOfStates::get_rotConsts(std::vector<double> &mmtsInt)
  {
    if (m_RC_chk <= -1){
      ErrorContext e(this->getHost()->getName());
      if (m_RC_chk == -1)
        //        cinfo << "Rotational constants were not defined but requested." << endl;
      --m_RC_chk;
      return UNDEFINED_TOP; // treat as a non-rotor
    }

    mmtsInt.clear();
    mmtsInt.push_back(m_RotCstA);
    mmtsInt.push_back(m_RotCstB);
    mmtsInt.push_back(m_RotCstC);
    ++m_RC_chk;

    // The classification of rotors is simplified to only three following types.
    // 3-D rotors may have other attributes, but in ILT they are treated as the same type. 

    if ((mmtsInt[0] + mmtsInt[1] + mmtsInt[2]) == 0.)
      return UNDEFINED_TOP; // not a rotor
    else if ((mmtsInt[0] * mmtsInt[1] * mmtsInt[2]) == 0.)
      return LINEAR;        // 2-D linear
    else
      return NONLINEAR;     // 3-D symmetric/asymmetric/spherical top
  }


  //
  // Calculate the rovibrational density of states.
  //
  bool gDensityOfStates::calcDensityOfStates()
  {
    bool recalc(false);
    const size_t MaximumCell = m_host->getEnv().MaxCell;
    if (MaximumCell > m_MaximumCell) {
      recalc = true;
      m_MaximumCell = MaximumCell;
    }

    if (recalc) {
      // Calculate density of states.
      bool ret(true);
      for (size_t i(0); ret && i < m_DOSCalculators.size(); ++i)
        ret = ret && m_DOSCalculators[i]->countCellDOS(this, m_host->getEnv());
      if (!ret)
        return false;
    }

    if (IsNan(m_ZPE))
    {
      //cinfo << "calculation of DOS cutailed because no ZPE" << endl;
    }
    else
    {
      const int cellOffset = get_cellOffset();
      std::vector<double> cellEne;
      getCellEnergies(MaximumCell, m_host->getEnv().CellSize, cellEne);
      calcGrainAverages(m_host->getEnv().MaxGrn, m_host->getEnv().cellPerGrain(), cellOffset, m_cellDOS, cellEne, m_grainDOS, m_grainEne);
    }

    if (recalc) {
      testDensityOfStates();
    }

    recalculateDOScompleted();

    return true;
  }

  // Calculate classical energy
  double gDensityOfStates::getClassicalEnergy(){
    //Basically use the frequencies to calculate the contribution of ZPE from harmonic oscillators approximation
    double ZC = 0.0;
    for (unsigned int i = 0; i < m_VibFreq.size(); ++i)
      ZC += m_VibFreq[i] / 2.0;
    return get_zpe() - ZC;
  }


  //
  // Test the rovibronic density of states
  //
  void gDensityOfStates::testDensityOfStates()
  {
    const int MaximumGrain = m_host->getEnv().MaxGrn;
    const int MaximumCell = m_host->getEnv().MaxCell;
    std::vector<double> cellEne;
    getCellEnergies(MaximumCell, m_host->getEnv().CellSize, cellEne);

    // Partition functions that are higher than the current simulation temperature will not be output.
    const double temperature = 1. / (boltzmann_RCpK * m_host->getEnv().beta);
    const int max_nplus1 = int(temperature / 100.);

    if (m_host->isMolType("modelled") || m_host->isMolType("transitionState")){
      string comment("Rovibronic partition function calculation at various temperatures. qtot : product of QM partition functions for vibrations (1-D harmonic oscillator) and classical partition functions for rotations.  sumc : cell based partition function. sumg : grain based partition function ");

      PersistPtr ppList = m_host->get_PersistentPointer()->XmlWriteMainElement("me:densityOfStatesList", comment);

      if (m_host->getFlags().testDOSEnabled) ctest << endl << "Test rovibronic density of states for: " << m_host->getName() << "\n{\n";
      if (m_host->getFlags().testDOSEnabled) ctest << "      T           qtot           sumc           sumg\n";


      //loop through predefined test temperatures
      for (int n = 0; n < max_nplus1; ++n) {
        double temp = 100.0*static_cast<double>(n + 2);
        double beta = 1.0 / (boltzmann_RCpK*temp);

        // Calculate rovibronic partition functions based on cells.
        double cellCanPrtnFn = canonicalPartitionFunction(m_cellDOS, cellEne, beta);

        // Calculate rovibronic partition functions based on grains.
        double grainCanPrtnFn = canonicalPartitionFunction(m_grainDOS, m_grainEne, beta);

        // Calculate rovibronic partition functions, using analytical formula where possible.
        double qtot(1.0);
        for (vector<DensityOfStatesCalculator*>::size_type j = 0; j < m_DOSCalculators.size(); ++j) {
          qtot *= m_DOSCalculators[j]->canPrtnFnCntrb(this, beta);
        }

        if (m_host->getFlags().testDOSEnabled) {
          formatFloat(ctest, temp, 6, 7);
          formatFloat(ctest, qtot, 6, 15);
          formatFloat(ctest, cellCanPrtnFn, 6, 15);
          formatFloat(ctest, grainCanPrtnFn, 6, 15);
          ctest << endl;
        }

        //Add to XML document
        PersistPtr ppItem = ppList->XmlWriteElement("me:densityOfStates");
        ppItem->XmlWriteValueElement("me:T", temp, 6);
        ppItem->XmlWriteValueElement("me:qtot", qtot, 6);
        ppItem->XmlWriteValueElement("me:sumc", cellCanPrtnFn, 6);
        ppItem->XmlWriteValueElement("me:sumg", grainCanPrtnFn, 6);
      }
      if (m_host->getFlags().testDOSEnabled) ctest << "}" << endl;
    }

    if (m_host->getFlags().cellDOSEnabled){
      ctest << endl << "Cell rovibronic density of states of " << m_host->getName() << endl << "{" << endl;
      for (int i(0); i < MaximumCell; ++i){
        formatFloat(ctest, cellEne[i], 6, 15);
        formatFloat(ctest, m_cellDOS[i], 6, 15);
        ctest << endl;
      }
      ctest << "}" << endl;
    }

    if (m_host->getFlags().grainDOSEnabled && (m_host->isMolType("modelled") || m_host->isMolType("transitionState"))){
      ctest << endl << "Grain rovibronic density of states of " << m_host->getName() << endl << "{" << endl;
      for (int i(0); i < MaximumGrain; ++i){
        formatFloat(ctest, m_grainEne[i], 6, 15);
        formatFloat(ctest, m_grainDOS[i], 6, 15);
        ctest << endl;
      }
      ctest << "}" << endl;
    }
  }

  double gDensityOfStates::get_zpe() {
    if (m_ZPE_chk == -1) {
      //      cinfo << "m_ZPE was not defined but requested in " << m_host->getName() << ". Default value " << m_ZPE << " is given." << endl;
      --m_ZPE_chk;
    }
    else if (m_ZPE_chk < -1) {
      --m_ZPE_chk;
    }
    else {
      ++m_ZPE_chk;
    }
    return double(m_ZPE);
  }

  double gDensityOfStates::get_scaleFactor() {
    if (m_scaleFactor_chk == -1){
      cinfo << "m_scaleFactor was not defined but requested in " << m_host->getName() << ". Default value " << m_scaleFactor << " is given." << endl;
      --m_scaleFactor_chk;
    }
    else if (m_scaleFactor_chk < -1) {
      --m_scaleFactor_chk;
    }
    else {
      ++m_scaleFactor_chk;
    }
    return m_scaleFactor;
  }

  double gDensityOfStates::get_Sym(void){
    if (m_Sym_chk == -1){
      cinfo << "m_Sym was not defined but requested in " << m_host->getName() << ". Default value " << m_Sym << " is given." << endl;
      --m_Sym_chk;
    }
    else if (m_Sym_chk < -1) {
      --m_Sym_chk;
    }
    else {
      ++m_Sym_chk;
    }
    return m_Sym;
  }

  void gDensityOfStates::get_VibFreq(std::vector<double>& vibFreq){
    const double scalefactor = get_scaleFactor();
    for (unsigned int i = 0; i < m_VibFreq.size(); ++i)
      vibFreq.push_back(m_VibFreq[i] * scalefactor);
  }

  bool gDensityOfStates::removeVibFreq(double freq) {
    vector<double>::iterator pos = find(m_VibFreq.begin(), m_VibFreq.end(), freq);
    if (pos == m_VibFreq.end())
      return false;
    m_VibFreq.erase(pos);
    return true;
  }

  int gDensityOfStates::getSpinMultiplicity(){
    if (m_SpinMultiplicity_chk >= 0) {
      ++m_SpinMultiplicity_chk;
    }
    else {
      cinfo << "m_SpinMultiplicity was not defined but requested in " << m_host->getName() << ". Default value " << m_SpinMultiplicity << " is given." << endl;
    }
    return m_SpinMultiplicity;
  }

  int gDensityOfStates::get_cellOffset(void) {
    double modulus = fmod(get_zpe() - m_host->getEnv().EMin, double(m_host->getEnv().GrainSize)) / m_host->getEnv().CellSize;
    return int(max(modulus, 0.0));
  };

  //
  // Get grain density of states.
  //
  void gDensityOfStates::getGrainDensityOfStates(vector<double> &grainDOS, const int startGrnIdx, const int ignoreCellNumber) {
    // If density of states have not already been calcualted then do so.
    if (!calcDensityOfStates()){
      throw (std::runtime_error("Failed calculating DOS."));
    }
    if (ignoreCellNumber == 0){ // If there is no cells ignored in this grain, the grain DOS dose not need to be recalculated.
      grainDOS = m_grainDOS;
    }
    else{ // Some cells are ignored in this grain, as they do not occur in this part of reaction.
      // first deal with the first grain.
      // const int MaximumCell = m_host->getEnv().MaxCell;
      const int gsz = m_host->getEnv().GrainSize;
      const int cellOffset = get_cellOffset();
      const int grnStartCell = startGrnIdx * gsz - cellOffset;
      double partialDOS(0.0);
      for (int i(ignoreCellNumber); i < gsz; ++i){
        partialDOS += m_cellDOS[i + grnStartCell];
      }
      grainDOS.clear();
      grainDOS.push_back(partialDOS);
      for (int i(startGrnIdx + 1); i < int(m_grainDOS.size()); ++i){
        grainDOS.push_back(m_grainDOS[i]);
      }
    }
  }

  //
  // Get grain energies.
  //
  void gDensityOfStates::getGrainEnergies(vector<double> &grainEne) {
    // If density of states have not already been calcualted then do so.
    if (!calcDensityOfStates())
      cerr << "Failed calculating DOS";
    grainEne = m_grainEne;
  }

  //
  // Get Grain canonical partition function for rotational, vibrational, and electronic contributions.
  //
  double gDensityOfStates::rovibronicGrnCanPrtnFn() {

    // If density of states have not already been calculated then do so.
    if (!calcDensityOfStates())
      cerr << "Failed calculating DOS";

    return canonicalPartitionFunction(m_grainDOS, m_grainEne, m_host->getEnv().beta);

  }

  //
  // Calculate standard thermodynamic quantities as a function of temperature.
  // The calculation is based on cell densities of states.
  //
  bool gDensityOfStates::thermodynamicsFunctions(double temp, double unitFctr, double& enthalpy, double& entropy, double& gibbsFreeEnergy) {

    std::vector<double> cellEne;
    getCellEnergies(m_host->getEnv().MaxCell, m_host->getEnv().CellSize, cellEne);

    calcDensityOfStates();

    double beta;
    if (temp > 0.0) {
      beta = 1.0 / (boltzmann_RCpK*temp);
    }
    else {
      return false;
    }

    // Calculate rovibronic partition functions based on cells.
    double cellCanPrtnFn = canonicalPartitionFunction(m_cellDOS, cellEne, beta);

    // The following calculates the mean internal molecular energy.
    double internalEnergy = canonicalMeanEnergy(m_cellDOS, cellEne, beta);

    // The rovibronic partition function must be corrected for translation 
    // and (assuming an ideal gas) molecular indistinguishability.
    double molarVol = 1.e+06*idealGasC / (boltzmann_RCpK*beta*atm_in_pascal);  // cm3
    gibbsFreeEnergy = unitFctr*(-log(cellCanPrtnFn)
      - log(tp_C * pow((m_host->getStruc().getMass() / beta), 1.5)*molarVol)
      + log(AvogadroC)) / beta;

    // The enthalpy must be corrected for translation by an additional 3kT/2.
    enthalpy = unitFctr*(internalEnergy + 5.0 / (2.0*beta));

    entropy = (enthalpy - gibbsFreeEnergy) / temp;

    return true;
  }

  // Calculate vibrational frequencies from molecular Hessian. This method 
  // projects out the overall translation and rotation vectors as defined 
  // in Wilson, Decius and Cross. Molecular Vibrations, Dover 1980. See
  // also Miller, Handy and Adams JCP, Vol, 72, 99 (1980).
  bool gDensityOfStates::FrqsFromHessian() {

    const size_t msize = m_Hessian->size();
    m_Modes = new dMatrix(msize, 0.0);

    gStructure& gs = m_host->getStruc();
    bool HasCoords = gs.ReadStructure();
    if (!HasCoords || 3 * size_t(gs.NumAtoms()) != msize) {
      cerr << "The dimension of the defined Hessian for " << getHost()->getName() << " does not match the specified number of atoms.";
      return false;
    }

    // Get atomic massess and coordinates.

    vector<double> atomicMasses, xx, yy, zz;
    gs.getAtomicMasses(atomicMasses);
    gs.getXCoords(xx);
    gs.getYCoords(yy);
    gs.getZCoords(zz);

    // Initialize mass weights.

    vector<double> massWeights;
    size_t i(0), j(0);
    for (; i < atomicMasses.size(); i++) {
      double weight = sqrt(atomicMasses[i]);
      for (j = 0; j < 3; j++) {
        massWeights.push_back(weight);
      }
    }

    // Mass weight Hessian.

    double convFactor(1.0);
    if (m_HessianUnits == "kJ/mol/Ang2") {
      // Nothing to do.
    }
    else if (m_HessianUnits == "kcal/mol/Ang2") {
      convFactor *= Calorie_in_Joule;
    }
    else if (m_HessianUnits == "Hartree/Bohr2") {
      convFactor *= Hartree_In_kJperMol / (bohr_in_angstrom * bohr_in_angstrom);
    }
    else {
      throw (std::runtime_error("Unknown Hessian units."));
    }

    for (i = 0; i < msize; i++) {
      for (j = i; j < msize; j++) {
        (*m_Hessian)[i][j] *= convFactor / (massWeights[i] * massWeights[j]);
        (*m_Hessian)[j][i] = (*m_Hessian)[i][j];
      }
    }

    // Rotate Hessian.

    dMatrix axisAlignment(3, 0.0);
    gs.getAlignmentMatrix(axisAlignment);
    dMatrix Transform(msize, 0.0);
    for (i = 0; i < msize; i += 3) {
      size_t ii(0), jj(0);
      for (ii = 0; ii < 3; ii++) {
        for (jj = 0; jj < 3; jj++) {
          Transform[i + ii][i + jj] = axisAlignment[ii][jj];
        }
      }
    }

    dMatrix tmpHessian = (*m_Hessian)*Transform;
    Transform.Transpose();
    *m_Hessian = Transform*tmpHessian;

    // X Translation projector.

    vector<double> mode(msize, 0.0);
    for (i = 0; i < msize; i += 3)
      mode[i] = massWeights[i];

    UpdateProjector(mode);

    // Y Translation projector.

    ShiftTransVector(mode);
    UpdateProjector(mode);

    // Z Translation projector.

    ShiftTransVector(mode);
    UpdateProjector(mode);

    // Rotational modes.
    RotationVector(yy, 2, 1.0, zz, 1, -1.0, massWeights, mode);
    UpdateProjector(mode);

    RotationVector(xx, 2, -1.0, zz, 0, 1.0, massWeights, mode);
    UpdateProjector(mode);

    RotationVector(xx, 1, 1.0, yy, 0, -1.0, massWeights, mode);
    UpdateProjector(mode);

    // Project out translational and rotational modes.

    vector<double> freqs(msize, 0.0);
    calculateFreqs(freqs, m_host->isMolType("transitionState"));

    // Save projected frequencies. Note need to check if configuration is a transtion state.

    size_t firstFrqIdx = (m_host->isMolType("transitionState")) ? 7 : 6;
    m_VibFreq.clear();
    for (i = firstFrqIdx; i < msize; i++)
      m_VibFreq.push_back(freqs[i]);

    return true;
  }

  // Function to calculate the vibrational frequencies from a projected Hessian matrix.
  bool gDensityOfStates::calculateFreqs(vector<double> &freqs, bool projectTransStateMode) {

    const size_t msize = m_Hessian->size();

    dMatrix tmp = *m_Modes;
    tmp.Transpose();
    dMatrix Projector = (*m_Modes)*tmp;
    for (size_t i = 0; i < msize; i++) {
      for (size_t j = 0; j < msize; j++) {
        Projector[i][j] *= -1.0;
      }
      Projector[i][i] += 1.0;
    }

    dMatrix tmpHessian = Projector*(*m_Hessian)*Projector;

    tmpHessian.diagonalize(&freqs[0]);

    double convFactor = conHess2Freq / (2.0*M_PI);
    for (size_t m(0); m < msize; m++) {
      if (freqs[m] > 0.0) {

        // Mostly vibration modes.
        freqs[m] = convFactor*sqrt(freqs[m]);

      }
      else if (projectTransStateMode) {

        // Add Transition state mode to the projected set after orthogonalization.
        // The magic number of "1.0" below is 1 cm-1 and used to filter out any 
        // small -ve frequencies associated with existing projected modes.
        double imFreq = convFactor*sqrt(fabs(freqs[m]));
        if (imFreq > 1.0) {

          vector<double> mode(msize, 0.0);
          for (size_t j(0); j < msize; j++) {
            mode[j] = tmpHessian[j][m];
          }

          orthogonalizeMode(mode);
        }
      }
    }

    return true;
  }

  // This method is used to project a mode from the stored Hessian and
  // re-calculate the remaining frequencies.
  bool gDensityOfStates::projectMode(vector<double> &mode) {

    bool status(true);

    const size_t msize = m_Hessian->size();

    gStructure& gs = m_host->getStruc();
    vector<double> atomicMasses;
    gs.getAtomicMasses(atomicMasses);

    size_t i, j, k;
    vector<double> massWeights(msize, 0.0);
    for (j = 0, i = 0; j < atomicMasses.size(); j++) {
      double weight = sqrt(atomicMasses[j]);
      for (k = 0; k < 3; k++, i++) {
        mode[i] *= weight;
      }
    }

    // Orthogonalize and project out mode.

    orthogonalizeMode(mode);

    vector<double> freqs(msize, 0.0);
    calculateFreqs(freqs);

    // Save projected frequencies.

    size_t nfreq = m_VibFreq.size();
    m_VibFreq.clear();
    for (i = msize - nfreq + 1; i < msize; i++)
      m_VibFreq.push_back(freqs[i]);

    return status;
  }

  // This method is used to orthogonalize a mode against existing
  // projected modes and then add it to the projected set.
  bool gDensityOfStates::orthogonalizeMode(vector<double> &mode) {

    const size_t msize = m_Hessian->size();

    // Orthogonalize against existing modes. 

    size_t i(0), j(0);
    for (i = 0; i < m_nModes; i++) {
      double sum(0.0);
      for (j = 0; j < msize; j++) {
        sum += mode[j] * (*m_Modes)[j][i];
      }
      for (j = 0; j < msize; j++) {
        mode[j] -= sum * (*m_Modes)[j][i];
      }
    }

    UpdateProjector(mode);

    return true;
  }

  // Helper function to shift translation projection vector.
  void gDensityOfStates::ShiftTransVector(vector<double> &mode) {
    const size_t msize = mode.size();
    for (size_t i = msize - 1; i > 0; i--)
      mode[i] = mode[i - 1];
    mode[0] = 0.0;
  }

  // Helper function to create projector.
  void gDensityOfStates::UpdateProjector(vector<double> &mode) {

    // Normalize mode.

    double NormFctr(0.0);
    size_t i(0);
    for (; i < mode.size(); i++) {
      NormFctr += mode[i] * mode[i];
    }

    NormFctr = 1.0 / sqrt(NormFctr);
    for (i = 0; i < mode.size(); i++) {
      mode[i] *= NormFctr;
    }

    // Add mode to existing set.

    for (i = 0; i < m_Modes->size(); i++)
      (*m_Modes)[i][m_nModes] = mode[i];
    m_nModes++;

  }

  // Function to calculate the rotational mode vectors.
  void gDensityOfStates::RotationVector(vector<double> &aa, size_t loca, double sgna, vector<double> &bb, size_t locb, double sgnb, vector<double> &massWeights, vector<double> &mode) {

    mode.clear();
    size_t ncoords = aa.size();
    size_t i(0);
    for (; i < ncoords; i++) {
      vector<double> rcross(3, 0.0);
      rcross[loca] = sgna*aa[i];
      rcross[locb] = sgnb*bb[i];
      for (size_t j(0); j < 3; j++) {
        mode.push_back(rcross[j]);
      }
    }

    // Mass weight vector ;
    for (i = 0; i < mode.size(); i++) {
      mode[i] *= massWeights[i];
    }

  }

  // This method tests if a rotor is heavy. It is a helper method
  // used to assess if a QM method will be expensive for calculating
  // the energy levels of an asymmetic top.
  bool gDensityOfStates::IsHeavyTop(size_t n) {

    gStructure& gs = m_host->getStruc();
    vector<double> atomicMasses;
    gs.getAtomicMasses(atomicMasses);

    size_t nHeavyAtoms(0) ;
    for (size_t j(0); j < atomicMasses.size(); j++) {
      if (atomicMasses[j] > atomMass("H")) {
	     nHeavyAtoms++ ;
	  }
    }
    return (nHeavyAtoms > n) ;
  }

}//namespace
