//-------------------------------------------------------------------------------------------
//
// ThermodynamicTable.cpp
//
// Author: Struan Robertson
// Date:   06/Mar/2011
//
// This file contains the declaration and implementation of the plug-in class that calculates
// the thermodynamics tables for all the molecules specified in the input file.
//
//-------------------------------------------------------------------------------------------

#include <functional>
#include "../System.h"
#include "../calcmethod.h"
#include "../gDensityOfStates.h"
#include "../gStructure.h"

namespace mesmer
{
  class ThermodynamicTable : public CalcMethod
  {
  public:

    ThermodynamicTable(const char* id) : m_id(id),
      m_nTemp(20),
      m_TempInterval(50.0),
      m_Unit("kJ/mol"),
      m_makeNasaPoly(true)
    {
      Register();
    }

    virtual ~ThermodynamicTable() {}
    virtual const char* getID()  { return m_id; }

    virtual bool DoesOwnParsing() { return true; }

    // Function to do the work
    virtual bool DoCalculation(System* pSys);

    double SdivR(vector<double>::iterator i, double T)const;

  private:

    // Read any data from XML and store in this instance. 
    bool ReadParameters(PersistPtr ppControl);
    virtual bool ParseData(PersistPtr pp); //preferred method

    string underlineText(const string& text) const;

    string writeTableHeader(const string& unit) const;

    void writeTableEntry(Molecule *pmol, double temp, double unitFctr, string & header) const;

    const char* m_id;

    string WriteNASAPoly(Molecule* pmol, vector<double> coeffs, double TLo, double TMid, double THi);

    // TODO This would be better (with T, U random-access iterators) as
    // U FitPoly(unsigned order, T xstart, T xend, U ystart);
    // so that it is more like STL algorithms and can be used with arrays, vectors, etc.
    vector<double> FitPoly(unsigned order,
      vector<double>::const_iterator xstart,
      vector<double>::const_iterator xend,
      vector<double>::const_iterator ystart)const;

    int m_nTemp;
    double m_TempInterval, m_Tmin, m_Tmax, m_Tmid;
    string m_Unit;
    double m_unitFctr;
    bool m_makeNasaPoly;
  } ;

  ////////////////////////////////////////////////
  //Global instance
  ThermodynamicTable theThermodynamicTable("ThermodynamicTable");
  ///////////////////////////////////////////////

  bool ThermodynamicTable::ParseData(PersistPtr pp)
  {
    ErrorContext("ThermodynamicTable");
    //Read in parameters in child elements of CalcMethod or use defaults.xml
    //called from ParseForPlugin in System::Parse()
    m_Unit = pp->XmlReadValue("units");
    m_Tmin = pp->XmlReadDouble("me:Tmin");
    m_Tmax = pp->XmlReadDouble("me:Tmax");
    m_TempInterval = pp->XmlReadDouble("me:Tstep");
    // If Tmid is non-zero the NASA polynomial has two temperature ranges.
    m_Tmid = pp->XmlReadDouble("me:Tmid");
    m_unitFctr = 1.0 / kJPerMol_in_RC;
    if (m_Unit == "kcal/mol")
      m_unitFctr = 1.0 / kCalPerMol_in_RC;

    System* pSys = getParent();
    MoleculeManager* pMoleculeManager = pSys->getMoleculeManager();
    if (pSys->getReactionManager()->size() == 0)
    {
      //No reactions specified, so read in all the molecules here.
      PersistPtr ppmol = pMoleculeManager->get_PersistPtr();
      while (ppmol = ppmol->XmlMoveTo("molecule"))
      {
        // Get the name of the molecule.
        const char* reftxt = ppmol->XmlReadValue("id");
		const char* roletxt = ppmol->XmlReadValue("role");
        if (reftxt) {
          //Use molType="forThermo" to activate DOS properties.
          //pMoleculeManager->addmol(string(reftxt), string("forThermo"), pSys->getEnv(), pSys->m_Flags);
		  pMoleculeManager->addmol(string(reftxt), string(roletxt), pSys->getEnv(), pSys->m_Flags);
        }
      }

      MesmerEnv& Env = pSys->getEnv();
      Env.GrainSize = 100;
      Env.MaxGrn = 1000;
      Env.MaxCell = Env.GrainSize * Env.MaxGrn;

      // Determine if DOS test information is to appear.
      PersistPtr ppControl = pSys->getPersistPtr()->XmlMoveTo("me:control");
      pSys->m_Flags.testDOSEnabled = ppControl->XmlReadBoolean("me:testDOS");
      if (pSys->m_Flags.testDOSEnabled)
        pSys->getEnv().beta = 1.0 / (boltzmann_RCpK*double(m_nTemp)*m_TempInterval);
    }
    if (MolecularComponent::getEnergyConvention() == "arbitrary")
      m_makeNasaPoly = false;

    if (pMoleculeManager->size() == 0 || MolecularComponent::getEnergyConvention().empty())
    {
      cerr << "No suitable molecules have been specified." << endl;
      return false;
    }
    return true;
  }

  bool ThermodynamicTable::DoCalculation(System* pSys)
  {
    // Make provision for the special case of T = 298.15.
    bool   tempLessThan298;
    double temp289(298.15);

    // Loop over all molecules producing a table for each molecule that
    // has an energy specified.
    // A NASA polynomial will also be produced if 7 or more temperatures
    // have been requested in both the upper and lower polynomials 
    // or, if Tmid=0, the single polynomial.

    int nTemps = static_cast<int>(std::floor((m_Tmax - m_Tmin) / m_TempInterval))+1;
    int nTempsLower = static_cast<int>(std::floor((m_Tmid - m_Tmin) / m_TempInterval)+1);
    bool enoughPoints =  !m_Tmid && nTemps > 6 
         || m_Tmid && nTempsLower> 6 && (nTemps-nTempsLower)>=6;
    if (!enoughPoints)
    {
      cinfo << "Too few data points to fit NASA polynomials." << endl;
      m_makeNasaPoly = false;
    }
    double R = boltzmann_C * AvogadroC;
    if (m_Unit == "kcal/mol") R *= Calorie_in_Joule;

    MoleculeManager* pMoleculeManager = pSys->getMoleculeManager();
    MoleculeManager::constMolIter molItr = pMoleculeManager->begin();
    MoleculeManager::constMolIter molItrEnd = pMoleculeManager->end();
    for (; molItr != molItrEnd; molItr++)
    {
      vector<double> temperature, Hf /* enthalpy of formation / R */;
      Molecule *pmol = molItr->second;

      // Restrict output for molecules without a specified energy.
      double Hf298local = NaN;
      if (m_makeNasaPoly)
      {
        Hf298local = pmol->getDOS().get_Hf298Thermo();
        if (IsNan(Hf298local))
          cinfo << "Restricted thermo output for " << pmol->getName()
          << " because it has no non-arbitrary energy data." << endl;
      }
      PersistPtr pp = pmol->get_PersistentPointer();
      pp = pp->XmlWriteMainElement("me:thermoTable", "", true); //will replace an existing element
      //pp = pp->XmlWriteElement("me:thermoTable");
      pp->XmlWriteAttribute("unitsT", "K");
      pp->XmlWriteAttribute("unitsH", m_Unit);
      pp->XmlWriteAttribute("unitsS", m_Unit.substr(1, m_Unit.length()) + "/K");
      pp->XmlWriteAttribute("unitsG", m_Unit);
      if (!IsNan(Hf298local))
        pp->XmlWriteAttribute("unitsHf", m_Unit);

      double S298;//Always calculated. NOTE kJ/mol/K

      double enthalpy298, dummy;
      pmol->getDOS().thermodynamicsFunctions(298.15, m_unitFctr,
        enthalpy298, S298, dummy);
      tempLessThan298 = true;
	  ctest << "\nthermodynamic data based on qtot begin:    " << pmol->getName() << endl;
	  ctest << "unit:[cal][mol][K]" << endl;
      for (double temp = m_Tmin; temp <= m_Tmax; temp += m_TempInterval)
      {
        double T = temp;
        if (tempLessThan298 && temp > temp289)
        {
          // Special case of T = 289.15
          tempLessThan298 = false;
          T = temp289;
          temp -= m_TempInterval;
        }
        temperature.push_back(T);

        double enthalpy(0.0), entropy(0.0), gibbsFreeEnergy(0.0);
        pmol->getDOS().thermodynamicsFunctions(T, m_unitFctr,
          enthalpy, entropy, gibbsFreeEnergy);
        PersistPtr ppVal = pp->XmlWriteElement("me:thermoValue");
        ppVal->XmlWriteAttribute("T", T, 2, true);
        ppVal->XmlWriteAttribute("H", enthalpy, 4, true);
        ppVal->XmlWriteAttribute("S", entropy*1000, 4, true);
        ppVal->XmlWriteAttribute("G", gibbsFreeEnergy, 4, true);
        if (!IsNan(Hf298local))
        {
          Hf.push_back((enthalpy - enthalpy298 + Hf298local) * 1000 / R); //e.g. J/mol
          ppVal->XmlWriteAttribute("Hf", Hf.back()*R / 1000, 4, true); //back to kJ/mol
        }
      }
	  ctest << "thermodynamic data based on qtot end:    " << pmol->getName() << endl << endl;

      if (m_makeNasaPoly && !IsNan(Hf298local))
      {
        //Fit NASA polynomial to enthalpy data
        //H/R =a6 + T*a1 + T^2*a2/2 + a3*T^3/3 + a4*T^4/4 + a5*T^5/5
        vector<double> fits1, fits2; //a1, a2/2, a3/3, etc
        if (m_Tmid == 0) // single range (duplicated)
        {
          fits1 = FitPoly(6, temperature.begin(), temperature.end(), Hf.begin());
          fits1[2] *= 2; fits1[3] *= 3; fits1[4] *= 4; fits1[5] *= 5;
          fits2 = fits1;
        }
        else //two ranges
        {
          vector<double>::iterator itermid = find(temperature.begin(), temperature.end(), m_Tmid);
          if (itermid == temperature.end())
          {
            cerr << "In NASA polynomial fits the middle temperature"
              "must be one of the specified temperatures." << endl;
            return false;
          }
          int nlowerrange = itermid - temperature.begin();
          fits1 = FitPoly(6, itermid, temperature.end(), Hf.begin() + nlowerrange); //upper range
          fits2 = FitPoly(6, temperature.begin(), itermid+1, Hf.begin()); //lower range
          fits1[2] *= 2; fits1[3] *= 3; fits1[4] *= 4; fits1[5] *= 5;
          fits2[2] *= 2; fits2[3] *= 3; fits2[4] *= 4; fits2[5] *= 5;
        }
        vector<double> coeffs(15);
        copy(fits1.begin() + 1, fits1.end(), coeffs.begin());
        copy(fits2.begin() + 1, fits2.end(), coeffs.begin() + 7);

        coeffs[5] = fits1[0];
        coeffs[12] = fits2[0];
        coeffs[14] = Hf298local/R;

        //Set a14 to match S at 298.15K
        coeffs[13] = 0.0;
        coeffs[13] = S298*1000/R - SdivR(coeffs.begin()+7, 298.15);
 
        //Set a7 to match a) S at 298K for one range; b) S at Tmid for two range;
        if(m_Tmid==0)
          coeffs[6] = coeffs[13];
        else
        {
          double Smid(0.0), dummy1(0.0), dummy2(0.0); //dummys must be diffferent variables!
          pmol->getDOS().thermodynamicsFunctions(m_Tmid, m_unitFctr,
            dummy1, Smid, dummy2);
          coeffs[6] = 0.0;
          coeffs[6] = Smid*1000/R - SdivR(coeffs.begin(), m_Tmid);
        }

        // Output to XML using a CML property for Nasa Polynomials
        // previously used in OpenBabel.
        PersistPtr pp = pmol->get_PersistentPointer();
        PersistPtr ppProp = pp->XmlMoveTo("propertyList");
        //ppProp = (ppProp ? ppProp : pp)->XmlWriteMainElement("property","",true);
        ppProp = (ppProp ? ppProp : pp)->XmlWriteElement("property");
        ppProp->XmlWriteAttribute("dictRef", "NasaPolynomial");

        stringstream ss;
        ss << temperature[0];
        PersistPtr ppScalar = ppProp->XmlWriteValueElement("scalar", ss.str());
        ppScalar->XmlWriteAttribute("dictRef", "NasaLowT");

        ss.str("");
        ss << temperature.back();
        ppScalar = ppProp->XmlWriteValueElement("scalar", ss.str());
        ppScalar->XmlWriteAttribute("dictRef", "NasaHighT");

        ss.str("");
        ss << m_Tmid ? m_Tmid : temperature.back();
        ppScalar = ppProp->XmlWriteValueElement("scalar", ss.str());
        ppScalar->XmlWriteAttribute("dictRef", "NasaMidT");

        ppScalar = ppProp->XmlWriteValueElement("scalar", "G");
        ppScalar->XmlWriteAttribute("dictRef", "Phase");

        stringstream vals;
        std::copy(coeffs.begin(), coeffs.end(), ostream_iterator<double>(vals, " "));
        ppScalar = ppProp->XmlWriteValueElement("array", vals.str());
        ppScalar->XmlWriteAttribute("dictRef", "NasaCoeffs");
        ppScalar->XmlWriteAttribute("size", "15");

        string poly = WriteNASAPoly(pmol, coeffs, temperature[0],
          m_Tmid ? m_Tmid : temperature.back(), temperature.back());
        ppScalar = ppProp->XmlWriteValueElement(
          "scalar", poly, true); //Output polynomial as CDATA
        ppScalar->XmlWriteAttribute("dictRef", "NasaPolynomial");
        ctest << poly << endl; // for QA test
      }
    }

    ////TEST FitPoly
    //vector<double> xdata = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 };
    //vector<double> ydata = { 1, 6, 17, 34, 57, 86, 121, 162, 209, 262, 321 };
    ////Solution is 3 x2 + 2 x + 1 //ok

    //vector<double> xdata = { 1000, 1050, 1100, 1150, 1200, 1250, 1300, 1350, 1400, 2450, 1500 };
    //vector<double> ydata = { 22.69882061, 24.44133150, 26.19332847, 27.95387862, 29.72215460,
    //  31.49742271, 33.27903189, 35.06640374, 36.85902371, 38.65643325, 40.45822296 };
    //vector<double> result = FitPoly(6, xdata.begin(), xdata.end(), ydata.begin());
    //double t = xdata[10];
    //double val(0);
    //for (int i = result.size()-1; i >= 0; --i)
    //  val = (val*t + result[i]); //ok

     return true;
  }

  //Return coefficients of x in a polynomial x^0 to x^order calculated
  //using a least squares fit data points (xdata,ydata).
  //Matrix elements from http://www.codecogs.com/library/maths/approximation/regression/discrete.php
  vector<double> ThermodynamicTable::FitPoly(unsigned order,
    vector<double>::const_iterator xstart,
    vector<double>::const_iterator xend,
    vector<double>::const_iterator ystart)const
  {
    TMatrix<double> matrix(order);
    unsigned n = xend - xstart;
    vector<double> rhs(order);
    if (!(order!=0 && n!=0)) //size of ystart not checked
      return rhs; //empty on error
    double sum;
    for (unsigned ir = 0; ir != order; ++ir) //each row
    {
      for (unsigned ic = 0; ic != order; ++ic) //each column
      {
        sum = 0.0;
        for (unsigned j = 0; j != n; ++j)
          sum += pow(*(xstart+j), int(ir+ic));
        matrix[ir][ic] = sum;
      }
      sum=0.0;
      for (unsigned j = 0; j != n; ++j)
        sum += pow(*(xstart+j), int(ir)) * *(ystart+j);
      rhs[ir] = sum;
    }
    matrix.solveLinearEquationSet(&rhs[0]);
    return rhs;
  }

  double ThermodynamicTable::SdivR(vector<double>::iterator i, double T) const
  {
    //S/R  = a1 lnT + a2 T + a3 T^2 /2 + a4 T^3 /3 + a5 T^4 /4 + a7
    //return *i*log(T) + *(i+1)*T + *(i+2)*T*T / 2 + *(i+3)*T*T*T / 3 + *(i+4)*T*T*T*T / 4 + *(i+6);
    return *i*log(T) +T*(*(i+1) + T*(*(i+2)/2 + T*(*(i+3)/3 + T*(*(i+4)/4)))) + *(i+6);
  }

  string ThermodynamicTable::WriteNASAPoly(Molecule* pmol, vector<double> coeffs,
                                           double TLo, double TMid, double THi)
  {
    stringstream ss;
    unsigned int i;
#ifdef _MSC_VER
    unsigned oldf = _set_output_format(_TWO_DIGIT_EXPONENT);
#endif
    ss << '\n';
    ss << left << setw(24) << pmol->getName().substr(0, 24);
    map<string, int> Comp = pmol->getStruc().GetElementalComposition();
    int npad = 4 - Comp.size();
    map<string, int>::const_iterator itr = Comp.begin() ;
    for (;  itr != Comp.end() ; itr++ )
      ss << left << setw(2) << itr->first << right << setw(3) << itr->second ;
    for (; npad; --npad)
      ss << "     ";
    ss << right << 'G' << fixed << setprecision(3) << setw(10) << TLo;
    ss << setw(10) << THi << setw(9) << TMid << "    01" << '\n';

    ss << scientific << setprecision(7);
    for (i = 0; i<5; ++i)
      ss << setw(15) << coeffs[i];
    ss << "    2\n";
    for (i = 5; i<10; ++i)
      ss << setw(15) << coeffs[i];
    ss << "    3\n";
    for (i = 10; i<15; ++i)
      ss << setw(15) << coeffs[i];
    ss << "    4" << endl;

#ifdef _MSC_VER
    _set_output_format(oldf);
#endif

    return ss.str();
  }
}//namespace

