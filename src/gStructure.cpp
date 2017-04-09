#include "Molecule.h"
#include "System.h"
#include "vector3.h"
#include "gDensityOfStates.h"
#include "gStructure.h"

using namespace OpenBabel;

namespace mesmer
{
  gStructure::gStructure(mesmer::Molecule *pMol) : m_MolecularWeight(-1),
    m_PrincipalMI(3, 0.0),
    m_AxisAlignment(NULL),
    Atoms(),
    Bonds(),
    m_atomicOrder(),
    m_HasCoords(false),
    m_verbose(false),
    m_RotBondIDs()
  {
    ErrorContext c(pMol->getName());
    m_host = pMol;
    PersistPtr pp = pMol->get_PersistentPointer();
    PersistPtr ppPropList = pp->XmlMoveTo("propertyList");
    if (!ppPropList)
      ppPropList = pp; //Be forgiving; we can get by without a propertyList element
    double MW = ppPropList->XmlReadPropertyDouble("me:MW", optional);
    ReadStructure();
	if (IsNan(MW))
    {
      ReadStructure();
      if (Atoms.empty())
        cinfo << "Neither chemical structure nor "
        "Molecular Weight as an XML property were provided." << endl;
      else
        MW = CalcMW();
    }
    setMass(MW);
  }

  map<string, int> gStructure::GetElementalComposition() const
  {
    map<string, int> Composition;
    if (!Atoms.empty())
      for (map<string, atom>::const_iterator it = Atoms.begin(); it != Atoms.end(); ++it)
        ++Composition[it->second.element];
    return Composition;
  }

  // Provide a function to define particular counts of the convolved DOS of two molecules.
  bool countDimerCellDOS
    (gDensityOfStates& pDOS1, gDensityOfStates& pDOS2, std::vector<double>& rctsCellDOS){
      std::vector<double> rct1CellDOS;
      std::vector<double> rct2CellDOS;
      if (!pDOS1.getCellDensityOfStates(rct1CellDOS) || !pDOS2.getCellDensityOfStates(rct2CellDOS))
        return false;
      std::vector<double> rotConsts;
      if (pDOS1.get_rotConsts(rotConsts) == UNDEFINED_TOP){
        rctsCellDOS = rct2CellDOS;
      }
      else if (pDOS2.get_rotConsts(rotConsts) == UNDEFINED_TOP){
        rctsCellDOS = rct1CellDOS;
      }
      else{
        FastLaplaceConvolution(rct1CellDOS, rct2CellDOS, rctsCellDOS);
      }
      return true;
  }

  //Returns true if atoms have coordinates
  bool gStructure::ReadStructure()
  {
    if (!Atoms.empty())
      return m_HasCoords;
    PersistPtr ppMol = getHost()->get_PersistentPointer();
    PersistPtr ppAtom = ppMol->XmlMoveTo("atomArray");
    if (!ppAtom) // there may not be an <atomArray> element
      ppAtom = ppMol;
    while (ppAtom = ppAtom->XmlMoveTo("atom"))
    {
      atom at;
      const char* el = ppAtom->XmlReadValue("elementType");
      if (!el)
      {
        cerr << "<atom> elements must have an elementType attribute" << endl;
        return false;
      }
      at.element = el;
      const char* pId = ppAtom->XmlReadValue("id", optional);
      if (pId)
        at.id = pId;
      else
      {
        //Atom has no id so give it one.
        stringstream ss;
        ss << "a" << (Atoms.size()+1) ;
        at.id = ss.str();
        ppAtom->XmlWriteAttribute("id", at.id);
      }
      double x3, y3, z3;
      x3 = ppAtom->XmlReadDouble("x3", optional);
      y3 = ppAtom->XmlReadDouble("y3", optional);
      z3 = ppAtom->XmlReadDouble("z3", optional);
      if (!IsNan(x3) && !IsNan(y3) && !IsNan(z3))
      {
        at.coords.Set(x3, y3, z3);
        if (x3 != 0 || y3 != 0 || z3 != 0)
          m_HasCoords = true; //at least one atom with non-zero coordinates
      }
      Atoms[at.id] = at;
      m_atomicOrder.push_back(at.id);
    }

    // If coordinates are defined, shift coordinates to the centre of mass/principal axis frame.
    if (m_HasCoords)
      AlignCoords();

    //Read all the bonds. For each bond add a connect to each atom
    PersistPtr ppBond = ppMol->XmlMoveTo("bondArray");
    int ibond = 1;
    if (!ppBond) // there may not be an <bondArray> element
      ppBond = ppMol;
    while (ppBond = ppBond->XmlMoveTo("bond"))
    {
      const char* pId = ppBond->XmlReadValue("id", optional);
      string id;
      if (pId)
        id = pId;
      else
      {
        //id is e.g. "bond3", if not provided
        stringstream ss;
        ss << " bond" << ibond;
        id = ss.str();
      }

      const char* pRefs = ppBond->XmlReadValue("atomRefs2");
      if (!pRefs) return false;
      string refs(pRefs);
      string::size_type pos = refs.find_first_of(" ,");
      string::size_type pos2 = refs.find_last_of(" ,");
      if (pos == string::npos) return false;
      string atomref1 = refs.substr(0, pos);
      string atomref2 = refs.substr(pos2 + 1);
      Bonds[id] = make_pair(atomref1, atomref2);
      Atoms[atomref1].connects.push_back(atomref2);
      Atoms[atomref2].connects.push_back(atomref1);
      ++ibond;
    }

    return m_HasCoords;
  }

  // Method to shift coordinates to the centre of mass/principal axis frame. 
  bool gStructure::AlignCoords() {

    // Calculate centre of mass and subtract from coordinates, and calculate mass weights.

    bool status(true);

    //Determine centre of mass
    map<string, atom>::iterator iter;
    vector3 centreOfMass;
    double mt = 0.0;
    for (iter = Atoms.begin(); iter != Atoms.end(); ++iter) {
      double mass = atomMass(iter->second.element);
      centreOfMass += iter->second.coords * mass;
      mt += mass;
    }
    centreOfMass /= mt;

    dMatrix MI(3);
    double sxx = 0.0, syy = 0.0, szz = 0.0, sxy = 0.0, sxz = 0.0, syz = 0.0;
    for (iter = Atoms.begin(); iter != Atoms.end(); ++iter)
    {
      vector3 c = iter->second.coords - centreOfMass;
      iter->second.coords.Set(c.x(), c.y(), c.z());
      double  m = atomMass(iter->second.element);
      sxx += m * c.x() * c.x();
      syy += m * c.y() * c.y();
      szz += m * c.z() * c.z();
      sxy += m * c.x() * c.y();
      sxz += m * c.x() * c.z();
      syz += m * c.y() * c.z();
    }

    if (NumAtoms() == 2)
      m_PrincipalMI[0] = szz;
    else {
      MI[0][0] = syy + szz;
      MI[1][1] = sxx + szz;
      MI[2][2] = sxx + syy;
      MI[0][1] = MI[1][0] = -sxy;
      MI[0][2] = MI[2][0] = -sxz;
      MI[1][2] = MI[2][1] = -syz;

      MI.diagonalize(&m_PrincipalMI[0]);

      // Save the Alignment matrix 

      m_AxisAlignment = new dMatrix(MI);

      // Rotate coordinates to principal axis frame.

      MI.Transpose();
      for (iter = Atoms.begin(); iter != Atoms.end(); ++iter) {
        vector<double> r(3, 0.0);
        iter->second.coords.Get(&r[0]);
        r *= MI;
        iter->second.coords.Set(&r[0]);
      }
    }

    return status;
  }

  double gStructure::CalcMW()
  {
    map<string, atom>::iterator iter;
    double MW = 0.0;
    for (iter = Atoms.begin(); iter != Atoms.end(); ++iter)
      MW += atomMass(iter->second.element);
    return MW;
  }

  // Returns in atomset the IDs of all the atoms attached to atomID via bonds,
  // but does not include any atoms already in atomset or atoms beyond them.
  // Handles rings. (Recursive function) 
  void gStructure::GetAttachedAtoms(vector<string>& atomset, const string& atomID)
  {
    atomset.push_back(atomID);
    vector<string>::iterator coniter;
    for (coniter = Atoms[atomID].connects.begin(); coniter != Atoms[atomID].connects.end(); ++coniter)
    {
      if (find(atomset.begin(), atomset.end(), *coniter) != atomset.end())
        continue;
      GetAttachedAtoms(atomset, *coniter);
    }
  }

  double gStructure::CalcMomentAboutAxis(vector<string> atomset, vector3 at1, vector3 at2)
  {
    double sumMoment = 0.0;
    vector<string>::iterator iter;
    for (iter = atomset.begin(); iter != atomset.end(); ++iter)
    {
      vector3 a = Atoms[*iter].coords;
      double d = Point2Line(a, at1, at2);
      sumMoment += atomMass(Atoms[*iter].element) * d * d;
    }
    if (sumMoment > 0.0){
      return sumMoment;
    }
    else {
      string errorMsg = "Problem with calculation of reduced moment of inertia. Either there is something wrong with the atom co-ordinates or the atomic symbol is not present in unitsConversion.cpp ";
      throw (std::runtime_error(errorMsg));
    }
  }

  // Calculates internal rotation eigenvector about an axis define by at1 and at2.
  bool gStructure::CalcInternalRotVec(vector<string> atomset, vector3 at1, vector3 at2, vector<double> &mode, bool ApplyMWeight)
  {
    vector3 diff = at1 - at2;
    diff.normalize();
    vector<string>::iterator iter;
    for (iter = atomset.begin(); iter != atomset.end(); ++iter)
    {
      double massWeight = (ApplyMWeight) ? sqrt(atomMass(Atoms[*iter].element)) : 1.0;
      vector3 a = Atoms[*iter].coords - at1;
      vector3 b = cross(a, diff);
      int atomicOrder = getAtomicOrder(*iter);
      if (atomicOrder >= 0) {
        size_t location = 3 * size_t(atomicOrder);
        for (size_t i(location), j(0); j < 3; i++, j++) {
          mode[i] = massWeight*b[j];
        }
      }
      else {
        string errorMsg = "Problem with calculation of internal rotation eigenvector. Atomic order is not correactly defined.";
        throw (std::runtime_error(errorMsg));
      }
    }

    return true;
  }

  // For a bond between atom at1 and atom at2 find all the atoms connected to at1
  // excluding those connect via at2.
  void gStructure::findRotorConnectedAtoms(vector<string> &atomset, const string at1, const string at2) {
    atomset.push_back(at2); // Will not look beyond this atom on the other side of the bond.
    GetAttachedAtoms(atomset, at1);
    atomset.erase(atomset.begin()); // The other side of the bond is not in this set.
  }

  // Calculate the reduce moment of inertia about axis defined by specifed atoms.
  double gStructure::reducedMomentInertia(pair<string, string>& bondats) {

    vector3 coords1 = GetAtomCoords(bondats.first);
    vector3 coords2 = GetAtomCoords(bondats.second);

    // Calculate moment of inertia about bond axis of atoms on one side of bond...
    vector<string> atomset;
    findRotorConnectedAtoms(atomset, bondats.first, bondats.second);
    double mm1 = CalcMomentAboutAxis(atomset, coords1, coords2);

    //...and the other side of the bond
    atomset.clear();
    findRotorConnectedAtoms(atomset, bondats.second, bondats.first);
    double mm2 = CalcMomentAboutAxis(atomset, coords1, coords2);

    /*
    Is the reduced moment of inertia needed about the bond axis or, separately for the set of
    atoms on each side of the bond, about a parallel axis through their centre of mass?
    See:
    http://www.ccl.net/chemistry/resources/messages/2001/03/21.005-dir/index.html
    http://www.ccl.net/chemistry/resources/messages/2001/03/31.002-dir/index.html
    The bond axis is used here.
    */

    return mm1 * mm2 / (mm1 + mm2); //units a.u.*Angstrom*Angstrom

  }

  // Calculate the reduce moment of inertia about axis defined by specifed atoms.
  void gStructure::reducedMomentInertiaAngular(string bondID, double phase, vector<double>& angles,
    vector<double>& redInvMOI, PersistPtr ppConfigData) {

      // Save coordinates.
      exportToXYZ("orig_coords", false, ppConfigData) ;

      map<string, atom>::iterator iter;
      vector<vector3> coordinates;
      for (iter = Atoms.begin(); iter != Atoms.end(); ++iter) {
        coordinates.push_back(iter->second.coords);
      }

  	//  printf("orig_coords:\n");
	  //for(int i=0 ; i<m_atomicOrder.size(); i++)
	  //{
		 // printf("%s\t%f\t%f\t%f\n", Atoms[m_atomicOrder[i]].id.c_str(), Atoms[m_atomicOrder[i]].coords.x(), Atoms[m_atomicOrder[i]].coords.y(), Atoms[m_atomicOrder[i]].coords.z());
	  //}

      pair<string, string> bondats = GetAtomsOfBond(bondID);
      atom &at1 = Atoms[bondats.first];
      atom &at2 = Atoms[bondats.second];

      // Move fragment so that at1 is at origin.

      vector3 origin = at1.coords;
      for (iter = Atoms.begin(); iter != Atoms.end(); ++iter) {
        iter->second.coords -= origin;
      }

      exportToXYZ("origin_at_at1") ;

      // Rotate molecule so that the at1-at2 bond is along z -axis.

      vector3 bond = at2.coords;
      dMatrix rotZ(3, 0.0);
      double radius = sqrt(bond.x()*bond.x() + bond.y()*bond.y());
      double sgn = (bond.x()*bond.y() > 0.0) ? 1.0 : -1.0 ;
      double cosPhi(0.0);
      double sinPhi(0.0);
      if (radius > 0.0) {
        cosPhi =     fabs(bond.x() / radius) ;
        sinPhi = sgn*fabs(bond.y() / radius) ;
      }
      rotZ[0][0] = rotZ[1][1] = cosPhi;
      rotZ[2][2] = 1.0;
      rotZ[0][1] = sinPhi;
      rotZ[1][0] = -rotZ[0][1];

      for (iter = Atoms.begin(); iter != Atoms.end(); ++iter) {
        vector<double> r(3, 0.0);
        iter->second.coords.Get(&r[0]);
        r *= rotZ;
        iter->second.coords.Set(&r[0]);
      }
	  bond = at2.coords;

      exportToXYZ("at2_inXZ");

  //	  printf("before along z-axis:\n");
  //    for (iter = Atoms.begin(); iter != Atoms.end(); ++iter) {
		//vector<double> r(3, 0.0);
  //      iter->second.coords.Get(&r[0]);
		//cout<< iter->second.element <<"\t";
		//vector<double>::iterator it;
		//for(it=r.begin(); it!=r.end();it++)
		//{
		//	printf("%f\t",*it);
		//}
		//printf("\n");
  //    }

      dMatrix rotY(3, 0.0);
      sgn = (bond.z()*bond.x() > 0.0) ? 1.0 : -1.0 ;
//      double cosTheta =     fabs(bond.z() / bond.length()) ;			
//      double sinTheta = sgn*fabs(bond.x() / bond.length());
      double cosTheta = bond.z() / bond.length();			
      double sinTheta = bond.x() / bond.length();
	  rotY[0][0] = rotY[2][2] = cosTheta;
      rotY[1][1] = 1.0;
      rotY[0][2] = -sinTheta ;
      rotY[2][0] = sinTheta ;

      for (iter = Atoms.begin(); iter != Atoms.end(); ++iter) {
        vector<double> r(3, 0.0);
        iter->second.coords.Get(&r[0]);
        r *= rotY;
        iter->second.coords.Set(&r[0]);
      }

      exportToXYZ("z_axis=at1-at2") ;

	 // printf("before rotation:\n");
  //    for (iter = Atoms.begin(); iter != Atoms.end(); ++iter) {
		//vector<double> r(3, 0.0);
  //      iter->second.coords.Get(&r[0]);
		//cout<< iter->second.element <<"\t";
		//vector<double>::iterator it;
		//for(it=r.begin(); it!=r.end();it++)
		//{
		//	printf("%f\t",*it);
		//}
		//printf("\n");
  //    }

      // Determine the content of one of the fragments so that it can moved relative to the other.

      vector<string> atomset;
      findRotorConnectedAtoms(atomset, bondats.first, bondats.second);

      // Rotate fragments so that they are in phase with the potential.

      const double Angle(fmod(phase, 360.0)*M_PI / 180.);
      dMatrix rot(3, 0.0);
      rot[0][0] = rot[1][1] = cos(Angle);
      rot[2][2] = 1.0;
      rot[0][1] = sin(Angle);
      rot[1][0] = -rot[0][1];

      for (size_t j(0); j < atomset.size(); j++) {
        atom &at = Atoms[atomset[j]];
        vector<double> r(3, 0.0);
        at.coords.Get(&r[0]);
        r *= rot;
        at.coords.Set(&r[0]);
      }

      // Rotate one fragment relative to the other.

      redInvMOI.clear();
      redInvMOI.push_back(conMntInt2RotCnt*getGRIT(bondID));
      const size_t nAngle(angles.size());
      const double dAngle = 2.0*M_PI / double(nAngle);
      angles[0] = 0.0;
      rot.reset(rot.size());
      rot[0][0] = rot[1][1] = cos(dAngle);
      rot[2][2] = 1.0;
      rot[0][1] = sin(dAngle);
      rot[1][0] = -rot[0][1];

      for (size_t i(1); i < nAngle; i++) {

        angles[i] = angles[i - 1] + dAngle;

        for (size_t j(0); j < atomset.size(); j++) {
          atom &at = Atoms[atomset[j]];
          vector<double> r(3, 0.0);
          at.coords.Get(&r[0]);
          r *= rot;
          at.coords.Set(&r[0]);
        }

        exportToXYZ(("rotate_fragment_" + toString(i)).c_str()) ;
        redInvMOI.push_back(conMntInt2RotCnt*getGRIT(bondID)) ;
  	  
		//printf("after rotation:\n");
		//for (iter = Atoms.begin(); iter != Atoms.end(); ++iter) {
		//	vector<double> r(3, 0.0);
	 //       iter->second.coords.Get(&r[0]);
		//	cout<< iter->second.element <<"\t";
		//	vector<double>::iterator it;
		//	for(it=r.begin(); it!=r.end();it++)
		//	{
		//		printf("%f\t",*it);
		//	}
		//	printf("\n");
	 //     }

      }

      // Restore original coordinates.

      size_t i(0);
      for (iter = Atoms.begin(); iter != Atoms.end(); ++iter, i++) {
        iter->second.coords = coordinates[i];
      }

      exportToXYZ("orig_coords", true) ;
  }

  // Calculate the internal rotation eigenvector. Based on the internal rotation 
  // mode vector as defined by Sharma, Raman and Green, J. Phys. Chem. (2010).
  // Typically this vector is used to project out an internal rotational mode
  // from a Hessian.
  void gStructure::internalRotationVector(string bondID, vector<double>& mode, bool ApplyMWeight) {

    pair<string, string> bondats = GetAtomsOfBond(bondID);

    vector3 coords1 = GetAtomCoords(bondats.first);
    vector3 coords2 = GetAtomCoords(bondats.second);

    // Calculate moment of inertia about bond axis of atoms on one side of bond...
    vector<string> atomset1;
    findRotorConnectedAtoms(atomset1, bondats.first, bondats.second);
    double mm1 = CalcMomentAboutAxis(atomset1, coords1, coords2);
    CalcInternalRotVec(atomset1, coords1, coords2, mode, ApplyMWeight);

    //...and the other side of the bond
    vector<string> atomset2;
    findRotorConnectedAtoms(atomset2, bondats.second, bondats.first);
    double mm2 = CalcMomentAboutAxis(atomset2, coords1, coords2);
    CalcInternalRotVec(atomset2, coords2, coords1, mode, ApplyMWeight);

    // In the following weights are applied to relative rotation of each 
    // fragment with respect to each other. In the limit that one fragment
    // is infinitely massive the rotation will be confined to the other 
    // fragment. 

    if (!ApplyMWeight) {
      double fctr1 = mm2 / (mm1 + mm2);
      ApplyInertiaWeighting(atomset1, mode, fctr1);

      double fctr2 = mm1 / (mm1 + mm2);
      ApplyInertiaWeighting(atomset2, mode, fctr2);
    }

  }

  // Calculates the GRIT for the current set of coodinates.
  double gStructure::getGRIT(string bondID) {

    size_t msize(m_RotBondIDs.size() + 3);
	////--------------------test begin-----------------------
	//size_t msize(1 + 3);
	////-------------------- test end -----------------------
    dMatrix GRIT(msize, 0.0);

    // Save coordinates.

    map<string, atom>::iterator iter;
    vector<vector3> coordinates;
    for (iter = Atoms.begin(); iter != Atoms.end(); ++iter) {
      coordinates.push_back(iter->second.coords);
    }

    // Determine centre of mass.

    vector3 centreOfMass;
    double sm(0.0);
    for (iter = Atoms.begin(); iter != Atoms.end(); ++iter) {
      double mass = atomMass(iter->second.element);
      centreOfMass += iter->second.coords * mass;
      sm += mass;
    }
    centreOfMass /= sm;

    // Determine ther inertial tensor for overall rotation.

    double sxx = 0.0, syy = 0.0, szz = 0.0, sxy = 0.0, sxz = 0.0, syz = 0.0;
    for (iter = Atoms.begin(); iter != Atoms.end(); ++iter) {
      vector3 c = iter->second.coords - centreOfMass;
      iter->second.coords = c;
      double  m = atomMass(iter->second.element);
      sxx += m * c.x() * c.x();
      syy += m * c.y() * c.y();
      szz += m * c.z() * c.z();
      sxy += m * c.x() * c.y();
      sxz += m * c.x() * c.z();
      syz += m * c.y() * c.z();
    }

    GRIT[0][0] = syy + szz;
    GRIT[1][1] = sxx + szz;
    GRIT[2][2] = sxx + syy;
    GRIT[0][1] = GRIT[1][0] = -sxy;
    GRIT[0][2] = GRIT[2][0] = -sxz;
    GRIT[1][2] = GRIT[2][1] = -syz;

	////--------------------test begin-----------------------
 //   // Calculate the velocity vectors (based on the Sharma, Raman and Green vector).

 //   vector<vector<double> > velocities(m_RotBondIDs.size(), vector<double>(3 * NumAtoms(), 0.0));
 //   for (size_t i(0); i < m_RotBondIDs.size(); i++) {

 //     vector<double> velocity(3 * NumAtoms(), 0.0);

 //     internalRotationVector(m_RotBondIDs[i], velocity, false);

 //     // Remove centre of mass velocity.
 //     vector3 centreOfMassVelocity;
 //     for (size_t j(0); j < m_atomicOrder.size(); j++){
 //       double mass = atomMass((Atoms.find(m_atomicOrder[j]))->second.element);
 //       vector3 vtmp;
 //       vtmp.Set(&velocity[3 * j]);
 //       centreOfMassVelocity += mass*vtmp;
 //     }
 //     centreOfMassVelocity /= sm;

 //     for (size_t j(0), idx(0); j < m_atomicOrder.size(); j++){
 //       for (size_t n(0); n < 3; idx++, n++) {
 //         velocity[idx] -= centreOfMassVelocity[n];
 //       }
 //     }

 //     velocities[i] = velocity;
 //   }

	////get the index of bondID
 //   size_t index_bond(3);
 //   for (index_bond = 0; index_bond < m_RotBondIDs.size() && m_RotBondIDs[index_bond] != bondID; index_bond++);

 //   // Calculate the internal kinetic energy terms.
	//vector<double> &vi = velocities[index_bond];
	//double smk(0.0);
	//for (size_t m(0), idx(0); m < m_atomicOrder.size(); m++)
	//{
	//	double mass = atomMass((Atoms.find(m_atomicOrder[m]))->second.element);
	//	for (size_t l(0); l < 3; l++, idx++) 
	//	{
	//		smk += mass*vi[idx] * vi[idx];
	//	}
	//}
	//GRIT[3][3] = smk;

 //   // Calculate the Coriolis terms.
 //   vector3 coriolis;
 //   for (iter = Atoms.begin(); iter != Atoms.end(); ++iter) 
	//{
	//	size_t ll = 3 * getAtomicOrder(iter->first);
	//	vector3 r = iter->second.coords;
	//	vector3 vtmp;
	//	vtmp.Set(&vi[ll]);
	//	double mass = atomMass(iter->second.element);
	//	coriolis += mass*cross(r, vtmp);
 //   }
 //   GRIT[0][3] = GRIT[3][0] = coriolis.x();
 //   GRIT[1][3] = GRIT[3][1] = coriolis.y();
 //   GRIT[2][3] = GRIT[3][2] = coriolis.z();


 //   // Restore original coordinates.

 //   size_t i(0);
 //   for (iter = Atoms.begin(); iter != Atoms.end(); ++iter, i++) {
 //     iter->second.coords = coordinates[i];
 //   }

 //   if (m_verbose) {
 //     string MatrixTitle("Generalized rotation inertia tensor:");
 //     GRIT.print(MatrixTitle, ctest);
 //     ctest << endl;
 //   }

 //   dMatrix invGRIT(GRIT);
 //   invGRIT.invertLUdecomposition();
 //   if (m_verbose) {
 //     string MatrixTitle = "Inverse of Generalized rotation inertia tensor:";
 //     invGRIT.print(MatrixTitle, ctest);
 //   }

 //   return invGRIT[3][3];
	////--------------------test end-----------------


    // Calculate the velocity vectors (based on the Sharma, Raman and Green vector).

    vector<vector<double> > velocities(m_RotBondIDs.size(), vector<double>(3 * NumAtoms(), 0.0));
    for (size_t i(0); i < m_RotBondIDs.size(); i++) {

      vector<double> velocity(3 * NumAtoms(), 0.0);

      internalRotationVector(m_RotBondIDs[i], velocity, false);

      // Remove centre of mass velocity.
      vector3 centreOfMassVelocity;
      for (size_t j(0); j < m_atomicOrder.size(); j++){
        double mass = atomMass((Atoms.find(m_atomicOrder[j]))->second.element);
        vector3 vtmp;
        vtmp.Set(&velocity[3 * j]);
        centreOfMassVelocity += mass*vtmp;
      }
      centreOfMassVelocity /= sm;

      for (size_t j(0), idx(0); j < m_atomicOrder.size(); j++){
        for (size_t n(0); n < 3; idx++, n++) {
          velocity[idx] -= centreOfMassVelocity[n];
        }
      }

      velocities[i] = velocity;
    }

    // Calculate the internal kinetic energy terms.

    for (size_t i(0), ii(3); i < m_RotBondIDs.size(); i++, ii++) {
      vector<double> &vi = velocities[i];
      for (size_t j(i), jj(ii); j < m_RotBondIDs.size(); j++, jj++) {
        vector<double> &vj = velocities[j];
        double smk(0.0);
        for (size_t m(0), idx(0); m < m_atomicOrder.size(); m++){
          double mass = atomMass((Atoms.find(m_atomicOrder[m]))->second.element);
          for (size_t l(0); l < 3; l++, idx++) {
            smk += mass*vi[idx] * vj[idx];
          }
        }
        GRIT[ii][jj] = smk;
        if (ii != jj)
          GRIT[jj][ii] = GRIT[ii][jj];
      }
    }

    // Calculate the Coriolis terms.

    for (size_t i(0), ii(3); i < m_RotBondIDs.size(); i++, ii++) {
      vector<double> &vi = velocities[i];
      vector3 coriolis;
      for (iter = Atoms.begin(); iter != Atoms.end(); ++iter) {
        size_t ll = 3 * getAtomicOrder(iter->first);
        vector3 r = iter->second.coords;
        vector3 vtmp;
        vtmp.Set(&vi[ll]);
        double mass = atomMass(iter->second.element);
        coriolis += mass*cross(r, vtmp);
      }
      GRIT[0][ii] = GRIT[ii][0] = coriolis.x();
      GRIT[1][ii] = GRIT[ii][1] = coriolis.y();
      GRIT[2][ii] = GRIT[ii][2] = coriolis.z();
    }

    // Restore original coordinates.

    size_t i(0);
    for (iter = Atoms.begin(); iter != Atoms.end(); ++iter, i++) {
      iter->second.coords = coordinates[i];
    }

    if (m_verbose) {
      string MatrixTitle("Generalized rotation inertia tensor:");
      GRIT.print(MatrixTitle, ctest);
      ctest << endl;
    }

    dMatrix invGRIT(GRIT);
    invGRIT.invertLUdecomposition();
    if (m_verbose) {
      string MatrixTitle = "Inverse of Generalized rotation inertia tensor:";
      invGRIT.print(MatrixTitle, ctest);
    }
    // Locate the reduced moment of inertia associated with the bond.

    size_t idx(3);
    for (size_t i(0); i < m_RotBondIDs.size() && m_RotBondIDs[i] != bondID; i++, idx++);

    return invGRIT[idx][idx];
  }

  // Apply inertia weighting to the raw internal rotation velocity vector.
  void gStructure::ApplyInertiaWeighting(vector<string> &atomset, vector<double> &velocity, double fctr) const
  {
    for (size_t j(0); j < atomset.size(); j++)
    {
      size_t i = 3 * size_t(getAtomicOrder(atomset[j]));
      for (size_t n(0); n < 3; i++, n++) {
        velocity[i] *= fctr;
      }
    }
  }

  // Returns the rotational constants (in cm-1) in a vector.
  vector<double> gStructure::CalcRotConsts()
  {
    vector<double> RotConsts(m_PrincipalMI.size(), 0.0); //cm-1
    if (NumAtoms() < 2)
      return RotConsts; //empty
    for (size_t i = 0; i < m_PrincipalMI.size(); ++i) {
      RotConsts[i] = (m_PrincipalMI[i] == 0.0) ? 0.0 : conMntInt2RotCnt / m_PrincipalMI[i];
    }

    return RotConsts;
  }

  bool gStructure::GetConstituentAtomThermo
    (double& ZPE, double&Hf0, double& Hf298, double& dH298, double& dStdH298)
  {
    //TODO check non-initialized var?
    map<string, int> Comp = GetElementalComposition();
    if (Comp.empty())
    {
      cerr << "To mix thermodynamic and computational energies "
        "the molecule needs chemical structure (an atomList at least)" << endl;
      return false;
    }
    ZPE = Hf0 = Hf298 = dH298 = dStdH298 = 0.0;
    // for (auto c : Comp)   C++11!
    for (map<string, int>::iterator it = Comp.begin(); it != Comp.end(); ++it)
    { 

      //get values from librarymols.xml
      PersistPtr ppMol = GetFromLibrary(it->first, PersistPtr());
      if (ppMol)
      {
        //All must be in kJ/mol - not checked
        ZPE      += it->second * ppMol->XmlReadPropertyDouble("me:ZPE", optional);//computational
        Hf0      += it->second * ppMol->XmlReadPropertyDouble("me:Hf0", optional);
        Hf298    += it->second * ppMol->XmlReadPropertyDouble("me:Hf298", optional);
        dH298    += it->second * ppMol->XmlReadPropertyDouble("me:H0-H298", optional);
        dStdH298 += it->second * ppMol->XmlReadPropertyDouble("me:stdH0-H298", optional);
      }
      else
      {
        cerr << "librarymols.xml does not have an entry for the atom " << it->first << endl;
        return false;
      }
    }
    return !(IsNan(ZPE) || IsNan(Hf298) || IsNan(dH298));
  }

  // Returns an ordered array of masses.
  void gStructure::getAtomicMasses(vector<double> &AtomicMasses) const {
    AtomicMasses.clear();
    for (size_t i(0); i < m_atomicOrder.size(); i++){
      double mass = atomMass((Atoms.find(m_atomicOrder[i]))->second.element);
      AtomicMasses.push_back(mass);
    }
  }

  // Returns an ordered array of coordinates.
  void gStructure::getAtomicCoords(vector<double> &coords, AxisLabel cartLabel) const {
    coords.clear();
    for (size_t i(0); i < m_atomicOrder.size(); i++){
      double coord = (Atoms.find(m_atomicOrder[i]))->second.coords[cartLabel];
      coords.push_back(coord);
    }
  }

  // Returns an ordered array of X coordinates.
  void gStructure::getXCoords(vector<double> &coords) const {
    getAtomicCoords(coords, X);
  }

  // Returns an ordered array of Y coordinates.
  void gStructure::getYCoords(vector<double> &coords) const {
    getAtomicCoords(coords, Y);
  }

  // Returns an ordered array of Z coordinates.
  void gStructure::getZCoords(vector<double> &coords) const {
    getAtomicCoords(coords, Z);
  }

	//convert the string contain 3n-6-t non-tortional modes to a (3n-6-t)*(3n) matrix, which is a two-dimensional vector for convenience
	void gStructure::nonTortionalVectors(std::string modesStr, vector<vector<double>> &modes)
	{
		int i,j;
		vector<vector<int>> atomsInModes;
		OpenBabel::vector3 coord_I, coord_J, coord_K, coord_L;
		double RIJ1,RIJ2,RIJ3,RJK1,RJK2,RJK3,RKL1,RKL2,RKL3,RIJ,RJK,RKL,PIJK,CIJK,SIJK,PJKL,CJKL,SJKL,CIJKL,SIJKL;

		//read in the modes from modesStr
		for(i = 0; i < modesStr.length(); i++)
		{
			if(modesStr[i] == '\t')
			{
				modesStr[i] = ' ';
			}
		}

		std::string::size_type pos1, pos2, pos3;
		pos1 = 0;
		pos2 = modesStr.find(' ');
		pos3 = modesStr.find('-');;
		if(std::string::npos != pos3)
		{
			atomsInModes.push_back(vector<int>());
		}
		while(std::string::npos != pos3)
		{
			while(std::string::npos != pos3 && pos3<=pos2)
			{
				(atomsInModes.back()).push_back(atoi((modesStr.substr(pos1, pos3-pos1)).c_str()));
				pos1 = pos3+1;
				pos3 = modesStr.find('-', pos1);
			}
			if(atomsInModes.back().size() > 0)
			{
				atomsInModes.back().push_back(atoi((modesStr.substr(pos1, pos2-pos1)).c_str()));
				if(std::string::npos != pos3)
				{
					atomsInModes.push_back(vector<int>());							
				}
			}

			if(std::string::npos != pos3)
			{
				pos1=pos2+1;
				pos2=modesStr.find(' ', pos1);
				if(std::string::npos == pos2)
				{
					pos2 = modesStr.length();
				}
			}
		}

		////print modes details for test
		//for(i=0; i < atomsInModes.size(); i++)
		//{
		//	for(j=0; j < atomsInModes[i].size(); j++)
		//	{
		//		printf("%d\t", atomsInModes[i][j]);
		//	}
		//	printf("\n");
		//}

		//convert the info about label of atoms in modes to vectors about vibration modes 
		modes = vector<vector<double>>(atomsInModes.size(), vector<double>(3*NumAtoms(), 0.0));
		//details follows to assign values in modes
		for(i=0; i < atomsInModes.size(); i++)
		{
			switch(atomsInModes[i].size())
			{
			case 2:
				//bond vibration for stretch
				coord_I = Atoms[m_atomicOrder[atomsInModes[i][0]-1]].coords;
				coord_J = Atoms[m_atomicOrder[atomsInModes[i][1]-1]].coords;
				RIJ1 = coord_J.x() - coord_I.x();
				RIJ2 = coord_J.y() - coord_I.y();
				RIJ3 = coord_J.z() - coord_I.z();
				RIJ = (coord_I - coord_J).length();

				modes[i][3*atomsInModes[i][0]-3] = -RIJ1/RIJ;
				modes[i][3*atomsInModes[i][0]-2] = -RIJ2/RIJ;			
				modes[i][3*atomsInModes[i][0]-1] = -RIJ3/RIJ;
				modes[i][3*atomsInModes[i][1]-3] = RIJ1/RIJ;
				modes[i][3*atomsInModes[i][1]-2] = RIJ2/RIJ;			
				modes[i][3*atomsInModes[i][1]-1] = RIJ3/RIJ;

				//printf("%d\t%d\t%d\n", i+1, atomsInModes[i][0], atomsInModes[i][1]);
				//printf("%f\t%f\t%f\n", coord_I.x(), coord_I.y(), coord_I.z());
				//printf("%f\t%f\t%f\n", coord_J.x(), coord_J.y(), coord_J.z());
				//printf("%f\t%f\t%f\t%f\n", RIJ1, RIJ2, RIJ3, RIJ);
				
				break;
			case 3:
				//linear bend
				coord_I = Atoms[m_atomicOrder[atomsInModes[i][0]-1]].coords;
				coord_J = Atoms[m_atomicOrder[atomsInModes[i][1]-1]].coords;
				coord_K = Atoms[m_atomicOrder[atomsInModes[i][2]-1]].coords;
				RIJ1 = coord_J.x() - coord_I.x();
				RIJ2 = coord_J.y() - coord_I.y();
				RIJ3 = coord_J.z() - coord_I.z();
				RJK1 = coord_K.x() - coord_J.x();
				RJK2 = coord_K.y() - coord_J.y();
				RJK3 = coord_K.z() - coord_J.z();
				RIJ = (coord_I - coord_J).length();
				RJK = (coord_J - coord_K).length();
				PIJK = OpenBabel::dot((coord_J - coord_I), (coord_J - coord_K))/RIJ/RJK;
				if(PIJK > 1)
				{
					PIJK =1.0;
				}
				else if(PIJK < -1)
				{
					PIJK = -1.0;
				}
				PIJK = acos(PIJK);
				CIJK = cos(PIJK);
				SIJK = sin(PIJK);

				modes[i][3*atomsInModes[i][0]-3] = (-CIJK*RIJ1*RJK-RIJ*RJK1)/(SIJK*pow(RIJ,2)*RJK);
				modes[i][3*atomsInModes[i][0]-2] = (-CIJK*RIJ2*RJK-RIJ*RJK2)/(SIJK*pow(RIJ,2)*RJK);
				modes[i][3*atomsInModes[i][0]-1] = (-CIJK*RIJ3*RJK-RIJ*RJK3)/(SIJK*pow(RIJ,2)*RJK);
				modes[i][3*atomsInModes[i][1]-3] = (-RIJ*RIJ1*RJK+CIJK*RIJ1*pow(RJK,2)-CIJK*pow(RIJ,2)*RJK1+RIJ*RJK*RJK1)/(SIJK*pow(RIJ,2)*pow(RJK,2));
				modes[i][3*atomsInModes[i][1]-2] = (-RIJ*RIJ2*RJK+CIJK*RIJ2*pow(RJK,2)-CIJK*pow(RIJ,2)*RJK2+RIJ*RJK*RJK2)/(SIJK*pow(RIJ,2)*pow(RJK,2));
				modes[i][3*atomsInModes[i][1]-1] = (-RIJ*RIJ3*RJK+CIJK*RIJ3*pow(RJK,2)-CIJK*pow(RIJ,2)*RJK3+RIJ*RJK*RJK3)/(SIJK*pow(RIJ,2)*pow(RJK,2));
				modes[i][3*atomsInModes[i][2]-3] = (RIJ1*RJK+CIJK*RIJ*RJK1)/(SIJK*RIJ*pow(RJK,2));
				modes[i][3*atomsInModes[i][2]-2] = (RIJ2*RJK+CIJK*RIJ*RJK2)/(SIJK*RIJ*pow(RJK,2));
				modes[i][3*atomsInModes[i][2]-1] = (RIJ3*RJK+CIJK*RIJ*RJK3)/(SIJK*RIJ*pow(RJK,2));

				//printf("%d\t%d\t%d\t%d\n", i+1, atomsInModes[i][0], atomsInModes[i][1], atomsInModes[i][2]);
				//printf("%f\t%f\t%f\n", coord_I.x(), coord_I.y(), coord_I.z());
				//printf("%f\t%f\t%f\n", coord_J.x(), coord_J.y(), coord_J.z());
				//printf("%f\t%f\t%f\n", coord_K.x(), coord_K.y(), coord_K.z());
				//printf("%f\t%f\t%f\t%f\n", RIJ1, RIJ2, RIJ3, RIJ);
				//printf("%f\t%f\t%f\t%f\n", RJK1, RJK2, RJK3, RJK);
				//printf("%f\t%f\t%f\n", PIJK, CIJK, SIJK);
				//non-linear bend should be added in the future
				break;
			case 4:
				//rotations that haven't been included in the tortional modes
				coord_I = Atoms[m_atomicOrder[atomsInModes[i][0]-1]].coords;
				coord_J = Atoms[m_atomicOrder[atomsInModes[i][1]-1]].coords;
				coord_K = Atoms[m_atomicOrder[atomsInModes[i][2]-1]].coords;
				coord_L = Atoms[m_atomicOrder[atomsInModes[i][3]-1]].coords;
				RIJ1 = coord_J.x() - coord_I.x();
				RIJ2 = coord_J.y() - coord_I.y();
				RIJ3 = coord_J.z() - coord_I.z();
				RJK1 = coord_K.x() - coord_J.x();
				RJK2 = coord_K.y() - coord_J.y();
				RJK3 = coord_K.z() - coord_J.z();
				RKL1 = coord_L.x() - coord_K.x();
				RKL2 = coord_L.y() - coord_K.y();
				RKL3 = coord_L.z() - coord_K.z();
				RIJ = (coord_I - coord_J).length();
				RJK = (coord_J - coord_K).length();
				RKL = (coord_K - coord_L).length();
				PIJK = OpenBabel::dot((coord_J - coord_I), (coord_J - coord_K))/RIJ/RJK;
				if(PIJK > 1)
				{
					PIJK =1.0;
				}
				else if(PIJK < -1)
				{
					PIJK = -1.0;
				}
				PIJK = acos(PIJK);
				CIJK = cos(PIJK);
				SIJK = sin(PIJK);

				PJKL = OpenBabel::dot((coord_K - coord_J), (coord_K - coord_L))/RJK/RKL;
				if(PJKL > 1)
				{
					PJKL =1.0;
				}
				else if(PJKL < -1)
				{
					PJKL = -1.0;
				}
				PJKL = acos(PJKL);
				CJKL = cos(PJKL);
				SJKL = sin(PJKL);
				CIJKL=((-RIJ2*RJK1+RIJ1*RJK2)*(-RJK2*RKL1+RJK1*RKL2)+(RIJ3*RJK1-RIJ1*RJK3)*(RJK3*RKL1-RJK1*RKL3)+(-RIJ3*RJK2+RIJ2*RJK3)*(-RJK3*RKL2+RJK2*RKL3))/(SIJK*SJKL*RIJ*RJK*RJK*RKL);
				SIJKL=((-RIJ3*RJK2+RIJ2*RJK3)*RKL1+(RIJ3*RJK1-RIJ1*RJK3)*RKL2+(-(RIJ2*RJK1)+RIJ1*RJK2)*RKL3)/(RIJ*RJK*RKL*SIJK*SJKL);

				modes[i][3*atomsInModes[i][0]-3] = SIJKL*RIJ1/(CIJKL*pow(SIJK,2)*pow(RIJ,2))+CIJK*SIJKL*RJK1/(CIJKL*pow(SIJK,2)*RIJ*RJK)+(RJK3*RKL2-RJK2*RKL3)/(CIJKL*SIJK*SJKL*RIJ*RJK*RKL);
				modes[i][3*atomsInModes[i][0]-2] = SIJKL*RIJ2/(CIJKL*pow(SIJK,2)*pow(RIJ,2))+CIJK*SIJKL*RJK2/(CIJKL*pow(SIJK,2)*RIJ*RJK)+(-RJK3*RKL1+RJK1*RKL3)/(CIJKL*SIJK*SJKL*RIJ*RJK*RKL);
				modes[i][3*atomsInModes[i][0]-1] = SIJKL*RIJ3/(CIJKL*pow(SIJK,2)*pow(RIJ,2))+CIJK*SIJKL*RJK3/(CIJKL*pow(SIJK,2)*RIJ*RJK)+(RJK2*RKL1-RJK1*RKL2)/(CIJKL*SIJK*SJKL*RIJ*RJK*RKL);
				modes[i][3*atomsInModes[i][1]-3] =-(SIJKL*RIJ1/(CIJKL*pow(SIJK,2)*pow(RIJ,2)))+CIJK*SIJKL*(RIJ1-RJK1)/(CIJKL*pow(SIJK,2)*RIJ*RJK)-SIJKL*RJK1/(CIJKL*pow(RJK,2))+SIJKL*RJK1/(CIJKL*pow(SIJK,2)*pow(RJK,2))+SIJKL*RJK1/(CIJKL*pow(SJKL,2)*pow(RJK,2))+CJKL*SIJKL*RKL1/(CIJKL*pow(SJKL,2)*RJK*RKL)+(-RIJ3*RKL2-RJK3*RKL2+RIJ2*RKL3+RJK2*RKL3)/(CIJKL*SIJK*SJKL*RIJ*RJK*RKL);
				modes[i][3*atomsInModes[i][1]-2] =-SIJKL*RIJ2/(CIJKL*pow(SIJK,2)*pow(RIJ,2))+CIJK*SIJKL*(RIJ2-RJK2)/(CIJKL*pow(SIJK,2)*RIJ*RJK)-SIJKL*RJK2/(CIJKL*pow(RJK,2))+SIJKL*RJK2/(CIJKL*pow(SIJK,2)*pow(RJK,2))+SIJKL*RJK2/(CIJKL*pow(SJKL,2)*pow(RJK,2))+CJKL*SIJKL*RKL2/(CIJKL*pow(SJKL,2)*RJK*RKL)+(RIJ3*RKL1+RJK3*RKL1-RIJ1*RKL3-RJK1*RKL3)/(CIJKL*SIJK*SJKL*RIJ*RJK*RKL); 
				modes[i][3*atomsInModes[i][1]-1] =-SIJKL*RIJ3/(CIJKL*pow(SIJK,2)*pow(RIJ,2))+CIJK*SIJKL*(RIJ3-RJK3)/(CIJKL*pow(SIJK,2)*RIJ*RJK)-SIJKL*RJK3/(CIJKL*pow(RJK,2))+SIJKL*RJK3/(CIJKL*pow(SIJK,2)*pow(RJK,2))+SIJKL*RJK3/(CIJKL*pow(SJKL,2)*pow(RJK,2))+(-RIJ2*RKL1-RJK2*RKL1+RIJ1*RKL2+RJK1*RKL2)/(CIJKL*SIJK*SJKL*RIJ*RJK*RKL)+CJKL*SIJKL*RKL3/(CIJKL*pow(SJKL,2)*RJK*RKL);   
				modes[i][3*atomsInModes[i][2]-3] =-CIJK*SIJKL*RIJ1/(CIJKL*pow(SIJK,2)*RIJ*RJK)+SIJKL*RJK1/(CIJKL*pow(RJK,2))-SIJKL*RJK1/(CIJKL*pow(SIJK,2)*pow(RJK,2))-SIJKL*RJK1/(CIJKL*pow(SJKL,2)*pow(RJK,2))+CJKL*SIJKL*(RJK1-RKL1)/(CIJKL*pow(SJKL,2)*RJK*RKL)+SIJKL*RKL1/(CIJKL*pow(SJKL,2)*pow(RKL,2))+(RIJ3*RJK2-RIJ2*RJK3+RIJ3*RKL2-RIJ2*RKL3)/(CIJKL*SIJK*SJKL*RIJ*RJK*RKL);
				modes[i][3*atomsInModes[i][2]-2] =-CIJK*SIJKL*RIJ2/(CIJKL*pow(SIJK,2)*RIJ*RJK)+SIJKL*RJK2/(CIJKL*pow(RJK,2))-SIJKL*RJK2/(CIJKL*pow(SIJK,2)*pow(RJK,2))-SIJKL*RJK2/(CIJKL*pow(SJKL,2)*pow(RJK,2))+CJKL*SIJKL*(RJK2-RKL2)/(CIJKL*pow(SJKL,2)*RJK*RKL)+SIJKL*RKL2/(CIJKL*pow(SJKL,2)*pow(RKL,2))+(-RIJ3*RJK1+RIJ1*RJK3-RIJ3*RKL1+RIJ1*RKL3)/(CIJKL*SIJK*SJKL*RIJ*RJK*RKL);
				modes[i][3*atomsInModes[i][2]-1] =-CIJK*SIJKL*RIJ3/(CIJKL*pow(SIJK,2)*RIJ*RJK)+SIJKL*RJK3/(CIJKL*pow(RJK,2))-SIJKL*RJK3/(CIJKL*pow(SIJK,2)*pow(RJK,2))-SIJKL*RJK3/(CIJKL*pow(SJKL,2)*pow(RJK,2))+(RIJ2*RJK1-RIJ1*RJK2+RIJ2*RKL1-RIJ1*RKL2)/(CIJKL*SIJK*SJKL*RIJ*RJK*RKL)+CJKL*SIJKL*(RJK3-RKL3)/(CIJKL*pow(SJKL,2)*RJK*RKL)+SIJKL*RKL3/(CIJKL*pow(SJKL,2)*pow(RKL,2));
				modes[i][3*atomsInModes[i][3]-3] =-CJKL*SIJKL*RJK1/(CIJKL*pow(SJKL,2)*RJK*RKL)+(-RIJ3*RJK2+RIJ2*RJK3)/(CIJKL*SIJK*SJKL*RIJ*RJK*RKL)-SIJKL*RKL1/(CIJKL*pow(SJKL,2)*pow(RKL,2)); 
				modes[i][3*atomsInModes[i][3]-2] =-CJKL*SIJKL*RJK2/(CIJKL*pow(SJKL,2)*RJK*RKL)+(RIJ3*RJK1-RIJ1*RJK3)/(CIJKL*SIJK*SJKL*RIJ*RJK*RKL)-SIJKL*RKL2/(CIJKL*pow(SJKL,2)*pow(RKL,2));
				modes[i][3*atomsInModes[i][3]-1] =(-RIJ2*RJK1+RIJ1*RJK2)/(CIJKL*SIJK*SJKL*RIJ*RJK*RKL)-CJKL*SIJKL*RJK3/(CIJKL*pow(SJKL,2)*RJK*RKL)-SIJKL*RKL3/(CIJKL*pow(SJKL,2)*pow(RKL,2));

				//printf("%d\t%d\t%d\t%d\t%d\n", i+1, atomsInModes[i][0], atomsInModes[i][1], atomsInModes[i][2], atomsInModes[i][3]);
				//printf("%f\t%f\t%f\n", coord_I.x(), coord_I.y(), coord_I.z());
				//printf("%f\t%f\t%f\n", coord_J.x(), coord_J.y(), coord_J.z());
				//printf("%f\t%f\t%f\n", coord_K.x(), coord_K.y(), coord_K.z());
				//printf("%f\t%f\t%f\n", coord_L.x(), coord_L.y(), coord_L.z());
				//printf("%f\t%f\t%f\t%f\n", RIJ1, RIJ2, RIJ3, RIJ);
				//printf("%f\t%f\t%f\t%f\n", RJK1, RJK2, RJK3, RJK);
				//printf("%f\t%f\t%f\t%f\n", RKL1, RKL2, RKL3, RKL);
				//printf("%f\t%f\t%f\n", PIJK, CIJK, SIJK);
				//printf("%f\t%f\t%f\n", PJKL, CJKL, SJKL);
				//printf("%f\t%f\n", CIJKL, SIJKL);

				break;
			default:
				cerr << "non-tortional modes cannot be recognized because the number of atoms is " << atomsInModes[i].size() << endl;
				break;
			}
		}
	}

  // Export to xmol format.
  void gStructure::exportToXYZ(const char* txt, bool last, PersistPtr ppConfigData){

    // Only write something if the verbosity flag has been set.
    if (!m_verbose)
      return;

    cinfo << Atoms.size() << endl ;
    cinfo << getHost()->getName();
    if(txt)
      cinfo << "-" << txt;
    cinfo << endl ;
    map<string, atom>::const_iterator iter;
    for (iter=Atoms.begin(); iter!=Atoms.end(); ++iter) {
      const atom &at = iter->second ;
      cinfo << setw(10) << at.element ;
      cinfo << formatFloat(at.coords.x(), 6, 15)  
        << formatFloat(at.coords.y(), 6, 15) 
        << formatFloat(at.coords.z(), 6, 15) << endl ;
    }
    exportToCML(txt, last, ppConfigData);
  }

  void gStructure::exportToCML(const char* txt, bool last, PersistPtr ppConfigData){
    static PersistPtr pp;
    if(ppConfigData)
      pp = ppConfigData;
    PersistPtr ppMol = pp->XmlWriteElement("molecule");
    string id = getHost()->getName();
    ppMol->XmlWriteAttribute("id", id  + '-' + txt);
    PersistPtr ppAtList = ppMol->XmlWriteElement("atomArray");
    map<string, atom>::const_iterator iter;
    for (iter=Atoms.begin(); iter!=Atoms.end(); ++iter) {
      const atom &at = iter->second ;
      PersistPtr ppAt = ppAtList->XmlWriteElement("atom");
      ppAt->XmlWriteAttribute("id", at.id);
      ppAt->XmlWriteAttribute("elementType", at.element);
      ppAt->XmlWriteAttribute("x3", at.coords.x(), 4, true);
      ppAt->XmlWriteAttribute("y3", at.coords.y(), 4, true);
      ppAt->XmlWriteAttribute("z3", at.coords.z(), 4, true);
    }
    PersistPtr ppBList = ppMol->XmlWriteElement("bondArray");
    for (map<string, pair<string, string> >::const_iterator iter=Bonds.begin();
      iter!=Bonds.end(); ++iter) {
        PersistPtr ppB = ppBList->XmlWriteElement("bond");
        ppB->XmlWriteAttribute("id", iter->first);
        ppB->XmlWriteAttribute("atomRefs2", iter->second.first + " " + iter->second.second);
        ppB->XmlWriteAttribute("order", "1");
        if(iter->first==m_RotBondIDs.back())
          ppB->XmlWriteAttribute("color", "pink");
    }
  }
}//namespace

//XML written as text; not used.
//void gStructure::exportToCML(const char* txt, bool last, PersistPtr pp){
//  static stringstream ss;
//  static PersistPtr ppConfigData;
//  if(pp)
//    ppConfigData = pp;
//  ss << " <molecule id=\"" << getHost()->getName();
//  if(txt)
//    ss <<'-' << txt;
//  ss << "\">\n" << "  <atomArray>\n";
//  map<string, atom>::const_iterator iter;
//  for (iter=Atoms.begin(); iter!=Atoms.end(); ++iter) {
//    const atom &at = iter->second ;
//    ss  << "   <atom id=" << at.id << " elementType=\"" << at.element << "\" ";
//    ss << setprecision(4) << fixed
//          << "x3=\"" << at.coords.x() <<"\" "
//          << "y3=\"" << at.coords.y() <<"\" "
//          << "z3=\"" << at.coords.z() <<"\"/>\n";
//  }
//  ss << "   </atomArray>\n   <bondArray>";
//  int i(1);
//  for (map<string, pair<string, string> >::const_iterator iter=Bonds.begin();
//         iter!=Bonds.end(); ++iter, ++i) {
//    if(i%2) ss << '\n';
//    ss << "   <bond atomRefs2=\"" << iter->second.first
//          <<' ' << iter->second.second << " order=\"1\"/>";
//  }
//  ss << "\n  </bondArray>\n </molecule>" << endl;
//  if(last)
//    // Written as text (rather than XML) because it is available.
//    // Should be seen ok as valid XML.
//    ppConfigData->XmlWrite(ss.str());
//  ss.str()="";
//}
