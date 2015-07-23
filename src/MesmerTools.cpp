//MesmerTools.cpp
#include "MesmerTools.h"

using namespace Constants;

namespace mesmer
{

  // translation contribution for the partition function of two molecules
  double translationalContribution(const double m1, const double m2, const double beta){
    // Translational contribution
    // 2.0593e19 = conversion factor,  1e-6*(((cm-1 -> j)/(h2*na)))^3/2
    // double tp_C = 2.0593e19 * pow(2. * M_PI ,1.5);

    return (tp_C * pow(m1 * m2 / ((m1 + m2) * beta), 1.5));
  }

  double canonicalPartitionFunction(const vector<double>& DOS, const vector<double>& Ene, const double beta){

    double CanPrtnFn(0.0) ;
    for (size_t i(0), j(DOS.size()-1); i < DOS.size(); i++, j--) {
      if (DOS[j] > 0.0)
        CanPrtnFn += exp( log(DOS[j]) - beta*Ene[j] ) ;
    }
    return CanPrtnFn;

  }

  double canonicalMeanEnergy(const vector<double>& DOS, const vector<double>& Ene, const double beta){

    double meanEnergy(0.0), CanPrtnFn(0.0) ;
    for (size_t i(0), j(DOS.size()-1); i < DOS.size(); i++, j--) {
      if (DOS[j] > 0.0) {
        double tmp  = exp( log(DOS[j]) - beta*Ene[j] ) ;
        CanPrtnFn  += tmp ;
        meanEnergy += Ene[j]*tmp ;
      }
    }
    return meanEnergy/CanPrtnFn ;

  }

  // Function to calculate the thermodynamic data and output to the .test file
  // prtnFn[0] is the partition function z1*z2*...*zj*...*zn
  // prtnFn[1] denotes for sum(z'[j]/z[j])
  // prtnFn[2] denotes for sum((z'[j]/z[j])')=sum(z''[j]/z[j]-(z'[j]/z[j])^2)
  // z'[j] is dz/d(1/T)
  void thermodynamicCalc(const double* prtnFn,const double beta, double MW, double * out_thermoData)
  {
	  double temp = 1.0/boltzmann_RCpK/beta;
	  double S, Cp, HmH0;

	  S = idealGasC/Calorie_in_Joule*(2.5+1.5*log(2*M_PI*MW/1000)-4*log(AvogadroC)-3*log(PlancksConstant_in_JouleSecond)-log(atm_in_pascal)+2.5*log(idealGasC*temp)+log(prtnFn[0])-prtnFn[1]/temp);
	  Cp = idealGasC/Calorie_in_Joule*(prtnFn[2]/temp/temp+2.5);
	  HmH0 = idealGasC/Calorie_in_Joule*(-prtnFn[1]+temp*2.5);
	  out_thermoData[0] = HmH0;
	  out_thermoData[1] = S;
	  out_thermoData[2] = Cp;
	  ctest << "temperature, Q, H(T)-H(0), S, and Cp ([cal][mol][K]):    " << temp << "    " << prtnFn[0] << "    " << HmH0 << "    " << S << "    " << Cp << endl;
  }

  //
  // Calculate the average grain energy and then number of states per grain.
  //
  void calcGrainAverages(const size_t &MaximumGrain, const size_t &cellPerGrain, const size_t &cellOffset, const vector<double>& CellDOS, const vector<double>& CellEne, vector<double>& grainDOS, vector<double>& grainEne)
  {
    grainEne.clear() ;
    grainDOS.clear() ;
    grainEne.resize(MaximumGrain, 0.) ;
    grainDOS.resize(MaximumGrain, 0.) ;

    // Check that there are enough cells.
    //if (GrainSize < 1) {
    //  throw (std::runtime_error("The number of Cells is insufficient to produce requested number of Grains.")); 
    //}

    size_t idx1 = 0 ;
    size_t idx2 = 0 ;
    for (size_t i(0) ; i < MaximumGrain ; ++i ) {

      size_t idx3(idx1);

	  // Account for the cell off set against PES grid by altering
	  // the range of first grain average.

	  const size_t cellRange = (i == 0) ? cellPerGrain - cellOffset : cellPerGrain ; 

      // Calculate the number of states in a grain.
      double gNOS = 0.0 ;
      for (size_t j(0) ; j < cellRange ; ++j, ++idx1 ){
        gNOS += CellDOS[idx1] ;
      }

      // Calculate average energy of the grain if it contains sum states.
       if ( gNOS > 0.0 ){
        double gSE = 0.0 ; // grain sum of state energy
        for (size_t j(0) ; j < cellRange ; ++j, ++idx3 ){
          gSE += CellEne[idx3] * CellDOS[idx3] ;
        }
        grainDOS[idx2] = gNOS ;
        grainEne[idx2] = gSE/gNOS ;
        idx2++ ;
      }

    }

    // Issue warning if number of grains produced is less that requested.

    if ( idx2 != MaximumGrain ) {
      cinfo << "Number of grains produced is not equal to that requested" << once << endl
        << "Number of grains requested: " << MaximumGrain << once << endl
        << "Number of grains produced : " << idx2 << once << endl;
    }
  }
}
