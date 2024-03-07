// ====================================================================
// This file is part of PhaseTracer

// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
// ====================================================================

#ifndef POTENTIAL_TEST_MODEL_HPP_INCLUDED
#define POTENTIAL_TEST_MODEL_HPP_INCLUDED

/**
   NZ2HDM T square corrections JMC
*/

#include <vector>
#include <complex>
#include <cmath>
#include <Eigen/Eigenvalues>

#include "one_loop_potential.hpp"
#include "pow.hpp"


namespace EffectivePotential {


class NZ2HDM_Tsq_4 : public Potential {


 public:
  
  
  void set_v(double v_) {v = v_;}
  void set_tb(double tb_) {tb = tb_;}
  void set_lS(double lS_) {lS = lS_;}
  void set_lS1(double tb_, double lbeta_, double lS2_, double lS3_) {lS1 = ( 1 + square(tb_) ) * lbeta_ - lS2_ * square(tb) - lS3_ * square(tb);}
  void set_lS2(double lS2_) {lS2 = lS2_;}
  void set_lS3_real(double lS3_) {lS3.real(lS3_);}
  void set_lS3_imag(double lS3_) {lS3.imag(lS3_);}
  void set_mS(double ms_, double lbeta_) {
  	mS2 = square(ms_) - lbeta_ * square(v); }
  void set_m3(double m3in_, double mH_, double tb) {
  	m32 = ( square(mH_) - square(m3in_) ) / ( tb + 1. / tb );}
  void set_m1(double mh_, double mH_, double tb, double cbma_) {
  	m12 = m32 * tb - .5 * ( square(mh_) + ( square(mH_) - square(mh_) ) * cbma_ * ( cbma_ + sqrt( 1. - square(cbma_) ) * tb ) ); }
  void set_m2(double mh_, double mH_, double tb, double cbma_) {
  	m22 = m32 / tb - .5 * ( square(mh_) + ( square(mH_) - square(mh_) ) * cbma_ * ( cbma_ - sqrt( 1. - square(cbma_) ) / tb ) ); }
  void set_l1(double mh_, double mH_, double tb, double cbma_) {
  	l1 = square(mh_) - tb * ( m32 - square(mH_) / ( tb + 1./tb ) ) * ( square(tb) + 1. ) + ( square(mh_) - square(mH_) ) * ( square(cbma_) * ( square(tb) - 1. ) - 2. * tb * cbma_ * sqrt(1.-square(cbma_)) ); l1/=square(v); }
  void set_l2(double mh_, double mH_, double tb, double cbma_) {
  	l2 = square(mh_) - ( m32 - square(mH_) / ( tb + 1./tb ) ) * ( 1. / square(tb) + 1. ) / tb + ( square(mh_) - square(mH_) ) * ( square(cbma_) * ( 1. / square(tb) - 1. ) - 2. / tb * cbma_ * sqrt(1.-square(cbma_)) ); l2/=square(v); }
  void set_l3(double mh_, double mH_, double mHc_, double tb, double cbma_) {
  	l3 = square(mh_) + 2 * square(mHc_) - 2 * square(mH_) - ( m32 - square(mH_) / ( tb + 1./tb ) ) * ( tb + 1./tb ) - ( square(mh_) - square(mH_) ) * ( 2. * square(cbma_) + cbma_ * sqrt(1.-square(cbma_)) * ( tb + 1./tb ) ); l3/=square(v); }
  void set_l4(double mH_, double mA_, double mHc_, double tb) {
  	l4 = square(mA_) - 2 * square(mHc_) + square(mH_) + ( m32 - square(mH_) / ( tb + 1./tb ) ) * ( tb + 1./tb ); l4/=square(v); }
  void set_l5(double mH_, double mA_, double tb) {
  	l5 = square(mH_) - square(mA_) + ( m32 - square(mH_) / ( tb + 1./tb ) ) * ( tb + 1./tb ); l5/=square(v); }
  void set_l6(double l6_) {l6 = l6_;}
  void set_l7(double l7_) {l7 = l7_;}
  
  
  void check_param() {
  	std::cout<<"m12 "<<m12<<" m22 "<<m22<<" mS2 "<<mS2<<" m32 "<<m32<<std::endl
  	<<"l1 "<<l1<<" l2 "<<l2<<" l3 "<<l3<<" l4 "<<l4<<" l5 "<<l5<<" l6 "<<l6<<" l7 "<<l7<<std::endl
  	<<"lS "<<lS<<" lS1 "<<lS1<<" lS2 "<<lS2<<" lS3 "<<lS3<<std::endl<<std::endl; }
  
  
  double V(Eigen::VectorXd phi, double T) const override {
  
  
  
//    double p1 = phi[0];
//    double p2 = phi[1];
//    double s = phi[2];
//    double delta = 0;
//    double delta = phi[3];
    
//    double mh2 = ( 3.*l1 + 2.*l3 + l4 + lS1 ) * square(p1) + (3.*l2 + 2.*l3 + l4 + lS2) * square(p2) + (3.*lS + 4.*(lS1 + lS2)) * square(s) + lS3 * p1 * p2 * sin(delta);
//    double mw2 = square(g)*(p1*p1+p2*p2)/4.;
//    double mz2 = (square(g)+square(gp))*(p1*p1+p2*p2)/4.;
//    double mt2 = square(mt/v) * (1+1./square(tb)) * (p2*p2);
//    
//    return .5 * ( m12*square(p1) + m22*square(p2) )
//           - m32*p1*p2*cos(delta)
//           + .125 * ( l1*pow_4(p1) + l2*pow_4(p2) )
//           + .25 * (l3+l4)*square(p1)*square(p2)
//           + .5 * mS2*square(s) + 0.25 * lS*pow_4(s)
//           + .5 * ( lS1*square(p1)*square(s) + lS2*square(p2)*square(s) + lS3*p1*p2*square(s)*sin(delta) )
//           + square(T)/24. * ( mh2 + 3. * ( 2. * ( mw2 + mt2 ) + mz2 ) );



//    double h1 = phi[0];
//    double h2 = phi[1];
//    double a2 = phi[2];
//    double s = phi[3];
//    
//    double mh2 = ( 3.*l1 + 2.*l3 + l4 + lS1 ) * (square(h1)+square(a1)) + (3.*l2 + 2.*l3 + l4 + lS2) * (square(h2)+square(a2)) + (3.*lS + 4.*(lS1 + lS2)) * square(s) + lS3 * (a1*h2-h1*a2);
//    double mw2 = square(g)*(h1*h1+a1*a1+h2*h2+a2*a2)/4.;
//    double mz2 = (square(g)+square(gp))*(h1*h1+a1*a1+h2*h2+a2*a2)/4.;
//    double mt2 = square(mt/v) * (1+1./square(tb)) * (h2*h2+a2*a2);
//    
//    return .5 * ( m12*(square(h1)+square(a1)) + m22*(square(h2)+square(a2)) )
////           - m32*h1*h2
//           - m32*(h1*h2+a1*a2)
//           + .125 * ( l1*square((square(h1)+square(a1))) + l2*square(square(h2)+square(a2)))
//           + .25 * (l3+l4)*(square(h1)+square(a1))*(square(h2)+square(a2))
//           + .5 * mS2*square(s) + 0.25 * lS*pow_4(s)
//           + .5 * ( lS1*(square(h1)+square(a1))*square(s) + lS2*(square(h2)+square(a2))*square(s) + lS3*(a1*h2-h1*a2)*square(s) )
//           + square(T)/24. * ( mh2 + 3. * ( 2. * ( mw2 + mt2 ) + mz2 ) );
  
  
  
    std::complex<double> phi1;
    phi1.real(phi[0]/sqrt(2));
    phi1.imag(0);
  
    std::complex<double> phi2;
    phi2.real(phi[1]/sqrt(2));
    phi2.imag(phi[2]/sqrt(2));
    std::complex<double> phi1C = std::conj(phi1);
    std::complex<double> phi2C = std::conj(phi2);
    double s = phi[3];
    
    double mh2 = real(2.*( 3.*l1 + 2.*l3 + l4 + lS1 ) * phi1C*phi1 + 2.*(3.*l2 + 2.*l3 + l4 + lS2) * phi2C*phi2 + (3.*lS + 4.*(lS1 + lS2)) * square(s) + lS3*phi1C*phi2 + std::conj(lS3*phi1C*phi2));
    double mw2 = real(square(g)*(phi1C*phi1 + phi2C*phi2)/2.);
    double mz2 = real((square(g)+square(gp))*(phi1C*phi1 + phi2C*phi2)/2.);
    double mt2 = real(2.*square(mt/v) * (1+1./square(tb)) * phi2C*phi2);
    
//    std::cout<<h1<<std::endl
//    	<<a1<<std::endl
//    	<<phi1<<std::endl
//    	<<phi1C<<std::endl
//    	<<real(m32)<<std::endl
//    	<<m32+std::conj(m32)<<std::endl
//    	<<real(m32)*(h1*h2+a1*a2)<<std::endl
//    	<<m32*phi1C*phi2+std::conj(m32*phi1C*phi2)<<std::endl<<std::endl;
    
    return real( .5 * ( l1*phi1C*phi1*phi1C*phi1 + l2*phi2C*phi2*phi2C*phi2 )
           + l3*phi1C*phi1*phi2C*phi2 + l4*phi1C*phi2*phi2C*phi1
           + .25*lS*s*s*s*s + lS1*phi1C*phi1*s*s + lS2*phi2C*phi2*s*s
           + m12*phi1C*phi1 + m22*phi2C*phi2 + .5*mS2*s*s
           + ( .5*l5*phi1C*phi2 + l6*phi1C*phi1 + l7*phi2C*phi2 ) * phi1C*phi2 + .5*lS3*phi1C*phi2*s*s - m32*phi1C*phi2
           + std::conj( ( .5*l5*phi1C*phi2 + l6*phi1C*phi1 + l7*phi2C*phi2 ) * phi1C*phi2 + .5*lS3*phi1C*phi2*s*s - m32*phi1C*phi2 ) )
           + square(T)/24. * ( mh2 + 3. * ( 2. * ( mw2 + mt2 ) + mz2 ) );
  }
  
  size_t get_n_scalars() const override { return 4; }


  std::vector<Eigen::VectorXd> apply_symmetry(Eigen::VectorXd phi) const override {
    
    auto S1 = phi;
    S1[0] = - phi[0];
    S1[1] = - phi[1];
    S1[2] = - phi[2];
    
    auto S2 = phi;
    S2[3] = - phi[3];
    
    return {S1, S2};
  }
  

 private:
  
//  double a1 = 0.;
  double g = .649368;
  double gp = .357880;
  double mt = 165.3;
  double v = 246.;
  double tb = 1.;
  double m12 = pow(100., 2);
  double m22 = pow(100., 2);
  double mS2 = pow(100., 2);
  double m32 = pow(50., 2);
  double l1 = .1;
  double l2 = .1;
  double l3 = .1;
  double l4 = .1;
  double l5 = .1;
  double l6 = 0.;
  double l7 = 0.;
  double lS = .1;
  double lS1 = .1;
  double lS2 = .1;
  std::complex<double> lS3;
  
  
};


}  // namespace EffectivePotential

#endif
