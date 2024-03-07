/**
  NZ2HDM_Tsq_4 JMC
*/

#include <iostream>

#include "models/NZ2HDM_Tsq_4.hpp" // Located in effective-potential/include/models
#include "transition_finder.hpp"
#include "phase_finder.hpp"
#include "phase_plotter.hpp"
#include "potential_plotter.hpp"
#include "potential_line_plotter.hpp"
#include "logger.hpp"


int main(int argc, char* argv[]) {

  const bool debug_mode = argc > 1 and strcmp(argv[1],"-d")==0;

  // Set level of screen  output
  if (debug_mode) {
      LOGGER(debug);
  } else {
      LOGGER(fatal);
  }

  // Construct our model
  EffectivePotential::NZ2HDM_Tsq_4 model;
  
  const double v = 246.;
  const double mh = 125.;
  const double mH = 400.;
  const double mA = 400.;
  const double mHc = 400.;
  const double m3in = 200.;
  const double tb = 1.;
  const double cbma = 0.;
  const double ms = 130.;
  const double lS = 1.;
  const double lbeta = 0.5;
  const double lS2 = 1.5*lbeta;
  const double lS3_real = 0.;
  const double lS3_imag = 0.115;
  const double l6 = 0.;
  const double l7 = 0.;
  
  model.set_v(v);
  model.set_tb(tb);
  model.set_lS(lS);
  model.set_lS1(tb, lbeta, lS2, lS3_real);
  model.set_lS2(lS2);
  model.set_lS3_real(lS3_real);
  model.set_lS3_imag(lS3_imag);
  model.set_mS(ms, lbeta);
  model.set_m3(m3in, mH, tb);
  model.set_m1(mh, mH, tb, cbma);
  model.set_m2(mh, mH, tb, cbma);
  model.set_l1(mh, mH, tb, cbma);
  model.set_l2(mh, mH, tb, cbma);
  model.set_l3(mh, mH, mHc, tb, cbma);
  model.set_l4(mH, mA, mHc, tb);
  model.set_l5(mH, mA, tb);
  model.set_l6(l6);
  model.set_l7(l7);
  
  model.check_param();

  // Make PhaseFinder object and find the phases
  PhaseTracer::PhaseFinder pf(model);
  pf.set_n_ew_scalars(3);
//  pf.set_dt_max_abs(1);
//  pf.set_dt_start_rel(0.001);
//  pf.set_dt_max_rel(0.001);
//  pf.set_x_abs_jump(0.01);
//  pf.set_x_rel_jump(0.001);
//  pf.set_x_abs_identical(0.01);
//  pf.set_x_rel_identical(0.001);
  pf.set_hessian_singular_rel_tol(5e-3);
  pf.find_phases();
  std::cout <<pf;

  // Make TransitionFinder object and find the transitions
  PhaseTracer::TransitionFinder tf(pf);
  tf.set_n_ew_scalars(3);
  tf.find_transitions();
  std::cout <<std::setprecision (15) <<tf;

  if (debug_mode) {
//    PhaseTracer::potential_plotter(model, tf.get_transitions().front().TC, "NZ2HDM_Tsq_4", 0., 2., 0.01, -2., 0., 0.01);
    PhaseTracer::potential_line_plotter(model, tf.get_transitions(), "NZ2HDM_Tsq_4");
    PhaseTracer::phase_plotter(tf, "NZ2HDM_Tsq_4");
  }
}
