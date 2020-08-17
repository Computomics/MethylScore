/*
  Copyright (C) 2008 Cold Spring Harbor Laboratory
  Authors: Andrew D. Smith
  modified by Joerg Hagmann, Max Planck Institute, Tuebingen

  This file is part of methpipe.

  methpipe is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  methpipe is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with methpipe; if not, write to the Free Software
  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
*/

#ifndef TWO_STATE_HMM_HPP
#define TWO_STATE_HMM_HPP

#include "smithlab_utils.hpp"
#include <memory>

struct betabin;

class TwoStateHMMB {
public:
  
  TwoStateHMMB(const double mp, const double tol,
	       const size_t max_itr, const bool v, bool d = false) : 
    MIN_PROB(mp), tolerance(tol), max_iterations(max_itr),    
    VERBOSE(v), DEBUG(d) {}
  
  double
  ViterbiDecoding(const std::vector<std::pair<double, double> > &values,
		  const std::vector<char> &contexts,
		  const std::vector<size_t> &reset_points,
		  const std::vector<double> &start_trans, 
		  const std::vector<std::vector<double> > &trans, 
		  const std::vector<double> &end_trans, 
		  std::vector<std::vector<double> > &alphas,
		  std::vector<std::vector<double> > &betas,
		  std::vector<bool> &ml_classes) const;
  
  
  double
  BaumWelchTraining(const std::vector<std::pair<double, double> > &values,
		    const std::vector<char> &contexts,
		    const std::vector<size_t> &reset_points,
		    std::vector<double> &start_trans,
		    std::vector<std::vector<double> > &trans, 
		    std::vector<double> &end_trans,
		    std::vector<std::vector<double> > &alphas,
		    std::vector<std::vector<double> > &betas) const;
  
  double
  PosteriorDecoding(const std::vector<std::pair<double, double> > &values,
		    const std::vector<char> &contexts,
		    const std::vector<size_t> &reset_points,
		    const std::vector<double> &start_trans, 
		    const std::vector<std::vector<double> > &trans, 
		    const std::vector<double> &end_trans, 
		    std::vector<std::vector<double> > &alphas,
		    std::vector<std::vector<double> > &betas,
		    std::vector<bool> &classes,
		    std::vector<double> &llr_scores) const;
/*  
  void
  PosteriorScores(const std::vector<std::pair<double, double> > &values,
		  const std::vector<size_t> &reset_points,
		  const std::vector<double> &start_trans, 
		  const std::vector<std::vector<double> > &trans, 
		  const std::vector<double> &end_trans, 
		  const double fg_alpha, const double fg_beta,
		  const double bg_alpha, const double bg_beta,
		  const std::vector<bool> &classes,
		  std::vector<double> &llr_scores) const;
  
  void
  PosteriorScores(const std::vector<std::pair<double, double> > &values,
		  const std::vector<size_t> &reset_points,
		  const std::vector<double> &start_trans, 
		  const std::vector<std::vector<double> > &trans, 
		  const std::vector<double> &end_trans, 
		  const double fg_alpha, const double fg_beta,
		  const double bg_alpha, const double bg_beta,
		  const bool class_id,
		  std::vector<double> &llr_scores) const;
  
  void
  TransitionPosteriors(const std::vector<std::pair<double, double> > &values,
		       const std::vector<size_t> &reset_points,
		       const std::vector<double> &start_trans, 
		       const std::vector<std::vector<double> > &trans, 
		       const std::vector<double> &end_trans, 
		       const double fg_alpha, const double fg_beta,
		       const double bg_alpha, const double bg_beta,
		       const size_t transition,
		       std::vector<double> &scores) const;
*/
  std::string
  error_log() const;
  
  static const size_t FG_TO_BG_TRANSITION = 1;
  static const size_t BG_TO_FG_TRANSITION = 2;
    
private:

  double
  ViterbiDecoding(const std::vector<std::pair<double, double> > &values,
		  const std::vector<char> &contexts,
		  const std::vector<size_t> &reset_points,
		  double p_sf, double p_sb,
		  double p_ff, double p_fb, double p_ft,
		  double p_bf, double p_bb, double p_bt,
		  const std::vector<std::vector<betabin> > &distributions,
		  std::vector<bool> &ml_classes) const;

  double
  BaumWelchTraining(const std::vector<std::pair<double, double> > &values,
                    const std::vector<char> &contexts,
		    const std::vector<size_t> &reset_points,
		    double &p_sf, double &p_sb,
		    double &p_ff, double &p_fb, double &p_ft,
		    double &p_bf, double &p_bb, double &p_bt,
		    std::vector<std::vector<betabin> > &distributions) const;

  double
  PosteriorDecoding(const std::vector<std::pair<double, double> > &values,
		    const std::vector<char> &contexts,
		    const std::vector<size_t> &reset_points,
		    double p_sf, double p_sb,
		    double p_ff, double p_fb, double p_ft,
		    double p_bf, double p_bb, double p_bt,
		    const std::vector<std::vector<betabin> > &distributions,
		    std::vector<bool> &classes,
		    std::vector<double> &llr_scores) const;
/*
  void
  PosteriorScores(const std::vector<std::pair<double, double> > &values,
		  const std::vector<size_t> &reset_points,
		  double p_sf, double p_sb,
		  double p_ff, double p_fb, double p_ft,
		  double p_bf, double p_bb, double p_bt,
		  const betabin &fg_distro, const betabin &bg_distro,
		  const std::vector<bool> &classes,
		  std::vector<double> &llr_scores) const;

  void
  PosteriorScores(const std::vector<std::pair<double, double> > &values,
		  const std::vector<size_t> &reset_points,
		  double p_sf, double p_sb,
		  double p_ff, double p_fb, double p_ft,
		  double p_bf, double p_bb, double p_bt,
		  const double fg_alpha, const double fg_beta,
		  const double bg_alpha, const double bg_beta,
		  const bool class_id,
		  std::vector<double> &llr_scores) const;

  void
  TransitionPosteriors(const std::vector<std::pair<double, double> > &values,
		       const std::vector<size_t> &reset_points,
		       double p_sf, double p_sb,
		       double p_ff, double p_fb, double p_ft,
		       double p_bf, double p_bb, double p_bt,
		       const betabin &fg_distro, const betabin &bg_distro,
		       const size_t transition,
		       std::vector<double> &scores) const;
*/
  double
  single_iteration(const std::vector<std::pair<double, double> > &values,
		   const std::vector<char> &contexts,
		   const std::vector<std::vector<double> > &vals_a,
		   const std::vector<std::vector<double> > &vals_b,
		   const std::vector<size_t> &reset_points,
		   std::vector<std::pair<double, double> > &forward,
		   std::vector<std::pair<double, double> > &backward,
		   double &p_sf, double &p_sb,
		   double &p_ff, double &p_fb, double &p_ft,
		   double &p_bf, double &p_bb, double &p_bt,
		   std::vector<std::vector<betabin> > &distributions) const;

  double 
  forward_algorithm(const std::vector<std::pair<double, double> > &vals,
		    const std::vector<char> &contexts,
		    const size_t start, const size_t end,
		    const double lp_sf, const double lp_sb,
		    const double lp_ff, const double lp_fb, const double lp_ft,
		    const double lp_bf, const double lp_bb, const double lp_bt,
		    const std::vector<std::vector<betabin> > &distributions,
		    std::vector<std::pair<double, double> > &f) const;
  double 
  backward_algorithm(const std::vector<std::pair<double, double> > &vals,
		     const std::vector<char> &contexts,
		     const size_t start, const size_t end,
		     const double lp_sf, const double lp_sb,
		     const double lp_ff, const double lp_fb, const double lp_ft,
		     const double lp_bf, const double lp_bb, const double lp_bt,
		     const std::vector<std::vector<betabin> > &distributions,
		     std::vector<std::pair<double, double> > &b) const;

  double
  log_sum_log_vec(const std::vector<double> &vals, size_t limit) const;
  
  void 
  estimate_emissions(const std::vector<std::pair<double, double> > &f,
		     const std::vector<std::pair<double, double> > &b,
		     std::vector<std::vector<std::vector<double> > > &em_probs,
		     const std::vector<char> &contexts) const;
  
  void
  estimate_transitions(const std::vector<std::pair<double, double> > &vals,
		       const std::vector<char> &contexts,
		       const size_t start, const size_t end,
		       const std::vector<std::pair<double, double> > &f,
		       const std::vector<std::pair<double, double> > &b,
		       const double total,
		       const std::vector<std::vector<betabin> > &fg_distributions,
		       const double lp_ff, const double lp_fb,
		       const double lp_bf, const double lp_bb,
		       const double lp_ft, const double lp_bt,
		       std::vector<double> &ff_vals,
		       std::vector<double> &fb_vals,
		       std::vector<double> &bf_vals,
		       std::vector<double> &bb_vals) const;
  
  double
  log_sum_log(const double p, const double q) const;
  
  double MIN_PROB;
  double tolerance;
  size_t max_iterations;
  bool VERBOSE;
  bool DEBUG;
  
  mutable size_t emission_correction_count;
};

#endif
