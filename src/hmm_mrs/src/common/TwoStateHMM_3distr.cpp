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

#include "TwoStateHMM_3distr.hpp"

#include <iomanip>
#include <numeric>
#include <limits>
#include <cmath>

#include <gsl/gsl_sf_psi.h>
#include <gsl/gsl_sf_gamma.h>

// #pragma omp <rest of pragma>

using std::vector;
using std::pair;
using std::setw;
using std::max;
using std::cerr;
using std::endl;
using std::string;
using std::setprecision;


struct betabin {
  betabin(const double a, const double b) : 
    alpha(a), beta(b), lnbeta_helper(gsl_sf_lnbeta(a, b)) {}
  double operator()(const pair<double, double> &val) const;
  void fit(const vector<double> &mrates,
	   const vector<double> &invrates,
	   const vector<double> &p);
  string tostring() const;
  double alpha;
  double beta;
  double lnbeta_helper;
};

string
betabin::tostring() const {
  std::ostringstream os;
  os << setprecision(4) << alpha << " " << setprecision(4) << beta;
  return os.str();
}

double 
betabin::operator()(const pair<double, double> &val) const {
  const size_t x = static_cast<size_t>(val.first);
  const size_t n = static_cast<size_t>(x + val.second);
  //return gsl_sf_lnchoose(n, x) + 
  //  gsl_sf_lnbeta(alpha + x, beta + val.second) - lnbeta_helper;
  double ret = gsl_sf_lnchoose(n, x) + 
    gsl_sf_lnbeta(alpha + x, beta + val.second) - lnbeta_helper;
  return ret;
}

inline static double 
sign(double x) {
  return (x >= 0) ? 1.0 : -1.0;
}
static const double tolerance = 1e-10;

inline static double
invpsi(const double tolerance, const double x) {
//std::cout << "# Calculating invpsi with x=" << x << " ...";
  double L = 1.0, Y = std::exp(x);
  while (L > tolerance) {
    Y += L*sign(x - gsl_sf_psi(Y));
    L /= 2.0;
  }
//std::cout << " done\n";
  return Y;
}

static double
movement(const double curr, const double prev) {
  return std::abs(curr - prev)/std::max(std::fabs(curr), std::fabs(prev));
}

void
betabin::fit(const vector<double> &mrates, const vector<double> &invrates,
	     const vector<double> &p) {
//std::cout << "# Fitting Beta Distribution ...";
  const double p_total = std::accumulate(p.begin(), p.end(), 0.0);
  const double alpha_rhs = inner_product(mrates.begin(), mrates.end(), 
					 p.begin(), 0.0)/p_total;
  const double beta_rhs = inner_product(invrates.begin(), invrates.end(), 
					p.begin(), 0.0)/p_total;

/*double alpha_rhs;
double beta_rhs;
//for (vector<double>::const_iterator mrates_i=mrates.begin(); mrates_i<mrates.end(); ++mrates_i) {
std::cout << "Size of mrates: " << mrates.size() << std::endl;
vector<double>::const_iterator mrates_i=mrates.begin()+int(mrates.size()/4);
alpha_rhs = inner_product(mrates.begin(), mrates_i, p.begin(), 0.0)/p_total;
std::cout << (mrates_i-mrates.begin()) << " - " << alpha_rhs << std::endl;
mrates_i=mrates.begin()+int(mrates.size()/2);
alpha_rhs = inner_product(mrates.begin(), mrates_i, p.begin(), 0.0)/p_total;
std::cout << (mrates_i-mrates.begin()) << " - " << alpha_rhs << std::endl;
mrates_i=mrates.begin()+int(mrates.size()*3/4);
alpha_rhs = inner_product(mrates.begin(), mrates_i, p.begin(), 0.0)/p_total;
std::cout << (mrates_i-mrates.begin()) << " - " << alpha_rhs << std::endl;
mrates_i=mrates.begin()+int(mrates.size()*7/8);
alpha_rhs = inner_product(mrates.begin(), mrates_i, p.begin(), 0.0)/p_total;
std::cout << (mrates_i-mrates.begin()) << " - " << alpha_rhs << std::endl;
mrates_i=mrates.begin()+int(mrates.size()*15/16);
alpha_rhs = inner_product(mrates.begin(), mrates_i, p.begin(), 0.0)/p_total;
std::cout << (mrates_i-mrates.begin()) << " - " << alpha_rhs << std::endl;
mrates_i=mrates.begin()+int(mrates.size()-1);
alpha_rhs = inner_product(mrates.begin(), mrates_i, p.begin(), 0.0)/p_total;
std::cout << (mrates_i-mrates.begin()) << " - " << alpha_rhs << std::endl;
//}

std::cout << "inner products calculated: ptotal: " << p_total << " alpha_rhs: " << alpha_rhs << " beta_rhs: " << beta_rhs << "\n";
return;*/

  double prev_alpha = 0.0, prev_beta = 0.0;
  alpha = beta = 0.01;
  while (movement(alpha, prev_alpha) > tolerance &&
	 movement(beta, prev_beta) > tolerance) {
    prev_alpha = alpha;
    prev_beta = beta;
//std::cout << "pa " << prev_alpha << " pb " << prev_beta << "\n";
    alpha = invpsi(tolerance, gsl_sf_psi(prev_alpha + prev_beta) + alpha_rhs);
    beta = invpsi(tolerance, gsl_sf_psi(prev_alpha + prev_beta) + beta_rhs);
//std::cout << "a " << alpha << " b " << beta << "\n";
  }
  lnbeta_helper = gsl_sf_lnbeta(alpha, beta);
//std::cout << "Fitting Done\n";
}


////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////


inline double
TwoStateHMMB::log_sum_log(const double p, const double q) const {
  if (p == 0) {return q;}
  else if (q == 0) {return p;}
  const double larger = (p > q) ? p : q;
  const double smaller = (p > q) ? q : p;
  return larger + log(1.0 + exp(smaller - larger));
}


double
TwoStateHMMB::log_sum_log_vec(const vector<double> &vals, size_t limit) const {
  const vector<double>::const_iterator x = 
    std::max_element(vals.begin(), vals.begin() + limit);
  const double max_val = *x;
  const size_t max_idx = x - vals.begin();
  double sum = 1.0;
  for (size_t i = 0; i < limit; ++i) {
    if (i != max_idx) {
      sum += exp(vals[i] - max_val);
#ifdef DEBUG
      assert(finite(sum));
#endif
    }
  }
  return max_val + log(sum);
}



double
TwoStateHMMB::forward_algorithm(const vector<pair<double, double> > &vals,
				const vector<char> &contexts,
				const size_t start, const size_t end,
				const double lp_sf, const double lp_sb,
				const double lp_ff, const double lp_fb, 
				const double lp_ft,
				const double lp_bf, const double lp_bb, 
				const double lp_bt,
				const vector<vector<betabin> > &distributions,
				vector<pair<double, double> > &f) const {
  f[start].first = distributions[0][contexts[start]](vals[start]) + lp_sf;
  f[start].second = distributions[1][contexts[start]](vals[start]) + lp_sb;
  for (size_t i = start + 1; i < end; ++i) {
    // assert(finite(fg_distro.log_likelihood(vals[i])));
    const size_t k = i - 1;
    f[i].first = (distributions[0][contexts[i]](vals[i]) +
		  log_sum_log(f[k].first + lp_ff, f[k].second + lp_bf));
    f[i].second = (distributions[1][contexts[i]](vals[i]) + 
		   log_sum_log(f[k].first + lp_fb, f[k].second + lp_bb));
  }
  return log_sum_log(f[end - 1].first + lp_ft, f[end - 1].second + lp_bt);

}

double
TwoStateHMMB::backward_algorithm(const vector<pair<double, double> > &vals,
				const vector<char> &contexts,
				const size_t start, const size_t end,
				const double lp_sf, const double lp_sb,
				const double lp_ff, const double lp_fb, 
				const double lp_ft,
				const double lp_bf, const double lp_bb, 
				const double lp_bt,
				const vector<vector<betabin> > &distributions,
				vector<pair<double, double> > &b) const {
  b[end - 1].first = lp_ft;
  b[end - 1].second = lp_bt;
  for (size_t k = end - 1; k > start; --k) {
    size_t i = k - 1;
    const double fg_a = distributions[0][contexts[k]](vals[k]) + b[k].first;
    const double bg_a = distributions[1][contexts[k]](vals[k]) + b[k].second;
    b[i].first = log_sum_log(fg_a + lp_ff, bg_a + lp_fb);
    b[i].second = log_sum_log(fg_a + lp_bf, bg_a + lp_bb);
  }
  return log_sum_log(b[start].first + distributions[0][contexts[start]](vals[start]) + lp_sf,
		     b[start].second + distributions[1][contexts[start]](vals[start]) + lp_sb);
}


void
TwoStateHMMB::estimate_emissions(const vector<pair<double, double> > &f,
				const vector<pair<double, double> > &b,
				vector<vector<vector<double> > > &em_probs,
				const vector<char> &contexts) const {
//std::cout << "# Estimate Emissions\n";
  for (size_t i = 0; i < b.size(); ++i) {
    const double fg = (f[i].first + b[i].first);
    const double bg = (f[i].second + b[i].second);
    const double denom = log_sum_log(fg, bg);
    em_probs[0][contexts[i]].push_back(exp(fg - denom));
    em_probs[1][contexts[i]].push_back(exp(bg - denom));
  }
//std::cout << "  Done\n";
}



void
TwoStateHMMB::estimate_transitions(const vector<pair<double, double> > &vals,
				  const vector<char> &contexts,
				  const size_t start, const size_t end,
				  const vector<pair<double, double> > &f,
				  const vector<pair<double, double> > &b,
				  const double total,
				  const vector<vector<betabin> > &distributions,
				  const double lp_ff, const double lp_fb,
				  const double lp_bf, const double lp_bb,
				  const double lp_ft, const double lp_bt,
				  vector<double> &ff_vals,
				  vector<double> &fb_vals,
				  vector<double> &bf_vals,
				  vector<double> &bb_vals) const {
//std::cout << "# Estimate Transitions\n";
  for (size_t i = start + 1; i < end; ++i) {
    const size_t k = i - 1;
    const double b_first = b[i].first + distributions[0][contexts[i]](vals[i]) - total;
    const double b_second = b[i].second + distributions[1][contexts[i]](vals[i]) - total;
    
    const double ff = f[k].first;
    const double bb = f[k].second;
    
    ff_vals[k] = ff + lp_ff + b_first;
    fb_vals[k] = ff + lp_fb + b_second;
    
    bf_vals[k] = bb + lp_bf + b_first;
    bb_vals[k] = bb + lp_bb + b_second;
  }
}



double
TwoStateHMMB::single_iteration(const vector<pair<double, double> > &values,
			       const vector<char> &contexts,
			       const vector<vector<double> > &mrates,
			       const vector<vector<double> > &invrates,
			       const vector<size_t> &reset_points,
			       vector<pair<double, double> > &forward,
			       vector<pair<double, double> > &backward,
			       double &p_sf, double &p_sb,
			       double &p_ff, double &p_fb, double &p_ft,
			       double &p_bf, double &p_bb, double &p_bt,
			       vector<vector<betabin> > &distributions) const {
//if (DEBUG) std::cout << "in single iteration: " << mrates.size() << " " << mrates[0].size() << std::endl;
  vector<double> log_fg_expected;
  vector<double> log_bg_expected;
  
  double total_score = 0;
  
  const double lp_sf = log(p_sf);
  const double lp_sb = log(p_sb);
  const double lp_ff = log(p_ff);
  const double lp_fb = log(p_fb);
  const double lp_ft = log(p_ft);
  const double lp_bf = log(p_bf);
  const double lp_bb = log(p_bb);
  const double lp_bt = log(p_bt);
  
  assert(finite(lp_sf) && finite(lp_sb) && 
	 finite(lp_ff) && finite(lp_fb) && finite(lp_ft) && 
	 finite(lp_bf) && finite(lp_bb) && finite(lp_bt));

  // for estimating transitions
  vector<double> ff_vals(values.size(), 0);
  vector<double> fb_vals(values.size(), 0);
  vector<double> bf_vals(values.size(), 0);
  vector<double> bb_vals(values.size(), 0);
  
  
  // #pragma omp parallel  
  for (size_t i = 0; i < reset_points.size() - 1; ++i) {
//if (DEBUG) std::cout << "forward " << reset_points[i] << " " << reset_points[i+1] << std::endl;

    const double score = forward_algorithm(values, contexts,
					   reset_points[i], 
					   reset_points[i + 1],
					   lp_sf, lp_sb,
					   lp_ff, lp_fb, lp_ft,
					   lp_bf, lp_bb, lp_bt,
					   distributions, forward);
//if (DEBUG) std::cout << "backward " << reset_points[i] << " " << reset_points[i+1] << std::endl;
    const double backward_score = 
      backward_algorithm(values, contexts,
			 reset_points[i], 
			 reset_points[i + 1],
			 lp_sf, lp_sb,
			 lp_ff, lp_fb, lp_ft,
			 lp_bf, lp_bb, lp_bt,
			 distributions, backward);
    
    if (DEBUG && (fabs(score - backward_score)/
		  max(score, backward_score)) > 1e-10)
      cerr << "fabs(score - backward_score)/"
	   << "max(score, backward_score) > 1e-10" << endl;
    
    estimate_transitions(values, contexts,
			 reset_points[i], 
			 reset_points[i + 1],
			 forward, backward,
			 score, 
			 distributions,
			 lp_ff, lp_fb, lp_bf, 
			 lp_bb, lp_ft, lp_bt,
			 ff_vals, fb_vals,
			 bf_vals, bb_vals);

    total_score += score;
  }

  // Subtracting 1 from the limit of the summation because the final
  // term has no meaning since there is no transition to be counted
  // from the final observation (they all must go to terminal state)
  const double p_ff_new_estimate = exp(log_sum_log_vec(ff_vals, values.size() - 1));
  const double p_fb_new_estimate = exp(log_sum_log_vec(fb_vals, values.size() - 1));
  const double p_bf_new_estimate = exp(log_sum_log_vec(bf_vals, values.size() - 1));
  const double p_bb_new_estimate = exp(log_sum_log_vec(bb_vals, values.size() - 1));
  
  double denom = (p_ff_new_estimate + p_fb_new_estimate);
  p_ff = p_ff_new_estimate/denom - p_ft/2.0;
  p_fb = p_fb_new_estimate/denom - p_ft/2.0;
  
  if (p_ff < MIN_PROB) {
    if (DEBUG)
      cerr << "p_ff < MIN_PROB" << endl;
    p_ff = MIN_PROB;
  }
  
  if (p_fb < MIN_PROB) {
    if (DEBUG)
      cerr << "p_fb < MIN_PROB" << endl;
    p_fb = MIN_PROB;
  }
  
  denom = (p_bf_new_estimate + p_bb_new_estimate);
  p_bf = p_bf_new_estimate/denom - p_bt/2.0;
  p_bb = p_bb_new_estimate/denom - p_bt/2.0;
  
  if (p_bf < MIN_PROB) {
    if (DEBUG)
      cerr << "p_bf < MIN_PROB" << endl;
    p_bf = MIN_PROB;
  }

  if (p_bb < MIN_PROB) {
    if (DEBUG)
      cerr << "p_bb < MIN_PROB" << endl;
    p_bb = MIN_PROB;
  }

  p_sb = (p_bb + p_fb)/2.0;
  p_sf = (p_bf + p_ff)/2.0;
  
  // for estimating emissions
//  vector<double> fg_probs(values.size());
  vector<vector<vector<double> > > em_probs(mrates.size(), *(new vector<vector<double> >(distributions[0].size())));
//  for (size_t i = 0; i < 3; ++i) {
//    em_probs[0][i] = *(new vector<double>(mrates[i].size()));
//    em_probs[1][i] = *(new vector<double>(mrates[i].size()));
//  }

  estimate_emissions(forward, backward, em_probs, contexts);

//std::cout << "before fit: " << mrates.size() << " " << mrates[0].size() << std::endl;
//std::cout << "distr size " << distributions.size() << " " << distributions[0].size() << std::endl;
  for (size_t i = 0; i < distributions.size(); ++i) {
    for (size_t j = 0; j < distributions[0].size(); ++j) {
      if (mrates[j].size() == 0) {
        cerr << "Cannot train distributions for context " << j << " in state " << i << " (no sites of this context in that state!)\n";
      } else { 
        distributions[i][j].fit(mrates[j], invrates[j], em_probs[i][j]);
//      distributions[i][j].fit(mrates, invrates, bg_probs);
      }
    }
  }

  return total_score;
}



double
TwoStateHMMB::BaumWelchTraining(const vector<pair<double, double> > &values,
				const vector<char> &contexts,
				const vector<size_t> &reset_points,
				vector<double> &start_trans,
				vector<vector<double> > &trans, 
				vector<double> &end_trans,
				vector<vector<double> > &alphas,
				vector<vector<double> > &betas) const {
  
//  betabin fg_distro(fg_alpha, fg_beta);
//  betabin bg_distro(bg_alpha, bg_beta);

  assert(start_trans.size() >= 2);
  assert(end_trans.size() >= 2);
  assert(trans.size() >= 2);
  for (size_t i = 0; i < trans.size(); ++i)
    assert(trans[i].size() >= 2);
  assert(alphas.size() >= 2); assert(betas.size() >= 2);
  assert(alphas[0].size() >= 3); assert(betas[0].size() >= 3);
  assert(alphas[1].size() >= 3); assert(betas[1].size() >= 3);
  
  vector<vector<betabin> > distributions;
  for (size_t i = 0; i < alphas.size(); ++i) {
    vector<betabin> tmp;
    for (size_t c = 0; c < 3; ++c) {
      betabin distr(alphas[i][c], betas[i][c]);
      tmp.push_back(distr);
    }
    distributions.push_back(tmp);
  }

  const double score = BaumWelchTraining(values, contexts, reset_points,
					 start_trans[0], start_trans[1],
					 trans[0][0], trans[0][1], end_trans[0],
					 trans[1][0], trans[1][1], end_trans[1],
					 distributions);
  
  // Update alphas and betas:
  for (size_t i = 0; i < alphas.size(); ++i) {
    for (size_t j = 0; j < alphas[0].size(); ++j) {
      alphas[i][j] = distributions[i][j].alpha;
      betas[i][j]  = distributions[i][j].beta;
    }
  }
  //fg_alpha = fg_distro.alpha;
  //fg_beta = fg_distro.beta;

  //bg_alpha = bg_distro.alpha;
  //bg_beta = bg_distro.beta;

  return score;
}



double
TwoStateHMMB::BaumWelchTraining(const vector<pair<double, double> > &values,
				const vector<char> &contexts,
				const vector<size_t> &reset_points,
				double &p_sf, double &p_sb,
				double &p_ff, double &p_fb, double &p_ft,
				double &p_bf, double &p_bb, double &p_bt,
				vector<vector<betabin> > &distributions) const {
  
  vector<pair<double, double> > forward(values.size(), pair<double, double>(0, 0));
  vector<pair<double, double> > backward(values.size(), pair<double, double>(0, 0));
  
  if (VERBOSE) {
    cerr << setw(4)  << "ITR"
    	 << setw(10) << "p_01"
    	 << setw(10) << "p_10";
    for (size_t i = 0; i < distributions.size(); ++i) {
      cerr << setw(8) << i << setw(6) << " (CG)"
	   << setw(8) << i << setw(6) << " (CHG)"
	   << setw(8) << i << setw(6) << " (CHH)";
    }
    cerr << setw(14) << "SCORE"
	 << setw(14) << "PREV_SCORE"
    	 << setw(14) << "DELTA"
    	 << endl;

    cerr << setw(4) << "0"
	 << setw(10) << p_fb
	 << setw(10) << p_bf;
    for (size_t i = 0; i < distributions.size(); ++i)
      for (size_t j = 0; j < distributions[0].size(); ++j)
	cerr << setw(14) << distributions[i][j].tostring();
    cerr << endl;
  }
  
  double prev_total = -std::numeric_limits<double>::max();

  vector<vector<double> > mrates(distributions[0].size());
  vector<vector<double> > invrates(distributions[0].size());
  for (size_t i = 0; i < values.size(); ++i) {
    mrates[contexts[i]].push_back(
      log(std::min(std::max(values[i].first/(values[i].first + values[i].second), 1e-2), 1.0 - 1e-2))
    );
//std::cout << "v " << values[i].first << "-" << values[i].second << " c " << contexts[i] << " " << mrates[contexts[i]][i] << "\n";
    invrates[contexts[i]].push_back(
      log(1 - std::min(std::max(values[i].first/(values[i].first + values[i].second), 1e-2), 1.0 - 1e-2))
    );
  }
//std::cout << "after mrate creation: values: " << values.size() << " contexts: " << contexts.size() << " mrates: " << mrates.size() << " " << mrates[0].size() << std::endl;

  
  for (size_t i = 0; i < max_iterations; ++i) {
    
    double p_sf_est = p_sf;
    double p_sb_est = p_sb;
    double p_ff_est = p_ff;
    double p_fb_est = p_fb;
    double p_bf_est = p_bf;
    double p_bb_est = p_bb;
    double p_ft_est = p_ft;
    double p_bt_est = p_bt;

    double total = single_iteration(values, contexts,
				    mrates, invrates,
				    reset_points,
				    forward, backward,
				    p_sf_est, p_sb_est,
				    p_ff_est, p_fb_est, p_ft_est,
				    p_bf_est, p_bb_est, p_bt_est,
				    distributions);
    
    if (VERBOSE) {
      cerr << setw(4) << i + 1
	   << setw(10) << p_fb_est
	   << setw(10) << p_bf_est;
      for (size_t i = 0; i < distributions.size(); ++i)
	for (size_t j = 0; j < distributions[0].size(); ++j)
	  cerr << setw(14) << distributions[i][j].tostring();
      cerr << setw(14) << total
	   << setw(14) << prev_total
	   << setw(14) << (total - prev_total)/std::fabs(total)
	   << endl;
    }
    if ((total - prev_total) < tolerance) {// && (total - prev_total) >= 0) {
      if (VERBOSE)
	cerr << "CONVERGED" << endl << endl;
      break;
    }
    
    p_sf = p_sf_est;
    p_sb = p_sb_est;
    p_ff = p_ff_est;
    p_fb = p_fb_est;
    p_bf = p_bf_est;
    p_bb = p_bb_est;
    p_ft = p_ft_est;
    p_bt = p_bt_est;

    prev_total = total;
  }
  return prev_total;
}


/*
void
TwoStateHMMB::PosteriorScores(const vector<pair<double, double> > &values,
			      const vector<size_t> &reset_points,
			      const vector<double> &start_trans,
			      const vector<vector<double> > &trans, 
			      const vector<double> &end_trans,
			      const double fg_alpha, const double fg_beta,
			      const double bg_alpha, const double bg_beta,
			      const vector<bool> &classes,
			      vector<double> &llr_scores) const {

  const betabin fg_distro(fg_alpha, fg_beta);
  const betabin bg_distro(bg_alpha, bg_beta);

  assert(start_trans.size() >= 2);
  assert(end_trans.size() >= 2);
  assert(trans.size() >= 2);
  for (size_t i = 0; i < trans.size(); ++i)
    assert(trans[i].size() >= 2);
  
  return PosteriorScores(values, reset_points,
			 start_trans[0], start_trans[1],
			 trans[0][0], trans[0][1], end_trans[0],
			 trans[1][0], trans[1][1], end_trans[1],
			 fg_distro, bg_distro, classes, llr_scores);
}



void
TwoStateHMMB::PosteriorScores(const vector<pair<double, double> > &values,
			      const vector<size_t> &reset_points,
			      double p_sf, double p_sb,
			      double p_ff, double p_fb, double p_ft,
			      double p_bf, double p_bb, double p_bt,
			      const betabin &fg_distro,
			      const betabin &bg_distro,
			      const vector<bool> &classes,
			      vector<double> &llr_scores) const {

  double total_score = 0;
  
  const double lp_sf = log(p_sf);
  const double lp_sb = log(p_sb);
  const double lp_ff = log(p_ff);
  const double lp_fb = log(p_fb);
  const double lp_ft = log(p_ft);
  const double lp_bf = log(p_bf);
  const double lp_bb = log(p_bb);
  const double lp_bt = log(p_bt);
  
  assert(finite(lp_sf) && finite(lp_sb) && 
	 finite(lp_ff) && finite(lp_fb) && finite(lp_ft) && 
	 finite(lp_bf) && finite(lp_bb) && finite(lp_bt));

  vector<pair<double, double> > forward(values.size(), pair<double, double>(0, 0));
  vector<pair<double, double> > backward(values.size(), pair<double, double>(0, 0));

  for (size_t i = 0; i < reset_points.size() - 1; ++i) {
    const double score = forward_algorithm(values, 
					   reset_points[i],
					   reset_points[i + 1],
					   lp_sf, lp_sb,
					   lp_ff, lp_fb, lp_ft,
					   lp_bf, lp_bb, lp_bt,
					   fg_distro, bg_distro, forward);
    
    const double backward_score = 
      backward_algorithm(values,
			 reset_points[i],
			 reset_points[i + 1],
			 lp_sf, lp_sb,
			 lp_ff, lp_fb, lp_ft,
			 lp_bf, lp_bb, lp_bt,
			 fg_distro, bg_distro, backward);
    
    if (DEBUG && (fabs(score - backward_score)/
		  max(score, backward_score)) > 1e-10)
      cerr << "fabs(score - backward_score)/"
	   << "max(score, backward_score) > 1e-10" << endl;

    total_score += score;
  }
  
  llr_scores.resize(values.size());
  for (size_t i = 0; i < values.size(); ++i) {
    const double fg_state = forward[i].first + backward[i].first;
    const double bg_state = forward[i].second + backward[i].second;
    if (classes[i])
      llr_scores[i] = (fg_state - bg_state);
    else 
      llr_scores[i] = (bg_state - fg_state);
  }
}



void
TwoStateHMMB::PosteriorScores(const vector<pair<double, double> > &values,
			      const vector<size_t> &reset_points,
			      const vector<double> &start_trans,
			      const vector<vector<double> > &trans, 
			      const vector<double> &end_trans,
			      const double fg_alpha, const double fg_beta,
			      const double bg_alpha, const double bg_beta,
			      const bool fg_class,
			      vector<double> &llr_scores) const {
  
  const betabin fg_distro(fg_alpha, fg_beta);
  const betabin bg_distro(bg_alpha, bg_beta);


  assert(start_trans.size() >= 2);
  assert(end_trans.size() >= 2);
  assert(trans.size() >= 2);
  for (size_t i = 0; i < trans.size(); ++i)
    assert(trans[i].size() >= 2);
  
  return PosteriorScores(values, reset_points,
			 start_trans[0], start_trans[1],
			 trans[0][0], trans[0][1], end_trans[0],
			 trans[1][0], trans[1][1], end_trans[1],
			 fg_distro, bg_distro, fg_class, llr_scores);
}




void
TwoStateHMMB::PosteriorScores(const vector<pair<double, double> > &values,
			     const vector<size_t> &reset_points,
			     double p_sf, double p_sb,
			     double p_ff, double p_fb, double p_ft,
			     double p_bf, double p_bb, double p_bt,
			     const betabin &fg_distro,
			     const betabin &bg_distro,
			     const bool fg_class,
			     vector<double> &llr_scores) const {
  
  double total_score = 0;
  
  const double lp_sf = log(p_sf);
  const double lp_sb = log(p_sb);
  const double lp_ff = log(p_ff);
  const double lp_fb = log(p_fb);
  const double lp_ft = log(p_ft);
  const double lp_bf = log(p_bf);
  const double lp_bb = log(p_bb);
  const double lp_bt = log(p_bt);
  
  assert(finite(lp_sf) && finite(lp_sb) && 
	 finite(lp_ff) && finite(lp_fb) && finite(lp_ft) && 
	 finite(lp_bf) && finite(lp_bb) && finite(lp_bt));
  
  vector<pair<double, double> > forward(values.size(), pair<double, double>(0, 0));
  vector<pair<double, double> > backward(values.size(), pair<double, double>(0, 0));
  
  for (size_t i = 0; i < reset_points.size() - 1; ++i) {
    const double score = 
      forward_algorithm(values, reset_points[i], reset_points[i + 1],
			lp_sf, lp_sb, lp_ff, lp_fb, lp_ft, lp_bf, lp_bb, lp_bt,
			fg_distro, bg_distro, forward);
    
    const double backward_score = 
      backward_algorithm(values, reset_points[i], reset_points[i + 1],
			 lp_sf, lp_sb, lp_ff, lp_fb, lp_ft,
			 lp_bf, lp_bb, lp_bt,
			 fg_distro, bg_distro, backward);
    
    if (DEBUG && (fabs(score - backward_score)/
		  max(score, backward_score)) > 1e-10)
      cerr << "fabs(score - backward_score)/"
	   << "max(score, backward_score) > 1e-10" << endl;
    total_score += score;
  }
  
  llr_scores.resize(values.size());
  for (size_t i = 0; i < values.size(); ++i) {
    const double fg_state = forward[i].first + backward[i].first;
    const double bg_state = forward[i].second + backward[i].second;
    if (fg_class)
      llr_scores[i] = exp(fg_state - log_sum_log(fg_state, bg_state));
    else
      llr_scores[i] = exp(bg_state - log_sum_log(fg_state, bg_state));
    //     if (fg_class)
    //       llr_scores[i] = (fg_state - bg_state);
    //     else
    //       llr_scores[i] = (bg_state - fg_state);
  }
}





void
TwoStateHMMB::TransitionPosteriors(const vector<pair<double, double> > &values,
				   const vector<size_t> &reset_points,
				   const vector<double> &start_trans,
				   const vector<vector<double> > &trans, 
				   const vector<double> &end_trans,
				   const double fg_alpha, const double fg_beta,
				   const double bg_alpha, const double bg_beta,
				   const size_t transition,
				   vector<double> &llr_scores) const {
  
  const betabin fg_distro(fg_alpha, fg_beta);
  const betabin bg_distro(bg_alpha, bg_beta);

  assert(start_trans.size() >= 2);
  assert(end_trans.size() >= 2);
  assert(trans.size() >= 2);
  for (size_t i = 0; i < trans.size(); ++i)
    assert(trans[i].size() >= 2);
  
  return TransitionPosteriors(values, reset_points,
			      start_trans[0], start_trans[1],
			      trans[0][0], trans[0][1], end_trans[0],
			      trans[1][0], trans[1][1], end_trans[1],
			      fg_distro, bg_distro, transition, llr_scores);
}


void
TwoStateHMMB::TransitionPosteriors(const vector<pair<double, double> > &values,
				   const vector<size_t> &reset_points,
				   double p_sf, double p_sb,
				   double p_ff, double p_fb, double p_ft,
				   double p_bf, double p_bb, double p_bt,
				   const betabin &fg_distro,
				   const betabin &bg_distro,
				   const size_t transition,
				   vector<double> &scores) const {
  

  double total_score = 0;
  
  const double lp_sf = log(p_sf);
  const double lp_sb = log(p_sb);
  const double lp_ff = log(p_ff);
  const double lp_fb = log(p_fb);
  const double lp_ft = log(p_ft);
  const double lp_bf = log(p_bf);
  const double lp_bb = log(p_bb);
  const double lp_bt = log(p_bt);
  
  assert(finite(lp_sf) && finite(lp_sb) && 
	 finite(lp_ff) && finite(lp_fb) && finite(lp_ft) && 
	 finite(lp_bf) && finite(lp_bb) && finite(lp_bt));
  
  vector<pair<double, double> > forward(values.size(), pair<double, double>(0, 0));
  vector<pair<double, double> > backward(values.size(), pair<double, double>(0, 0));
  
  for (size_t i = 0; i < reset_points.size() - 1; ++i) {
    const double score = forward_algorithm(values, 
					   reset_points[i],
					   reset_points[i + 1],
					   lp_sf, lp_sb,
					   lp_ff, lp_fb, lp_ft,
					   lp_bf, lp_bb, lp_bt,
					   fg_distro, bg_distro, forward);
    
    const double backward_score = 
      backward_algorithm(values, reset_points[i], reset_points[i + 1],
			 lp_sf, lp_sb, lp_ff, lp_fb, lp_ft,
			 lp_bf, lp_bb, lp_bt,
			 fg_distro, bg_distro, backward);
    
    if (DEBUG && (fabs(score - backward_score)/
		  max(score, backward_score)) > 1e-10)
      cerr << "fabs(score - backward_score)/"
	   << "max(score, backward_score) > 1e-10" << endl;

    total_score += score;
  }
  
  scores.resize(values.size());
  size_t j = 0;
  for (size_t i = 0; i < values.size(); ++i) {
    if (i == reset_points[j]) {
      ++j;
      scores[i] = 0;
    }
    else {
      const double fg_to_fg_state = forward[i - 1].first + lp_ff + // transition
	// emission for value i + 1
	fg_distro(values[i]) + backward[i].first;
      const double fg_to_bg_state = forward[i - 1].first + lp_fb + 
	bg_distro(values[i]) + backward[i].second;
      const double bg_to_fg_state = forward[i - 1].second + lp_bf + 
	fg_distro(values[i]) + backward[i].first;
      const double bg_to_bg_state = forward[i - 1].second + lp_bb + 
	bg_distro(values[i]) + backward[i].second;
      const double denom = log_sum_log(log_sum_log(fg_to_fg_state, fg_to_bg_state),
				       log_sum_log(bg_to_fg_state, bg_to_bg_state));
      double numerator = fg_to_fg_state;
      if (transition == 1)
	numerator = fg_to_bg_state;
      if (transition == 2)
	numerator = bg_to_fg_state;
      if (transition == 3)
	numerator = bg_to_bg_state;
      scores[i] = exp(numerator - denom);
    }
  }
}
*/

double
TwoStateHMMB::PosteriorDecoding(const vector<pair<double, double> > &values,
				const vector<char> &contexts,
				const vector<size_t> &reset_points,
				const vector<double> &start_trans,
				const vector<vector<double> > &trans, 
				const vector<double> &end_trans,
				vector<vector<double> > &alphas,
				vector<vector<double> > &betas,
				vector<bool> &classes,
				vector<double> &llr_scores) const {

//  const betabin fg_distro(fg_alpha, fg_beta);
//  const betabin bg_distro(bg_alpha, bg_beta);
  vector<vector<betabin> > distributions;
  for (size_t i = 0; i < alphas.size(); ++i) {
    vector<betabin> tmp;
    for (size_t c = 0; c < 3; ++c) {
      betabin distr(alphas[i][c], betas[i][c]);
      tmp.push_back(distr);
    }
    distributions.push_back(tmp);
  }



  assert(start_trans.size() >= 2);
  assert(end_trans.size() >= 2);
  assert(trans.size() >= 2);
  for (size_t i = 0; i < trans.size(); ++i)
    assert(trans[i].size() >= 2);
  
  return PosteriorDecoding(values, contexts, reset_points,
			   start_trans[0], start_trans[1],
			   trans[0][0], trans[0][1], end_trans[0],
			   trans[1][0], trans[1][1], end_trans[1],
			   distributions, classes, llr_scores);
}


double
TwoStateHMMB::PosteriorDecoding(const vector<pair<double, double> > &values,
			       const vector<char> &contexts,
			       const vector<size_t> &reset_points,
			       double p_sf, double p_sb,
			       double p_ff, double p_fb, double p_ft,
			       double p_bf, double p_bb, double p_bt,
			       const vector<vector<betabin> > &distributions,
			       vector<bool> &classes,
			       vector<double> &llr_scores) const {
  
  double total_score = 0;
  
  const double lp_sf = log(p_sf);
  const double lp_sb = log(p_sb);
  const double lp_ff = log(p_ff);
  const double lp_fb = log(p_fb);
  const double lp_ft = log(p_ft);
  const double lp_bf = log(p_bf);
  const double lp_bb = log(p_bb);
  const double lp_bt = log(p_bt);
  
  assert(finite(lp_sf) && finite(lp_sb) && 
	 finite(lp_ff) && finite(lp_fb) && finite(lp_ft) && 
	 finite(lp_bf) && finite(lp_bb) && finite(lp_bt));

  vector<pair<double, double> > forward(values.size(), pair<double, double>(0, 0));
  vector<pair<double, double> > backward(values.size(), pair<double, double>(0, 0));

  for (size_t i = 0; i < reset_points.size() - 1; ++i) {
    const double score = forward_algorithm(values, contexts,
					   reset_points[i],
					   reset_points[i + 1],
					   lp_sf, lp_sb,
					   lp_ff, lp_fb, lp_ft,
					   lp_bf, lp_bb, lp_bt,
					   distributions, forward);
    
    const double backward_score = 
      backward_algorithm(values, contexts,
			 reset_points[i],
			 reset_points[i + 1],
			 lp_sf, lp_sb,
			 lp_ff, lp_fb, lp_ft,
			 lp_bf, lp_bb, lp_bt,
			 distributions, backward);
    
    if (DEBUG && (fabs(score - backward_score)/
		  max(score, backward_score)) > 1e-10)
      cerr << "fabs(score - backward_score)/"
	   << "max(score, backward_score) > 1e-10" << endl;

    total_score += score;
  }
  
  classes.resize(values.size());
//   llr_scores.resize(values.size());
//   for (size_t i = 0; i < values.size(); ++i) {
//     const double fg_state = forward[i].first + backward[i].first;
//     const double bg_state = forward[i].second + backward[i].second;
//     classes[i] = static_cast<bool>(fg_state > bg_state);
//     llr_scores[i] = (fg_state - bg_state);
    

  llr_scores.resize(values.size());
  for (size_t i = 0; i < values.size(); ++i) {
    const double fg_state = forward[i].first + backward[i].first;
    const double bg_state = forward[i].second + backward[i].second;
    
    // const double bg_state = forward[i].second + backward[i].second;
    classes[i] = static_cast<bool>(fg_state > bg_state);
//std::cout << classes[i];
    
    // if (fg_class)
    llr_scores[i] = exp(fg_state - log_sum_log(fg_state, bg_state));
    //     else
    //       llr_scores[i] = exp(bg_state - log_sum_log(fg_state, bg_state));
    //     if (fg_class)
    //       llr_scores[i] = (fg_state - bg_state);
    //     else
    //       llr_scores[i] = (bg_state - fg_state);
  }
//std::cout << std::endl;
  
  
  return total_score;
}


/*************************************************************
 *
 * Functions for Viterbi training and decoding.
 *
 *************************************************************/


double
TwoStateHMMB::ViterbiDecoding(const vector<pair<double, double> > &values,
			      const vector<char> &contexts,
			      const vector<size_t> &reset_points,
			      const vector<double> &start_trans,
			      const vector<vector<double> > &trans, 
			      const vector<double> &end_trans,
			      vector<vector<double> > &alphas,
			      vector<vector<double> > &betas,
			      vector<bool> &classes) const {
  
  vector<vector<betabin> > distributions;
  for (size_t i = 0; i < alphas.size(); ++i) {
    vector<betabin> tmp;
    for (size_t c = 0; c < 3; ++c) {
      betabin distr(alphas[i][c], betas[i][c]);
      tmp.push_back(distr);
    }
    distributions.push_back(tmp);
  }
  
  assert(start_trans.size() >= 2);
  assert(end_trans.size() >= 2);
  assert(trans.size() >= 2);
  for (size_t i = 0; i < trans.size(); ++i)
    assert(trans[i].size() >= 2);
  
  return ViterbiDecoding(values, contexts, reset_points,
			 start_trans[0], start_trans[1],
			 trans[0][0], trans[0][1], end_trans[0],
			 trans[1][0], trans[1][1], end_trans[1],
			 distributions, classes);
}


double
TwoStateHMMB::ViterbiDecoding(const vector<pair<double, double> > &values,
			     const vector<char> &contexts,
			     const vector<size_t> &reset_points,
			     double p_sf, double p_sb,
			     double p_ff, double p_fb, double p_ft,
			     double p_bf, double p_bb, double p_bt,
			     const vector<vector<betabin> > &distributions,
			     vector<bool> &ml_classes) const {
  const double lp_sf = log(p_sf);
  const double lp_sb = log(p_sb);
  const double lp_ff = log(p_ff);
  const double lp_fb = log(p_fb);
  const double lp_ft = log(p_ft);
  const double lp_bf = log(p_bf);
  const double lp_bb = log(p_bb);
  const double lp_bt = log(p_bt);
  
  // ml_classes = vector<bool>(values.size());
  double total = 0;
  for (size_t i = 0; i < reset_points.size() - 1; ++i) {

    const size_t lim = reset_points[i + 1] - reset_points[i];
    
    vector<pair<double, double> > v(lim, pair<double, double>(0, 0));
    vector<pair<size_t, size_t> > trace(lim, pair<size_t, size_t>(0, 0));
    
    v.front().first = distributions[0][contexts[reset_points[i]]](values[reset_points[i]]) + lp_sf;
    v.front().second = distributions[1][contexts[reset_points[i]]](values[reset_points[i]]) + lp_sb;

    for (size_t j = 1; j < lim; ++j) {
      
      const double ff = v[j - 1].first + lp_ff;
      const double bf = v[j - 1].second + lp_bf;
      const double fg_log_emmit =
	distributions[0][contexts[reset_points[i]+j]](values[reset_points[i] + j]);
      if (ff > bf) {
	v[j].first = fg_log_emmit + ff;
	trace[j].first = 0;
      }
      else {
	v[j].first = fg_log_emmit + bf;
	trace[j].first = 1;
      }
    
      const double fb = v[j - 1].first + lp_fb;
      const double bb = v[j - 1].second + lp_bb;
      const double bg_log_emmit = 
	distributions[1][contexts[reset_points[i]+j]](values[reset_points[i] + j]);
      if (fb > bb) {
	v[j].second = bg_log_emmit + fb;
	trace[j].second = 0;
      }
      else {
	v[j].second = bg_log_emmit + bb;
	trace[j].second = 1;
      }
    }
    v.back().first += lp_ft;
    v.back().second += lp_bt;
    
    vector<bool> inner_ml_classes;
    
    // do the traceback
    size_t prev = 0;
    if (v.back().first > v.back().second) {
      inner_ml_classes.push_back(true);
      prev = trace.back().first;
    }
    else {
      inner_ml_classes.push_back(false);
      prev = trace.back().second;
    }
    
    for (size_t j = trace.size() - 1; j > 0; --j) {
      const size_t k = j - 1;
      if (prev == 0) {
	inner_ml_classes.push_back(true);
	prev = trace[k].first;
      }
      else {
	inner_ml_classes.push_back(false);
	prev = trace[k].second;
      }
    }
    reverse(inner_ml_classes.begin(), inner_ml_classes.end());
    ml_classes.insert(ml_classes.end(), inner_ml_classes.begin(), 
		      inner_ml_classes.end());
    
    total += max(v.back().first, v.back().second);
  }
  return total;
}
