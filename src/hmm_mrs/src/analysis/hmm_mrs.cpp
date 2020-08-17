/* Copyright (C) 2009-2012 University of Southern California
 *                         Andrew D Smith
 * Author: Andrew D. Smith
 * modified by Joerg Hagmann, Max Planck Institute, Tuebingen
 *
 * This is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this software; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA
 * 02110-1301 USA
 */

#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"
#include "GenomicRegion.hpp"
#include "OptionParser.hpp"
#include "TwoStateHMM_3distr.hpp"

#include <numeric>
#include <cmath>
#include <fstream>
#include <iomanip>

#include <unistd.h>

using std::string;
using std::vector;
using std::cout;
using std::endl;
using std::cerr;
using std::numeric_limits;
using std::max;
using std::min;
using std::pair;

double
get_fdr_cutoff(const vector<double> &scores, const double fdr) {
  if (fdr <= 0)
    return numeric_limits<double>::max();
  else if (fdr > 1)
    return numeric_limits<double>::min();
  vector<double> local(scores);
  std::sort(local.begin(), local.end());
  size_t i = 0;
  for (; i < local.size() - 1 &&
	 local[i+1] < fdr*static_cast<double>(i+1)/local.size(); ++i);
  return local[i];
}


/*static void
get_domain_scores(const vector<bool> &classes,
		  const vector<pair<double, double> > &meth,
		  const vector<size_t> &reset_points,
		  vector<double> &scores) {
  static const bool CLASS_ID = true;
  size_t n_cpgs = 0, reset_idx = 1;
  bool in_domain = false;
  double score = 0;
  for (size_t i = 0; i < classes.size(); ++i) {
    if (reset_points[reset_idx] == i) {
      if (in_domain) {
	in_domain = false;
	scores.push_back(score);
	score = 0;
      }
      ++reset_idx;
    }
    if (classes[i] == CLASS_ID) {
      in_domain = true;
      score += 1.0 - (meth[i].first/(meth[i].first + meth[i].second));
      ++n_cpgs;
    }
    else if (in_domain) {
      in_domain = false;
      scores.push_back(score);
      score = 0;
    }
  }
}*/

size_t
trim_end_of_region(int region_end, int region_len,
		   const vector<SimpleGenomicRegion> &cpgs,
		   const vector<pair<double, double> > &meth,
		   double trimrate, vector<GenomicRegion> &domains,
		   const bool human) {
  double methrate;
  for (int j=region_end; j>=region_end-region_len; --j) {
    methrate = meth[j].second/(meth[j].first + meth[j].second);
	if (human == true) methrate = 1 - methrate;
    if (methrate >= trimrate) {
      domains.back().set_end(cpgs[j].get_start());
      //domains.back().pop_site();
      return (region_end-j);
    }
  }
  return region_len;
}

static void
build_and_score_domains(const bool VERBOSE,
			const bool HUMAN,
			const vector<SimpleGenomicRegion> &cpgs,
			const vector<pair<double, double> > &meth,
			const vector<size_t> &reset_points,
			const vector<bool> &classes,
			double trimrate, int merge_regions,
			double lower_bound_rate,
			double upper_bound_rate,
			vector<GenomicRegion> &domains,
			vector<double> &scores) {

// IMPORTANT:
// - domains.get_score() are just the number of sites within a domain
// - scores are the sum of meth.rates within a domain, and this is used for
//   subsequent p-value and FDR calculations!

  static const bool CLASS_ID = true;
  size_t n_domains = 0, reset_idx = 1, prev_end_idx = 0, prev_start_idx = 0;
  bool in_domain = false, first_methyl = false;
  double score_in = 0, score_out = 0, methrate = 0;
  GenomicRegion last_domain;
  int FP = 0, FN = 0, FP_trimmed = 0, FN_trimmed = 0;
  for (size_t i = 0; i < classes.size(); ++i) {
//std::cout<<cpgs[i].get_start()<<"\t"<<(meth[i].second/(meth[i].first + meth[i].second))<<"\t"<<classes[i]<<std::endl;
    if (reset_points[reset_idx] == i || classes[i] != CLASS_ID) {
      if (in_domain) {
        // END OF DOMAIN
	in_domain = false;
	if (domains.back().get_score() > 0) { // trimming function sets region end
	  domains.back().set_score(
			domains.back().get_score() -
			trim_end_of_region(prev_end_idx, domains.back().get_score(),
					   cpgs, meth, trimrate, domains, HUMAN)
	  );
	  prev_end_idx = prev_start_idx + domains.back().get_score() - 1;
	}
	// n_cpgs can now be 0 after trimming:
	if (domains.back().get_score() == 0) domains.pop_back(); //
	else {
	  scores.push_back(score_in);
	  // calculate avg and .75 percentile coverage:
	  domains.back().calc_average_cov(meth, prev_start_idx, domains.back().get_score());
	  domains.back().calc_percentile_cov(meth, prev_start_idx, domains.back().get_score(), 0.75);
	}
	//n_cpgs = 0;
	first_methyl = false;
	score_in = 0;
      }

      if (classes[i] != CLASS_ID) {
	// NOT IN DOMAIN
	methrate = meth[i].second/(meth[i].first + meth[i].second);
	if (HUMAN == true) methrate = 1 - methrate;
	score_out += methrate;
      }
      if (reset_points[reset_idx] == i && reset_idx < reset_points.size()) ++reset_idx;

      // Store quality criterion
      if (methrate > upper_bound_rate) FN++;
    }

    if (classes[i] == CLASS_ID) {
      // IN METHYLATED REGION
      if (!in_domain) {
        // POTENTIAL START OF DOMAIN
	in_domain = true;
	if (!domains.empty()) last_domain = domains.back();
	domains.push_back(GenomicRegion(cpgs[i]));
	domains.back().set_name(toa(n_domains++));
      }
      methrate = meth[i].second/(meth[i].first + meth[i].second);
	  if (HUMAN == true) methrate = 1 - methrate;
      score_in += methrate;

      // Trim low methylated sites at the front of a region
      if (first_methyl != CLASS_ID) {
	if (methrate >= trimrate) {
	  // Merge this region with the previous if they lie
	  // less than 'merge_regions' bases apart
	  if (n_domains>1 && domains.back().same_chrom(last_domain) &&
	      abs(int(cpgs[i].get_start()) - int(last_domain.get_end()))<=merge_regions) {
	    // MERGE WITH PREVIOUS
	    scores.pop_back();
	    domains.pop_back();
	    domains.back().set_score(domains.back().get_score()+i-prev_end_idx-1);
	    //for (int pos=int(prev_end_idx+1); pos<int(i); pos++)
		//domains.back().add_site(meth[pos].first + meth[pos].second);
	    score_in = scores.back() + score_out;
	    //n_cpgs += i - prev_end_idx;
	  }
	  else {
	    // REAL START OF DOMAIN (no merging, no more front trimming)
	    domains.back().set_start(cpgs[i].get_start());
	    domains.back().set_end(cpgs[i].get_start());
	    //domains.back().add_site(meth[i].first + meth[i].second);
	    prev_start_idx = i;
	  }
	  first_methyl = true;
	  score_out = 0;
	}
	else {
	  if (methrate > upper_bound_rate) FN_trimmed++;
	  score_out += methrate;
	}
      }

      if (first_methyl == CLASS_ID) {
	domains.back().set_score(domains.back().get_score()+1); //
	//domains.back().add_site(meth[i].first + meth[i].second);
	prev_end_idx = i;
	if (methrate < lower_bound_rate) FP_trimmed++;
      }

      // Store quality criterion
      if (methrate < lower_bound_rate) FP++;
    }
  }

  if (upper_bound_rate >= 0) {
    std::cout << std::endl << "*** QUALITY CRITERIA ***\n";
    std::cout << "  Low methylated positions (< " << int(lower_bound_rate*100) << "\%) in DMRs     : " << FP;
    if (trimrate > 0) { std::cout << " (after trimming: " << FP_trimmed << ")"; }
    std::cout << "\n  High methylated positions (> " << int(upper_bound_rate*100) << "\%) not in DMRs: " << FN;
    if (trimrate > 0) { std::cout << " (after trimming: " << FN+FN_trimmed << ")"; }
    std::cout << std::endl << std::endl;
  }
}

template <class T, class U> static void
separate_regions(const bool VERBOSE, const size_t desert_size,
		 vector<SimpleGenomicRegion> &cpgs,
		 vector<T> &meth, vector<U> &reads,
		 vector<size_t> &reset_points) {
  if (VERBOSE)
    cerr << "[SEPARATING BY CPG DESERT]" << endl;
/*  // eliminate the zero-read cpgs
  size_t j = 0;
  for (size_t i = 0; i < cpgs.size(); ++i)
    if (reads[i] > 0) {
      cpgs[j] = cpgs[i];
      meth[j] = meth[i];
      reads[j] = reads[i];
      ++j;
    }
  cpgs.erase(cpgs.begin() + j, cpgs.end());
  meth.erase(meth.begin() + j, meth.end());
  reads.erase(reads.begin() + j, reads.end());
*/
  // segregate cpgs
  std::tr1::unordered_map<std::string, bool> hash;
  size_t prev_cpg = 0;
  for (size_t i = 0; i < cpgs.size(); ++i) {
    if (i>0 && cpgs[i].same_chrom(cpgs[i-1]) == false) {
      if (hash[cpgs[i].get_chrom()] == false)
        hash.insert( std::pair<std::string, bool>(cpgs[i].get_chrom(), 1) );
      else {
        std::ostringstream coords;
        coords << cpgs[i].get_chrom() << ":" << cpgs[i].get_start();
        throw SMITHLABException("CpGs not sorted! "+coords.str());
      }
    }
    else if (i>0 && cpgs[i].same_chrom(cpgs[i-1]) && cpgs[i-1].get_start() >= cpgs[i].get_start()) {
      std::ostringstream coords;
      coords << cpgs[i].get_chrom() << ":" << cpgs[i].get_start();
      throw SMITHLABException("CpGs not sorted! "+coords.str());
    }

//    if (i>0 && cpgs[i].same_chrom(cpgs[i - 1]) == false)
//      reset_points.push_back(0); // 0 is the marker for chromosome/scaffold breakpoint

    const size_t dist = (i > 0 && cpgs[i].same_chrom(cpgs[i - 1])) ?
      cpgs[i].get_start() - prev_cpg : numeric_limits<size_t>::max();
    if (dist > desert_size)
      reset_points.push_back(i);
    prev_cpg = cpgs[i].get_start();
  }
  reset_points.push_back(cpgs.size()-1);
  if (VERBOSE)
    cerr << "CPGS RETAINED: " << cpgs.size() << endl
	 << "DESERTS REMOVED: " << reset_points.size() - 2 << endl << endl;
}


/*
 * FUNCTION TO "FOLD" THE METHYLATION PROFILE SO THAT THE MIDDLE
 * METHYLATION BECOMES LOWER METHYLATION, AND BOTH THE LOW AND HIGH
 * METHYLATION BECOME HIGH. THIS METHOD ACTUALLY SEEMS TO WORK.
 */
static void
make_partial_meth(vector<GenomicRegion> &cpgs) {
  for (size_t i = 0; i < cpgs.size(); ++i) {
    double cpg_meth = cpgs[i].get_score();
    if (cpg_meth > 0.5)
      cpg_meth = 1.0 - cpg_meth;
    cpgs[i].set_score(1.0 - (2*cpg_meth));
  }
}


static void
load_cpgs(const bool VERBOSE, const bool PARTIAL_METH,
	  string cpgs_file, vector<SimpleGenomicRegion> &cpgs,
	  vector<pair<double, double> > &meth,
	  vector<char> &contexts, vector<size_t> &reads,
          int min_cov, size_t matrix_idx) {
  vector<GenomicRegion> cpgs_in;

  if (matrix_idx == 0)
    ReadBEDFile(cpgs_file, cpgs_in);
  else {
    //if (!checkGenomeMatrixFormat(cpgs_file))
    //  throw SMITHLABException("Wrong genome matrix file format: "+cpgs_file);
    ReadGenomeMatrixFile(cpgs_file, min_cov, matrix_idx, cpgs_in);
  }

  //if (!check_sorted(cpgs_in))
    //throw SMITHLABException("CpGs not sorted in file \"" + cpgs_file + "\"");

  if (PARTIAL_METH)
    make_partial_meth(cpgs_in);

  int cov;
  string r;
  for (size_t i = 0; i < cpgs_in.size(); ++i) {
    r = cpgs_in[i].get_name();

    cov = atoi(r.substr(r.find_first_of(":") + 1).c_str());

    if (cov >= min_cov) {
      cpgs.push_back(SimpleGenomicRegion(cpgs_in[i]));
      if (cpgs_in[i].get_context() < 1)
        throw SMITHLABException("Invalid context in line: " + cpgs_in[i].tostring());
      contexts.push_back(cpgs_in[i].get_context()-1);
      meth.push_back(std::make_pair(cpgs_in[i].get_score(), 0.0));

      reads.push_back(cov);

      meth.back().first = floor(meth.back().first * cov + 0.5);
      meth.back().second = floor(cov - meth.back().first + 0.5);
    }
  }
  if (VERBOSE)
    cerr << "TOTAL CPGS: " << cpgs.size() << " " << meth.size() << endl
	 << "MEAN COVERAGE: "
	 << accumulate(reads.begin(), reads.end(), 0.0)/reads.size()
	 << endl << endl;
}



static void
shuffle_cpgs(const TwoStateHMMB &hmm,
		 const bool human,
	     vector<pair<double, double> > meth,
	     const vector<char> &contexts,
	     vector<size_t> &reset_points,
	     const vector<double> &start_trans,
	     const vector<vector<double> > &trans,
	     const vector<double> &end_trans,
	     vector<vector<double> > &alphas,
	     vector<vector<double> > &betas,
	     vector<double> &domain_scores,
	     const vector<SimpleGenomicRegion> &cpgs,
	     double trimrate, int merge_regions, bool viterbi) {
  srand(time(0) + getpid());
  random_shuffle(meth.begin(), meth.end());
  vector<bool> classes;
  vector<double> scores;
  if (viterbi == true) {
    hmm.ViterbiDecoding(meth, contexts, reset_points, start_trans, trans,
			end_trans, alphas, betas, classes);
  }
  else {
    hmm.PosteriorDecoding(meth, contexts, reset_points, start_trans, trans,
			  end_trans, alphas, betas, classes, scores);
  }
  //get_domain_scores(classes, meth, reset_points, domain_scores);
  vector<GenomicRegion> domains;
  build_and_score_domains(false, human, cpgs, meth, reset_points, classes,
			  trimrate, merge_regions, -1, -1, domains, domain_scores);
  sort(domain_scores.begin(), domain_scores.end());
}


static void
assign_p_values(const vector<double> &random_scores,
		const vector<double> &observed_scores,
		vector<double> &p_values) {
  const double n_randoms =
      random_scores.size() == 0 ? 1 : random_scores.size();
  for (size_t i = 0; i < observed_scores.size(); ++i)
    p_values.push_back((random_scores.end() -
			upper_bound(random_scores.begin(),
				    random_scores.end(),
				    observed_scores[i]))/n_randoms);
}

static void
read_params_file(const bool VERBOSE,
		 const string &params_file,
		 vector<double> &start_trans,
		 vector<vector<double> > &trans,
		 vector<double> &end_trans,
		 double &fdr_cutoff,
		 vector<vector<double> > &alphas,
		 vector<vector<double> > &betas) {
  if (VERBOSE) cerr << "[PARAMETERS:]\n";
  string jnk;
  size_t states, distros;
  std::ifstream in(params_file.c_str());
  in >> jnk >> start_trans[0]
     >> jnk >> start_trans[1]
     >> jnk >> trans[0][0]
     >> jnk >> trans[0][1]
     >> jnk >> trans[1][0]
     >> jnk >> trans[1][1]
     >> jnk >> end_trans[0]
     >> jnk >> end_trans[1]
     >> jnk >> states
     >> jnk >> distros;
  for (size_t i=0; i<states; ++i) {
    alphas.push_back(*(new vector<double>(distros)));
    betas.push_back(*(new vector<double>(distros)));
    for (size_t j=0; j<distros; ++j)
      in >> jnk >> alphas[i][j]
         >> jnk >> betas[i][j];
  }
  in >> jnk >> fdr_cutoff;

  if (VERBOSE) {
    cerr << "S_0\t" << start_trans[0] << endl
	 << "S_1\t" << start_trans[1] << endl
	 << "p_00\t" << trans[0][0] << endl
	 << "p_01\t" << trans[0][1] << endl
	 << "p_10\t" << trans[1][0] << endl
	 << "p_11\t" << trans[1][1] << endl
	 << "0_E\t" << end_trans[0] << endl
	 << "1_E\t" << end_trans[1] << endl
	 << "FDR_CUTOFF\t" << fdr_cutoff << endl;
    for (size_t i=0; i<states; ++i)
      for (size_t j=0; j<distros; ++j)
        cerr << "ALPHA_" << i << j << "\t" << alphas[i][j] << endl
             << "BETA_"  << i << j << "\t" << betas[i][j]  << endl;
  }
}

static void
write_params_file(const string &outfile,
		  const vector<double> &start_trans,
		  const vector<vector<double> > &trans,
		  const vector<double> &end_trans,
		  vector<vector<double> > &alphas,
		  vector<vector<double> > &betas) {

  std::ofstream of;
  if (!outfile.empty()) of.open(outfile.c_str());
  std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());

  out.precision(30);
  out << "S_1\t" << start_trans[0] << endl
      << "S_2\t" << start_trans[1] << endl
      << "1_1\t" << trans[0][0] << endl
      << "1_2\t" << trans[0][1] << endl
      << "2_1\t" << trans[1][0] << endl
      << "2_2\t" << trans[1][1] << endl
      << "1_E\t" << end_trans[0] << endl
      << "2_E\t" << end_trans[1] << endl
      << "STATES\t" << alphas.size() << endl
      << "DISTROS\t" << alphas[0].size() << endl;
  for (size_t i = 0; i < alphas.size(); ++i)
    for (size_t j = 0; j < alphas[0].size(); ++j)
      out << "ALPHA" << i << "_" << j << "\t" << alphas[i][j] << endl
          << "BETA"  << i << "_" << j << "\t" << betas[i][j]  << endl;
}


int
main(int argc, const char **argv) {

  bool VERBOSE = true;
  size_t good_hmr_count = 0;
  try {

    size_t matrix_idx = 0;
    string sampleid = "";

    string outfile;

    int min_cov = 1;
    size_t desert_size = 1000;
    size_t max_iterations = 10;
    double trimrate = 0;
    int merge_regions = 0;
    double lower_bound_rate = -1;
    double upper_bound_rate = -1;
    int min_nr_c = -1;

    // run mode flags
    //bool VERBOSE = true;
    bool PARTIAL_METH = false;
    bool all_regions = false;
    double fdr = 0.01;
    bool viterbi = false;
	bool HUMAN = false;

    // corrections for small values (not parameters):
    double tolerance = 1e-10;
    double min_prob  = 1e-10;

    string params_in_file;
    string params_out_file;

    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]), "Program for identifying "
			   "MRs in methylation data", "<genome_matrix-file/cpg-BED-file>");
    opt_parse.add_opt("idx", 'x', "sample index, if input is genome matrix (specify '1' "
		      "for sample 1 in column 5, ...); do not specify if input is BED", false, matrix_idx);
    opt_parse.add_opt("sampleid", 'y', "sample name, arbitrary, will appear as last column in output",
          false, sampleid);
    opt_parse.add_opt("cov", 'c', "minimum coverage at each cytosine (default: 1)", false, min_cov);
    opt_parse.add_opt("out", 'o', "output file (default: stdout)",
		      false, outfile);
    opt_parse.add_opt("desert", 'd', "maximal distance between covered cytosines in MRs (default: 1000)",
		      false, desert_size);
    opt_parse.add_opt("itr", 'i', "max number of iterations (default: 10)", false, max_iterations);
//    opt_parse.add_opt("partial", '\0', "identify PMRs instead of HMRs",
//		      false, PARTIAL_METH);
	opt_parse.add_opt("human", 'H', "mode to run on human methylation (default: false)", false, HUMAN);
    opt_parse.add_opt("trim", 't', "trim sites with methylation rate <= t (0<=t<=1)"
		      " at start and end of an MR (default: no trimming)", false, trimrate);
    opt_parse.add_opt("merge", 'm', "merge regions not more than m bps apart (default: no merging)",
		      false, merge_regions);
    opt_parse.add_opt("params-in", 'P', "HMM parameters file (no training)",
		      false, params_in_file);
    opt_parse.add_opt("params-out", 'p', "write HMM parameters to this file",
		      false, params_out_file);
    opt_parse.add_opt("fdr", 'f', "FDR cutoff for permutation test (default: 0.01)", false, fdr);
    opt_parse.add_opt("nr-c", 'n', "minimum number of covered cytosines in a MR. Implies skipping "
              " the permutation test", false, min_nr_c);
    opt_parse.add_opt("all", 'a', "report all regions, do not perform the "
		      "permutation test", false, all_regions);
    opt_parse.add_opt("viterbi", 'V', "use Viterbi Decoding instead of "
		      "Posterior Decoding", false, viterbi);
    opt_parse.add_opt("quiet", 'q', "quiet mode; no output on console during runtime", false, VERBOSE);
    opt_parse.add_opt("lb", 'l', "lower bound methylation rate for statistics", false,
		      lower_bound_rate);
    opt_parse.add_opt("ub", 'u', "upper bound methylation rate for statistics", false,
		      upper_bound_rate);

    vector<string> leftover_args;
    opt_parse.parse(argc, argv, leftover_args);
    if (argc == 1 || opt_parse.help_requested()) {
      cerr << opt_parse.help_message() << endl
	   << opt_parse.about_message() << endl;
      return EXIT_SUCCESS;
    }
    if (opt_parse.about_requested()) {
      cerr << opt_parse.about_message() << endl;
      return EXIT_SUCCESS;
    }
    if (opt_parse.option_missing()) {
      cerr << opt_parse.option_missing_message() << endl;
      return EXIT_SUCCESS;
    }
    if (leftover_args.empty()) {
      cerr << opt_parse.help_message() << endl;
      return EXIT_SUCCESS;
    }
    if (trimrate < 0 || trimrate > 1) {
      cerr << "ERROR: Trim rate must be between 0 and 1!\n";
      return EXIT_SUCCESS;
    }
    if (merge_regions >= int(desert_size)) {
      cerr << "ERROR: Distance between regions for merging must be"
	      " smaller than desert size!\n";
      return EXIT_SUCCESS;
    }
    const string cpgs_file = leftover_args.front();

    if (VERBOSE) {
      for (int i=0; i<argc; ++i)
        cerr << argv[i] << " ";
      cerr << endl;
    }
    /****************** END COMMAND LINE OPTIONS *****************/

    // separate the regions by chrom and by desert
    vector<SimpleGenomicRegion> cpgs;
    vector<pair<double, double> > meth;
    vector<char> contexts;

    vector<size_t> reads;
    if (VERBOSE) cerr << "[LOAD CYTOSINES]\n";
    load_cpgs(VERBOSE, PARTIAL_METH, cpgs_file, cpgs, meth, contexts, reads, min_cov, matrix_idx);

    // separate the regions by chrom and by desert, and eliminate
    // those isolated CpGs
    vector<size_t> reset_points;
    separate_regions(VERBOSE, desert_size, cpgs, meth, reads, reset_points);

    // Initialization of HMM:
    vector<double> start_trans(2, 0.5), end_trans(2, 1e-10);
    vector<vector<double> > trans(2, vector<double>(2, 0.25));
    trans[0][0] = trans[1][1] = 0.75;

    double fdr_cutoff = std::numeric_limits<double>::max();
    vector<vector<double> > alphas;
    vector<vector<double> > betas;

    if (!params_in_file.empty()) {
      // READ THE PARAMETERS FILE
      read_params_file(VERBOSE, params_in_file,
		       start_trans, trans, end_trans, fdr_cutoff,
		       alphas, betas);
      max_iterations = 0;
    }
    else {
      const double n_reads =
	accumulate(reads.begin(), reads.end(), 0.0)/reads.size();
      //fg_alpha = 0.33*n_reads;
      //fg_beta = 0.67*n_reads;
      //bg_alpha = 0.67*n_reads;
      //bg_beta = 0.33*n_reads;

      int states = 2;
      for (int s=0; s<states; ++s) {
	alphas.push_back(*(new vector<double>()));
	betas.push_back(*(new vector<double>()));
        for (int c=0; c<3; ++c) {
          alphas[s].push_back(fabs((s==0)-0.8)*n_reads);
          betas[s].push_back(fabs((s==0)-0.2)*n_reads);
        }
      }
    }
//for (size_t i=0; i<alphas.size(); ++i)
//  for (size_t j=0; j<alphas[0].size(); ++j)
//    cout << i << " " << j << ": " << alphas[i][j] << " " << betas[i][j] << endl;

    // Create HMM
    const TwoStateHMMB hmm(min_prob, tolerance, max_iterations, VERBOSE, 0);

    if (max_iterations > 0) {
      if (VERBOSE) cerr << "[TRAINING]" << endl;
      hmm.BaumWelchTraining(meth, contexts, reset_points, start_trans, trans,
			    end_trans, alphas, betas);
    }

    if (!params_out_file.empty()) {
      // WRITE ALL THE HMM PARAMETERS:
      write_params_file(params_out_file, start_trans, trans, end_trans,
			alphas, betas);
    }

    /***********************************
     * STEP 5: DECODE THE DOMAINS
     */
    vector<bool> classes;
    if (VERBOSE) cerr << endl << "[FINAL DECODING]\n";
    if (viterbi == false) {
      vector<double> scores;
      hmm.PosteriorDecoding(meth, contexts, reset_points, start_trans, trans,
			    end_trans, alphas, betas, classes, scores);
    }
    else {
      hmm.ViterbiDecoding(meth, contexts, reset_points, start_trans, trans,
			  end_trans, alphas, betas, classes);
    }

//for (size_t i=0; i<1000; ++i)
//  cout << classes[i];
//cout << endl;

    if (VERBOSE) cerr << "[BUILD AND ASSESS REGIONS]\n";
    vector<GenomicRegion> domains;
    vector<double> domain_scores;
    build_and_score_domains(VERBOSE, HUMAN, cpgs, meth, reset_points, classes,
			    trimrate, merge_regions, lower_bound_rate,
			    upper_bound_rate, domains, domain_scores);

    vector<double> p_values;
    if (all_regions == false && min_nr_c == -1) {
//      vector<double> domain_scores;
//      get_domain_scores(classes, meth, reset_points, domain_scores);

      vector<double> random_scores;
      shuffle_cpgs(hmm, HUMAN, meth, contexts, reset_points, start_trans, trans, end_trans,
  		   alphas, betas, random_scores, cpgs, trimrate, merge_regions, viterbi);

      assign_p_values(random_scores, domain_scores, p_values);

      if (fdr_cutoff == numeric_limits<double>::max())
        fdr_cutoff = get_fdr_cutoff(p_values, fdr);

      if (!params_out_file.empty()) {
        std::ofstream out(params_out_file.c_str(), std::ios::app);
        out << "FDR_CUTOFF\t"
  	  << std::setprecision(30) << fdr_cutoff << endl;
        out.close();
      }
    }

    //vector<GenomicRegion> domains;
    //build_domains(VERBOSE, cpgs, reset_points, classes, domains);

    std::ofstream of;
    if (!outfile.empty()) of.open(outfile.c_str());
    std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());

    //size_t good_hmr_count = 0;
    for (size_t i = 0; i < domains.size(); ++i) {
      if (all_regions == true || (min_nr_c > -1 && domains[i].get_score() >= min_nr_c) || (min_nr_c == -1 && p_values[i] < fdr_cutoff)) {
        if (sampleid == "")
      	  domains[i].set_name("MR" + smithlab::toa(good_hmr_count++));
        else {
          domains[i].set_name(sampleid);
          good_hmr_count++;
        }
      	out << domains[i] << '\n';
      }
    }
  }
  catch (SMITHLABException &e) {
    cerr << "ERROR:\t" << e.what() << endl;
    return EXIT_FAILURE;
  }
  catch (std::bad_alloc &ba) {
    cerr << "ERROR: could not allocate memory" << endl;
    return EXIT_FAILURE;
  }
  if (VERBOSE) cerr << "FOUND " << good_hmr_count << " REGIONS.\n[DONE]\n";
  return EXIT_SUCCESS;
}
