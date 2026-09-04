#include <Rcpp.h>
#include <string>
#include <vector>
#include <map>
#include "bgzf.h"
#include <fstream>
#include <random>

#include <iostream>
#include <cmath>

using namespace Rcpp;

int BgzfGetLine(BGZF* fp, std::string& line);

//' Build Index File for a BGZF Reference Panel
//'
//' Reads a BGZF-compressed reference panel genotype file and writes an index
//' recording the byte offset (\code{fpos}) and line length for each SNP.
//' The resulting index file is required by most other functions in this package.
//'
//' @param reference_data_file Character. Path to the BGZF-compressed reference
//'   genotype file (e.g. \code{"/33KG/33kg_geno.gz"}).
//' @param output_file Character. Path for the output index file.
//'
//' @return Invisibly returns \code{NULL}. Writes the index to \code{output_file}.
//'   Each line contains: \code{fpos line_length}.
//'
//' @examples
//' \dontrun{
//' indexer(
//'   reference_data_file = "/33KG/33kg_geno.gz",
//'   output_file = "33kg_index.txt"
//' )
//' }
//' @export
// [[Rcpp::export]]
void indexer(std::string reference_data_file,
             std::string output_file){
  
  BGZF* fp = bgzf_open(reference_data_file.c_str(), "r");
  std::ofstream outfile;
  outfile.open(output_file.c_str());
  
  if(!fp){
    
    Rcpp::stop("ERROR: can't open reference data file '"+reference_data_file+"'");
  }
  int last_char;
  std::string line;
  long long int fpos;
  
  while(true){
    fpos = bgzf_tell(fp);
    last_char=BgzfGetLine(fp,line);
    if(last_char == -1)
      break;
    outfile<<fpos<<" "<<line.length()<<std::endl;
  }
  
  outfile.close();
  bgzf_close(fp);
}

//' Compute Overall Reference Allele Frequency (AF1) per SNP
//'
//' Reads the BGZF-compressed reference genotype file and computes the
//' alternate allele frequency (AF1) across all populations combined.
//' AF1 is rounded up to 5 decimal places.
//'
//' @param reference_data_file Character. Path to the BGZF-compressed reference
//'   genotype file.
//' @param num_pops Integer. Total number of populations in the reference panel
//'   (e.g. \code{29} for 33KG).
//' @param output_file Character. Path for the output file. Each line contains
//'   the AF1 value for the corresponding SNP.
//'
//' @return Invisibly returns \code{NULL}. Writes one AF1 value per line to
//'   \code{output_file}.
//'
//' @examples
//' \dontrun{
//' cal_af1ref(
//'   reference_data_file = "/33KG/33kg_geno.gz",
//'   num_pops = 29,
//'   output_file = "33kg_af1ref.txt"
//' )
//' }
//' @export
// [[Rcpp::export]]
void cal_af1ref(std::string reference_data_file,
                 int num_pops,
                 std::string output_file){
  BGZF* fp = bgzf_open(reference_data_file.c_str(), "r");
  std::ofstream outfile;
  outfile.open(output_file.c_str());
  if(!fp){
    
    Rcpp::stop("ERROR: can't open reference data file '"+reference_data_file+"'");
  }
  int last_char;
  std::string line;
  double af1ref;
  int num_subj;
  //int cc=0;
  while(true){
    num_subj=0;
    last_char = BgzfGetLine(fp, line);
    if(last_char == -1)
      break;
    std::istringstream buffer(line);
    double allele_counter=0;
    for(int k=0; k<num_pops; k++){
      std::string geno_str;
      buffer >> geno_str;
      num_subj += geno_str.length();
      for(int i=0; i<geno_str.length(); i++){
        allele_counter += (double)(geno_str[i]-'0');
      }
    }
    af1ref = allele_counter/(2*num_subj);
    af1ref = std::ceil(af1ref*100000.0)/100000.0;  //round up to 5 decimal places
    outfile<<std::setprecision(5)<<std::fixed<<af1ref<<std::endl;
    //cc++;
    //if(cc>100)
    //  break;
  }
  outfile.close();
  bgzf_close(fp);
}

//' Extract Genotype Data for a Single Chromosome
//'
//' Extracts all SNP genotype records for a specified chromosome from the
//' BGZF-compressed reference panel, using the index file to seek efficiently.
//'
//' @param chr_num Integer. Chromosome number (1--22).
//' @param num_pops Integer. Total number of populations in the reference panel.
//' @param index_data_file Character. Path to the reference panel index file
//'   (output of \code{\link{indexer}}).
//' @param reference_data_file Character. Path to the BGZF-compressed reference
//'   genotype file.
//' @param ref_out_file Character. Output file path. Each line is a raw genotype
//'   record for one SNP.
//'
//' @return Invisibly returns \code{NULL}. Writes genotype records to
//'   \code{ref_out_file}.
//'
//' @examples
//' \dontrun{
//' extract_chr_data(
//'   chr_num = 22,
//'   num_pops = 29,
//'   index_data_file = "/33KG/33kg_index.gz",
//'   reference_data_file = "/33KG/33kg_geno.gz",
//'   ref_out_file = "33kg_chr22_geno.txt"
//' )
//' }
//' @export
// [[Rcpp::export]]
void extract_chr_data(int chr_num,
                      int num_pops,
                      std::string index_data_file,
                      std::string reference_data_file,
                      std::string ref_out_file){
  
  BGZF* fpi = bgzf_open(index_data_file.c_str(), "r");
  if(!fpi){
    
    Rcpp::stop("ERROR: can't open index data file '"+index_data_file+"'");
  }
  BGZF* fpd = bgzf_open(reference_data_file.c_str(), "r");
  if(!fpd){
    
    Rcpp::stop("ERROR: can't open reference data file '"+reference_data_file+"'");
  }
  
  std::ofstream data_out;
  data_out.open(ref_out_file.c_str());
  
  int last_char;
  std::string index_line, data_line;
  std::string rsid, a1, a2;
  int chr;
  double af1ref;
  long long int bp, fpos;
  
  while(true){
    last_char = BgzfGetLine(fpi, index_line);
    if(last_char == -1) //EOF
      break;
    
    std::istringstream buffer(index_line);
    buffer >> rsid >> chr >> bp >> a1 >> a2 >> af1ref >> fpos;
    
    if(chr==chr_num){
      bgzf_seek(fpd, fpos, SEEK_SET);
      BgzfGetLine(fpd, data_line);
      //write data info
      data_out<<data_line<<std::endl;
    }
  }
  data_out.close();
  bgzf_close(fpi);
  bgzf_close(fpd);
}

//' Extract Genotype Data for Specific Populations on One Chromosome
//'
//' Extracts SNP genotype data for a user-specified subset of populations
//' on a single chromosome. Population names are matched case-insensitively
//' against the reference population description file.
//'
//' @param chr_num Integer. Chromosome number (1--22).
//' @param pop_vec Character vector. Population abbreviations to include
//'   (e.g. \code{c("CEU", "YRI")}). Case-insensitive.
//' @param index_data_file Character. Path to the reference panel index file.
//' @param reference_data_file Character. Path to the BGZF-compressed reference
//'   genotype file.
//' @param reference_pop_desc_file Character. Path to the population description
//'   file (columns: \code{pop_abb num_subj sup_pop_abb}, with a header row).
//' @param ref_out_file Character. Output file path. Each line contains genotype
//'   strings and AF1 values for the selected populations only.
//'
//' @return Invisibly returns \code{NULL}. Writes filtered genotype records to
//'   \code{ref_out_file}.
//'
//' @examples
//' \dontrun{
//' extract_chr_pop_data(
//'   chr_num = 22,
//'   pop_vec = c("CEU", "YRI"),
//'   index_data_file = "/33KG/33kg_index.gz",
//'   reference_data_file = "/33KG/33kg_geno.gz",
//'   reference_pop_desc_file = "/33KG/33kg_pop_desc.txt",
//'   ref_out_file = "33kg_chr22_CEU_YRI_geno.txt"
//' )
//' }
//' @export
// [[Rcpp::export]]
void extract_chr_pop_data(int chr_num,
                      std::vector<std::string> pop_vec,
                      std::string index_data_file,
                      std::string reference_data_file,
                      std::string reference_pop_desc_file,
                      std::string ref_out_file){
  
  // Read pop_vec and convert pops uppercase
  std::vector<std::string> pop_vec_input;
  for(int i=0; i<pop_vec.size(); i++){
    std::string pop = pop_vec[i];
    std::transform(pop.begin(), pop.end(), pop.begin(), ::toupper); //make capital
    pop_vec_input.push_back(pop);
  }
  
  // Read reference_pop_desc_file 
  std::string ref_desc_file = reference_pop_desc_file;
  std::ifstream in_ref_desc(ref_desc_file.c_str());
  
  if(!in_ref_desc){
    Rcpp::Rcout<<std::endl;
    Rcpp::stop("ERROR: can't open reference population description file '"+ref_desc_file+"'");
  }
  std::string line;
  std::string pop_abb, sup_pop_abb;
  int pop_num_subj;
  std::vector<std::string> ref_pop_vec;
  std::vector<int> ref_pop_size_vec;
  std::vector<std::string> ref_sup_pop_vec;
  
  std::getline(in_ref_desc, line); //read header of input file.  
  while(std::getline(in_ref_desc, line)){
    std::istringstream buffer(line);
    buffer >> pop_abb >> pop_num_subj >> sup_pop_abb;
    ref_pop_vec.push_back(pop_abb);
    ref_pop_size_vec.push_back(pop_num_subj);
    ref_sup_pop_vec.push_back(sup_pop_abb);
  }//while
  int num_pops;
  num_pops=ref_pop_vec.size();
  in_ref_desc.close();
  
  // init pop_flag vector
  std::vector<int> pop_flag_vec;
  for(int i=0; i<num_pops; i++){
    std::string pop = ref_pop_vec[i];
    if(std::find(pop_vec_input.begin(), pop_vec_input.end(), pop)!=pop_vec_input.end()){ //if pop is found in pop_vec_input
      pop_flag_vec.push_back(1);
    } else {
      pop_flag_vec.push_back(0);
    }
  }
  
  BGZF* fpi = bgzf_open(index_data_file.c_str(), "r");
  if(!fpi){
    
    Rcpp::stop("ERROR: can't open index data file '"+index_data_file+"'");
  }
  BGZF* fpd = bgzf_open(reference_data_file.c_str(), "r");
  if(!fpd){
    
    Rcpp::stop("ERROR: can't open reference data file '"+reference_data_file+"'");
  }
    
  std::ofstream data_out;
  data_out.open(ref_out_file.c_str());
  
  int last_char;
  std::string index_line, data_line;
  std::string rsid, a1, a2;
  int chr;
  double af1ref;
  long long int bp, fpos;
  
  while(true){
    last_char = BgzfGetLine(fpi, index_line);
    if(last_char == -1) //EOF
      break;
    
    std::istringstream buffer(index_line);
    buffer >> rsid >> chr >> bp >> a1 >> a2 >> af1ref >> fpos;
    
    if(chr==chr_num){
      bgzf_seek(fpd, fpos, SEEK_SET);
      last_char = BgzfGetLine(fpd, data_line);      
      if(last_char == -1) //EOF
        break;
      
      std::istringstream data_buffer(data_line);
      for(int k=0; k<num_pops; k++){
        std::string geno_str;
        data_buffer >> geno_str;
        if(pop_flag_vec[k])
          data_out<<geno_str<<" "; // write genotype string
      }
      for(int k=0; k<num_pops; k++){
        std::string af1_pop;
        data_buffer >> af1_pop;
        if(pop_flag_vec[k])
          data_out<<af1_pop<<" ";  // write af1 of each pop
      }
      //write data info
      data_out<<std::endl;
    }
  }
  data_out.close();
  bgzf_close(fpi);
  bgzf_close(fpd);
  
}

//' Extract Population-Level AF1 for All SNPs on One Chromosome
//'
//' For each SNP on the specified chromosome, reads the per-population allele
//' frequency (AF1) fields from the reference genotype file and writes them
//' to the output — one line per SNP, space-separated AF1 values in the order
//' populations appear in the population description file.
//'
//' @param chr_num Integer. Chromosome number (1--22).
//' @param index_data_file Character. Path to the reference panel index file.
//' @param reference_data_file Character. Path to the BGZF-compressed reference
//'   genotype file.
//' @param reference_pop_desc_file Character. Path to the population description
//'   file (columns: \code{pop_abb num_subj sup_pop_abb}, with a header row).
//' @param ref_out_file Character. Output file path. Each line contains
//'   space-separated AF1 values for all populations for one SNP.
//'
//' @return Invisibly returns \code{NULL}. Writes AF1 matrix rows to
//'   \code{ref_out_file}.
//'
//' @examples
//' \dontrun{
//' extract_all_af1(
//'   chr_num = 22,
//'   index_data_file = "/33KG/33kg_index.gz",
//'   reference_data_file = "/33KG/33kg_geno.gz",
//'   reference_pop_desc_file = "/33KG/33kg_pop_desc.txt",
//'   ref_out_file = "33kg_chr22_af1.txt"
//' )
//' }
//' @export
// [[Rcpp::export]]
void extract_all_af1(int chr_num,
                     std::string index_data_file,
                     std::string reference_data_file,
                     std::string reference_pop_desc_file,
                     std::string ref_out_file){
  
  // Read reference_pop_desc_file 
  std::string ref_desc_file = reference_pop_desc_file;
  std::ifstream in_ref_desc(ref_desc_file.c_str());
  
  if(!in_ref_desc){
    Rcpp::Rcout<<std::endl;
    Rcpp::stop("ERROR: can't open reference population description file '"+ref_desc_file+"'");
  }
  std::string line;
  std::string pop_abb, sup_pop_abb;
  int pop_num_subj;
  std::vector<std::string> ref_pop_vec;
  std::vector<int> ref_pop_size_vec;
  std::vector<std::string> ref_sup_pop_vec;
  
  std::getline(in_ref_desc, line); //read header of input file.  
  while(std::getline(in_ref_desc, line)){
    std::istringstream buffer(line);
    buffer >> pop_abb >> pop_num_subj >> sup_pop_abb;
    ref_pop_vec.push_back(pop_abb);
    ref_pop_size_vec.push_back(pop_num_subj);
    ref_sup_pop_vec.push_back(sup_pop_abb);
  }//while
  int num_pops;
  num_pops=ref_pop_vec.size();
  in_ref_desc.close();
  
  BGZF* fpi = bgzf_open(index_data_file.c_str(), "r");
  if(!fpi){
    
    Rcpp::stop("ERROR: can't open index data file '"+index_data_file+"'");
  }
  BGZF* fpd = bgzf_open(reference_data_file.c_str(), "r");
  if(!fpd){
    
    Rcpp::stop("ERROR: can't open reference data file '"+reference_data_file+"'");
  }
  
  std::ofstream data_out;
  data_out.open(ref_out_file.c_str());
  
  int last_char;
  std::string index_line, data_line;
  std::string rsid, a1, a2;
  int chr;
  double af1ref;
  long long int bp, fpos;
  
  while(true){
    last_char = BgzfGetLine(fpi, index_line);
    if(last_char == -1) //EOF
      break;
    
    std::istringstream buffer(index_line);
    buffer >> rsid >> chr >> bp >> a1 >> a2 >> af1ref >> fpos;
    
    if(chr==chr_num){
      bgzf_seek(fpd, fpos, SEEK_SET);
      last_char = BgzfGetLine(fpd, data_line);      
      if(last_char == -1) //EOF
        break;
      
      std::istringstream data_buffer(data_line);
      for(int k=0; k<num_pops; k++){
        std::string geno_str;
        data_buffer >> geno_str; //skip genotypes
      }
      for(int k=0; k<num_pops; k++){
        std::string af1_pop;
        data_buffer >> af1_pop;
        data_out<<af1_pop<<" ";  // write af1 of each pop
      }
      //write data info
      data_out<<std::endl;
    }
  }
  data_out.close();
  bgzf_close(fpi);
  bgzf_close(fpd);
}



// Internal helper: chr_num == -1 means process all chromosomes
static void simulate_af1_z_impl(
    int chr_num,
    std::vector<std::string> pop_vec,
    std::vector<int> num_sim_vec,
    std::string index_data_file,
    std::string reference_data_file,
    std::string reference_pop_desc_file,
    std::string ref_out_file)
{
  // Normalise population names to uppercase
  std::vector<std::string> pop_vec_input;
  for(int i = 0; i < (int)pop_vec.size(); i++){
    std::string pop = pop_vec[i];
    std::transform(pop.begin(), pop.end(), pop.begin(), ::toupper);
    pop_vec_input.push_back(pop);
  }

  // Parse population description file
  std::ifstream in_ref_desc(reference_pop_desc_file.c_str());
  if(!in_ref_desc)
    Rcpp::stop("ERROR: can't open reference population description file '" +
               reference_pop_desc_file + "'");

  std::string line, pop_abb, sup_pop_abb;
  int pop_num_subj;
  std::vector<std::string> ref_pop_vec;
  std::vector<int>         ref_pop_size_vec;

  std::getline(in_ref_desc, line); // skip header
  while(std::getline(in_ref_desc, line)){
    std::istringstream buf(line);
    buf >> pop_abb >> pop_num_subj >> sup_pop_abb;
    ref_pop_vec.push_back(pop_abb);
    ref_pop_size_vec.push_back(pop_num_subj);
  }
  int num_pops = ref_pop_vec.size();
  in_ref_desc.close();

  // Validate inputs
  if(pop_vec_input.size() != num_sim_vec.size())
    Rcpp::stop("ERROR: pop_vec and num_sim_vec must be the same length.");
  if(pop_vec_input.empty())
    Rcpp::stop("ERROR: pop_vec must not be empty.");

  // Map each requested population onto its position in the reference panel.
  // num_sim_by_ref is indexed by reference-panel population order, so that
  // pop_vec may be given in any order without permuting the sample sizes.
  std::vector<int> pop_flag_vec(num_pops, 0);
  std::vector<int> num_sim_by_ref(num_pops, 0);
  for(int i = 0; i < (int)pop_vec_input.size(); i++){
    std::vector<std::string>::iterator it =
      std::find(ref_pop_vec.begin(), ref_pop_vec.end(), pop_vec_input[i]);
    if(it == ref_pop_vec.end())
      Rcpp::stop("ERROR: population '" + pop_vec_input[i] +
                 "' not found in reference population description file.");
    int ref_idx = (int)(it - ref_pop_vec.begin());
    if(pop_flag_vec[ref_idx])
      Rcpp::stop("ERROR: population '" + pop_vec_input[i] +
                 "' is listed more than once in pop_vec.");
    if(num_sim_vec[i] <= 0)
      Rcpp::stop("ERROR: num_sim_vec must contain positive values (population '" +
                 pop_vec_input[i] + "').");
    pop_flag_vec[ref_idx]   = 1;
    num_sim_by_ref[ref_idx] = num_sim_vec[i];
  }

  int total_num_subj = 0;
  for(int i = 0; i < num_pops; i++) total_num_subj += num_sim_by_ref[i];
  if(total_num_subj < 3)
    Rcpp::stop("ERROR: the total number of simulated subjects must be at "
               "least 3 to compute an association Z-score.");

  // Open BGZF files
  BGZF* fpi = bgzf_open(index_data_file.c_str(), "r");
  if(!fpi)
    Rcpp::stop("ERROR: can't open index data file '" + index_data_file + "'");
  BGZF* fpd = bgzf_open(reference_data_file.c_str(), "r");
  if(!fpd)
    Rcpp::stop("ERROR: can't open reference data file '" + reference_data_file + "'");

  std::ofstream data_out(ref_out_file.c_str());
  data_out << "rsid chr bp a1 a2 sim_af1 sim_z" << std::endl;

  // Draw bootstrap sample indices and null phenotype once, reused across SNPs
  std::random_device rd;
  std::mt19937 gen(rd());
  std::normal_distribution<> norm_dist(0, 1);

  std::vector<double> response(total_num_subj);
  for(int i = 0; i < total_num_subj; i++) response[i] = norm_dist(gen);

  std::vector<int> geno_index_vec;
  for(int k = 0; k < num_pops; k++){
    if(pop_flag_vec[k]){
      std::uniform_int_distribution<> dis(0, ref_pop_size_vec[k] - 1);
      for(int j = 0; j < num_sim_by_ref[k]; j++)
        geno_index_vec.push_back(dis(gen));
    }
  }

  // Main loop over index file
  int last_char;
  std::string index_line, data_line, rsid, a1, a2;
  int chr;
  double af1ref;
  long long int bp, fpos;

  while(true){
    last_char = BgzfGetLine(fpi, index_line);
    if(last_char == -1) break;

    std::istringstream ibuf(index_line);
    ibuf >> rsid >> chr >> bp >> a1 >> a2 >> af1ref >> fpos;

    // chr_num == -1 means all chromosomes; otherwise filter
    if(chr_num != -1 && chr != chr_num) continue;

    bgzf_seek(fpd, fpos, SEEK_SET);
    last_char = BgzfGetLine(fpd, data_line);
    if(last_char == -1) break;

    // Extract bootstrap genotypes for selected populations
    std::istringstream dbuf(data_line);
    std::vector<double> geno_vec;
    double allele_counter = 0.0;
    int subj_counter = 0;
    for(int k = 0; k < num_pops; k++){
      std::string geno_str;
      dbuf >> geno_str;
      if(pop_flag_vec[k]){
        for(int j = 0; j < num_sim_by_ref[k]; j++){
          double geno = (double)(geno_str[geno_index_vec[j + subj_counter]] - '0');
          allele_counter += geno;
          geno_vec.push_back(geno);
        }
        subj_counter += num_sim_by_ref[k];
      }
    }

    // Simulated AF1
    double sim_af1 = std::ceil(allele_counter / (2.0 * total_num_subj) * 100000.0) / 100000.0;

    // OLS association Z-score under null
    double x_mean = std::accumulate(geno_vec.begin(), geno_vec.end(), 0.0) / geno_vec.size();
    double y_mean = std::accumulate(response.begin(), response.end(), 0.0) / response.size();
    double sum_xy = std::inner_product(geno_vec.begin(), geno_vec.end(), response.begin(), 0.0);
    double sum_xx = std::inner_product(geno_vec.begin(), geno_vec.end(), geno_vec.begin(), 0.0);
    double sxy = sum_xy - (int)geno_vec.size() * x_mean * y_mean;
    double sxx = sum_xx - (int)geno_vec.size() * std::pow(x_mean, 2);
    double beta1 = sxy / sxx;
    double beta0 = y_mean - beta1 * x_mean;

    double sse = 0.0;
    for(int i = 0; i < (int)geno_vec.size(); i++){
      double r = response[i] - (beta0 + beta1 * geno_vec[i]);
      sse += r * r;
    }
    double std_err = std::sqrt((sse / (geno_vec.size() - 2)) / sxx);
    double sim_z = std::ceil((beta1 / std_err) * 100000.0) / 100000.0;

    data_out << rsid << " " << chr << " " << bp << " " << a1 << " " << a2 << " "
             << std::setprecision(5) << std::fixed << sim_af1 << " "
             << std::setprecision(5) << std::fixed << sim_z << std::endl;
  }

  data_out.close();
  bgzf_close(fpi);
  bgzf_close(fpd);
}


//' Simulate AF1 and Association Z-scores Across All Chromosomes
//'
//' For each SNP across all chromosomes, draws a bootstrap sample of subjects
//' from specified populations (sampling with replacement), simulates a null
//' phenotype (standard normal), and computes the OLS association Z-score.
//' Use this to generate genome-wide null Z-score distributions.
//'
//' @param pop_vec Character vector. Population abbreviations to sample from
//'   (e.g. `c("CEU", "YRI")`). Case-insensitive. Each population may appear
//'   at most once.
//' @param num_sim_vec Integer vector. Number of subjects to sample from each
//'   population, where `num_sim_vec[i]` is the sample size for `pop_vec[i]`.
//'   Must be the same length as `pop_vec` and contain positive values. The
//'   pairing is positional, so `pop_vec` may be given in any order; it need
//'   not follow the order of the population description file.
//' @param index_data_file Character. Path to the reference panel index file.
//' @param reference_data_file Character. Path to the BGZF-compressed reference
//'   genotype file.
//' @param reference_pop_desc_file Character. Path to the population description
//'   file.
//' @param ref_out_file Character. Output file path. Space-separated columns:
//'   `rsid chr bp a1 a2 sim_af1 sim_z`.
//'
//' @return Invisibly returns `NULL`. Writes simulation results to
//'   `ref_out_file`.
//'
//' @seealso [simulate_af1_z()] for single-chromosome simulation.
//'
//' @examples
//' \dontrun{
//' simulate_af1_z_allchr(
//'   pop_vec = c("CEU", "YRI"),
//'   num_sim_vec = c(100, 100),
//'   index_data_file = "/33KG/33kg_index.gz",
//'   reference_data_file = "/33KG/33kg_geno.gz",
//'   reference_pop_desc_file = "/33KG/33kg_pop_desc.txt",
//'   ref_out_file = "sim_allchr_results.txt"
//' )
//' }
//' @export
// [[Rcpp::export]]
void simulate_af1_z_allchr(std::vector<std::string> pop_vec,
                            std::vector<int> num_sim_vec,
                            std::string index_data_file,
                            std::string reference_data_file,
                            std::string reference_pop_desc_file,
                            std::string ref_out_file){
  simulate_af1_z_impl(-1, pop_vec, num_sim_vec, index_data_file,
                      reference_data_file, reference_pop_desc_file, ref_out_file);
}


//' Simulate AF1 and Association Z-scores for One Chromosome
//'
//' Same simulation procedure as [simulate_af1_z_allchr()] but restricted to a
//' single chromosome. Bootstrap subject indices are drawn once and reused
//' across all SNPs on that chromosome (consistent sampling).
//'
//' @param chr_num Integer. Chromosome number (1--22).
//' @param pop_vec Character vector. Population abbreviations to sample from.
//'   Case-insensitive. Each population may appear at most once.
//' @param num_sim_vec Integer vector. Number of subjects to sample per
//'   population, where `num_sim_vec[i]` is the sample size for `pop_vec[i]`.
//'   Must be the same length as `pop_vec` and contain positive values. The
//'   pairing is positional, so `pop_vec` may be given in any order; it need
//'   not follow the order of the population description file.
//' @param index_data_file Character. Path to the reference panel index file.
//' @param reference_data_file Character. Path to the BGZF-compressed reference
//'   genotype file.
//' @param reference_pop_desc_file Character. Path to the population description
//'   file.
//' @param ref_out_file Character. Output file path. Space-separated columns:
//'   `rsid chr bp a1 a2 sim_af1 sim_z`.
//'
//' @return Invisibly returns `NULL`. Writes simulation results to
//'   `ref_out_file`.
//'
//' @seealso [simulate_af1_z_allchr()] for genome-wide simulation.
//'
//' @examples
//' \dontrun{
//' simulate_af1_z(
//'   chr_num = 22,
//'   pop_vec = c("CEU", "YRI"),
//'   num_sim_vec = c(100, 100),
//'   index_data_file = "/33KG/33kg_index.gz",
//'   reference_data_file = "/33KG/33kg_geno.gz",
//'   reference_pop_desc_file = "/33KG/33kg_pop_desc.txt",
//'   ref_out_file = "sim_chr22_results.txt"
//' )
//' }
//' @export
// [[Rcpp::export]]
void simulate_af1_z(int chr_num,
                    std::vector<std::string> pop_vec,
                    std::vector<int> num_sim_vec,
                    std::string index_data_file,
                    std::string reference_data_file,
                    std::string reference_pop_desc_file,
                    std::string ref_out_file){
  simulate_af1_z_impl(chr_num, pop_vec, num_sim_vec, index_data_file,
                      reference_data_file, reference_pop_desc_file, ref_out_file);
}



//' Extract Genotype Data for a Genomic Region
//'
//' Extracts all SNP genotype records within a specified base-pair range
//' on a chromosome. Uses the index file for efficient random access into
//' the BGZF-compressed reference panel.
//'
//' @param chr_num Integer. Chromosome number (1--22).
//' @param start_bp Integer. Start base pair position (inclusive).
//' @param end_bp Integer. End base pair position (inclusive).
//' @param num_pops Integer. Total number of populations in the reference panel.
//' @param index_data_file Character. Path to the reference panel index file.
//' @param reference_data_file Character. Path to the BGZF-compressed reference
//'   genotype file.
//' @param ref_out_file Character. Output file path. Each line is a raw genotype
//'   record for one SNP in the region.
//'
//' @return Invisibly returns \code{NULL}. Writes genotype records to
//'   \code{ref_out_file}.
//'
//' @examples
//' \dontrun{
//' extract_reg_data(
//'   chr_num = 14,
//'   start_bp = 104000000,
//'   end_bp = 104200000,
//'   num_pops = 29,
//'   index_data_file = "/33KG/33kg_index.gz",
//'   reference_data_file = "/33KG/33kg_geno.gz",
//'   ref_out_file = "33kg_chr14_region_geno.txt"
//' )
//' }
//' @export
// [[Rcpp::export]]
void extract_reg_data(int chr_num,
                      int start_bp,
                      int end_bp,
                      int num_pops,
                      std::string index_data_file,
                      std::string reference_data_file,
                      std::string ref_out_file){
  
  BGZF* fpi = bgzf_open(index_data_file.c_str(), "r");
  if(!fpi){
    
    Rcpp::stop("ERROR: can't open index data file '"+index_data_file+"'");
  }
  BGZF* fpd = bgzf_open(reference_data_file.c_str(), "r");
  if(!fpd){
    
    Rcpp::stop("ERROR: can't open reference data file '"+reference_data_file+"'");
  }
  
  std::ofstream data_out;
  data_out.open(ref_out_file.c_str());
  
  int last_char;
  std::string index_line, data_line;
  std::string rsid, a1, a2;
  int chr;
  double af1ref;
  long long int bp, fpos;
  
  while(true){
    last_char = BgzfGetLine(fpi, index_line);
    if(last_char == -1) //EOF
      break;
    
    std::istringstream buffer(index_line);
    buffer >> rsid >> chr >> bp >> a1 >> a2 >> af1ref >> fpos;
    
    if(chr==chr_num && (bp >= start_bp && bp <= end_bp)){
      bgzf_seek(fpd, fpos, SEEK_SET);
      BgzfGetLine(fpd, data_line);
      //write data info
      data_out<<data_line<<std::endl;
    }
  }
  data_out.close();
  bgzf_close(fpi);
  bgzf_close(fpd);
}



//' Test Whether a BGZF File Can Be Opened
//'
//' Opens the specified BGZF-compressed file, reads the first line, prints it
//' to the console, and closes the file. Use this to verify that a reference
//' panel file is readable before running longer extraction jobs.
//'
//' @param gz_file Character. Path to the BGZF-compressed file to test.
//'
//' @return Invisibly returns \code{NULL}. Prints the first line of the file
//'   to the console.
//'
//' @examples
//' \dontrun{
//' test_gz_file("/33KG/33kg_geno.gz")
//' }
//' @export
// [[Rcpp::export]]
void test_gz_file(std::string gz_file){
  BGZF* fp = bgzf_open(gz_file.c_str(), "r");
  if(!fp){
    
    Rcpp::stop("ERROR: can't open index data file '"+gz_file+"'");
  }
  std::string line;
  BgzfGetLine(fp, line);
  Rcpp::Rcout << line << std::endl;
  bgzf_close(fp);
}


//' @name get_geno_info
//' @keywords internal
// [[Rcpp::export]]
std::string get_geno_info(int64_t fpos,
                   std::string reference_data_file){

  BGZF* fp = bgzf_open(reference_data_file.c_str(), "r");
  if(!fp){
    
    Rcpp::stop("ERROR: can't open reference data file '"+reference_data_file+"'");
  }
  Rcpp::Rcout << "fpos: " << fpos << std::endl;
  std::string line;
  bgzf_seek(fp, fpos, SEEK_SET);
  BgzfGetLine(fp, line);
  
  bgzf_close(fp);
  return line;
}






//' @name largeval
//' @keywords internal
// [[Rcpp::export]]
Rcpp::NumericVector largeval(int64_t val) {
  //int64_t val = 9223372036854775807LL - 1;
  Rcpp::Rcout << "C++ value: " << val << "\n";
  Rcpp::NumericVector dbl(1);
  std::memcpy(&(dbl[0]), &val, sizeof(double));
  dbl.attr("class") = "integer64";
  return dbl;
}



int BgzfGetLine(BGZF* fp, std::string& line){
  line.erase();
  int i=0;
  int c;
  while(true){
    i++;
    c = bgzf_getc(fp);
    if(c == -2){
      Rcpp::stop("ERROR: can't read " + std::to_string(i) + "-th character from BGZF file");
    }
    if(c == -1){ // end of file                                                                                 
      break;
    }
    if(c == 10){ // end of line                                                                                 
      break;
    }
    line += static_cast<char>(c);
  }
  return c;
}

void LoadProgressBar( int percent ){
  std::string bar;
  for(int i = 0; i < 50; i++){
    if( i < (percent/2)){
      bar.replace(i,1,"=");
    }else if( i == (percent/2)){
      bar.replace(i,1,">");
    }else{
      bar.replace(i,1," ");
    }
  }
  Rcpp::Rcout << "\r" << percent << "%  " << "[" << bar << "] " << std::flush;
}


//' @name simulate_zscore
//' @keywords internal
// [[Rcpp::export]]
void simulate_zscore(){
  // Initialize vectors
  std::vector<double> response(100);
  std::vector<int> predictor(100);
  
  // Fill response vector with standard normal values
  std::random_device rd;
  std::mt19937 rng(rd());
  std::normal_distribution<double> dist(0.0, 1.0);
  
  for (int i = 0; i < 100; i++) {
    response[i] = dist(rng);
  }
  
  // Fill predictor vector with (0, 1, 2) values
  for (int i = 0; i < 100; i++) {
    predictor[i] = i % 3;
  }
  
  // Calculate the least square estimators of beta0 and beta1
  double x_mean = accumulate(predictor.begin(), predictor.end(), 0.0) / predictor.size();
  double y_mean = accumulate(response.begin(), response.end(), 0.0) / response.size();
  double sum_xy = inner_product(predictor.begin(), predictor.end(), response.begin(), 0.0);
  double sum_xx = inner_product(predictor.begin(), predictor.end(), predictor.begin(), 0.0);
  double beta1 = (sum_xy - predictor.size() * x_mean * y_mean) / (sum_xx - predictor.size() * pow(x_mean, 2));
  double beta0 = y_mean - beta1 * x_mean;
  
  // Calculate the standard error of beta1
  double sum_residuals_squared = 0;
  for (int i = 0; i < predictor.size(); i++) {
    double residual = response[i] - (beta0 + beta1 * predictor[i]);
    sum_residuals_squared += pow(residual, 2);
  }
  double mse = sum_residuals_squared / (predictor.size() - 2);
  double std_err = sqrt(mse / sum_xx);
  
  // Calculate two-sided z-score for the regression coefficient
  double z_c1 = beta1 / std_err;
  //double p_c1 = 2 * (1 - gsl_cdf_tdist(fabs(t_c1), predictor.size() - 2));
  //double z_c1 = gsl_cdf_gaussian_Pinv(p_c1 / 2, 1);
  
  // Print two-sided z-score for the regression coefficient
  Rcpp::Rcout << "Two-sided z-score for the regression coefficient: " << z_c1 << std::endl;
}
  
