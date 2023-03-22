#include <Rcpp.h>
#include <string>
#include <vector>
#include <map>
#include "bgzf.h"
#include <fstream>
#include <random>

#include <iostream>
#include <cmath>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

using namespace Rcpp;

int BgzfGetLine(BGZF* fp, std::string& line);

// [[Rcpp::export]]
void indexer(std::string reference_data_file,
             std::string output_file){
  
  BGZF* fp = bgzf_open(reference_data_file.c_str(), "r");
  std::ofstream outfile;
  outfile.open(output_file.c_str());
  
  if(!fp){
    std::cout<<std::endl;
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

// [[Rcpp::export]]
void cal_af1ref(std::string reference_data_file,
                 int num_pops,
                 std::string output_file){
  BGZF* fp = bgzf_open(reference_data_file.c_str(), "r");
  std::ofstream outfile;
  outfile.open(output_file.c_str());
  if(!fp){
    std::cout<<std::endl;
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
    //std::cout<<std::setprecision(5)<<std::fixed<<af1ref<<std::endl;
    //cc++;
    //if(cc>100)
    //  break;
  }
  outfile.close();
  bgzf_close(fp);
}

// [[Rcpp::export]]
void extract_chr_data(int chr_num,
                      int num_pops,
                      std::string index_data_file,
                      std::string reference_data_file,
                      std::string ref_out_file){
  
  BGZF* fpi = bgzf_open(index_data_file.c_str(), "r");
  if(!fpi){
    std::cout<<std::endl;
    Rcpp::stop("ERROR: can't open index data file '"+index_data_file+"'");
  }
  BGZF* fpd = bgzf_open(reference_data_file.c_str(), "r");
  if(!fpd){
    std::cout<<std::endl;
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
    //std::cout<<pop<<std::endl;
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
    std::cout<<std::endl;
    Rcpp::stop("ERROR: can't open index data file '"+index_data_file+"'");
  }
  BGZF* fpd = bgzf_open(reference_data_file.c_str(), "r");
  if(!fpd){
    std::cout<<std::endl;
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
    std::cout<<std::endl;
    Rcpp::stop("ERROR: can't open index data file '"+index_data_file+"'");
  }
  BGZF* fpd = bgzf_open(reference_data_file.c_str(), "r");
  if(!fpd){
    std::cout<<std::endl;
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



// [[Rcpp::export]]
void simulate_af1_z(int chr_num,
                  std::vector<std::string> pop_vec,
                  std::vector<int> num_sim_vec,
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
  
  // total number of bootstrap samples
  int total_num_subj=0;
  for(int i=0; i<num_sim_vec.size();i++){
    total_num_subj += num_sim_vec[i];
  }
  
  // open reference panel index data
  BGZF* fpi = bgzf_open(index_data_file.c_str(), "r");
  if(!fpi){
    std::cout<<std::endl;
    Rcpp::stop("ERROR: can't open index data file '"+index_data_file+"'");
  }
  // open reference panel genotype data
  BGZF* fpd = bgzf_open(reference_data_file.c_str(), "r");
  if(!fpd){
    std::cout<<std::endl;
    Rcpp::stop("ERROR: can't open reference data file '"+reference_data_file+"'");
  }
  
  // open output file
  std::ofstream data_out;
  data_out.open(ref_out_file.c_str());
  // write header
  data_out<<"rsid chr bp a1 a2 sim_af1 sim_z"<<std::endl;
  //data_out<<"rsid chr bp a1 a2 sim_af1 beta1 beta0 mse sxy sxx std_err sim_z"<<std::endl;
  
  int last_char;
  std::string index_line, data_line;
  std::string rsid, a1, a2;
  int chr;
  double af1ref;
  long long int bp, fpos;
  
  // Simulate response variable (standard normal variates)
  std::random_device rd;
  std::mt19937 gen(rd());
  std::normal_distribution<> dist(0, 1);
  std::vector<double> response; // Used to store simulated response variable values.
  for (int i = 0; i < total_num_subj; i++) {
    double z = dist(gen);
    response.push_back(z);
  }  
  
  /*
  const gsl_rng_type* T;
  gsl_rng* rng;
  gsl_rng_env_setup();
  T = gsl_rng_default;
  rng = gsl_rng_alloc(T);
  for (int i = 0; i < total_num_subj; i++) {
    double z = gsl_ran_ugaussian(rng);
    response.push_back(z);
  }
  gsl_rng_free(rng);
  */
  
  while(true){
    last_char = BgzfGetLine(fpi, index_line);
    if(last_char == -1) //EOF
      break;
    
    std::istringstream buffer(index_line);
    buffer >> rsid >> chr >> bp >> a1 >> a2 >> af1ref >> fpos;
    
    double allele_counter=0.0;
    if(chr==chr_num){
      
      bgzf_seek(fpd, fpos, SEEK_SET);
      last_char = BgzfGetLine(fpd, data_line);      
      if(last_char == -1) //EOF
        break;
      
      std::istringstream data_buffer(data_line);
      std::vector<double> geno_vec;
      int pop_counter=0;
      for(int k=0; k<num_pops; k++){
        std::string geno_str;
        data_buffer >> geno_str;
        if(pop_flag_vec[k]) {
          for(int j=0; j<num_sim_vec[pop_counter]; j++){
            std::uniform_int_distribution<> dis(0, ref_pop_size_vec[k] - 1);
            int ran_index = dis(gen);
            double geno = (double)(geno_str[ran_index] - '0');
            allele_counter += geno;
            geno_vec.push_back(geno);
          }
          pop_counter++;
        }
      }
      // compute af1 of simulated genotype
      double sim_af1 = allele_counter/(2*total_num_subj);
      sim_af1 = std::ceil(sim_af1*100000.0)/100000.0;  //round up to 5 decimal places
      
      // Compute association Z-score using simulated genotype and phenotype under null
      // Calculate the least square estimators of beta0 and beta1
      double x_mean = accumulate(geno_vec.begin(), geno_vec.end(), 0.0) / geno_vec.size();
      double y_mean = accumulate(response.begin(), response.end(), 0.0) / response.size();
      double sum_xy = inner_product(geno_vec.begin(), geno_vec.end(), response.begin(), 0.0);
      double sum_xx = inner_product(geno_vec.begin(), geno_vec.end(), geno_vec.begin(), 0.0);
      double sxy = sum_xy - geno_vec.size() * x_mean * y_mean;
      double sxx = sum_xx - geno_vec.size() * pow(x_mean, 2);
      double beta1 = sxy/sxx;
      double beta0 = y_mean - beta1 * x_mean;
      
      // Calculate the standard error of beta1
      double sum_residuals_squared = 0;
      for (int i = 0; i < geno_vec.size(); i++) {
        double residual = response[i] - (beta0 + beta1 * geno_vec[i]);
        sum_residuals_squared += pow(residual, 2);
      }
      double mse = sum_residuals_squared / (geno_vec.size() - 2);
      double std_err = sqrt(mse/sxx);
      
      // Calculate two-sided z-score for the regression coefficient
      double sim_z = beta1 / std_err;
      sim_z = std::ceil(sim_z*100000.0)/100000.0;  //round up to 5 decimal places  
      // write results in file
      data_out<<rsid<<" "<<chr<<" "<<bp<<" "<<a1<<" "<<a2<<" ";
      data_out<<std::setprecision(5)<<std::fixed<<sim_af1<<" ";
      //data_out<<std::setprecision(5)<<std::fixed<<beta1<<" ";
      //data_out<<std::setprecision(5)<<std::fixed<<beta0<<" ";
      //data_out<<std::setprecision(5)<<std::fixed<<mse<<" ";
      //data_out<<std::setprecision(5)<<std::fixed<<sxy<<" ";
      //data_out<<std::setprecision(5)<<std::fixed<<sxx<<" ";
      //data_out<<std::setprecision(5)<<std::fixed<<std_err<<" ";
      data_out<<std::setprecision(5)<<std::fixed<<sim_z<<std::endl;
    }
  }
  data_out.close(); //close output filestream
  bgzf_close(fpi);  //close reference index file
  bgzf_close(fpd);  //close reference data file
  
}

// [[Rcpp::export]]
void simulate_af1_z2(int chr_num,
                    std::vector<std::string> pop_vec,
                    std::vector<int> num_sim_vec,
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
  
  // total number of bootstrap samples
  int total_num_subj=0;
  for(int i=0; i<num_sim_vec.size();i++){
    total_num_subj += num_sim_vec[i];
  }
  
  // open reference panel index data
  BGZF* fpi = bgzf_open(index_data_file.c_str(), "r");
  if(!fpi){
    std::cout<<std::endl;
    Rcpp::stop("ERROR: can't open index data file '"+index_data_file+"'");
  }
  // open reference panel genotype data
  BGZF* fpd = bgzf_open(reference_data_file.c_str(), "r");
  if(!fpd){
    std::cout<<std::endl;
    Rcpp::stop("ERROR: can't open reference data file '"+reference_data_file+"'");
  }
  
  // open output file
  std::ofstream data_out;
  data_out.open(ref_out_file.c_str());
  // write header
  data_out<<"rsid chr bp a1 a2 sim_af1 sim_z"<<std::endl;
  //data_out<<"rsid chr bp a1 a2 sim_af1 beta1 beta0 mse sxy sxx std_err sim_z"<<std::endl;
  
  int last_char;
  std::string index_line, data_line;
  std::string rsid, a1, a2;
  int chr;
  double af1ref;
  long long int bp, fpos;
  
  // Simulate response variable (standard normal variates)
  std::random_device rd;
  std::mt19937 gen(rd());
  std::normal_distribution<> dist(0, 1);
  std::vector<double> response; // Used to store simulated response variable values.
  for (int i = 0; i < total_num_subj; i++) {
    double z = dist(gen);
    response.push_back(z);
  }
  
  // To randomly draw with replacement a subject's genotypes from 
  // each combination of ethnicities, randomly select indexes 
  // of subjects in each ethnic group and store the information 
  // in geno_index_vec. The geno_index_vec will be used to simulate
  // genotypes
  std::vector<int> geno_index_vec;
  int pop_counter=0;
  for(int k=0; k<num_pops; k++){
    if(pop_flag_vec[k]) {
      for(int j=0; j<num_sim_vec[pop_counter]; j++){
        std::uniform_int_distribution<> dis(0, ref_pop_size_vec[k] - 1);
        int ran_index = dis(gen);
        geno_index_vec.push_back(ran_index);
      }
      pop_counter++;
    }
  }
  
  while(true){
    last_char = BgzfGetLine(fpi, index_line);
    if(last_char == -1) //EOF
      break;
    
    std::istringstream buffer(index_line);
    buffer >> rsid >> chr >> bp >> a1 >> a2 >> af1ref >> fpos;
    
    double allele_counter=0.0;
    if(chr==chr_num){
      
      bgzf_seek(fpd, fpos, SEEK_SET);
      last_char = BgzfGetLine(fpd, data_line);      
      if(last_char == -1) //EOF
        break;
      
      std::istringstream data_buffer(data_line);
      std::vector<double> geno_vec;
      int pop_counter=0;
      int subj_counter=0;
      for(int k=0; k<num_pops; k++){
        std::string geno_str;
        data_buffer >> geno_str;
        if(pop_flag_vec[k]) {
          for(int j=0; j<num_sim_vec[pop_counter]; j++){
            //std::uniform_int_distribution<> dis(0, ref_pop_size_vec[k] - 1);
            int ran_index = geno_index_vec[j+subj_counter];
            double geno = (double)(geno_str[ran_index] - '0');
            allele_counter += geno;
            geno_vec.push_back(geno);
          }
          subj_counter = subj_counter + num_sim_vec[pop_counter];
          pop_counter++;
        }
      }
      // compute af1 of simulated genotype
      double sim_af1 = allele_counter/(2*total_num_subj);
      sim_af1 = std::ceil(sim_af1*100000.0)/100000.0;  //round up to 5 decimal places
      
      // Compute association Z-score using simulated genotype and phenotype under null
      // Calculate the least square estimators of beta0 and beta1
      double x_mean = accumulate(geno_vec.begin(), geno_vec.end(), 0.0) / geno_vec.size();
      double y_mean = accumulate(response.begin(), response.end(), 0.0) / response.size();
      double sum_xy = inner_product(geno_vec.begin(), geno_vec.end(), response.begin(), 0.0);
      double sum_xx = inner_product(geno_vec.begin(), geno_vec.end(), geno_vec.begin(), 0.0);
      double sxy = sum_xy - geno_vec.size() * x_mean * y_mean;
      double sxx = sum_xx - geno_vec.size() * pow(x_mean, 2);
      double beta1 = sxy/sxx;
      double beta0 = y_mean - beta1 * x_mean;
      
      // Calculate the standard error of beta1
      double sum_residuals_squared = 0;
      for (int i = 0; i < geno_vec.size(); i++) {
        double residual = response[i] - (beta0 + beta1 * geno_vec[i]);
        sum_residuals_squared += pow(residual, 2);
      }
      double mse = sum_residuals_squared / (geno_vec.size() - 2);
      double std_err = sqrt(mse/sxx);
      
      // Calculate two-sided z-score for the regression coefficient
      double sim_z = beta1 / std_err;
      sim_z = std::ceil(sim_z*100000.0)/100000.0;  //round up to 5 decimal places  
      // write results in file
      data_out<<rsid<<" "<<chr<<" "<<bp<<" "<<a1<<" "<<a2<<" ";
      data_out<<std::setprecision(5)<<std::fixed<<sim_af1<<" ";
      //data_out<<std::setprecision(5)<<std::fixed<<beta1<<" ";
      //data_out<<std::setprecision(5)<<std::fixed<<beta0<<" ";
      //data_out<<std::setprecision(5)<<std::fixed<<mse<<" ";
      //data_out<<std::setprecision(5)<<std::fixed<<sxy<<" ";
      //data_out<<std::setprecision(5)<<std::fixed<<sxx<<" ";
      //data_out<<std::setprecision(5)<<std::fixed<<std_err<<" ";
      data_out<<std::setprecision(5)<<std::fixed<<sim_z<<std::endl;
    }
  }
  data_out.close(); //close output filestream
  bgzf_close(fpi);  //close reference index file
  bgzf_close(fpd);  //close reference data file
  
}



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
    std::cout<<std::endl;
    Rcpp::stop("ERROR: can't open index data file '"+index_data_file+"'");
  }
  BGZF* fpd = bgzf_open(reference_data_file.c_str(), "r");
  if(!fpd){
    std::cout<<std::endl;
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



// [[Rcpp::export]]
void test_gz_file(std::string gz_file){
  BGZF* fp = bgzf_open(gz_file.c_str(), "r");
  if(!fp){
    std::cout<<std::endl;
    Rcpp::stop("ERROR: can't open index data file '"+gz_file+"'");
  }
  std::string line;
  int last_char = BgzfGetLine(fp, line);
  std::cout<<line<<std::endl;
  bgzf_close(fp);
}


// [[Rcpp::export]]
std::string get_geno_info(int64_t fpos, 
                   std::string reference_data_file){

  BGZF* fp = bgzf_open(reference_data_file.c_str(), "r");
  if(!fp){
    std::cout<<std::endl;
    Rcpp::stop("ERROR: can't open reference data file '"+reference_data_file+"'");
  }
  std::cout<<"fpos: "<<fpos<<std::endl;
  std::string line;
  bgzf_seek(fp, fpos, SEEK_SET);
  BgzfGetLine(fp, line);
  
  bgzf_close(fp);
  return line;
}






// [[Rcpp::export]]
Rcpp::NumericVector largeval (int64_t val) {
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
      std::cout<<"Error: can't read "<<i<<"-th character"<<std::endl;
      exit(EXIT_FAILURE);
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
  std::cout<< "\r" <<percent << "%  " <<"[" << bar << "] " <<std::flush;
}


// [[Rcpp::export]]
void simulate_zscore(){
  // Initialize vectors
  std::vector<double> response(100);
  std::vector<int> predictor(100);
  
  // Fill response vector with standard normal values
  const gsl_rng_type* T;
  gsl_rng* rng;

  gsl_rng_env_setup();
  T = gsl_rng_default;
  rng = gsl_rng_alloc(T);
  
  for (int i = 0; i < 100; i++) {
    double z = gsl_ran_ugaussian(rng);
    response[i] = z;
  }
  gsl_rng_free(rng);
  
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
  std::cout << "Two-sided z-score for the regression coefficient: " << z_c1 << std::endl;
}
  