#include <Rcpp.h>
#include <string>
#include <vector>
#include <map>
#include "bgzf.h"
#include <fstream>
#include <random>

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
  //for(int i=0; i<pop_vec_input.size(); i++) {std::cout<<pop_vec_input[i]<<std::endl;}
  
  
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
  
  //std::cout<<num_pops<<std::endl;
  //for(int i=0; i<ref_pop_vec.size(); i++) {std::cout<<ref_pop_vec[i]<<std::endl;}
  //for(int i=0; i<ref_pop_vec.size(); i++) {std::cout<<ref_pop_size_vec[i]<<std::endl;}
  //for(int i=0; i<ref_pop_vec.size(); i++) {std::cout<<ref_sup_pop_vec[i]<<std::endl;}
  
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
  
  //for(int i=0; i<pop_flag_vec.size(); i++) {std::cout<<pop_flag_vec[i]<<std::endl;}
  

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
void simulate_af1(int chr_num,
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
    //std::cout<<pop<<std::endl;
  }
  //for(int i=0; i<pop_vec_input.size(); i++) {std::cout<<pop_vec_input[i]<<std::endl;}
  
  
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
  
  //std::cout<<num_pops<<std::endl;
  //for(int i=0; i<ref_pop_vec.size(); i++) {std::cout<<ref_pop_vec[i]<<std::endl;}
  //for(int i=0; i<ref_pop_vec.size(); i++) {std::cout<<ref_pop_size_vec[i]<<std::endl;}
  //for(int i=0; i<ref_pop_vec.size(); i++) {std::cout<<ref_sup_pop_vec[i]<<std::endl;}
  
  // init pop_flag vector
  std::vector<int> pop_flag_vec;
  int total_num_subj=0;
  for(int i=0; i<num_pops; i++){
    std::string pop = ref_pop_vec[i];
    if(std::find(pop_vec_input.begin(), pop_vec_input.end(), pop)!=pop_vec_input.end()){ //if pop is found in pop_vec_input
      pop_flag_vec.push_back(1);
      total_num_subj += ref_pop_size_vec[i];
    } else {
      pop_flag_vec.push_back(0);
    }
  }
  
  //for(int i=0; i<pop_flag_vec.size(); i++) {std::cout<<pop_flag_vec[i]<<std::endl;}
  
  
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
  
  std::random_device rd;
  std::mt19937 gen(rd());
  
  while(true){
    last_char = BgzfGetLine(fpi, index_line);
    if(last_char == -1) //EOF
      break;
    
    std::istringstream buffer(index_line);
    buffer >> rsid >> chr >> bp >> a1 >> a2 >> af1ref >> fpos;
    
    double sim_af1=0.0;
    double allele_counter=0.0;
    if(chr==chr_num){
      bgzf_seek(fpd, fpos, SEEK_SET);
      last_char = BgzfGetLine(fpd, data_line);      
      if(last_char == -1) //EOF
        break;
      
      std::istringstream data_buffer(data_line);
      for(int k=0; k<num_pops; k++){
        std::string geno_str;
        data_buffer >> geno_str;
        if(pop_flag_vec[k]) {
          for(int j=0; j<num_sim_vec.size(); j++){
            std::uniform_int_distribution<> dis(0, ref_pop_size_vec[k] - 1);
            int ran_index = dis(gen);
            allele_counter += (double)(geno_str[ran_index]-'0');
          }
        }
      }
      // compute af1 of simulated genotype
      sim_af1 = allele_counter/(2*total_num_subj);
      sim_af1 = std::ceil(sim_af1*100000.0)/100000.0;  //round up to 5 decimal places
      // write results in file
      data_out<<rsid<<" "<<chr<<" "<<bp<<" "<<a1<<" "<<a2<<" ";
      data_out<<std::setprecision(5)<<std::fixed<<af1ref<<std::endl;
      data_out<<std::endl;
    }
  }
  data_out.close();
  bgzf_close(fpi);
  bgzf_close(fpd);
  
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
  