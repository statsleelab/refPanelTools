#include <Rcpp.h>
#include <string>
#include <vector>
#include <map>
#include "bgzf.h"
#include <fstream>

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
  