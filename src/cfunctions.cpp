// Copyright (C) 2017-2018  Arne A.N. Bruyneel
//
// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation; either version 2
// of the License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.


#include <Rcpp.h>
#include <iostream>
#include <string>
#include <regex>
#include <cstdio>
#include <memory>
#include <stdexcept>
#include <array>
#include <cstdarg>

using namespace Rcpp;
#include "chelper.h"
#include "bioio.h" // fasta file lib

#include <boost/algorithm/string.hpp>

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(BH)]]

// functions

// [[Rcpp::export]]
std::string binassignC(int pos, int posend, NumericVector y) {
  // extract vars
  int lower = y[0];
  int upper = y[1];

  // output bins for each condition
  if (pos < lower && posend > upper) {
    return "D";
  } else if (pos < lower && posend > lower) {
    return "B";
  } else if (pos < upper && posend > upper) {
    return "C";
  } else {
    return "A";
  }
}


// [[Rcpp::export]]
DataFrame binassignChimericC(DataFrame df, NumericVector y) {
  // extract vars for range determination
  int lower = y[0];
  int upper = y[1];
  int pcrlower = y[2];
  int pcrupper = y[3];

  // access the columns of the data.frame
  IntegerVector pos = df["pos"];
  IntegerVector posend = df["posend"];
  CharacterVector qn = df["qname"];
  CharacterVector bin (qn.size());
  LogicalVector out_lower (qn.size());
  LogicalVector out_upper (qn.size());

  // preliminary determination
  LogicalVector select = (pos < lower) & (posend > upper);

  // unique names of reads, defined chimeras
  CharacterVector temp_unique_chimer = unique(qn);

  // loop through each pair of chimeric reads
  for (int j=0; j<temp_unique_chimer.size(); j++) {
    LogicalVector chimrows = whichC(qn, temp_unique_chimer[j]);
    // subset
    LogicalVector temp_data_select = select[chimrows];
    IntegerVector temp_data_pos = pos[chimrows];
    IntegerVector temp_data_posend = posend[chimrows];

    if (sum(temp_data_select) == temp_data_select.size()) {
      // select these chimeras
      select[chimrows] = true;
      bin[chimrows] = "D";

    } else {
      // more complex cenarios
      if (temp_data_select.size() == 2) {
        if ((temp_data_pos[0] < lower) & (temp_data_posend[1] > upper) ) {
          select[chimrows] = true;
          out_lower[chimrows] = LogicalVector {true, false};
          out_upper[chimrows] = LogicalVector {false, true};
          bin[chimrows] = "D";

        } else if ((temp_data_pos[1] < lower) & (temp_data_posend[0] > upper) ) {
          select[chimrows] = true;
          out_lower[chimrows] = LogicalVector {false, true};
          out_upper[chimrows] = LogicalVector {true, false};
          bin[chimrows] = "D";

        } else {
          // simplify this by looking at the min of the pos and the max of the posends
          int temp_chim_start = vecmin(temp_data_pos);
          int temp_chim_stop = vecmax(temp_data_posend);

          // bins
          if (temp_chim_start < lower & temp_chim_stop > lower) {
            bin[chimrows] = "B";
          } else if (temp_chim_start < upper & temp_chim_stop > upper) {
            bin[chimrows] = "C";
          } else if (temp_chim_start > pcrlower & temp_chim_stop < pcrupper) {
            bin[chimrows] = "A";
          }
        }

      }

    }


  } // loop unique

  // return a new data frame
  return DataFrame::create(_["select"]= select, _["bin"]= bin, _["lower"]= out_lower, _["upper"]= out_upper);
}

// [[Rcpp::export]]
std::vector< std::string > inssitesC(std::string allele) {

  std::vector< std::string > tmp;
  std::vector< std::string > tmpout;

  std::regex allele_regex("(-?[[:digit:]]+)_(-?[[:digit:]]+)ins([[:alpha:]]+)");
  std::regex suballele_regex("(-?[[:digit:]]+)|ins([[:alpha:]]+)");
  auto allele_begin = std::sregex_iterator(allele.begin(), allele.end(), allele_regex);
  auto allele_end = std::sregex_iterator();

  for (std::sregex_iterator i = allele_begin; i != allele_end; ++i) {
    std::smatch match = *i;
    std::string match_str = match.str();
    tmp.push_back(match_str);
  }

  if (tmp.size() > 0) {
    for (int j=0; j < tmp.size(); j++) {
      auto suballele_begin = std::sregex_iterator(tmp[j].begin(), tmp[j].end(), suballele_regex);
      auto suballele_end = std::sregex_iterator();

      for (std::sregex_iterator k = suballele_begin; k != suballele_end; ++k) {
        std::smatch match = *k;
        std::string match_str;

        if (match.str().length() > 3) {
          if (match.str().substr(0,3) == "ins") {
            match_str = match.str().substr(3,match.str().length()-3);
          } else {
            match_str = match.str();
          }

        } else {
          match_str = match.str();
        }

        //
        tmpout.push_back(match_str);
      }

    }

  }

  return tmpout;
}

// [[Rcpp::export]]
IntegerVector readposendC(DataFrame df) {
  // parse df
  IntegerVector pos = df["pos"];
  CharacterVector cigar = df["cigar"];
  // prep output
  IntegerVector output(pos.size());

  //
  output = pos + readposendintermediateC(cigar);

  return(output);
}

// [[Rcpp::export]]
std::string seqinsertedC(std::string seq, std::vector< std::string > cigarsubs, int cigarsubpos) {
  //initialize variables
  int seq_start = 0;
  cigarsubpos = cigarsubpos-1; // R vs C++ index

  //loop thorugh the cigar string, up to the sub before the start of the one of interest.
  for (int i=0; i < cigarsubpos; i++) {
    if (cigarsubs.at(i) == "H") {
      //the hard clipped seq is not present in the final read
    } else if (cigarsubs.at(i) == "S") {
      seq_start = seq_start + std::stoi(cigarsubs.at(i-1));
    } else if (cigarsubs.at(i) == "D") {
      //deletion hence the next matched bases will have the same position as the deleted ones, if they would have not been deleted.
    } else if (cigarsubs.at(i) == "I") {
      seq_start = seq_start + std::stoi(cigarsubs.at(i-1));
    } else if (cigarsubs.at(i) == "M") {
      seq_start = seq_start + std::stoi(cigarsubs.at(i-1));
    }
  }

  // subtring of seq with the start and length
  return(seq.substr(seq_start, std::stoi(cigarsubs.at(cigarsubpos-1))));
}

// [[Rcpp::export]]
std::vector< bool > bamQCC(std::vector< std::int64_t > flag) {
  // output
  std::vector< bool > output;

  // loop through the vector of flags
  // bit math for flag determination
  for (int j=0; j<flag.size(); j++) {
    output.push_back(((flag.at(j) & ( 1 << 2 )) >> 2) || ((flag.at(j) & ( 1 << 9 )) >> 9) || ((flag.at(j) & ( 1 << 10 )) >> 10));
  }

  return output;
}


// [[Rcpp::export]]
bool mergedSanityCheckC(DataFrame df) {
  // extract params
  CharacterVector cigars = df["cigar"];
  CharacterVector seqs = df["seq"];

  // set output
  bool out = false;

  // there should only be one row in the df
  // double check, just in case
  if (seqs.size() != 1) {
    out = false;

  } else {
    // ok to proceed
    std::string cigar = std::string(cigars[0]);
    std::string seq = std::string(seqs[0]);
    std::vector < std::string > cigarsubs = cigarsubsC(cigar);
    int seqlength = seq.length();

    // set defaults
    bool noerror = true;
    int templength = 0;

    //loop thorugh the cigar string, and determine what the seq length should be based on the cigar string.
    for (int c = 0; c < cigarsubs.size(); c++) {
      if (cigarsubs.at(c)== "H") {
        //the hard clipped seq is not present in the final read

      } else if (cigarsubs.at(c) == "S") {
        templength = templength + std::stoi(cigarsubs[c-1]);

      } else if (cigarsubs.at(c) == "D") {
        //deletion not present in seq
        if (std::stoi(cigarsubs.at(c-1)) == 0) {
          noerror = false;
        }

      } else if (cigarsubs.at(c) == "I") {
        templength = templength + std::stoi(cigarsubs.at(c-1));

      } else if (cigarsubs.at(c) == "M") {
        templength = templength + std::stoi(cigarsubs.at(c-1));

      }
    }

    out = ((seqlength == templength) && noerror);

  } // if rows in df

  //return
  return(out);

}


// [[Rcpp::export]]
std::string getreferenceC(std::string file, std::string chr, int lower, int upper) {
  //
  const auto index = bioio::read_fasta_index(file+".fai");
  const std::string contig {chr};

  size_t chunk_begin {size_t(lower-1)};
  size_t chunk_size  {size_t(upper-lower+1)};

  const auto chunk = bioio::read_fasta_contig(file, index.at(contig), chunk_begin, chunk_size);

  return(chunk);
}

// [[Rcpp::export]]
std::vector < int > checkatgC(DataFrame df, DataFrame gene) {
  //get the ATG position with respect to the start of the alignment
  //are there any indels before the putative start codon?

  //init
  IntegerVector atgv = gene["ATG"];
  CharacterVector starttypev = gene["StartType"];
  IntegerVector posv = df["pos"];
  CharacterVector cigarv = df["cigar"];
  CharacterVector seqv = df["seq"];

  int atgloc = atgv[0];
  std::string type = std::string(starttypev[0]);
  bool output_atg = false;
  int patgpos = atgloc - posv[0] + 1; //initial guess: this assumes that there no indels before the start
  int shift = 0;
  int current = 1;
  std::string cigar = std::string (cigarv[0]);
  std::vector< std::string > cigarsubs = cigarsubsC(cigar);
  bool found = false;
  int atgpos = 0;
  std::string seq = std::string(seqv[0]);

  if (patgpos < 1) {
    // if it's less than the first position, then basically the start codon preseeds the read.
    // we can determine this negative value.
    found = false;

  } else {
    if (cigarsubs.size() > 2) {
      int posexplained = 0; // position upto which the indels etc are explained

      while ((found == false) && (current < cigarsubs.size())) {
        if (cigarsubs[current] == "D") {
          posexplained += std::stoi(cigarsubs[current-1]);

        } else if (cigarsubs[current] == "I") {
          shift += std::stoi(cigarsubs[current-1]);

        } else if (cigarsubs[current] == "M") {
          posexplained += std::stoi(cigarsubs[current-1]);

        }

        if (posexplained >= patgpos) {
          found = true;
        }

        // update
        current = current+2;
      }

    } else if ((cigarsubs.size() == 2) && (cigarsubs[1] == "M")) {
      if (patgpos < std::stoi(cigarsubs[0])) {
        found = true;
      }
    }
  }

  if (found) {
    // atg position
    atgpos = patgpos + shift;

    // is there a deletion or insertion
    if (insordelC(cigar)) {
      // so there is at least one deletion/insertion
      // where is it located with respect to the ATG location
      for (int d=3; d < cigarsubs.size(); d++) {
        if ((cigarsubs[d] == "D") | (cigarsubs[d] == "I")) {
          // tmp
          int posstart, posend;

          // calc start
          IntegerVector rangestart = mathrangeC(0,d-3,2);
          int rangecigarsubssumstart = 0;
          for (int f=0; f < rangestart.size(); f++) {
            rangecigarsubssumstart += std::stoi(cigarsubs[rangestart[f]]);
          }
          posstart = rangecigarsubssumstart + 1 - sum(cigarsubsinsertionsC(cigarsubs, d-2, 0, "cigar"));

          if (cigarsubs[d] == "I") {
            posend = posstart;

          } else {
            // calc stop
            IntegerVector rangestop = mathrangeC(0,d-1,2);
            int rangecigarsubssumstop = 0;
            for (int f=0; f < rangestop.size(); f++) {
              rangecigarsubssumstop += std::stoi(cigarsubs[rangestop[f]]);
            }
            //
            posend = rangecigarsubssumstop - sum(cigarsubsinsertionsC(cigarsubs, d-2, 0, "cigar"));
          }

          if (mathintersectC(posstart, posend, atgpos, atgpos+2)) {
            output_atg = true;
          }
        }
      }
    }

  } else {
    // return default
    atgpos = patgpos + shift;
    // this works if the ATG is withing the gap or if it is before
  }

  // prep output
  std::vector < int > output(3);
  output[0] = output_atg;
  output[1] = atgpos;
  output[2] = found;

  return(output);
}


// [[Rcpp::export]]
std::vector < int > checkstopC(DataFrame df, DataFrame gene) {
  // get the stop position with respect to the start of the alignment
  // are there any indels before the putative stop codon?

  //init
  IntegerVector stopcv = gene["StopCodon"];
  CharacterVector starttypev = gene["StartType"];
  IntegerVector posv = df["pos"];
  CharacterVector cigarv = df["cigar"];
  CharacterVector seqv = df["seq"];

  int stoploc = stopcv[0];
  bool output_stop = false;
  int pstoppos = stoploc - posv[0] + 1; // initial guess: this assumes that there no indels before the stop codon
  int shift = 0;
  int current = 1;
  int stoppos = 0;
  std::string cigar = std::string (cigarv[0]);
  std::vector< std::string > cigarsubs = cigarsubsC(cigar);
  bool found = false;
  std::string seq = std::string(seqv[0]);


  if (cigarsubs.size() > 2) {
    // position upto which the indels etc are explained
    int posexplained = 0;

    while ((found == false) && (current < cigarsubs.size())) {

      if (cigarsubs[current] == "D") {
        posexplained += std::stoi(cigarsubs[current-1]);

      } else if (cigarsubs[current] == "I") {
        shift += std::stoi(cigarsubs[current-1]);

      } else if (cigarsubs[current] == "M") {
        posexplained += std::stoi(cigarsubs[current-1]);

      }

      if (posexplained >= pstoppos) {
        found = true;
      }

      // update
      current = current + 2;
    }
  } else if ((cigarsubs.size() == 2) && (cigarsubs[1] == "M")) {
    if (pstoppos < std::stoi(cigarsubs[0])) {
      found = true;
    }
  }

  // if found, in case of leading and trail or no stop in reach, there is no need to find one
  if (found) {
    // stop position
    stoppos = pstoppos + shift;

    if (insordelC(cigar)) {
      // so there is at least one deletion/insertion
      // where is it located with respect to the stop location
      for (int d=3; d < cigarsubs.size(); d++) {
        if ((cigarsubs[d] == "D") | (cigarsubs[d] == "I")) {
          // tmp
          int posstart, posend;

          // calc posstart
          IntegerVector rangestart = mathrangeC(0,d-3,2);
          int rangecigarsubssumstart = 0;
          for (int f=0; f < rangestart.size(); f++) {
            rangecigarsubssumstart += std::stoi(cigarsubs[rangestart[f]]);
          }
          posstart = rangecigarsubssumstart + 1 - sum(cigarsubsinsertionsC(cigarsubs, d-2, 0, "cigar"));

          if (cigarsubs[d] == "I") {
            posend = posstart;

          } else {
            // calc stop
            IntegerVector rangestop = mathrangeC(0,d-1,2);
            int rangecigarsubssumstop = 0;
            for (int f=0; f < rangestop.size(); f++) {
              rangecigarsubssumstop += std::stoi(cigarsubs[rangestop[f]]);
            }
            //
            posend = rangecigarsubssumstop - sum(cigarsubsinsertionsC(cigarsubs, d-2, 0, "cigar"));
          }

          if (mathintersectC(posstart,posend,stoppos,stoppos+2)) {
            output_stop = true;
          }
        }
      }
    }
  } else {
    // return NULL or NA
    stoppos = 0;
    output_stop = false;
  }

  // prep output
  std::vector < int > output(3);
  output[0] = output_stop;
  output[1] = stoppos;
  output[2] = found;

  return(output);
}

// [[Rcpp::export]]
std::string mutreferenceC(DataFrame df, std::vector < std::string > cigarsubs, DataFrame gene, std::string index) {
  //init
  IntegerVector posv = df["pos"];
  IntegerVector posendv = df["posend"];
  CharacterVector cigarv = df["cigar"];
  CharacterVector seqv = df["seq"];
  IntegerVector atgv = gene["ATG"];
  IntegerVector stopcv = gene["StopCodon"];
  CharacterVector starttypev = gene["StartType"];
  CharacterVector stoptypev = gene["StopType"];
  CharacterVector chrv = gene["Chr"];
  IntegerVector pcrstartv = gene["PCRStart"];
  IntegerVector pcrstopv = gene["PCRStop"];

  // output
  std::string cds;

  //extract info
  int posendobj = posendv[0];
  int pos = posv[0];
  std::string chr = std::string (chrv[0]);
  std::string starttype = std::string (starttypev[0]);
  std::string stoptype = std::string (stoptypev[0]);
  std::string cigar = std::string (cigarv[0]);

  // the sequence may have some clipping at the start of end, check the cigar and make certain this is removed, when doing stop/stop position inputation.
  std::string prelimseq = std::string (seqv[0]);
  int start, seqlength;

  // retrieve the read sequence
  // sequence clipped of S
  // prelimseq might start with 'S' in the cigar
  if (cigarsubs[1] == "S") {
    start = std::stoi(cigarsubs[0]);
    cigarsubs.erase(cigarsubs.begin());
    cigarsubs.erase(cigarsubs.begin());
  } else {
    start = 0;
  }

  if (cigarsubs[1] == "H") {
    cigarsubs.erase(cigarsubs.begin());
    cigarsubs.erase(cigarsubs.begin());
  }

  // the prelimseq might also end in a 'S'
  if (cigarsubs.back() == "S") {
    seqlength = prelimseq.length() - std::stoi(cigarsubs[cigarsubs.size()-2]) - start;
    cigarsubs.pop_back();
    cigarsubs.pop_back();
  } else {
    seqlength = prelimseq.length();
  }

  if (cigarsubs.back() == "H") {
    cigarsubs.pop_back();
    cigarsubs.pop_back();
  }

  // seq of the read
  std::string seq = prelimseq.substr(start, seqlength);

  // assess whether there are any indels:
  if (insordelC(cigar)) {
    std::vector < int > positions;
    std::vector < std::string > insertions;
    // so there is at least one deletion/insertion
    // where is it located with respect to refernce
    for (int d = 3; d < cigarsubs.size(); d++) {
      if ((cigarsubs[d] == "D") | (cigarsubs[d] == "I")) {
        // locations are common
        int posstart, posend;
        // calc posstart
        IntegerVector rangestart = mathrangeC(0,d-3,2);
        int rangecigarsubssumstart = 0;
        for (int f=0; f < rangestart.size(); f++) {
          rangecigarsubssumstart += std::stoi(cigarsubs[rangestart[f]]);
        }
        posstart = posv[0] - 1 + rangecigarsubssumstart + 1 - sum(cigarsubsinsertionsC(cigarsubs, d-2, 0, "cigar"));

        if ((positions.size() == 0) && (posstart > posv[0])) {
          positions.push_back(posv[0]);
        }

        // identify inserted bases
        if (cigarsubs[d] == "D") {
          // deleted bases should be excluded from the reference
          // a deletion lacks an insertion of bases
          insertions.push_back("");
          // calc stop
          IntegerVector rangestop = mathrangeC(0,d-1,2);
          int rangecigarsubssumstop = 0;
          for (int f=0; f < rangestop.size(); f++) {
            rangecigarsubssumstop += std::stoi(cigarsubs[rangestop[f]]);
          }
          //
          posend = posv[0] - 1 + rangecigarsubssumstop - sum(cigarsubsinsertionsC(cigarsubs, d-2, 0, "cigar"));

        } else if (cigarsubs[d] == "I") {
          // inserted bases need to be extracted from the obj and introduced into the reference
          insertions.push_back(seqinsertedC(seq, cigarsubs, d+1)); // +1 due to R C++ difference in indices
          posend = posstart-1;

        }

        // save start and stop positions
        positions.push_back(posstart-1);
        if (posendobj < posend + 1) {
          positions.push_back(posendobj);
        } else {
          positions.push_back(posend + 1);
        }

      }
    }

    if ((positions.back() < posendobj) || ((positions.back() <= posendobj) && (cigarsubstail(cigarsubs) == "1M"))) {
      positions.push_back(posendobj);
      insertions.push_back("");
    }

    // parse the mutated reference coding sequence
    IntegerVector parseids = mathrangeC(1, positions.size(), 2);

    for (int i=0; i < parseids.size(); i++) {
      std::string additionalseq = getreferenceC(index, chr, positions[parseids[i]-1], positions[parseids[i]]);
      cds = cds+additionalseq+insertions[parseids[i]/2];
    }

  } else {
    // matches only
    // the unmutate cds is:
    cds = getreferenceC(index, chr, pos, posendobj);
  }

  // return
  return(cds);

}

// [[Rcpp::export]]
List snpsC(DataFrame df, DataFrame gene, std::string index) {
  //init
  IntegerVector posv = df["pos"];
  IntegerVector posendv = df["posend"];
  CharacterVector cigarv = df["cigar"];
  CharacterVector seqv = df["seq"];
  IntegerVector atgv = gene["ATG"];
  std::string cigar = std::string (cigarv[0]);
  int pos = posv[0];
  int posendobj = posendv[0];
  std::vector <std::string> cigarsubs = cigarsubsC(cigar);
  std::string prelimseq = std::string (seqv[0]);

  // temps
  int start, seqlength;
  std::string seq;
  std::vector <int> prelimsnppositions;

  // output
  std::vector <int> positions;
  std::vector <int> shifts;
  std::vector <char> altbase;
  std::vector <std::string> snppositions;
  int snpcount = 0;

  // retrieve the read sequence
  // sequence clipped of S
  // prelimseq might start with 'S' in the cigar
  if (cigarsubs[1] == "S") {
    start = std::stoi(cigarsubs[0]);
    cigarsubs.erase(cigarsubs.begin());
    cigarsubs.erase(cigarsubs.begin());
  } else {
    start = 0;
  }

  if (cigarsubs[1] == "H") {
    cigarsubs.erase(cigarsubs.begin());
    cigarsubs.erase(cigarsubs.begin());
  }

  // the prelimseq might also end in a 'S'
  if (cigarsubs.back() == "S") {
    seqlength = prelimseq.length() - std::stoi(cigarsubs[cigarsubs.size()-2]) - start;
    cigarsubs.pop_back();
    cigarsubs.pop_back();
  } else {
    seqlength = prelimseq.length();
  }

  if (cigarsubs.back() == "H") {
    cigarsubs.pop_back();
    cigarsubs.pop_back();
  }

  // seq of the read
  seq = prelimseq.substr(start, seqlength);

  // mutated reference
  std::string mergedcigarsubs = generalpastelongC(cigarsubs, "");
  DataFrame dummyobj = DataFrame::create(_["pos"] = pos, _["seq"] = seq, _["cigar"] = mergedcigarsubs);
  IntegerVector dummyobjposend = readposendC(dummyobj);
  dummyobj = DataFrame::create(_["pos"] = pos, _["seq"] = seq, _["cigar"] = mergedcigarsubs, _["posend"] = dummyobjposend);
  std::string mutref = mutreferenceC(dummyobj, cigarsubs, gene, index);

  // error handling
  if (mutref.length() != seq.length()) {
    Rcout << "mutref : " << mutref.length() << " : " << mergedcigarsubs << std::endl;
    Rcout << "mutref : " << mutref << std::endl;
    Rcout << "seqref : " << seq.length() << std::endl;
    Rcout << "seqref : " << seq << std::endl;
    Rcout << "prelimseq : " << prelimseq << std::endl;
    throw "Fatal error: SNP determination failure: unequal lengths!";
  }

  // these positions will have to reflect the reference, at this stage it ignores deletions/insertions
  for (int k=0; k < mutref.length(); k++) {
    if (mutref[k] != seq[k]) {
      prelimsnppositions.push_back(pos - 1 + k);
      altbase.push_back(seq[k]);
    }
  }

  // determine correction factors
  for (int d = 1; d < cigarsubs.size(); d++) {
    if (cigarsubs[d] == "M") {
      if (positions.size() == 0) {
        // so basically if empty
        positions.push_back(pos);
        shifts.push_back(0);
      }

    } else if ((cigarsubs[d] == "D") | (cigarsubs[d] == "I")) {
      // locations are common
      int posstart, posend;
      int indellength;

      // calc posstart
      IntegerVector rangestart = mathrangeC(0,d-3,2);
      int rangecigarsubssumstart = 0;
      for (int f=0; f < rangestart.size(); f++) {
        rangecigarsubssumstart += std::stoi(cigarsubs[rangestart[f]]);
      }
      posstart = pos + rangecigarsubssumstart - sum(cigarsubsdeletionsC(cigarsubs, d-2, 0, "cigar"));
      // don't substract insertions, but deletions because they are included in the mut reference
      indellength = std::stoi(cigarsubs[d-1]);

      // identify the shifts
      if (cigarsubs[d] == "D") {
        // deleted bases should be excluded from the reference
        // a deletion lacks an insertion of bases
        shifts.push_back(shifts.back() + indellength);

      } else if (cigarsubs[d] == "I") {
        shifts.push_back(shifts.back() - indellength);

      }

      // save start position of current indel
      if (posendobj < posstart) {
        positions.push_back(posendobj);
      } else {
        positions.push_back(posstart);
      }

    }
  }

  //
  if (positions.back() < posendobj) {
    positions.push_back(posendobj);
  }

  // correct positions
  for (int i=1; i < positions.size(); i++) {
    for (int l=0; l < prelimsnppositions.size(); l++) {
      if ((prelimsnppositions[l] >= positions[i-1]) && (prelimsnppositions[l] < positions[i])) {
        prelimsnppositions[l] = prelimsnppositions[l] + shifts[i-1];
      }
    }
  }

  // convert to relative positions
  // make strings from snppositions
  for (int j=0; j < prelimsnppositions.size(); j++) {
    snppositions.push_back(std::to_string(checknegCsingle(prelimsnppositions[j] + 2 - atgv[0])));
  }

  // if no snps
  if (snppositions.size() == 0) {
    snpcount = 0;
    snppositions.push_back("");
    altbase.push_back('\0');
  } else {
    snpcount = snppositions.size();
  }

  // return list(count, snppositions, altbase)
  return Rcpp::List::create(Rcpp::Named("count") = snpcount,
                            Rcpp::Named("snppos") = snppositions,
                            Rcpp::Named("altbase") = altbase);
}

// [[Rcpp::export]]
std::string variantnameC(List dels, List ins, List snp) {
  // deal with List inputs
  bool bool_dels = dels[0];
  std::vector <std::string> delspos = dels[1];
  bool bool_ins = ins[0];
  std::vector <std::string> inspos = ins[1];
  std::vector <std::string> insseq = ins[2];
  int no_snps = snp[0];
  std::vector <std::string> snpspos = snp[1];
  std::vector <std::string> snpsseq = snp[2];

  // output
  std::string out;

  if (bool_dels | bool_ins) {
    // process output
    out += "exon";

    if (bool_dels) {
      for (int i=0; i < delspos.size(); i++) {
        if (delspos[i].length() > 0) {
          out += "."+delspos[i]+"del";
        }
      }
    }

    if (bool_ins) {
      for (int i=0; i < inspos.size(); i++) {
        if (inspos[i].length() > 0) {
          out += "."+inspos[i]+"ins"+insseq[i];
        }
      }
    }

  } else {
    out = "WT";
  }

  // snps
  if (no_snps > 0) {
    for (int i=0; i < snpspos.size(); i++) {
      if (snpsseq[i].length() == 1) {
        out += ".SNP("+snpspos[i]+":"+snpsseq[i]+")";
      }
    }
  }

  return(out);
}



