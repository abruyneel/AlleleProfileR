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

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(BH)]]

Rcpp::LogicalVector whichC (Rcpp::CharacterVector vectorin, Rcpp::String match) {
  Rcpp::LogicalVector out(vectorin.size());
  for (int i=0; i<vectorin.size(); i++) {
    if (vectorin[i] == match) {
      out[i] = true;
    } else {
      out[i] = false;
    }
  }

  return(out);
}

int whichCint (Rcpp::IntegerVector vectorin, int match) {
  // only for unique matches
  int out = -1;
  bool found = false;
  int temp = 0;
  while (!found && temp < vectorin.size()) {
    if (vectorin[temp] == match) {
      out = temp;
      found = true;
    } else {
      temp += 1;
    }
  }

  return(out);
}


int vecmin(IntegerVector x) {
  IntegerVector::iterator it = std::min_element(x.begin(), x.end());
  return *it;
}

int vecmax(IntegerVector x) {
  IntegerVector::iterator it = std::max_element(x.begin(), x.end());
  return *it;
}

// [[Rcpp::export]]
std::string generalpastelongC(std::vector <std::string> x, std::string sep = "") {
  std::string c;
  // concatenate inputs
  for (int i = 0; i < x.size(); i++) {
    if (i > 0) {
      c += sep;
      c += x[i];
    } else {
      c += x[i];
    }
  }

  return c;
}

// [[Rcpp::export]]
bool mathandC(LogicalVector x) {
  bool out = false;
  if (x.size() > 0) {
      out = sum(x) == x.size();
  }
  return(out);
}

// [[Rcpp::export]]
bool mathorC(LogicalVector x) {
  bool out = false;
  if (x.size() > 0) {
      out = sum(x) > 0;
  }
  return(out);
}

// [[Rcpp::export]]
bool mathintersectC(int a, int b, int x, int y) {
  bool out = false;
  if ((a <= b) & (x <= y)) {
    if ((a < x) & (b < x)) {
      out = false;
    } else if ((x < a) & (y < a)) {
      out = false;
    } else if ((b >= x) & (a <= x)) {
      out = true;
    } else if ((y >= a) & (x <= a)) {
      out = true;
    } else if ((a >= x) & (b <= y)) {
      out = true;
    } else if ((x >= a) & (y <= b)) {
      out = true;
    }
  }

  return(out);
}

// [[Rcpp::export]]
IntegerVector mathrangeC(int start, int stop, int step) {
  IntegerVector out(1);
  out[0] = start;

  while(out[out.size()-1] + step <= stop) {
    out.push_back(out[out.size()-1] + step);
  }

  return(out);
}

// [[Rcpp::export]]
std::vector < int > checknegC(std::vector < int > val) {
  // since there is no zero, deal with the negative values: add -1 to all neg or zero values
  if (val.size() > 0) {
    std::vector < int > valout (val.size());
    for (int i = 0; i < val.size(); i++) {
      if (val.at(i) <= 0) {
        valout.at(i) = val.at(i)-1;
      } else {
        valout.at(i) = val.at(i);
      }
    }
    return(valout);
  } else {
    return(val);
  }
}

int checknegCsingle(int val) {
  // since there is no zero, deal with the negative values: add -1 to all neg or zero values
  if (val <= 0) {
    return(val-1);
  } else {
    return(val);
  }
}

// [[Rcpp::export]]
bool iswithincutrangeC(std::vector < int > cutpos, int range, int pos, int length, int atgloc) {
  // parse inputs
  IntegerVector tempcutpos (cutpos.size());
  for (int i = 0; i < cutpos.size(); i++) {
    tempcutpos[i] = cutpos[i] - atgloc + 1;
  }
  bool inrange = false;

  // min/max of the cutposton
  inrange = mathintersectC(pos, pos+length-1, vecmin(tempcutpos)-range+1, vecmax(tempcutpos)+range-1);

  return(inrange);
}

// [[Rcpp::export]]
bool iswithincdsC(bool atgfound, bool stopfound, int pos, int length, int seqalignpos, int atgloc, int stoploc) {
  // parse inputs
  int stoppos = stoploc - atgloc + 1;

  // check all possible outcomes
  if (atgfound && stopfound) {
    // ideal condition, is our indel after the start and before the stop?
    return(mathintersectC(pos, pos+length-1, 0, stoppos));

  } else if (!atgfound && stopfound) {
    if (pos < stoppos) {
      return(true);
    } else {
      return(false);
    }

  } else if (atgfound && !stopfound) {
    if (pos > 0) {
      return(true);
    } else {
      return(false);
    }

  } else {
    // neither stop nor start are in the read.
    // we could be either in a UTR region
    // or somewhere in the middle or a CDS
    if ((atgloc < seqalignpos+pos) && (seqalignpos+pos < stoploc)) {
      return(true);
    } else {
      return(false);
    }

  }

}

// [[Rcpp::export]]
std::vector< std::string > cigarsubsC(std::string cigar) {

  std::vector< std::string > tmp;
  std::regex cigar_regex("([[:digit:]]+)|([[:alpha:]])");
  auto cigar_begin = std::sregex_iterator(cigar.begin(), cigar.end(), cigar_regex);
  auto cigar_end = std::sregex_iterator();

  for (std::sregex_iterator i = cigar_begin; i != cigar_end; ++i) {
    std::smatch match = *i;
    std::string match_str = match.str();
    tmp.push_back(match_str);
  }

  return tmp;
}

// [[Rcpp::export]]
std::string cigarsubstail(std::vector < std::string> cigarsubs) {
  std::string out;
  int l = cigarsubs.size();
  out = cigarsubs[l-2]+cigarsubs[l-1];
  return(out);
}

IntegerVector readposendintermediateC(CharacterVector cigar) {
  // output
  IntegerVector output;

  //
  std::regex cigar_regex("([[:digit:]]+)|([[:alpha:]])");

  // loop through the vector of cigar strings
  for (int j=0; j<cigar.size(); j++) {
    //
    std::vector< std::string > tmp;
    std::string cigart =  std::string(cigar[j]);
    auto cigar_begin = std::sregex_iterator(cigart.begin(), cigart.end(), cigar_regex);
    auto cigar_end = std::sregex_iterator();

    for (std::sregex_iterator i = cigar_begin; i != cigar_end; ++i) {
      std::smatch match = *i;
      std::string match_str = match.str();
      tmp.push_back(match_str);
    }

    if (tmp.size() == 2 && tmp.at(1) == "M") {
      output.push_back(std::stoi( tmp.at(0) ) - 1);

    } else {
      std::int32_t tempout = -1;
      int i;
      for (i=0; i<tmp.size(); i++) {
        // matches or deletions only
        if (tmp.at(i+1) == "M" || tmp.at(i+1) == "D") {
          tempout = tempout + std::stoi( tmp.at(i) );
        }
        //
        i = i+1; // odd integers only
      }

      output.push_back(tempout);

    }

  }

  return output;
}

// [[Rcpp::export]]
std::vector< std::string > codonsC(std::string seq, int start, int stop = 0) {
  // init
  std::vector< std::string > out;
  int stopv;
  if (stop == 0) {
    stopv = seq.length();
  } else {
    stopv = stop;
  }

  std::string text = seq.substr(start-1, stopv-start+1);

  if (text.length() > 2) {
    //
    int c = 0;
    while (c < text.length()) {
      out.push_back(text.substr(c, 3));
      c += 3;
    }
    //
    return(out);
  } else {
    out.push_back(text);
    return(out);
  }
}

// [[Rcpp::export]]
std::string codonselectC(std::string seq, int def, int select = 0, int shift = 0) {
  int startposselectcodon = def + shift + 3 * select - 1;
  return(seq.substr(startposselectcodon, 3));
}


bool insordelC(std::string cigar) {
  /// regex cigar
  std::vector< std::string > tmp;
  std::regex cigar_regex("(-?[[:digit:]]+)D|(-?[[:digit:]]+)I");
  auto cigar_begin = std::sregex_iterator(cigar.begin(), cigar.end(), cigar_regex);
  auto cigar_end = std::sregex_iterator();

  for (std::sregex_iterator i = cigar_begin; i != cigar_end; ++i) {
    std::smatch match = *i;
    std::string match_str = match.str();
    tmp.push_back(match_str);
  }

  return(tmp.size() > 0);
}

// [[Rcpp::export]]
NumericVector cigarsubsdeletionsC(std::vector< std::string > cigarsubs, int end, int start = 0, std::string type = "cigar") {
  // output
  NumericVector out;

  if (type == "cigar") {
    IntegerVector range = mathrangeC(1,end,2);
    for (int i=0; i < range.size(); i++) {
      if (cigarsubs.at(range.at(i)) == "D") {
        out.push_back(std::stoi(cigarsubs.at(range.at(i)-1)));
      }
    }

  } else if (type == "base") {
    //loop thorugh the cigar string, and determine what max position is.
    int maxpos = 0;
    for (int c=0; c < cigarsubs.size(); c++) {
      if (cigarsubs[c]== "H") {
        //the hard clipped seq is not present in the final read
      } else if (cigarsubs[c] == "S") {
        maxpos = maxpos + std::stoi(cigarsubs[c-1]);
      } else if (cigarsubs[c] == "D") {
        //deletion not present in seq
      } else if (cigarsubs[c] == "I") {
        maxpos = maxpos + std::stoi(cigarsubs[c-1]);
      } else if (cigarsubs[c] == "M") {
        maxpos = maxpos + std::stoi(cigarsubs[c-1]);
      }
    }

    // base pair position is c
    int c = start;
    //correct the i parameter for cigarsubs as well
    int tempstart = start;
    int itemp = -1;

    while (tempstart > 0) {
      itemp = itemp+2;
      if (cigarsubs[itemp]== "H") {
        //the hard clipped seq is not present in the final read
      } else if (cigarsubs[itemp] == "S") {
        tempstart = tempstart - std::stoi(cigarsubs[itemp-1]);
      } else if (cigarsubs[itemp] == "D") {
        //
      } else if (cigarsubs[itemp] == "I") {
        tempstart = tempstart - std::stoi(cigarsubs[itemp-1]);
      } else if (cigarsubs[itemp] == "M") {
        tempstart = tempstart - std::stoi(cigarsubs[itemp-1]);
      }
    }

    tempstart = -tempstart;
    int i = itemp;
    //note that we may be half way in the next one!

    if (i+1 == cigarsubs.size()) {
      if (cigarsubs[i] == "D") {
        out.push_back(std::stoi(cigarsubs.at(i-1)));
      } else {
        c = c + std::stoi(cigarsubs[i-1])+tempstart;
      }

    } else {
      while ((c < end) && (c < maxpos) && (i+1 < cigarsubs.size())) {
        //update
        i = i + 2;

        if (cigarsubs[i] == "D") {
          out.push_back(std::stoi(cigarsubs.at(i-1)));
        } else {
          c = c + std::stoi(cigarsubs[i-1])+tempstart;
        }

        if ((tempstart > 0) && (cigarsubs[i] != "D")) {
          //basically we're somewhere halfway in the next one.
          //we'll substract the value, but only the first time
          tempstart = 0;
        }
      }

    }

  }

  return(out);
}

// [[Rcpp::export]]
NumericVector cigarsubsinsertionsC(std::vector< std::string > cigarsubs, int end, int start = 0, std::string type = "cigar") {
  // output
  NumericVector out;

  if (type == "cigar") {
    IntegerVector range = mathrangeC(1,end,2);
    for (int i=0; i < range.size(); i++) {
      if (cigarsubs.at(range.at(i)) == "I") {
        out.push_back(std::stoi(cigarsubs.at(range.at(i)-1)));
      }
    }

  } else if (type == "base") {
    // loop through the cigar string, and determine what max position is.
    int maxpos = 0;
    for (int c=0; c < cigarsubs.size(); c++) {
      if (cigarsubs[c] == "H") {
        // the hard clipped seq is not present in the final read
      } else if (cigarsubs[c] == "S") {
        maxpos = maxpos + std::stoi(cigarsubs[c-1]);
      } else if (cigarsubs[c] == "D") {
        // deletion not present in seq
      } else if (cigarsubs[c] == "I") {
        maxpos = maxpos + std::stoi(cigarsubs[c-1]);
      } else if (cigarsubs[c] == "M") {
        maxpos = maxpos + std::stoi(cigarsubs[c-1]);
      }
    }

    // base pair position is c
    int c = start;
    // correct the i parameter for cigarsubs as well
    int tempstart = start;
    int itemp = -1;

    while (tempstart > 0) {
      itemp = itemp+2;
      if (cigarsubs[itemp]== "H") {
        // the hard clipped seq is not present in the final read
      } else if (cigarsubs[itemp] == "S") {
        tempstart = tempstart - std::stoi(cigarsubs[itemp-1]);
      } else if (cigarsubs[itemp] == "D") {
        //
      } else if (cigarsubs[itemp] == "I") {
        tempstart = tempstart - std::stoi(cigarsubs[itemp-1]);
      } else if (cigarsubs[itemp] == "M") {
        tempstart = tempstart - std::stoi(cigarsubs[itemp-1]);
      }
    }

    tempstart = -tempstart;
    int i = itemp;
    // note that we may be half way in the next one!

    // check first whether we've reached the last one already
    if (i+1 == cigarsubs.size()) {
      if (cigarsubs[i] == "I") {
        out.push_back(std::stoi(cigarsubs.at(i-1)));
        c = c + std::stoi(cigarsubs[i-1])+tempstart;

      } else if (cigarsubs[i] == "D") {
        // do not update c

      } else {
        c = c + std::stoi(cigarsubs[i-1])+tempstart;

      }

    } else {
      while ((c < end) && (c < maxpos) && (i+1 < cigarsubs.size())) {
        // update
        i = i + 2;

        if (cigarsubs[i] == "I") {
          out.push_back(std::stoi(cigarsubs.at(i-1)));
          c = c + std::stoi(cigarsubs[i-1])+tempstart;

        } else if (cigarsubs[i] == "D") {
          // do not update c

        } else {
          c = c + std::stoi(cigarsubs[i-1])+tempstart;

        }

        if ((tempstart > 0) && (cigarsubs[i] != "D")) {
          // basically we're somewhere halfway in the next one.
          // we'll substract the value, but only the first time
          tempstart = 0;
        }
      }
    }

  }

  return(out);
}

// [[Rcpp::export]]
NumericVector cigarsubsSHC(std::vector< std::string > cigarsubs, int end, int start) {
  // output
  NumericVector out;

  IntegerVector range = mathrangeC(1,end,2);
  for (int i=0; i < range.size(); i++) {
    if (cigarsubs.at(range.at(i)) == "S" || cigarsubs.at(range.at(i)) == "H" ) {
      out.push_back(std::stoi(cigarsubs.at(range.at(i)-1)));
    }
  }

  return(out);
}


// [[Rcpp::export]]
bool iswithininsertionC(std::vector< std::string > cigarsubs, int pos) {
  // output
  bool out = false;

  // loop through the cigar string
  int currentpos = 0;
  for (int c=0; c < cigarsubs.size(); c++) {
    if (cigarsubs[c] == "S" || cigarsubs[c] == "M") {
      // update pos
      currentpos = currentpos + std::stoi(cigarsubs[c-1]);
    } else if ( cigarsubs[c] == "I") {
      if ( currentpos < pos && currentpos + std::stoi(cigarsubs[c-1]) > pos) {
        out = true;
      }

      // update pos
      currentpos = currentpos + std::stoi(cigarsubs[c-1]);
    }
  }

  return(out);
}


