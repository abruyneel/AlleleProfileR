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


// [[Rcpp::export]]
DataFrame readdeptc(DataFrame df, IntegerVector range) {
  //init
  CharacterVector cigarv = df["cigar"];
  IntegerVector posv = df["pos"];
  IntegerVector posendv = df["posend"];

  // output init
  IntegerVector counts (range.size());

  // for each read
  for (int i=0; i < cigarv.size(); i++) {
    std::vector< std::string > tempcigarsubs = cigarsubsC(std::string(cigarv[i]));

    if ((tempcigarsubs.size() == 2) && (tempcigarsubs[1] == "M")) {
      // matches only
      for (int j=posv[i]; j < posendv[i]+1; j++) {
        counts[whichCint(range, j)] += 1;
      }

    } else {
      // more complex, read contains insertions and deletions
      for (int g=1; g < tempcigarsubs.size(); g++) {
        // current is a match
        if ((g == 1) && (tempcigarsubs[1] == "M")) {
          // matches only
          for (int j=posv[i]; j < posv[i]+std::stoi(tempcigarsubs[0]); j++) {
            counts[whichCint(range, j)] += 1;
          }

        } else if (tempcigarsubs[g] == "M"){
          int posstart = 0;
          int posend = 0;
          // start position
          IntegerVector rangestart = mathrangeC(0,g-2,2);
          int rangecigarsubssumstart = 0;
          for (int f=0; f < rangestart.size(); f++) {
            rangecigarsubssumstart += std::stoi(tempcigarsubs[rangestart[f]]);
          }
          posstart = rangecigarsubssumstart - sum(cigarsubsinsertionsC(tempcigarsubs, g-2, 0, "cigar")) - sum(cigarsubsSHC(tempcigarsubs,g-2));
          // end position
          posend = posstart + std::stoi(tempcigarsubs[g-1]);
          // counts
          for (int j=posv[i]+posstart; j < (posv[i]+posend); j++) {
            counts[whichCint(range, j)] += 1;
          }

        }

      }
    }

  }

  // return
  return(DataFrame::create(_["pos"] = range, _["counts"] = counts));
}

// [[Rcpp::export]]
DataFrame seglogoC(DataFrame df, DataFrame frame) {
  // init bamtable
  CharacterVector cigarv = df["unifiedcigar"];
  IntegerVector unifiedposv = df["unifiedpos"];
  CharacterVector unifiedseqv = df["unifiedseq"];

  // init (output)frame
  IntegerVector absposv = frame["abspos"];
  IntegerVector relposv = frame["relpos"];
  CharacterVector refv = frame["ref"];
  IntegerVector Av = frame["A"];
  IntegerVector Tv = frame["T"];
  IntegerVector Gv = frame["G"];
  IntegerVector Cv = frame["C"];
  IntegerVector Nv = frame["N"];

  // process
  for (int i=0; i < cigarv.size(); i++) {
    // for each read
    std::vector< std::string > tempcigarsubs = cigarsubsC(std::string(cigarv[i]));
    std::string currentseq = std::string(unifiedseqv[i]);

    // LOOP through positions and cigar
    int currentpos = unifiedposv[i];
    int base = 0;

    // loop thorugh the cigar string
    for (int c=1; c < tempcigarsubs.size(); c++) {
      if (tempcigarsubs[c]== "H") {
        // the hard clipped seq is not present in the final read, hence ignore

      } else if (tempcigarsubs[c] == "S") {
        // # Starting 'S' is removed in the main processing code, hence do nothing, as it will be at the end, if it isn't at the end: terminate
        if (c != (tempcigarsubs.size()-1)) {
          throw "Fatal error: error S";
        }

      } else if (tempcigarsubs[c] == "D") {
        // deletion not present in seq
        // set these as 'N' in the outout
        int currentlength = std::stoi(tempcigarsubs[c-1]);
        if (currentlength >= 1) {
          for (int ll=0; ll < currentlength; ll++) {
            // is within range?
            int posupdate = whichCint(absposv, currentpos);
            if (posupdate != -1) {
              // update output
              Nv[whichCint(absposv, currentpos)] += 1;
            }

            // update currentpos and but NOT base
            currentpos += 1;
          }
        }

      } else if (tempcigarsubs[c] == "I") {
        // these positions do not occur in the reference, hence ignore
        // update base but NOT currentpos
        base += std::stoi(tempcigarsubs[c-1]);

      } else if (tempcigarsubs[c] == "M") {
        int currentlength = std::stoi(tempcigarsubs[c-1]);
        if (currentlength >= 1) {
          for (int ll=0; ll < currentlength; ll++) {
            // is within range?
            int posupdate = whichCint(absposv, currentpos);

            if (posupdate != -1) {
              // update output
              switch (currentseq[base]) {
              case 'A':
                Av[posupdate] += 1;
                break;
              case 'T':
                Tv[posupdate] += 1;
                break;
              case 'G':
                Gv[posupdate] += 1;
                break;
              case 'C':
                Cv[posupdate] += 1;
                break;
              }
            }

            // update currentpos
            currentpos += 1;
            base += 1;
          }
        }
      }
    }
  }

  // return
  return(DataFrame::create(_["abspos"] = absposv, _["relpos"] = relposv, _["ref"] = refv,
                           _["A"] = Av, _["T"] = Tv, _["G"] = Gv, _["C"] = Cv, _["N"] = Nv));
}

