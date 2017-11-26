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

#ifndef CHELPER_H
#define CHELPER_H

#include <Rcpp.h>
using namespace Rcpp;

Rcpp::LogicalVector whichC (Rcpp::CharacterVector vectorin, Rcpp::String match);
int vecmin(Rcpp::IntegerVector x);
int vecmax(Rcpp::IntegerVector x);
bool mathintersectC(int a, int b, int x, int y);
Rcpp::IntegerVector mathrangeC(int start, int stop, int step);
Rcpp::IntegerVector readposendintermediateC(CharacterVector cigar);
std::vector< std::string > codonsC(std::string seq, int start, int stop);
std::string codonselectC(std::string seq, int def, int select, int shift);
bool insordelC(std::string cigar);
std::string cigarsubstail(std::vector < std::string> cigarsubs);
std::string generalpastelongC(std::vector <std::string> x, std::string sep);
std::vector < int > checknegC(std::vector < int > val);
int checknegCsingle(int val);
std::vector< std::string > cigarsubsC(std::string cigar);
int whichCint (Rcpp::IntegerVector vectorin, int match);
Rcpp::NumericVector cigarsubsinsertionsC(std::vector< std::string > cigarsubs, int end, int start, std::string type);
Rcpp::NumericVector cigarsubsdeletionsC(std::vector< std::string > cigarsubs, int end, int start, std::string type);
Rcpp::NumericVector cigarsubsSHC(std::vector< std::string > cigarsubs, int end, int start = 0);

#endif
