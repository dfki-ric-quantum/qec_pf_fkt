// COPYRIGHT NOTICE
//  This file is part of isingZ, a program to compute partition
//      functions and free energies for domain walls in 2D Ising models.
//  Copyright 2012 by Creighton K. Thomas (creightonthomas@gmail.com),
//                    A. Alan Middleton (aam@syr.edu)
//  The development of this software was supported in part by the National Science
//  Foundation under grant DMR-1006731.
//
//  The program isingZ is free software; you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation; either version 3 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with this program; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
//  See http://www.gnu.org/licenses/gpl.txt for more details.
//
//

// main.cc
// See usage comment in main routine.
// Computes and reports partition functions and free energies for combinations
// of boundary conditions.

#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
#include <iomanip>
#include "Sample.h"
#include "FINDmatrix.h"
#include <cstdlib>
#include "exp_log.h"

void createDirectory(const std::string &path) {
    std::string command = "mkdir -p " + path;
    int status = system(command.c_str());
    if (status != 0) {
        std::cerr << "Error creating directory: " << path << std::endl;
    }
}

void findPartition(Sample &S, const std::string &outputFile, const int precision) {
  FINDmatrix X(&S);

  FINDmatrix Ypls1(X);
  FINDmatrix Yneg1(X);

  Ypls1.wrapHorz(1);
  Yneg1.wrapHorz(-1);

  FINDmatrix Ypls2(Ypls1);
  FINDmatrix Yneg2(Yneg1);

  dataType y1 = Ypls1.Zvert(1);
  dataType y2 = Yneg1.Zvert(1);
  dataType y3 = Ypls2.Zvert(-1);
  dataType y4 = Yneg2.Zvert(-1);

  dataType prefactor = S.get_Z_prefactor();

  dataType ZPP = abs(prefactor*0.5*( y1+y2+y3+y4));
  dataType ZPA = abs(prefactor*0.5*(-y1-y2+y3+y4));
  dataType ZAP = abs(prefactor*0.5*(-y1+y2-y3+y4));
  dataType ZAA = abs(prefactor*0.5*(-y1+y2+y3-y4));

  std::ofstream outFile(outputFile.c_str());

  // Set precision based on the input precision parameter
  outFile.precision(int(precision * 0.301)); // Convert bits to decimal digits
  outFile << std::scientific;               // Use scientific notation

  outFile << ZPP << "\t"
          << ZPA << "\t"
          << ZAP << "\t"
          << ZAA << "\t";
  outFile.close();
}

int main(int argc, char* argv[])
{
  if (argc < 8 || argc > 9)
  {
    std::cout << "FIND2DIsing: computes partition function of 2D Ising model on a square lattice\n";
    std::cout << "usage: " << argv[0] << " bitsOfPrecision Lx Ly seed probability temperature directory [std dev] \n";
    return 1;
  }

  int prec = atoi(argv[1]);
  mpf_set_default_prec(prec);
  std::cout.precision(int(prec*0.301));

  int x   = atoi(argv[2]);
  int y   = atoi(argv[3]);
  int seed = atoi(argv[4]);
  double prob = atof(argv[5]);
  double T_frac = atof(argv[6]);

  bool useGaussian = (argc == 9);
  double stddev = 0.0;
  dataType T_nish = 1.0; // Nishimori temperature, default is 1.0
  if (prob!=0){
    T_nish = 2/std::log((1-prob)/prob); // if use Gaussian noise is false, couplings are normalized to one and Nishimori temperature is p dependent
  }
  if (useGaussian) {
    T_nish = 1.0; // For Gaussian noise, couplings are not normalized and p dependent and Nishimori temperature normalized to one.
    stddev = std::atof(argv[8]);
    if (stddev <= 0) {
      std::cerr << "Error: Std dev must be positive.\n";
      return 1;
    }
  }
  dataType T = T_frac*T_nish;

  std::string directory = argv[7];

  std::string input =   directory + "/interactionsGaussian/" +
                            std::to_string(prob) + "/" +
                            std::to_string(x) + "/" +
                            std::to_string(y) + "/" +
                            std::to_string(stddev) + "/" +
                            std::to_string(seed) + "/interaction_lattice.txt";

  Sample S(input, T);

  std::string outputDir =   directory + "/resultsGaussian/" +
                            std::to_string(prob) + "/" +
                            std::to_string(stddev) + "/" +
                            std::to_string(x) + "/" +
                            std::to_string(y) + "/" +
                            std::to_string(T_frac) + "/" +
                            std::to_string(prec) + "/" +
                            std::to_string(seed);

  createDirectory(outputDir);

  std::string outputFile = outputDir + "/Z.txt";

  findPartition(S, outputFile, prec);
  std::cout << "Z results written to: " << outputDir << std::endl;
  return 0;
}
