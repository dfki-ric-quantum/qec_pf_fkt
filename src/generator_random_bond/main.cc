#include <cmath>
#include <fstream>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <iostream>

void createDirectory(const std::string &path) {
  std::string command = "mkdir -p " + path;
  int status = system(command.c_str());
  if (status != 0) {
      std::cerr << "Error creating directory: " << path << std::endl;
  }
}

int main(int argc, char *argv[]) {
  if (argc < 6 || argc > 7) {
    std::cout << "Usage: " << argv[0]
              << " Lx Ly seed probability directory [std deviation]\n";
    return 1;
  }

  int Lx = atoi(argv[1]);
  int Ly = atoi(argv[2]);
  int seed = atoi(argv[3]);
  double prob = atof(argv[4]);
  std::string directory = argv[5];
  bool useGaussian = (argc == 7);
  double stddev = 0.0;
  if (useGaussian) {
    stddev = std::atof(argv[6]);
    if (stddev <= 0) {
      std::cerr << "Error: Std dev must be positive.\n";
      return 1;
    }
  }

  gsl_rng *rng;
  rng = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set(rng, seed);

  std::string outputDir = directory + "/interactionsGaussian/" + std::to_string(prob) +
                          "/" + std::to_string(Lx) + "/" + std::to_string(Ly) +
                          "/" + std::to_string(stddev) + "/" + std::to_string(seed);

  createDirectory(outputDir);

  std::ofstream outFile(outputDir + "/interaction_lattice.txt");

  outFile << Lx << " " << Ly << "\n";
  for (int j = 0; j < Ly; j++) {
    for (int i = 0; i < Lx; i++) {
      double valueE, valueS;

      if (useGaussian) {
        double probE = gsl_ran_gaussian(rng, stddev) + prob;
        double probS = gsl_ran_gaussian(rng, stddev) + prob;

        probE = std::min(std::max(probE, 1e-4), 0.5 - 1e-10); // ensures physicality of error probabilities
        probS = std::min(std::max(probS, 1e-4), 0.5 - 1e-10);

        int flip_interactionE = (gsl_rng_uniform(rng) < probE) ? -1 : 1;
        int flip_interactionS = (gsl_rng_uniform(rng) < probS) ? -1 : 1;

        valueE = flip_interactionE * 0.5 * std::log((1.0 - probE) / probE);
        valueS = flip_interactionS * 0.5 * std::log((1.0 - probS) / probS);
      } else {
        valueE = (gsl_rng_uniform(rng) < prob) ? -1 : 1;
        valueS = (gsl_rng_uniform(rng) < prob) ? -1 : 1;
      }

      outFile << i << "\t" << j << "\tE\t" << valueE << "\n";
      outFile << i << "\t" << j << "\tS\t" << valueS << "\n";
    }
  }

  outFile.close();
  gsl_rng_free(rng);
  return 0;
}
