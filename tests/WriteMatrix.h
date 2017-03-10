#ifndef INCLUDE_TESTING_WRITE_EIGEN_MATRIX_H
#define INCLUDE_TESTING_WRITE_EIGEN_MATRIX_H

#include <Eigen/Core>
#include <fstream>
#include <string>

using namespace std::string_literals;

void writeMatrix(
  const std::string& name,
  const Eigen::MatrixXd& matrix
) {
  std::ofstream outStream(name + ".csv"s);

  for(unsigned i = 0; i < matrix.rows(); i++) {
    for(unsigned j = 0; j < matrix.cols(); j++) {
      outStream << matrix(i, j);
      if(j != matrix.cols() - 1) outStream << ",";
    }
    
    outStream << std::endl;
  }

  outStream.close();
}

#endif
