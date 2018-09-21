project(spectra_lib)

set(spectra_VERSION 0.6.2)

set(HEADERS_All
  ${CMAKE_CURRENT_SOURCE_DIR}/Spectra/Util/CompInfo.h
  ${CMAKE_CURRENT_SOURCE_DIR}/Spectra/Util/TypeTraits.h
  ${CMAKE_CURRENT_SOURCE_DIR}/Spectra/Util/SimpleRandom.h
  ${CMAKE_CURRENT_SOURCE_DIR}/Spectra/Util/SelectionRule.h
  ${CMAKE_CURRENT_SOURCE_DIR}/Spectra/Util/GEigsMode.h
  ${CMAKE_CURRENT_SOURCE_DIR}/Spectra/contrib/PartialSVDSolver.h
  ${CMAKE_CURRENT_SOURCE_DIR}/Spectra/contrib/LOBPCGSolver.h
  ${CMAKE_CURRENT_SOURCE_DIR}/Spectra/MatOp/DenseSymShiftSolve.h
  ${CMAKE_CURRENT_SOURCE_DIR}/Spectra/MatOp/SparseSymShiftSolve.h
  ${CMAKE_CURRENT_SOURCE_DIR}/Spectra/MatOp/SparseGenMatProd.h
  ${CMAKE_CURRENT_SOURCE_DIR}/Spectra/MatOp/internal/SymGEigsCholeskyOp.h
  ${CMAKE_CURRENT_SOURCE_DIR}/Spectra/MatOp/internal/SymGEigsRegInvOp.h
  ${CMAKE_CURRENT_SOURCE_DIR}/Spectra/MatOp/DenseSymMatProd.h
  ${CMAKE_CURRENT_SOURCE_DIR}/Spectra/MatOp/SparseSymMatProd.h
  ${CMAKE_CURRENT_SOURCE_DIR}/Spectra/MatOp/SparseRegularInverse.h
  ${CMAKE_CURRENT_SOURCE_DIR}/Spectra/MatOp/DenseGenComplexShiftSolve.h
  ${CMAKE_CURRENT_SOURCE_DIR}/Spectra/MatOp/DenseCholesky.h
  ${CMAKE_CURRENT_SOURCE_DIR}/Spectra/MatOp/DenseGenMatProd.h
  ${CMAKE_CURRENT_SOURCE_DIR}/Spectra/MatOp/SparseCholesky.h
  ${CMAKE_CURRENT_SOURCE_DIR}/Spectra/MatOp/DenseGenRealShiftSolve.h
  ${CMAKE_CURRENT_SOURCE_DIR}/Spectra/MatOp/SparseGenRealShiftSolve.h
  ${CMAKE_CURRENT_SOURCE_DIR}/Spectra/SymEigsSolver.h
  ${CMAKE_CURRENT_SOURCE_DIR}/Spectra/GenEigsComplexShiftSolver.h
  ${CMAKE_CURRENT_SOURCE_DIR}/Spectra/LinAlg/UpperHessenbergEigen.h
  ${CMAKE_CURRENT_SOURCE_DIR}/Spectra/LinAlg/DoubleShiftQR.h
  ${CMAKE_CURRENT_SOURCE_DIR}/Spectra/LinAlg/UpperHessenbergQR.h
  ${CMAKE_CURRENT_SOURCE_DIR}/Spectra/LinAlg/TridiagEigen.h
  ${CMAKE_CURRENT_SOURCE_DIR}/Spectra/GenEigsSolver.h
  ${CMAKE_CURRENT_SOURCE_DIR}/Spectra/GenEigsRealShiftSolver.h
  ${CMAKE_CURRENT_SOURCE_DIR}/Spectra/SymEigsShiftSolver.h
  ${CMAKE_CURRENT_SOURCE_DIR}/Spectra/SymGEigsSolver.h
)

add_library(spectra INTERFACE)
target_sources(spectra INTERFACE $<BUILD_INTERFACE:${HEADERS_All}>)
target_include_directories(
  spectra INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:$<INSTALL_PREFIX>/include>
)

source_group("Headers" FILES ${HEADERS_All})
