export FFTW_UTIL_INSTALL_DIR=$ANITA_UTIL_INSTALL_DIR
export FFTW_UTIL_INC_DIR=${FFTW_UTIL_INSTALL_DIR}/include
export LD_LIBRARY_PATH=${FFTW_UTIL_INSTALL_DIR}/lib:${LD_LIBRARY_PATH}
export PATH=${FFTW_UTIL_INSTALL_DIR}/bin:${PATH}
