# module load openmpi/4.1.3

. $HOME/zen2_dec/dec_spack/share/spack/setup-env.sh

spack load libconfig libxml2 jansson@2.9 # texinfo@5.0 jansson

HDF5_DIR="/qfs/people/tang584/install/hdf5"
HDF5_BIN="$HDF5_DIR/bin"
HDF5_LIB="$HDF5_DIR/lib"
HDF5_INCLUDE="$HDF5_DIR/include"
PATH="${HDF5_BIN}:${PATH}"
LD_LIBRARY_PATH="${HDF5_LIB}:${LD_LIBRARY_PATH}"
C_INCLUDE_PATH="${HDF5_INCLUDE}:${C_INCLUDE_PATH}"
HDF5_PKG_CONFIG_PATH="$HDF5_DIR/lib/pkgconfig"
PKG_CONFIG_PATH="${HDF5_PKG_CONFIG_PATH}:${PKG_CONFIG_PATH}"

NETCDFC_DIR="/qfs/people/tang584/install/netcdf-c"
NETCDFC_BIN="$NETCDFC_DIR/bin"
NETCDFC_LIB="$NETCDFC_DIR/lib"
NETCDFC_INCLUDE="$NETCDFC_DIR/include"
PATH="${NETCDFC_BIN}:${PATH}"
LD_LIBRARY_PATH="${NETCDFC_LIB}:${LD_LIBRARY_PATH}"
C_INCLUDE_PATH="${NETCDFC_INCLUDE}:${C_INCLUDE_PATH}"
NETCDFC_PKG_CONFIG_PATH="$NETCDFC_DIR/lib/pkgconfig"
PKG_CONFIG_PATH="${NETCDFC_PKG_CONFIG_PATH}:${PKG_CONFIG_PATH}"

LIBCONFIG_DIR="/qfs/people/tang584/zen2_dec/dec_spack/opt/spack/linux-centos7-zen2/gcc-9.1.0/libconfig-1.7.2-htzxtlcxpvlcyb3cti7sr7ttnzuiwdjb"
LD_LIBRARY_PATH="${LIBCONFIG_DIR}:${LD_LIBRARY_PATH}"

JANSSON_DIR=/qfs/people/tang584/zen2_dec/dec_spack/opt/spack/linux-centos7-zen2/gcc-9.1.0/jansson-2.9-6mxyawnx3ftbm325dyjlumi3m2bw2c3a
LD_LIBRARY_PATH="${JANSSON_DIR}:${LD_LIBRARY_PATH}"