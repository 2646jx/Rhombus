USE_HEXL=OFF    # Use AVX512 to accelerate HE operations
USE_ZSTD=OFF    # Use zstd to compress the size of HE ciphertexts

WORK_DIR=`pwd`
BUILD_DIR=$WORK_DIR/build
DEPS_DIR=$WORK_DIR/deps

# SEAL
git clone https://github.com/microsoft/SEAL.git $DEPS_DIR/SEAL

# build SEAL
target=SEAL
cd $DEPS_DIR/$target
git checkout 119dc32    # v4.1.2
mkdir -p $BUILD_DIR/deps/$target
cd $BUILD_DIR/deps/$target
cmake $DEPS_DIR/$target -DCMAKE_INSTALL_PREFIX=$BUILD_DIR\
                        -DCMAKE_PREFIX_PATH=$BUILD_DIR\
                        -DSEAL_USE_MSGSL=OFF\
                        -DSEAL_USE_ZLIB=OFF\
                        -DSEAL_USE_ZSTD=$USE_ZSTD\
                        -DSEAL_USE_INTEL_HEXL=$USE_HEXL\
                        -DCMAKE_BUILD_TYPE=Release\
                        -DSEAL_THROW_ON_TRANSPARENT_CIPHERTEXT=OFF
make install -j4
