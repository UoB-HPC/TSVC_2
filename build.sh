SCRIPT_PATH=$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")
export SRC_DIR=$SCRIPT_PATH/src

[ -d build_dir ] || mkdir build_dir
cd build_dir

# Configure the build with cmake and build it.
cmake ..
cmake --build . --parallel 4
