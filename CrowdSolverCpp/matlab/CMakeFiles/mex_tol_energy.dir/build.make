# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.17

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Disable VCS-based implicit rules.
% : %,v


# Disable VCS-based implicit rules.
% : RCS/%


# Disable VCS-based implicit rules.
% : RCS/%,v


# Disable VCS-based implicit rules.
% : SCCS/s.%


# Disable VCS-based implicit rules.
% : s.%


.SUFFIXES: .hpux_make_needs_suffix_list


# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/local/Cellar/cmake/3.17.1/bin/cmake

# The command to remove a file.
RM = /usr/local/Cellar/cmake/3.17.1/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/vismay/recode/crowds/CrowdSolverCpp/matlab

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/vismay/recode/crowds/CrowdSolverCpp/matlab

# Include any dependencies generated for this target.
include CMakeFiles/mex_tol_energy.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/mex_tol_energy.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/mex_tol_energy.dir/flags.make

CMakeFiles/mex_tol_energy.dir/mex_tol_energy.cpp.o: CMakeFiles/mex_tol_energy.dir/flags.make
CMakeFiles/mex_tol_energy.dir/mex_tol_energy.cpp.o: mex_tol_energy.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/vismay/recode/crowds/CrowdSolverCpp/matlab/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/mex_tol_energy.dir/mex_tol_energy.cpp.o"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/mex_tol_energy.dir/mex_tol_energy.cpp.o -c /Users/vismay/recode/crowds/CrowdSolverCpp/matlab/mex_tol_energy.cpp

CMakeFiles/mex_tol_energy.dir/mex_tol_energy.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/mex_tol_energy.dir/mex_tol_energy.cpp.i"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/vismay/recode/crowds/CrowdSolverCpp/matlab/mex_tol_energy.cpp > CMakeFiles/mex_tol_energy.dir/mex_tol_energy.cpp.i

CMakeFiles/mex_tol_energy.dir/mex_tol_energy.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/mex_tol_energy.dir/mex_tol_energy.cpp.s"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/vismay/recode/crowds/CrowdSolverCpp/matlab/mex_tol_energy.cpp -o CMakeFiles/mex_tol_energy.dir/mex_tol_energy.cpp.s

# Object files for target mex_tol_energy
mex_tol_energy_OBJECTS = \
"CMakeFiles/mex_tol_energy.dir/mex_tol_energy.cpp.o"

# External object files for target mex_tol_energy
mex_tol_energy_EXTERNAL_OBJECTS =

mex_tol_energy.mexmaci64: CMakeFiles/mex_tol_energy.dir/mex_tol_energy.cpp.o
mex_tol_energy.mexmaci64: CMakeFiles/mex_tol_energy.dir/build.make
mex_tol_energy.mexmaci64: /Applications/MATLAB_R2019b.app/bin/maci64/libmex.dylib
mex_tol_energy.mexmaci64: /Applications/MATLAB_R2019b.app/bin/maci64/libmx.dylib
mex_tol_energy.mexmaci64: /Applications/MATLAB_R2019b.app/bin/maci64/libeng.dylib
mex_tol_energy.mexmaci64: /Applications/MATLAB_R2019b.app/bin/maci64/libmat.dylib
mex_tol_energy.mexmaci64: /Applications/MATLAB_R2019b.app/bin/maci64/libmex.dylib
mex_tol_energy.mexmaci64: /Applications/MATLAB_R2019b.app/bin/maci64/libmx.dylib
mex_tol_energy.mexmaci64: /Applications/MATLAB_R2019b.app/bin/maci64/libeng.dylib
mex_tol_energy.mexmaci64: /Applications/MATLAB_R2019b.app/bin/maci64/libmat.dylib
mex_tol_energy.mexmaci64: CMakeFiles/mex_tol_energy.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/vismay/recode/crowds/CrowdSolverCpp/matlab/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX shared library mex_tol_energy.mexmaci64"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/mex_tol_energy.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/mex_tol_energy.dir/build: mex_tol_energy.mexmaci64

.PHONY : CMakeFiles/mex_tol_energy.dir/build

CMakeFiles/mex_tol_energy.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/mex_tol_energy.dir/cmake_clean.cmake
.PHONY : CMakeFiles/mex_tol_energy.dir/clean

CMakeFiles/mex_tol_energy.dir/depend:
	cd /Users/vismay/recode/crowds/CrowdSolverCpp/matlab && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/vismay/recode/crowds/CrowdSolverCpp/matlab /Users/vismay/recode/crowds/CrowdSolverCpp/matlab /Users/vismay/recode/crowds/CrowdSolverCpp/matlab /Users/vismay/recode/crowds/CrowdSolverCpp/matlab /Users/vismay/recode/crowds/CrowdSolverCpp/matlab/CMakeFiles/mex_tol_energy.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/mex_tol_energy.dir/depend

