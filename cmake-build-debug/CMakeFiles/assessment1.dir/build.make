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
CMAKE_COMMAND = "/Users/zhiyongliu/Library/Application Support/JetBrains/Toolbox/apps/CLion/ch-0/203.6682.181/CLion.app/Contents/bin/cmake/mac/bin/cmake"

# The command to remove a file.
RM = "/Users/zhiyongliu/Library/Application Support/JetBrains/Toolbox/apps/CLion/ch-0/203.6682.181/CLion.app/Contents/bin/cmake/mac/bin/cmake" -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/zhiyongliu/Documents/Imperial_Github/ACSE-5/assessment1

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/zhiyongliu/Documents/Imperial_Github/ACSE-5/assessment1/cmake-build-debug

# Include any dependencies generated for this target.
include CMakeFiles/assessment1.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/assessment1.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/assessment1.dir/flags.make

CMakeFiles/assessment1.dir/main.cpp.o: CMakeFiles/assessment1.dir/flags.make
CMakeFiles/assessment1.dir/main.cpp.o: ../main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/zhiyongliu/Documents/Imperial_Github/ACSE-5/assessment1/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/assessment1.dir/main.cpp.o"
	/Library/Developer/CommandLineTools/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/assessment1.dir/main.cpp.o -c /Users/zhiyongliu/Documents/Imperial_Github/ACSE-5/assessment1/main.cpp

CMakeFiles/assessment1.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/assessment1.dir/main.cpp.i"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/zhiyongliu/Documents/Imperial_Github/ACSE-5/assessment1/main.cpp > CMakeFiles/assessment1.dir/main.cpp.i

CMakeFiles/assessment1.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/assessment1.dir/main.cpp.s"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/zhiyongliu/Documents/Imperial_Github/ACSE-5/assessment1/main.cpp -o CMakeFiles/assessment1.dir/main.cpp.s

CMakeFiles/assessment1.dir/Matrix.cpp.o: CMakeFiles/assessment1.dir/flags.make
CMakeFiles/assessment1.dir/Matrix.cpp.o: ../Matrix.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/zhiyongliu/Documents/Imperial_Github/ACSE-5/assessment1/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/assessment1.dir/Matrix.cpp.o"
	/Library/Developer/CommandLineTools/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/assessment1.dir/Matrix.cpp.o -c /Users/zhiyongliu/Documents/Imperial_Github/ACSE-5/assessment1/Matrix.cpp

CMakeFiles/assessment1.dir/Matrix.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/assessment1.dir/Matrix.cpp.i"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/zhiyongliu/Documents/Imperial_Github/ACSE-5/assessment1/Matrix.cpp > CMakeFiles/assessment1.dir/Matrix.cpp.i

CMakeFiles/assessment1.dir/Matrix.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/assessment1.dir/Matrix.cpp.s"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/zhiyongliu/Documents/Imperial_Github/ACSE-5/assessment1/Matrix.cpp -o CMakeFiles/assessment1.dir/Matrix.cpp.s

CMakeFiles/assessment1.dir/Solver.cpp.o: CMakeFiles/assessment1.dir/flags.make
CMakeFiles/assessment1.dir/Solver.cpp.o: ../Solver.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/zhiyongliu/Documents/Imperial_Github/ACSE-5/assessment1/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object CMakeFiles/assessment1.dir/Solver.cpp.o"
	/Library/Developer/CommandLineTools/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/assessment1.dir/Solver.cpp.o -c /Users/zhiyongliu/Documents/Imperial_Github/ACSE-5/assessment1/Solver.cpp

CMakeFiles/assessment1.dir/Solver.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/assessment1.dir/Solver.cpp.i"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/zhiyongliu/Documents/Imperial_Github/ACSE-5/assessment1/Solver.cpp > CMakeFiles/assessment1.dir/Solver.cpp.i

CMakeFiles/assessment1.dir/Solver.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/assessment1.dir/Solver.cpp.s"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/zhiyongliu/Documents/Imperial_Github/ACSE-5/assessment1/Solver.cpp -o CMakeFiles/assessment1.dir/Solver.cpp.s

CMakeFiles/assessment1.dir/CSRMatrix.cpp.o: CMakeFiles/assessment1.dir/flags.make
CMakeFiles/assessment1.dir/CSRMatrix.cpp.o: ../CSRMatrix.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/zhiyongliu/Documents/Imperial_Github/ACSE-5/assessment1/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object CMakeFiles/assessment1.dir/CSRMatrix.cpp.o"
	/Library/Developer/CommandLineTools/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/assessment1.dir/CSRMatrix.cpp.o -c /Users/zhiyongliu/Documents/Imperial_Github/ACSE-5/assessment1/CSRMatrix.cpp

CMakeFiles/assessment1.dir/CSRMatrix.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/assessment1.dir/CSRMatrix.cpp.i"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/zhiyongliu/Documents/Imperial_Github/ACSE-5/assessment1/CSRMatrix.cpp > CMakeFiles/assessment1.dir/CSRMatrix.cpp.i

CMakeFiles/assessment1.dir/CSRMatrix.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/assessment1.dir/CSRMatrix.cpp.s"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/zhiyongliu/Documents/Imperial_Github/ACSE-5/assessment1/CSRMatrix.cpp -o CMakeFiles/assessment1.dir/CSRMatrix.cpp.s

CMakeFiles/assessment1.dir/Test.cpp.o: CMakeFiles/assessment1.dir/flags.make
CMakeFiles/assessment1.dir/Test.cpp.o: ../Test.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/zhiyongliu/Documents/Imperial_Github/ACSE-5/assessment1/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building CXX object CMakeFiles/assessment1.dir/Test.cpp.o"
	/Library/Developer/CommandLineTools/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/assessment1.dir/Test.cpp.o -c /Users/zhiyongliu/Documents/Imperial_Github/ACSE-5/assessment1/Test.cpp

CMakeFiles/assessment1.dir/Test.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/assessment1.dir/Test.cpp.i"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/zhiyongliu/Documents/Imperial_Github/ACSE-5/assessment1/Test.cpp > CMakeFiles/assessment1.dir/Test.cpp.i

CMakeFiles/assessment1.dir/Test.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/assessment1.dir/Test.cpp.s"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/zhiyongliu/Documents/Imperial_Github/ACSE-5/assessment1/Test.cpp -o CMakeFiles/assessment1.dir/Test.cpp.s

# Object files for target assessment1
assessment1_OBJECTS = \
"CMakeFiles/assessment1.dir/main.cpp.o" \
"CMakeFiles/assessment1.dir/Matrix.cpp.o" \
"CMakeFiles/assessment1.dir/Solver.cpp.o" \
"CMakeFiles/assessment1.dir/CSRMatrix.cpp.o" \
"CMakeFiles/assessment1.dir/Test.cpp.o"

# External object files for target assessment1
assessment1_EXTERNAL_OBJECTS =

assessment1: CMakeFiles/assessment1.dir/main.cpp.o
assessment1: CMakeFiles/assessment1.dir/Matrix.cpp.o
assessment1: CMakeFiles/assessment1.dir/Solver.cpp.o
assessment1: CMakeFiles/assessment1.dir/CSRMatrix.cpp.o
assessment1: CMakeFiles/assessment1.dir/Test.cpp.o
assessment1: CMakeFiles/assessment1.dir/build.make
assessment1: CMakeFiles/assessment1.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/zhiyongliu/Documents/Imperial_Github/ACSE-5/assessment1/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Linking CXX executable assessment1"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/assessment1.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/assessment1.dir/build: assessment1

.PHONY : CMakeFiles/assessment1.dir/build

CMakeFiles/assessment1.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/assessment1.dir/cmake_clean.cmake
.PHONY : CMakeFiles/assessment1.dir/clean

CMakeFiles/assessment1.dir/depend:
	cd /Users/zhiyongliu/Documents/Imperial_Github/ACSE-5/assessment1/cmake-build-debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/zhiyongliu/Documents/Imperial_Github/ACSE-5/assessment1 /Users/zhiyongliu/Documents/Imperial_Github/ACSE-5/assessment1 /Users/zhiyongliu/Documents/Imperial_Github/ACSE-5/assessment1/cmake-build-debug /Users/zhiyongliu/Documents/Imperial_Github/ACSE-5/assessment1/cmake-build-debug /Users/zhiyongliu/Documents/Imperial_Github/ACSE-5/assessment1/cmake-build-debug/CMakeFiles/assessment1.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/assessment1.dir/depend

