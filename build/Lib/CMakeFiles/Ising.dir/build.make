# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.16

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


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
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/wenqingxie/Ising-Wolff/Code

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/wenqingxie/Ising-Wolff/build

# Include any dependencies generated for this target.
include Lib/CMakeFiles/Ising.dir/depend.make

# Include the progress variables for this target.
include Lib/CMakeFiles/Ising.dir/progress.make

# Include the compile flags for this target's objects.
include Lib/CMakeFiles/Ising.dir/flags.make

Lib/CMakeFiles/Ising.dir/gridState.cpp.o: Lib/CMakeFiles/Ising.dir/flags.make
Lib/CMakeFiles/Ising.dir/gridState.cpp.o: /home/wenqingxie/Ising-Wolff/Code/Lib/gridState.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/wenqingxie/Ising-Wolff/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object Lib/CMakeFiles/Ising.dir/gridState.cpp.o"
	cd /home/wenqingxie/Ising-Wolff/build/Lib && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/Ising.dir/gridState.cpp.o -c /home/wenqingxie/Ising-Wolff/Code/Lib/gridState.cpp

Lib/CMakeFiles/Ising.dir/gridState.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/Ising.dir/gridState.cpp.i"
	cd /home/wenqingxie/Ising-Wolff/build/Lib && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/wenqingxie/Ising-Wolff/Code/Lib/gridState.cpp > CMakeFiles/Ising.dir/gridState.cpp.i

Lib/CMakeFiles/Ising.dir/gridState.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/Ising.dir/gridState.cpp.s"
	cd /home/wenqingxie/Ising-Wolff/build/Lib && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/wenqingxie/Ising-Wolff/Code/Lib/gridState.cpp -o CMakeFiles/Ising.dir/gridState.cpp.s

Lib/CMakeFiles/Ising.dir/isingState.cpp.o: Lib/CMakeFiles/Ising.dir/flags.make
Lib/CMakeFiles/Ising.dir/isingState.cpp.o: /home/wenqingxie/Ising-Wolff/Code/Lib/isingState.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/wenqingxie/Ising-Wolff/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object Lib/CMakeFiles/Ising.dir/isingState.cpp.o"
	cd /home/wenqingxie/Ising-Wolff/build/Lib && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/Ising.dir/isingState.cpp.o -c /home/wenqingxie/Ising-Wolff/Code/Lib/isingState.cpp

Lib/CMakeFiles/Ising.dir/isingState.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/Ising.dir/isingState.cpp.i"
	cd /home/wenqingxie/Ising-Wolff/build/Lib && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/wenqingxie/Ising-Wolff/Code/Lib/isingState.cpp > CMakeFiles/Ising.dir/isingState.cpp.i

Lib/CMakeFiles/Ising.dir/isingState.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/Ising.dir/isingState.cpp.s"
	cd /home/wenqingxie/Ising-Wolff/build/Lib && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/wenqingxie/Ising-Wolff/Code/Lib/isingState.cpp -o CMakeFiles/Ising.dir/isingState.cpp.s

Lib/CMakeFiles/Ising.dir/monteCarlo.cpp.o: Lib/CMakeFiles/Ising.dir/flags.make
Lib/CMakeFiles/Ising.dir/monteCarlo.cpp.o: /home/wenqingxie/Ising-Wolff/Code/Lib/monteCarlo.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/wenqingxie/Ising-Wolff/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object Lib/CMakeFiles/Ising.dir/monteCarlo.cpp.o"
	cd /home/wenqingxie/Ising-Wolff/build/Lib && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/Ising.dir/monteCarlo.cpp.o -c /home/wenqingxie/Ising-Wolff/Code/Lib/monteCarlo.cpp

Lib/CMakeFiles/Ising.dir/monteCarlo.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/Ising.dir/monteCarlo.cpp.i"
	cd /home/wenqingxie/Ising-Wolff/build/Lib && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/wenqingxie/Ising-Wolff/Code/Lib/monteCarlo.cpp > CMakeFiles/Ising.dir/monteCarlo.cpp.i

Lib/CMakeFiles/Ising.dir/monteCarlo.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/Ising.dir/monteCarlo.cpp.s"
	cd /home/wenqingxie/Ising-Wolff/build/Lib && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/wenqingxie/Ising-Wolff/Code/Lib/monteCarlo.cpp -o CMakeFiles/Ising.dir/monteCarlo.cpp.s

# Object files for target Ising
Ising_OBJECTS = \
"CMakeFiles/Ising.dir/gridState.cpp.o" \
"CMakeFiles/Ising.dir/isingState.cpp.o" \
"CMakeFiles/Ising.dir/monteCarlo.cpp.o"

# External object files for target Ising
Ising_EXTERNAL_OBJECTS =

Lib/libIsing.a: Lib/CMakeFiles/Ising.dir/gridState.cpp.o
Lib/libIsing.a: Lib/CMakeFiles/Ising.dir/isingState.cpp.o
Lib/libIsing.a: Lib/CMakeFiles/Ising.dir/monteCarlo.cpp.o
Lib/libIsing.a: Lib/CMakeFiles/Ising.dir/build.make
Lib/libIsing.a: Lib/CMakeFiles/Ising.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/wenqingxie/Ising-Wolff/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Linking CXX static library libIsing.a"
	cd /home/wenqingxie/Ising-Wolff/build/Lib && $(CMAKE_COMMAND) -P CMakeFiles/Ising.dir/cmake_clean_target.cmake
	cd /home/wenqingxie/Ising-Wolff/build/Lib && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/Ising.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
Lib/CMakeFiles/Ising.dir/build: Lib/libIsing.a

.PHONY : Lib/CMakeFiles/Ising.dir/build

Lib/CMakeFiles/Ising.dir/clean:
	cd /home/wenqingxie/Ising-Wolff/build/Lib && $(CMAKE_COMMAND) -P CMakeFiles/Ising.dir/cmake_clean.cmake
.PHONY : Lib/CMakeFiles/Ising.dir/clean

Lib/CMakeFiles/Ising.dir/depend:
	cd /home/wenqingxie/Ising-Wolff/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/wenqingxie/Ising-Wolff/Code /home/wenqingxie/Ising-Wolff/Code/Lib /home/wenqingxie/Ising-Wolff/build /home/wenqingxie/Ising-Wolff/build/Lib /home/wenqingxie/Ising-Wolff/build/Lib/CMakeFiles/Ising.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : Lib/CMakeFiles/Ising.dir/depend

