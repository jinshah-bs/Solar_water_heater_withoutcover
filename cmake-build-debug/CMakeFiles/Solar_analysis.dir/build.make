# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.20

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

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /var/lib/snapd/snap/clion/169/bin/cmake/linux/bin/cmake

# The command to remove a file.
RM = /var/lib/snapd/snap/clion/169/bin/cmake/linux/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/thinkcomputational/CLionProjects/Solar_water_heater_withoutcover

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/thinkcomputational/CLionProjects/Solar_water_heater_withoutcover/cmake-build-debug

# Include any dependencies generated for this target.
include CMakeFiles/Solar_analysis.dir/depend.make
# Include the progress variables for this target.
include CMakeFiles/Solar_analysis.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/Solar_analysis.dir/flags.make

CMakeFiles/Solar_analysis.dir/main.cpp.o: CMakeFiles/Solar_analysis.dir/flags.make
CMakeFiles/Solar_analysis.dir/main.cpp.o: ../main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/thinkcomputational/CLionProjects/Solar_water_heater_withoutcover/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/Solar_analysis.dir/main.cpp.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/Solar_analysis.dir/main.cpp.o -c /home/thinkcomputational/CLionProjects/Solar_water_heater_withoutcover/main.cpp

CMakeFiles/Solar_analysis.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/Solar_analysis.dir/main.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/thinkcomputational/CLionProjects/Solar_water_heater_withoutcover/main.cpp > CMakeFiles/Solar_analysis.dir/main.cpp.i

CMakeFiles/Solar_analysis.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/Solar_analysis.dir/main.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/thinkcomputational/CLionProjects/Solar_water_heater_withoutcover/main.cpp -o CMakeFiles/Solar_analysis.dir/main.cpp.s

CMakeFiles/Solar_analysis.dir/read_data.cpp.o: CMakeFiles/Solar_analysis.dir/flags.make
CMakeFiles/Solar_analysis.dir/read_data.cpp.o: ../read_data.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/thinkcomputational/CLionProjects/Solar_water_heater_withoutcover/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/Solar_analysis.dir/read_data.cpp.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/Solar_analysis.dir/read_data.cpp.o -c /home/thinkcomputational/CLionProjects/Solar_water_heater_withoutcover/read_data.cpp

CMakeFiles/Solar_analysis.dir/read_data.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/Solar_analysis.dir/read_data.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/thinkcomputational/CLionProjects/Solar_water_heater_withoutcover/read_data.cpp > CMakeFiles/Solar_analysis.dir/read_data.cpp.i

CMakeFiles/Solar_analysis.dir/read_data.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/Solar_analysis.dir/read_data.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/thinkcomputational/CLionProjects/Solar_water_heater_withoutcover/read_data.cpp -o CMakeFiles/Solar_analysis.dir/read_data.cpp.s

CMakeFiles/Solar_analysis.dir/mat_prop.cpp.o: CMakeFiles/Solar_analysis.dir/flags.make
CMakeFiles/Solar_analysis.dir/mat_prop.cpp.o: ../mat_prop.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/thinkcomputational/CLionProjects/Solar_water_heater_withoutcover/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object CMakeFiles/Solar_analysis.dir/mat_prop.cpp.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/Solar_analysis.dir/mat_prop.cpp.o -c /home/thinkcomputational/CLionProjects/Solar_water_heater_withoutcover/mat_prop.cpp

CMakeFiles/Solar_analysis.dir/mat_prop.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/Solar_analysis.dir/mat_prop.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/thinkcomputational/CLionProjects/Solar_water_heater_withoutcover/mat_prop.cpp > CMakeFiles/Solar_analysis.dir/mat_prop.cpp.i

CMakeFiles/Solar_analysis.dir/mat_prop.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/Solar_analysis.dir/mat_prop.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/thinkcomputational/CLionProjects/Solar_water_heater_withoutcover/mat_prop.cpp -o CMakeFiles/Solar_analysis.dir/mat_prop.cpp.s

CMakeFiles/Solar_analysis.dir/TDMA.cpp.o: CMakeFiles/Solar_analysis.dir/flags.make
CMakeFiles/Solar_analysis.dir/TDMA.cpp.o: ../TDMA.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/thinkcomputational/CLionProjects/Solar_water_heater_withoutcover/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object CMakeFiles/Solar_analysis.dir/TDMA.cpp.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/Solar_analysis.dir/TDMA.cpp.o -c /home/thinkcomputational/CLionProjects/Solar_water_heater_withoutcover/TDMA.cpp

CMakeFiles/Solar_analysis.dir/TDMA.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/Solar_analysis.dir/TDMA.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/thinkcomputational/CLionProjects/Solar_water_heater_withoutcover/TDMA.cpp > CMakeFiles/Solar_analysis.dir/TDMA.cpp.i

CMakeFiles/Solar_analysis.dir/TDMA.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/Solar_analysis.dir/TDMA.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/thinkcomputational/CLionProjects/Solar_water_heater_withoutcover/TDMA.cpp -o CMakeFiles/Solar_analysis.dir/TDMA.cpp.s

CMakeFiles/Solar_analysis.dir/HTC.cpp.o: CMakeFiles/Solar_analysis.dir/flags.make
CMakeFiles/Solar_analysis.dir/HTC.cpp.o: ../HTC.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/thinkcomputational/CLionProjects/Solar_water_heater_withoutcover/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building CXX object CMakeFiles/Solar_analysis.dir/HTC.cpp.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/Solar_analysis.dir/HTC.cpp.o -c /home/thinkcomputational/CLionProjects/Solar_water_heater_withoutcover/HTC.cpp

CMakeFiles/Solar_analysis.dir/HTC.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/Solar_analysis.dir/HTC.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/thinkcomputational/CLionProjects/Solar_water_heater_withoutcover/HTC.cpp > CMakeFiles/Solar_analysis.dir/HTC.cpp.i

CMakeFiles/Solar_analysis.dir/HTC.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/Solar_analysis.dir/HTC.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/thinkcomputational/CLionProjects/Solar_water_heater_withoutcover/HTC.cpp -o CMakeFiles/Solar_analysis.dir/HTC.cpp.s

CMakeFiles/Solar_analysis.dir/Write_data.cpp.o: CMakeFiles/Solar_analysis.dir/flags.make
CMakeFiles/Solar_analysis.dir/Write_data.cpp.o: ../Write_data.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/thinkcomputational/CLionProjects/Solar_water_heater_withoutcover/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Building CXX object CMakeFiles/Solar_analysis.dir/Write_data.cpp.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/Solar_analysis.dir/Write_data.cpp.o -c /home/thinkcomputational/CLionProjects/Solar_water_heater_withoutcover/Write_data.cpp

CMakeFiles/Solar_analysis.dir/Write_data.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/Solar_analysis.dir/Write_data.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/thinkcomputational/CLionProjects/Solar_water_heater_withoutcover/Write_data.cpp > CMakeFiles/Solar_analysis.dir/Write_data.cpp.i

CMakeFiles/Solar_analysis.dir/Write_data.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/Solar_analysis.dir/Write_data.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/thinkcomputational/CLionProjects/Solar_water_heater_withoutcover/Write_data.cpp -o CMakeFiles/Solar_analysis.dir/Write_data.cpp.s

# Object files for target Solar_analysis
Solar_analysis_OBJECTS = \
"CMakeFiles/Solar_analysis.dir/main.cpp.o" \
"CMakeFiles/Solar_analysis.dir/read_data.cpp.o" \
"CMakeFiles/Solar_analysis.dir/mat_prop.cpp.o" \
"CMakeFiles/Solar_analysis.dir/TDMA.cpp.o" \
"CMakeFiles/Solar_analysis.dir/HTC.cpp.o" \
"CMakeFiles/Solar_analysis.dir/Write_data.cpp.o"

# External object files for target Solar_analysis
Solar_analysis_EXTERNAL_OBJECTS =

Solar_analysis: CMakeFiles/Solar_analysis.dir/main.cpp.o
Solar_analysis: CMakeFiles/Solar_analysis.dir/read_data.cpp.o
Solar_analysis: CMakeFiles/Solar_analysis.dir/mat_prop.cpp.o
Solar_analysis: CMakeFiles/Solar_analysis.dir/TDMA.cpp.o
Solar_analysis: CMakeFiles/Solar_analysis.dir/HTC.cpp.o
Solar_analysis: CMakeFiles/Solar_analysis.dir/Write_data.cpp.o
Solar_analysis: CMakeFiles/Solar_analysis.dir/build.make
Solar_analysis: CMakeFiles/Solar_analysis.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/thinkcomputational/CLionProjects/Solar_water_heater_withoutcover/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_7) "Linking CXX executable Solar_analysis"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/Solar_analysis.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/Solar_analysis.dir/build: Solar_analysis
.PHONY : CMakeFiles/Solar_analysis.dir/build

CMakeFiles/Solar_analysis.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/Solar_analysis.dir/cmake_clean.cmake
.PHONY : CMakeFiles/Solar_analysis.dir/clean

CMakeFiles/Solar_analysis.dir/depend:
	cd /home/thinkcomputational/CLionProjects/Solar_water_heater_withoutcover/cmake-build-debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/thinkcomputational/CLionProjects/Solar_water_heater_withoutcover /home/thinkcomputational/CLionProjects/Solar_water_heater_withoutcover /home/thinkcomputational/CLionProjects/Solar_water_heater_withoutcover/cmake-build-debug /home/thinkcomputational/CLionProjects/Solar_water_heater_withoutcover/cmake-build-debug /home/thinkcomputational/CLionProjects/Solar_water_heater_withoutcover/cmake-build-debug/CMakeFiles/Solar_analysis.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/Solar_analysis.dir/depend

