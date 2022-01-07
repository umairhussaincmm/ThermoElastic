# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.10

# Default target executed when no arguments are given to make.
default_target: all

.PHONY : default_target

# Allow only one "make -f Makefile2" at a time, but pass parallelism.
.NOTPARALLEL:


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
CMAKE_SOURCE_DIR = /home/ubuntu/CLionProjects/Thermo-Elastic

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/ubuntu/CLionProjects/Thermo-Elastic

#=============================================================================
# Targets provided globally by CMake.

# Special rule for the target rebuild_cache
rebuild_cache:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Running CMake to regenerate build system..."
	/usr/bin/cmake -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR)
.PHONY : rebuild_cache

# Special rule for the target rebuild_cache
rebuild_cache/fast: rebuild_cache

.PHONY : rebuild_cache/fast

# Special rule for the target edit_cache
edit_cache:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Running CMake cache editor..."
	/usr/bin/ccmake -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR)
.PHONY : edit_cache

# Special rule for the target edit_cache
edit_cache/fast: edit_cache

.PHONY : edit_cache/fast

# The main all target
all: cmake_check_build_system
	$(CMAKE_COMMAND) -E cmake_progress_start /home/ubuntu/CLionProjects/Thermo-Elastic/CMakeFiles /home/ubuntu/CLionProjects/Thermo-Elastic/CMakeFiles/progress.marks
	$(MAKE) -f CMakeFiles/Makefile2 all
	$(CMAKE_COMMAND) -E cmake_progress_start /home/ubuntu/CLionProjects/Thermo-Elastic/CMakeFiles 0
.PHONY : all

# The main clean target
clean:
	$(MAKE) -f CMakeFiles/Makefile2 clean
.PHONY : clean

# The main clean target
clean/fast: clean

.PHONY : clean/fast

# Prepare targets for installation.
preinstall: all
	$(MAKE) -f CMakeFiles/Makefile2 preinstall
.PHONY : preinstall

# Prepare targets for installation.
preinstall/fast:
	$(MAKE) -f CMakeFiles/Makefile2 preinstall
.PHONY : preinstall/fast

# clear depends
depend:
	$(CMAKE_COMMAND) -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR) --check-build-system CMakeFiles/Makefile.cmake 1
.PHONY : depend

#=============================================================================
# Target rules for targets named info

# Build rule for target.
info: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 info
.PHONY : info

# fast build rule for target.
info/fast:
	$(MAKE) -f CMakeFiles/info.dir/build.make CMakeFiles/info.dir/build
.PHONY : info/fast

#=============================================================================
# Target rules for targets named strip_comments

# Build rule for target.
strip_comments: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 strip_comments
.PHONY : strip_comments

# fast build rule for target.
strip_comments/fast:
	$(MAKE) -f CMakeFiles/strip_comments.dir/build.make CMakeFiles/strip_comments.dir/build
.PHONY : strip_comments/fast

#=============================================================================
# Target rules for targets named run

# Build rule for target.
run: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 run
.PHONY : run

# fast build rule for target.
run/fast:
	$(MAKE) -f CMakeFiles/run.dir/build.make CMakeFiles/run.dir/build
.PHONY : run/fast

#=============================================================================
# Target rules for targets named Thermo-Elastic

# Build rule for target.
Thermo-Elastic: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 Thermo-Elastic
.PHONY : Thermo-Elastic

# fast build rule for target.
Thermo-Elastic/fast:
	$(MAKE) -f CMakeFiles/Thermo-Elastic.dir/build.make CMakeFiles/Thermo-Elastic.dir/build
.PHONY : Thermo-Elastic/fast

#=============================================================================
# Target rules for targets named debug

# Build rule for target.
debug: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 debug
.PHONY : debug

# fast build rule for target.
debug/fast:
	$(MAKE) -f CMakeFiles/debug.dir/build.make CMakeFiles/debug.dir/build
.PHONY : debug/fast

#=============================================================================
# Target rules for targets named release

# Build rule for target.
release: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 release
.PHONY : release

# fast build rule for target.
release/fast:
	$(MAKE) -f CMakeFiles/release.dir/build.make CMakeFiles/release.dir/build
.PHONY : release/fast

#=============================================================================
# Target rules for targets named runclean

# Build rule for target.
runclean: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 runclean
.PHONY : runclean

# fast build rule for target.
runclean/fast:
	$(MAKE) -f CMakeFiles/runclean.dir/build.make CMakeFiles/runclean.dir/build
.PHONY : runclean/fast

#=============================================================================
# Target rules for targets named distclean

# Build rule for target.
distclean: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 distclean
.PHONY : distclean

# fast build rule for target.
distclean/fast:
	$(MAKE) -f CMakeFiles/distclean.dir/build.make CMakeFiles/distclean.dir/build
.PHONY : distclean/fast

Initial_values.o: Initial_values.cpp.o

.PHONY : Initial_values.o

# target to build an object file
Initial_values.cpp.o:
	$(MAKE) -f CMakeFiles/Thermo-Elastic.dir/build.make CMakeFiles/Thermo-Elastic.dir/Initial_values.cpp.o
.PHONY : Initial_values.cpp.o

Initial_values.i: Initial_values.cpp.i

.PHONY : Initial_values.i

# target to preprocess a source file
Initial_values.cpp.i:
	$(MAKE) -f CMakeFiles/Thermo-Elastic.dir/build.make CMakeFiles/Thermo-Elastic.dir/Initial_values.cpp.i
.PHONY : Initial_values.cpp.i

Initial_values.s: Initial_values.cpp.s

.PHONY : Initial_values.s

# target to generate assembly for a file
Initial_values.cpp.s:
	$(MAKE) -f CMakeFiles/Thermo-Elastic.dir/build.make CMakeFiles/Thermo-Elastic.dir/Initial_values.cpp.s
.PHONY : Initial_values.cpp.s

ThermoElastic.o: ThermoElastic.cpp.o

.PHONY : ThermoElastic.o

# target to build an object file
ThermoElastic.cpp.o:
	$(MAKE) -f CMakeFiles/Thermo-Elastic.dir/build.make CMakeFiles/Thermo-Elastic.dir/ThermoElastic.cpp.o
.PHONY : ThermoElastic.cpp.o

ThermoElastic.i: ThermoElastic.cpp.i

.PHONY : ThermoElastic.i

# target to preprocess a source file
ThermoElastic.cpp.i:
	$(MAKE) -f CMakeFiles/Thermo-Elastic.dir/build.make CMakeFiles/Thermo-Elastic.dir/ThermoElastic.cpp.i
.PHONY : ThermoElastic.cpp.i

ThermoElastic.s: ThermoElastic.cpp.s

.PHONY : ThermoElastic.s

# target to generate assembly for a file
ThermoElastic.cpp.s:
	$(MAKE) -f CMakeFiles/Thermo-Elastic.dir/build.make CMakeFiles/Thermo-Elastic.dir/ThermoElastic.cpp.s
.PHONY : ThermoElastic.cpp.s

assemble_system.o: assemble_system.cpp.o

.PHONY : assemble_system.o

# target to build an object file
assemble_system.cpp.o:
	$(MAKE) -f CMakeFiles/Thermo-Elastic.dir/build.make CMakeFiles/Thermo-Elastic.dir/assemble_system.cpp.o
.PHONY : assemble_system.cpp.o

assemble_system.i: assemble_system.cpp.i

.PHONY : assemble_system.i

# target to preprocess a source file
assemble_system.cpp.i:
	$(MAKE) -f CMakeFiles/Thermo-Elastic.dir/build.make CMakeFiles/Thermo-Elastic.dir/assemble_system.cpp.i
.PHONY : assemble_system.cpp.i

assemble_system.s: assemble_system.cpp.s

.PHONY : assemble_system.s

# target to generate assembly for a file
assemble_system.cpp.s:
	$(MAKE) -f CMakeFiles/Thermo-Elastic.dir/build.make CMakeFiles/Thermo-Elastic.dir/assemble_system.cpp.s
.PHONY : assemble_system.cpp.s

main.o: main.cpp.o

.PHONY : main.o

# target to build an object file
main.cpp.o:
	$(MAKE) -f CMakeFiles/Thermo-Elastic.dir/build.make CMakeFiles/Thermo-Elastic.dir/main.cpp.o
.PHONY : main.cpp.o

main.i: main.cpp.i

.PHONY : main.i

# target to preprocess a source file
main.cpp.i:
	$(MAKE) -f CMakeFiles/Thermo-Elastic.dir/build.make CMakeFiles/Thermo-Elastic.dir/main.cpp.i
.PHONY : main.cpp.i

main.s: main.cpp.s

.PHONY : main.s

# target to generate assembly for a file
main.cpp.s:
	$(MAKE) -f CMakeFiles/Thermo-Elastic.dir/build.make CMakeFiles/Thermo-Elastic.dir/main.cpp.s
.PHONY : main.cpp.s

output_results.o: output_results.cpp.o

.PHONY : output_results.o

# target to build an object file
output_results.cpp.o:
	$(MAKE) -f CMakeFiles/Thermo-Elastic.dir/build.make CMakeFiles/Thermo-Elastic.dir/output_results.cpp.o
.PHONY : output_results.cpp.o

output_results.i: output_results.cpp.i

.PHONY : output_results.i

# target to preprocess a source file
output_results.cpp.i:
	$(MAKE) -f CMakeFiles/Thermo-Elastic.dir/build.make CMakeFiles/Thermo-Elastic.dir/output_results.cpp.i
.PHONY : output_results.cpp.i

output_results.s: output_results.cpp.s

.PHONY : output_results.s

# target to generate assembly for a file
output_results.cpp.s:
	$(MAKE) -f CMakeFiles/Thermo-Elastic.dir/build.make CMakeFiles/Thermo-Elastic.dir/output_results.cpp.s
.PHONY : output_results.cpp.s

run.o: run.cpp.o

.PHONY : run.o

# target to build an object file
run.cpp.o:
	$(MAKE) -f CMakeFiles/Thermo-Elastic.dir/build.make CMakeFiles/Thermo-Elastic.dir/run.cpp.o
.PHONY : run.cpp.o

run.i: run.cpp.i

.PHONY : run.i

# target to preprocess a source file
run.cpp.i:
	$(MAKE) -f CMakeFiles/Thermo-Elastic.dir/build.make CMakeFiles/Thermo-Elastic.dir/run.cpp.i
.PHONY : run.cpp.i

run.s: run.cpp.s

.PHONY : run.s

# target to generate assembly for a file
run.cpp.s:
	$(MAKE) -f CMakeFiles/Thermo-Elastic.dir/build.make CMakeFiles/Thermo-Elastic.dir/run.cpp.s
.PHONY : run.cpp.s

setup_system.o: setup_system.cpp.o

.PHONY : setup_system.o

# target to build an object file
setup_system.cpp.o:
	$(MAKE) -f CMakeFiles/Thermo-Elastic.dir/build.make CMakeFiles/Thermo-Elastic.dir/setup_system.cpp.o
.PHONY : setup_system.cpp.o

setup_system.i: setup_system.cpp.i

.PHONY : setup_system.i

# target to preprocess a source file
setup_system.cpp.i:
	$(MAKE) -f CMakeFiles/Thermo-Elastic.dir/build.make CMakeFiles/Thermo-Elastic.dir/setup_system.cpp.i
.PHONY : setup_system.cpp.i

setup_system.s: setup_system.cpp.s

.PHONY : setup_system.s

# target to generate assembly for a file
setup_system.cpp.s:
	$(MAKE) -f CMakeFiles/Thermo-Elastic.dir/build.make CMakeFiles/Thermo-Elastic.dir/setup_system.cpp.s
.PHONY : setup_system.cpp.s

solve.o: solve.cpp.o

.PHONY : solve.o

# target to build an object file
solve.cpp.o:
	$(MAKE) -f CMakeFiles/Thermo-Elastic.dir/build.make CMakeFiles/Thermo-Elastic.dir/solve.cpp.o
.PHONY : solve.cpp.o

solve.i: solve.cpp.i

.PHONY : solve.i

# target to preprocess a source file
solve.cpp.i:
	$(MAKE) -f CMakeFiles/Thermo-Elastic.dir/build.make CMakeFiles/Thermo-Elastic.dir/solve.cpp.i
.PHONY : solve.cpp.i

solve.s: solve.cpp.s

.PHONY : solve.s

# target to generate assembly for a file
solve.cpp.s:
	$(MAKE) -f CMakeFiles/Thermo-Elastic.dir/build.make CMakeFiles/Thermo-Elastic.dir/solve.cpp.s
.PHONY : solve.cpp.s

# Help Target
help:
	@echo "The following are some of the valid targets for this Makefile:"
	@echo "... all (the default if no target is provided)"
	@echo "... clean"
	@echo "... depend"
	@echo "... rebuild_cache"
	@echo "... edit_cache"
	@echo "... info"
	@echo "... strip_comments"
	@echo "... run"
	@echo "... Thermo-Elastic"
	@echo "... debug"
	@echo "... release"
	@echo "... runclean"
	@echo "... distclean"
	@echo "... Initial_values.o"
	@echo "... Initial_values.i"
	@echo "... Initial_values.s"
	@echo "... ThermoElastic.o"
	@echo "... ThermoElastic.i"
	@echo "... ThermoElastic.s"
	@echo "... assemble_system.o"
	@echo "... assemble_system.i"
	@echo "... assemble_system.s"
	@echo "... main.o"
	@echo "... main.i"
	@echo "... main.s"
	@echo "... output_results.o"
	@echo "... output_results.i"
	@echo "... output_results.s"
	@echo "... run.o"
	@echo "... run.i"
	@echo "... run.s"
	@echo "... setup_system.o"
	@echo "... setup_system.i"
	@echo "... setup_system.s"
	@echo "... solve.o"
	@echo "... solve.i"
	@echo "... solve.s"
.PHONY : help



#=============================================================================
# Special targets to cleanup operation of make.

# Special rule to run CMake to check the build system integrity.
# No rule that depends on this can have commands that come from listfiles
# because they might be regenerated.
cmake_check_build_system:
	$(CMAKE_COMMAND) -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR) --check-build-system CMakeFiles/Makefile.cmake 0
.PHONY : cmake_check_build_system

