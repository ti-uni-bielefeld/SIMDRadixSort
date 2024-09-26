#===========================================================================
#
# Makefile
#
# This file is part of the following software:
# 
#    - the low-level C++ template SIMD library
#    - the SIMD implementation of the MinWarping and the 2D-Warping methods 
#      for local visual homing.
# 
# The software is provided based on the accompanying license agreement in the
# file LICENSE.md.
# The software is provided "as is" without any warranty by the licensor and
# without any liability of the licensor, and the software may not be
# distributed by the licensee; see the license agreement for details.
# 
# (C) Ralf Möller
#     Computer Engineering
#     Faculty of Technology
#     Bielefeld University
#     www.ti.uni-bielefeld.de
#
#===========================================================================

binaries = simdRadixSortGeneric

CXX ?= g++

# flags:
# - save assembly code: -masm=intel -save-temps -fverbose-asm
# - profiling: -pg for gprof

# warnings recommended by:
# https://github.com/cpp-best-practices/cppbestpractices/blob/master/02-Use_the_Tools_Available.md#gcc--clang
# (also has explanations for the flags)
warning_flags = -Wall -Wextra -pedantic -Wpedantic \
	-Wnon-virtual-dtor -Wcast-align -Wunused -Woverloaded-virtual \
	-Wsign-compare -Wmisleading-indentation -Wformat=2 \
	-Wimplicit-fallthrough \
	-Wno-format-nonliteral -Wno-cast-align

flags_avx512 = -mavx512f -mavx512bw -mavx512dq -mavx512vl -mpopcnt

# optimization flags
# -funroll-loops is faster in both min-warping phases
optflags ?= -O3 -funroll-loops

flags = $(warning_flags) -fno-var-tracking $(flags_avx512) -pthread -ggdb -std=c++11 $(optflags)

.PHONY: all
all: $(binaries)

%.o: %.C
	$(CXX) -MMD -MP $(flags) -c $< -o $@
$(binaries): %: %.o
	$(CXX) $(flags) -o $@ $< $(libraries)

# os dependent definitions
ifeq ($(OS),Windows_NT) 
RM = del /Q /F
NULL = NUL
else
RM = rm -f
NULL = /dev/null
endif

.PHONY: clean
clean:
	$(RM) $(binaries) *.o *.d *~ *.tmp >$(NULL) 2>&1

.PHONY: info
info:
	@echo "vector extensions:"
	@echo | $(CXX) -xc++ -dM -E - $(flags) | grep -e SSE -e AVX -e POPCNT | sort
	@echo "compiler:" $(CXX)
	@echo "flags:   " $(flags)
	@echo "binaries:" $(binaries)

-include *.d
