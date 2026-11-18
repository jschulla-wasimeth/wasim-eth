# WaSiM makefile for Linux and macOS
#
#   make                      OpenMP build             -> ./wasim
#   make MPI=1                MPI + OpenMP build       -> ./wasim-mpi
#   make NCXX=system          use the system netCDF C++ library (ncxx4-config) instead of src/netcdfcpp_Linux
#   make ARCH=-march=native   CPU-specific tuning (changes results slightly, e.g. through FMA; see tests/golden)
#   make DEBUG=1              -O0 -g, no stripping
#   make clean                remove build products (also those of the legacy makefile_* files)
#
# Requirements: a C++11 compiler with OpenMP, netCDF-C (nc-config), for MPI=1 an MPI compiler wrapper (mpicxx).
# On macOS: Homebrew (brew install netcdf llvm [open-mpi] [netcdf-cxx]); uses Homebrew's LLVM clang if installed,
# otherwise (or with APPLE_CLANG=1) Apple clang with brew install libomp.
# Variables can be overridden on the command line, e.g. make CXX=g++-14 MPICXX=mpiicpc NETCDF_LIBS=...

MPI   ?= 0
NCXX  ?= bundled
DEBUG ?= 0
ARCH  ?=
OS    := $(shell uname -s)

# --- compiler and OpenMP -------------------------------------------------------------------------
ifeq ($(OS),Darwin)
  PLATFORM_DEF = -D MACOS64
  BREW ?= $(shell brew --prefix 2>/dev/null || echo /opt/homebrew)
  # Default: Homebrew's LLVM clang, which ships its own OpenMP runtime, so the flags are the same
  # as for clang on Linux. Without it (or with APPLE_CLANG=1): Apple clang + Homebrew libomp.
  LLVM ?= $(shell brew --prefix llvm 2>/dev/null || echo $(BREW)/opt/llvm)
  APPLE_CLANG ?= $(if $(wildcard $(LLVM)/bin/clang++),0,1)
  ifeq ($(APPLE_CLANG),0)
    ifeq ($(origin CXX),default)
      CXX = $(LLVM)/bin/clang++
      CXXTAG = brew-llvm
    endif
    OMPFLAGS ?= -fopenmp
    # use LLVM's libomp and libc++ at run time, not the system's (see `brew info llvm`)
    OMPLIBS  ?= -fopenmp -L$(LLVM)/lib -Wl,-rpath,$(LLVM)/lib -L$(LLVM)/lib/c++ -Wl,-rpath,$(LLVM)/lib/c++
  else
    LIBOMP ?= $(shell brew --prefix libomp 2>/dev/null || echo $(BREW)/opt/libomp)
    ifeq ($(origin CXX),default)
      CXX = clang++
      CXXTAG = apple-clang
    endif
    # Apple clang has no OpenMP runtime of its own
    OMPFLAGS ?= -Xpreprocessor -fopenmp -I$(LIBOMP)/include
    OMPLIBS  ?= -L$(LIBOMP)/lib -lomp
  endif
else
  PLATFORM_DEF = -D LINUX64
  ifeq ($(origin CXX),default)
    CXX = g++
  endif
  OMPFLAGS ?= -fopenmp
  OMPLIBS  ?= -fopenmp
endif

ifeq ($(MPI),1)
  MPICXX ?= mpicxx
  # Open MPI and MPICH wrappers honour these to pick the underlying compiler
  export OMPI_CXX = $(CXX)
  export MPICH_CXX = $(CXX)
  COMPILER = $(MPICXX)
  MPIDEFS  = -D _MPI -D OMPI_SKIP_MPICXX
  TARGET  ?= wasim-mpi
  VARIANT  = mpi
else
  COMPILER = $(CXX)
  MPIDEFS  =
  TARGET  ?= wasim
  VARIANT  = omp
endif

# --- netCDF --------------------------------------------------------------------------------------
NETCDF_CFLAGS ?= $(shell nc-config --cflags 2>/dev/null)
NETCDF_LIBS   ?= $(shell nc-config --libs 2>/dev/null || echo -lnetcdf)

ifeq ($(NCXX),system)
  NCXX_CFLAGS ?= -I$(shell ncxx4-config --includedir 2>/dev/null)
  NCXX_LIBS   ?= $(shell ncxx4-config --libs 2>/dev/null || echo -lnetcdf_c++4)
  NCXX_SRC     =
else
  # bundled copy of netcdf-cxx4 (an older release, see TODO.md)
  NCXX_DIR     = netcdfcpp_Linux
  NCXX_CFLAGS  = -I$(NCXX_DIR)
  NCXX_LIBS    =
  NCXX_SRC     = $(wildcard $(NCXX_DIR)/*.cc)
endif

# --- flags ---------------------------------------------------------------------------------------
ifeq ($(DEBUG),1)
  OPT = -O0 -g
else
  OPT = -O3
endif

DEFS     = $(PLATFORM_DEF) $(MPIDEFS) -D _OMP3 -D _NETCDF -D USE_NETCDF4
CXXFLAGS = -std=c++11 $(OPT) $(ARCH) $(DEFS) $(OMPFLAGS) $(NETCDF_CFLAGS) $(NCXX_CFLAGS) -Wno-unused-result
LDLIBS   = $(NCXX_LIBS) $(NETCDF_LIBS) $(OMPLIBS) -lm

# --- sources -------------------------------------------------------------------------------------
SRC = abfluss.cc cspline.cc dataassim.cc defltuzr.cc etpetr.cc expokorr.cc gifreader.cc gwflow.cc \
      interzep.cc irrigate.cc meteoint.cc modapuzr.cc modefuzr.cc netcdfhlp.cc parser.cc siltingup.cc \
      snow_canopy.cc snowmod.cc stringp.cc surfrout.cc teilmod.cc unsatzon.cc wasimuzr.cc wsstream.cc \
      $(NCXX_SRC)

CXXTAG  ?= $(notdir $(firstword $(CXX)))
BUILDDIR = build/$(VARIANT)-$(NCXX)-$(CXXTAG)
OBJ      = $(SRC:%.cc=$(BUILDDIR)/%.o)
DEP      = $(OBJ:.o=.d)

.PHONY: all clean help

all: $(TARGET)

$(TARGET): $(OBJ)
	$(COMPILER) $(CXXFLAGS) $(OBJ) -o $@ $(LDLIBS)
ifneq ($(DEBUG),1)
	strip $@
endif

$(BUILDDIR)/%.o: %.cc
	@mkdir -p $(dir $@)
	$(COMPILER) $(CXXFLAGS) -MMD -MP -c $< -o $@

-include $(DEP)

help:
	@sed -n '2,14p' Makefile | sed 's/^# \{0,1\}//'

clean:
	rm -rf build
	rm -f wasim wasim-mpi
	rm -f *.o netcdfcpp*/*.o *.a
	rm -f wasimuzr-* wasimtop-*
