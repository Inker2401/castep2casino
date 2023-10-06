# Usage: make <target> where <target> indicates any one of the following
#   casino2castep : compile the program
#   clean   : delete executable and object (compiled source) files to prepare for rebuilding from scratch


# --------------------------------------------------------
# Set Fortran Compiler - one of gfortran, ifort
F90 := gfortran
# Set FFT library - currently only fftw3 supported
FFTLIB:=fftw3
# Set optimisation - fast/debug
BUILD := debug
# --------------------------------------------------------
# ************* END OF EDITABLE SECTION ******************

# Set the source and object directories
SRC_DIR:=./src
OBJ_DIR:=./obj
# Set the modules and dependencies
MAIN_FILE = $(SRC_DIR)/main.f90
MODULES = constants.f90 math.f90 io.f90 latt.f90 basis.f90 density.f90 casino.f90

# Strip whitespaces
F90 := $(strip $(F90))
BUILD := $(strip $(BUILD))

# Check if valid parameters have been specified
VALID_F90=gfortran ifort
VALID_BUILD=fast debug
VALID_FFT=fftw3

$(if $(filter-out $(VALID_F90), $(F90)), \
	$(error F90 should be one of the following supported compilers:  $(VALID_F90), not $(F90))\
)
$(if $(filter-out $(VALID_BUILD), $(BUILD)), \
	$(error BUILD should be one of: $(VALID_BUILD), not $(BUILD) )\
)
$(if $(filter-out $(VALID_FFT), $(FFTLIB)), \
	$(error FFTLIB should be one of: $(VALID_FFT), not $(FFTLIB) )\
)

# Set compiler optimisation flags
ifeq ($(F90), gfortran)
    FFLAGS_FAST = -O2 -fimplicit-none -std=f2018 -pedantic -fcheck=do,bounds -Wline-truncation
    FFLAGS_DEBUG = -O0 -g -fimplicit-none -std=f2018 -Wall -Wconversion -Wcharacter-truncation -Wdo-subscript -Wsurprising -Wunused -fcheck=all -fbacktrace
    MOD_FLAGS = -I$(OBJ_DIR) -J$(OBJ_DIR)
else ifeq ($(F90), ifort)
    FFLAGS_FAST = -O2 -standard-semantics
    FFLAGS_DEBUG = -O0 -g -standard-semantics -debug extended -traceback
    MOD_FLAGS = -I$(OBJ_DIR) -module $(OBJ_DIR)
endif

ifeq ($(BUILD), fast)
FFLAGS = $(FFLAGS_FAST)
else ifeq ($(BUILD), debug)
FFLAGS = $(FFLAGS_DEBUG)
endif

# Set FFT flags
FFTFLAG = -l$(FFTLIB)

export F90
export BUILD
export FFLAGS_FAST
export FFLAGS_DEBUG
export FFTFLAG

# Set location of module files
MOD_FILES = $(addprefix $(OBJ_DIR)/, $(MODULES))
OBJ_FILES = $(addprefix $(OBJ_DIR)/, $(MODULES:.f90=.o))
EXECUTABLE = casino2castep.e

# Set the actual compile command - Do not modify directly
COMPILE.F90 = $(F90) $(FFLAGS) $(FFTFLAG) $(MOD_FLAGS)

# Specify targets
.phony : all casino2castep clean

all : casino2castep

casino2castep : $(OBJ_FILES) $(MAIN_FILE)
	$(COMPILE.F90) $^ -o $(EXECUTABLE)

clean :
	rm -rf $(EXECUTABLE) $(OBJ_DIR)

# Make directory if it does not already exist
$(OBJ_DIR) :
	mkdir -p $(OBJ_DIR)

# Build the modules
$(OBJ_DIR)/%.o : $(SRC_DIR)/%.f90 | $(OBJ_DIR)
	$(COMPILE.F90) -c $< -o $@
