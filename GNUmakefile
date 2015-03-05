CASTRO_DIR ?= /path/to/Castro

PRECISION        = DOUBLE
PROFILE          = FALSE
DEBUG            = FALSE
DIM              = 2
COMP	         = g++
FCOMP	         = gfortran

USE_MPI          = TRUE
USE_OMP          = FALSE

USE_GRAV         = TRUE
USE_REACT        = FALSE
USE_MODELPARSER  = TRUE

CASTRO_DIR = ../..

# This sets the EOS directory in $(CASTRO_DIR)/EOS
EOS_dir     := gamma_law_general

# This sets the network directory in $(NETWORK_HOME)
Network_dir := general_null
GENERAL_NET_INPUTS := ./planet.net

Bpack   := ./Make.package
Blocs   := .

include ../Make.Castro
