#------------------------------------------------------------------------------
#                  GEOS-Chem Global Chemical Transport Model                  #
#------------------------------------------------------------------------------
#BOP
#
# !MODULE: Makefile (in the MOSAIC/ subdirectory)
#
# !DESCRIPTION: This makefile compiles the MOSAICIA code.  
#  Object files (*.o) are bundled into the libMOSAIC.a library 
#  (located in the LIB directory).  Module files (*.mod) are copied to 
#  the MOD directory. 
#\\
#\\
# !REMARKS:
# To build the programs, call "make" with the following syntax:
#                                                                             .
#   make -jN TARGET REQUIRED-FLAGS [ OPTIONAL-FLAGS ]
#                                                                             .
# To display a complete list of options, type "make help".
#                                                                             .
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%% NOTE: Normally you will not have to call this Makefile directly,     %%%
# %%% it will be called automatically from the Makefile in the directory   %%%
# %%% just above this one!                                                 %%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#                                                                             .
# Makefile uses the following variables:
#                                                                             .
# Variable   Description
# --------   -----------
# SHELL      Specifies the shell for "make" to use (usually SHELL=/bin/sh)
# ROOTDIR    Specifies the root-level directory of the GEOS-Chem code
# HDR        Specifies the directory where GEOS-Chem include files are found
# LIB        Specifies the directory where library files (*.a) are stored
# MOD        Specifies the directory where module files (*.mod) are stored
# AR         Sys var w/ name of library creator program (i.e., "ar", "ranlib")
# MAKE       Sys var w/ name of Make command (i.e, "make" or "gmake")
# R8         Specifies the command to treat "REAL" as "REAL*8"
#
# !REVISION HISTORY: 
#  14 Nov 2014 - S. D. Eastham-Initial version
#EOP
#------------------------------------------------------------------------------
#BOC

###############################################################################
###                                                                         ###
###  Initialization section                                                 ###
###                                                                         ###
###############################################################################

# Define variables
SHELL :=/bin/bash
ROOT  :=..
HELP  :=$(ROOT)/help
LIB   :=$(ROOT)/lib
MOD   :=$(ROOT)/mod

# Include header file.  This returns variables CC, F90, FREEFORM, LD, R8,
# as well as the default Makefile compilation rules for source code files.
include $(ROOT)/Makefile_header.mk

# List of source files
#SRC   :=$(wildcard *.F) $(wildcard *.F90)

# We need to specialize this - only certain files are needed


# Replace .F and .F90 extensions with *.o
TMP   :=$(SRC:.F=.o)
OBJ   :=$(TMP:.F90=.o)

###############################################################################
###                                                                         ###
###  Makefile targets: type "make help" for a complete listing!             ###
###                                                                         ###
###############################################################################

.PHONY: clean help

lib: $(OBJ)
	$(AR) crs libMOSAIC.a $(OBJ)
	mv libMOSAIC.a $(LIB)

clean:
	rm -f *.o *.mod

help:
	@$(MAKE) -C $(HELP)

###############################################################################
###                                                                         ###
###  Dependencies listing                                                   ###
###  (grep "USE " to get the list of module references!)                    ###
###                                                                         ###
###  From this list of dependencies, the "make" utility will figure out     ###
###  correct order of compilation (so we don't have to do that ourselves).  ###
###  This also allows us to compile on multiple processors with "make -j".  ###
###                                                                         ###
###  NOTES:                                                                 ###
###  (1) Only specify object-file dependencies that are within this         ###
###       directory.  Object files in other directories will be referenced  ### 
###       at link-time.                                                     ###
###  (2) For "make -jN" (i.e. compile N files simultaneously), all files    ###
###       in this directory must have a listed dependency.                  ###
###                                                                         ###
###############################################################################

#
# Object lists
#
DATAMODU =  \
	module_data_mosaic_kind.o  \
	module_data_mosaic_constants.o  \
	module_data_mosaic_main.o  \
	module_data_mosaic_aero.o  \
	module_data_mosaic_asect.o \
	module_data_mosaic_gas.o   \
	module_data_mosaic_cloud.o \
	module_data_mosaic_pmcmos.o

PMCMOS =  \
	module_pmcmos_aer.o   \
	module_pmcmos_init.o  \
	module_pmcmos_subaa.o

AEROSOL =  \
	aerchemistry.o  \
	aeroptical.o  \
	ASTEM.o  \
	cputime.o  \
	init_aerosol.o  \
	miedriver.o  \
	module_coag.o  \
	module_coag_3d.o  \
	module_movesect.o  \
	module_movesect_3d.o  \
	module_newnuc.o  \
	module_peg_util.o  \
	module_sect_iface.o  \
	mosaic_box.25.o  \
	module_mosaic_ext.o \
	module_mosaic_lsode.o

CLOUD =  \
	cldchemistry.o

GAS =  \
	function_arr.o  \
	function_troe.o  \
	gaschemistry.o  \
	gasintegrator.o  \
	gasrateconstants_bio.o  \
	gasrateconstants_com.o  \
	gasrateconstants.o  \
	gasrateconstants_het.o  \
	gasrateconstants_mar.o  \
	gasrateconstants_urb.o  \
	gasrates_bio.o  \
	gasrates_com.o  \
	gasrates.o  \
	gasrates_het.o  \
	gasrates_mar.o  \
	gasrates_urb.o  \
	jac.o  \
	loadperoxyparameters.o  \
	mapgas_bio.o  \
	mapgas_com.o  \
	mapgas_mar.o  \
	mapgasspecies.o  \
	mapgas_urb.o  \
	ode_bio.o  \
	ode_com.o  \
	ode_gas.o  \
	ode_het.o  \
	ode_mar.o  \
	ode_urb.o  \
	peroxyrateconstants.o  \
	photoconstants_fixed.o  \
	photoconstants_solar.o  \
	photoparam1.o  \
	photoparam2.o  \
	selectgasregime.o  \
	setgas_bio.o  \
	setgas_com.o  \
	setgasindices.o  \
	setgas_mar.o  \
	setgas_urb.o  \
	solarzenithangle.o

SOLVER =  \
	dlsode.o  \
	dlsodes.o \
	diprep.o  \
	dprjs.o   \
	dsolss.o  \
	cdrv.o    \
	dprep.o   \
	adjlr.o   \
	cntnzu.o  \
	jgroup.o  \
	nnfc.o    \
	nnsc.o    \
	nntc.o    \
	nroc.o    \
	nsfc.o    \
	odrv.o    \
	md.o      \
	sro.o     \
	mdi.o     \
	mdm.o     \
	mdp.o     \
	mdu.o 

module_data_mosaic_constants.o:   module_data_mosaic_kind.o

module_data_mosaic_main.o:   module_data_mosaic_kind.o  module_data_mosaic_constants.o

module_data_mosaic_aero.o:   module_data_mosaic_kind.o  module_data_mosaic_main.o

module_data_mosaic_asect.o:  module_data_mosaic_kind.o  module_data_mosaic_aero.o

module_data_mosaic_cloud.o:  module_data_mosaic_kind.o

module_data_mosaic_gas.o:    module_data_mosaic_kind.o  module_data_mosaic_main.o

module_data_mosaic_pmcmos.o: module_data_mosaic_kind.o  module_data_mosaic_main.o  	module_data_mosaic_aero.o  

init_aerosol.o:              module_pmcmos_aer.o  module_mosaic_init.o $(DATAMODU)

module_coag.o:               module_peg_util.o  $(DATAMODU)

module_coag_3d.o:            module_peg_util.o  $(DATAMODU)

module_movesect.o:           module_peg_util.o  $(DATAMODU)

module_movesect_3d.o:        module_peg_util.o  $(DATAMODU)

module_newnuc.o:             module_peg_util.o  $(DATAMODU)

module_sect_iface.o:         module_coag.o  module_coag_3d.o  module_movesect.o  module_movesect_3d.o  module_newnuc.o  $(DATAMODU)

aerchemistry.o:              module_sect_iface.o  mosaic_box.25.o $(DATAMODU)

module_pmcmos_init.o:        module_pmcmos_subaa.o module_mosaic_init.o $(DATAMODU)

readinputfile.o:             module_pmcmos_init.o  $(DATAMODU)

ASTEM.o:	             module_mosaic_ext.o module_mosaic_support.o module_print_aer.o $(DATAMODU)

mosaic_box.25.o:             ASTEM.o module_mosaic_ext.o module_mosaic_lsode.o $(DATAMODU)

module_mosaic_lsode.o:       module_mosaic_ext.o $(DATAMODU)

module_mosaic_init.o:        module_data_mosaic_kind.o module_data_mosaic_aero.o

$(PMCMOS):  $(DATAMODU)

$(AEROSOL):  $(DATAMODU)
#EOC

