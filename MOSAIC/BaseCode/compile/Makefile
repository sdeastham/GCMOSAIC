#######################################################################
#
# Makefile for mosaic.x
#
#######################################################################


#COMPILER = $(COMP)
 COMPILER = pgf90
#foptions = $(FLAGS)
 foptions = -g -C -byteswapio -Ktrap=fp -O0
 COMP_FLAGS_FIXED =

RM = /bin/rm


#
# Suffix rules
#
.SUFFIXES : .f90 .f .o .x

.f90.o:
	$(COMPILER)  -c $(foptions) $(COMP_FLAGS_FREE)  $<
.f.o:
	$(COMPILER)  -c $(foptions) $(COMP_FLAGS_FIXED) $<


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

MAIN =  \
	domassbalance.o  \
	function_nbllen.o  \
	function_watervapor.o  \
	init_data_modules.o  \
	integratechemistry.o  \
	main.o  \
	parametric_analysis_mode.o  \
	module_print_aer.o  \
	print_gas.o  \
	readinputfile.o  \
	setaircomposition.o  \
	setiofiles.o  \
	setrunparameters.o  \
	time_integration_mode.o  \
	updateemissions.o  \
	updatemetfields.o  \
	updatetime.o \
	mosaic_allocate_memory.o \
	module_mosaic_support.o \
	module_mosaic_init.o


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

$(CLOUD):  $(DATAMODU)

$(GAS):  $(DATAMODU)

#$(SOLVER):  $(DATAMODU)

$(MAIN):  $(DATAMODU)


mosaic.x:  $(DATAMODU)  $(PMCMOS)  $(AEROSOL)  $(CLOUD)  $(GAS)  $(SOLVER)  $(MAIN)
	$(COMPILER)  $(foptions)  -o $@  $(DATAMODU)  $(PMCMOS)  $(AEROSOL)  $(CLOUD)  $(GAS)  $(SOLVER)  $(MAIN)


clean:
	$(RM)  mosaic.x  *.o  *.mod



