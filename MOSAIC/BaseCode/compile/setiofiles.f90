subroutine SetIOfiles

  use module_data_mosaic_main,   only:lun_sect_170, lun_sect_171, &
       lun_sect_172, lun_sect_180, lun_sect_183, lun_sect_184,    &
       lun_sect_185, lun_sect_186, lun_sect_188, lun_sect_190,    &
       inputfile, lun_inp, lun_gas, lun_drydist, lun_wetdist,     &
       lun_species, lun_aeroptic, lun_aer_status, lun_aer,        &
       lun_fullout,                                               &
       iwrite_aer_bin, mgas, gas_output, maer, drydist_output,    &
       wetdist_output, species_output, aeroptic_output, maeroptic,&
       aer_output, fullout_fname
  use module_data_mosaic_gas,    only:
  use module_data_mosaic_aero,   only: mYes, nbin_a
  use module_data_mosaic_cloud,  only:
  use module_data_mosaic_pmcmos, only: lun_pmcmos_in1,            &
       lun_pmcmos_in2
  use module_data_mosaic_asect,  only: lunout

  implicit none

  !Local Variables
  character(len=7) ::  bb

  integer :: ibin
  integer :: nchar, nbllen, nbb


  ! INPUT FILES
  !BSINGH - temporarily read case name from a file
  !open(100,file = 'b_inp.txt',status='old') !BSINGH Uncomment this line BSINGH
  write(6,*)'Enter gas input filename. Example: case1.inp'
  read(5,*)inputfile  !BSINGH - Commented out this line
  !read(100,*)inputfile !BSINGH - reading from the file  -!BSINGH Uncomment this line BSINGH
  !      inputfile = 'test.inp'

  lun_pmcmos_in1 = 16 ; lun_pmcmos_in2 = 17

  lun_sect_170 = 170
  lun_sect_171 = 171
  lun_sect_172 = 172
  lun_sect_180 = 180
  lun_sect_183 = 183
  lun_sect_184 = 184
  lun_sect_185 = 185
  lun_sect_186 = 186
  lun_sect_188 = 188
  lun_sect_190 = 190

  lunout = lun_sect_170

  lun_inp = 10
  open(lun_inp, file = inputfile)
  call ReadInputFile
  close(lun_inp)


  !------------------------------------------------------------------
  ! OUTPUT FILES
  !
  nchar = nbllen(inputfile) - 4

  ! GAS output file
  if(mgas .eq. mYES)then

     lun_gas = 19
     gas_output = 'output/'//inputfile(1:nchar)//'.gas.txt'
     open(lun_gas, file = gas_output)

  endif




  ! AEROSOL output files
  if(maer .eq. mYES)then        ! UNCOMMENT THIS LINE
     !      if(1 .eq. 0)then            ! COMMENT THIS LINE

     ! dry aerosol number, area, and mass distribution output file
     lun_drydist = 20
     drydist_output='output/'//inputfile(1:nchar)//'.dist.dry.txt'
     open(lun_drydist, file = drydist_output)

     ! wet aerosol number, area, and mass distribution output file
     lun_wetdist = 21
     wetdist_output='output/'//inputfile(1:nchar)//'.dist.wet.txt'
     open(lun_wetdist, file = wetdist_output)

     ! species distribution output file
     lun_species = 22
     species_output=   &
          'output/'//inputfile(1:nchar)//'.dist.species.txt'
     open(lun_species, file = species_output)

     ! aerosol optical info output file
     lun_aeroptic = 23
     aeroptic_output=   &
          'output/'//inputfile(1:nchar)//'.aeroptic.txt'
     if (maeroptic > 0) then
        open(lun_aeroptic, file = aeroptic_output)
     end if

     ! "fullout" output file with all gas and aerosol species
     lun_fullout = 24
     fullout_fname=   &
          'output/'//inputfile(1:nchar)//'.fullout.txt'
     open(lun_fullout, file = fullout_fname)


     ! aerosol bin output files
     do ibin = 1, nbin_a
        write(bb, '(a,i3.3)')'.bin',ibin
        aer_output(ibin) =   &
             'output/'//inputfile(1:nchar)//bb//'.txt'
        lun_aer_status(ibin) = 0
        if (nbin_a <= 20) then
           lun_aer(ibin) = 30 + (ibin-1)
           if (iwrite_aer_bin > 0) &
                open( lun_aer(ibin), file=aer_output(ibin) )
        else
           ! open/write/close these files on each entry to avoid
           ! problems with lun_aer being too large
           lun_aer(ibin) = 30
        endif
     enddo

  endif

  return
end subroutine SetIOfiles


!---------------------------------------------------------------------
subroutine lun_aer_open(ibin)
  use module_data_mosaic_aero
  use module_data_mosaic_main
  implicit none
  integer, intent(in) :: ibin

  if (nbin_a <= 20) return
  if ((ibin < 1) .or. (ibin > nbin_a)) then
     write(*,*) '*** lun_aer_open -- bad ibin = ', ibin
     stop
  end if

  ! open/write/close these files on each entry to avoid
  ! problems with lun_aer being too large
  if (lun_aer_status(ibin) <= 0) then
     open( unit=lun_aer(ibin), &
          file=aer_output(ibin), status='unknown' )
     lun_aer_status(ibin) = 1
  else
     open( unit=lun_aer(ibin), &
          file=aer_output(ibin), status='old', position='append' )
  endif
  return
end subroutine lun_aer_open


!---------------------------------------------------------------------
subroutine lun_aer_close(ibin)
  use module_data_mosaic_main
  use module_data_mosaic_aero
  implicit none
  integer, intent(in) :: ibin

  if (nbin_a <= 20) return
  if ((ibin < 1) .or. (ibin > nbin_a)) then
     write(*,*) '*** lun_aer_close -- bad ibin = ', ibin
     stop
  end if

  ! open/write/close these files on each entry to avoid
  ! problems with lun_aer being too large
  close( unit=lun_aer(ibin) )
  return
end subroutine lun_aer_close


