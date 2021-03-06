module module_print_aer
! Aerosol-phase Species
! subroutine for printing output at iprint time steps
!
!--------------------------------------------------------------------
!
  implicit none
  private

  public:: print_aer
contains
      subroutine print_aer(idum, jaerosolstate,isteps_ASTEM,  &
           iter_MESA,aer,gas,electrolyte,mc,num_a,Dp_dry_a,Dp_wet_a, &
           area_dry_a,area_wet_a,mass_wet_a,mass_dry_a,water_a)

      use module_data_mosaic_main
      use module_data_mosaic_aero

      implicit none

      !Subroutine Arguments
      integer, intent(in) :: idum, isteps_ASTEM
      integer, intent(in), dimension(nbin_a_max) :: jaerosolstate,iter_MESA

      real(r8), intent(in), dimension(nbin_a_max) :: num_a,Dp_dry_a,Dp_wet_a,area_dry_a,water_a
      real(r8), intent(in), dimension(nbin_a_max) :: area_wet_a,mass_wet_a,mass_dry_a
      real(r8), intent(in), dimension(ngas_volatile) :: gas
      real(r8), intent(in), dimension(Ncation,nbin_a_max) :: mc
      real(r8), intent(in), dimension(naer,3,nbin_a_max) :: aer
      real(r8), intent(in), dimension(nelectrolyte,3,nbin_a_max) :: electrolyte

      !Local Variables
      integer l, ibin, lbin, je, iv, ifix
      real(r8) :: timeofday, time_dum, cppb(ntot_max), conv1, conv2, conv3,   &
           conv4,   &
           tot_SO4, tot_NO3, tot_Cl, tot_NH4, tot_Na, tot_Ca, tot_H
      real(r8) :: Dpfix_lo(25), Dpfix_hi(25), Dpmean_fix(25), alpha,   &
           num_fix(25), area_fix(25), mass_fix(25)
      real(r8) :: dlogDpdry(nbin_a), dlogDpwet(nbin_a)
      real(r8), allocatable, dimension(:) :: pH
      character*7 diam(nbin_a),numb(nbin_a),area(nbin_a),mass(nbin_a)




      if(.not.allocated(pH)) allocate(pH(nbin_a_max))

      !BSINGH - Following two variables were temporarily initialized, please modify if required
      Dpfix_lo(:) = 0.0
      Dpfix_hi(:) = 0.0

      conv1 = 1.e15/avogad	 ! converts (molec/cc) to (nmol/m^3)
      conv2 = 1./conv1		 ! converts (nmol/m^3) to (molec/cc)
      conv3 = conv2/cair_mlc*ppb ! converts (nmol/m^3) to (ppbv)
!      conv4 = 1.e-3		 ! converts (nmol/m^3) to (umol/m^3)
      conv4 = 1.0		 ! converts (nmol/m^3) to (nmol/m^3)

!      if(idum .eq. 0)then
!        tot_SO4 = aer(iso4_a,jtotal,1) + gas(ih2so4_g)
!        tot_NO3 = aer(ino3_a,jtotal,1) + gas(ihno3_g)
!        tot_Cl  = aer(icl_a ,jtotal,1) + gas(ihcl_g)
!        tot_NH4 = aer(inh4_a,jtotal,1) + gas(inh3_g)
!        tot_Na  = aer(ina_a, jtotal,1)
!        tot_Ca  = aer(ica_a, jtotal,1)
!        tot_H   = (2.*tot_SO4 + tot_NO3 + tot_Cl) - (tot_NH4 + tot_Na)
!
!        write(18,*)'Moles per m^3'
!        write(18,*)'H   = ', tot_H*1.e-9
!        write(18,*)'NH4 = ', tot_NH4*1.e-9
!        write(18,*)'SO4 = ', tot_SO4*1.e-9
!        write(18,*)'NO3 = ', tot_NO3*1.e-9
!        write(18,*)'Cl  = ', tot_Cl*1.e-9
!        write(18,*)'Na  = ', tot_Na*1.e-9
!      endif


      if(idum .eq. 0)then
        time_dum = 0.0
        timeofday = time_UTC_beg
      else
        time_dum = time_hrs
        timeofday = time_UTC
      endif


! output number, surface, and mass distributions
      if(idum .eq. 0)then

! load Dpfix grid
        alpha = 10.**0.1
        Dpfix_lo(1) = 0.005 ! um
        Dpfix_hi(1) = 0.005*alpha
        do ifix = 2, 25
          Dpfix_lo(ifix) = Dpfix_hi(ifix-1)
          Dpfix_hi(ifix) = Dpfix_lo(ifix)*alpha
        enddo

        do ifix = 1, 25
          Dpmean_fix(ifix) = (Dpfix_lo(ifix) + Dpfix_hi(ifix))/2.
!          write(6,660)ifix, Dpfix_lo(ifix), Dpfix_hi(ifix)
        enddo
660     format(i4, f10.5, 2x, f10.5)

! DpgN Lagrangian grid (header)
        do ibin = 1, nbin_a
          write(diam(ibin),'(a,i3.3,a)')'Dp(',ibin,')'
          write(numb(ibin),'(a,i3.3,a)')'Np(',ibin,')'
          write(area(ibin),'(a,i3.3,a)')'Sp(',ibin,')'
          write(mass(ibin),'(a,i3.3,a)')'Mp(',ibin,')'
        enddo
        if (iwrite_aer_dist > 0)   &
        write(lun_drydist,200)   &
                              (diam(ibin), ibin =1, nbin_a),   &
                              (numb(ibin), ibin =1, nbin_a),   &
                              (area(ibin), ibin =1, nbin_a),   &
                              (mass(ibin), ibin =1, nbin_a)


! Dpfix grid (header)
!        write(lun_drydist,200)
!     &                        (diam(ibin), ibin =1, 25),
!     &                        (numb(ibin), ibin =1, 25),
!     &                        (area(ibin), ibin =1, 25),
!     &                        (mass(ibin), ibin =1, 25)

! Lagrangian grid (header)
        if (iwrite_aer_dist > 0)   &
        write(lun_wetdist,200)(diam(ibin), ibin =1, nbin_a),   &
                              (numb(ibin), ibin =1, nbin_a),   &
                              (area(ibin), ibin =1, nbin_a),   &
                              (mass(ibin), ibin =1, nbin_a)

      endif

200   format(' UTC(hr)    t(hr)    ', 10000(a7,6x))





! do rebinning to a fixed Dp grid
      do ifix = 1, 25

        num_fix(ifix) = 0.
        area_fix(ifix) = 0.
        mass_fix(ifix) = 0.

        do ibin = 1, nbin_a
          if(Dp_dry_a(ibin)*1.e4 .ge. Dpfix_lo(ifix) .and.   &
             Dp_dry_a(ibin)*1.e4 .lt. Dpfix_hi(ifix))then

            num_fix(ifix) = num_fix(ifix) + num_a(ibin)
            area_fix(ifix) = area_fix(ifix) + area_dry_a(ibin)
            mass_fix(ifix) = mass_fix(ifix) + mass_dry_a(ibin)

          endif
        enddo

            num_fix(ifix) = num_fix(ifix)/0.1  !  = dN/dLogDp
            area_fix(ifix)= area_fix(ifix)/0.1 !  = dS/dLogDp
            mass_fix(ifix)= mass_fix(ifix)/0.1 !  = dM/dLogDp

      enddo

! DpgN
      if (iwrite_aer_dist > 0)   &
      write(lun_drydist,201)   &
        timeofday,				   &  ! time of day in UTC
        time_dum, 				   &  ! hours since start of simulation
        (Dp_dry_a(ibin)*1.e4,ibin=1, nbin_a),	   &  ! Dp_dry_a is in cm. output in um
        (num_a(ibin), ibin=1, nbin_a),		   &  ! num_a is in #/cc. output in #/cc
        (area_dry_a(ibin)*1.e8,ibin=1, nbin_a),	   &  ! area_dry_a is in cm^2/cc. output in um^2/cc
        (mass_dry_a(ibin)*1.e12,ibin=1, nbin_a)	! mass_dry_a is in g/cc. output in ug/m^3


! Dpfix grid
!      write(lun_drydist,201)
!     &  timeofday,				! time of day in UTC
!     &  time_dum, 				! hours since start of simulation
!     &  (Dpmean_fix(ifix),ifix=1, 25),		! Dpmean_fix is in um
!     &  (num_fix(ifix), ifix=1, 25),		! num_a is in #/cc. output in #/cc
!     &  (area_fix(ifix)*1.e8,ifix=1, 25),	! area_dry_a is in cm^2/cc. output in um^2/cc
!     &  (mass_fix(ifix)*1.e12,ifix=1, 25)	! mass_dry_a is in g/cc. output in ug/m^3



      if (iwrite_aer_dist > 0)   &
      write(lun_wetdist,201)   &
        timeofday,				   &  ! time of day in UTC
        time_dum, 				   &  ! hours since start of simulation
        (Dp_wet_a(ibin)*1.e4,ibin=1, nbin_a),	   &  ! Dp_wet_a is in cm. output in um
        (num_a(ibin), ibin=1, nbin_a),		   &  ! num_a is in #/cc. output in #/cc
        (area_wet_a(ibin)*1.e8,ibin=1, nbin_a),	   &  ! area_wet_a is in cm^2/cc. output in um^2/cc
        (mass_wet_a(ibin)*1.e12,ibin=1, nbin_a)	! mass_wet_a is in g/cc. output in ug/m^3

201   format(f7.3,2x,f8.3,10000(2x,e11.5))









! calculate dlogDpdry and dlogDpwet
      if(nbin_a .le. 2)then

        do ibin = 1, nbin_a
          dlogDpdry(ibin) = -1.0
          dlogDpwet(ibin) = -1.0
        enddo

      else

        dlogDpdry(1) = log10(Dp_dry_a(2)/Dp_dry_a(1))
        dlogDpwet(1) = log10(Dp_wet_a(2)/Dp_wet_a(1))

        dlogDpdry(nbin_a) = log10(Dp_dry_a(nbin_a)/Dp_dry_a(nbin_a-1))
        dlogDpwet(nbin_a) = log10(Dp_wet_a(nbin_a)/Dp_wet_a(nbin_a-1))

        do ibin = 2, nbin_a-1
          dlogDpdry(ibin)=log10( (Dp_dry_a(ibin+1) + Dp_dry_a(ibin))/   &
                                 (Dp_dry_a(ibin-1) + Dp_dry_a(ibin)) )

          dlogDpwet(ibin)=log10( (Dp_wet_a(ibin+1) + Dp_wet_a(ibin))/   &
                                 (Dp_wet_a(ibin-1) + Dp_wet_a(ibin)) )
        enddo

      endif





! output species distribution
!      if(1 == 1)goto 102

      if(idum .eq. 0) then
         if (iwrite_aer_species > 0)   &
         write(lun_species,100)
      end if

!101   format(f7.3,2x,f8.3,2x,i4,2x,a6,1x,i4,2x,i4,3x,f6.3,72(2x,e12.6))
101   format(f7.3,2x,f8.3,2x,i4,2x,a6,1x,i4,2x,i4,3x,f6.3,72(2x,e28.20))

      do ibin = 1, nbin_a

        if(mc(jc_h,ibin) .gt. 0.0)then
          pH(ibin) = -log10(mc(jc_h,ibin))
        else
          pH(ibin) = 0.0
        endif

        if (iwrite_aer_species > 0)   &
        write(lun_species,101)   &
        timeofday,				   &  ! time of day in UTC
        time_dum, 				   &  ! hours since start of simulation
        ibin,					   &  ! integer
        phasestate(jaerosolstate(ibin)),	   &  ! 1=solid, 2=liquid, 3=mixed
        iter_mesa(ibin),			   &  ! number of mesa iterations
        isteps_ASTEM,				   &  ! number of astem steps
        pH(ibin),                                  &  ! pH
        Dp_dry_a(ibin)*1.e4,			   &  ! microns (same as Dp_dry_a)
        Dp_wet_a(ibin)*1.e4,			   &  ! microns
        num_a(ibin),				   &  ! #/cc(air)
        gas(ih2so4_g)*conv4,			   &  ! nmol/m^3
        gas(ihno3_g)*conv4,			   &  ! nmol/m^3
        gas(ihcl_g)*conv4, 			   &  ! nmol/m^3
        gas(inh3_g)*conv4,			   &  ! nmol/m^3
        gas(imsa_g)*conv4,			   &  ! nmol/m^3
        gas(iaro1_g)*conv4,			   &  ! nmol/m^3
        gas(iaro2_g)*conv4,			   &  ! nmol/m^3
        gas(ialk1_g)*conv4,			   &  ! nmol/m^3
        gas(iole1_g)*conv4,			   &  ! nmol/m^3
        gas(iapi1_g)*conv4,			   &  ! nmol/m^3
        gas(iapi2_g)*conv4,			   &  ! nmol/m^3
        gas(ilim1_g)*conv4,			   &  ! nmol/m^3
        gas(ilim2_g)*conv4,			   &  ! nmol/m^3
        aer(iso4_a,jtotal,ibin)*conv4,		   &  ! nmol/m^3
        aer(ino3_a,jtotal,ibin)*conv4,		   &  ! nmol/m^3
        aer(icl_a, jtotal,ibin)*conv4,		   &  ! nmol/m^3
        aer(imsa_a,jtotal,ibin)*conv4,		   &  ! nmol/m^3
        aer(inh4_a,jtotal,ibin)*conv4,		   &  ! nmol/m^3
        aer(ina_a, jtotal,ibin)*conv4,		   &  ! nmol/m^3
        aer(ica_a, jtotal,ibin)*conv4,		   &  ! nmol/m^3
        aer(ioin_a,jtotal,ibin)*conv4,		   &  ! ng/m^3
        aer(ibc_a,jtotal,ibin)*conv4,		   &  ! ng/m^3
        aer(ioc_a,jtotal,ibin)*conv4,		   &  ! ng/m^3
        aer(iaro1_a,jtotal,ibin)*conv4,		   &  ! nmol/m^3
        aer(iaro2_a,jtotal,ibin)*conv4,		   &  ! nmol/m^3
        aer(ialk1_a,jtotal,ibin)*conv4,		   &  ! nmol/m^3
        aer(iole1_a,jtotal,ibin)*conv4,		   &  ! nmol/m^3
        aer(iapi1_a,jtotal,ibin)*conv4,		   &  ! nmol/m^3
        aer(iapi2_a,jtotal,ibin)*conv4,		   &  ! nmol/m^3
        aer(ilim1_a,jtotal,ibin)*conv4,		   &  ! nmol/m^3
        aer(ilim2_a,jtotal,ibin)*conv4,		   &  ! nmol/m^3
        water_a(ibin)*55.555*1.e9*conv4,	   &  ! nmol/m^3(air)
        aer(iso4_a,jliquid,ibin)*conv4,		   &  ! nmol/m^3
        aer(ino3_a,jliquid,ibin)*conv4,		   &  ! nmol/m^3
        aer(icl_a, jliquid,ibin)*conv4,		   &  ! nmol/m^3
        aer(imsa_a,jliquid,ibin)*conv4,		   &  ! nmol/m^3
        mc(jc_h,ibin)*1.e9*water_a(ibin)*conv4,	   &  ! nmol/m^3
        aer(inh4_a,jliquid,ibin)*conv4,		   &  ! nmol/m^3
        aer(ina_a, jliquid,ibin)*conv4,		   &  ! nmol/m^3
        aer(ica_a, jliquid,ibin)*conv4,		   &  ! nmol/m^3
        (electrolyte(je,jsolid,ibin)*conv4, je=1,nsalt),    &  ! 15 nmol/m^3
        electrolyte(jcaso4,jsolid,ibin)*conv4,     &  ! nmol/m^3
        electrolyte(jcaco3,jsolid,ibin)*conv4,     &  ! nmol/m^3
        area_dry_a(ibin)*1.e8,			   &  ! area_dry_a is in cm^2/cc. output in um^2/cc
        mass_dry_a(ibin)*1.e12,			   &  ! mass_dry_a is in g/cc. output in ug/m^3
        area_wet_a(ibin)*1.e8,			   &  ! area_wet_a is in cm^2/cc. output in um^2/cc
        mass_wet_a(ibin)*1.e12,			   &  ! mass_wet_a is in g/cc. output in ug/m^3
        num_a(ibin)/dlogDpdry(ibin),		   &  ! dN/dlogDpdry (1/cc)
        num_a(ibin)/dlogDpwet(ibin),		   &  ! dN/dlogDpwet (1/cc)
        area_dry_a(ibin)*1.e8/dlogDpdry(ibin),	   &  ! dS/dlogDpdry in um^2/cc
        mass_dry_a(ibin)*1.e12/dlogDpdry(ibin),	   &  ! dM/dlogDpdry in ug/m^3
        area_wet_a(ibin)*1.e8/dlogDpwet(ibin),	   &  ! dS/dlogDpwet in um^2/cc
        mass_wet_a(ibin)*1.e12/dlogDpwet(ibin),	   &  ! dM/dlogDpwet in ug/m^3
        dlogDpdry(ibin),   &
        dlogDpwet(ibin)

     end do


        if (iwrite_aer_species > 0)   &
        write(lun_species,*)


100   format(   &
       ' UTC(hr)    t(hr)  ibin  jstate iMESA iASTEEM  pH     ',   &
       ' Dpdry(um)     Dpwet(um)     N(#/cc)       H2SO4(g)     ',   &
       ' HNO3(g)       HCl(g)        NH3(g)        MSA(g)       ',   &
       ' ARO1(g)       ARO2(g)       ALK1(g)       OLE1(g)      ',   &
       ' API1(g)       API2(g)       LIM1(g)       LIM2(g)      ',   &
       ' SO4(t)        NO3(t)        Cl(t)         MSA(t)       ',   &
       ' NH4(t)        Na(t)         Ca(t)         OIN(t)       ',   &
       ' BC(t)         OC(t)         ARO1(t)       ARO2(t)      ',   &
       ' ALK1(t)       OLE1(t)       API1(t)       API2(t)      ',   &
       ' LIM1(t)       LIM2(t)       H2O(l)        SO4=(l)      ',   &
       ' NO3-(l)       Cl-(l)        MSA-(l)       H+(l)        ',   &
       ' NH4+(l)       Na+(l)        Ca++(l)       AmSO4(s)     ',   &
       ' Lvcite(s)     NH4HSO4(s)    NH4MSA(s)     NH4NO3(s)    ',   &
       ' NH4Cl(s)      Na2SO4(s)     Na3HSO4(s)    NaHSO4(s)    ',   &
       ' NaMSA(s)      NaNO3(s)      NaCl(s)       CaNO3(s)     ',   &
       ' CaCl2(s)      CaMSA2(s)     CaSO4(s)      CaCO3(s)     ',   &
       ' Sdry          Mdry          Swet          Mwet       ',   &
       ' dN/dlogDpdry  dN/dlogDpwet  dS/dlogDpdry  dM/dlogDpdry ',   &
       ' dS/dlogDpwet  dM/dlogDpwet  dlogDpdry     dlogDpwet')



! 102     continue














! output composition for each bin
      do ibin = 1, nbin_a
        if(jaerosolstate(ibin) .eq. no_aerosol)cycle

        if (iwrite_aer_bin > 0)   &
        call lun_aer_open(ibin)

        if(idum .eq. 0) then
           if (iwrite_aer_bin > 0)   &
           write(lun_aer(ibin),300)
        end if

        if(mc(jc_h,ibin) .gt. 0.0)then
          pH(ibin) = -log10(mc(jc_h,ibin))
        else
          pH(ibin) = 0.0
        endif

        if (iwrite_aer_bin > 0)   &
        write(lun_aer(ibin),301)   &
        timeofday,				   &  ! time of day in UTC
        RH,      				   &  ! time_dum hours since start of simulation
        ibin,					   &  ! integer
        phasestate(jaerosolstate(ibin)),	   &  ! 1=solid, 2=liquid, 3=mixed
        iter_mesa(ibin),			   &  ! number of mesa iterations
        isteps_ASTEM,				   &  ! number of astem steps
        Dp_dry_a(ibin)*1.e4,			   &  ! microns
        num_a(ibin),				   &  ! #/cc(air)
        gas(ih2so4_g)*conv4,			   &  ! nmol/m^3
        gas(ihno3_g)*conv4,			   &  ! nmol/m^3
        gas(ihcl_g)*conv4, 			   &  ! nmol/m^3
        gas(inh3_g)*conv4,			   &  ! nmol/m^3
        gas(imsa_g)*conv4,			   &  ! nmol/m^3
        gas(iaro1_g)*conv4,			   &  ! nmol/m^3
        gas(iaro2_g)*conv4,			   &  ! nmol/m^3
        gas(ialk1_g)*conv4,			   &  ! nmol/m^3
        gas(iole1_g)*conv4,			   &  ! nmol/m^3
        gas(iapi1_g)*conv4,			   &  ! nmol/m^3
        gas(iapi2_g)*conv4,			   &  ! nmol/m^3
        gas(ilim1_g)*conv4,			   &  ! nmol/m^3
        gas(ilim2_g)*conv4,			   &  ! nmol/m^3
        aer(iso4_a,jliquid,ibin)*conv4,		   &  ! nmol/m^3
        aer(ino3_a,jliquid,ibin)*conv4,		   &  ! nmol/m^3
        aer(icl_a, jliquid,ibin)*conv4,		   &  ! nmol/m^3
        aer(imsa_a,jliquid,ibin)*conv4,		   &  ! nmol/m^3
        mc(jc_h,ibin)*1.e9*water_a(ibin)*conv4,	   &  ! nmol/m^3
        aer(inh4_a,jliquid,ibin)*conv4,		   &  ! nmol/m^3
        aer(ina_a, jliquid,ibin)*conv4,		   &  ! nmol/m^3
        aer(ica_a, jliquid,ibin)*conv4,		   &  ! nmol/m^3
        water_a(ibin)*55.555*1.e9*conv4,	   &  ! nmol/m^3(air)
        pH(ibin),				   &  ! pH
        (electrolyte(je,jsolid,ibin)*conv4, je=1,nsalt),    &  ! 15 umol/m^3
        electrolyte(jcaso4,jsolid,ibin)*conv4,     &  ! nmol/m^3
        electrolyte(jcaco3,jsolid,ibin)*conv4,     &  ! nmol/m^3
        aer(ioin_a,jsolid,ibin)*conv4,		   &  ! ng/m^3
        aer(ibc_a,jsolid,ibin)*conv4,		   &  ! ng/m^3
        aer(ioc_a,jsolid,ibin)*conv4,	   	   &  ! ng/m^3
        aer(iaro1_a,jtotal,ibin)*conv4,		   &  ! nmol/m^3
        aer(iaro2_a,jtotal,ibin)*conv4,		   &  ! nmol/m^3
        aer(ialk1_a,jtotal,ibin)*conv4,		   &  ! nmol/m^3
        aer(iole1_a,jtotal,ibin)*conv4,		   &  ! nmol/m^3
        aer(iapi1_a,jtotal,ibin)*conv4,		   &  ! nmol/m^3
        aer(iapi2_a,jtotal,ibin)*conv4,		   &  ! nmol/m^3
        aer(ilim1_a,jtotal,ibin)*conv4,		   &  ! nmol/m^3
        aer(ilim2_a,jtotal,ibin)*conv4		      ! nmol/m^3	! razx

        if (iwrite_aer_bin > 0)   &
        call lun_aer_close(ibin)

     enddo



300   format(   &
       ' UTC(hr)    RH(%)  ibin  jstate iMESA iASTEEM   Dp(um)',   &
       '       N(#/cc)      H2SO4(g)     HNO3(g)      HCl(g)',   &
       '       NH3(g)       MSA(g)       ARO1(g)      ARO2(g)',   &
       '      ALK1(g)      OLE1(g)      API1(g)      API2(g)',   &
       '      LIM1(g)      LIM2(g)      SO4=(l)      NO3-(l)',   &
       '      Cl-(l)        MSA-(l)      H+(l)        NH4+(l)',   &
       '      Na+(l)       Ca++(l)     Water        pH       AmSO4(s)',   &
       '      Lvcite(s)     NH4HSO4(s)    NH4MSA(s)     NH4NO3(s)',   &
       '     NH4Cl(s)      Na2SO4(s)    Na3HSO4(s)',   &
       '     NaHSO4(s)     NaMSA(s)      NaNO3(s)      NaCl(s)',   &
       '       CaNO3(s)      CaCl2(s)      CaMSA2(s)      CaSO4(s)',   &
       '     CaCO3(s)      OIN(s)        BC(s)         OC(o)',   &
       '         ARO1(o)       ARO2(o)       ALK1(o)       OLE1(o)',   &
       '       API1(o)       API2(o)       LIM1(o)       LIM2(o)')

!301   format(f7.3,2x,f8.3,2x,i4,2x,a6,1x,i4,2x,i4,2x,24(2x,e11.5),   &
!             2x,f6.3,28(3x,e11.5))

301   format(f7.3,2x,f8.3,2x,i4,2x,a6,1x,i4,2x,i4,2x,24(2x,e28.20),   &
             2x,f6.3,28(3x,e28.20))

      return
      end subroutine print_aer

    end module module_print_aer
