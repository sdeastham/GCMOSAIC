subroutine IntegrateChemistry(it,                          & !intent-ins
     jaerosolstate, dp_wet_a,                              & !intent-inouts
     cair_mol_m3, cair_mol_cc, aH2O_a, gam_ratio,iter_mesa ) !intent-outs

  use module_data_mosaic_kind,  only: r8
  use module_data_mosaic_aero,  only: nbin_a_max, nsalt
  use module_data_mosaic_main,  only: told_sec, tcur_sec, mgas, maer, mcld
  implicit none

  !Subroutine arguments
  integer,  intent(in)    :: it
  integer,  intent(inout), dimension(nbin_a_max) :: jaerosolstate

  integer,  intent(out),   dimension(nbin_a_max) :: iter_MESA

  real(r8), intent(inout), dimension(nbin_a_max) :: dp_wet_a

  real(r8), intent(out) :: cair_mol_m3,cair_mol_cc
  real(r8), intent(out),   dimension(nbin_a_max) :: aH2O_a,gam_ratio

  !Local Variables
  real(r8) :: t_in, t_out

  t_in = told_sec
  t_out= tcur_sec

  if(mgas.eq.1)then
     call GasChemistry(t_in, t_out)
  endif

  if(maer.eq.1)then
     call AerChemistry(it, t_out, t_in,                         & !intent-ins
          jaerosolstate, dp_wet_a,                              & !intent-inouts
          cair_mol_m3, cair_mol_cc, aH2O_a, gam_ratio,iter_mesa ) !intent-outs
     
  endif

  if(mcld.eq.1)then
     call CldChemistry(t_in, t_out)
  endif


  return
end subroutine IntegrateChemistry
