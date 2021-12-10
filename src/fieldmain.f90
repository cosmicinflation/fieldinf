!   This file is part of FieldInf
!
!   Copyright (C) 2005-2021 C. Ringeval
!   
!   FieldInf is free software: you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation, either version 3 of the License, or
!   (at your option) any later version.
!
!   FieldInf is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with FieldInf.  If not, see <https://www.gnu.org/licenses/>.





!example test program using fieldinf for the large field models to get
!everything

program fieldmain

  use fieldprec, only : kp
  use infbgmodel, only : infbgparam, set_infbg_param
  use infbg, only : infbgphys, infbgdata, print_infbgphys
  use infbg, only : set_infbg_ini, bg_field_evol
  use inftorad, only : inftoradcosmo, print_inftoradcosmo
  use inftorad, only : set_inftorad_cosmo
  use inftorad, only : hubble_splinexit, infhubblexit
  use infpert, only : set_infpert_ini
  use infpert, only : scalNum 
!power_spectrum_scal() returns an array (scalnum,scalnum) with all cross correlators
!between fields and zeta (for single field scalnum=2); the last entry
!is Pzeta
  use infpert, only : power_spectrum_scal
  use infpert, only : power_spectrum_tens

  implicit none

!all inflationary model parameters
  type(infbgparam) :: infParam

!background physical quantities
  type(infbgphys) :: bgIni, bgEnd

!to store the background trajectory
  type(infbgdata), pointer :: ptrBgData => null()

!from inflation to radiation, reheating quantities
  type(inftoradcosmo) :: infCosmo

!to get some quantities at Hubble exit for any mode (such as N*)
  type(infhubblexit) :: atHkExit



!true if the model parameter are ok
  logical :: setupDone

!the reheating parameter R (rescaled)
  real(kp) :: lnReheat

!perturbation mode in Mpc^-1
  real(kp) :: kmpc, lnkmpc, lnKmpcMin, lnKmpcMax
  real(kp) :: Nkstar

!power spectra
  real(kp) :: powerTens, powerZeta
  real(kp), dimension(scalNum,scalNum) :: powerScal

  integer :: i, nkmpc



  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!specifying the model, its parameters, and initial field value;
!checkout infbgmodel.f90 or man libaspic
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!initialization
  infParam%consts = 0.

!model name
  infParam%name = 'largef'

!large field has 2 params, V=M^4 phi^p  
!M
  infparam%consts(1) = 1e-4
!p
  infparam%consts(2) = 2

!initial field values (in Planck unit). If zero, triggers an automatic
!guess from slow-roll integration using libaspic
  infParam%matters = 0._kp

!update parameters + consistency checks
  setupDone = set_infbg_param(infParam)
  if (.not.setupDone) stop 'model parameters incorrect'

!set and return all initial physical quantities in bgIni

  bgIni = set_infbg_ini(infParam)

  write(*,*)'Initial Conditions'
  call print_infbgphys(bgIni,'bgIni=')
  write(*,*)

  write(*,*)'INITIAL CONDITIONS SET!'




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!background integration
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!integrate the background till the end of inflation

!simplest usage for getting the background only
!  bgEnd = bg_field_evol(bgIni)

!simplest usage to also get the perturbations. Trajectory stored in
!  ptrBgData
  bgEnd = bg_field_evol(bgIni,ptrStart=ptrBgData)

  write(*,*)'End of Inflation'
  call print_infbgphys(bgEnd,'bgEnd=')
  write(*,*)

  write(*,*)'BACKGROUND INTEGRATED'




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Embedding inflation into FLRW cosmological model: from inflation to
!radiation from the reheating parameter
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!R, the reheating parameter (zero means radiation-like reheating)
  lnReheat = 0._kp

  infCosmo = set_inftorad_cosmo(infParam,bgIni,bgEnd,lnReheat)

  write(*,*)'Reheating settings'
  call print_inftoradcosmo(infCosmo,'infCosmo= ')
  write(*,*)
  
  if (.not.associated(ptrBgData)) stop 




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Linear perturbations
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!prepare integrating the perturbations
  call set_infpert_ini(infCosmo,ptrBgData)


!ln(k) varies from min to max (k in Mpc^-1)
  lnKmpcMin = -18._kp  
  lnKmpcMax = 0._kp

!number of modes in between
  nkmpc = 100

  do i=1,nkmpc-1
     lnkmpc =  lnkmpcMin + real(i-1)*(lnkmpcMax - lnkmpcMin)/real(nkmpc - 1)
     kmpc = exp(lnkmpc)

!get quantities at Hubble exit (optional)
     atHkExit = hubble_splinexit(infCosmo,kmpc)
     NkStar = atHkExit%bfold

     powerScal = power_spectrum_scal(infCosmo,kmpc)
     powerZeta = powerScal(scalNum,scalNum)

     powerTens = power_spectrum_tens(infCosmo,kmpc)

     write(*,*)'kmpc= Nkstar= Pzeta= Ph=',kmpc,Nkstar,powerZeta,powerTens

  enddo

  write(*,*)'PERTURBATIONS INTEGRATED.'
  
end program fieldmain
