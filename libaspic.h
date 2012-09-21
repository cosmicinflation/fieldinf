!import required potential and derivatives from libaspic

!large field  
  use lfisr, only : lfi_norm_potential, lfi_norm_deriv_potential
  use lfisr, only : lfi_norm_deriv_second_potential

!mixed large field
  use mlfisr, only : mlfi_norm_potential, mlfi_norm_deriv_potential
  use mlfisr, only : mlfi_norm_deriv_second_potential

!radiative corrected massive
  use rcmisr, only : rcmi_norm_potential, rcmi_norm_deriv_potential
  use rcmisr, only : rcmi_norm_deriv_second_potential

!radiative corrected quadratic
  use rcqisr, only : rcqi_norm_potential, rcqi_norm_deriv_potential
  use rcqisr, only : rcqi_norm_deriv_second_potential

!natural with a plus sign
  use pnisr, only : pni_norm_potential, pni_norm_deriv_potential
  use pnisr, only : pni_norm_deriv_second_potential

!natural with a minus sign
  use mnisr, only : mni_norm_potential, mni_norm_deriv_potential
  use mnisr, only : mni_norm_deriv_second_potential

!exponential susy
  use esisr, only : esi_norm_potential, esi_norm_deriv_potential
  use esisr, only : esi_norm_deriv_second_potential

!power law
  use plisr, only : pli_norm_potential, pli_norm_deriv_potential
  use plisr, only : pli_norm_deriv_second_potential

!khaler moduli first type
  use kmiisr, only : kmii_norm_potential, kmii_norm_deriv_potential
  use kmiisr, only : kmii_norm_deriv_second_potential

!horizon flow linear
  use hf1isr, only : hf1i_norm_potential, hf1i_norm_deriv_potential
  use hf1isr, only : hf1i_norm_deriv_second_potential

!small field
  use sfisr, only : sfi_norm_potential, sfi_norm_deriv_potential
  use sfisr, only : sfi_norm_deriv_second_potential

!global susy with loop inflation
  use lisr, only : li_norm_potential, li_norm_deriv_potential
  use lisr, only : li_norm_deriv_second_potential

!intermediate inflation
  use iisr, only : ii_norm_potential, ii_norm_deriv_potential
  use iisr, only : ii_norm_deriv_second_potential

!khaler moduli second type
  use kmiiisr, only : kmiii_norm_potential, kmiii_norm_deriv_potential
  use kmiiisr, only : kmiii_norm_deriv_second_potential

!coleman-weinberg inflation
  use cwisr, only : cwi_norm_potential, cwi_norm_deriv_potential
  use cwisr, only : cwi_norm_deriv_second_potential

!higgs inflation
  use hisr, only : hi_norm_potential, hi_norm_deriv_potential
  use hisr, only : hi_norm_deriv_second_potential

!twisted inflation
  use twisr, only : twi_norm_potential, twi_norm_deriv_potential
  use twisr, only : twi_norm_deriv_second_potential
  
!topological double well inflation
  use dwisr, only : dwi_norm_potential, dwi_norm_deriv_potential
  use dwisr, only : dwi_norm_deriv_second_potential

!mutated hilltop inflation
  use mhisr, only : mhi_norm_potential, mhi_norm_deriv_potential
  use mhisr, only : mhi_norm_deriv_second_potential

!logamediate inflation
  use lmi1sr, only : lmi1_norm_potential, lmi1_norm_deriv_potential
  use lmi1sr, only : lmi1_norm_deriv_second_potential

  use lmi2sr, only : lmi2_norm_potential, lmi2_norm_deriv_potential
  use lmi2sr, only : lmi2_norm_deriv_second_potential

!R + R^p inflation
  use rpi1sr, only : rpi1_norm_potential, rpi1_norm_deriv_potential
  use rpi1sr, only : rpi1_norm_deriv_second_potential

  use rpi2sr, only : rpi2_norm_potential, rpi2_norm_deriv_potential
  use rpi2sr, only : rpi2_norm_deriv_second_potential

!beta exponential inflation
  use beisr, only : bei_norm_potential, bei_norm_deriv_potential
  use beisr, only : bei_norm_deriv_second_potential

!radion assisted gauge inflation
  use rgisr, only : rgi_norm_potential, rgi_norm_deriv_potential
  use rgisr, only : rgi_norm_deriv_second_potential
