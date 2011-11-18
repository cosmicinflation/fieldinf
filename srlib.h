!import required potential and derivated from libsrmodels

!large field  
  use lfisr, only : lfi_norm_potential, lfi_norm_deriv_potential
  use lfisr, only : lfi_norm_deriv_second_potential

!mixer large field
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
