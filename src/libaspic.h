!import required potential and derivatives from libaspic

!large field  
  use lfisr, only : lfi_norm_potential, lfi_norm_deriv_potential
  use lfisr, only : lfi_norm_deriv_second_potential

!mixed large field
  use gmlfisr, only : gmlfi_norm_potential, gmlfi_norm_deriv_potential
  use gmlfisr, only : gmlfi_norm_deriv_second_potential

!radiative corrected massive
  use rcmisr, only : rcmi_norm_potential, rcmi_norm_deriv_potential
  use rcmisr, only : rcmi_norm_deriv_second_potential

!radiative corrected quadratic
  use rcqisr, only : rcqi_norm_potential, rcqi_norm_deriv_potential
  use rcqisr, only : rcqi_norm_deriv_second_potential

!natural inflation
  use nisr, only : ni_norm_potential, ni_norm_deriv_potential
  use nisr, only : ni_norm_deriv_second_potential

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
  use lmicommon, only : lmi_norm_potential, lmi_norm_deriv_potential
  use lmicommon, only : lmi_norm_deriv_second_potential

!R + R^2p inflation
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

!mssm inflation
  use mssmisr, only : mssmi_norm_potential, mssmi_norm_deriv_potential
  use mssmisr, only : mssmi_norm_deriv_second_potential
  
!generalized mssm inflation
  use gmssmisr, only : gmssmi_norm_potential, gmssmi_norm_deriv_potential
  use gmssmisr, only : gmssmi_norm_deriv_second_potential

!original renormalizable inflection point inflation
  use oripisr, only : oripi_norm_potential, oripi_norm_deriv_potential
  use oripisr, only : oripi_norm_deriv_second_potential

!brane susy breaking inflation
  use bsusybisr, only : bsusybi_norm_potential, bsusybi_norm_deriv_potential
  use bsusybisr, only : bsusybi_norm_deriv_second_potential

!tip inflation
  use tisr, only : ti_norm_potential, ti_norm_deriv_potential
  use tisr, only : ti_norm_deriv_second_potential

!pseudo natural inflation
  use psnisr, only : psni_norm_potential, psni_norm_deriv_potential
  use psnisr, only : psni_norm_deriv_second_potential

!non-canonical Kahler inflation
  use nckisr, only : ncki_norm_potential, ncki_norm_deriv_potential
  use nckisr, only : ncki_norm_deriv_second_potential

!valley hybrid inflation
  use vhisr, only : vhi_norm_potential, vhi_norm_deriv_potential
  use vhisr, only : vhi_norm_deriv_second_potential

!dynamical susy inflation
  use dsisr, only : dsi_norm_potential, dsi_norm_deriv_potential
  use dsisr, only : dsi_norm_deriv_second_potential

!arctan inflation
  use aisr, only : ai_norm_potential, ai_norm_deriv_potential
  use aisr, only : ai_norm_deriv_second_potential

!constant ns inflation a, b, c and d
  use cnaisr, only : cnai_norm_potential, cnai_norm_deriv_potential
  use cnaisr, only : cnai_norm_deriv_second_potential

  use cnbisr, only : cnbi_norm_potential, cnbi_norm_deriv_potential
  use cnbisr, only : cnbi_norm_deriv_second_potential

  use cncisr, only : cnci_norm_potential, cnci_norm_deriv_potential
  use cncisr, only : cnci_norm_deriv_second_potential

  use cndisr, only : cndi_norm_potential, cndi_norm_deriv_potential
  use cndisr, only : cndi_norm_deriv_second_potential

!constant spectrum inflation
  use csisr, only : csi_norm_potential, csi_norm_deriv_potential
  use csisr, only : csi_norm_deriv_second_potential

!orientifold inflation
  use oisr, only : oi_norm_potential, oi_norm_deriv_potential
  use oisr, only : oi_norm_deriv_second_potential

!supergravity brane inflation
  use sbisr, only : sbi_norm_potential, sbi_norm_deriv_potential
  use sbisr, only : sbi_norm_deriv_second_potential

!spontaneous symmetry breaking inflation
  use ssbicommon, only : ssbi_norm_potential, ssbi_norm_deriv_potential
  use ssbicommon, only : ssbi_norm_deriv_second_potential

!running mass inflation
  use rmicommon, only : rmi_norm_potential, rmi_norm_deriv_potential
  use rmicommon, only : rmi_norm_deriv_second_potential

!logarithmic potential inflation
  use lpicommon, only : lpi_norm_potential, lpi_norm_deriv_potential
  use lpicommon, only : lpi_norm_deriv_second_potential

!open string tachyonic inflation
  use ostisr, only : osti_norm_potential, osti_norm_deriv_potential
  use ostisr, only : osti_norm_deriv_second_potential

!Witten O Raifeartaigh Higgs inflation
  use wrisr, only : wri_norm_potential, wri_norm_deriv_potential
  use wrisr, only : wri_norm_deriv_second_potential

!Inverse monomial inflation
  use imisr, only : imi_norm_potential, imi_norm_deriv_potential
  use imisr, only : imi_norm_deriv_second_potential

!renormalizable inflection point inflation
  use ripisr, only : ripi_norm_potential, ripi_norm_deriv_potential
  use ripisr, only : ripi_norm_deriv_second_potential

!generalized renormalizable inflection point inflatoin
  use gripisr, only : gripi_norm_potential, gripi_norm_deriv_potential
  use gripisr, only : gripi_norm_deriv_second_potential

!brane inflation potential in the large x limit
  use bisr, only : bi_norm_potential, bi_norm_deriv_potential
  use bisr, only : bi_norm_deriv_second_potential

!kklt brane inflation potential (witout assuming the large x limit)
  use kkltisr, only : kklti_norm_potential, kklti_norm_deriv_potential
  use kkltisr, only : kklti_norm_deriv_second_potential

!N-formalism inflation
  use nficommon, only : nfi_norm_potential, nfi_norm_deriv_potential
  use nficommon, only : nfi_norm_deriv_second_potential

!Cublicly corrected Starobinksi inflation
  use ccsicommon, only : ccsi_norm_potential, ccsi_norm_deriv_potential
  use ccsicommon, only : ccsi_norm_deriv_second_potential

!Dual inflation
  use disr, only : di_norm_potential, di_norm_deriv_potential
  use disr, only : di_norm_deriv_second_potential

!Hybrid Natural inflation
  use hnisr, only : hni_norm_potential, hni_norm_deriv_potential
  use hnisr, only : hni_norm_deriv_second_potential

!Non-renormalizable corrected loop inflation
  use nclisr, only : ncli_norm_potential, ncli_norm_deriv_potential
  use nclisr, only : ncli_norm_deriv_second_potential

!Viatcheslav Fyodorovich Mukhanov inflation
  use vfmi, only : vfmi_norm_potential, vfmi_norm_deriv_potential
  use vfmi, only : vfmi_norm_deriv_second_potential
