MODULE MODULE_GHG_CHEM

CONTAINS

  subroutine ghg_chem_driver(dt,config_flags,            &
         t_phy,chem,rho_phy,ghg_oh,                       &
         ids,ide, jds,jde, kds,kde,                                        &
         ims,ime, jms,jme, kms,kme,                                        &
         its,ite, jts,jte, kts,kte                                         )
  USE module_configure
  USE module_state_description
  IMPLICIT NONE
   TYPE(grid_config_rec_type),  INTENT(IN   )    :: config_flags

   REAL, INTENT(IN   ) ::    dt

   INTEGER,      INTENT(IN   ) ::                                   &
                                  ids,ide, jds,jde, kds,kde,               &
                                  ims,ime, jms,jme, kms,kme,               &
                                  its,ite, jts,jte, kts,kte
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme, num_chem ),                 &
         INTENT(INOUT ) ::                                   chem
   REAL,  DIMENSION( ims:ime , kms:kme , jms:jme ),                        &
          INTENT(IN   ) ::                     t_phy,rho_phy
   REAL,  DIMENSION( ims:ime , kms:kme, jms:jme ),           &
          INTENT(IN   ) ::   ghg_oh
  integer :: nmx,cmx,tmx,i,j,k,imx,jmx,lmx
  real*8  :: tmp_oh, koh_ch4_cnst, koh_c13_cnst, &
          koh_co_cnst, koh_ch4, koh_c13, koh_co
  real*8, DIMENSION (1,1,1) :: tmp, oh_in, airden
  real*8, DIMENSION (1,1,1,9) :: tc1, tc2  !ch4 and c13
  real*8, DIMENSION (1,1,1,6) :: tc4   ! c2h6

!( 36)  CH4 + OH ->  CH3O2 + H2O + LCH4                                      rate = 2.45E-12*exp(  -1775./t)   ( 79)
!(336)  CTOT + OH ->  OH + LCTOT                                             rate = 2.44E-12*exp(  -1775./t)   (379)

! C2H6+OH=C2H5O2{+H2O}           : ARR2(8.7e-12, 1070., TEMP) ;
! CH4+OH=CH3O2+H2O               : ARR2(2.45e-12, 1775.0, TEMP) ;
! CO+HO=HO2+CO2                  : k57(TEMP,C_M) ;  from racm

  imx=1
  jmx=1
  lmx=1
  nmx=9
  cmx=5
  tmx=6

  do j=jts,jte
    do i=its,ite
      do k=kts,kte-1

          ! All source tagged tracers go through same chemical loss
          tc1(1,1,1,1)=chem(i,k,j,p_ch4_ant)*1.d-6
          tc1(1,1,1,2)=chem(i,k,j,p_ch4_bio)*1.d-6
          tc1(1,1,1,3)=chem(i,k,j,p_ch4_bbu)*1.d-6
          tc1(1,1,1,4)=chem(i,k,j,p_ch4_tst)*1.d-6
          tc1(1,1,1,5)=chem(i,k,j,p_ch4_ong)*1.d-6
          tc1(1,1,1,6)=chem(i,k,j,p_ch4_agr)*1.d-6
          tc1(1,1,1,7)=chem(i,k,j,p_ch4_wst)*1.d-6
          tc1(1,1,1,8)=chem(i,k,j,p_ch4_wet)*1.d-6
          tc1(1,1,1,9)=chem(i,k,j,p_ch4_bck)*1.d-6

          tc2(1,1,1,1)=chem(i,k,j,p_c13_ant)*1.d-6
          tc2(1,1,1,2)=chem(i,k,j,p_c13_bio)*1.d-6
          tc2(1,1,1,3)=chem(i,k,j,p_c13_bbu)*1.d-6
          tc2(1,1,1,4)=chem(i,k,j,p_c13_tst)*1.d-6
          tc2(1,1,1,5)=chem(i,k,j,p_c13_ong)*1.d-6
          tc2(1,1,1,6)=chem(i,k,j,p_c13_agr)*1.d-6
          tc2(1,1,1,7)=chem(i,k,j,p_c13_wst)*1.d-6
          tc2(1,1,1,8)=chem(i,k,j,p_c13_wet)*1.d-6
          tc2(1,1,1,9)=chem(i,k,j,p_c13_bck)*1.d-6

          tc4(1,1,1,1)=chem(i,k,j,p_eth_ant)*1.d-6
          tc4(1,1,1,2)=chem(i,k,j,p_eth_bck)*1.d-6
          tc4(1,1,1,3)=chem(i,k,j,p_eth_bbu)*1.d-6
          tc4(1,1,1,4)=chem(i,k,j,p_eth_tst)*1.d-6
          tc4(1,1,1,5)=chem(i,k,j,p_eth_bio)*1.d-6
          tc4(1,1,1,6)=chem(i,k,j,p_eth_ong)*1.d-6

          airden(1,1,1)=rho_phy(i,k,j)
          tmp(1,1,1)=t_phy(i,k,j)
          oh_in(1,1,1)=ghg_oh(i,k,j)*1.d-6  ! mol/mol

          call chmdrv_sink( imx,jmx,lmx,nmx,cmx,tmx,&
               dt, tmp, airden, oh_in, tc1, tc2, tc4)

          !Currently we do not include fire emissions for all the tracers
          !tracer_bbu will be same as tracer_bck but without chemical loss
          
          chem(i,k,j,p_ch4_ant)=tc1(1,1,1,1)*1.e6  ! ppm
          chem(i,k,j,p_ch4_bio)=tc1(1,1,1,2)*1.e6
          !chem(i,k,j,p_ch4_bbu)=tc1(1,1,1,3)*1.e6
          chem(i,k,j,p_ch4_tst)=tc1(1,1,1,4)*1.e6
          chem(i,k,j,p_ch4_ong)=tc1(1,1,1,5)*1.e6
          chem(i,k,j,p_ch4_agr)=tc1(1,1,1,6)*1.e6
          chem(i,k,j,p_ch4_wst)=tc1(1,1,1,7)*1.e6
          chem(i,k,j,p_ch4_wet)=tc1(1,1,1,8)*1.e6
          chem(i,k,j,p_ch4_bck)=tc1(1,1,1,9)*1.e6
  
          chem(i,k,j,p_c13_ant)=tc2(1,1,1,1)*1.e6  ! ppm
          chem(i,k,j,p_c13_bio)=tc2(1,1,1,2)*1.e6
          !chem(i,k,j,p_c13_bbu)=tc2(1,1,1,3)*1.e6
          chem(i,k,j,p_c13_tst)=tc2(1,1,1,4)*1.e6
          chem(i,k,j,p_c13_ong)=tc2(1,1,1,5)*1.e6
          chem(i,k,j,p_c13_agr)=tc2(1,1,1,6)*1.e6
          chem(i,k,j,p_c13_wst)=tc2(1,1,1,7)*1.e6
          chem(i,k,j,p_c13_wet)=tc2(1,1,1,8)*1.e6
          chem(i,k,j,p_c13_bck)=tc2(1,1,1,9)*1.e6

          chem(i,k,j,p_eth_ant)=tc4(1,1,1,1)*1.e6  ! ppm
          chem(i,k,j,p_eth_bck)=tc4(1,1,1,2)*1.e6
          !chem(i,k,j,p_eth_bbu)=tc4(1,1,1,3)*1.e6
          chem(i,k,j,p_eth_tst)=tc4(1,1,1,4)*1.e6
          chem(i,k,j,p_eth_bio)=tc4(1,1,1,5)*1.e6
          chem(i,k,j,p_eth_ong)=tc4(1,1,1,6)*1.e6

       enddo
     enddo
   enddo
end subroutine ghg_chem_driver

SUBROUTINE chmdrv_sink( imx,jmx,lmx,nmx,cmx,tmx, &
               dt1,tmp,airden,oh,tc1,tc2,tc4)

! ****************************************************************************
! **                                                                        **
! **  Chemistry subroutine to calculate methane loss due to OH              **
! **                                                                        **
! ****************************************************************************
  
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: nmx,imx,jmx,lmx,cmx,tmx
  integer :: ndt1
  real, intent(in) :: dt1
  REAL*8, DIMENSION(imx,jmx,lmx), INTENT(IN) :: tmp, airden
  REAL*8, DIMENSION(imx,jmx,lmx), INTENT(IN) :: oh
  REAL*8, INTENT(INOUT) :: tc1(imx,jmx,lmx,nmx),tc2(imx,jmx,lmx,nmx)
  REAL*8, INTENT(INOUT) :: tc4(imx,jmx,lmx,tmx)

  ! executable statements
  ndt1=ifix(dt1)
  if(ndt1.le.0)stop

     CALL chem_ch4_loss(imx,jmx,lmx,nmx,ndt1, tmp, airden, oh, tc1, tc2)
     !CALL chem_co_loss(imx,jmx,lmx,cmx,ndt1, tmp, airden, oh, tc3)
     CALL chem_eth_loss(imx,jmx,lmx,tmx,ndt1, tmp, airden, oh, tc4)
  
END SUBROUTINE chmdrv_sink

!=============================================================================
SUBROUTINE chem_ch4_loss( imx,jmx,lmx,nmx,ndt1, tmp, airden, oh, tc1, tc2)

! ****************************************************************************
! *                                                                          *
! *  This is CH4 chemistry subroutine.                                       *
! *                                                                          *
! *  R1:   CH4 + OH ->  CH3O2 + H2O + LCH4                                   *
! *        rate = 2.45E-12*exp(  -1775./t)                                   *
! *                                                                          *
! *  R2:   C13 + OH ->                                                       *
! *        rate = 2.43866E-12*exp(  -1775./t)                                *
! **************************************************************************** 

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: nmx, ndt1,imx,jmx,lmx
  REAL*8, DIMENSION(imx,jmx,lmx), INTENT(IN) :: tmp, airden
  REAL*8, DIMENSION(imx,jmx,lmx), INTENT(IN) :: oh
  REAL*8, DIMENSION(imx,jmx,lmx,nmx), INTENT(INOUT) :: tc1,tc2

  ! -- Factor to convert AIRDEN from kgair/m3 to molecules/cm3: 
  REAL*8, PARAMETER :: f = 1000.0 / 28.97 * 6.022D23 * 1.0D-6
  INTEGER :: i, j, l
  REAL(KIND=8) :: tk, rk1, rk2, rk3,  xoh,  xx, o2
  REAL*8, DIMENSION(nmx):: chem1_0, chem1_oh, chem2_0, chem2_oh
  
  ! executable statements
 
  rk1 = 0.0d0
  rk2 = 0.0d0
 
  DO l = 1,lmx
!CMIC$ doall autoscope
     DO j = 1,jmx
        DO i = 1,imx
           
           tk = tmp(i,j,l)
           o2 = airden(i,j,l) * f * 0.21
           chem1_0 = tc1(i,j,l,:)
           chem2_0 = tc2(i,j,l,:)
! ****************************************************************************
! *  (1) CH4 + OH:  RK1                                                      *
! ****************************************************************************
           IF (oh(i,j,l) > 0.0) THEN
              !rate constant in cm^3 molec^-1 s^-1 for 2nd order
              !oh is in mol/mol, convert to molecules/cm3
               rk1 = 2.45D-12*EXP(-1775.0/tk) * oh(i,j,l)*airden(i,j,l)*f
           END IF
           
           chem1_oh = chem1_0 * EXP( -rk1  * REAL(ndt1) )
           tc1(i,j,l,:) = MAX(chem1_oh, 1.0D-32)

! ****************************************************************************
! *  (2) C13 + OH:  RK2                                                      *
! ****************************************************************************
           IF (oh(i,j,l) > 0.0) THEN
           ! JianHe: mean fractionation based on Saueressig et al. (2001), 
           ! (medsink: (1.0039+1.0054)/2), hisink: 1.0039)
                 rk2 = 2.45D-12/1.00465*EXP(-1775.0/tk) * oh(i,j,l)*airden(i,j,l)*f
           END IF

           chem2_oh = chem2_0 * EXP( -rk2  * REAL(ndt1) )
           tc2(i,j,l,:) = MAX(chem2_oh, 1.0D-32)

           
        END DO
     END DO
  END DO
  
END SUBROUTINE chem_ch4_loss
 
!=============================================================================
SUBROUTINE chem_co_loss( imx,jmx,lmx,cmx,ndt1, tmp, airden, oh, tc3)

! ****************************************************************************
! *                                                                          *
! *  This is CO chemistry subroutine.                                        *
! *                 
! *  R1: ! CO+HO=HO2+CO2                : k57(TEMP,C_M) ;  from racm         *
! ****************************************************************************

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: cmx, ndt1,imx,jmx,lmx
  REAL*8, DIMENSION(imx,jmx,lmx), INTENT(IN) :: tmp, airden
  REAL*8, DIMENSION(imx,jmx,lmx), INTENT(IN) :: oh
  REAL*8, DIMENSION(imx,jmx,lmx,cmx), INTENT(INOUT) :: tc3

  ! -- Factor to convert AIRDEN from kgair/m3 to molecules/cm3: 
  REAL*8, PARAMETER :: f = 1000.0 / 28.97 * 6.022D23 * 1.0D-6

  INTEGER :: i, j, l
  REAL(KIND=8) :: tk, rk,rk1, rk2, rk3,  xoh,  xx
  REAL*8, DIMENSION(cmx):: chem_0, chem_oh

  REAL*8 :: temp       ! temperature [K]
  REAL*8 :: c_m        ! air concentration [molecules/cm3]
  REAL*8 :: k0_300Kn   ! low pressure limit at 300 K
  REAL*8 :: nn         ! exponent for low pressure limit
  REAL*8 :: kinf_300Kn ! high pressure limit at 300 K
  REAL*8 :: mn         ! exponent for high pressure limit
  REAL*8 :: zt_help, k0_T, kinf_T, k_ratio
  REAL*8 :: k57troe, k57cact, k57

  rk = 0.0d0

  DO l = 1,lmx
     DO j = 1,jmx
        DO i = 1,imx

          k0_300Kn = 5.9D-33
          nn = 1.4
          kinf_300Kn = 1.1D-12
          mn = -1.3

          temp = tmp(i,j,l)
          c_m = airden(i,j,l) * f
          chem_0 = tc3(i,j,l,:)

          zt_help = 300./temp
          k0_T    = k0_300Kn   * zt_help**(nn) * c_m ! k_0   at current T
          kinf_T  = kinf_300Kn * zt_help**(mn)       ! k_inf at current T
          k_ratio = k0_T/kinf_T
          k57troe   = k0_T/(1.+k_ratio)*0.6**(1./(1.+LOG10(k_ratio)**2))

          k0_300Kn = 1.5D-13
          nn = -0.6
          kinf_300Kn = 2.9D9
          mn = -6.1

          k0_T    = k0_300Kn   * zt_help**(nn)! k_0   at current T
          kinf_T  = kinf_300Kn * zt_help**(mn) / c_m ! k_inf at current T
          k_ratio = k0_T/kinf_T
          k57cact = k0_T/(1.+k_ratio)*0.6**(1./(1.+LOG10(k_ratio)**2))

          k57 = k57troe + k57cact   ! total rate constant

          rk = k57*oh(i,j,l)*c_m

          chem_oh = chem_0 * EXP( -rk  * REAL(ndt1) )
          tc3(i,j,l,:) = MAX(chem_oh, 1.0D-32)

        END DO
      END DO
    END DO

END SUBROUTINE chem_co_loss


!=============================================================================
SUBROUTINE chem_eth_loss( imx,jmx,lmx,tmx,ndt1, tmp, airden, oh, tc4)

! ****************************************************************************
! *                                                                          *
! *  This is C2H6 chemistry subroutine.                                      *
! *                 
! *  R1: ! C2H6+OH=C2H5O2{+H2O}           : ARR2(8.7e-12, 1070., TEMP)       *
! ****************************************************************************

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: tmx, ndt1,imx,jmx,lmx
  REAL*8, DIMENSION(imx,jmx,lmx), INTENT(IN) :: tmp, airden
  REAL*8, DIMENSION(imx,jmx,lmx), INTENT(IN) :: oh
  REAL*8, DIMENSION(imx,jmx,lmx,tmx), INTENT(INOUT) :: tc4

  ! -- Factor to convert AIRDEN from kgair/m3 to molecules/cm3: 
  REAL*8, PARAMETER :: f = 1000.0 / 28.97 * 6.022D23 * 1.0D-6

  INTEGER :: i, j, l
  REAL(KIND=8) :: tk, rk,rk1, rk2, rk3,  xoh,  xx
  REAL*8, DIMENSION(tmx):: chem_0, chem_oh

  rk = 0.0d0

  DO l = 1,lmx
     DO j = 1,jmx
        DO i = 1,imx

           tk = tmp(i,j,l)
           chem_0 = tc4(i,j,l,:)

           IF (oh(i,j,l) > 0.0) THEN
              !rate constant in cm^3 molec^-1 s^-1 for 2nd order
              !oh is in mol/mol, convert to molecules/cm3

               rk = 8.7D-12*EXP(-1070.0/tk) * oh(i,j,l)*airden(i,j,l)*f
           END IF

           chem_oh = chem_0 * EXP( -rk  * REAL(ndt1) )
           tc4(i,j,l,:) = MAX(chem_oh, 1.0D-32)

        END DO
      END DO
    END DO

END SUBROUTINE chem_eth_loss

!=============================================================================

END MODULE MODULE_GHG_CHEM
