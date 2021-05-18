!***********************SUBROUTINE****************************
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! 太陽位置を求めるサブルーチン
! PHI: 計算対象地点の緯度 ( 度 ),
! LON: 計算対象地点の経度 ( 度 )
! TM: 標準時 ( 時 ) ，
! LONS: 標準時の地点の経度 ( 度 ) ( 日本では 135.0)
! SINH: 太陽高度のサイン , COSH: 太陽高度のコサイン
! SINA: 太陽方位角のサイン， COSA: 太陽方位角のコサイン
! ET: 均時差 ( 度 ), T: 時角 ( 度 )
! SINDLT: 赤緯のサイン , COSDLT: 赤緯のコサイン
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

SUBROUTINE SUNLHA(PHI,LON,TM,LONS,SINDLT,COSDLT, ET,SINH,COSH,SINA,COSA)
    implicit none
    integer :: TM
    real:: T, PHIRAD, TRAD
    real:: PHI, LON, LONS, SINDLT, COSDLT
    real:: ET, SINH, COSH, SINA, COSA
    real, PARAMETER:: RAD=3.141592*2./360.
  
    LONS=135.0
    T = 15. * (TM - 12.) + (LON - LONS) + ET
    PHIRAD = PHI * RAD
    TRAD = T * RAD
  
  SINH = SIN(PHIRAD) * SINDLT+COS(PHIRAD) * COSDLT * COS(TRAD)
  COSH = SQRT(1.0 - SINH**2.0)
  SINA = COSDLT*SIN(TRAD)/COSH
  COSA = (SINH*SIN(PHIRAD)-SINDLT) / (COSH*COS(PHIRAD))
  
  END SUBROUTINE SUNLHA
  
  !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  !  赤坂の式により，赤緯と均時差を求めるサブルーチン
  !  YEAR: 西暦年 , NDAY: 年間通日 day number, start from 1.1
  !  SINDLT: 赤緯のサイン 
  !  COSDLT: 赤緯のコサイン
  !  ET: 均時差 ( 度 )
  !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  
  SUBROUTINE SUNLD(YEAR,NDAY,SINDLT,COSDLT,ET)
      implicit none
  real:: NDAY, SINDLT, COSDLT, ET
  real:: M, D0, N, EPS, V, VEPS, VE2
  integer :: YEAR
  real, PARAMETER ::RAD=3.141592*2./360.
  real, PARAMETER :: DLT0=-23.4393*RAD
  
  N = YEAR - 1968.
  
  D0 = 3.71 + 0.2596*N - INT((N+3.)/4.)
  
  M = 0.9856 * (NDAY -D0)
  EPS = 12.3901 + 0.0172 * (N + M/360.)
  V = M + 1.914 * SIN(M*RAD) + 0.02 * SIN(2.*M*RAD)
  VEPS=(V+EPS)*RAD
  VE2 = 2. * VEPS
  
  ET = (M - V) - ATAN(0.043*SIN(VE2)/(1.-0.043*COS(VE2))) /RAD
  
  SINDLT = COS(VEPS) * SIN(DLT0)
  COSDLT = SQRT(ABS(1. - SINDLT**2.))
  
  END SUBROUTINE SUNLD