!------------------------------------------------------
!Program transform weather into HASP style
!              2020.9.3  Z.J
!------------------------------------------------------
program HASP
    implicit none
    character(len=80) :: data     !AMeDaS data
    integer :: temp(24), humi(24), irra(24), irrad(24)
    integer :: wind(24), windi(24), windv(24)
    integer :: month(12)
    integer :: TM, mon, day, week, i, days_of_year
    integer :: YEAR, first_day, iyr, num
    real    :: NDAY, SINDLT, COSDLT, ET
    real    :: PHI, LON, LONS, ALT
    real    :: SINH, COSH, SINA, COSA, sinD
    real    :: SNH(24)           !Solar Elevation Angle for 24 hours'  [-]
    real    :: SNH1, SNH2, SNH3  !Solar Elevation Angle for 3 hours'   [-]
    real    :: IG1, IG2, IG3     !Global Irradiance for 3 hours'       [0.01MJ/m2]
    real    :: temp1, humi1
    real    :: IG
    real    :: nora_r
    real    :: YY, MM, DD, HH    !for Matsumoto sunlocation function
    integer :: IB, ID ,nora      !ID_R：水平面天空日射量 IB_R：法線面直達日射量 
    character(len=80) :: fin     !! len=80  very important!
    character(len=80) :: fout    !! len=80  very important!  
    character(len=3)  :: stnum   !station number
    character(len=4)  :: cyr     !year number
    character(len=80)  :: WD(7)  !weather data in HASP style 

    integer :: PHI_D
    integer :: LON_D
    real :: PHI_N(842), PHI_M                           !計算対象地点の緯度 ( degree,minitue )
    real :: LON_N(842), LON_M                           !計算対象地点の経度 ( degree,minitue )  
    real :: ALT_N(842)                                  !ALTITUDE OF SITE  [m]
    character(len=80) :: St_info                   !file name of station information
    character(len=6) ::p, l, a, b, c, d, e, f

!************read longitude, latitude, and altitude***************
    St_info = '../StationInformation/StnInfo1981.dat' 
    open(3,file=St_info)
    do  num = 1,842
        read(3,*)a,b,c,d,e,f,p,l,ALT_N(num)
        read(p(1:2),*)PHI_D
        read(p(3:5),*)PHI_M  
        read(l(1:3),*)LON_D
        read(l(4:6),*)LON_M
        PHI_N(num) = PHI_D + PHI_M/10/60
        LON_N(num) = LON_D + LON_M/10/60
        !write(*,*)num, PHI(num), LON(num), ALT(num)
    end do
    close(3)


    do YEAR=1986,1986                                   !1st do : number of a year
        do num=565,565                                      !2nd do : number of a station
            LON=LON_N(num)
            PHI=PHI_N(num)
            ALT=ALT_N(num)
            write(stnum,'(i3.3)')num
            write(cyr,'(i4)')YEAR
            
            
    fin = '../ASCiidata/'//cyr//'/'//trim(stnum)//'0'// cyr //'.dat'
    fout = '../HASPdata/' //cyr//'/'//trim(stnum)//'0'// cyr //'.has' !trim(adjustl(stnum))      
    open(9,file=fin)
    open(2,file=fout)
    
    call day_of_month (year, days_of_year, month, first_day)
    NDAY=0                                              !年間通日 , for Erbs model
    week=first_day + 1                                  !day of a week at the first in a year (day = day + 1 in HASP)
    call day_of_week (week,week)                        !judge the day of a week
    
    do mon=1,12                                         !3rd do : number of a month
        do day=1,month(mon)                             !4th do : number of a day in a month
            NDAY=NDAY+1                                 !年間通日

!read AMeDaS data of 24 hours'
            do TM= 1,24                                 !5th do : number of a hour in a day
              read(9,'(A70)')data                       !AMeDaS data
              read(data(9:12),*)temp(TM)                !気温
              read(data(17:20),*)humi(TM)               !絶対湿度
              read(data(25:28),*)irra(TM)               !全天日射量
              read(data(33:36),*)irrad(TM)              !大気放射量
              read(data(43:44),*)windi(TM)              !風向
              read(data(49:52),*)windv(TM)              !風速
!--------------------------------1.  Akasaka sunlocation function------------------------------------------------------
              call SUNLD(YEAR,NDAY,SINDLT,COSDLT,ET)
              call SUNLHA(PHI,LON,TM,LONS,SINDLT,COSDLT, ET,SINH,COSH,SINA,COSA)
!--------------------------------2.  Matsumoto sunlocation function----------------------------------------------------
              !YY = real(YEAR)                  !Matsumoto function
              !MM = real(mon)                   !Matsumoto function
              !DD = real(day)                   !Matsumoto function
              !HH = real(TM)                    !Matsumoto function
              !call Sunlocation_Matsumoto(YY,MM,DD,HH,days_of_year,PHI,LON,SINH,COSH,Et,sinD)
!----------------------------------------------------------------------------------------------------------------------
              SNH(TM)=SINH                              !Solar Elevation Angle for 24 hours'  [-]
            end do                                      !end 5th do 

        do TM= 1,24                                     !6th do : number of a hour in a day
!calculate normal direct irradiation and horizontal diffuse irradiation---------直散分離計算
!assignment of Global Irradiance for 3 hours'   [0.01MJ/m2]
            if ( TM == 1 )then
              IG1=-999.
              IG2=irra(TM)
              IG3=irra(TM+1)
            else if ( TM == 24 )then
              IG1=irra(TM-1)
              IG2=irra(TM)
              IG3=-999.
            else
              IG1=irra(TM-1)
              IG2=irra(TM)
              IG3=irra(TM+1)
            end if
!assignment of Solar Elevation Angle for 3 hours'   [-]
            if ( TM == 1 )then
                SNH1=-999.
                SNH2=SNH(TM)
                SNH3=SNH(TM+1)
            else if ( TM == 24 )then
                SNH1=SNH(TM-1)
                SNH2=SNH(TM)
                SNH3=-999.
            else
                SNH1=SNH(TM-1)
                SNH2=SNH(TM)
                SNH3=SNH(TM+1)
            end if

            humi1=humi(TM)*0.0001                        !absolute humidity  from [0.1g/kg] to [kg/kg] 
            temp1=temp(TM)*0.1                           !temperature        from [0.1℃] to [℃]
!----------------------------------------1.  Perez model---------------------------------------------------------
            call Perez_model(IG1,IG2,IG3,SNH1,SNH2,SNH3,humi1,temp1,ALT,IB,ID) !calculate by Perez model
!----------------------------------------2.  Erbs model----------------------------------------------------------
            !call Erbs_model(IG2,SINH,IB,ID)  !calculate by Erbs model
!----------------------------------------------------------------------------------------------------------------

!calculate nocturnal radiation----------------------夜間放射
            nora_r=((temp(TM)*0.1+273.15)**4*5.67/10**8)*0.86 - irrad(TM)*2.3889    !地面放射　－　大気放射
            nora=nint(nora_r)

!write weater data in HASP style            
            write(WD(1)(3*TM-2:3*TM),'(i3)')(temp(TM)+500) !temperature * 24 hours
            write(WD(2)(3*TM-2:3*TM),'(i3)')humi(TM)       !absolute humidity * 24 hours
            write(WD(3)(3*TM-2:3*TM),'(i3)')IB             !normal direct irradiation * 24 hours
            write(WD(4)(3*TM-2:3*TM),'(i3)')ID             !horizontal diffuse irradiation * 24 hours
            write(WD(5)(3*TM-2:3*TM),'(i3)')nora           !nocturnal radiation * 24 hours
            write(WD(6)(3*TM-2:3*TM),'(i3)')windi(TM)      !wind direction * 24 hours
            write(WD(7)(3*TM-2:3*TM),'(i3)')windv(TM)      !wind velocity * 24 hours
        end do                                             !end 6th do : end 24 hour

!output the 7 kinds of data in HASP
        do i = 1,7
            WD(i)(73:74)=cyr(3:4)                          !add year
            write(WD(i)(75:76),'(i2)')mon                  !add month
            write(WD(i)(77:78),'(i2)')day                  !add day
            write(WD(i)(79:79),'(i1)')week                 !add day of a week
            write(WD(i)(80:80),'(i1)')i                    !add number of data
            write(2,'(a)')WD(i)                            !output the 7 kinds of data in HASP style
        end do

        week=week+1                                      !calculate the next day of a week
        call day_of_week (week,week)                     !judge the day of a week
        end do                                           !end 4th do :day
    end do                                               !end 3rd do :month
  close(2)
  close(9)
  end do                                                 !end 2nd do :station number
end do                                                   !end 1st do :year

end program HASP
!*********************************END PROGRAM****************************************




!***********************SUBROUTINE****************************
subroutine day_of_month (year, days_of_year, month, first_day)   
    implicit none
    integer, intent(in)  :: year                  !year number
    integer, intent(out)    :: days_of_year          !the number of days in a year
    integer :: n
    integer, intent(out) :: month(1:12)
    integer, intent(out) :: first_day         !the first day of a year
    integer              ::D
    real                 ::d1
    
!judging if the year is leap year.    
    if (mod(year,100)==0) then                          
        if(mod(year,400)==0) then
           days_of_year=366
        else
         days_of_year=365
        end if
    else if (mod(year,4)==0) then
         days_of_year=366
    else
         days_of_year=365
    end if

!calculate the number of days in a month
     do n=1, 12, 1
        select case (n)
        case(1,3,5,7,8,10,12)
            month(n)=31
        case(2)
            month(n)=28
        case(4,6,9,11)
            month(n)=30
        end select
     
        if (days_of_year==366) then
            month(2)=29
        else
        end if    
    end do

!calculate the first day of a year
    d1=(year-1)/4                                       ! int(x) 整数化 (小数点以下切り捨て) mod(x,y) x を y で割ったあまり 
    D=aint(d1)*(365*3+366)+mod((year-1),4)*365          ! day betwween 0001/01/01 to year/01/01 
    first_day=mod(D,7)                                  ! which day of year/01/01 
        
end subroutine day_of_month
!-----------------------------------------

subroutine day_of_week (week,week_real)   
    implicit none
    integer, intent(in)  :: week
    integer, intent(out) :: week_real
    if (week==8) then
        week_real = 1
    else
        week_real = week
    end if
end subroutine day_of_week
!----------------------------------------


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
!-----------------Erbs model-------------------------------
  SUBROUTINE Erbs_model(IG2,SINH,IB,ID)
    implicit none
    real, intent(in) :: IG2       !Global Irradiance                           [0.01MJ/m2h] 
    real, intent(in) :: SINH     !Solar Elevation Angle sinh                  [-] 
    real :: IG 
    real :: Kt                   !clearness index Kt                          [-]
    real :: IB_R                 !normal direct irradiation                             [MJ/m2]
    real :: ID_R                 !DIFFUSE IRRADIANCE                          [MJ/m2]
    integer :: IB                !normal direct irradiation                             [kcal/m2h] 
    integer :: ID                !DIFFUSE IRRADIANCE                          [kcal/m2h] 
    real, parameter :: IO=4.921  !Solar constant                              [MJ/m2]

    IG=IG2*0.01       ![MJ/m2]
    Kt=IG/IO/SINH    !clearness index Kt [-]
    if ( Kt <= 0.22 ) then
        ID_R=IG*(1.0 - 0.09*Kt)
      else if ( Kt <= 0.8 )then
        ID_R=IG*(0.9511 - 0.1604*Kt + 4.388*Kt**2 - 16.638*Kt**3 + 12.336*Kt**4)
      else if (Kt > 0.8)then 
        ID_R=IG*0.165        
    end if
    
    if ( SINH < 0.0694 )then
        IB_R=0
      else if( SINH >= 0.0694 )then    !sin(3.98)=0.0694 
        IB_R=(IG-ID_R)/SINH
        ID_R=ID_R*238.89               ![kcal/m2h]
        IB_R=IB_R*238.89               ![kcal/m2h]
    end if

    IB=nint(IB_R) !integer part of normal direct irradiation [kcal/m2h]
    ID=nint(ID_R) !integer part of Diffuse Irrandiance [kcal/m2h]

end SUBROUTINE Erbs_model


!*****************************PEREZ MODEL*******************************
subroutine Perez_model(IG1,IG2,IG3,SNH1,SNH2,SNH3,humi,temp,ALT,IB,ID)
    implicit none
!---------------ARGUMENTS-------------------------
    real, intent(in) :: ALT        !altitude of site                            [m]
    real, intent(in) :: humi       !absolute humidity                           [kg/kg] 
    real, intent(in) :: temp       !temperature                                 [℃]
    real, intent(in) :: SNH1, SNH2, SNH3  !Solar Elevation Angle for 3 hours'   [-]
    real, intent(in) :: IG1, IG2, IG3     !Global Irradiance for 3 hours'       [0.01MJ/m2]
    real :: IO                     !Solar constant                              [MJ/m2]
    real :: TD                     !Dew point temperature                       [℃]
    real :: RTOD                   !1 rad = 57.29577951°                        [°]
    real :: IG(3)                  !Global Irradiance                           [MJ/m2]
    real :: SNH(3)                 !Solar Elevation Angle sinh                  [-] 
    real :: Z(3)                   !Zenith Angle                                [RADIANS]
    real :: Zd(3)                  !Zenith Angle                                [°]
    real :: EP                     !vapor pressure-Dalton's law                 [hPa]
    real :: EE                     !saturated vapor pressure-Tetens equation    [hPa]
    real :: E                      !vapor pressure-dew point temperature        [Pa]
    real :: EW                     !saturated vapor pressure-Wagner  equation   [Pa]
    real :: RH                     !relative humidity                           [%]
    real :: y                      !argument for dew point temperature          [-]
    real :: T                      !temperature                                 [K]
    real :: W                      !precipitable water                          [cm]
    real :: AM(3)                  !relative air mass                           [–]
    real :: CALT                   !coefficient                                 [-]
    real :: KT1(3)                 !clearness index Kt'                         [-]
    real :: DKT1                   !stability index ∆Kt'                        [-] 
    real :: A, B, C                !coefficient                                 [-]
    real :: Knc, BMAX              !coefficient                                 [-]
!--------------RETURNS--------------------------
    integer, intent(out) :: IB     !normal direct irradiation                   [kcal/m2h] 
    integer, intent(out) :: ID     !DIFFUSE IRRADIANCE                          [kcal/m2h] 
    real :: IB_R                   !normal direct irradiation                   [MJ/m2]
    real :: ID_R                   !DIFFUSE IRRADIANCE                          [MJ/m2]
    real :: CM(6,6,7,5)            !CORRECTION MATRIX                           [-]
    real :: Kt(3)                  !clearness index Kt                          [-]
    real :: KTBIN(6)               !clearness index Kt'                         [-]
    real :: ZBIN(6)                !solar zenith angle Z                        [°]
    real :: WBIN(4)                !precipitable water W                        [cm]
    real :: DKTBIN(6)              !stability index ∆Kt'                        [-]
    
    real, PARAMETER :: RAD=3.141592*2./360.      
    integer :: I, J, K, L, I1, J1, K1, L1
    DATA RTOD / 57.295779513082316 /                                           ![1 rad]
    DATA KTBIN /0.24, 0.4, 0.56, 0.7, 0.8, 1.0/                 
    DATA ZBIN / 25.0, 40.0, 55.0, 70.0, 80.0, 90.0/
    DATA DKTBIN / 0.015, 0.035, 0.07, 0.15, 0.3, 1.0/
    DATA WBIN / 1.0, 2.0, 3.0,1.E+38 /

DATA ((CM(1,1,I,J),J=1,5),I=1,7) / &
0.385230, 0.385230, 0.385230, 0.462880, 0.317440, &
0.338390, 0.338390, 0.221270, 0.316730, 0.503650, &
0.235680, 0.235680, 0.241280, 0.157830, 0.269440, &
0.830130, 0.830130, 0.171970, 0.841070, 0.457370, &
0.548010, 0.548010, 0.478000, 0.966880, 1.036370, &
0.548010, 0.548010, 1.000000, 3.012370, 1.976540, &
0.582690, 0.582690, 0.229720, 0.892710, 0.569950  /
DATA ((CM(1,2,I,J),J=1,5),I=1,7) / &
0.131280, 0.131280, 0.385460, 0.511070, 0.127940, &
0.223710, 0.223710, 0.193560, 0.304560, 0.193940, &
0.229970, 0.229970, 0.275020, 0.312730, 0.244610, &
0.090100, 0.184580, 0.260500, 0.687480, 0.579440, &
0.131530, 0.131530, 0.370190, 1.380350, 1.052270, &
1.116250, 1.116250, 0.928030, 3.525490, 2.316920, &
0.090100, 0.237000, 0.300040, 0.812470, 0.664970 /
DATA ((CM(1,3,I,J),J=1,5),I=1,7) / &
0.587510, 0.130000, 0.400000, 0.537210, 0.832490, &
0.306210, 0.129830, 0.204460, 0.500000, 0.681640, &
0.224020, 0.260620, 0.334080, 0.501040, 0.350470, &
0.421540, 0.753970, 0.750660, 3.706840, 0.983790, &
0.706680, 0.373530, 1.245670, 0.864860, 1.992630, &
4.864400, 0.117390, 0.265180, 0.359180, 3.310820, &
0.392080, 0.493290, 0.651560, 1.932780, 0.898730 /
DATA ((CM(1,4,I,J),J=1,5),I=1,7) / &
0.126970, 0.126970, 0.126970, 0.126970, 0.126970, &
0.810820, 0.810820, 0.810820, 0.810820, 0.810820, &
3.241680, 2.500000, 2.291440, 2.291440, 2.291440, &
4.000000, 3.000000, 2.000000, 0.975430, 1.965570, &
12.494170, 12.494170, 8.000000, 5.083520, 8.792390, &
21.744240, 21.744240, 21.744240, 21.744240, 21.744240, &
3.241680, 12.494170, 1.620760, 1.375250, 2.331620 /
DATA ((CM(1,5,I,J),J=1,5),I=1,7) / &
0.126970, 0.126970, 0.126970, 0.126970, 0.126970, &
0.810820, 0.810820, 0.810820, 0.810820, 0.810820, &
3.241680, 2.500000, 2.291440, 2.291440, 2.291440, &
4.000000, 3.000000, 2.000000, 0.975430, 1.965570, &
12.494170, 12.494170, 8.000000, 5.083520, 8.792390, &
21.744240, 21.744240, 21.744240, 21.744240, 21.744240, &
3.241680, 12.494170, 1.620760, 1.375250, 2.331620 /
DATA ((CM(1,6,I,J),J=1,5),I=1,7) / &
0.126970, 0.126970, 0.126970, 0.126970, 0.126970, &
0.810820, 0.810820, 0.810820, 0.810820, 0.810820, &
3.241680, 2.500000, 2.291440, 2.291440, 2.291440, &
4.000000, 3.000000, 2.000000, 0.975430, 1.965570, &
12.494170, 12.494170, 8.000000, 5.083520, 8.792390, &
21.744240, 21.744240, 21.744240, 21.744240, 21.744240, &
3.241680, 12.494170, 1.620760, 1.375250, 2.331620 /
DATA ((CM(2,1,I,J),J=1,5),I=1,7) / &
0.337440, 0.337440, 0.969110, 1.097190, 1.116080, &
0.337440, 0.337440, 0.969110, 1.116030, 0.623900, &
0.337440, 0.337440, 1.530590, 1.024420, 0.908480, &
0.584040, 0.584040, 0.847250, 0.914940, 1.289300, &
0.337440, 0.337440, 0.310240, 1.435020, 1.852830, &
0.337440, 0.337440, 1.015010, 1.097190, 2.117230, &
0.337440, 0.337440, 0.969110, 1.145730, 1.476400 /
DATA ((CM(2,2,I,J),J=1,5),I=1,7) / &
0.300000, 0.300000, 0.700000, 1.100000, 0.796940, &
0.219870, 0.219870, 0.526530, 0.809610, 0.649300, &
0.386650, 0.386650, 0.119320, 0.576120, 0.685460, &
0.746730, 0.399830, 0.470970, 0.986530, 0.785370, &
0.575420, 0.936700, 1.649200, 1.495840, 1.335590, &
1.319670, 4.002570, 1.276390, 2.644550, 2.518670, &
0.665190, 0.678910, 1.012360, 1.199940, 0.986580 /
DATA ((CM(2,3,I,J),J=1,5),I=1,7) / &
0.378870, 0.974060, 0.500000, 0.491880, 0.665290, &
0.105210, 0.263470, 0.407040, 0.553460, 0.582590, &
0.312900, 0.345240, 1.144180, 0.854790, 0.612280, &
0.119070, 0.365120, 0.560520, 0.793720, 0.802600, &
0.781610, 0.837390, 1.270420, 1.537980, 1.292950, &
1.152290, 1.152290, 1.492080, 1.245370, 2.177100, &
0.424660, 0.529550, 0.966910, 1.033460, 0.958730 /
DATA ((CM(2,4,I,J),J=1,5),I=1,7) / &
0.310590, 0.714410, 0.252450, 0.500000, 0.607600, &
0.975190, 0.363420, 0.500000, 0.400000, 0.502800, &
0.175580, 0.196250, 0.476360, 1.072470, 0.490510, &
0.719280, 0.698620, 0.657770, 1.190840, 0.681110, &
0.426240, 1.464840, 0.678550, 1.157730, 0.978430, &
2.501120, 1.789130, 1.387090, 2.394180, 2.394180, &
0.491640, 0.677610, 0.685610, 1.082400, 0.735410 /
DATA ((CM(2,5,I,J),J=1,5),I=1,7) / &
0.597000, 0.500000, 0.300000, 0.310050, 0.413510, &
0.314790, 0.336310, 0.400000, 0.400000, 0.442460, &
0.166510, 0.460440, 0.552570, 1.000000, 0.461610, &
0.401020, 0.559110, 0.403630, 1.016710, 0.671490, &
0.400360, 0.750830, 0.842640, 1.802600, 1.023830, &
3.315300, 1.510380, 2.443650, 1.638820, 2.133990, &
0.530790, 0.745850, 0.693050, 1.458040, 0.804500 /
DATA ((CM(2,6,I,J),J=1,5),I=1,7) / &
0.597000, 0.500000, 0.300000, 0.310050, 0.800920, &
0.314790, 0.336310, 0.400000, 0.400000, 0.237040, &
0.166510, 0.460440, 0.552570, 1.000000, 0.581990, &
0.401020, 0.559110, 0.403630, 1.016710, 0.898570, &
0.400360, 0.750830, 0.842640, 1.802600, 3.400390, &
3.315300, 1.510380, 2.443650, 1.638820, 2.508780, &
0.204340, 1.157740, 2.003080, 2.622080, 1.409380 /
DATA ((CM(3,1,I,J),J=1,5),I=1,7) / &
1.242210, 1.242210, 1.242210, 1.242210, 1.242210, &
0.056980, 0.056980, 0.656990, 0.656990, 0.925160, &
0.089090, 0.089090, 1.040430, 1.232480, 1.205300, &
1.053850, 1.053850, 1.399690, 1.084640, 1.233340, &
1.151540, 1.151540, 1.118290, 1.531640, 1.411840, &
1.494980, 1.494980, 1.700000, 1.800810, 1.671600, &
1.018450, 1.018450, 1.153600, 1.321890, 1.294670 /
DATA ((CM(3,2,I,J),J=1,5),I=1,7) / &
0.700000, 0.700000, 1.023460, 0.700000, 0.945830, &
0.886300, 0.886300, 1.333620, 0.800000, 1.066620, &
0.902180, 0.902180, 0.954330, 1.126690, 1.097310, &
1.095300, 1.075060, 1.176490, 1.139470, 1.096110, &
1.201660, 1.201660, 1.438200, 1.256280, 1.198060, &
1.525850, 1.525850, 1.869160, 1.985410, 1.911590, &
1.288220, 1.082810, 1.286370, 1.166170, 1.119330 /
DATA ((CM(3,3,I,J),J=1,5),I=1,7) / &
0.600000, 1.029910, 0.859890, 0.550000, 0.813600, &
0.604450, 1.029910, 0.859890, 0.656700, 0.928840, &
0.455850, 0.750580, 0.804930, 0.823000, 0.911000, &
0.526580, 0.932310, 0.908620, 0.983520, 0.988090, &
1.036110, 1.100690, 0.848380, 1.035270, 1.042380, &
1.048440, 1.652720, 0.900000, 2.350410, 1.082950, &
0.817410, 0.976160, 0.861300, 0.974780, 1.004580 /
DATA ((CM(3,4,I,J),J=1,5),I=1,7) / &
0.782110, 0.564280, 0.600000, 0.600000, 0.665740, &
0.894480, 0.680730, 0.541990, 0.800000, 0.669140, &
0.487460, 0.818950, 0.841830, 0.872540, 0.709040, &
0.709310, 0.872780, 0.908480, 0.953290, 0.844350, &
0.863920, 0.947770, 0.876220, 1.078750, 0.936910, &
1.280350, 0.866720, 0.769790, 1.078750, 0.975130, &
0.725420, 0.869970, 0.868810, 0.951190, 0.829220 /
DATA ((CM(3,5,I,J),J=1,5),I=1,7) / &
0.791750, 0.654040, 0.483170, 0.409000, 0.597180, &
0.566140, 0.948990, 0.971820, 0.653570, 0.718550, &
0.648710, 0.637730, 0.870510, 0.860600, 0.694300, &
0.637630, 0.767610, 0.925670, 0.990310, 0.847670, &
0.736380, 0.946060, 1.117590, 1.029340, 0.947020, &
1.180970, 0.850000, 1.050000, 0.950000, 0.888580, &
0.700560, 0.801440, 0.961970, 0.906140, 0.823880 /
DATA ((CM(3,6,I,J),J=1,5),I=1,7) / &
0.500000, 0.500000, 0.586770, 0.470550, 0.629790, &
0.500000, 0.500000, 1.056220, 1.260140, 0.658140, &
0.500000, 0.500000, 0.631830, 0.842620, 0.582780, &
0.554710, 0.734730, 0.985820, 0.915640, 0.898260, &
0.712510, 1.205990, 0.909510, 1.078260, 0.885610, &
1.899260, 1.559710, 1.000000, 1.150000, 1.120390, &
0.653880, 0.793120, 0.903320, 0.944070, 0.796130 /
DATA ((CM(4,1,I,J),J=1,5),I=1,7) / &
1.000000, 1.000000, 1.050000, 1.170380, 1.178090, &
0.960580, 0.960580, 1.059530, 1.179030, 1.131690, &
0.871470, 0.871470, 0.995860, 1.141910, 1.114600, &
1.201590, 1.201590, 0.993610, 1.109380, 1.126320, &
1.065010, 1.065010, 0.828660, 0.939970, 1.017930, &
1.065010, 1.065010, 0.623690, 1.119620, 1.132260, &
1.071570, 1.071570, 0.958070, 1.114130, 1.127110 /
DATA ((CM(4,2,I,J),J=1,5),I=1,7) / &
0.950000, 0.973390, 0.852520, 1.092200, 1.096590, &
0.804120, 0.913870, 0.980990, 1.094580, 1.042420, &
0.737540, 0.935970, 0.999940, 1.056490, 1.050060, &
1.032980, 1.034540, 0.968460, 1.032080, 1.015780, &
0.900000, 0.977210, 0.945960, 1.008840, 0.969960, &
0.600000, 0.750000, 0.750000, 0.844710, 0.899100, &
0.926800, 0.965030, 0.968520, 1.044910, 1.032310 /
DATA ((CM(4,3,I,J),J=1,5),I=1,7) / &
0.850000, 1.029710, 0.961100, 1.055670, 1.009700, &
0.818530, 0.960010, 0.996450, 1.081970, 1.036470, &
0.765380, 0.953500, 0.948260, 1.052110, 1.000140, &
0.775610, 0.909610, 0.927800, 0.987800, 0.952100, &
1.000990, 0.881880, 0.875950, 0.949100, 0.893690, &
0.902370, 0.875960, 0.807990, 0.942410, 0.917920, &
0.856580, 0.928270, 0.946820, 1.032260, 0.972990 /
DATA ((CM(4,4,I,J),J=1,5),I=1,7) / &
0.750000, 0.857930, 0.983800, 1.056540, 0.980240, &
0.750000, 0.987010, 1.013730, 1.133780, 1.038250, &
0.800000, 0.947380, 1.012380, 1.091270, 0.999840, &
0.800000, 0.914550, 0.908570, 0.999190, 0.915230, &
0.778540, 0.800590, 0.799070, 0.902180, 0.851560, &
0.680190, 0.317410, 0.507680, 0.388910, 0.646710, &
0.794920, 0.912780, 0.960830, 1.057110, 0.947950 /
DATA ((CM(4,5,I,J),J=1,5),I=1,7) / &
0.750000, 0.833890, 0.867530, 1.059890, 0.932840, &
0.979700, 0.971470, 0.995510, 1.068490, 1.030150, &
0.858850, 0.987920, 1.043220, 1.108700, 1.044900, &
0.802400, 0.955110, 0.911660, 1.045070, 0.944470, &
0.884890, 0.766210, 0.885390, 0.859070, 0.818190, &
0.615680, 0.700000, 0.850000, 0.624620, 0.669300, &
0.835570, 0.946150, 0.977090, 1.049350, 0.979970 /
DATA ((CM(4,6,I,J),J=1,5),I=1,7) / &
0.689220, 0.809600, 0.900000, 0.789500, 0.853990, &
0.854660, 0.852840, 0.938200, 0.923110, 0.955010, &
0.938600, 0.932980, 1.010390, 1.043950, 1.041640, &
0.843620, 0.981300, 0.951590, 0.946100, 0.966330, &
0.694740, 0.814690, 0.572650, 0.400000, 0.726830, &
0.211370, 0.671780, 0.416340, 0.297290, 0.498050, &
0.843540, 0.882330, 0.911760, 0.898420, 0.960210 /
DATA ((CM(5,1,I,J),J=1,5),I=1,7) / &
1.054880, 1.075210, 1.068460, 1.153370, 1.069220, &
1.000000, 1.062220, 1.013470, 1.088170, 1.046200, &
0.885090, 0.993530, 0.942590, 1.054990, 1.012740, &
0.920000, 0.950000, 0.978720, 1.020280, 0.984440, &
0.850000, 0.908500, 0.839940, 0.985570, 0.962180, &
0.800000, 0.800000, 0.810080, 0.950000, 0.961550, &
1.038590, 1.063200, 1.034440, 1.112780, 1.037800 /
DATA ((CM(5,2,I,J),J=1,5),I=1,7) / &
1.017610, 1.028360, 1.058960, 1.133180, 1.045620, &
0.920000, 0.998970, 1.033590, 1.089030, 1.022060, &
0.912370, 0.949930, 0.979770, 1.020420, 0.981770, &
0.847160, 0.935300, 0.930540, 0.955050, 0.946560, &
0.880260, 0.867110, 0.874130, 0.972650, 0.883420, &
0.627150, 0.627150, 0.700000, 0.774070, 0.845130, &
0.973700, 1.006240, 1.026190, 1.071960, 1.017240 /
DATA ((CM(5,3,I,J),J=1,5),I=1,7) / &
1.028710, 1.017570, 1.025900, 1.081790, 1.024240, &
0.924980, 0.985500, 1.014100, 1.092210, 0.999610, &
0.828570, 0.934920, 0.994950, 1.024590, 0.949710, &
0.900810, 0.901330, 0.928830, 0.979570, 0.913100, &
0.761030, 0.845150, 0.805360, 0.936790, 0.853460, &
0.626400, 0.546750, 0.730500, 0.850000, 0.689050, &
0.957630, 0.985480, 0.991790, 1.050220, 0.987900 /
DATA ((CM(5,4,I,J),J=1,5),I=1,7) / &
0.992730, 0.993880, 1.017150, 1.059120, 1.017450, &
0.975610, 0.987160, 1.026820, 1.075440, 1.007250, &
0.871090, 0.933190, 0.974690, 0.979840, 0.952730, &
0.828750, 0.868090, 0.834920, 0.905510, 0.871530, &
0.781540, 0.782470, 0.767910, 0.764140, 0.795890, &
0.743460, 0.693390, 0.514870, 0.630150, 0.715660, &
0.934760, 0.957870, 0.959640, 0.972510, 0.981640 /
DATA ((CM(5,5,I,J),J=1,5),I=1,7) / &
0.965840, 0.941240, 0.987100, 1.022540, 1.011160, &
0.988630, 0.994770, 0.976590, 0.950000, 1.034840, &
0.958200, 1.018080, 0.974480, 0.920000, 0.989870, &
0.811720, 0.869090, 0.812020, 0.850000, 0.821050, &
0.682030, 0.679480, 0.632450, 0.746580, 0.738550, &
0.668290, 0.445860, 0.500000, 0.678920, 0.696510, &
0.926940, 0.953350, 0.959050, 0.876210, 0.991490 /
DATA ((CM(5,6,I,J),J=1,5),I=1,7) / &
0.948940, 0.997760, 0.850000, 0.826520, 0.998470, &
1.017860, 0.970000, 0.850000, 0.700000, 0.988560, &
1.000000, 0.950000, 0.850000, 0.606240, 0.947260, &
1.000000, 0.746140, 0.751740, 0.598390, 0.725230, &
0.922210, 0.500000, 0.376800, 0.517110, 0.548630, &
0.500000, 0.450000, 0.429970, 0.404490, 0.539940, &
0.960430, 0.881630, 0.775640, 0.596350, 0.937680 /
DATA ((CM(6,1,I,J),J=1,5),I=1,7) / &
1.030000, 1.040000, 1.000000, 1.000000, 1.049510, &
1.050000, 0.990000, 0.990000, 0.950000, 0.996530, &
1.050000, 0.990000, 0.990000, 0.820000, 0.971940, &
1.050000, 0.790000, 0.880000, 0.820000, 0.951840, &
1.000000, 0.530000, 0.440000, 0.710000, 0.928730, &
0.540000, 0.470000, 0.500000, 0.550000, 0.773950, &
1.038270, 0.920180, 0.910930, 0.821140, 1.034560 /
DATA ((CM(6,2,I,J),J=1,5),I=1,7) / &
1.041020, 0.997520, 0.961600, 1.000000, 1.035780, &
0.948030, 0.980000, 0.900000, 0.950360, 0.977460, &
0.950000, 0.977250, 0.869270, 0.800000, 0.951680, &
0.951870, 0.850000, 0.748770, 0.700000, 0.883850, &
0.900000, 0.823190, 0.727450, 0.600000, 0.839870, &
0.850000, 0.805020, 0.692310, 0.500000, 0.788410, &
1.010090, 0.895270, 0.773030, 0.816280, 1.011680 /
DATA ((CM(6,3,I,J),J=1,5),I=1,7) / &
1.022450, 1.004600, 0.983650, 1.000000, 1.032940, &
0.943960, 0.999240, 0.983920, 0.905990, 0.978150, &
0.936240, 0.946480, 0.850000, 0.850000, 0.930320, &
0.816420, 0.885000, 0.644950, 0.817650, 0.865310, &
0.742960, 0.765690, 0.561520, 0.700000, 0.827140, &
0.643870, 0.596710, 0.474460, 0.600000, 0.651200, &
0.971740, 0.940560, 0.714880, 0.864380, 1.001650 /
DATA ((CM(6,4,I,J),J=1,5),I=1,7) / &
0.995260, 0.977010, 1.000000, 1.000000, 1.035250, &
0.939810, 0.975250, 0.939980, 0.950000, 0.982550, &
0.876870, 0.879440, 0.850000, 0.900000, 0.917810, &
0.873480, 0.873450, 0.751470, 0.850000, 0.863040, &
0.761470, 0.702360, 0.638770, 0.750000, 0.783120, &
0.734080, 0.650000, 0.600000, 0.650000, 0.715660, &
0.942160, 0.919100, 0.770340, 0.731170, 0.995180 /
DATA ((CM(6,5,I,J),J=1,5),I=1,7) / &
0.952560, 0.916780, 0.920000, 0.900000, 1.005880, &
0.928620, 0.994420, 0.900000, 0.900000, 0.983720, &
0.913070, 0.850000, 0.850000, 0.800000, 0.924280, &
0.868090, 0.807170, 0.823550, 0.600000, 0.844520, &
0.769570, 0.719870, 0.650000, 0.550000, 0.733500, &
0.580250, 0.650000, 0.600000, 0.500000, 0.628850, &
0.904770, 0.852650, 0.708370, 0.493730, 0.949030 /
DATA ((CM(6,6,I,J),J=1,5),I=1,7) / &
0.911970, 0.800000, 0.800000, 0.800000, 0.956320, &
0.912620, 0.682610, 0.750000, 0.700000, 0.950110, &
0.653450, 0.659330, 0.700000, 0.600000, 0.856110, &
0.648440, 0.600000, 0.641120, 0.500000, 0.695780, &
0.570000, 0.550000, 0.598800, 0.400000, 0.560150, &
0.475230, 0.500000, 0.518640, 0.339970, 0.520230, &
0.743440, 0.592190, 0.603060, 0.316930, 0.794390 /

IG(1)=IG1*0.01  !Global Irradiance at n-1 o'clock   [MJ/m2]
IG(2)=IG2*0.01  !Global Irradiance at n   o'clock   [MJ/m2]
IG(3)=IG3*0.01  !Global Irradiance at n+1 o'clock   [MJ/m2]
SNH(1)=SNH1     !Solar Elevation Angle at n-1 o'clock
SNH(2)=SNH2     !Solar Elevation Angle at n   o'clock
SNH(3)=SNH3     !Solar Elevation Angle at n+1 o'clock

IO=4.921        !Solar constant 
CALT = EXP(-0.0001184*ALT)

!calculate Dew-point temperature TD
EP = (1013.25-ALT/9)*(28.9398373*humi/(28.9398373*humi+18.))    
EE = 6.11*10.**(7.5*temp/(237.3+temp))      
RH = EP/EE*100.                            
T = temp+273.15   
EW = exp(-6096.9385*(T**(-1.))+21.2409642-2.711193*(0.01)*T &
& +1.673952*(10.**(-5.))*(T**2.)+2.433502*log(T))
E = RH/100.*EW
y = log(E/611.213)
if (y>=0.)then
  TD = 13.715*y+8.4262*0.1*y**2.+ 1.9048*0.01*y**3.+ 7.8158*0.001*y**4.
else if (y<0.)then
  TD = 13.7204*y+7.36631*0.1*y**2.+ 3.32136*0.01*y**3.+ 7.78591*0.0001*y**4.
end if
!-------------------------------------------------------------------------
J = 1 
K = 3

!missing data
if ((IG(1) == -999.0) .OR. ( SNH(1) == -999.0 )) then
  J = 2
  KT1(1) = -999.0
end if
if ((IG(3) == -999.0) .OR. ( SNH(3) == -999.0 )) then
  K = 2
  KT1(3) = -999.0
end if

!calculate clearness index Kt' and solar zenith angle Z
do I=J,K,1
    if ( SNH(I) < 0.0 ) then
      KT1(I) = -999.0
     else 
      Zd(I)= acos(SNH(I))* RTOD
      Kt(I) = IG(I) / (IO * MAX(0.052336,SNH(I)))
      AM(I) = AMIN1(15.25,1. / (SNH(I)+0.15 * (93.885 - Zd(I))**(-1.253)))
      KT1(I) = Kt(I) / ( 1.031 * EXP(-1.4 / (0.9 + 9.4 / (AM(I) * CALT))) + 0.1)
    end if
end do

Knc=0.866 - 0.122*AM(2) + 0.0121*AM(2)**2 - 0.000653*AM(2)**3 + 0.000014*AM(2)**4

!calculate coefficient A, B, C
if (Kt(2) <= 0.6) then
    A=0.512-1.56*Kt(2)+ 2.286*Kt(2)**2 - 2.222*Kt(2)**3
    B= 0.37 + 0.962*Kt(2)
    C=-0.28 + Kt(2)*0.932 - 2.048*Kt(2)**2
else if (Kt(2) > 0.6) then
    A=-5.743+Kt(2)*21.77 - 27.49*Kt(2)**2 + 11.56*kt(2)**3
    B= 41.4 - 118.5*kt(2) +66.05*Kt(2)**2 +31.9*kt(2)**3
    C=-47.01+Kt(2)*184.2 - 222*Kt(2)**2 + 73.81*kt(2)**3
end if

BMAX = Knc - (A + B * EXP(C * AM(2)))

!calculate stability index ∆Kt'
if (( KT1(1) == -999.0 ) .AND. ( KT1(3) == -999.0 )) then
    K = 7
else if ( KT1(1) == -999.0 ) then
    DKT1 = ABS(KT1(3) - KT1(2))
else if ( KT1(3) == -999.0 ) then
    DKT1 = ABS(KT1(2) - KT1(1))
else
    DKT1 = (ABS(KT1(2) - KT1(1)) + ABS(KT1(3) - KT1(2))) / 2.0
end if

!select index of stability index ∆Kt'
do K1=1,5,1
    if (  (DKT1 > DKTBIN(K1)) ) then   
        K=K1+1
    else if (  (DKT1 < DKTBIN(1)) ) then  
        K=1
    end if
end do
!select index of clearness index Kt'
DO I1=1,5,1
    if (  (KT1(2) > KTBIN(I1)) ) then
      I=I1+1
    else if (  (KT1(2) < KTBIN(1)) ) then  
      I=1
    end if
end do   

!select index of solar zenith angle Zd
DO J1=1,5,1
    if ( (Zd(2) > ZBIN(J1)) ) then
        J=J1+1
    else if (  (Zd(2) < ZBIN(1)) ) then 
        J=1
    end if
end do

!calculate precipitable water W
if ( TD == -999.0 ) then
    L = 5
  else
    W = EXP( -0.075 + 0.07 * TD )
end if

!select index of precipitable water W
DO L1=1,3,1
    if ((W > WBIN(L1)) ) then
        L=L1+1
    else if ( (W < WBIN(1)) ) then
        L=1
    end if
end do

IB_R = IO * BMAX * CM(I,J,K,L)              !normal direct irradiation [kcal/m2h]
ID_R = IG(2) - IB_R*MAX(0.052336,SNH(2))    !Diffuse Irrandiance [kcal/m2h]

!abnormal data
if ((IG(2) < 0.01) .or. (SNH(2) <= 0)) then 
    IB_R = 0.0
    ID_R = IG(2)
end if

if ( BMAX < 0.0 ) then
    IB_R = 0.0
    ID_R = IG(2)
end if

if (ID_R<0.) then
    ID_R = 0.
end if

if (ID_R>IG(2)) then
    ID_R = IG(2)
end if

ID_R=ID_R*238.85     !normal direct irradiation [kcal/m2h]
IB_R=IB_R*238.85     !Diffuse Irrandiance [kcal/m2h]

IB=nint(IB_R)    !integer part of normal direct irradiation [kcal/m2h]
ID=nint(ID_R)    !integer part of Diffuse Irrandiance [kcal/m2h]

end subroutine Perez_model



!*****************************************************************************
subroutine Sunlocation_Matsumoto(YY,MM,DD,HH,days_of_year,PHI,LON,SINH,COSH,Et,sinD)
    implicit none
    real, intent(in) :: YY, MM, DD, HH
    integer, intent(in) :: days_of_year
    real(8) :: JDu      !Julian day of UTC
    real :: Tu       !Julian century of UTC
    real :: T        !Julian century of TCG
    real :: asec     !average right ascension(平均赤経)
    real :: psi      !celestial longitude(視黄経) [deg.]
    real :: eps      !章動を含み黄道傾斜角         [deg.]
    real :: r        !動径 (地心距離) (返値は AU 単位)
    real, intent(out) :: Et       !均時差 (deg.)
    real :: hAngle   !時角 (deg.)
    real :: Decl     !(視) 赤緯 Decl(deg.)
    real :: PHI, LON  !Latitude緯度, Longitude経度
    real :: TM
    real, intent(out) :: SinH, sinD, CosH ! 太陽高度角の正弦，余弦
    real :: SinA, CosA ! 太陽方位角の正弦，余弦
    real :: CosD ! 視赤緯の正弦，余弦 
    real, PARAMETER   :: rpd = 3.1415926535898/180
    real, PARAMETER   :: ZERO = 0.0000001
    real, PARAMETER   :: LONS = 135.0    !135°E  [deg.]
    
    TM = HH

  call Julian_UTC(days_of_year,YY,MM,DD,HH,JDu,Tu)
  call Julian_TCG(days_of_year,YY,MM,DD,HH,Tu,T) 
  call right_ascension (Tu,asec)
  call celestial_longitude(T,asec,psi,eps,r,Et,Decl)
    SinD = sin( Decl * rpd )
    CosD = cos( Decl * rpd )
    hAngle = 15. * (TM - 12.) + (LON - LONS) + Et
    SinH = sin( PHI * rpd ) * SinD + cos( PHI * rpd ) * CosD * Cos(hAngle * rpd)
    CosH = sqrt( 1.0 - SinH * SinH )
    if ( ZERO > abs( CosH ) )then   ! H = PI/2(90deg.), A = 0 とみなす
        SinA = 0.0
        CosA = 1.0
        SinH = 1.0
        CosH = 0.0
    else 
        SinA = CosD * sin(hAngle * rpd) / CosH;
        CosA = sqrt( 1.0 - SinA * SinA )
    end if
end subroutine Sunlocation_Matsumoto
!--------------------------write real YY MM DD HH-------------------------
 subroutine write_date(days_of_year,YY,MM,DD,HH,YYu,MMu,DDu,HHu)
    implicit none
    integer, intent(in)   :: days_of_year   !days_of_year
    real, intent(in)   :: YY   !year  
    real, intent(in)   :: MM   !month
    real, intent(in)   :: DD   !day           can be 1.23 day
    real, intent(in)   :: HH   !hour 
    real, intent(out)  :: YYu  !real year
    real, intent(out)  :: MMu  !real month
    real, intent(out)  :: DDu  !real day
    real, intent(out)  :: HHu  !real hour
    if (HH < 0) then
        HHu = HH + 24 
        if (DD == 1) then 
            if (MM==2.or.MM==4.or.MM==6.or.MM==8.or.MM==9.or.MM==11)then
                DDu = 31
                MMu = MM - 1
                YYu = YY
            else if (MM==3) then
              if (days_of_year==366) then
                DDu = 29
                MMu = MM - 1
                YYu = YY
              else if (days_of_year==365) then
                DDu = 28
                MMu = MM - 1
                YYu = YY
              end if
            else if (MM==5.or.MM==7.or.MM==10.or.MM==12) then
                DDu = 30
                MMu = MM - 1
                YYu = YY
            else if (MM == 1) then
                DDu = 31
                MMu = 12
                YYu = YY -1
            end if
        else if (DD > 1) then
            DDu = DD -1
            MMu = MM
            YYu = YY
        end if
    else if (HH > 24) then
        HHu = HH - 24 
        if (DD == 31) then 
            if (MM==1.or.MM==3.or.MM==5.or.MM==7.or.MM==8.or.MM==10)then
                DDu = 1
                MMu = MM + 1
                YYu = YY
            else if (MM==12) then
                DDu = 1
                MMu = 1
                YYu = YY + 1
            end if
        else if (DD == 30) then
            if (MM==4.or.MM==6.or.MM==9.or.MM==11) then
                DDu = 1
                MMu = MM + 1
                YYu = YY
            end if
        else if (DD == 28 ) then
            if (days_of_year==366) then
                DDu = 29
                MMu = 2
                YYu = YY
            else if (days_of_year==356) then
                DDu = 1
                MMu = 3
                YYu = YY
            end if
        else if (DD == 29 ) then
                DDu = 1
                MMu = 3
                YYu = YY
        else 
            DDu = DD + 1
            MMu = MM
            YYu = YY
        end if
    else if (HH <= 24) then
        HHu = HH  
        DDu = DD
        MMu = MM
        YYu = YY
    end if
 end subroutine write_date
!---------------------------calculate Julian day/century ---------------------------------
 subroutine Julian_day(YY,MM,DD,HH,JD,JC)
    implicit none
    real, intent(in)   :: YY   !year  AFTER 1582.10.15
    real, intent(in)   :: MM   !month
    real, intent(in)   :: DD   !day 
    real, intent(in)   :: HH   !hour
    real(8), intent(out)  :: JD   !Julian day
    real, intent(out)  :: JC   !Julian century
    real(8)               :: A, B, Y, M, D
    Y = YY
    M = MM
    if (MM <= 2) then
        Y = Y - 1.0
        M = M + 12.0
    end if
    A  = nint(Y/100)
    B  = nint(A/4) - nint(Y/100) + 2
    D  = (HH - 12) / 24
    JD = floor(365.25*Y) + floor(30.6001*(M+1)) + DD +HH/12 +B + 1720994.5  !12:00 + 0.5
    !write(*,*)JD,d
    JC = (JD - 2451545) / 36535
end subroutine Julian_day
!-----------------and calculate the Julian century of TCG-----------------------------
subroutine Julian_UTC(days_of_year,YY,MM,DD,HH,JDu,Tu)
    implicit none
    integer, intent(in)   :: days_of_year   !days_of_year
    real, intent(in)   :: YY   !year  AFTER 1582.10.15
    real, intent(in)   :: MM   !month
    real, intent(in)   :: DD   !day    can be 1.23 day
    real, intent(in)   :: HH   !hour
    real, intent(out)  :: Tu   !Julian century of 
    real(8), intent(out)  :: JDu
    real               :: YYt   !year of TCG
    real               :: MMt   !month of TCG
    real               :: DDt   !day of TCG
    real               :: HHt   !hour of TCG
    real(8)               :: JDt   !Julian day of TCG
    real               :: T1   !time difference between TCG - UTC   [s]

    call write_date(days_of_year,YY,MM,DD,HH-9,YYt,MMt,DDt,HHt)
    call Julian_day(YYt,MMt,DDt,HHt,JDu,Tu)
    !write(*,*)JDu
 end subroutine Julian_UTC
!-----------------------change UTC(JST) to TCG time---------------------------------
!-----------------and calculate the Julian century of TCG-----------------------------
 subroutine Julian_TCG(days_of_year,YY,MM,DD,HH,Tu,T)
    implicit none
    integer, intent(in)   :: days_of_year   !days_of_year
    real, intent(in)   :: YY   !year  AFTER 1582.10.15
    real, intent(in)   :: MM   !month
    real, intent(in)   :: DD   !day    can be 1.23 day
    real, intent(in)   :: HH   !hour
    real, intent(in)   :: Tu   !Julian century of UTC
    real, intent(out)  :: T    !Julian century of TCG
    real               :: YYt   !year of TCG
    real               :: MMt   !month of TCG
    real               :: DDt   !day of TCG
    real               :: HHt   !hour of TCG
    real(8)               :: JDt   !Julian day of TCG
    real               :: T1   !time difference between TCG - UTC   [s]

    T1 = (80.84308/(1+0.2605601 * exp( -4.423790 * Tu )) - 0.311) / 3600  ![h]
    call write_date(days_of_year,YY,MM,DD,HH-9+T1,YYt,MMt,DDt,HHt)
    call Julian_day(YYt,MMt,DDt,HHt,JDt,T)
 end subroutine Julian_TCG
!------------------------calculate average right ascension(平均赤経)------------------------
 subroutine right_ascension(Tu,asec)
    implicit none
    real, intent(in)   :: Tu      !Julian century of UTC
    real, intent(out)  :: asec    !average right ascension [deg.]
    asec = (((0.093104 - 0.0000062*Tu)*Tu + 8640184.812866)*Tu + 67310.54841)/240 !15/3600=1/240[deg./s]
    !write(*,*)asec,tu
end subroutine right_ascension
!-----------------------calculate celestial longitude(視黄経)------------------------- 
!----------------------and 章動を含み黄道傾斜角 ; 動径 (地心距離)---------------------- 
 subroutine celestial_longitude(T,asec,psi,eps,r,Et,Decl)
    implicit none
    real, intent(in)  :: T                       !Julian century of TCG
    real, intent(in)  :: asec                    !average right ascension [deg.]
    real, intent(out) :: psi                     !celestial longitude(視黄経) [o]
    real, intent(out) :: eps                     !章動を含み黄道傾斜角
    real, intent(out) :: r                       !動径 (地心距離) (返値は AU 単位)
    real, intent(out) :: Et                      !均時差 (deg.)
    real, intent(out) :: Decl                    !(視) 赤緯 Decl(deg.)
    real              :: Pi(18), Qi(18), Ri(18)  !coefficient
    real              :: P1(9), Q1(9), R1(9)     !coefficient
    real              :: tanam, tpce, dAsc
    real              :: sind, cosd
    integer           :: I
    real, PARAMETER   :: rpd = 3.1415926535898/180
    real, PARAMETER   :: E1 = -89381.448000
    real, PARAMETER   :: E2 = 46.815000
    real, PARAMETER   :: E3 = 0.000590
    real, PARAMETER   :: E4 = -0.001813
    DATA Pi / 1.9147, 0.0200, 0.0020, 0.0018,  0.0018, 0.0015, &
              0.0013, 0.0007, 0.0007, 0.0007,  0.0006, 0.0005, &
              0.0005, 0.0004, 0.0004, -0.0048, 0.0048, -0.0004 /
    DATA Qi / 35999.05, 71998.10, 32964.00, 19.00,    445267.00, 45038.00, &
              22519.00, 65929.00, 3035.00,  9038.00,  33718.00,  155.00,   &
              2281.00,  29930.00, 31557.00, 35999.00, 1934.00,   72002.00  /
    DATA Ri / 267.52, 265.10, 158.00, 159.00, 208.00, 254.00, &
              352.00, 45.00,  110.00, 64.00,  316.00, 118.00, &
              221.00, 48.00,  161.00, 268.00, 145.00, 111.00  /
    DATA P1 / 1.000140, 0.016706, 0.000139, 0.000031,  0.000016, &
              0.000016, 0.000005, 0.000005, -0.000042  /
    DATA Q1 / 0.00, 35999.05, 71998.00, 445267.00, 32964.00, &
              45038.00, 22519.00, 33718.00, 35999.00  /
    DATA R1 / 0.00, 177.53, 175.00, 298.00,  68.00, &
              164.00, 233.00, 226.00, 178.00  /
    !!---章動を含み黄道傾斜角
    eps = -0.00256 *cos( mod( 1934.0 * T + 235.0, 360.0 ) * rpd )
    eps = eps - 0.00015 * cos( mod( 72002.0 * T + 201.0, 360.0 ) * rpd )
    eps = eps + ( ((E4 * T + E3) * T + E2) * T + E1 ) / 3600    !in [deg.]  
    eps = -eps  !!!松本の式(13)　逆にしたと考えます。
    !!---視黄経  
    psi = 0.0
    do I=1,18,1
        psi = Pi(I) * cos( mod( (Qi(I)*T + Ri(I)), 360.0 )*rpd ) + psi
        if (I == 16) then 
            psi = Pi(I) * T * cos( mod( (Qi(I)*T + Ri(I)), 360.0 )*rpd ) + psi
        end if
    end do
    psi = psi + 36000.7695 * T + 280.4602  -0.00569 !?
    !!---視黄経動径 (地心距離)
    r = 0
    do I=2,8,1
        r = P1(I) * cos( mod( (Q1(I)*T + R1(I)), 360.0 )*rpd ) + r
    end do
    r = P1(1) + P1(9) * T *cos( mod( (Q1(9)*T + R1(9)), 360.0 )*rpd ) + r
    !!---均時差 (deg.)
    Et = 0
    do I=17,18
        Et = Pi(I) * cos( mod( (Qi(I)*T + Ri(I)), 360.0 )*rpd ) + Et
    end do
    Et = ( Et - 0.00569 ) * cos( mod( eps, 360.0 ) * rpd )
    tanam = tan( asec * rpd )
    tpce  = tan( mod( psi, 360.0 ) * rpd ) * cos( mod( eps, 360.0 ) * rpd )
    dAsc = atan( (tanam - tpce) / (1.0 + tanam * tpce) ) / rpd    !rad.->deg.
    Et = mod( Et + dAsc, 360.0 )
    !!---(視) 赤緯 Decl(deg.)
    sind = sin( mod( psi, 360.0 ) * rpd ) * sin( mod( eps, 360.0 ) * rpd )
    cosd = sqrt( 1.0 - sind * sind )
    Decl = atan( sind / cosd ) / rpd     !rad.->deg.
    
 end subroutine celestial_longitude
!----------------------------------------------------------------------------------- 
  