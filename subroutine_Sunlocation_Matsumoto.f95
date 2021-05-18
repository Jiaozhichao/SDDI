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