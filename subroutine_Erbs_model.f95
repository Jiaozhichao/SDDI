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