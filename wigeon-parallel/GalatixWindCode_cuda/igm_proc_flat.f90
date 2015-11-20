module IGM_heating_cooling
  use datablock
  use readin_out
  
contains
  subroutine extQ2(fq,Q)
    real,dimension(5)::fq(5)
    integer,parameter::Nx=4
    real::Q,n_eff,den,cden,vx,vy,vz,eng,vk,pres,tem_H,tem10
    real::CHrate,CR
    real,dimension(Nx+1) ::temx,coolx
    data tem_min,tem_max/4.00,8.50/
    n_eff=x_hydrogen+y_helium/4.
    den=fq(1)
    cden=1./den
    vx=fq(2)*cden
    vy=fq(3)*cden
    vz=fq(4)*cden
    eng=fq(5)
    vk=0.5*(vx*vx+vy*vy+vz*vz)
    pres=gm1*(eng*cden-vk)
    tem_H=pres*tem_unit/n_eff   !need to be checked  ,seems that should divided by n_eff
    if(tem_H<0.) then
       Q=0
       return
    endif
    
    tem10=log10(tem_H)
    
    i_tem=int(float(n_cool-1)*(tem10-tem_min)/(tem_max-tem_min))
    
    Nx2=Nx/2
    im=i_tem-Nx2
    ip=i_tem+Nx2
    
    if(im>=1.and.ip<=(n_cool-1)) then
       do i=1,Nx+1  
          temx(i)=tem_min+float(i+im-1)*(tem_max-tem_min)/float(n_cool-1)
       enddo
       coolx(1:Nx+1)=coolrate(im:ip)
       call polint(temx,coolx,Nx+1,tem10,CR,df)
       CHrate=CR
       Q=cooling_coff*den*den*10**(CHrate)
    else
       CHrate=0.
       Q=0
    endif
  
    return
  end subroutine extQ2
  
       
  subroutine extQ(zf,xJ,fq,Q)
    real,dimension(5)::fq(5)
    real ::zf,n_eff,den,cden,vx,vy,vz,eng,vk,pres,den_H,Q,xJ
    
    if (zf>z_reion) then
       n_eff = x_hydrogen+y_helium/4.
    else
       n_eff = 2.*x_hydrogen +3. * y_helium /4.
    endif
    
    den=fq(1)
    cden=1./den
    vx=fq(2)*cden
    vy=fq(3)*cden
    vz=fq(4)*cden
    eng=fq(5)
    vk=0.5*(vx*vx+vy*vy+vz*vz)
    pres=gm1*(eng*cden-vk)*cf*cf
    tem_H = pres * T_unit * n_eff
    den_H = n_unit * den * cf * cf * cf 
    call GasModeling(zf,xJ,tem_H,den_H,CH_rate)  
    Q = cooling_coff * af * den * den * CH_rate
    
    return
  end subroutine extQ
  
! -----------------------------------------------------------------------
  subroutine GasModeling(z_redshift,xJ21,tem_igm,den_igm,CHrate)
    integer,parameter::Ns=1000,Nx=4
    real, parameter::uv_index=1. ! pow index of the UV spectrum 
    real :: tem_igm,den_igm,CR,CHrate
    data tem1,tem2,temlog_max,denlog_min,denlog_max/3.95,0.,8.,-8.,2.0/
    
    real,dimension(Ns)::tem,den
    real,dimension(Ns,Ns)::ch_table
    real,dimension(Nx+1)::temx,denx
    real,dimension(Nx+1,Nx+1)::chx
    
    integer :: kd1,kd2
    save tem,den,ch_table,kd1,kd2
    data kd1,kd2/0,0/
    
    if (xJ21==0.) tem_min=tem1
    if (xJ21>0.) tem_min=tem2
    
    if ((kd1==0.and.xJ21==0.).or.(kd2==0.and.xJ21>0.)) then  
       do j=1,Ns
          do i=1,Ns
             tem(i)=tem_min+(tem_max-tem_min)*float(i)/float(Ns)
             den(j)=denlog_min+(denlog_max-denlog_min)*float(j)/float(Ns)
             tem_i=10.**tem(i)
             den_j=10.**den(j)
             call IGMprocess(z_redshift,tem_i,den_j,xJ21,uv_index,CR)         
             ch_table(i,j)=CR
          enddo
       enddo
       if (kd1==0.and.xJ21==0.) then
          write(*,'(''generate cooling-heating table without UV '')')
          kd1=1
       endif
       if (kd2==0.and.xJ21>0.) then
          write(*,'(''generate cooling-heating table with UV'')')    
          kd2=1
       endif
    endif
    
    if (den_igm<=0..or.tem_igm<=tem_min) then
       CHrate=0.
       return
    endif
    
    tem10=log10(tem_igm)
    den10=log10(den_igm)
    
    i_tem=int(float(Ns)*(tem10-tem_min)/(tem_max-tem_min))
    j_den=int(float(Ns)*(den10-den_min)/(den_max-den_min))
    
    Nx2=Nx/2
    im=i_tem-Nx2
    ip=i_tem+Nx2
    jm=j_den-Nx2
    jp=j_den+Nx2
    
    if (im>=1.and.ip<=Ns.and.jm>=1.and.jp<=Ns) then
       temx(1:Nx+1)=tem(im:ip)
       denx(1:Nx+1)=den(jm:jp)
       chx(1:Nx+1,1:Nx+1)=ch_table(im:ip,jm:jp)
       CHrate=CR      
    else 
       CHrate=0.
    endif
    
    return
  end subroutine GasModeling
  
! ---------------------------- IGMprocess -------------------------------------

  subroutine IGMprocess(z_redshift,T,den_H,xJ21,uv_index,CR)
    real :: alpha_H1,alpha_He1,alpha_He2,alpha_D
    real :: gamma_eH0,gamma_eHe0,gamma_eHe1,gamma_ph0,gamma_pHe0,gamma_pHe1
    real :: xH0,xH1,xHe0,xHe1,xHe2,xe
    real :: z_redshift,T,den_H,xJ21,uv_index,CR
    

! computing the ionizing cofficients
    call ionizing_coff

! computing the fraction of each component in den_H
    call fraction_den

! computing the net cooling rate 
    call CoolingRate

    return
    
  contains
    
! ---------------------------- sub-subroutine ---------------------------------

! computing the fraction of each components 
    subroutine fraction_den
      real, parameter::eps=1.E-10
      
      y = 0.25 * Y_helium / X_helium
      
      if (xJ21 == 0.) then
         
         xH0 = alpha_H1 / ( alpha_H1 + gamma_eH0 )
         xH1 = 1. - xH0
         xHe1 = y /(1. + ( alpha_He1 + alpha_d ) /gamma_eHe0 + gamma_eHe1 /alpha_He2 )
         xHe0 = xHe1 * ( alpha_He1 + alpha_D) /gamma_eHe0
         xHe2 = xHe1 * gamma_eHe1 /alpha_He2
         xe = xH1 + xHe1 + 2. * xHe2
         
      else if (xJ21 > 0.) then 
         xe0 = 1.D0
10       xH0 = alpha_H1 /( alpha_H1 + gamma_eH0 + Gamma_pH0 /Xe0 /Den_H )
         xH1 = 1. - xH0
         xHe1 = y /( 1.+ ( alpha_He2 + alpha_d ) / ( gamma_eHe0 + gamma_pHe0 /Xe0 /den_H ) &
&      + ( gamma_eHe1 + gamma_pHe1 /xe0 /Den_H ) /alpha_He2 )
         xHe0 = xHe1 * ( alpha_He1 + alpha_D ) /( gamma_eHe0 + gamma_pHe0 /xe0 /den_H )
         xHe2 = xHe1 * ( gamma_eHe1 + gamma_pHe1 /xe0 /den_H ) /alpha_He2
         xe = xH1 + xHe1 + 2. * xHe2
         
         if (abs(xe-xe0)/xe0<eps) then
            return
         else 
            xe0 = xe
            goto 10
         endif
         
      endif
      
    end subroutine fraction_den

! ionizing rate in unit of 10^-10 cm^3 s^-1  
    subroutine ionizing_coff
      real:: T3,T5,T6,Ts,T5s
      
      T3=T/1000.
      T5=T3/100.
      T6=T5/10.
      Ts=sqrt(T)
      T5s=sqrt(T5s)
      
      alpha_H1 = 0.84 /Ts * T3**(-0.2) /(1.+T6**0.7)
      alpha_He1 = 1.5 * T**(-0.6353)
      alpha_He2 = 3.36 /Ts * T3**(-0.2) /(1.+ T6**0.7)
      alpha_d = (1.9E+7) * T**(-1.5) * exp(-4.7/T5) * (1.+0.3 * exp(-0.94/T5))
      
      gamma_eH0 = 0.585 * Ts * exp(-157809.1/T)/(1.+T5s)
      gamma_eHe0 = 0.238 * Ts * exp(-285335.4/T)/(1.+T5s)
      gamma_eHe1 = 0.0568 * Ts * exp(-631515.0/T)/(1.+T5s)
      
      if (xJ21>0.) then
         xn = uv_index
         gamma_pH0 = 0.126* xJ21 /(3.+xn)
         gamma_pHe0 = 0.148 * xJ21 * 0.553**xn * (1.66/(xn+2.05)-0.66/(xn+3.05))
         gamma_pHe1 = 0.0334 * xJ21 * 0.249**xn /(3.+xn)
      endif
      
      return
    end subroutine ionizing_coff

! cooling rate in unit of 10^-21 ergs cm^-3 s^-1   
    subroutine CoolingRate
      real :: T3,T5,T6,Ts,T5s,TT
      real :: CIC_H0,CIC_He0,CIC_He1,CIC_He2S,DRC_He1,FFC_ALL
      real :: RCC_H0,RCC_He1,RCC_He2,CEC_H0,CEC_He1,CMB
      
      T3=T/1000.
      T5=T3/100.
      T6=T5/10.
      Ts=sqrt(T)
      T5s=sqrt(T5)
      TT=Ts/(1.+T5s)
      
! collisional ionization     
      CIC_H0 = 1.27 * TT * exp(-157809.1/T) * xe * xH0
      CIC_He0 = 0.938 * TT * exp(-285335.4/T) * xe * xHe0
      CIC_He1 = 0.495 * TT * exp(-631515.0/T) * xe * xHe1
      CIC_He2S = (5.01D-6) * T**(-0.1687) / (1.+T5s) * exp(-55338.0/T) * xe * xe * xHe1 * den_H
      CIC = CIC_H0 + CIC_He0 + CIC_He1 + CIC_He2S

! recombination 	
      RCC_H0 = (8.70D-6) * Ts * T3**(-0.2) / (1.+T6**0.7) * xe * xH1
      RCC_He1 = (1.55D-5) * T**0.3647 * xe * xHe1
      RCC_He2 = (3.48D-5) * Ts * T3**(-0.2) / (1.+T6**0.7) * xe * xHe2
      RCC = RCC_H0 + RCC_He1 + RCC_He2

! collisional excitation         
      CEC_H0 = 750. * exp(-118348.0/T)/(1.+T5s)*xe*xH0
      CEC_He1 = (5.54D+4) * T**(-0.397) * exp(-473638.0/T) / (1.+T5s) * xe * xHe1
      CEC = CEC_H0 + CEC_He1

! dielectric recombination	
      DRC_He1 = (1.24D+8) * T**(-1.5) * exp(-4.7/T5) * (1.+0.3 * exp(-0.94/T5)) * xe * xHe1 
      DRC = DRC_He1

! free-free emission 
      gff = 1.1 + 0.34 * exp(-(5.5-log10(T))**2/3.0)
      FFC_ALL=(1.42E-6) * gff * Ts * ( xH1 + xHe1 + 4. * xHe2 ) * Xe
      FFC = FFC_ALL

! CMB + electron Compton scattering 
      CMB = (5.406D-10) * (T5-2.7E-5/(1. + z_redshift))*(1.+ z_redshift)**4 * xe / den_H

! total cooling rate
      CR = CIC + RCC + CEC + DRC + FFC + CMB

! UV vackground heating
      if (xJ21 > 0.) then
         xn = uv_index
         Heat_H0 = 0.291 * xJ21 /(2.+xn) /(3.+xn)
         Heat_He0 = 0.584 * xJ21 * 0.553**xn * (1.66/(xn+1.05)-2.32/(xn+2.05)+0.66/(xn+3.05))
         Heat_He1 = 0.292 * xJ21 * 0.249**xn /(2.+xn) /(3.+xn)
         HR = (xH0 * Heat_H0 + xHe1 * Heat_He1 + xHe0 * Heat_He0) /den_H
      endif

! the net cooling rate
      CR = CR - HR
           
      return
    end subroutine CoolingRate
    
  end subroutine IGMprocess
  
end module IGM_heating_cooling
