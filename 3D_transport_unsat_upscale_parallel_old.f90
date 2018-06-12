!!!**************** discription ***************
!!! A FVM parallel pseudo-implicit code for unsaturated advection-diffusion-reaction transport
!!! double-precision, wall is set in x-,y-directions, 1st and 2nd B.C. is used in z-direction
!!! stagger arrangment, biofilm donot move
!!! advection and diffusion terms are implicit
!!! reaction term of aqueous organic carbon (lc) and oxygen (lo) is semi-implicit to increase numerical stability and avoid iteration
!!! otherwise, reaction terms implement explicit method
!!!
!!! change alpha_c with theta
!!! combined gas phase with aqueous phase
!!! use a standard form of equation to code. other equations can be transformed to this standard equation, but be careful the minor difference
!!!
!!! all carbon-related mass, including CO2, are represented by carbon (c)
!!! this code is modified based on previous codes where \partial(sc)/\partial(t)=-km0*(Ka*sc-lc)*A. To avoid change this code, we modifed 
!!! \partial(sc)/\partial(t)=-km*(sc-Kcc*lc)*A=-km*Kcc*(1/Kcc*sc-lc)*A, in the upscaling paper where alpha=km*A=a1*Dlc*A. this change aims to match the published eqn.
!!! ignore the effects of water potential on microbial activity
!!! ignore the effects of water and gas percolation on diffusion
!!!
!!! developed by Zhifeng Yan
!!! 2015-07-28
!!!********************************************



Module GLOBAL
!!! for grid 
  Integer::i,j,k,Nx,Ny,Nz,Nt,LX,LY,LZ,ilen

!!! for 
  Double precision,Dimension(:,:,:),Allocatable::ne,tmp_ne,Se,S,Sm,theta,tmp_theta,psi,epsi,tmp_epsi
  Double precision,Dimension(:,:,:),Allocatable::u,v,w,c_o,c,c0,lc_o,lc,lcc,sc_o,sc,b_o,b,ba,k_m,alpha_c
  Double precision,Dimension(:,:,:),Allocatable::beta1_o,beta2_o,alpha1_o,alpha2_o,beta1_d,beta2_d,alpha1_d,alpha2_d
  Double precision,Dimension(:,:,:),Allocatable::to_o,to,lo_o,lo,go_o,go,d_o,d,d0,ld_o,ld,gd_o,gd,td_o,td
  Double precision,Dimension(:,:,:),Allocatable::tmp_lc_o,tmp_to_o,tmp_td_o,tmp_ld_o,tmp_gd_o
  Double precision,Dimension(:,:,:),Allocatable::tmp_lc,tmp_to,tmp_td
  Double precision,Dimension(:,:,:),Allocatable::Dlc_eff,Db_eff,Dlo_eff,Dgo_eff,Dld_eff,Dgd_eff
  Double precision,Dimension(:,:,:),Allocatable::tmp_Dlc_eff,tmp_Dlo_eff,tmp_Dgo_eff,tmp_Dld_eff,tmp_Dgd_eff
  Double precision,Dimension(:,:,:),Allocatable::tmp_beta1_o,tmp_beta2_o,tmp_alpha1_o,tmp_alpha2_o
  Double precision,Dimension(:,:,:),Allocatable::tmp_beta1_d,tmp_beta2_d,tmp_alpha1_d,tmp_alpha2_d
  logical, allocatable, dimension(:,:,:):: wall,tmp_wall
  Double Precision::tmo,tfld,tfgd,tmd,psi_opt,alpha2
  Double Precision::Ka1,Ka2,KH,R_igl,T_abs,PH,Khd,KPH,w1,Ka
  Double precision,Dimension(:,:),Allocatable::flo,fgo,mlo,mgo,fld,fgd,mld,mgd,to_d,to_t,td_d,td_t
  Double precision,Dimension(:,:,:),Allocatable::q,aa,bb,cc,dd,ee,ff,gg
  Double precision,Dimension(:),Allocatable::ma,mb,mc,md,me,mf,mg,HOHO,HH,QQ
  logical, allocatable, dimension(:):: WW


!!! for input and output
  Integer::lcbc_d,lcbc_t,bbc_d,bbc_t,tobc_d,tobc_t,tdbc_d,tdbc_t
  Integer::gdbc_d,gdbc_t,Max_step,Print_step,File_step,t0
  Double precision::dx,dy,dz,dt,Dlc,Db,Dld,Dgd,Dlo,Dgo,qc_0,qo_0,qd_0,Kc,Ko,Dk,Yd
  Double precision::lambda1,tau1,m1,mu,mmm_l,nnn_l,mmm_g,nnn_g,Kho,Kcc,rhos,OM,a1,b1,SA
  Double precision::lc_d,lc_t,b_d,b_t,lo_d,lo_t,go_d,go_t,ld_d,ld_t,gd_d,gd_t,omega
  Character(64)::PMname,OUT_DIR,IN_DIR,Rootname

!!! for computation
  Integer::istep,nxx,nyy,nzz,nn,pp
  Double Precision::starttime,endtime,totaltime,HH_eps,HH_err,HH_err_local
  Double precision::mean_por,Porlocal,Porsum,clocalmax,cmax,lclocalmax,lcmax,sclocalmax,scmax
  Double precision::blocalmax,bmax,tolocalmax,tomax,lolocalmax,lomax,golocalmax,gomax
  Double precision::tdlocalmax,tdmax,ldlocalmax,ldmax,gdlocalmax,gdmax
  Double precision::clocalsum,csum,lclocalsum,lcsum,sclocalsum,scsum,blocalsum,bsum,ldlocalsum,ldsum,gdlocalsum,gdsum
  Double precision::tolocalsum,tosum,lolocalsum,losum,golocalsum,gosum,tdlocalsum,tdsum
  Double precision,PARAMETER::Pi=3.1415927
  Double precision,PARAMETER::SMALL=1.e-20

!!! for parallel
  Integer,PARAMETER::MEMO=8
  INTEGER::Nstart,error,mype,npes,idex,iprev,inext,nallgrp


End Module GLOBAL




PROGRAM USMSSolver

  Use GLOBAL
  Implicit None
  include 'mpif.h'
  integer::stat(MPI_STATUS_SIZE)

  Call input

  nxx = int(0.5*Nx)+1
  nyy = int(0.5*Ny)+1
  nzz = int(0.5*Nz)+1
  Nt = Nx*Ny*Nz
  m1 = 1.-1./tau1
  Ka = 1./Kcc

  ilen = Nx*Ny

!!! memory issue

  Allocate (u(Nx+1,Ny,Nz),v(Nx,Ny+1,Nz),w(Nx,Ny,Nz+1))
  Allocate (c_o(Nx,Ny,Nz),c(Nx,Ny,Nz),c0(Nx,Ny,Nz),lc_o(Nx,Ny,Nz),lc(Nx,Ny,Nz),sc_o(Nx,Ny,Nz),sc(Nx,Ny,Nz))
  Allocate (b_o(Nx,Ny,Nz),b(Nx,Ny,Nz))
  Allocate (to_o(Nx,Ny,Nz),to(Nx,Ny,Nz),td_o(Nx,Ny,Nz),td(Nx,Ny,Nz),lo_o(Nx,Ny,Nz),lo(Nx,Ny,Nz),go_o(Nx,Ny,Nz),go(Nx,Ny,Nz))
  Allocate (ld_o(Nx,Ny,Nz),ld(Nx,Ny,Nz),gd_o(Nx,Ny,Nz),gd(Nx,Ny,Nz))
  Allocate (ne(Nx,Ny,Nz),Se(Nx,Ny,Nz),S(Nx,Ny,Nz),Sm(Nx,Ny,Nz),theta(Nx,Ny,Nz),epsi(Nx,Ny,Nz))
  Allocate (psi(Nx,Ny,Nz),ba(Nx,Ny,Nz),k_m(Nx,Ny,Nz),alpha_c(Nx,Ny,Nz))
  Allocate (beta1_o(Nx,Ny,Nz),beta2_o(Nx,Ny,Nz),alpha1_o(Nx,Ny,Nz),alpha2_o(Nx,Ny,Nz))
  Allocate (beta1_d(Nx,Ny,Nz),beta2_d(Nx,Ny,Nz),alpha1_d(Nx,Ny,Nz),alpha2_d(Nx,Ny,Nz))
  Allocate (Dlc_eff(Nx,Ny,Nz),Db_eff(Nx,Ny,Nz),Dlo_eff(Nx,Ny,Nz),Dgo_eff(Nx,Ny,Nz),Dld_eff(Nx,Ny,Nz),Dgd_eff(Nx,Ny,Nz))

  Allocate (tmp_lc_o(Nx,Ny,2),tmp_lc(Nx,Ny,2),tmp_to_o(Nx,Ny,2),tmp_to(Nx,Ny,2),tmp_td_o(Nx,Ny,2),tmp_td(Nx,Ny,2))
  Allocate (tmp_Dlc_eff(Nx,Ny,2),tmp_Dlo_eff(Nx,Ny,2),tmp_Dgo_eff(Nx,Ny,2),tmp_Dld_eff(Nx,Ny,2),tmp_Dgd_eff(Nx,Ny,2))
!  Allocate (tmp_ld_o(Nx,Ny,2),tmp_ld(Nx,Ny,2),tmp_gd_o(Nx,Ny,2),tmp_gd(Nx,Ny,2))
  Allocate (tmp_beta1_o(Nx,Ny,2),tmp_beta2_o(Nx,Ny,2),tmp_alpha1_o(Nx,Ny,2),tmp_alpha2_o(Nx,Ny,2))
  Allocate (tmp_beta1_d(Nx,Ny,2),tmp_beta2_d(Nx,Ny,2),tmp_alpha1_d(Nx,Ny,2),tmp_alpha2_d(Nx,Ny,2))
  Allocate (wall(Nx,Ny,Nz))
  Allocate (tmp_wall(Nx,Ny,2),tmp_ne(Nx,Ny,2),tmp_theta(Nx,Ny,2),tmp_epsi(Nx,Ny,2))
  Allocate (flo(Nx,Ny),fgo(Nx,Ny),fld(Nx,Ny),fgd(Nx,Ny),mlo(Nx,Ny),mgo(Nx,Ny),mld(Nx,Ny),mgd(Nx,Ny))
  Allocate (to_d(Nx,Ny),to_t(Nx,Ny),td_d(Nx,Ny),td_t(Nx,Ny))

  Allocate (q(Nx,Ny,Nz),aa(Nx,Ny,Nz),bb(Nx,Ny,Nz),cc(Nx,Ny,Nz),dd(Nx,Ny,Nz),ee(Nx,Ny,Nz),ff(Nx,Ny,Nz),gg(Nx,Ny,Nz))
  Allocate (ma(Nt),mb(Nt),mc(Nt),md(Nt),me(Nt),mf(Nt),mg(Nt),HOHO(Nt),HH(Nt),QQ(Nt),WW(Nt))




  starttime = MPI_WTIME()


!!! read porosity
  Call ReadDensity(ne(:,:,:),Nz,Nx,Ny,Nstart,MEMO,&
       Trim(IN_DIR)//'/'//Trim(PMname))
  
  
  !notice
  ne(1,:,:) = 0.
  ne(Nx,:,:) = 0.
  ne(:,1,:) = 0.
  ne(:,Ny,:) = 0.


  wall = .false. 
  do i=1,Nx
     do j=1,Ny
        do k=1,Nz  
           if (ne(i,j,k) .eq. 0.) then
              wall(i,j,k) = .true. 
           endif
        enddo
     enddo
  enddo


!!! exchange porosity and wall at interface

  call mpi_sendrecv(ne(:,:,Nz),ilen,MPI_Double_precision,inext,810, &
       tmp_ne(:,:,1),ilen,MPI_Double_precision,iprev,810,nallgrp,stat,error)
  call mpi_sendrecv(ne(:,:,1),ilen,MPI_Double_precision,iprev,820, &
       tmp_ne(:,:,2),ilen,MPI_Double_precision,inext,820,nallgrp,stat,error)

  call mpi_sendrecv(wall(:,:,Nz),ilen,MPI_LOGICAL,inext,310, &
       tmp_wall(:,:,1),ilen,MPI_LOGICAL,iprev,310,nallgrp,stat,error)
  call mpi_sendrecv(wall(:,:,1),ilen,MPI_LOGICAL,iprev,320, &
       tmp_wall(:,:,2),ilen,MPI_LOGICAL,inext,320,nallgrp,stat,error)

  CALL MPI_BARRIER(MPI_COMM_WORLD,error)



  Porlocal = sum(ne)
  CALL MPI_REDUCE(Porlocal,Porsum,1,MPI_Double_precision,MPI_SUM,0,MPI_COMM_WORLD,error)

  mean_por = Porsum/(LX-2)/(LY-2)/LZ



  IF (mype .eq. 0 ) THEN
     print*,'##'
     print*,'average porosity  =', mean_por
     Print*,'porosity---------','nxx=',nxx,'nyy=',nyy
     Print*,ne(nxx,nyy,:)
  Endif




!!! read initial velocity and pressure


  CALL ReadDensity(u(:,:,:),Nz,Nx+1,Ny,Nstart,MEMO,&
       trim(IN_DIR)//'/'//trim(PMname)//'_u')
  CALL ReadDensity(v(:,:,:),Nz,Nx,Ny+1,Nstart,MEMO,&
       trim(IN_DIR)//'/'//trim(PMname)//'_v')
  CALL ReadDensity(w(:,:,:),Nz+1,Nx,Ny,Nstart,MEMO,&
       trim(IN_DIR)//'/'//trim(PMname)//'_w')
  CALL ReadDensity(S(:,:,:),Nz,Nx,Ny,Nstart,MEMO,&
       trim(IN_DIR)//'/'//trim(PMname)//'_S')


  if(t0==0) then
!!$     CALL ReadDensity(u(:,:,:),Nz,Nx,Ny,Nstart,MEMO,&
!!$          trim(IN_DIR)//'/'//trim(PMname)//'_u')
!!$     CALL ReadDensity(v(:,:,:),Nz,Nx,Ny,Nstart,MEMO,&
!!$          trim(IN_DIR)//'/'//trim(PMname)//'_v')
!!$     CALL ReadDensity(w(:,:,:),Nz,Nx,Ny,Nstart,MEMO,&
!!$          trim(IN_DIR)//'/'//trim(PMname)//'_w')
     CALL ReadDensity(sc_o(:,:,:),Nz,Nx,Ny,Nstart,MEMO,&
          trim(IN_DIR)//'/'//trim(PMname)//'_sc')
     CALL ReadDensity(lc_o(:,:,:),Nz,Nx,Ny,Nstart,MEMO,&
          trim(IN_DIR)//'/'//trim(PMname)//'_lc')
     CALL ReadDensity(b_o(:,:,:),Nz,Nx,Ny,Nstart,MEMO,&
          trim(IN_DIR)//'/'//trim(PMname)//'_b')
     CALL ReadDensity(lo_o(:,:,:),Nz,Nx,Ny,Nstart,MEMO,&
          trim(IN_DIR)//'/'//trim(PMname)//'_lo')
     CALL ReadDensity(go_o(:,:,:),Nz,Nx,Ny,Nstart,MEMO,&
          trim(IN_DIR)//'/'//trim(PMname)//'_go')
     CALL ReadDensity(ld_o(:,:,:),Nz,Nx,Ny,Nstart,MEMO,&
          trim(IN_DIR)//'/'//trim(PMname)//'_ld')
     CALL ReadDensity(gd_o(:,:,:),Nz,Nx,Ny,Nstart,MEMO,&
          trim(IN_DIR)//'/'//trim(PMname)//'_gd')
  else
!!$
!!$     Rootname = trim(OUT_DIR)//'/'//trim(PMname)//'_u_'// &
!!$          't'// &   
!!$          CHAR(mod(t0,10000000)/1000000+48)//&
!!$          CHAR(mod(t0,1000000)/100000+48)//&
!!$          CHAR(mod(t0,100000)/10000+48)//&
!!$          CHAR(mod(t0,10000)/1000+48)//&
!!$          CHAR(mod(t0,1000)/100+48)//&
!!$          CHAR(mod(t0,100)/10+48)//&
!!$          CHAR(mod(t0,10)+48)
!!$     CALL ReadDensity(u(:,:,:),Nz,Nx,Ny,Nstart,MEMO,trim(Rootname))
!!$
!!$     Rootname = trim(OUT_DIR)//'/'//trim(PMname)//'_v_'// &
!!$          't'// &   
!!$          CHAR(mod(t0,10000000)/1000000+48)//&
!!$          CHAR(mod(t0,1000000)/100000+48)//&
!!$          CHAR(mod(t0,100000)/10000+48)//&
!!$          CHAR(mod(t0,10000)/1000+48)//&
!!$          CHAR(mod(t0,1000)/100+48)//&
!!$          CHAR(mod(t0,100)/10+48)//&
!!$          CHAR(mod(t0,10)+48)
!!$     CALL ReadDensity(v(:,:,:),Nz,Nx,Ny,Nstart,MEMO,trim(Rootname))
!!$
!!$     Rootname = trim(OUT_DIR)//'/'//trim(PMname)//'_w_'// &
!!$          't'// &   
!!$          CHAR(mod(t0,10000000)/1000000+48)//&
!!$          CHAR(mod(t0,1000000)/100000+48)//&
!!$          CHAR(mod(t0,100000)/10000+48)//&
!!$          CHAR(mod(t0,10000)/1000+48)//&
!!$          CHAR(mod(t0,1000)/100+48)//&
!!$          CHAR(mod(t0,100)/10+48)//&
!!$          CHAR(mod(t0,10)+48)
!!$     CALL ReadDensity(w(:,:,:),Nz,Nx,Ny,Nstart,MEMO,trim(Rootname))

     Rootname = trim(OUT_DIR)//'/'//trim(PMname)//'_sc_'// &
          't'// &   
          CHAR(mod(t0,10000000)/1000000+48)//&
          CHAR(mod(t0,1000000)/100000+48)//&
          CHAR(mod(t0,100000)/10000+48)//&
          CHAR(mod(t0,10000)/1000+48)//&
          CHAR(mod(t0,1000)/100+48)//&
          CHAR(mod(t0,100)/10+48)//&
          CHAR(mod(t0,10)+48)
     CALL ReadDensity(sc_o(:,:,:),Nz,Nx,Ny,Nstart,MEMO,trim(Rootname))

     Rootname = trim(OUT_DIR)//'/'//trim(PMname)//'_lc_'// &
          't'// &   
          CHAR(mod(t0,10000000)/1000000+48)//&
          CHAR(mod(t0,1000000)/100000+48)//&
          CHAR(mod(t0,100000)/10000+48)//&
          CHAR(mod(t0,10000)/1000+48)//&
          CHAR(mod(t0,1000)/100+48)//&
          CHAR(mod(t0,100)/10+48)//&
          CHAR(mod(t0,10)+48)
     CALL ReadDensity(lc_o(:,:,:),Nz,Nx,Ny,Nstart,MEMO,trim(Rootname))


     Rootname = trim(OUT_DIR)//'/'//trim(PMname)//'_b_'// &
          't'// &   
          CHAR(mod(t0,10000000)/1000000+48)//&
          CHAR(mod(t0,1000000)/100000+48)//&
          CHAR(mod(t0,100000)/10000+48)//&
          CHAR(mod(t0,10000)/1000+48)//&
          CHAR(mod(t0,1000)/100+48)//&
          CHAR(mod(t0,100)/10+48)//&
          CHAR(mod(t0,10)+48)
     CALL ReadDensity(b_o(:,:,:),Nz,Nx,Ny,Nstart,MEMO,trim(Rootname))

     Rootname = trim(OUT_DIR)//'/'//trim(PMname)//'_lo_'// &
          't'// &   
          CHAR(mod(t0,10000000)/1000000+48)//&
          CHAR(mod(t0,1000000)/100000+48)//&
          CHAR(mod(t0,100000)/10000+48)//&
          CHAR(mod(t0,10000)/1000+48)//&
          CHAR(mod(t0,1000)/100+48)//&
          CHAR(mod(t0,100)/10+48)//&
          CHAR(mod(t0,10)+48)
     CALL ReadDensity(lo_o(:,:,:),Nz,Nx,Ny,Nstart,MEMO,trim(Rootname))

     Rootname = trim(OUT_DIR)//'/'//trim(PMname)//'_go_'// &
          't'// &   
          CHAR(mod(t0,10000000)/1000000+48)//&
          CHAR(mod(t0,1000000)/100000+48)//&
          CHAR(mod(t0,100000)/10000+48)//&
          CHAR(mod(t0,10000)/1000+48)//&
          CHAR(mod(t0,1000)/100+48)//&
          CHAR(mod(t0,100)/10+48)//&
          CHAR(mod(t0,10)+48)
     CALL ReadDensity(go_o(:,:,:),Nz,Nx,Ny,Nstart,MEMO,trim(Rootname))

     Rootname = trim(OUT_DIR)//'/'//trim(PMname)//'_ld_'// &
          't'// &   
          CHAR(mod(t0,10000000)/1000000+48)//&
          CHAR(mod(t0,1000000)/100000+48)//&
          CHAR(mod(t0,100000)/10000+48)//&
          CHAR(mod(t0,10000)/1000+48)//&
          CHAR(mod(t0,1000)/100+48)//&
          CHAR(mod(t0,100)/10+48)//&
          CHAR(mod(t0,10)+48)
     CALL ReadDensity(ld_o(:,:,:),Nz,Nx,Ny,Nstart,MEMO,trim(Rootname))

     Rootname = trim(OUT_DIR)//'/'//trim(PMname)//'_gd_'// &
          't'// &   
          CHAR(mod(t0,10000000)/1000000+48)//&
          CHAR(mod(t0,1000000)/100000+48)//&
          CHAR(mod(t0,100000)/10000+48)//&
          CHAR(mod(t0,10000)/1000+48)//&
          CHAR(mod(t0,1000)/100+48)//&
          CHAR(mod(t0,100)/10+48)//&
          CHAR(mod(t0,10)+48)
     CALL ReadDensity(gd_o(:,:,:),Nz,Nx,Ny,Nstart,MEMO,trim(Rootname))

!!$
!!$     Rootname = trim(OUT_DIR)//'/'//trim(PMname)//'_S_'// &
!!$          't'// &   
!!$          CHAR(mod(t0,10000000)/1000000+48)//&
!!$          CHAR(mod(t0,1000000)/100000+48)//&
!!$          CHAR(mod(t0,100000)/10000+48)//&
!!$          CHAR(mod(t0,10000)/1000+48)//&
!!$          CHAR(mod(t0,1000)/100+48)//&
!!$          CHAR(mod(t0,100)/10+48)//&
!!$          CHAR(mod(t0,10)+48)
!!$     CALL ReadDensity(S(:,:,:),Nz,Nx,Ny,Nstart,MEMO,trim(Rootname))


  endif






  IF (mype .eq. 0 ) THEN
     print*,'##'
     print*,'initilization done'
  Endif 




!!$!!! calculate co2 henry's constant according to equilibrium reaction 
!!$!!! consider effects of PH on CO2 flux and production
!!$  
!!$  Ka1 = 10**(-6.3)
!!$  Ka2 = 10**(-10.25)
!!$  KH = 10**(-1.5)      ! equilibrium constant between H2CO3* and CO2(g)
!!$  R_igl = 8.314
!!$  T_abs = 298
!!$  PH = 6.8
!!$
!!$  Khd = 12./4400.*R_igl*T_abs*KH
!!$  KPH = 1.+Ka1/10**(-PH)+Ka1*Ka2/(10**(-PH))**2
!!$
!!$ !notice
!!$  IF (mype .eq. 0 ) THEN
!!$     print*,'Khd=',Khd,'KPH=',KPH
!!$  Endif






!!! calculate hydrodynamic interaction

  do i=1,Nx
     do j=1,Ny
        do k=1,Nz  
           if(.not.wall(i,j,k)) then
              theta(i,j,k) = ne(i,j,k)*S(i,j,k)
              epsi(i,j,k) = ne(i,j,k)*(1.0-S(i,j,k))

              beta1_o(i,j,k) = Kho*ne(i,j,k)/(theta(i,j,k)*Kho+epsi(i,j,k))
              beta2_o(i,j,k) = ne(i,j,k)/(theta(i,j,k)*Kho+epsi(i,j,k))
              alpha1_o(i,j,k) = theta(i,j,k)*beta1_o(i,j,k)/ne(i,j,k)
              alpha2_o(i,j,k) = epsi(i,j,k)*beta2_o(i,j,k)/ne(i,j,k)

              beta1_d(i,j,k) = KPH*Khd*ne(i,j,k)/(theta(i,j,k)*KPH*Khd+epsi(i,j,k))
              beta2_d(i,j,k) = ne(i,j,k)/(theta(i,j,k)*KPH*Khd+epsi(i,j,k))
              alpha1_d(i,j,k) = theta(i,j,k)*beta1_d(i,j,k)/ne(i,j,k)
              alpha2_d(i,j,k) = epsi(i,j,k)*beta2_d(i,j,k)/ne(i,j,k)
           endif
        enddo
     enddo
  enddo






!!! exchange water content at interface (need to put in cycle if theta changes)

  call mpi_sendrecv(theta(:,:,Nz),ilen,MPI_Double_precision,inext,910, &
       tmp_theta(:,:,1),ilen,MPI_Double_precision,iprev,910,nallgrp,stat,error)
  call mpi_sendrecv(theta(:,:,1),ilen,MPI_Double_precision,iprev,920, &
       tmp_theta(:,:,2),ilen,MPI_Double_precision,inext,920,nallgrp,stat,error)

  call mpi_sendrecv(epsi(:,:,Nz),ilen,MPI_Double_precision,inext,930, &
       tmp_epsi(:,:,1),ilen,MPI_Double_precision,iprev,930,nallgrp,stat,error)
  call mpi_sendrecv(epsi(:,:,1),ilen,MPI_Double_precision,iprev,940, &
       tmp_epsi(:,:,2),ilen,MPI_Double_precision,inext,940,nallgrp,stat,error)




  call mpi_sendrecv(beta1_o(:,:,Nz),ilen,MPI_double_precision,inext,950, &
       tmp_beta1_o(:,:,1),ilen,MPI_double_precision,iprev,950,nallgrp,stat,error)
  call mpi_sendrecv(beta1_o(:,:,1),ilen,MPI_double_precision,iprev,960, &
       tmp_beta1_o(:,:,2),ilen,MPI_double_precision,inext,960,nallgrp,stat,error)

  call mpi_sendrecv(beta2_o(:,:,Nz),ilen,MPI_double_precision,inext,970, &
       tmp_beta2_o(:,:,1),ilen,MPI_double_precision,iprev,970,nallgrp,stat,error)
  call mpi_sendrecv(beta2_o(:,:,1),ilen,MPI_double_precision,iprev,980, &
       tmp_beta2_o(:,:,2),ilen,MPI_double_precision,inext,980,nallgrp,stat,error)

  call mpi_sendrecv(alpha1_o(:,:,Nz),ilen,MPI_double_precision,inext,150, &
       tmp_alpha1_o(:,:,1),ilen,MPI_double_precision,iprev,150,nallgrp,stat,error)
  call mpi_sendrecv(alpha1_o(:,:,1),ilen,MPI_double_precision,iprev,160, &
       tmp_alpha1_o(:,:,2),ilen,MPI_double_precision,inext,160,nallgrp,stat,error)

  call mpi_sendrecv(alpha2_o(:,:,Nz),ilen,MPI_double_precision,inext,170, &
       tmp_alpha2_o(:,:,1),ilen,MPI_double_precision,iprev,170,nallgrp,stat,error)
  call mpi_sendrecv(alpha2_o(:,:,1),ilen,MPI_double_precision,iprev,180, &
       tmp_alpha2_o(:,:,2),ilen,MPI_double_precision,inext,180,nallgrp,stat,error)




  call mpi_sendrecv(beta1_d(:,:,Nz),ilen,MPI_double_precision,inext,1950, &
       tmp_beta1_d(:,:,1),ilen,MPI_double_precision,iprev,1950,nallgrp,stat,error)
  call mpi_sendrecv(beta1_d(:,:,1),ilen,MPI_double_precision,iprev,1960, &
       tmp_beta1_d(:,:,2),ilen,MPI_double_precision,inext,1960,nallgrp,stat,error)

  call mpi_sendrecv(beta2_d(:,:,Nz),ilen,MPI_double_precision,inext,1970, &
       tmp_beta2_d(:,:,1),ilen,MPI_double_precision,iprev,1970,nallgrp,stat,error)
  call mpi_sendrecv(beta2_d(:,:,1),ilen,MPI_double_precision,iprev,1980, &
       tmp_beta2_d(:,:,2),ilen,MPI_double_precision,inext,1980,nallgrp,stat,error)

  call mpi_sendrecv(alpha1_d(:,:,Nz),ilen,MPI_double_precision,inext,1150, &
       tmp_alpha1_d(:,:,1),ilen,MPI_double_precision,iprev,1150,nallgrp,stat,error)
  call mpi_sendrecv(alpha1_d(:,:,1),ilen,MPI_double_precision,iprev,1160, &
       tmp_alpha1_d(:,:,2),ilen,MPI_double_precision,inext,1160,nallgrp,stat,error)

  call mpi_sendrecv(alpha2_d(:,:,Nz),ilen,MPI_double_precision,inext,1170, &
       tmp_alpha2_d(:,:,1),ilen,MPI_double_precision,iprev,1170,nallgrp,stat,error)
  call mpi_sendrecv(alpha2_d(:,:,1),ilen,MPI_double_precision,iprev,1180, &
       tmp_alpha2_d(:,:,2),ilen,MPI_double_precision,inext,1180,nallgrp,stat,error)

  CALL MPI_BARRIER(MPI_COMM_WORLD,error)






!!! calculate initial total oxygen and inorganic carbon conc.

  Do k=1,Nz
     Do j=1,Ny
        Do i=1,Nx
           if(.not.wall(i,j,k)) then
              to_o(i,j,k) = (theta(i,j,k)*lo_o(i,j,k) + epsi(i,j,k)*go_o(i,j,k))/ne(i,j,k)
              td_o(i,j,k) = (theta(i,j,k)*ld_o(i,j,k) + epsi(i,j,k)*gd_o(i,j,k))/ne(i,j,k)
           endif
        Enddo
     Enddo
  Enddo



  
!!! calculate boundary value for to and ?
  
  if(mype.eq.0) then
     Do j=1,Ny
        Do i=1,Nx
           if(.not.wall(i,j,1)) then
              to_d(i,j) = (theta(i,j,1)*lo_d + epsi(i,j,1)*go_d)/ne(i,j,1)
              td_d(i,j) = (theta(i,j,1)*ld_d + epsi(i,j,1)*gd_d)/ne(i,j,1)
           endif
        Enddo
     Enddo
  endif
  
  if(mype.eq.npes-1) then
     Do j=1,Ny
        Do i=1,Nx
           if(.not.wall(i,j,Nz)) then
              to_t(i,j) = (theta(i,j,Nz)*lo_t + epsi(i,j,Nz)*go_t)/ne(i,j,Nz)
              td_t(i,j) = (theta(i,j,Nz)*ld_t + epsi(i,j,Nz)*gd_t)/ne(i,j,Nz)
           endif
        Enddo
     Enddo
  endif




!!$!!! calculate initial total carbon conc.
!!$
!!$  Do k=1,Nz
!!$     Do j=1,Ny
!!$        Do i=1,Nx
!!$           if(wall(i,j,k)) then
!!$              c0(i,j,k) = 0.0
!!$              d0(i,j,k) = 0.0
!!$           else
!!$              c0(i,j,k) = lc_o(i,j,k)*theta(i,j,k) + sc_o(i,j,k)*rhos*(1.-ne(i,j,k))
!!$              d0(i,j,k) = ld_o(i,j,k)*theta(i,j,k) + gd_o(i,j,k)*epsi(i,j,k)*12./44.
!!$           endif
!!$        Enddo
!!$     Enddo
!!$  Enddo





!notde
  print*,'S =',S(nxx,nyy,nzz)





!!! calculate other effective diffusion coeff 
  do i=1,Nx
     do j=1,Ny
        do k=1,Nz   
           if(wall(i,j,k)) then
              Dlc_eff(i,j,k) = 0.
              Db_eff(i,j,k) = 0.
              Dlo_eff(i,j,k) = 0.
              Dgo_eff(i,j,k) = 0.
              Dld_eff(i,j,k) = 0.
              Dgd_eff(i,j,k) = 0.
           else
              Dlc_eff(i,j,k) = ne(i,j,k)**mmm_l * &
                   (theta(i,j,k)/ne(i,j,k))**nnn_l * Dlc
              Db_eff(i,j,k) = (ne(i,j,k))**mmm_l * &
                   ((theta(i,j,k))/(ne(i,j,k)))**nnn_l * Db
              Dlo_eff(i,j,k) = (ne(i,j,k))**mmm_l * &
                   ((theta(i,j,k))/(ne(i,j,k)))**nnn_l * Dlo
              Dgo_eff(i,j,k) = ne(i,j,k)**mmm_g*((epsi(i,j,k))/(ne(i,j,k)))**nnn_g * Dgo
              Dld_eff(i,j,k) = (ne(i,j,k))**mmm_l * &
                   ((theta(i,j,k))/(ne(i,j,k)))**nnn_l * Dld
              Dgd_eff(i,j,k) = ne(i,j,k)**mmm_g*((epsi(i,j,k))/(ne(i,j,k)))**nnn_g * Dgd
           endif
        enddo
     enddo
  enddo



     call mpi_sendrecv(Dlc_eff(:,:,Nz),ilen,MPI_double_precision,inext,130, &
          tmp_Dlc_eff(:,:,1),ilen,MPI_double_precision,iprev,130,nallgrp,stat,error)
     call mpi_sendrecv(Dlc_eff(:,:,1),ilen,MPI_double_precision,iprev,140, &
          tmp_Dlc_eff(:,:,2),ilen,MPI_double_precision,inext,140,nallgrp,stat,error)
     


     call mpi_sendrecv(Dlo_eff(:,:,Nz),ilen,MPI_double_precision,inext,130, &
          tmp_Dlo_eff(:,:,1),ilen,MPI_double_precision,iprev,130,nallgrp,stat,error)
     call mpi_sendrecv(Dlo_eff(:,:,1),ilen,MPI_double_precision,iprev,140, &
          tmp_Dlo_eff(:,:,2),ilen,MPI_double_precision,inext,140,nallgrp,stat,error)
     
     call mpi_sendrecv(Dgo_eff(:,:,Nz),ilen,MPI_double_precision,inext,190, &
          tmp_Dgo_eff(:,:,1),ilen,MPI_double_precision,iprev,190,nallgrp,stat,error)
     call mpi_sendrecv(Dgo_eff(:,:,1),ilen,MPI_double_precision,iprev,200, &
          tmp_Dgo_eff(:,:,2),ilen,MPI_double_precision,inext,200,nallgrp,stat,error)
     


     call mpi_sendrecv(Dld_eff(:,:,Nz),ilen,MPI_double_precision,inext,130, &
          tmp_Dld_eff(:,:,1),ilen,MPI_double_precision,iprev,130,nallgrp,stat,error)
     call mpi_sendrecv(Dld_eff(:,:,1),ilen,MPI_double_precision,iprev,140, &
          tmp_Dld_eff(:,:,2),ilen,MPI_double_precision,inext,140,nallgrp,stat,error)
     
     call mpi_sendrecv(Dgd_eff(:,:,Nz),ilen,MPI_double_precision,inext,190, &
          tmp_Dgd_eff(:,:,1),ilen,MPI_double_precision,iprev,190,nallgrp,stat,error)
     call mpi_sendrecv(Dgd_eff(:,:,1),ilen,MPI_double_precision,iprev,200, &
          tmp_Dgd_eff(:,:,2),ilen,MPI_double_precision,inext,200,nallgrp,stat,error)
     





  do i=1,Nx
     do j=1,Ny
        do k=1,Nz  
           if(.not.wall(i,j,k)) then
              k_m(i,j,k) = a1*Dlc*theta(i,j,k)/(theta(i,j,k)+b1)   
              alpha_c(i,j,k) = rhos*(1.-ne(i,j,k))*k_m(i,j,k)/theta(i,j,k)*SA*Kcc  !!! alpha_c=rhos*(1-ne)*k_m/theta*SA
           endif
        enddo
     enddo
  enddo

!!$!notde 
!!$alpha_1=1.
!!$alpha_2=0.


!notde
!print*,ne(nxx,nyy,nzz),theta(nxx,nyy,nzz),'alpha_c=',alpha_c(nxx,nyy,nzz)





  
  
!!! set values for walls
  do i=1,Nx
     do j=1,Ny
        do k=1,Nz  
           if(wall(i,j,k)) then
              sc_o(i,j,k) = 0.
              lc_o(i,j,k) = 0.
              b_o(i,j,k) = 0.
              lo_o(i,j,k) = 0.
              ld_o(i,j,k) = 0.
              go_o(i,j,k) = 0.
              gd_o(i,j,k) = 0.
              to_o(i,j,k) = 0.
              td_o(i,j,k) = 0.
           endif
        enddo
     enddo
  enddo




  tmo = 0.
  tmd = 0.



  IF (mype.eq.0) THEN
     CALL WriteMatlab()
     print*,'begin computation'
  ENDIF


!!!************* begin to compute ************

  Do istep=t0+1,Max_step


     IF(mod(istep,Print_step)==0) THEN
        IF (mype .eq. 0 ) THEN
           print*,'-------------------------'
           print*,'istep=', istep
        Endif
     endif





!!!! solve lc semi-implicitly
!!!! the advection and diffusion terms are completely implicit
!!!! the reaction term is dealt with semi-implicit to increase numerical stability
!!!! but avoid iteration (lineralization)



!!! boundary condition setting and
!!! communicate info between bottom and top planes


     call mpi_sendrecv(lc_o(:,:,Nz),ilen,MPI_double_precision,inext,110, &
          tmp_lc_o(:,:,1),ilen,MPI_double_precision,iprev,110,nallgrp,stat,error)
     call mpi_sendrecv(lc_o(:,:,1),ilen,MPI_double_precision,iprev,120, &
          tmp_lc_o(:,:,2),ilen,MPI_double_precision,inext,120,nallgrp,stat,error)

     
     CALL MPI_BARRIER(MPI_COMM_WORLD,error)

     


     
     
     if(npes.eq.1) then

        Do i=2,Nx-1
           Do j=2,Ny-1
              if(wall(i,j,1)) then
                 aa(i,j,1) = 0.
                 bb(i,j,1) = 0.
                 cc(i,j,1) = 0.
                 dd(i,j,1) = 0.
                 ee(i,j,1) = 0.
                 ff(i,j,1) = 0.
                 gg(i,j,1) = 0.
                 q(i,j,1) = 0.
              else
                 bb(i,j,1) = -dt*(max(u(i,j,1),0.)+(Dlc_eff(i-1,j,1)+Dlc_eff(i,j,1))/2./dx) * &
                      (theta(i-1,j,1)+theta(i,j,1))/2./dx
                 cc(i,j,1) = dt*(min(u(i+1,j,1),0.)-(Dlc_eff(i+1,j,1)+Dlc_eff(i,j,1))/2./dx) * &
                      (theta(i+1,j,1)+theta(i,j,1))/2./dx
                 dd(i,j,1) = -dt*(max(v(i,j,1),0.)+(Dlc_eff(i,j-1,1)+Dlc_eff(i,j,1))/2./dy) * &
                      (theta(i,j-1,1)+theta(i,j,1))/2./dy
                 ee(i,j,1) = dt*(min(v(i,j+1,1),0.)-(Dlc_eff(i,j+1,1)+Dlc_eff(i,j,1))/2./dy) * &
                      (theta(i,j+1,1)+theta(i,j,1))/2./dy
                 ff(i,j,1) = 0.
                 gg(i,j,1) = dt*(min(w(i,j,1+1),0.)-(Dlc_eff(i,j,1+1)+Dlc_eff(i,j,1))/2./dz) * &
                      (theta(i,j,1+1)+theta(i,j,1))/2./dz

                 if(lcbc_d.eq.1) then 
                    aa(i,j,1) = theta(i,j,1) + dt*(-min(u(i,j,1),0.)*(theta(i-1,j,1)+theta(i,j,1))/2./dx + &
                         (Dlc_eff(i-1,j,1)+Dlc_eff(i,j,1))/2./dx*(theta(i-1,j,1)+theta(i,j,1))/2./dx + &
                         max(u(i+1,j,1),0.)*(theta(i+1,j,1)+theta(i,j,1))/2./dx + &
                         (Dlc_eff(i+1,j,1)+Dlc_eff(i,j,1))/2./dx*(theta(i+1,j,1)+theta(i,j,1))/2./dx - &
                         min(v(i,j,1),0.)*(theta(i,j-1,1)+theta(i,j,1))/2./dy + &
                         (Dlc_eff(i,j-1,1)+Dlc_eff(i,j,1))/2./dy*(theta(i,j-1,1)+theta(i,j,1))/2./dy + &
                         max(v(i,j+1,1),0.)*(theta(i,j+1,1)+theta(i,j,1))/2./dy + &
                         (Dlc_eff(i,j+1,1)+Dlc_eff(i,j,1))/2./dy*(theta(i,j+1,1)+theta(i,j,1))/2./dy + &
                         Dlc_eff(i,j,1)/dz*2.*theta(i,j,1)/dz + &
                         max(w(i,j,1+1),0.)*(theta(i,j,1+1)+theta(i,j,1))/2./dz + &
                         (Dlc_eff(i,j,1+1)+Dlc_eff(i,j,1))/2./dz*(theta(i,j,1+1)+theta(i,j,1))/2./dz) - &
                         dt*theta(i,j,1)*(-qc_0*b_o(i,j,1)/(Kc+lc_o(i,j,1))*lo_o(i,j,1)/(Ko+lo_o(i,j,1))) + &
                         dt*theta(i,j,1)*alpha_c(i,j,1)
                    q(i,j,1) = theta(i,j,1)*lc_o(i,j,1) + &
                         dt*theta(i,j,1)*alpha_c(i,j,1)*Ka*sc_o(i,j,1) + &
                         dt*(w(i,j,1)+2.*Dlc_eff(i,j,1)/dz)*theta(i,j,1)/dz*lc_d
                 endif
                 if(lcbc_d.eq.2) then 

                    aa(i,j,1) = theta(i,j,1) + dt*(-min(u(i,j,1),0.)*(theta(i-1,j,1)+theta(i,j,1))/2./dx + &
                         (Dlc_eff(i-1,j,1)+Dlc_eff(i,j,1))/2./dx*(theta(i-1,j,1)+theta(i,j,1))/2./dx + &
                         max(u(i+1,j,1),0.)*(theta(i+1,j,1)+theta(i,j,1))/2./dx + &
                         (Dlc_eff(i+1,j,1)+Dlc_eff(i,j,1))/2./dx*(theta(i+1,j,1)+theta(i,j,1))/2./dx - &
                         min(v(i,j,1),0.)*(theta(i,j-1,1)+theta(i,j,1))/2./dy + &
                         (Dlc_eff(i,j-1,1)+Dlc_eff(i,j,1))/2./dy*(theta(i,j-1,1)+theta(i,j,1))/2./dy + &
                         max(v(i,j+1,1),0.)*(theta(i,j+1,1)+theta(i,j,1))/2./dy + &
                         (Dlc_eff(i,j+1,1)+Dlc_eff(i,j,1))/2./dy*(theta(i,j+1,1)+theta(i,j,1))/2./dy - &
                         w(i,j,1)*theta(i,j,1)/dz + &
                         max(w(i,j,1+1),0.)*(theta(i,j,1+1)+theta(i,j,1))/2./dz + &
                         (Dlc_eff(i,j,1+1)+Dlc_eff(i,j,1))/2./dz*(theta(i,j,1+1)+theta(i,j,1))/2./dz) - &
                         dt*theta(i,j,1)*(-qc_0*b_o(i,j,1)/(Kc+lc_o(i,j,1))*lo_o(i,j,1)/(Ko+lo_o(i,j,1))) + &
                         dt*theta(i,j,1)*alpha_c(i,j,1)
                    q(i,j,1) = theta(i,j,1)*lc_o(i,j,1)  + &
                         dt*theta(i,j,1)*alpha_c(i,j,1)*Ka*sc_o(i,j,1)
                 endif

              endif
           enddo
        enddo


        Do i=2,Nx-1
           Do j=2,Ny-1
              if(wall(i,j,Nz)) then
                 aa(i,j,Nz) = 0.
                 bb(i,j,Nz) = 0.
                 cc(i,j,Nz) = 0.
                 dd(i,j,Nz) = 0.
                 ee(i,j,Nz) = 0.
                 ff(i,j,Nz) = 0.
                 gg(i,j,Nz) = 0.
                 q(i,j,Nz) = 0.
              else     
                 bb(i,j,Nz) = -dt*(max(u(i,j,Nz),0.)+(Dlc_eff(i-1,j,Nz)+Dlc_eff(i,j,Nz))/2./dx) * &
                      (theta(i-1,j,Nz)+theta(i,j,Nz))/2./dx
                 cc(i,j,Nz) = dt*(min(u(i+1,j,Nz),0.)-(Dlc_eff(i+1,j,Nz)+Dlc_eff(i,j,Nz))/2./dx) * &
                      (theta(i+1,j,Nz)+theta(i,j,Nz))/2./dx
                 dd(i,j,Nz) = -dt*(max(v(i,j,Nz),0.)+(Dlc_eff(i,j-1,Nz)+Dlc_eff(i,j,Nz))/2./dy) * &
                      (theta(i,j-1,Nz)+theta(i,j,Nz))/2./dy
                 ee(i,j,Nz) = dt*(min(v(i,j+1,Nz),0.)-(Dlc_eff(i,j+1,Nz)+Dlc_eff(i,j,Nz))/2./dy) * &
                      (theta(i,j+1,Nz)+theta(i,j,Nz))/2./dy
                 ff(i,j,Nz) = -dt*(max(w(i,j,Nz),0.)+(Dlc_eff(i,j,Nz-1)+Dlc_eff(i,j,Nz))/2./dz) * &
                      (theta(i,j,Nz-1)+theta(i,j,Nz))/2./dz
                 gg(i,j,Nz) = 0.

                 if(lcbc_t.eq.1) then
                    aa(i,j,Nz) = theta(i,j,Nz) + dt*(-min(u(i,j,Nz),0.)*(theta(i-1,j,Nz)+theta(i,j,Nz))/2./dx + &
                         (Dlc_eff(i-1,j,Nz)+Dlc_eff(i,j,Nz))/2./dx*(theta(i-1,j,Nz)+theta(i,j,Nz))/2./dx + &
                         max(u(i+1,j,Nz),0.)*(theta(i+1,j,Nz)+theta(i,j,Nz))/2./dx + &
                         (Dlc_eff(i+1,j,Nz)+Dlc_eff(i,j,Nz))/2./dx*(theta(i+1,j,Nz)+theta(i,j,Nz))/2./dx - &
                         min(v(i,j,Nz),0.)*(theta(i,j-1,Nz)+theta(i,j,Nz))/2./dy + &
                         (Dlc_eff(i,j-1,Nz)+Dlc_eff(i,j,Nz))/2./dy*(theta(i,j-1,Nz)+theta(i,j,Nz))/2./dy + &
                         max(v(i,j+1,Nz),0.)*(theta(i,j+1,Nz)+theta(i,j,Nz))/2./dy + &
                         (Dlc_eff(i,j+1,Nz)+Dlc_eff(i,j,Nz))/2./dy*(theta(i,j+1,Nz)+theta(i,j,Nz))/2./dy - &
                         min(w(i,j,Nz),0.)*(theta(i,j,Nz-1)+theta(i,j,Nz))/2./dz + &
                         (Dlc_eff(i,j,Nz-1)+Dlc_eff(i,j,Nz))/2./dz*(theta(i,j,Nz-1)+theta(i,j,Nz))/2./dz + &
                         Dlc_eff(i,j,Nz)/dz*2.*theta(i,j,Nz)/dz) - &
                         dt*theta(i,j,Nz)*(-qc_0*b_o(i,j,Nz)/(kc+lc_o(i,j,Nz))*lo_o(i,j,Nz)/(ko+lo_o(i,j,Nz))) + &
                         dt*theta(i,j,Nz)*alpha_c(i,j,Nz)
                    q(i,j,Nz) = theta(i,j,Nz)*lc_o(i,j,Nz)  + &
                         dt*theta(i,j,Nz)*alpha_c(i,j,Nz)*Ka*sc_o(i,j,Nz) - &
                         dt*(w(i,j,Nz+1)-2.*Dlc_eff(i,j,Nz)/dz)*theta(i,j,Nz)/dz*lc_t
                 endif
                 if(lcbc_t.eq.2) then
                    aa(i,j,Nz) = theta(i,j,Nz) + dt*(-min(u(i,j,Nz),0.)*(theta(i-1,j,Nz)+theta(i,j,Nz))/2./dx + &
                         (Dlc_eff(i-1,j,Nz)+Dlc_eff(i,j,Nz))/2./dx*(theta(i-1,j,Nz)+theta(i,j,Nz))/2./dx + &
                         max(u(i+1,j,Nz),0.)*(theta(i+1,j,Nz)+theta(i,j,Nz))/2./dx + &
                         (Dlc_eff(i+1,j,Nz)+Dlc_eff(i,j,Nz))/2./dx*(theta(i+1,j,Nz)+theta(i,j,Nz))/2./dx - &
                         min(v(i,j,Nz),0.)*(theta(i,j-1,Nz)+theta(i,j,Nz))/2./dy + &
                         (Dlc_eff(i,j-1,Nz)+Dlc_eff(i,j,Nz))/2./dy*(theta(i,j-1,Nz)+theta(i,j,Nz))/2./dy + &
                         max(v(i,j+1,Nz),0.)*(theta(i,j+1,Nz)+theta(i,j,Nz))/2./dy + &
                         (Dlc_eff(i,j+1,Nz)+Dlc_eff(i,j,Nz))/2./dy*(theta(i,j+1,Nz)+theta(i,j,Nz))/2./dy - &
                         min(w(i,j,Nz),0.)*(theta(i,j,Nz-1)+theta(i,j,Nz))/2./dz + &
                         (Dlc_eff(i,j,Nz-1)+Dlc_eff(i,j,Nz))/2./dz*(theta(i,j,Nz-1)+theta(i,j,Nz))/2./dz + &
                         w(i,j,Nz+1)*theta(i,j,Nz)/dz) - &
                         dt*theta(i,j,Nz)*(-qc_0*b_o(i,j,Nz)/(kc+lc_o(i,j,Nz))*lo_o(i,j,Nz)/(ko+lo_o(i,j,Nz))) + &
                         dt*theta(i,j,Nz)*alpha_c(i,j,Nz)
                    q(i,j,Nz) = theta(i,j,Nz)*lc_o(i,j,Nz)  + &
                         dt*theta(i,j,Nz)*alpha_c(i,j,Nz)*Ka*sc_o(i,j,Nz)
                 endif
              endif
           enddo
        enddo



        !! for inside

        Do i=2,Nx-1
           Do j=2,Ny-1
              Do k=2,Nz-1
                 if(wall(i,j,k)) then
                    aa(i,j,k) = 0.
                    bb(i,j,k) = 0.
                    cc(i,j,k) = 0.
                    dd(i,j,k) = 0.
                    ee(i,j,k) = 0.
                    ff(i,j,k) = 0.
                    gg(i,j,k) = 0.
                    q(i,j,k) = 0.
                 else
                    aa(i,j,k) = theta(i,j,k) + dt*(-min(u(i,j,k),0.)*(theta(i-1,j,k)+theta(i,j,k))/2./dx + &
                         (Dlc_eff(i-1,j,k)+Dlc_eff(i,j,k))/2./dx*(theta(i-1,j,k)+theta(i,j,k))/2./dx + &
                         max(u(i+1,j,k),0.)*(theta(i+1,j,k)+theta(i,j,k))/2./dx + &
                         (Dlc_eff(i+1,j,k)+Dlc_eff(i,j,k))/2./dx*(theta(i+1,j,k)+theta(i,j,k))/2./dx - &
                         min(v(i,j,k),0.)*(theta(i,j-1,k)+theta(i,j,k))/2./dy + &
                         (Dlc_eff(i,j-1,k)+Dlc_eff(i,j,k))/2./dy*(theta(i,j-1,k)+theta(i,j,k))/2./dy + &
                         max(v(i,j+1,k),0.)*(theta(i,j+1,k)+theta(i,j,k))/2./dy + &
                         (Dlc_eff(i,j+1,k)+Dlc_eff(i,j,k))/2./dy*(theta(i,j+1,k)+theta(i,j,k))/2./dy - &
                         min(w(i,j,k),0.)*(theta(i,j,k-1)+theta(i,j,k))/2./dz + &
                         (Dlc_eff(i,j,k-1)+Dlc_eff(i,j,k))/2./dz*(theta(i,j,k-1)+theta(i,j,k))/2./dz + &
                         max(w(i,j,k+1),0.)*(theta(i,j,k+1)+theta(i,j,k))/2./dz + &
                         (Dlc_eff(i,j,k+1)+Dlc_eff(i,j,k))/2./dz*(theta(i,j,k+1)+theta(i,j,k))/2./dz) - &
                         dt*theta(i,j,k)*(-qc_0*b_o(i,j,k)/(Kc+lc_o(i,j,k))*lo_o(i,j,k)/(Ko+lo_o(i,j,k))) + &
                         dt*theta(i,j,k)*alpha_c(i,j,k)

                    bb(i,j,k) = -dt*(max(u(i,j,k),0.)+(Dlc_eff(i-1,j,k)+Dlc_eff(i,j,k))/2./dx) * &
                         (theta(i-1,j,k)+theta(i,j,k))/2./dx
                    cc(i,j,k) = dt*(min(u(i+1,j,k),0.)-(Dlc_eff(i+1,j,k)+Dlc_eff(i,j,k))/2./dx) * &
                         (theta(i+1,j,k)+theta(i,j,k))/2./dx
                    dd(i,j,k) = -dt*(max(v(i,j,k),0.)+(Dlc_eff(i,j-1,k)+Dlc_eff(i,j,k))/2./dy) * &
                         (theta(i,j-1,k)+theta(i,j,k))/2./dy
                    ee(i,j,k) = dt*(min(v(i,j+1,k),0.)-(Dlc_eff(i,j+1,k)+Dlc_eff(i,j,k))/2./dy) * &
                         (theta(i,j+1,k)+theta(i,j,k))/2./dy
                    ff(i,j,k) = -dt*(max(w(i,j,k),0.)+(Dlc_eff(i,j,k-1)+Dlc_eff(i,j,k))/2./dz) * &
                         (theta(i,j,k-1)+theta(i,j,k))/2./dz
                    gg(i,j,k) = dt*(min(w(i,j,k+1),0.)-(Dlc_eff(i,j,k+1)+Dlc_eff(i,j,k))/2./dz) * &
                         (theta(i,j,k+1)+theta(i,j,k))/2./dz

                    q(i,j,k) = theta(i,j,k)*lc_o(i,j,k)  + &
                         dt*theta(i,j,k)*alpha_c(i,j,k)*Ka*sc_o(i,j,k)
                 endif
              Enddo
           Enddo
        Enddo


     else
        if(mype.eq.0) then

           Do i=2,Nx-1
              Do j=2,Ny-1
                 if(wall(i,j,1)) then
                    aa(i,j,1) = 0.
                    bb(i,j,1) = 0.
                    cc(i,j,1) = 0.
                    dd(i,j,1) = 0.
                    ee(i,j,1) = 0.
                    ff(i,j,1) = 0.
                    gg(i,j,1) = 0.
                    q(i,j,1) = 0.
                 else
                    bb(i,j,1) = -dt*(max(u(i,j,1),0.)+(Dlc_eff(i-1,j,1)+Dlc_eff(i,j,1))/2./dx) * &
                         (theta(i-1,j,1)+theta(i,j,1))/2./dx
                    cc(i,j,1) = dt*(min(u(i+1,j,1),0.)-(Dlc_eff(i+1,j,1)+Dlc_eff(i,j,1))/2./dx) * &
                         (theta(i+1,j,1)+theta(i,j,1))/2./dx
                    dd(i,j,1) = -dt*(max(v(i,j,1),0.)+(Dlc_eff(i,j-1,1)+Dlc_eff(i,j,1))/2./dy) * &
                         (theta(i,j-1,1)+theta(i,j,1))/2./dy
                    ee(i,j,1) = dt*(min(v(i,j+1,1),0.)-(Dlc_eff(i,j+1,1)+Dlc_eff(i,j,1))/2./dy) * &
                         (theta(i,j+1,1)+theta(i,j,1))/2./dy
                    ff(i,j,1) = 0.
                    gg(i,j,1) = dt*(min(w(i,j,1+1),0.)-(Dlc_eff(i,j,1+1)+Dlc_eff(i,j,1))/2./dz) * &
                         (theta(i,j,1+1)+theta(i,j,1))/2./dz

                    if(lcbc_d.eq.1) then 
                       aa(i,j,1) = theta(i,j,1) + dt*(-min(u(i,j,1),0.)*(theta(i-1,j,1)+theta(i,j,1))/2./dx + &
                            (Dlc_eff(i-1,j,1)+Dlc_eff(i,j,1))/2./dx*(theta(i-1,j,1)+theta(i,j,1))/2./dx + &
                            max(u(i+1,j,1),0.)*(theta(i+1,j,1)+theta(i,j,1))/2./dx + &
                            (Dlc_eff(i+1,j,1)+Dlc_eff(i,j,1))/2./dx*(theta(i+1,j,1)+theta(i,j,1))/2./dx - &
                            min(v(i,j,1),0.)*(theta(i,j-1,1)+theta(i,j,1))/2./dy + &
                            (Dlc_eff(i,j-1,1)+Dlc_eff(i,j,1))/2./dy*(theta(i,j-1,1)+theta(i,j,1))/2./dy + &
                            max(v(i,j+1,1),0.)*(theta(i,j+1,1)+theta(i,j,1))/2./dy + &
                            (Dlc_eff(i,j+1,1)+Dlc_eff(i,j,1))/2./dy*(theta(i,j+1,1)+theta(i,j,1))/2./dy + &
                            Dlc_eff(i,j,1)/dz*2.*theta(i,j,1)/dz + &
                            max(w(i,j,1+1),0.)*(theta(i,j,1+1)+theta(i,j,1))/2./dz + &
                            (Dlc_eff(i,j,1+1)+Dlc_eff(i,j,1))/2./dz*(theta(i,j,1+1)+theta(i,j,1))/2./dz) - &
                            dt*theta(i,j,1)*(-qc_0*b_o(i,j,1)/(Kc+lc_o(i,j,1))*lo_o(i,j,1)/(Ko+lo_o(i,j,1))) + &
                            dt*theta(i,j,1)*alpha_c(i,j,1)
                       q(i,j,1) = theta(i,j,1)*lc_o(i,j,1) + &
                            dt*theta(i,j,1)*alpha_c(i,j,1)*Ka*sc_o(i,j,1) + &
                            dt*(w(i,j,1)+2.*Dlc_eff(i,j,1)/dz)*theta(i,j,1)/dz*lc_d
                    endif
                    if(lcbc_d.eq.2) then 

                       aa(i,j,1) = theta(i,j,1) + dt*(-min(u(i,j,1),0.)*(theta(i-1,j,1)+theta(i,j,1))/2./dx + &
                            (Dlc_eff(i-1,j,1)+Dlc_eff(i,j,1))/2./dx*(theta(i-1,j,1)+theta(i,j,1))/2./dx + &
                            max(u(i+1,j,1),0.)*(theta(i+1,j,1)+theta(i,j,1))/2./dx + &
                            (Dlc_eff(i+1,j,1)+Dlc_eff(i,j,1))/2./dx*(theta(i+1,j,1)+theta(i,j,1))/2./dx - &
                            min(v(i,j,1),0.)*(theta(i,j-1,1)+theta(i,j,1))/2./dy + &
                            (Dlc_eff(i,j-1,1)+Dlc_eff(i,j,1))/2./dy*(theta(i,j-1,1)+theta(i,j,1))/2./dy + &
                            max(v(i,j+1,1),0.)*(theta(i,j+1,1)+theta(i,j,1))/2./dy + &
                            (Dlc_eff(i,j+1,1)+Dlc_eff(i,j,1))/2./dy*(theta(i,j+1,1)+theta(i,j,1))/2./dy - &
                            w(i,j,1)*theta(i,j,1)/dz + &
                            max(w(i,j,1+1),0.)*(theta(i,j,1+1)+theta(i,j,1))/2./dz + &
                            (Dlc_eff(i,j,1+1)+Dlc_eff(i,j,1))/2./dz*(theta(i,j,1+1)+theta(i,j,1))/2./dz) - &
                            dt*theta(i,j,1)*(-qc_0*b_o(i,j,1)/(Kc+lc_o(i,j,1))*lo_o(i,j,1)/(Ko+lo_o(i,j,1))) + &
                            dt*theta(i,j,1)*alpha_c(i,j,1)
                       q(i,j,1) = theta(i,j,1)*lc_o(i,j,1)  + &
                            dt*theta(i,j,1)*alpha_c(i,j,1)*Ka*sc_o(i,j,1)
                    endif

                 endif
              enddo
           enddo


           Do i=2,Nx-1
              Do j=2,Ny-1
                 if(wall(i,j,Nz)) then
                    aa(i,j,Nz) = 0.
                    bb(i,j,Nz) = 0.
                    cc(i,j,Nz) = 0.
                    dd(i,j,Nz) = 0.
                    ee(i,j,Nz) = 0.
                    ff(i,j,Nz) = 0.
                    gg(i,j,Nz) = 0.
                    q(i,j,Nz) = 0.
                 else
                    aa(i,j,Nz) = theta(i,j,Nz) + dt*(-min(u(i,j,Nz),0.)*(theta(i-1,j,Nz)+theta(i,j,Nz))/2./dx + &
                         (Dlc_eff(i-1,j,Nz)+Dlc_eff(i,j,Nz))/2./dx*(theta(i-1,j,Nz)+theta(i,j,Nz))/2./dx + &
                         max(u(i+1,j,Nz),0.)*(theta(i+1,j,Nz)+theta(i,j,Nz))/2./dx + &
                         (Dlc_eff(i+1,j,Nz)+Dlc_eff(i,j,Nz))/2./dx*(theta(i+1,j,Nz)+theta(i,j,Nz))/2./dx - &
                         min(v(i,j,Nz),0.)*(theta(i,j-1,Nz)+theta(i,j,Nz))/2./dy + &
                         (Dlc_eff(i,j-1,Nz)+Dlc_eff(i,j,Nz))/2./dy*(theta(i,j-1,Nz)+theta(i,j,Nz))/2./dy + &
                         max(v(i,j+1,Nz),0.)*(theta(i,j+1,Nz)+theta(i,j,Nz))/2./dy + &
                         (Dlc_eff(i,j+1,Nz)+Dlc_eff(i,j,Nz))/2./dy*(theta(i,j+1,Nz)+theta(i,j,Nz))/2./dy - &
                         min(w(i,j,Nz),0.)*(theta(i,j,Nz-1)+theta(i,j,Nz))/2./dz + &
                         (Dlc_eff(i,j,Nz-1)+Dlc_eff(i,j,Nz))/2./dz*(theta(i,j,Nz-1)+theta(i,j,Nz))/2./dz + &
                         max(w(i,j,Nz+1),0.)*(tmp_theta(i,j,2)+theta(i,j,Nz))/2./dz + &
                         (tmp_Dlc_eff(i,j,2)+Dlc_eff(i,j,Nz))/2./dz*(tmp_theta(i,j,2)+theta(i,j,Nz))/2./dz) - &
                         dt*theta(i,j,Nz)*(-qc_0*b_o(i,j,Nz)/(Kc+lc_o(i,j,Nz))*lo_o(i,j,Nz)/(Ko+lo_o(i,j,Nz))) + &
                         dt*theta(i,j,Nz)*alpha_c(i,j,Nz)

                    bb(i,j,Nz) = -dt*(max(u(i,j,Nz),0.)+(Dlc_eff(i-1,j,Nz)+Dlc_eff(i,j,Nz))/2./dx) * &
                         (theta(i-1,j,Nz)+theta(i,j,Nz))/2./dx
                    cc(i,j,Nz) = dt*(min(u(i+1,j,Nz),0.)-(Dlc_eff(i+1,j,Nz)+Dlc_eff(i,j,Nz))/2./dx) * &
                         (theta(i+1,j,Nz)+theta(i,j,Nz))/2./dx
                    dd(i,j,Nz) = -dt*(max(v(i,j,Nz),0.)+(Dlc_eff(i,j-1,Nz)+Dlc_eff(i,j,Nz))/2./dy) * &
                         (theta(i,j-1,Nz)+theta(i,j,Nz))/2./dy
                    ee(i,j,Nz) = dt*(min(v(i,j+1,Nz),0.)-(Dlc_eff(i,j+1,Nz)+Dlc_eff(i,j,Nz))/2./dy) * &
                         (theta(i,j+1,Nz)+theta(i,j,Nz))/2./dy
                    ff(i,j,Nz) = -dt*(max(w(i,j,Nz),0.)+(Dlc_eff(i,j,Nz-1)+Dlc_eff(i,j,Nz))/2./dz) * &
                         (theta(i,j,Nz-1)+theta(i,j,Nz))/2./dz
                    gg(i,j,Nz) = 0.

                    q(i,j,Nz) = theta(i,j,Nz)*lc_o(i,j,Nz)  + &
                         dt*theta(i,j,Nz)*alpha_c(i,j,Nz)*Ka*sc_o(i,j,Nz) - &
                         (dt*(min(w(i,j,Nz+1),0.)-(tmp_Dlc_eff(i,j,2)+Dlc_eff(i,j,Nz))/2./dz) * &
                         (tmp_theta(i,j,2)+theta(i,j,Nz))/2./dz)*tmp_lc_o(i,j,2)
                 endif
              Enddo
           Enddo



           !! for inside

           Do i=2,Nx-1
              Do j=2,Ny-1
                 Do k=2,Nz-1
                    if(wall(i,j,k)) then
                       aa(i,j,k) = 0.
                       bb(i,j,k) = 0.
                       cc(i,j,k) = 0.
                       dd(i,j,k) = 0.
                       ee(i,j,k) = 0.
                       ff(i,j,k) = 0.
                       gg(i,j,k) = 0.
                       q(i,j,k) = 0.
                    else
                       aa(i,j,k) = theta(i,j,k) + dt*(-min(u(i,j,k),0.)*(theta(i-1,j,k)+theta(i,j,k))/2./dx + &
                            (Dlc_eff(i-1,j,k)+Dlc_eff(i,j,k))/2./dx*(theta(i-1,j,k)+theta(i,j,k))/2./dx + &
                            max(u(i+1,j,k),0.)*(theta(i+1,j,k)+theta(i,j,k))/2./dx + &
                            (Dlc_eff(i+1,j,k)+Dlc_eff(i,j,k))/2./dx*(theta(i+1,j,k)+theta(i,j,k))/2./dx - &
                            min(v(i,j,k),0.)*(theta(i,j-1,k)+theta(i,j,k))/2./dy + &
                            (Dlc_eff(i,j-1,k)+Dlc_eff(i,j,k))/2./dy*(theta(i,j-1,k)+theta(i,j,k))/2./dy + &
                            max(v(i,j+1,k),0.)*(theta(i,j+1,k)+theta(i,j,k))/2./dy + &
                            (Dlc_eff(i,j+1,k)+Dlc_eff(i,j,k))/2./dy*(theta(i,j+1,k)+theta(i,j,k))/2./dy - &
                            min(w(i,j,k),0.)*(theta(i,j,k-1)+theta(i,j,k))/2./dz + &
                            (Dlc_eff(i,j,k-1)+Dlc_eff(i,j,k))/2./dz*(theta(i,j,k-1)+theta(i,j,k))/2./dz + &
                            max(w(i,j,k+1),0.)*(theta(i,j,k+1)+theta(i,j,k))/2./dz + &
                            (Dlc_eff(i,j,k+1)+Dlc_eff(i,j,k))/2./dz*(theta(i,j,k+1)+theta(i,j,k))/2./dz) - &
                            dt*theta(i,j,k)*(-qc_0*b_o(i,j,k)/(Kc+lc_o(i,j,k))*lo_o(i,j,k)/(Ko+lo_o(i,j,k))) + &
                            dt*theta(i,j,k)*alpha_c(i,j,k)

                       bb(i,j,k) = -dt*(max(u(i,j,k),0.)+(Dlc_eff(i-1,j,k)+Dlc_eff(i,j,k))/2./dx) * &
                            (theta(i-1,j,k)+theta(i,j,k))/2./dx
                       cc(i,j,k) = dt*(min(u(i+1,j,k),0.)-(Dlc_eff(i+1,j,k)+Dlc_eff(i,j,k))/2./dx) * &
                            (theta(i+1,j,k)+theta(i,j,k))/2./dx
                       dd(i,j,k) = -dt*(max(v(i,j,k),0.)+(Dlc_eff(i,j-1,k)+Dlc_eff(i,j,k))/2./dy) * &
                            (theta(i,j-1,k)+theta(i,j,k))/2./dy
                       ee(i,j,k) = dt*(min(v(i,j+1,k),0.)-(Dlc_eff(i,j+1,k)+Dlc_eff(i,j,k))/2./dy) * &
                            (theta(i,j+1,k)+theta(i,j,k))/2./dy
                       ff(i,j,k) = -dt*(max(w(i,j,k),0.)+(Dlc_eff(i,j,k-1)+Dlc_eff(i,j,k))/2./dz) * &
                            (theta(i,j,k-1)+theta(i,j,k))/2./dz
                       gg(i,j,k) = dt*(min(w(i,j,k+1),0.)-(Dlc_eff(i,j,k+1)+Dlc_eff(i,j,k))/2./dz) * &
                            (theta(i,j,k+1)+theta(i,j,k))/2./dz

                       q(i,j,k) = theta(i,j,k)*lc_o(i,j,k)  + &
                            dt*theta(i,j,k)*alpha_c(i,j,k)*Ka*sc_o(i,j,k)
                    endif
                 Enddo
              Enddo
           Enddo



        else if(mype.eq.npes-1) then

           Do i=2,Nx-1
              Do j=2,Ny-1
                 if(wall(i,j,1)) then
                    aa(i,j,1) = 0.
                    bb(i,j,1) = 0.
                    cc(i,j,1) = 0.
                    dd(i,j,1) = 0.
                    ee(i,j,1) = 0.
                    ff(i,j,1) = 0.
                    gg(i,j,1) = 0.
                    q(i,j,1) = 0.
                 else
                    aa(i,j,1) = theta(i,j,1) + dt*(-min(u(i,j,1),0.)*(theta(i-1,j,1)+theta(i,j,1))/2./dx + &
                         (Dlc_eff(i-1,j,1)+Dlc_eff(i,j,1))/2./dx*(theta(i-1,j,1)+theta(i,j,1))/2./dx + &
                         max(u(i+1,j,1),0.)*(theta(i+1,j,1)+theta(i,j,1))/2./dx + &
                         (Dlc_eff(i+1,j,1)+Dlc_eff(i,j,1))/2./dx*(theta(i+1,j,1)+theta(i,j,1))/2./dx - &
                         min(v(i,j,1),0.)*(theta(i,j-1,1)+theta(i,j,1))/2./dy + &
                         (Dlc_eff(i,j-1,1)+Dlc_eff(i,j,1))/2./dy*(theta(i,j-1,1)+theta(i,j,1))/2./dy + &
                         max(v(i,j+1,1),0.)*(theta(i,j+1,1)+theta(i,j,1))/2./dy + &
                         (Dlc_eff(i,j+1,1)+Dlc_eff(i,j,1))/2./dy*(theta(i,j+1,1)+theta(i,j,1))/2./dy - &
                         min(w(i,j,1),0.)*(tmp_theta(i,j,1)+theta(i,j,1))/2./dz + &
                         (tmp_Dlc_eff(i,j,1)+Dlc_eff(i,j,1))/2./dz*(tmp_theta(i,j,1)+theta(i,j,1))/2./dz + &
                         max(w(i,j,1+1),0.)*(theta(i,j,1+1)+theta(i,j,1))/2./dz + &
                         (Dlc_eff(i,j,1+1)+Dlc_eff(i,j,1))/2./dz*(theta(i,j,1+1)+theta(i,j,1))/2./dz) - &
                         dt*theta(i,j,1)*(-qc_0*b_o(i,j,1)/(Kc+lc_o(i,j,1))*lo_o(i,j,1)/(Ko+lo_o(i,j,1))) + &
                         dt*theta(i,j,1)*alpha_c(i,j,1)

                    bb(i,j,1) = -dt*(max(u(i,j,1),0.)+(Dlc_eff(i-1,j,1)+Dlc_eff(i,j,1))/2./dx) * &
                         (theta(i-1,j,1)+theta(i,j,1))/2./dx
                    cc(i,j,1) = dt*(min(u(i+1,j,1),0.)-(Dlc_eff(i+1,j,1)+Dlc_eff(i,j,1))/2./dx) * &
                         (theta(i+1,j,1)+theta(i,j,1))/2./dx
                    dd(i,j,1) = -dt*(max(v(i,j,1),0.)+(Dlc_eff(i,j-1,1)+Dlc_eff(i,j,1))/2./dy) * &
                         (theta(i,j-1,1)+theta(i,j,1))/2./dy
                    ee(i,j,1) = dt*(min(v(i,j+1,1),0.)-(Dlc_eff(i,j+1,1)+Dlc_eff(i,j,1))/2./dy) * &
                         (theta(i,j+1,1)+theta(i,j,1))/2./dy
                    ff(i,j,1) = 0.
                    gg(i,j,1) = dt*(min(w(i,j,1+1),0.)-(Dlc_eff(i,j,1+1)+Dlc_eff(i,j,1))/2./dz) * &
                         (theta(i,j,1+1)+theta(i,j,1))/2./dz

                    q(i,j,1) = theta(i,j,1)*lc_o(i,j,1)  + &
                         dt*theta(i,j,1)*alpha_c(i,j,1)*Ka*sc_o(i,j,1) - &
                         (-dt*(max(w(i,j,1),0.)+(tmp_Dlc_eff(i,j,1)+Dlc_eff(i,j,1))/2./dz) * &
                         (tmp_theta(i,j,1)+theta(i,j,1))/2./dz)*tmp_lc_o(i,j,1)
                 endif
              Enddo
           Enddo


           Do i=2,Nx-1
              Do j=2,Ny-1
                 if(wall(i,j,Nz)) then
                    aa(i,j,Nz) = 0.
                    bb(i,j,Nz) = 0.
                    cc(i,j,Nz) = 0.
                    dd(i,j,Nz) = 0.
                    ee(i,j,Nz) = 0.
                    ff(i,j,Nz) = 0.
                    gg(i,j,Nz) = 0.
                    q(i,j,Nz) = 0.
                 else     
                    bb(i,j,Nz) = -dt*(max(u(i,j,Nz),0.)+(Dlc_eff(i-1,j,Nz)+Dlc_eff(i,j,Nz))/2./dx) * &
                         (theta(i-1,j,Nz)+theta(i,j,Nz))/2./dx
                    cc(i,j,Nz) = dt*(min(u(i+1,j,Nz),0.)-(Dlc_eff(i+1,j,Nz)+Dlc_eff(i,j,Nz))/2./dx) * &
                         (theta(i+1,j,Nz)+theta(i,j,Nz))/2./dx
                    dd(i,j,Nz) = -dt*(max(v(i,j,Nz),0.)+(Dlc_eff(i,j-1,Nz)+Dlc_eff(i,j,Nz))/2./dy) * &
                         (theta(i,j-1,Nz)+theta(i,j,Nz))/2./dy
                    ee(i,j,Nz) = dt*(min(v(i,j+1,Nz),0.)-(Dlc_eff(i,j+1,Nz)+Dlc_eff(i,j,Nz))/2./dy) * &
                         (theta(i,j+1,Nz)+theta(i,j,Nz))/2./dy
                    ff(i,j,Nz) = -dt*(max(w(i,j,Nz),0.)+(Dlc_eff(i,j,Nz-1)+Dlc_eff(i,j,Nz))/2./dz) * &
                         (theta(i,j,Nz-1)+theta(i,j,Nz))/2./dz
                    gg(i,j,Nz) = 0.

                    if(lcbc_t.eq.1) then
                       aa(i,j,Nz) = theta(i,j,Nz) + dt*(-min(u(i,j,Nz),0.)*(theta(i-1,j,Nz)+theta(i,j,Nz))/2./dx + &
                            (Dlc_eff(i-1,j,Nz)+Dlc_eff(i,j,Nz))/2./dx*(theta(i-1,j,Nz)+theta(i,j,Nz))/2./dx + &
                            max(u(i+1,j,Nz),0.)*(theta(i+1,j,Nz)+theta(i,j,Nz))/2./dx + &
                            (Dlc_eff(i+1,j,Nz)+Dlc_eff(i,j,Nz))/2./dx*(theta(i+1,j,Nz)+theta(i,j,Nz))/2./dx - &
                            min(v(i,j,Nz),0.)*(theta(i,j-1,Nz)+theta(i,j,Nz))/2./dy + &
                            (Dlc_eff(i,j-1,Nz)+Dlc_eff(i,j,Nz))/2./dy*(theta(i,j-1,Nz)+theta(i,j,Nz))/2./dy + &
                            max(v(i,j+1,Nz),0.)*(theta(i,j+1,Nz)+theta(i,j,Nz))/2./dy + &
                            (Dlc_eff(i,j+1,Nz)+Dlc_eff(i,j,Nz))/2./dy*(theta(i,j+1,Nz)+theta(i,j,Nz))/2./dy - &
                            min(w(i,j,Nz),0.)*(theta(i,j,Nz-1)+theta(i,j,Nz))/2./dz + &
                            (Dlc_eff(i,j,Nz-1)+Dlc_eff(i,j,Nz))/2./dz*(theta(i,j,Nz-1)+theta(i,j,Nz))/2./dz + &
                            Dlc_eff(i,j,Nz)/dz*2.*theta(i,j,Nz)/dz) - &
                            dt*theta(i,j,Nz)*(-qc_0*b_o(i,j,Nz)/(kc+lc_o(i,j,Nz))*lo_o(i,j,Nz)/(ko+lo_o(i,j,Nz))) + &
                            dt*theta(i,j,Nz)*alpha_c(i,j,Nz)
                       q(i,j,Nz) = theta(i,j,Nz)*lc_o(i,j,Nz)  + &
                            dt*theta(i,j,Nz)*alpha_c(i,j,Nz)*Ka*sc_o(i,j,Nz) - &
                            dt*(w(i,j,Nz+1)-2.*Dlc_eff(i,j,Nz)/dz)*theta(i,j,Nz)/dz*lc_t
                    endif
                    if(lcbc_t.eq.2) then
                       aa(i,j,Nz) = theta(i,j,Nz) + dt*(-min(u(i,j,Nz),0.)*(theta(i-1,j,Nz)+theta(i,j,Nz))/2./dx + &
                            (Dlc_eff(i-1,j,Nz)+Dlc_eff(i,j,Nz))/2./dx*(theta(i-1,j,Nz)+theta(i,j,Nz))/2./dx + &
                            max(u(i+1,j,Nz),0.)*(theta(i+1,j,Nz)+theta(i,j,Nz))/2./dx + &
                            (Dlc_eff(i+1,j,Nz)+Dlc_eff(i,j,Nz))/2./dx*(theta(i+1,j,Nz)+theta(i,j,Nz))/2./dx - &
                            min(v(i,j,Nz),0.)*(theta(i,j-1,Nz)+theta(i,j,Nz))/2./dy + &
                            (Dlc_eff(i,j-1,Nz)+Dlc_eff(i,j,Nz))/2./dy*(theta(i,j-1,Nz)+theta(i,j,Nz))/2./dy + &
                            max(v(i,j+1,Nz),0.)*(theta(i,j+1,Nz)+theta(i,j,Nz))/2./dy + &
                            (Dlc_eff(i,j+1,Nz)+Dlc_eff(i,j,Nz))/2./dy*(theta(i,j+1,Nz)+theta(i,j,Nz))/2./dy - &
                            min(w(i,j,Nz),0.)*(theta(i,j,Nz-1)+theta(i,j,Nz))/2./dz + &
                            (Dlc_eff(i,j,Nz-1)+Dlc_eff(i,j,Nz))/2./dz*(theta(i,j,Nz-1)+theta(i,j,Nz))/2./dz + &
                            w(i,j,Nz+1)*theta(i,j,Nz)/dz) - &
                            dt*theta(i,j,Nz)*(-qc_0*b_o(i,j,Nz)/(kc+lc_o(i,j,Nz))*lo_o(i,j,Nz)/(ko+lo_o(i,j,Nz))) + &
                            dt*theta(i,j,Nz)*alpha_c(i,j,Nz)
                       q(i,j,Nz) = theta(i,j,Nz)*lc_o(i,j,Nz)  + &
                            dt*theta(i,j,Nz)*alpha_c(i,j,Nz)*Ka*sc_o(i,j,Nz)
                    endif
                 endif
              enddo
           enddo


           !! for inside

           Do i=2,Nx-1
              Do j=2,Ny-1
                 Do k=2,Nz-1
                    if(wall(i,j,k)) then
                       aa(i,j,k) = 0.
                       bb(i,j,k) = 0.
                       cc(i,j,k) = 0.
                       dd(i,j,k) = 0.
                       ee(i,j,k) = 0.
                       ff(i,j,k) = 0.
                       gg(i,j,k) = 0.
                       q(i,j,k) = 0.
                    else
                       aa(i,j,k) = theta(i,j,k) + dt*(-min(u(i,j,k),0.)*(theta(i-1,j,k)+theta(i,j,k))/2./dx + &
                            (Dlc_eff(i-1,j,k)+Dlc_eff(i,j,k))/2./dx*(theta(i-1,j,k)+theta(i,j,k))/2./dx + &
                            max(u(i+1,j,k),0.)*(theta(i+1,j,k)+theta(i,j,k))/2./dx + &
                            (Dlc_eff(i+1,j,k)+Dlc_eff(i,j,k))/2./dx*(theta(i+1,j,k)+theta(i,j,k))/2./dx - &
                            min(v(i,j,k),0.)*(theta(i,j-1,k)+theta(i,j,k))/2./dy + &
                            (Dlc_eff(i,j-1,k)+Dlc_eff(i,j,k))/2./dy*(theta(i,j-1,k)+theta(i,j,k))/2./dy + &
                            max(v(i,j+1,k),0.)*(theta(i,j+1,k)+theta(i,j,k))/2./dy + &
                            (Dlc_eff(i,j+1,k)+Dlc_eff(i,j,k))/2./dy*(theta(i,j+1,k)+theta(i,j,k))/2./dy - &
                            min(w(i,j,k),0.)*(theta(i,j,k-1)+theta(i,j,k))/2./dz + &
                            (Dlc_eff(i,j,k-1)+Dlc_eff(i,j,k))/2./dz*(theta(i,j,k-1)+theta(i,j,k))/2./dz + &
                            max(w(i,j,k+1),0.)*(theta(i,j,k+1)+theta(i,j,k))/2./dz + &
                            (Dlc_eff(i,j,k+1)+Dlc_eff(i,j,k))/2./dz*(theta(i,j,k+1)+theta(i,j,k))/2./dz) - &
                            dt*theta(i,j,k)*(-qc_0*b_o(i,j,k)/(Kc+lc_o(i,j,k))*lo_o(i,j,k)/(Ko+lo_o(i,j,k))) + &
                            dt*theta(i,j,k)*alpha_c(i,j,k)

                       bb(i,j,k) = -dt*(max(u(i,j,k),0.)+(Dlc_eff(i-1,j,k)+Dlc_eff(i,j,k))/2./dx) * &
                            (theta(i-1,j,k)+theta(i,j,k))/2./dx
                       cc(i,j,k) = dt*(min(u(i+1,j,k),0.)-(Dlc_eff(i+1,j,k)+Dlc_eff(i,j,k))/2./dx) * &
                            (theta(i+1,j,k)+theta(i,j,k))/2./dx
                       dd(i,j,k) = -dt*(max(v(i,j,k),0.)+(Dlc_eff(i,j-1,k)+Dlc_eff(i,j,k))/2./dy) * &
                            (theta(i,j-1,k)+theta(i,j,k))/2./dy
                       ee(i,j,k) = dt*(min(v(i,j+1,k),0.)-(Dlc_eff(i,j+1,k)+Dlc_eff(i,j,k))/2./dy) * &
                            (theta(i,j+1,k)+theta(i,j,k))/2./dy
                       ff(i,j,k) = -dt*(max(w(i,j,k),0.)+(Dlc_eff(i,j,k-1)+Dlc_eff(i,j,k))/2./dz) * &
                            (theta(i,j,k-1)+theta(i,j,k))/2./dz
                       gg(i,j,k) = dt*(min(w(i,j,k+1),0.)-(Dlc_eff(i,j,k+1)+Dlc_eff(i,j,k))/2./dz) * &
                            (theta(i,j,k+1)+theta(i,j,k))/2./dz

                       q(i,j,k) = theta(i,j,k)*lc_o(i,j,k)  + &
                            dt*theta(i,j,k)*alpha_c(i,j,k)*Ka*sc_o(i,j,k)
                    endif
                 Enddo
              Enddo
           Enddo


        else 

           Do i=2,Nx-1
              Do j=2,Ny-1
                 if(wall(i,j,1)) then
                    aa(i,j,1) = 0.
                    bb(i,j,1) = 0.
                    cc(i,j,1) = 0.
                    dd(i,j,1) = 0.
                    ee(i,j,1) = 0.
                    ff(i,j,1) = 0.
                    gg(i,j,1) = 0.
                    q(i,j,1) = 0.
                 else
                    aa(i,j,1) = theta(i,j,1) + dt*(-min(u(i,j,1),0.)*(theta(i-1,j,1)+theta(i,j,1))/2./dx + &
                         (Dlc_eff(i-1,j,1)+Dlc_eff(i,j,1))/2./dx*(theta(i-1,j,1)+theta(i,j,1))/2./dx + &
                         max(u(i+1,j,1),0.)*(theta(i+1,j,1)+theta(i,j,1))/2./dx + &
                         (Dlc_eff(i+1,j,1)+Dlc_eff(i,j,1))/2./dx*(theta(i+1,j,1)+theta(i,j,1))/2./dx - &
                         min(v(i,j,1),0.)*(theta(i,j-1,1)+theta(i,j,1))/2./dy + &
                         (Dlc_eff(i,j-1,1)+Dlc_eff(i,j,1))/2./dy*(theta(i,j-1,1)+theta(i,j,1))/2./dy + &
                         max(v(i,j+1,1),0.)*(theta(i,j+1,1)+theta(i,j,1))/2./dy + &
                         (Dlc_eff(i,j+1,1)+Dlc_eff(i,j,1))/2./dy*(theta(i,j+1,1)+theta(i,j,1))/2./dy - &
                         min(w(i,j,1),0.)*(tmp_theta(i,j,1)+theta(i,j,1))/2./dz + &
                         (tmp_Dlc_eff(i,j,1)+Dlc_eff(i,j,1))/2./dz*(tmp_theta(i,j,1)+theta(i,j,1))/2./dz + &
                         max(w(i,j,1+1),0.)*(theta(i,j,1+1)+theta(i,j,1))/2./dz + &
                         (Dlc_eff(i,j,1+1)+Dlc_eff(i,j,1))/2./dz*(theta(i,j,1+1)+theta(i,j,1))/2./dz) - &
                         dt*theta(i,j,1)*(-qc_0*b_o(i,j,1)/(Kc+lc_o(i,j,1))*lo_o(i,j,1)/(Ko+lo_o(i,j,1))) + &
                         dt*theta(i,j,1)*alpha_c(i,j,1)

                    bb(i,j,1) = -dt*(max(u(i,j,1),0.)+(Dlc_eff(i-1,j,1)+Dlc_eff(i,j,1))/2./dx) * &
                         (theta(i-1,j,1)+theta(i,j,1))/2./dx
                    cc(i,j,1) = dt*(min(u(i+1,j,1),0.)-(Dlc_eff(i+1,j,1)+Dlc_eff(i,j,1))/2./dx) * &
                         (theta(i+1,j,1)+theta(i,j,1))/2./dx
                    dd(i,j,1) = -dt*(max(v(i,j,1),0.)+(Dlc_eff(i,j-1,1)+Dlc_eff(i,j,1))/2./dy) * &
                         (theta(i,j-1,1)+theta(i,j,1))/2./dy
                    ee(i,j,1) = dt*(min(v(i,j+1,1),0.)-(Dlc_eff(i,j+1,1)+Dlc_eff(i,j,1))/2./dy) * &
                         (theta(i,j+1,1)+theta(i,j,1))/2./dy
                    ff(i,j,1) = 0.
                    gg(i,j,1) = dt*(min(w(i,j,1+1),0.)-(Dlc_eff(i,j,1+1)+Dlc_eff(i,j,1))/2./dz) * &
                         (theta(i,j,1+1)+theta(i,j,1))/2./dz

                    q(i,j,1) = theta(i,j,1)*lc_o(i,j,1)  + &
                         dt*theta(i,j,1)*alpha_c(i,j,1)*Ka*sc_o(i,j,1) - &
                         (-dt*(max(w(i,j,1),0.)+(tmp_Dlc_eff(i,j,1)+Dlc_eff(i,j,1))/2./dz) * &
                         (tmp_theta(i,j,1)+theta(i,j,1))/2./dz)*tmp_lc_o(i,j,1)
                 endif
              Enddo
           Enddo


           Do i=2,Nx-1
              Do j=2,Ny-1
                 if(wall(i,j,Nz)) then
                    aa(i,j,Nz) = 0.
                    bb(i,j,Nz) = 0.
                    cc(i,j,Nz) = 0.
                    dd(i,j,Nz) = 0.
                    ee(i,j,Nz) = 0.
                    ff(i,j,Nz) = 0.
                    gg(i,j,Nz) = 0.
                    q(i,j,Nz) = 0.
                 else
                    aa(i,j,Nz) = theta(i,j,Nz) + dt*(-min(u(i,j,Nz),0.)*(theta(i-1,j,Nz)+theta(i,j,Nz))/2./dx + &
                         (Dlc_eff(i-1,j,Nz)+Dlc_eff(i,j,Nz))/2./dx*(theta(i-1,j,Nz)+theta(i,j,Nz))/2./dx + &
                         max(u(i+1,j,Nz),0.)*(theta(i+1,j,Nz)+theta(i,j,Nz))/2./dx + &
                         (Dlc_eff(i+1,j,Nz)+Dlc_eff(i,j,Nz))/2./dx*(theta(i+1,j,Nz)+theta(i,j,Nz))/2./dx - &
                         min(v(i,j,Nz),0.)*(theta(i,j-1,Nz)+theta(i,j,Nz))/2./dy + &
                         (Dlc_eff(i,j-1,Nz)+Dlc_eff(i,j,Nz))/2./dy*(theta(i,j-1,Nz)+theta(i,j,Nz))/2./dy + &
                         max(v(i,j+1,Nz),0.)*(theta(i,j+1,Nz)+theta(i,j,Nz))/2./dy + &
                         (Dlc_eff(i,j+1,Nz)+Dlc_eff(i,j,Nz))/2./dy*(theta(i,j+1,Nz)+theta(i,j,Nz))/2./dy - &
                         min(w(i,j,Nz),0.)*(theta(i,j,Nz-1)+theta(i,j,Nz))/2./dz + &
                         (Dlc_eff(i,j,Nz-1)+Dlc_eff(i,j,Nz))/2./dz*(theta(i,j,Nz-1)+theta(i,j,Nz))/2./dz + &
                         max(w(i,j,Nz+1),0.)*(tmp_theta(i,j,2)+theta(i,j,Nz))/2./dz + &
                         (tmp_Dlc_eff(i,j,2)+Dlc_eff(i,j,Nz))/2./dz*(tmp_theta(i,j,2)+theta(i,j,Nz))/2./dz) - &
                         dt*theta(i,j,Nz)*(-qc_0*b_o(i,j,Nz)/(Kc+lc_o(i,j,Nz))*lo_o(i,j,Nz)/(Ko+lo_o(i,j,Nz))) + &
                         dt*theta(i,j,Nz)*alpha_c(i,j,Nz)

                    bb(i,j,Nz) = -dt*(max(u(i,j,Nz),0.)+(Dlc_eff(i-1,j,Nz)+Dlc_eff(i,j,Nz))/2./dx) * &
                         (theta(i-1,j,Nz)+theta(i,j,Nz))/2./dx
                    cc(i,j,Nz) = dt*(min(u(i+1,j,Nz),0.)-(Dlc_eff(i+1,j,Nz)+Dlc_eff(i,j,Nz))/2./dx) * &
                         (theta(i+1,j,Nz)+theta(i,j,Nz))/2./dx
                    dd(i,j,Nz) = -dt*(max(v(i,j,Nz),0.)+(Dlc_eff(i,j-1,Nz)+Dlc_eff(i,j,Nz))/2./dy) * &
                         (theta(i,j-1,Nz)+theta(i,j,Nz))/2./dy
                    ee(i,j,Nz) = dt*(min(v(i,j+1,Nz),0.)-(Dlc_eff(i,j+1,Nz)+Dlc_eff(i,j,Nz))/2./dy) * &
                         (theta(i,j+1,Nz)+theta(i,j,Nz))/2./dy
                    ff(i,j,Nz) = -dt*(max(w(i,j,Nz),0.)+(Dlc_eff(i,j,Nz-1)+Dlc_eff(i,j,Nz))/2./dz) * &
                         (theta(i,j,Nz-1)+theta(i,j,Nz))/2./dz
                    gg(i,j,Nz) = 0.

                    q(i,j,Nz) = theta(i,j,Nz)*lc_o(i,j,Nz)  + &
                         dt*theta(i,j,Nz)*alpha_c(i,j,Nz)*Ka*sc_o(i,j,Nz) - &
                         (dt*(min(w(i,j,Nz+1),0.)-(tmp_Dlc_eff(i,j,2)+Dlc_eff(i,j,Nz))/2./dz) * &
                         (tmp_theta(i,j,2)+theta(i,j,Nz))/2./dz)*tmp_lc_o(i,j,2)
                 endif
              Enddo
           Enddo



           !! for inside

           Do i=2,Nx-1
              Do j=2,Ny-1
                 Do k=2,Nz-1
                    if(wall(i,j,k)) then
                       aa(i,j,k) = 0.
                       bb(i,j,k) = 0.
                       cc(i,j,k) = 0.
                       dd(i,j,k) = 0.
                       ee(i,j,k) = 0.
                       ff(i,j,k) = 0.
                       gg(i,j,k) = 0.
                       q(i,j,k) = 0.
                    else
                       aa(i,j,k) = theta(i,j,k) + dt*(-min(u(i,j,k),0.)*(theta(i-1,j,k)+theta(i,j,k))/2./dx + &
                            (Dlc_eff(i-1,j,k)+Dlc_eff(i,j,k))/2./dx*(theta(i-1,j,k)+theta(i,j,k))/2./dx + &
                            max(u(i+1,j,k),0.)*(theta(i+1,j,k)+theta(i,j,k))/2./dx + &
                            (Dlc_eff(i+1,j,k)+Dlc_eff(i,j,k))/2./dx*(theta(i+1,j,k)+theta(i,j,k))/2./dx - &
                            min(v(i,j,k),0.)*(theta(i,j-1,k)+theta(i,j,k))/2./dy + &
                            (Dlc_eff(i,j-1,k)+Dlc_eff(i,j,k))/2./dy*(theta(i,j-1,k)+theta(i,j,k))/2./dy + &
                            max(v(i,j+1,k),0.)*(theta(i,j+1,k)+theta(i,j,k))/2./dy + &
                            (Dlc_eff(i,j+1,k)+Dlc_eff(i,j,k))/2./dy*(theta(i,j+1,k)+theta(i,j,k))/2./dy - &
                            min(w(i,j,k),0.)*(theta(i,j,k-1)+theta(i,j,k))/2./dz + &
                            (Dlc_eff(i,j,k-1)+Dlc_eff(i,j,k))/2./dz*(theta(i,j,k-1)+theta(i,j,k))/2./dz + &
                            max(w(i,j,k+1),0.)*(theta(i,j,k+1)+theta(i,j,k))/2./dz + &
                            (Dlc_eff(i,j,k+1)+Dlc_eff(i,j,k))/2./dz*(theta(i,j,k+1)+theta(i,j,k))/2./dz) - &
                            dt*theta(i,j,k)*(-qc_0*b_o(i,j,k)/(Kc+lc_o(i,j,k))*lo_o(i,j,k)/(Ko+lo_o(i,j,k))) + &
                            dt*theta(i,j,k)*alpha_c(i,j,k)

                       bb(i,j,k) = -dt*(max(u(i,j,k),0.)+(Dlc_eff(i-1,j,k)+Dlc_eff(i,j,k))/2./dx) * &
                            (theta(i-1,j,k)+theta(i,j,k))/2./dx
                       cc(i,j,k) = dt*(min(u(i+1,j,k),0.)-(Dlc_eff(i+1,j,k)+Dlc_eff(i,j,k))/2./dx) * &
                            (theta(i+1,j,k)+theta(i,j,k))/2./dx
                       dd(i,j,k) = -dt*(max(v(i,j,k),0.)+(Dlc_eff(i,j-1,k)+Dlc_eff(i,j,k))/2./dy) * &
                            (theta(i,j-1,k)+theta(i,j,k))/2./dy
                       ee(i,j,k) = dt*(min(v(i,j+1,k),0.)-(Dlc_eff(i,j+1,k)+Dlc_eff(i,j,k))/2./dy) * &
                            (theta(i,j+1,k)+theta(i,j,k))/2./dy
                       ff(i,j,k) = -dt*(max(w(i,j,k),0.)+(Dlc_eff(i,j,k-1)+Dlc_eff(i,j,k))/2./dz) * &
                            (theta(i,j,k-1)+theta(i,j,k))/2./dz
                       gg(i,j,k) = dt*(min(w(i,j,k+1),0.)-(Dlc_eff(i,j,k+1)+Dlc_eff(i,j,k))/2./dz) * &
                            (theta(i,j,k+1)+theta(i,j,k))/2./dz

                       q(i,j,k) = theta(i,j,k)*lc_o(i,j,k)  + &
                            dt*theta(i,j,k)*alpha_c(i,j,k)*Ka*sc_o(i,j,k)
                    endif
                 Enddo
              Enddo
           Enddo


        endif
     endif






     !! cabulate seven-diagonal elements near wall


     Do i=2,Nx-1
        Do j=2,Ny-1
           Do k=1,Nz
              if(wall(i-1,j,k)) then
                 aa(i,j,k) = aa(i,j,k) - dt*(Dlc_eff(i-1,j,k)+Dlc_eff(i,j,k))/2./dx*(theta(i-1,j,k)+theta(i,j,k))/2./dx
                 bb(i,j,k) = 0.
              endif
              if(wall(i+1,j,k)) then
                 aa(i,j,k) = aa(i,j,k) - dt*(Dlc_eff(i+1,j,k)+Dlc_eff(i,j,k))/2./dx*(theta(i+1,j,k)+theta(i,j,k))/2./dx
                 cc(i,j,k) = 0.
              endif
              if(wall(i,j-1,k)) then
                 aa(i,j,k) = aa(i,j,k) - dt*(Dlc_eff(i,j-1,k)+Dlc_eff(i,j,k))/2./dy*(theta(i,j-1,k)+theta(i,j,k))/2./dy
                 dd(i,j,k) = 0.
              endif
              if(wall(i,j+1,k)) then
                 aa(i,j,k) = aa(i,j,k) - dt*(Dlc_eff(i,j+1,k)+Dlc_eff(i,j,k))/2./dy*(theta(i,j+1,k)+theta(i,j,k))/2./dy
                 ee(i,j,k) = 0.
              endif
              
              if(npes.eq.1) then
                 if(k.ne.1 .and. wall(i,j,k-1)) then
                    aa(i,j,k) = aa(i,j,k) - dt*(Dlc_eff(i,j,k-1)+Dlc_eff(i,j,k))/2./dz*(theta(i,j,k-1)+theta(i,j,k))/2./dz
                    ff(i,j,k) = 0.
                 endif
                 if(k.ne.Nz .and. wall(i,j,k+1)) then
                    aa(i,j,k) = aa(i,j,k) - dt*(Dlc_eff(i,j,k+1)+Dlc_eff(i,j,k))/2./dz*(theta(i,j,k+1)+theta(i,j,k))/2./dz
                    gg(i,j,k) = 0.
                 endif
              else
                 if(mype.eq.0) then
                    if(k.ne.1 .and. wall(i,j,k-1)) then
                       aa(i,j,k) = aa(i,j,k) - dt*(Dlc_eff(i,j,k-1)+Dlc_eff(i,j,k))/2./dz*(theta(i,j,k-1)+theta(i,j,k))/2./dz
                       ff(i,j,k) = 0.
                    endif
                    if(k.eq.Nz) then
                       if(tmp_wall(i,j,2)) then
                          aa(i,j,Nz) = aa(i,j,Nz) - &
                               dt*(tmp_Dlc_eff(i,j,2)+Dlc_eff(i,j,Nz))/2./dz*(tmp_theta(i,j,2)+theta(i,j,Nz))/2./dz
                          gg(i,j,Nz) = 0.
                       endif
                    else
                       if(wall(i,j,k+1)) then
                          aa(i,j,k) = aa(i,j,k) - dt*(Dlc_eff(i,j,k+1)+Dlc_eff(i,j,k))/2./dz*(theta(i,j,k+1)+theta(i,j,k))/2./dz
                          gg(i,j,k) = 0.
                       endif
                    endif
                 else if(mype.eq.npes-1) then
                    if(k.eq.1) then
                       if(tmp_wall(i,j,1)) then
                          aa(i,j,1) = aa(i,j,1) - dt*(tmp_Dlc_eff(i,j,1)+Dlc_eff(i,j,1))/2./dz*(tmp_theta(i,j,1)+theta(i,j,1))/2./dz
                          ff(i,j,1) = 0.
                       endif
                    else
                       if(wall(i,j,k-1)) then
                          aa(i,j,k) = aa(i,j,k) - dt*(Dlc_eff(i,j,k-1)+Dlc_eff(i,j,k))/2./dz*(theta(i,j,k-1)+theta(i,j,k))/2./dz
                          ff(i,j,k) = 0.
                       endif
                    endif
                    if(k.ne.Nz .and. wall(i,j,k+1)) then
                       aa(i,j,k) = aa(i,j,k) - dt*(Dlc_eff(i,j,k+1)+Dlc_eff(i,j,k))/2./dz*(theta(i,j,k+1)+theta(i,j,k))/2./dz
                       gg(i,j,k) = 0.
                    endif
                 else 
                    if(k.eq.1) then
                       if(tmp_wall(i,j,1)) then
                          aa(i,j,1) = aa(i,j,1) - dt*(tmp_Dlc_eff(i,j,1)+Dlc_eff(i,j,1))/2./dz*(tmp_theta(i,j,1)+theta(i,j,1))/2./dz
                          ff(i,j,1) = 0.
                       endif
                    else
                       if(wall(i,j,k-1)) then
                          aa(i,j,k) = aa(i,j,k) - dt*(Dlc_eff(i,j,k-1)+Dlc_eff(i,j,k))/2./dz*(theta(i,j,k-1)+theta(i,j,k))/2./dz
                          ff(i,j,k) = 0.
                       endif
                    endif
                    if(k.eq.Nz) then
                       if(tmp_wall(i,j,2)) then
                          aa(i,j,Nz) = aa(i,j,Nz) - &
                               dt*(tmp_Dlc_eff(i,j,2)+Dlc_eff(i,j,Nz))/2./dz*(tmp_theta(i,j,2)+theta(i,j,Nz))/2./dz
                          gg(i,j,Nz) = 0.
                       endif
                    else
                       if(wall(i,j,k+1)) then
                          aa(i,j,k) = aa(i,j,k) - dt*(Dlc_eff(i,j,k+1)+Dlc_eff(i,j,k))/2./dz*(theta(i,j,k+1)+theta(i,j,k))/2./dz
                          gg(i,j,k) = 0.
                       endif
                    endif
                 endif
              endif

           Enddo
        Enddo
     Enddo




!!! calculate matrix elements for all side planes 




     !! calculate elements of matrix MM 


     Do k=1,Nz
        Do j=1,Ny
           Do i=1,Nx
              ma((k-1)*Nx*Ny+(j-1)*Nx+i) = aa(i,j,k)
              mb((k-1)*Nx*Ny+(j-1)*Nx+i) = bb(i,j,k)
              mc((k-1)*Nx*Ny+(j-1)*Nx+i) = cc(i,j,k)
              md((k-1)*Nx*Ny+(j-1)*Nx+i) = dd(i,j,k)
              me((k-1)*Nx*Ny+(j-1)*Nx+i) = ee(i,j,k)
              mf((k-1)*Nx*Ny+(j-1)*Nx+i) = ff(i,j,k)
              mg((k-1)*Nx*Ny+(j-1)*Nx+i) = gg(i,j,k)

              QQ((k-1)*Nx*Ny+(j-1)*Nx+i) = q(i,j,k)
              HOHO((k-1)*Nx*Ny+(j-1)*Nx+i) = lc_o(i,j,k)

              WW((k-1)*Nx*Ny+(j-1)*Nx+i) = wall(i,j,k)
           Enddo
        Enddo
     Enddo



!!$!notice 
!!$     if(mod(istep,Print_step).eq.0) then
!!$        Print*,'mype=',mype,'q'
!!$        print*,q(nxx,nyy,:)
!!$     endif





!!! solve lc by successive overrelaxation (SOR)



        nn = 1

        Do while(nn.lt.1e5)

           if(WW(1)) then
              HH(1) = 0.
           else
              HH(1) = omega*(QQ(1)-mc(1)*HOHO(2)-me(1)*HOHO(1+Nx)-mg(1)*HOHO(1+Nx*Ny))/ &
                   (ma(1)+SMALL) + (1-omega)*HOHO(1)
           endif
           Do pp=2,Nx
              if(WW(pp)) then
                 HH(pp) = 0.
              else
                 HH(pp) = omega*(QQ(pp)-mb(pp)*HH(pp-1)-mc(pp)*HOHO(pp+1)-me(pp)*HOHO(pp+Nx)- &
                      mg(pp)*HOHO(pp+Nx*Ny))/(ma(pp)+SMALL) + (1-omega)*HOHO(pp)
              endif
           Enddo
           Do pp=1+Nx,Nx*Ny
              if(WW(pp)) then
                 HH(pp) = 0.
              else
                 HH(pp) = omega*(QQ(pp)-md(pp)*HH(pp-Nx)-mb(pp)*HH(pp-1)-mc(pp)*HOHO(pp+1)- &
                      me(pp)*HOHO(pp+Nx)-mg(pp)*HOHO(pp+Nx*Ny))/(ma(pp)+SMALL) + &
                      (1-omega)*HOHO(pp)
              endif
           Enddo
           Do pp=1+Nx*Ny,Nt-Nx*Ny
              if(WW(pp)) then
                 HH(pp) = 0.
              else
                 HH(pp) = omega*(QQ(pp)-mf(pp)*HH(pp-Nx*Ny)-md(pp)*HH(pp-Nx)- &
                      mb(pp)*HH(pp-1)-mc(pp)*HOHO(pp+1)-me(pp)*HOHO(pp+Nx)- &
                      mg(pp)*HOHO(pp+Nx*Ny))/(ma(pp)+SMALL) + (1-omega)*HOHO(pp)
              endif
           Enddo
           Do pp=Nt-Nx*Ny+1,Nt-Nx
              if(WW(pp)) then
                 HH(pp) = 0.
              else
                 HH(pp) = omega*(QQ(pp)-mf(pp)*HH(pp-Nx*Ny)-md(pp)*HH(pp-Nx)- &
                      mb(pp)*HH(pp-1)-mc(pp)*HOHO(pp+1)-me(pp)*HOHO(pp+Nx))/(ma(pp)+SMALL) + &
                      (1-omega)*HOHO(pp)
              endif
           Enddo
           Do pp=Nt-Nx+1,Nt-1
              if(WW(pp)) then
                 HH(pp) = 0.
              else
                 HH(pp) = omega*(QQ(pp)-mf(pp)*HH(pp-Nx*Ny)-md(pp)*HH(pp-Nx)- &
                      mb(pp)*HH(pp-1)-mc(pp)*HOHO(pp+1))/(ma(pp)+SMALL) + (1-omega)*HOHO(pp)
              endif
           Enddo
           if(WW(Nt)) then
              HH(Nt) = 0.
           else
              HH(Nt) = omega*(QQ(Nt)-mf(Nt)*HH(Nt-Nx*Ny)-md(Nt)*HH(Nt-Nx)- &
                   mb(Nt)*HH(Nt-1))/(ma(Nt)+SMALL) +(1-omega)*HOHO(Nt)
           endif


!!! calculate error
           HH_err_local = Sum(Abs(HH-HOHO))
           CALL MPI_ALLREDUCE(HH_err_local,HH_err,1,MPI_double_precision,MPI_SUM,MPI_COMM_WORLD,error)
           HH_err = HH_err/LX/LY/LZ


!!$!notice 
!!$     if(mod(istep,Print_step).eq.0) then
!!$        Print*,'mype=',mype,'HH_err_local=',HH_err_local,'HH_err=',HH_err
!!$     endif


           If (HH_err.lt.HH_eps) exit


!!! update interface q

!!$        Do k=1,Nz
!!$           Do j=1,Ny
!!$              Do i=1,Nx
!!$                 lc(i,j,k) = HH((k-1)*Nx*Ny+(j-1)*Nx+i)
!!$              Enddo
!!$           Enddo
!!$        Enddo

           Do j=1,Ny
              Do i=1,Nx
                 lc(i,j,1) = HH((j-1)*Nx+i)
              Enddo
           Enddo

           Do j=1,Ny
              Do i=1,Nx
                 lc(i,j,Nz) = HH((Nz-1)*Nx*Ny+(j-1)*Nx+i)
              Enddo
           Enddo








!!! exchange h during iteration

           call mpi_sendrecv(lc(:,:,Nz),ilen,MPI_double_precision,inext,610, &
                tmp_lc(:,:,1),ilen,MPI_double_precision,iprev,610,nallgrp,stat,error)

           call mpi_sendrecv(lc(:,:,1),ilen,MPI_double_precision,iprev,620, &
                tmp_lc(:,:,2),ilen,MPI_double_precision,inext,620,nallgrp,stat,error)


           CALL MPI_BARRIER(MPI_COMM_WORLD,error)



           
           if(mype.ne.0) then
              Do j=1,Ny
                 Do i=1,Nx
                    q(i,j,1) = theta(i,j,1)*lc_o(i,j,1)  + &
                         dt*theta(i,j,1)*alpha_c(i,j,1)*Ka*sc_o(i,j,1) - &
                         (-dt*(max(w(i,j,1),0.)+(tmp_Dlc_eff(i,j,1)+Dlc_eff(i,j,1))/2./dz) * &
                         (tmp_theta(i,j,1)+theta(i,j,1))/2./dz)*tmp_lc(i,j,1)
                    QQ((j-1)*Nx+i) = q(i,j,1)
                 Enddo
              Enddo
           endif

           if(mype.ne.npes-1) then
              Do j=1,Ny
                 Do i=1,Nx
                    q(i,j,Nz) = theta(i,j,Nz)*lc_o(i,j,Nz)  + &
                         dt*theta(i,j,Nz)*alpha_c(i,j,Nz)*Ka*sc_o(i,j,Nz) - &
                         (dt*(min(w(i,j,Nz+1),0.)-(tmp_Dlc_eff(i,j,2)+Dlc_eff(i,j,Nz))/2./dz) * &
                         (tmp_theta(i,j,2)+theta(i,j,Nz))/2./dz)*tmp_lc(i,j,2)
                    QQ((Nz-1)*Nx*Ny+(j-1)*Nx+i) = q(i,j,Nz)
                 Enddo
              Enddo
           endif



           HOHO = HH
           nn = nn+1

        ENDDO




        !! update h

        Do k=1,Nz
           Do j=1,Ny
              Do i=1,Nx
                 lc(i,j,k) = HH((k-1)*Nx*Ny+(j-1)*Nx+i)
              Enddo
           Enddo
        Enddo

















!!!! solve sc semi-implicitl



     Do k=1,Nz
        Do j=1,Ny
           Do i=1,Nx
              if(ne(i,j,k).ne.1.0 .and. .not.wall(i,j,k)) then
                 sc(i,j,k) = 1./(1.+dt*a1*theta(i,j,k)*Dlc/(theta(i,j,k)+b1)*SA) * &
                      (sc_o(i,j,k) + dt*a1*theta(i,j,k)*Dlc/(theta(i,j,k)+b1)*SA*Kcc*lc_o(i,j,k))
              endif
           Enddo
        Enddo
     Enddo







!!! total carbon conc.


     Do k=1,Nz
        Do j=1,Ny
           Do i=1,Nx
              if(.not.wall(i,j,k)) then
                 c(i,j,k) = lc(i,j,k)*theta(i,j,k) + sc(i,j,k)*rhos*(1.-ne(i,j,k))
              endif
           Enddo
        Enddo
     Enddo





!!$!!! update aqueous carbon conc due to equilibrium with soil carbon 
!!$
!!$
!!$
!!$     Do k=1,Nz
!!$        Do j=1,Ny
!!$           Do i=1,Nx
!!$              if(.not.wall(i,j,k)) then
!!$                 c(i,j,k) = lc(i,j,k)*ne(i,j,k) + sc_o(i,j,k)*rhos*(1.-ne(i,j,k))
!!$                 sc(i,j,k) = c(i,j,k)/(Ka*ne(i,j,k)+rhos*(1.-ne(i,j,k)))
!!$                 lc(i,j,k) = Ka*c(i,j,k)/(Ka*ne(i,j,k)+rhos*(1.-ne(i,j,k)))
!!$              endif
!!$           Enddo
!!$        Enddo
!!$     Enddo





!!!!!! solve b explicitly



     Do i=1,Nx
        Do j=1,Ny
           Do k=1,Nz
              if(wall(i,j,k)) then
                 b(i,j,k) = 0.
              else
                 b(i,j,k) = b_o(i,j,k) + &
                      dt*(Yd*qc_0*b_o(i,j,k)*lc_o(i,j,k)/(Kc+lc_o(i,j,k))*lo_o(i,j,k)/(Ko+lo_o(i,j,k)) - Dk*b_o(i,j,k))
              endif
           Enddo
        Enddo
     Enddo









!!!! solve to (total oxygen) semi-implicitly



!!! boundary condition setting and
!!! communicate info between bottom and top planes


     call mpi_sendrecv(to_o(:,:,Nz),ilen,MPI_double_precision,inext,110, &
          tmp_to_o(:,:,1),ilen,MPI_double_precision,iprev,110,nallgrp,stat,error)
     call mpi_sendrecv(to_o(:,:,1),ilen,MPI_double_precision,iprev,120, &
          tmp_to_o(:,:,2),ilen,MPI_double_precision,inext,120,nallgrp,stat,error)


     CALL MPI_BARRIER(MPI_COMM_WORLD,error)

     

     
     if(npes.eq.1) then

        Do i=2,Nx-1
           Do j=2,Ny-1
              if(wall(i,j,1)) then
                 aa(i,j,1) = 0.
                 bb(i,j,1) = 0.
                 cc(i,j,1) = 0.
                 dd(i,j,1) = 0.
                 ee(i,j,1) = 0.
                 ff(i,j,1) = 0.
                 gg(i,j,1) = 0.
                 q(i,j,1) = 0.
              else
                 bb(i,j,1) = -dt*((alpha1_o(i-1,j,1)+alpha1_o(i,j,1))/2.*max(u(i,j,1),0.) + &
                      ((alpha1_o(i-1,j,1)+alpha1_o(i,j,1))/2.*(Dlo_eff(i-1,j,1)+Dlo_eff(i,j,1))/2. + &
                      (alpha2_o(i-1,j,1)+alpha2_o(i,j,1))/2.*(Dgo_eff(i-1,j,1)+Dgo_eff(i,j,1))/2.)/dx)* &
                      (ne(i-1,j,1)+ne(i,j,1))/2./dx
                 cc(i,j,1) = dt*((alpha1_o(i+1,j,1)+alpha1_o(i,j,1))/2.*min(u(i+1,j,1),0.) - &
                      ((alpha1_o(i+1,j,1)+alpha1_o(i,j,1))/2.*(Dlo_eff(i+1,j,1)+Dlo_eff(i,j,1))/2. + &
                      (alpha2_o(i+1,j,1)+alpha2_o(i,j,1))/2.*(Dgo_eff(i+1,j,1)+Dgo_eff(i,j,1))/2.)/dx)* &
                      (ne(i+1,j,1)+ne(i,j,1))/2./dx
                 dd(i,j,1) = -dt*((alpha1_o(i,j-1,1)+alpha1_o(i,j,1))/2.*max(v(i,j,1),0.) + &
                      ((alpha1_o(i,j-1,1)+alpha1_o(i,j,1))/2.*(Dlo_eff(i,j-1,1)+Dlo_eff(i,j,1))/2. + &
                      (alpha2_o(i,j-1,1)+alpha2_o(i,j,1))/2.*(Dgo_eff(i,j-1,1)+Dgo_eff(i,j,1))/2.)/dy)* &
                      (ne(i,j-1,1)+ne(i,j,1))/2./dy
                 ee(i,j,1) = dt*((alpha1_o(i,j+1,1)+alpha1_o(i,j,1))/2.*min(v(i,j+1,1),0.) - &
                      ((alpha1_o(i,j+1,1)+alpha1_o(i,j,1))/2.*(Dlo_eff(i,j+1,1)+Dlo_eff(i,j,1))/2. + &
                      (alpha2_o(i,j+1,1)+alpha2_o(i,j,1))/2.*(Dgo_eff(i,j+1,1)+Dgo_eff(i,j,1))/2.)/dy)* &
                      (ne(i,j+1,1)+ne(i,j,1))/2./dy
                 ff(i,j,1) = 0.
                 gg(i,j,1) = dt*((alpha1_o(i,j,1+1)+alpha1_o(i,j,1))/2.*min(w(i,j,1+1),0.) - &
                      ((alpha1_o(i,j,1+1)+alpha1_o(i,j,1))/2.*(Dlo_eff(i,j,1+1)+Dlo_eff(i,j,1))/2. + &
                      (alpha2_o(i,j,1+1)+alpha2_o(i,j,1))/2.*(Dgo_eff(i,j,1+1)+Dgo_eff(i,j,1))/2.)/dz)* &
                      (ne(i,j,1+1)+ne(i,j,1))/2./dz


                 if(tobc_d.eq.1) then 
                    aa(i,j,1) = ne(i,j,1) + dt*(-(alpha1_o(i-1,j,1)+alpha1_o(i,j,1))/2.*min(u(i,j,1),0.)* &
                         (ne(i-1,j,1)+ne(i,j,1))/2./dx + &
                         ((alpha1_o(i-1,j,1)+alpha1_o(i,j,1))/2.*(Dlo_eff(i-1,j,1)+Dlo_eff(i,j,1))/2. + &
                         (alpha2_o(i-1,j,1)+alpha2_o(i,j,1))/2.*(Dgo_eff(i-1,j,1)+Dgo_eff(i,j,1))/2.)/dx* &
                         (ne(i-1,j,1)+ne(i,j,1))/2./dx + &
                         (alpha1_o(i+1,j,1)+alpha1_o(i,j,1))/2.*max(u(i+1,j,1),0.)* &
                         (ne(i+1,j,1)+ne(i,j,1))/2./dx + &
                         ((alpha1_o(i+1,j,1)+alpha1_o(i,j,1))/2.*(Dlo_eff(i+1,j,1)+Dlo_eff(i,j,1))/2. + &
                         (alpha2_o(i+1,j,1)+alpha2_o(i,j,1))/2.*(Dgo_eff(i+1,j,1)+Dgo_eff(i,j,1))/2.)/dx* &
                         (ne(i+1,j,1)+ne(i,j,1))/2./dx - &
                         (alpha1_o(i,j-1,1)+alpha1_o(i,j,1))/2.*min(v(i,j,1),0.)*(ne(i,j-1,1)+ne(i,j,1))/2./dy + &
                         ((alpha1_o(i,j-1,1)+alpha1_o(i,j,1))/2.*(Dlo_eff(i,j-1,1)+Dlo_eff(i,j,1))/2. + &
                         (alpha2_o(i,j-1,1)+alpha2_o(i,j,1))/2.*(Dgo_eff(i,j-1,1)+Dgo_eff(i,j,1))/2.)/dy* &
                         (ne(i,j-1,1)+ne(i,j,1))/2./dy + &
                         (alpha1_o(i,j+1,1)+alpha1_o(i,j,1))/2.*max(v(i,j+1,1),0.)*(ne(i,j+1,1)+ne(i,j,1))/2./dy + &
                         ((alpha1_o(i,j+1,1)+alpha1_o(i,j,1))/2.*(Dlo_eff(i,j+1,1)+Dlo_eff(i,j,1))/2. + &
                         (alpha2_o(i,j+1,1)+alpha2_o(i,j,1))/2.*(Dgo_eff(i,j+1,1)+Dgo_eff(i,j,1))/2.)/dy* &
                         (ne(i,j+1,1)+ne(i,j,1))/2./dy + &
                         (alpha1_o(i,j,1)*Dlo_eff(i,j,1)+alpha2_o(i,j,1)*Dgo_eff(i,j,1))/dz*2.*ne(i,j,1)/dz + &
                         (alpha1_o(i,j,1+1)+alpha1_o(i,j,1))/2.*max(w(i,j,1+1),0.)*(ne(i,j,1+1)+ne(i,j,1))/2./dz + &
                         ((alpha1_o(i,j,1+1)+alpha1_o(i,j,1))/2.*(Dlo_eff(i,j,1+1)+Dlo_eff(i,j,1))/2. + &
                         (alpha2_o(i,j,1+1)+alpha2_o(i,j,1))/2.*(Dgo_eff(i,j,1+1)+Dgo_eff(i,j,1))/2.)/dz* &
                         (ne(i,j,1+1)+ne(i,j,1))/2./dz) - &
                         dt*ne(i,j,1)*(-alpha1_o(i,j,1)*qo_0*b_o(i,j,1)* &
                         lc_o(i,j,1)/(Kc+lc_o(i,j,1))/(Ko+lo_o(i,j,1)))

                    q(i,j,1) = ne(i,j,1)*to_o(i,j,1) - &
                            (-dt*(alpha1_o(i,j,1)*max(w(i,j,1),0.) + &
                            2.*(alpha1_o(i,j,1)*Dlo_eff(i,j,1)+alpha2_o(i,j,1)*Dgo_eff(i,j,1))/dz)* &
                            ne(i,j,1)/dz)*to_d(i,j)
                 endif
                 if(tobc_d.eq.2) then 
                    aa(i,j,1) = ne(i,j,1) + dt*(-(alpha1_o(i-1,j,1)+alpha1_o(i,j,1))/2.*min(u(i,j,1),0.)* &
                         (ne(i-1,j,1)+ne(i,j,1))/2./dx + &
                         ((alpha1_o(i-1,j,1)+alpha1_o(i,j,1))/2.*(Dlo_eff(i-1,j,1)+Dlo_eff(i,j,1))/2. + &
                         (alpha2_o(i-1,j,1)+alpha2_o(i,j,1))/2.*(Dgo_eff(i-1,j,1)+Dgo_eff(i,j,1))/2.)/dx* &
                         (ne(i-1,j,1)+ne(i,j,1))/2./dx + &
                         (alpha1_o(i+1,j,1)+alpha1_o(i,j,1))/2.*max(u(i+1,j,1),0.)* &
                         (ne(i+1,j,1)+ne(i,j,1))/2./dx + &
                         ((alpha1_o(i+1,j,1)+alpha1_o(i,j,1))/2.*(Dlo_eff(i+1,j,1)+Dlo_eff(i,j,1))/2. + &
                         (alpha2_o(i+1,j,1)+alpha2_o(i,j,1))/2.*(Dgo_eff(i+1,j,1)+Dgo_eff(i,j,1))/2.)/dx* &
                         (ne(i+1,j,1)+ne(i,j,1))/2./dx - &
                         (alpha1_o(i,j-1,1)+alpha1_o(i,j,1))/2.*min(v(i,j,1),0.)*(ne(i,j-1,1)+ne(i,j,1))/2./dy + &
                         ((alpha1_o(i,j-1,1)+alpha1_o(i,j,1))/2.*(Dlo_eff(i,j-1,1)+Dlo_eff(i,j,1))/2. + &
                         (alpha2_o(i,j-1,1)+alpha2_o(i,j,1))/2.*(Dgo_eff(i,j-1,1)+Dgo_eff(i,j,1))/2.)/dy* &
                         (ne(i,j-1,1)+ne(i,j,1))/2./dy + &
                         (alpha1_o(i,j+1,1)+alpha1_o(i,j,1))/2.*max(v(i,j+1,1),0.)*(ne(i,j+1,1)+ne(i,j,1))/2./dy + &
                         ((alpha1_o(i,j+1,1)+alpha1_o(i,j,1))/2.*(Dlo_eff(i,j+1,1)+Dlo_eff(i,j,1))/2. + &
                         (alpha2_o(i,j+1,1)+alpha2_o(i,j,1))/2.*(Dgo_eff(i,j+1,1)+Dgo_eff(i,j,1))/2.)/dy* &
                         (ne(i,j+1,1)+ne(i,j,1))/2./dy - &
                         alpha1_o(i,j,1)*min(w(i,j,1),0.)*ne(i,j,1)/dz + &
                         (alpha1_o(i,j,1+1)+alpha1_o(i,j,1))/2.*max(w(i,j,1+1),0.)*(ne(i,j,1+1)+ne(i,j,1))/2./dz + &
                         ((alpha1_o(i,j,1+1)+alpha1_o(i,j,1))/2.*(Dlo_eff(i,j,1+1)+Dlo_eff(i,j,1))/2. + &
                         (alpha2_o(i,j,1+1)+alpha2_o(i,j,1))/2.*(Dgo_eff(i,j,1+1)+Dgo_eff(i,j,1))/2.)/dz* &
                         (ne(i,j,1+1)+ne(i,j,1))/2./dz) - &
                         dt*ne(i,j,1)*(-alpha1_o(i,j,1)*qo_0*b_o(i,j,1)* &
                         lc_o(i,j,1)/(Kc+lc_o(i,j,1))/(Ko+lo_o(i,j,1)))

                    q(i,j,1) = ne(i,j,1)*to_o(i,j,1)
                 endif
              endif
           enddo
        enddo


        Do i=2,Nx-1
           Do j=2,Ny-1
              if(wall(i,j,Nz)) then
                 aa(i,j,Nz) = 0.
                 bb(i,j,Nz) = 0.
                 cc(i,j,Nz) = 0.
                 dd(i,j,Nz) = 0.
                 ee(i,j,Nz) = 0.
                 ff(i,j,Nz) = 0.
                 gg(i,j,Nz) = 0.
                 q(i,j,Nz) = 0.
              else     
                 bb(i,j,Nz) = -dt*((alpha1_o(i-1,j,Nz)+alpha1_o(i,j,Nz))/2.*max(u(i,j,Nz),0.) + &
                      ((alpha1_o(i-1,j,Nz)+alpha1_o(i,j,Nz))/2.*(Dlo_eff(i-1,j,Nz)+Dlo_eff(i,j,Nz))/2. + &
                      (alpha2_o(i-1,j,Nz)+alpha2_o(i,j,Nz))/2.*(Dgo_eff(i-1,j,Nz)+Dgo_eff(i,j,Nz))/2.)/dx)* &
                      (ne(i-1,j,Nz)+ne(i,j,Nz))/2./dx
                 cc(i,j,Nz) = dt*((alpha1_o(i+1,j,Nz)+alpha1_o(i,j,Nz))/2.*min(u(i+1,j,Nz),0.) - &
                      ((alpha1_o(i+1,j,Nz)+alpha1_o(i,j,Nz))/2.*(Dlo_eff(i+1,j,Nz)+Dlo_eff(i,j,Nz))/2. + &
                      (alpha2_o(i+1,j,Nz)+alpha2_o(i,j,Nz))/2.*(Dgo_eff(i+1,j,Nz)+Dgo_eff(i,j,Nz))/2.)/dx)* &
                      (ne(i+1,j,Nz)+ne(i,j,Nz))/2./dx
                 dd(i,j,Nz) = -dt*((alpha1_o(i,j-1,Nz)+alpha1_o(i,j,Nz))/2.*max(v(i,j,Nz),0.) + &
                      ((alpha1_o(i,j-1,Nz)+alpha1_o(i,j,Nz))/2.*(Dlo_eff(i,j-1,Nz)+Dlo_eff(i,j,Nz))/2. + &
                      (alpha2_o(i,j-1,Nz)+alpha2_o(i,j,Nz))/2.*(Dgo_eff(i,j-1,Nz)+Dgo_eff(i,j,Nz))/2.)/dy)* &
                      (ne(i,j-1,Nz)+ne(i,j,Nz))/2./dy
                 ee(i,j,Nz) = dt*((alpha1_o(i,j+1,Nz)+alpha1_o(i,j,Nz))/2.*min(v(i,j+1,Nz),0.) - &
                      ((alpha1_o(i,j+1,Nz)+alpha1_o(i,j,Nz))/2.*(Dlo_eff(i,j+1,Nz)+Dlo_eff(i,j,Nz))/2. + &
                      (alpha2_o(i,j+1,Nz)+alpha2_o(i,j,Nz))/2.*(Dgo_eff(i,j+1,Nz)+Dgo_eff(i,j,Nz))/2.)/dy)* &
                      (ne(i,j+1,Nz)+ne(i,j,Nz))/2./dy
                 ff(i,j,Nz) = -dt*((alpha1_o(i,j,Nz-1)+alpha1_o(i,j,Nz))/2.*max(w(i,j,Nz),0.) + &
                      ((alpha1_o(i,j,Nz-1)+alpha1_o(i,j,Nz))/2.*(Dlo_eff(i,j,Nz-1)+Dlo_eff(i,j,Nz))/2. + &
                      (alpha2_o(i,j,Nz-1)+alpha2_o(i,j,Nz))/2.*(Dgo_eff(i,j,Nz-1)+Dgo_eff(i,j,Nz))/2.)/dz)* &
                      (ne(i,j,Nz-1)+ne(i,j,Nz))/2./dz
                 gg(i,j,Nz) = 0.

                 if(tobc_t.eq.1) then
                    aa(i,j,Nz) = ne(i,j,Nz) + dt*(-(alpha1_o(i-1,j,Nz)+alpha1_o(i,j,Nz))/2.*min(u(i,j,Nz),0.)* &
                         (ne(i-1,j,Nz)+ne(i,j,Nz))/2./dx + &
                         ((alpha1_o(i-1,j,Nz)+alpha1_o(i,j,Nz))/2.*(Dlo_eff(i-1,j,Nz)+Dlo_eff(i,j,Nz))/2. + &
                         (alpha2_o(i-1,j,Nz)+alpha2_o(i,j,Nz))/2.*(Dgo_eff(i-1,j,Nz)+Dgo_eff(i,j,Nz))/2.)/dx* &
                         (ne(i-1,j,Nz)+ne(i,j,Nz))/2./dx + &
                         (alpha1_o(i+1,j,Nz)+alpha1_o(i,j,Nz))/2.*max(u(i+1,j,Nz),0.)* &
                         (ne(i+1,j,Nz)+ne(i,j,Nz))/2./dx + &
                         ((alpha1_o(i+1,j,Nz)+alpha1_o(i,j,Nz))/2.*(Dlo_eff(i+1,j,Nz)+Dlo_eff(i,j,Nz))/2. + &
                         (alpha2_o(i+1,j,Nz)+alpha2_o(i,j,Nz))/2.*(Dgo_eff(i+1,j,Nz)+Dgo_eff(i,j,Nz))/2.)/dx* &
                         (ne(i+1,j,Nz)+ne(i,j,Nz))/2./dx - &
                         (alpha1_o(i,j-1,Nz)+alpha1_o(i,j,Nz))/2.*min(v(i,j,Nz),0.)*(ne(i,j-1,Nz)+ne(i,j,Nz))/2./dy + &
                         ((alpha1_o(i,j-1,Nz)+alpha1_o(i,j,Nz))/2.*(Dlo_eff(i,j-1,Nz)+Dlo_eff(i,j,Nz))/2. + &
                         (alpha2_o(i,j-1,Nz)+alpha2_o(i,j,Nz))/2.*(Dgo_eff(i,j-1,Nz)+Dgo_eff(i,j,Nz))/2.)/dy* &
                         (ne(i,j-1,Nz)+ne(i,j,Nz))/2./dy + &
                         (alpha1_o(i,j+1,Nz)+alpha1_o(i,j,Nz))/2.*max(v(i,j+1,Nz),0.)*(ne(i,j+1,Nz)+ne(i,j,Nz))/2./dy + &
                         ((alpha1_o(i,j+1,Nz)+alpha1_o(i,j,Nz))/2.*(Dlo_eff(i,j+1,Nz)+Dlo_eff(i,j,Nz))/2. + &
                         (alpha2_o(i,j+1,Nz)+alpha2_o(i,j,Nz))/2.*(Dgo_eff(i,j+1,Nz)+Dgo_eff(i,j,Nz))/2.)/dy* &
                         (ne(i,j+1,Nz)+ne(i,j,Nz))/2./dy - &
                         (alpha1_o(i,j,Nz-1)+alpha1_o(i,j,Nz))/2.*min(w(i,j,Nz),0.)*(ne(i,j,Nz-1)+ne(i,j,Nz))/2./dz + &
                         ((alpha1_o(i,j,Nz-1)+alpha1_o(i,j,Nz))/2.*(Dlo_eff(i,j,Nz-1)+Dlo_eff(i,j,Nz))/2. + &
                         (alpha2_o(i,j,Nz-1)+alpha2_o(i,j,Nz))/2.*(Dgo_eff(i,j,Nz-1)+Dgo_eff(i,j,Nz))/2.)/dz* &
                         (ne(i,j,Nz-1)+ne(i,j,Nz))/2./dz + &
                         (alpha1_o(i,j,Nz)*Dlo_eff(i,j,Nz) + alpha2_o(i,j,Nz)*Dgo_eff(i,j,Nz))/dz*2.*ne(i,j,Nz)/dz) - &
                         dt*ne(i,j,Nz)*(-alpha1_o(i,j,Nz)*qo_0*b_o(i,j,Nz)* &
                         lc_o(i,j,Nz)/(Kc+lc_o(i,j,Nz))/(Ko+lo_o(i,j,Nz)))

                    q(i,j,Nz) = ne(i,j,Nz)*to_o(i,j,Nz) - dt * (alpha1_o(i,j,Nz)*min(w(i,j,Nz+1),0.) - &
                            2.*(alpha1_o(i,j,Nz)*Dlo_eff(i,j,Nz)+alpha2_o(i,j,Nz)*Dgo_eff(i,j,Nz))/dz)*ne(i,j,Nz)/dz*to_t(i,j)
                 endif
                 if(tobc_t.eq.2) then
                    aa(i,j,Nz) = ne(i,j,Nz) + dt*(-(alpha1_o(i-1,j,Nz)+alpha1_o(i,j,Nz))/2.*min(u(i,j,Nz),0.)* &
                         (ne(i-1,j,Nz)+ne(i,j,Nz))/2./dx + &
                         ((alpha1_o(i-1,j,Nz)+alpha1_o(i,j,Nz))/2.*(Dlo_eff(i-1,j,Nz)+Dlo_eff(i,j,Nz))/2. + &
                         (alpha2_o(i-1,j,Nz)+alpha2_o(i,j,Nz))/2.*(Dgo_eff(i-1,j,Nz)+Dgo_eff(i,j,Nz))/2.)/dx* &
                         (ne(i-1,j,Nz)+ne(i,j,Nz))/2./dx + &
                         (alpha1_o(i+1,j,Nz)+alpha1_o(i,j,Nz))/2.*max(u(i+1,j,Nz),0.)* &
                         (ne(i+1,j,Nz)+ne(i,j,Nz))/2./dx + &
                         ((alpha1_o(i+1,j,Nz)+alpha1_o(i,j,Nz))/2.*(Dlo_eff(i+1,j,Nz)+Dlo_eff(i,j,Nz))/2. + &
                         (alpha2_o(i+1,j,Nz)+alpha2_o(i,j,Nz))/2.*(Dgo_eff(i+1,j,Nz)+Dgo_eff(i,j,Nz))/2.)/dx* &
                         (ne(i+1,j,Nz)+ne(i,j,Nz))/2./dx - &
                         (alpha1_o(i,j-1,Nz)+alpha1_o(i,j,Nz))/2.*min(v(i,j,Nz),0.)*(ne(i,j-1,Nz)+ne(i,j,Nz))/2./dy + &
                         ((alpha1_o(i,j-1,Nz)+alpha1_o(i,j,Nz))/2.*(Dlo_eff(i,j-1,Nz)+Dlo_eff(i,j,Nz))/2. + &
                         (alpha2_o(i,j-1,Nz)+alpha2_o(i,j,Nz))/2.*(Dgo_eff(i,j-1,Nz)+Dgo_eff(i,j,Nz))/2.)/dy* &
                         (ne(i,j-1,Nz)+ne(i,j,Nz))/2./dy + &
                         (alpha1_o(i,j+1,Nz)+alpha1_o(i,j,Nz))/2.*max(v(i,j+1,Nz),0.)*(ne(i,j+1,Nz)+ne(i,j,Nz))/2./dy + &
                         ((alpha1_o(i,j+1,Nz)+alpha1_o(i,j,Nz))/2.*(Dlo_eff(i,j+1,Nz)+Dlo_eff(i,j,Nz))/2. + &
                         (alpha2_o(i,j+1,Nz)+alpha2_o(i,j,Nz))/2.*(Dgo_eff(i,j+1,Nz)+Dgo_eff(i,j,Nz))/2.)/dy* &
                         (ne(i,j+1,Nz)+ne(i,j,Nz))/2./dy - &
                         (alpha1_o(i,j,Nz-1)+alpha1_o(i,j,Nz))/2.*min(w(i,j,Nz),0.)*(ne(i,j,Nz-1)+ne(i,j,Nz))/2./dz + &
                         ((alpha1_o(i,j,Nz-1)+alpha1_o(i,j,Nz))/2.*(Dlo_eff(i,j,Nz-1)+Dlo_eff(i,j,Nz))/2. + &
                         (alpha2_o(i,j,Nz-1)+alpha2_o(i,j,Nz))/2.*(Dgo_eff(i,j,Nz-1)+Dgo_eff(i,j,Nz))/2.)/dz* &
                         (ne(i,j,Nz-1)+ne(i,j,Nz))/2./dz + &
                         alpha1_o(i,j,Nz)*max(w(i,j,Nz+1),0.)*ne(i,j,Nz)/dz) - &
                         dt*ne(i,j,Nz)*(-alpha1_o(i,j,Nz)*qo_0*b_o(i,j,Nz)* &
                         lc_o(i,j,Nz)/(Kc+lc_o(i,j,Nz))/(Ko+lo_o(i,j,Nz)))

                    q(i,j,Nz) = ne(i,j,Nz)*to_o(i,j,Nz)
                 endif


              endif
           enddo
        enddo


        !! for inside

        Do i=2,Nx-1
           Do j=2,Ny-1
              Do k=2,Nz-1
                 if(wall(i,j,k)) then
                    aa(i,j,k) = 0.
                    bb(i,j,k) = 0.
                    cc(i,j,k) = 0.
                    dd(i,j,k) = 0.
                    ee(i,j,k) = 0.
                    ff(i,j,k) = 0.
                    gg(i,j,k) = 0.
                    q(i,j,k) = 0.
                 else
                    aa(i,j,k) = ne(i,j,k) + dt*(-(alpha1_o(i-1,j,k)+alpha1_o(i,j,k))/2.*min(u(i,j,k),0.)* &
                         (ne(i-1,j,k)+ne(i,j,k))/2./dx + &
                         ((alpha1_o(i-1,j,k)+alpha1_o(i,j,k))/2.*(Dlo_eff(i-1,j,k)+Dlo_eff(i,j,k))/2. + &
                         (alpha2_o(i-1,j,k)+alpha2_o(i,j,k))/2.*(Dgo_eff(i-1,j,k)+Dgo_eff(i,j,k))/2.)/dx* &
                         (ne(i-1,j,k)+ne(i,j,k))/2./dx + &
                         (alpha1_o(i+1,j,k)+alpha1_o(i,j,k))/2.*max(u(i+1,j,k),0.)* &
                         (ne(i+1,j,k)+ne(i,j,k))/2./dx + &
                         ((alpha1_o(i+1,j,k)+alpha1_o(i,j,k))/2.*(Dlo_eff(i+1,j,k)+Dlo_eff(i,j,k))/2. + &
                         (alpha2_o(i+1,j,k)+alpha2_o(i,j,k))/2.*(Dgo_eff(i+1,j,k)+Dgo_eff(i,j,k))/2.)/dx* &
                         (ne(i+1,j,k)+ne(i,j,k))/2./dx - &
                         (alpha1_o(i,j-1,k)+alpha1_o(i,j,k))/2.*min(v(i,j,k),0.)*(ne(i,j-1,k)+ne(i,j,k))/2./dy + &
                         ((alpha1_o(i,j-1,k)+alpha1_o(i,j,k))/2.*(Dlo_eff(i,j-1,k)+Dlo_eff(i,j,k))/2. + &
                         (alpha2_o(i,j-1,k)+alpha2_o(i,j,k))/2.*(Dgo_eff(i,j-1,k)+Dgo_eff(i,j,k))/2.)/dy* &
                         (ne(i,j-1,k)+ne(i,j,k))/2./dy + &
                         (alpha1_o(i,j+1,k)+alpha1_o(i,j,k))/2.*max(v(i,j+1,k),0.)*(ne(i,j+1,k)+ne(i,j,k))/2./dy + &
                         ((alpha1_o(i,j+1,k)+alpha1_o(i,j,k))/2.*(Dlo_eff(i,j+1,k)+Dlo_eff(i,j,k))/2. + &
                         (alpha2_o(i,j+1,k)+alpha2_o(i,j,k))/2.*(Dgo_eff(i,j+1,k)+Dgo_eff(i,j,k))/2.)/dy* &
                         (ne(i,j+1,k)+ne(i,j,k))/2./dy - &
                         (alpha1_o(i,j,k-1)+alpha1_o(i,j,k))/2.*min(w(i,j,k),0.)*(ne(i,j,k-1)+ne(i,j,k))/2./dz + &
                         ((alpha1_o(i,j,k-1)+alpha1_o(i,j,k))/2.*(Dlo_eff(i,j,k-1)+Dlo_eff(i,j,k))/2. + &
                         (alpha2_o(i,j,k-1)+alpha2_o(i,j,k))/2.*(Dgo_eff(i,j,k-1)+Dgo_eff(i,j,k))/2.)/dz* &
                         (ne(i,j,k-1)+ne(i,j,k))/2./dz + &
                         (alpha1_o(i,j,k+1)+alpha1_o(i,j,k))/2.*max(w(i,j,k+1),0.)*(ne(i,j,k+1)+ne(i,j,k))/2./dz + &
                         ((alpha1_o(i,j,k+1)+alpha1_o(i,j,k))/2.*(Dlo_eff(i,j,k+1)+Dlo_eff(i,j,k))/2. + &
                         (alpha2_o(i,j,k+1)+alpha2_o(i,j,k))/2.*(Dgo_eff(i,j,k+1)+Dgo_eff(i,j,k))/2.)/dz* &
                         (ne(i,j,k+1)+ne(i,j,k))/2./dz) - &
                         dt*ne(i,j,k)*(-alpha1_o(i,j,k)*qo_0*b_o(i,j,k)* &
                         lc_o(i,j,k)/(Kc+lc_o(i,j,k))/(Ko+lo_o(i,j,k)))

                    bb(i,j,k) = -dt*((alpha1_o(i-1,j,k)+alpha1_o(i,j,k))/2.*max(u(i,j,k),0.) + &
                         ((alpha1_o(i-1,j,k)+alpha1_o(i,j,k))/2.*(Dlo_eff(i-1,j,k)+Dlo_eff(i,j,k))/2. + &
                         (alpha2_o(i-1,j,k)+alpha2_o(i,j,k))/2.*(Dgo_eff(i-1,j,k)+Dgo_eff(i,j,k))/2.)/dx)* &
                         (ne(i-1,j,k)+ne(i,j,k))/2./dx
                    cc(i,j,k) = dt*((alpha1_o(i+1,j,k)+alpha1_o(i,j,k))/2.*min(u(i+1,j,k),0.) - &
                         ((alpha1_o(i+1,j,k)+alpha1_o(i,j,k))/2.*(Dlo_eff(i+1,j,k)+Dlo_eff(i,j,k))/2. + &
                         (alpha2_o(i+1,j,k)+alpha2_o(i,j,k))/2.*(Dgo_eff(i+1,j,k)+Dgo_eff(i,j,k))/2.)/dx)* &
                         (ne(i+1,j,k)+ne(i,j,k))/2./dx
                    dd(i,j,k) = -dt*((alpha1_o(i,j-1,k)+alpha1_o(i,j,k))/2.*max(v(i,j,k),0.) + &
                         ((alpha1_o(i,j-1,k)+alpha1_o(i,j,k))/2.*(Dlo_eff(i,j-1,k)+Dlo_eff(i,j,k))/2. + &
                         (alpha2_o(i,j-1,k)+alpha2_o(i,j,k))/2.*(Dgo_eff(i,j-1,k)+Dgo_eff(i,j,k))/2.)/dy)* &
                         (ne(i,j-1,k)+ne(i,j,k))/2./dy
                    ee(i,j,k) = dt*((alpha1_o(i,j+1,k)+alpha1_o(i,j,k))/2.*min(v(i,j+1,k),0.) - &
                         ((alpha1_o(i,j+1,k)+alpha1_o(i,j,k))/2.*(Dlo_eff(i,j+1,k)+Dlo_eff(i,j,k))/2. + &
                         (alpha2_o(i,j+1,k)+alpha2_o(i,j,k))/2.*(Dgo_eff(i,j+1,k)+Dgo_eff(i,j,k))/2.)/dy)* &
                         (ne(i,j+1,k)+ne(i,j,k))/2./dy
                    ff(i,j,k) = -dt*((alpha1_o(i,j,k-1)+alpha1_o(i,j,k))/2.*max(w(i,j,k),0.) + &
                         ((alpha1_o(i,j,k-1)+alpha1_o(i,j,k))/2.*(Dlo_eff(i,j,k-1)+Dlo_eff(i,j,k))/2. + &
                         (alpha2_o(i,j,k-1)+alpha2_o(i,j,k))/2.*(Dgo_eff(i,j,k-1)+Dgo_eff(i,j,k))/2.)/dz)* &
                         (ne(i,j,k-1)+ne(i,j,k))/2./dz
                    gg(i,j,k) = dt*((alpha1_o(i,j,k+1)+alpha1_o(i,j,k))/2.*min(w(i,j,k+1),0.) - &
                         ((alpha1_o(i,j,k+1)+alpha1_o(i,j,k))/2.*(Dlo_eff(i,j,k+1)+Dlo_eff(i,j,k))/2. + &
                         (alpha2_o(i,j,k+1)+alpha2_o(i,j,k))/2.*(Dgo_eff(i,j,k+1)+Dgo_eff(i,j,k))/2.)/dz)* &
                         (ne(i,j,k+1)+ne(i,j,k))/2./dz

                    q(i,j,k) = ne(i,j,k)*to_o(i,j,k)
                 endif
              Enddo
           Enddo
        Enddo


     else 
        if(mype.eq.0) then

           Do i=2,Nx-1
              Do j=2,Ny-1
                 if(wall(i,j,1)) then
                    aa(i,j,1) = 0.
                    bb(i,j,1) = 0.
                    cc(i,j,1) = 0.
                    dd(i,j,1) = 0.
                    ee(i,j,1) = 0.
                    ff(i,j,1) = 0.
                    gg(i,j,1) = 0.
                    q(i,j,1) = 0.
                 else
                 bb(i,j,1) = -dt*((alpha1_o(i-1,j,1)+alpha1_o(i,j,1))/2.*max(u(i,j,1),0.) + &
                      ((alpha1_o(i-1,j,1)+alpha1_o(i,j,1))/2.*(Dlo_eff(i-1,j,1)+Dlo_eff(i,j,1))/2. + &
                      (alpha2_o(i-1,j,1)+alpha2_o(i,j,1))/2.*(Dgo_eff(i-1,j,1)+Dgo_eff(i,j,1))/2.)/dx)* &
                      (ne(i-1,j,1)+ne(i,j,1))/2./dx
                 cc(i,j,1) = dt*((alpha1_o(i+1,j,1)+alpha1_o(i,j,1))/2.*min(u(i+1,j,1),0.) - &
                      ((alpha1_o(i+1,j,1)+alpha1_o(i,j,1))/2.*(Dlo_eff(i+1,j,1)+Dlo_eff(i,j,1))/2. + &
                      (alpha2_o(i+1,j,1)+alpha2_o(i,j,1))/2.*(Dgo_eff(i+1,j,1)+Dgo_eff(i,j,1))/2.)/dx)* &
                      (ne(i+1,j,1)+ne(i,j,1))/2./dx
                 dd(i,j,1) = -dt*((alpha1_o(i,j-1,1)+alpha1_o(i,j,1))/2.*max(v(i,j,1),0.) + &
                      ((alpha1_o(i,j-1,1)+alpha1_o(i,j,1))/2.*(Dlo_eff(i,j-1,1)+Dlo_eff(i,j,1))/2. + &
                      (alpha2_o(i,j-1,1)+alpha2_o(i,j,1))/2.*(Dgo_eff(i,j-1,1)+Dgo_eff(i,j,1))/2.)/dy)* &
                      (ne(i,j-1,1)+ne(i,j,1))/2./dy
                 ee(i,j,1) = dt*((alpha1_o(i,j+1,1)+alpha1_o(i,j,1))/2.*min(v(i,j+1,1),0.) - &
                      ((alpha1_o(i,j+1,1)+alpha1_o(i,j,1))/2.*(Dlo_eff(i,j+1,1)+Dlo_eff(i,j,1))/2. + &
                      (alpha2_o(i,j+1,1)+alpha2_o(i,j,1))/2.*(Dgo_eff(i,j+1,1)+Dgo_eff(i,j,1))/2.)/dy)* &
                      (ne(i,j+1,1)+ne(i,j,1))/2./dy
                 ff(i,j,1) = 0.
                 gg(i,j,1) = dt*((alpha1_o(i,j,1+1)+alpha1_o(i,j,1))/2.*min(w(i,j,1+1),0.) - &
                      ((alpha1_o(i,j,1+1)+alpha1_o(i,j,1))/2.*(Dlo_eff(i,j,1+1)+Dlo_eff(i,j,1))/2. + &
                      (alpha2_o(i,j,1+1)+alpha2_o(i,j,1))/2.*(Dgo_eff(i,j,1+1)+Dgo_eff(i,j,1))/2.)/dz)* &
                      (ne(i,j,1+1)+ne(i,j,1))/2./dz


                 if(tobc_d.eq.1) then 
                    aa(i,j,1) = ne(i,j,1) + dt*(-(alpha1_o(i-1,j,1)+alpha1_o(i,j,1))/2.*min(u(i,j,1),0.)* &
                         (ne(i-1,j,1)+ne(i,j,1))/2./dx + &
                         ((alpha1_o(i-1,j,1)+alpha1_o(i,j,1))/2.*(Dlo_eff(i-1,j,1)+Dlo_eff(i,j,1))/2. + &
                         (alpha2_o(i-1,j,1)+alpha2_o(i,j,1))/2.*(Dgo_eff(i-1,j,1)+Dgo_eff(i,j,1))/2.)/dx* &
                         (ne(i-1,j,1)+ne(i,j,1))/2./dx + &
                         (alpha1_o(i+1,j,1)+alpha1_o(i,j,1))/2.*max(u(i+1,j,1),0.)* &
                         (ne(i+1,j,1)+ne(i,j,1))/2./dx + &
                         ((alpha1_o(i+1,j,1)+alpha1_o(i,j,1))/2.*(Dlo_eff(i+1,j,1)+Dlo_eff(i,j,1))/2. + &
                         (alpha2_o(i+1,j,1)+alpha2_o(i,j,1))/2.*(Dgo_eff(i+1,j,1)+Dgo_eff(i,j,1))/2.)/dx* &
                         (ne(i+1,j,1)+ne(i,j,1))/2./dx - &
                         (alpha1_o(i,j-1,1)+alpha1_o(i,j,1))/2.*min(v(i,j,1),0.)*(ne(i,j-1,1)+ne(i,j,1))/2./dy + &
                         ((alpha1_o(i,j-1,1)+alpha1_o(i,j,1))/2.*(Dlo_eff(i,j-1,1)+Dlo_eff(i,j,1))/2. + &
                         (alpha2_o(i,j-1,1)+alpha2_o(i,j,1))/2.*(Dgo_eff(i,j-1,1)+Dgo_eff(i,j,1))/2.)/dy* &
                         (ne(i,j-1,1)+ne(i,j,1))/2./dy + &
                         (alpha1_o(i,j+1,1)+alpha1_o(i,j,1))/2.*max(v(i,j+1,1),0.)*(ne(i,j+1,1)+ne(i,j,1))/2./dy + &
                         ((alpha1_o(i,j+1,1)+alpha1_o(i,j,1))/2.*(Dlo_eff(i,j+1,1)+Dlo_eff(i,j,1))/2. + &
                         (alpha2_o(i,j+1,1)+alpha2_o(i,j,1))/2.*(Dgo_eff(i,j+1,1)+Dgo_eff(i,j,1))/2.)/dy* &
                         (ne(i,j+1,1)+ne(i,j,1))/2./dy + &
                         (alpha1_o(i,j,1)*Dlo_eff(i,j,1)+alpha2_o(i,j,1)*Dgo_eff(i,j,1))/dz*2.*ne(i,j,1)/dz + &
                         (alpha1_o(i,j,1+1)+alpha1_o(i,j,1))/2.*max(w(i,j,1+1),0.)*(ne(i,j,1+1)+ne(i,j,1))/2./dz + &
                         ((alpha1_o(i,j,1+1)+alpha1_o(i,j,1))/2.*(Dlo_eff(i,j,1+1)+Dlo_eff(i,j,1))/2. + &
                         (alpha2_o(i,j,1+1)+alpha2_o(i,j,1))/2.*(Dgo_eff(i,j,1+1)+Dgo_eff(i,j,1))/2.)/dz* &
                         (ne(i,j,1+1)+ne(i,j,1))/2./dz) - &
                         dt*ne(i,j,1)*(-alpha1_o(i,j,1)*qo_0*b_o(i,j,1)* &
                         lc_o(i,j,1)/(Kc+lc_o(i,j,1))/(Ko+lo_o(i,j,1)))

                    q(i,j,1) = ne(i,j,1)*to_o(i,j,1) - &
                            (-dt*(alpha1_o(i,j,1)*max(w(i,j,1),0.) + &
                            2.*(alpha1_o(i,j,1)*Dlo_eff(i,j,1)+alpha2_o(i,j,1)*Dgo_eff(i,j,1))/dz)* &
                            ne(i,j,1)/dz)*to_d(i,j)
                 endif
                 if(tobc_d.eq.2) then 
                    aa(i,j,1) = ne(i,j,1) + dt*(-(alpha1_o(i-1,j,1)+alpha1_o(i,j,1))/2.*min(u(i,j,1),0.)* &
                         (ne(i-1,j,1)+ne(i,j,1))/2./dx + &
                         ((alpha1_o(i-1,j,1)+alpha1_o(i,j,1))/2.*(Dlo_eff(i-1,j,1)+Dlo_eff(i,j,1))/2. + &
                         (alpha2_o(i-1,j,1)+alpha2_o(i,j,1))/2.*(Dgo_eff(i-1,j,1)+Dgo_eff(i,j,1))/2.)/dx* &
                         (ne(i-1,j,1)+ne(i,j,1))/2./dx + &
                         (alpha1_o(i+1,j,1)+alpha1_o(i,j,1))/2.*max(u(i+1,j,1),0.)* &
                         (ne(i+1,j,1)+ne(i,j,1))/2./dx + &
                         ((alpha1_o(i+1,j,1)+alpha1_o(i,j,1))/2.*(Dlo_eff(i+1,j,1)+Dlo_eff(i,j,1))/2. + &
                         (alpha2_o(i+1,j,1)+alpha2_o(i,j,1))/2.*(Dgo_eff(i+1,j,1)+Dgo_eff(i,j,1))/2.)/dx* &
                         (ne(i+1,j,1)+ne(i,j,1))/2./dx - &
                         (alpha1_o(i,j-1,1)+alpha1_o(i,j,1))/2.*min(v(i,j,1),0.)*(ne(i,j-1,1)+ne(i,j,1))/2./dy + &
                         ((alpha1_o(i,j-1,1)+alpha1_o(i,j,1))/2.*(Dlo_eff(i,j-1,1)+Dlo_eff(i,j,1))/2. + &
                         (alpha2_o(i,j-1,1)+alpha2_o(i,j,1))/2.*(Dgo_eff(i,j-1,1)+Dgo_eff(i,j,1))/2.)/dy* &
                         (ne(i,j-1,1)+ne(i,j,1))/2./dy + &
                         (alpha1_o(i,j+1,1)+alpha1_o(i,j,1))/2.*max(v(i,j+1,1),0.)*(ne(i,j+1,1)+ne(i,j,1))/2./dy + &
                         ((alpha1_o(i,j+1,1)+alpha1_o(i,j,1))/2.*(Dlo_eff(i,j+1,1)+Dlo_eff(i,j,1))/2. + &
                         (alpha2_o(i,j+1,1)+alpha2_o(i,j,1))/2.*(Dgo_eff(i,j+1,1)+Dgo_eff(i,j,1))/2.)/dy* &
                         (ne(i,j+1,1)+ne(i,j,1))/2./dy - &
                         alpha1_o(i,j,1)*min(w(i,j,1),0.)*ne(i,j,1)/dz + &
                         (alpha1_o(i,j,1+1)+alpha1_o(i,j,1))/2.*max(w(i,j,1+1),0.)*(ne(i,j,1+1)+ne(i,j,1))/2./dz + &
                         ((alpha1_o(i,j,1+1)+alpha1_o(i,j,1))/2.*(Dlo_eff(i,j,1+1)+Dlo_eff(i,j,1))/2. + &
                         (alpha2_o(i,j,1+1)+alpha2_o(i,j,1))/2.*(Dgo_eff(i,j,1+1)+Dgo_eff(i,j,1))/2.)/dz* &
                         (ne(i,j,1+1)+ne(i,j,1))/2./dz) - &
                         dt*ne(i,j,1)*(-alpha1_o(i,j,1)*qo_0*b_o(i,j,1)* &
                         lc_o(i,j,1)/(Kc+lc_o(i,j,1))/(Ko+lo_o(i,j,1)))

                    q(i,j,1) = ne(i,j,1)*to_o(i,j,1)
                    endif

                 endif
              enddo
           enddo


           Do i=2,Nx-1
              Do j=2,Ny-1
                    if(wall(i,j,Nz)) then
                       aa(i,j,Nz) = 0.
                       bb(i,j,Nz) = 0.
                       cc(i,j,Nz) = 0.
                       dd(i,j,Nz) = 0.
                       ee(i,j,Nz) = 0.
                       ff(i,j,Nz) = 0.
                       gg(i,j,Nz) = 0.
                       q(i,j,Nz) = 0.
                    else
                       aa(i,j,Nz) = ne(i,j,Nz) + dt*(-(alpha1_o(i-1,j,Nz)+alpha1_o(i,j,Nz))/2.*min(u(i,j,Nz),0.)* &
                            (ne(i-1,j,Nz)+ne(i,j,Nz))/2./dx + &
                            ((alpha1_o(i-1,j,Nz)+alpha1_o(i,j,Nz))/2.*(Dlo_eff(i-1,j,Nz)+Dlo_eff(i,j,Nz))/2. + &
                            (alpha2_o(i-1,j,Nz)+alpha2_o(i,j,Nz))/2.*(Dgo_eff(i-1,j,Nz)+Dgo_eff(i,j,Nz))/2.)/dx* &
                            (ne(i-1,j,Nz)+ne(i,j,Nz))/2./dx + &
                            (alpha1_o(i+1,j,Nz)+alpha1_o(i,j,Nz))/2.*max(u(i+1,j,Nz),0.)* &
                            (ne(i+1,j,Nz)+ne(i,j,Nz))/2./dx + &
                            ((alpha1_o(i+1,j,Nz)+alpha1_o(i,j,Nz))/2.*(Dlo_eff(i+1,j,Nz)+Dlo_eff(i,j,Nz))/2. + &
                            (alpha2_o(i+1,j,Nz)+alpha2_o(i,j,Nz))/2.*(Dgo_eff(i+1,j,Nz)+Dgo_eff(i,j,Nz))/2.)/dx* &
                            (ne(i+1,j,Nz)+ne(i,j,Nz))/2./dx - &
                            (alpha1_o(i,j-1,Nz)+alpha1_o(i,j,Nz))/2.*min(v(i,j,Nz),0.)*(ne(i,j-1,Nz)+ne(i,j,Nz))/2./dy + &
                            ((alpha1_o(i,j-1,Nz)+alpha1_o(i,j,Nz))/2.*(Dlo_eff(i,j-1,Nz)+Dlo_eff(i,j,Nz))/2. + &
                            (alpha2_o(i,j-1,Nz)+alpha2_o(i,j,Nz))/2.*(Dgo_eff(i,j-1,Nz)+Dgo_eff(i,j,Nz))/2.)/dy* &
                            (ne(i,j-1,Nz)+ne(i,j,Nz))/2./dy + &
                            (alpha1_o(i,j+1,Nz)+alpha1_o(i,j,Nz))/2.*max(v(i,j+1,Nz),0.)*(ne(i,j+1,Nz)+ne(i,j,Nz))/2./dy + &
                            ((alpha1_o(i,j+1,Nz)+alpha1_o(i,j,Nz))/2.*(Dlo_eff(i,j+1,Nz)+Dlo_eff(i,j,Nz))/2. + &
                            (alpha2_o(i,j+1,Nz)+alpha2_o(i,j,Nz))/2.*(Dgo_eff(i,j+1,Nz)+Dgo_eff(i,j,Nz))/2.)/dy* &
                            (ne(i,j+1,Nz)+ne(i,j,Nz))/2./dy - &
                            (alpha1_o(i,j,Nz-1)+alpha1_o(i,j,Nz))/2.*min(w(i,j,Nz),0.)*(ne(i,j,Nz-1)+ne(i,j,Nz))/2./dz + &
                            ((alpha1_o(i,j,Nz-1)+alpha1_o(i,j,Nz))/2.*(Dlo_eff(i,j,Nz-1)+Dlo_eff(i,j,Nz))/2. + &
                            (alpha2_o(i,j,Nz-1)+alpha2_o(i,j,Nz))/2.*(Dgo_eff(i,j,Nz-1)+Dgo_eff(i,j,Nz))/2.)/dz* &
                            (ne(i,j,Nz-1)+ne(i,j,Nz))/2./dz + &
                            (tmp_alpha1_o(i,j,2)+alpha1_o(i,j,Nz))/2.*max(w(i,j,Nz+1),0.)*(tmp_ne(i,j,2)+ne(i,j,Nz))/2./dz + &
                            ((tmp_alpha1_o(i,j,2)+alpha1_o(i,j,Nz))/2.*(tmp_Dlo_eff(i,j,2)+Dlo_eff(i,j,Nz))/2. + &
                            (tmp_alpha2_o(i,j,2)+alpha2_o(i,j,Nz))/2.*(tmp_Dgo_eff(i,j,2)+Dgo_eff(i,j,Nz))/2.)/dz* &
                            (tmp_ne(i,j,2)+ne(i,j,Nz))/2./dz) - &
                            dt*ne(i,j,Nz)* &
                            (-alpha1_o(i,j,Nz)*qo_0*b_o(i,j,Nz)* &
                            lc_o(i,j,Nz)/(Kc+lc_o(i,j,Nz))/(Ko+lo_o(i,j,Nz)))

                       bb(i,j,Nz) = -dt*((alpha1_o(i-1,j,Nz)+alpha1_o(i,j,Nz))/2.*max(u(i,j,Nz),0.) + &
                            ((alpha1_o(i-1,j,Nz)+alpha1_o(i,j,Nz))/2.*(Dlo_eff(i-1,j,Nz)+Dlo_eff(i,j,Nz))/2. + &
                            (alpha2_o(i-1,j,Nz)+alpha2_o(i,j,Nz))/2.*(Dgo_eff(i-1,j,Nz)+Dgo_eff(i,j,Nz))/2.)/dx)* &
                            (ne(i-1,j,Nz)+ne(i,j,Nz))/2./dx
                       cc(i,j,Nz) = dt*((alpha1_o(i+1,j,Nz)+alpha1_o(i,j,Nz))/2.*min(u(i+1,j,Nz),0.) - &
                            ((alpha1_o(i+1,j,Nz)+alpha1_o(i,j,Nz))/2.*(Dlo_eff(i+1,j,Nz)+Dlo_eff(i,j,Nz))/2. + &
                            (alpha2_o(i+1,j,Nz)+alpha2_o(i,j,Nz))/2.*(Dgo_eff(i+1,j,Nz)+Dgo_eff(i,j,Nz))/2.)/dx)* &
                            (ne(i+1,j,Nz)+ne(i,j,Nz))/2./dx
                       dd(i,j,Nz) = -dt*((alpha1_o(i,j-1,Nz)+alpha1_o(i,j,Nz))/2.*max(v(i,j,Nz),0.) + &
                            ((alpha1_o(i,j-1,Nz)+alpha1_o(i,j,Nz))/2.*(Dlo_eff(i,j-1,Nz)+Dlo_eff(i,j,Nz))/2. + &
                            (alpha2_o(i,j-1,Nz)+alpha2_o(i,j,Nz))/2.*(Dgo_eff(i,j-1,Nz)+Dgo_eff(i,j,Nz))/2.)/dy)* &
                            (ne(i,j-1,Nz)+ne(i,j,Nz))/2./dy
                       ee(i,j,Nz) = dt*((alpha1_o(i,j+1,Nz)+alpha1_o(i,j,Nz))/2.*min(v(i,j+1,Nz),0.) - &
                            ((alpha1_o(i,j+1,Nz)+alpha1_o(i,j,Nz))/2.*(Dlo_eff(i,j+1,Nz)+Dlo_eff(i,j,Nz))/2. + &
                            (alpha2_o(i,j+1,Nz)+alpha2_o(i,j,Nz))/2.*(Dgo_eff(i,j+1,Nz)+Dgo_eff(i,j,Nz))/2.)/dy)* &
                            (ne(i,j+1,Nz)+ne(i,j,Nz))/2./dy
                       ff(i,j,Nz) = -dt*((alpha1_o(i,j,Nz-1)+alpha1_o(i,j,Nz))/2.*max(w(i,j,Nz),0.) + &
                             ((alpha1_o(i,j,Nz-1)+alpha1_o(i,j,Nz))/2.*(Dlo_eff(i,j,Nz-1)+Dlo_eff(i,j,Nz))/2. + &
                            (alpha2_o(i,j,Nz-1)+alpha2_o(i,j,Nz))/2.*(Dgo_eff(i,j,Nz-1)+Dgo_eff(i,j,Nz))/2.)/dz)* &
                            (ne(i,j,Nz-1)+ne(i,j,Nz))/2./dz
                       gg(i,j,Nz) = 0

                       q(i,j,Nz) = ne(i,j,Nz)*to_o(i,j,Nz) - &
                            (dt*((tmp_alpha1_o(i,j,2)+alpha1_o(i,j,Nz))/2.*min(w(i,j,Nz+1),0.) - &
                            ((tmp_alpha1_o(i,j,2)+alpha1_o(i,j,Nz))/2.*(tmp_Dlo_eff(i,j,2)+Dlo_eff(i,j,Nz))/2. + &
                            (tmp_alpha2_o(i,j,2)+alpha2_o(i,j,Nz))/2.*(tmp_Dgo_eff(i,j,2)+Dgo_eff(i,j,Nz))/2.)/dz)* &
                            (tmp_ne(i,j,2)+ne(i,j,Nz))/2./dz)*tmp_to_o(i,j,2)
                    endif
                 Enddo
              Enddo


           !! for inside

           Do i=2,Nx-1
              Do j=2,Ny-1
                 Do k=2,Nz-1
                    if(wall(i,j,k)) then
                       aa(i,j,k) = 0.
                       bb(i,j,k) = 0.
                       cc(i,j,k) = 0.
                       dd(i,j,k) = 0.
                       ee(i,j,k) = 0.
                       ff(i,j,k) = 0.
                       gg(i,j,k) = 0.
                       q(i,j,k) = 0.
                    else
                       aa(i,j,k) = ne(i,j,k) + dt*(-(alpha1_o(i-1,j,k)+alpha1_o(i,j,k))/2.*min(u(i,j,k),0.)* &
                            (ne(i-1,j,k)+ne(i,j,k))/2./dx + &
                            ((alpha1_o(i-1,j,k)+alpha1_o(i,j,k))/2.*(Dlo_eff(i-1,j,k)+Dlo_eff(i,j,k))/2. + &
                            (alpha2_o(i-1,j,k)+alpha2_o(i,j,k))/2.*(Dgo_eff(i-1,j,k)+Dgo_eff(i,j,k))/2.)/dx* &
                            (ne(i-1,j,k)+ne(i,j,k))/2./dx + &
                            (alpha1_o(i+1,j,k)+alpha1_o(i,j,k))/2.*max(u(i+1,j,k),0.)* &
                            (ne(i+1,j,k)+ne(i,j,k))/2./dx + &
                            ((alpha1_o(i+1,j,k)+alpha1_o(i,j,k))/2.*(Dlo_eff(i+1,j,k)+Dlo_eff(i,j,k))/2. + &
                            (alpha2_o(i+1,j,k)+alpha2_o(i,j,k))/2.*(Dgo_eff(i+1,j,k)+Dgo_eff(i,j,k))/2.)/dx* &
                            (ne(i+1,j,k)+ne(i,j,k))/2./dx - &
                            (alpha1_o(i,j-1,k)+alpha1_o(i,j,k))/2.*min(v(i,j,k),0.)*(ne(i,j-1,k)+ne(i,j,k))/2./dy + &
                            ((alpha1_o(i,j-1,k)+alpha1_o(i,j,k))/2.*(Dlo_eff(i,j-1,k)+Dlo_eff(i,j,k))/2. + &
                            (alpha2_o(i,j-1,k)+alpha2_o(i,j,k))/2.*(Dgo_eff(i,j-1,k)+Dgo_eff(i,j,k))/2.)/dy* &
                            (ne(i,j-1,k)+ne(i,j,k))/2./dy + &
                            (alpha1_o(i,j+1,k)+alpha1_o(i,j,k))/2.*max(v(i,j+1,k),0.)*(ne(i,j+1,k)+ne(i,j,k))/2./dy + &
                            ((alpha1_o(i,j+1,k)+alpha1_o(i,j,k))/2.*(Dlo_eff(i,j+1,k)+Dlo_eff(i,j,k))/2. + &
                            (alpha2_o(i,j+1,k)+alpha2_o(i,j,k))/2.*(Dgo_eff(i,j+1,k)+Dgo_eff(i,j,k))/2.)/dy* &
                            (ne(i,j+1,k)+ne(i,j,k))/2./dy - &
                            (alpha1_o(i,j,k-1)+alpha1_o(i,j,k))/2.*min(w(i,j,k),0.)*(ne(i,j,k-1)+ne(i,j,k))/2./dz + &
                            ((alpha1_o(i,j,k-1)+alpha1_o(i,j,k))/2.*(Dlo_eff(i,j,k-1)+Dlo_eff(i,j,k))/2. + &
                            (alpha2_o(i,j,k-1)+alpha2_o(i,j,k))/2.*(Dgo_eff(i,j,k-1)+Dgo_eff(i,j,k))/2.)/dz* &
                            (ne(i,j,k-1)+ne(i,j,k))/2./dz + &
                            (alpha1_o(i,j,k+1)+alpha1_o(i,j,k))/2.*max(w(i,j,k+1),0.)*(ne(i,j,k+1)+ne(i,j,k))/2./dz + &
                            ((alpha1_o(i,j,k+1)+alpha1_o(i,j,k))/2.*(Dlo_eff(i,j,k+1)+Dlo_eff(i,j,k))/2. + &
                            (alpha2_o(i,j,k+1)+alpha2_o(i,j,k))/2.*(Dgo_eff(i,j,k+1)+Dgo_eff(i,j,k))/2.)/dz* &
                            (ne(i,j,k+1)+ne(i,j,k))/2./dz) - &
                            dt*ne(i,j,k)*(-alpha1_o(i,j,k)*qo_0*b_o(i,j,k)* &
                            lc_o(i,j,k)/(Kc+lc_o(i,j,k))/(Ko+lo_o(i,j,k)))

                       bb(i,j,k) = -dt*((alpha1_o(i-1,j,k)+alpha1_o(i,j,k))/2.*max(u(i,j,k),0.) + &
                            ((alpha1_o(i-1,j,k)+alpha1_o(i,j,k))/2.*(Dlo_eff(i-1,j,k)+Dlo_eff(i,j,k))/2. + &
                            (alpha2_o(i-1,j,k)+alpha2_o(i,j,k))/2.*(Dgo_eff(i-1,j,k)+Dgo_eff(i,j,k))/2.)/dx)* &
                            (ne(i-1,j,k)+ne(i,j,k))/2./dx
                       cc(i,j,k) = dt*((alpha1_o(i+1,j,k)+alpha1_o(i,j,k))/2.*min(u(i+1,j,k),0.) - &
                            ((alpha1_o(i+1,j,k)+alpha1_o(i,j,k))/2.*(Dlo_eff(i+1,j,k)+Dlo_eff(i,j,k))/2. + &
                            (alpha2_o(i+1,j,k)+alpha2_o(i,j,k))/2.*(Dgo_eff(i+1,j,k)+Dgo_eff(i,j,k))/2.)/dx)* &
                            (ne(i+1,j,k)+ne(i,j,k))/2./dx
                       dd(i,j,k) = -dt*((alpha1_o(i,j-1,k)+alpha1_o(i,j,k))/2.*max(v(i,j,k),0.) + &
                            ((alpha1_o(i,j-1,k)+alpha1_o(i,j,k))/2.*(Dlo_eff(i,j-1,k)+Dlo_eff(i,j,k))/2. + &
                            (alpha2_o(i,j-1,k)+alpha2_o(i,j,k))/2.*(Dgo_eff(i,j-1,k)+Dgo_eff(i,j,k))/2.)/dy)* &
                            (ne(i,j-1,k)+ne(i,j,k))/2./dy
                       ee(i,j,k) = dt*((alpha1_o(i,j+1,k)+alpha1_o(i,j,k))/2.*min(v(i,j+1,k),0.) - &
                            ((alpha1_o(i,j+1,k)+alpha1_o(i,j,k))/2.*(Dlo_eff(i,j+1,k)+Dlo_eff(i,j,k))/2. + &
                            (alpha2_o(i,j+1,k)+alpha2_o(i,j,k))/2.*(Dgo_eff(i,j+1,k)+Dgo_eff(i,j,k))/2.)/dy)* &
                            (ne(i,j+1,k)+ne(i,j,k))/2./dy
                       ff(i,j,k) = -dt*((alpha1_o(i,j,k-1)+alpha1_o(i,j,k))/2.*max(w(i,j,k),0.) + &
                             ((alpha1_o(i,j,k-1)+alpha1_o(i,j,k))/2.*(Dlo_eff(i,j,k-1)+Dlo_eff(i,j,k))/2. + &
                            (alpha2_o(i,j,k-1)+alpha2_o(i,j,k))/2.*(Dgo_eff(i,j,k-1)+Dgo_eff(i,j,k))/2.)/dz)* &
                            (ne(i,j,k-1)+ne(i,j,k))/2./dz
                       gg(i,j,k) = dt*((alpha1_o(i,j,k+1)+alpha1_o(i,j,k))/2.*min(w(i,j,k+1),0.) - &
                            ((alpha1_o(i,j,k+1)+alpha1_o(i,j,k))/2.*(Dlo_eff(i,j,k+1)+Dlo_eff(i,j,k))/2. + &
                            (alpha2_o(i,j,k+1)+alpha2_o(i,j,k))/2.*(Dgo_eff(i,j,k+1)+Dgo_eff(i,j,k))/2.)/dz)* &
                            (ne(i,j,k+1)+ne(i,j,k))/2./dz

                       q(i,j,k) = ne(i,j,k)*to_o(i,j,k)
                    endif
                 Enddo
              Enddo
           Enddo



        else if(mype.eq.npes-1) then


           Do i=2,Nx-1
              Do j=2,Ny-1
                    if(wall(i,j,1)) then
                       aa(i,j,1) = 0.
                       bb(i,j,1) = 0.
                       cc(i,j,1) = 0.
                       dd(i,j,1) = 0.
                       ee(i,j,1) = 0.
                       ff(i,j,1) = 0.
                       gg(i,j,1) = 0.
                       q(i,j,1) = 0.
                    else
                       aa(i,j,1) = ne(i,j,1) + dt*(-(alpha1_o(i-1,j,1)+alpha1_o(i,j,1))/2.*min(u(i,j,1),0.)* &
                            (ne(i-1,j,1)+ne(i,j,1))/2./dx + &
                            ((alpha1_o(i-1,j,1)+alpha1_o(i,j,1))/2.*(Dlo_eff(i-1,j,1)+Dlo_eff(i,j,1))/2. + &
                            (alpha2_o(i-1,j,1)+alpha2_o(i,j,1))/2.*(Dgo_eff(i-1,j,1)+Dgo_eff(i,j,1))/2.)/dx* &
                            (ne(i-1,j,1)+ne(i,j,1))/2./dx + &
                            (alpha1_o(i+1,j,1)+alpha1_o(i,j,1))/2.*max(u(i+1,j,1),0.)* &
                            (ne(i+1,j,1)+ne(i,j,1))/2./dx + &
                            ((alpha1_o(i+1,j,1)+alpha1_o(i,j,1))/2.*(Dlo_eff(i+1,j,1)+Dlo_eff(i,j,1))/2. + &
                            (alpha2_o(i+1,j,1)+alpha2_o(i,j,1))/2.*(Dgo_eff(i+1,j,1)+Dgo_eff(i,j,1))/2.)/dx* &
                            (ne(i+1,j,1)+ne(i,j,1))/2./dx - &
                            (alpha1_o(i,j-1,1)+alpha1_o(i,j,1))/2.*min(v(i,j,1),0.)*(ne(i,j-1,1)+ne(i,j,1))/2./dy + &
                            ((alpha1_o(i,j-1,1)+alpha1_o(i,j,1))/2.*(Dlo_eff(i,j-1,1)+Dlo_eff(i,j,1))/2. + &
                            (alpha2_o(i,j-1,1)+alpha2_o(i,j,1))/2.*(Dgo_eff(i,j-1,1)+Dgo_eff(i,j,1))/2.)/dy* &
                            (ne(i,j-1,1)+ne(i,j,1))/2./dy + &
                            (alpha1_o(i,j+1,1)+alpha1_o(i,j,1))/2.*max(v(i,j+1,1),0.)*(ne(i,j+1,1)+ne(i,j,1))/2./dy + &
                            ((alpha1_o(i,j+1,1)+alpha1_o(i,j,1))/2.*(Dlo_eff(i,j+1,1)+Dlo_eff(i,j,1))/2. + &
                            (alpha2_o(i,j+1,1)+alpha2_o(i,j,1))/2.*(Dgo_eff(i,j+1,1)+Dgo_eff(i,j,1))/2.)/dy* &
                            (ne(i,j+1,1)+ne(i,j,1))/2./dy - &
                            (tmp_alpha1_o(i,j,1)+alpha1_o(i,j,1))/2.*min(w(i,j,1),0.)*(tmp_ne(i,j,1)+ne(i,j,1))/2./dz + &
                            ((tmp_alpha1_o(i,j,1)+alpha1_o(i,j,1))/2.*(tmp_Dlo_eff(i,j,1)+Dlo_eff(i,j,1))/2. + &
                            (tmp_alpha2_o(i,j,1)+alpha2_o(i,j,1))/2.*(tmp_Dgo_eff(i,j,1)+Dgo_eff(i,j,1))/2.)/dz* &
                            (tmp_ne(i,j,1)+ne(i,j,1))/2./dz + &
                            (alpha1_o(i,j,1+1)+alpha1_o(i,j,1))/2.*max(w(i,j,1+1),0.)*(ne(i,j,1+1)+ne(i,j,1))/2./dz + &
                            ((alpha1_o(i,j,1+1)+alpha1_o(i,j,1))/2.*(Dlo_eff(i,j,1+1)+Dlo_eff(i,j,1))/2. + &
                            (alpha2_o(i,j,1+1)+alpha2_o(i,j,1))/2.*(Dgo_eff(i,j,1+1)+Dgo_eff(i,j,1))/2.)/dz* &
                            (ne(i,j,1+1)+ne(i,j,1))/2./dz) - &
                            dt*ne(i,j,1)*(-alpha1_o(i,j,1)*qo_0*b_o(i,j,1)* &
                            lc_o(i,j,1)/(Kc+lc_o(i,j,1))/(Ko+lo_o(i,j,1)))

                       bb(i,j,1) = -dt*((alpha1_o(i-1,j,1)+alpha1_o(i,j,1))/2.*max(u(i,j,1),0.) + &
                            ((alpha1_o(i-1,j,1)+alpha1_o(i,j,1))/2.*(Dlo_eff(i-1,j,1)+Dlo_eff(i,j,1))/2. + &
                            (alpha2_o(i-1,j,1)+alpha2_o(i,j,1))/2.*(Dgo_eff(i-1,j,1)+Dgo_eff(i,j,1))/2.)/dx)* &
                            (ne(i-1,j,1)+ne(i,j,1))/2./dx
                       cc(i,j,1) = dt*((alpha1_o(i+1,j,1)+alpha1_o(i,j,1))/2.*min(u(i+1,j,1),0.) - &
                            ((alpha1_o(i+1,j,1)+alpha1_o(i,j,1))/2.*(Dlo_eff(i+1,j,1)+Dlo_eff(i,j,1))/2. + &
                            (alpha2_o(i+1,j,1)+alpha2_o(i,j,1))/2.*(Dgo_eff(i+1,j,1)+Dgo_eff(i,j,1))/2.)/dx)* &
                            (ne(i+1,j,1)+ne(i,j,1))/2./dx
                       dd(i,j,1) = -dt*((alpha1_o(i,j-1,1)+alpha1_o(i,j,1))/2.*max(v(i,j,1),0.) + &
                            ((alpha1_o(i,j-1,1)+alpha1_o(i,j,1))/2.*(Dlo_eff(i,j-1,1)+Dlo_eff(i,j,1))/2. + &
                            (alpha2_o(i,j-1,1)+alpha2_o(i,j,1))/2.*(Dgo_eff(i,j-1,1)+Dgo_eff(i,j,1))/2.)/dy)* &
                            (ne(i,j-1,1)+ne(i,j,1))/2./dy
                       ee(i,j,1) = dt*((alpha1_o(i,j+1,1)+alpha1_o(i,j,1))/2.*min(v(i,j+1,1),0.) - &
                            ((alpha1_o(i,j+1,1)+alpha1_o(i,j,1))/2.*(Dlo_eff(i,j+1,1)+Dlo_eff(i,j,1))/2. + &
                            (alpha2_o(i,j+1,1)+alpha2_o(i,j,1))/2.*(Dgo_eff(i,j+1,1)+Dgo_eff(i,j,1))/2.)/dy)* &
                            (ne(i,j+1,1)+ne(i,j,1))/2./dy
                       ff(i,j,1) = 0.
                       gg(i,j,1) = dt*((alpha1_o(i,j,1+1)+alpha1_o(i,j,1))/2.*min(w(i,j,1+1),0.) - &
                            ((alpha1_o(i,j,1+1)+alpha1_o(i,j,1))/2.*(Dlo_eff(i,j,1+1)+Dlo_eff(i,j,1))/2. + &
                            (alpha2_o(i,j,1+1)+alpha2_o(i,j,1))/2.*(Dgo_eff(i,j,1+1)+Dgo_eff(i,j,1))/2.)/dz)* &
                            (ne(i,j,1+1)+ne(i,j,1))/2./dz

                       q(i,j,1) = ne(i,j,1)*to_o(i,j,1) - &
                            (-dt*((tmp_alpha1_o(i,j,1)+alpha1_o(i,j,1))/2.*max(w(i,j,1),0.) + &
                             ((tmp_alpha1_o(i,j,1)+alpha1_o(i,j,1))/2.*(tmp_Dlo_eff(i,j,1)+Dlo_eff(i,j,1))/2. + &
                            (tmp_alpha2_o(i,j,1)+alpha2_o(i,j,1))/2.*(tmp_Dgo_eff(i,j,1)+Dgo_eff(i,j,1))/2.)/dz)* &
                            (tmp_ne(i,j,1)+ne(i,j,1))/2./dz)*tmp_to_o(i,j,1)
                    endif
              Enddo
           Enddo
           
           
           
           Do i=2,Nx-1
              Do j=2,Ny-1
                 if(wall(i,j,Nz)) then
                    aa(i,j,Nz) = 0.
                    bb(i,j,Nz) = 0.
                    cc(i,j,Nz) = 0.
                    dd(i,j,Nz) = 0.
                    ee(i,j,Nz) = 0.
                    ff(i,j,Nz) = 0.
                    gg(i,j,Nz) = 0.
                    q(i,j,Nz) = 0.
                 else     
                 bb(i,j,Nz) = -dt*((alpha1_o(i-1,j,Nz)+alpha1_o(i,j,Nz))/2.*max(u(i,j,Nz),0.) + &
                      ((alpha1_o(i-1,j,Nz)+alpha1_o(i,j,Nz))/2.*(Dlo_eff(i-1,j,Nz)+Dlo_eff(i,j,Nz))/2. + &
                      (alpha2_o(i-1,j,Nz)+alpha2_o(i,j,Nz))/2.*(Dgo_eff(i-1,j,Nz)+Dgo_eff(i,j,Nz))/2.)/dx)* &
                      (ne(i-1,j,Nz)+ne(i,j,Nz))/2./dx
                 cc(i,j,Nz) = dt*((alpha1_o(i+1,j,Nz)+alpha1_o(i,j,Nz))/2.*min(u(i+1,j,Nz),0.) - &
                      ((alpha1_o(i+1,j,Nz)+alpha1_o(i,j,Nz))/2.*(Dlo_eff(i+1,j,Nz)+Dlo_eff(i,j,Nz))/2. + &
                      (alpha2_o(i+1,j,Nz)+alpha2_o(i,j,Nz))/2.*(Dgo_eff(i+1,j,Nz)+Dgo_eff(i,j,Nz))/2.)/dx)* &
                      (ne(i+1,j,Nz)+ne(i,j,Nz))/2./dx
                 dd(i,j,Nz) = -dt*((alpha1_o(i,j-1,Nz)+alpha1_o(i,j,Nz))/2.*max(v(i,j,Nz),0.) + &
                      ((alpha1_o(i,j-1,Nz)+alpha1_o(i,j,Nz))/2.*(Dlo_eff(i,j-1,Nz)+Dlo_eff(i,j,Nz))/2. + &
                      (alpha2_o(i,j-1,Nz)+alpha2_o(i,j,Nz))/2.*(Dgo_eff(i,j-1,Nz)+Dgo_eff(i,j,Nz))/2.)/dy)* &
                      (ne(i,j-1,Nz)+ne(i,j,Nz))/2./dy
                 ee(i,j,Nz) = dt*((alpha1_o(i,j+1,Nz)+alpha1_o(i,j,Nz))/2.*min(v(i,j+1,Nz),0.) - &
                      ((alpha1_o(i,j+1,Nz)+alpha1_o(i,j,Nz))/2.*(Dlo_eff(i,j+1,Nz)+Dlo_eff(i,j,Nz))/2. + &
                      (alpha2_o(i,j+1,Nz)+alpha2_o(i,j,Nz))/2.*(Dgo_eff(i,j+1,Nz)+Dgo_eff(i,j,Nz))/2.)/dy)* &
                      (ne(i,j+1,Nz)+ne(i,j,Nz))/2./dy
                 ff(i,j,Nz) = -dt*((alpha1_o(i,j,Nz-1)+alpha1_o(i,j,Nz))/2.*max(w(i,j,Nz),0.) + &
                      ((alpha1_o(i,j,Nz-1)+alpha1_o(i,j,Nz))/2.*(Dlo_eff(i,j,Nz-1)+Dlo_eff(i,j,Nz))/2. + &
                      (alpha2_o(i,j,Nz-1)+alpha2_o(i,j,Nz))/2.*(Dgo_eff(i,j,Nz-1)+Dgo_eff(i,j,Nz))/2.)/dz)* &
                      (ne(i,j,Nz-1)+ne(i,j,Nz))/2./dz
                 gg(i,j,Nz) = 0.

                 if(tobc_t.eq.1) then
                    aa(i,j,Nz) = ne(i,j,Nz) + dt*(-(alpha1_o(i-1,j,Nz)+alpha1_o(i,j,Nz))/2.*min(u(i,j,Nz),0.)* &
                         (ne(i-1,j,Nz)+ne(i,j,Nz))/2./dx + &
                         ((alpha1_o(i-1,j,Nz)+alpha1_o(i,j,Nz))/2.*(Dlo_eff(i-1,j,Nz)+Dlo_eff(i,j,Nz))/2. + &
                         (alpha2_o(i-1,j,Nz)+alpha2_o(i,j,Nz))/2.*(Dgo_eff(i-1,j,Nz)+Dgo_eff(i,j,Nz))/2.)/dx* &
                         (ne(i-1,j,Nz)+ne(i,j,Nz))/2./dx + &
                         (alpha1_o(i+1,j,Nz)+alpha1_o(i,j,Nz))/2.*max(u(i+1,j,Nz),0.)* &
                         (ne(i+1,j,Nz)+ne(i,j,Nz))/2./dx + &
                         ((alpha1_o(i+1,j,Nz)+alpha1_o(i,j,Nz))/2.*(Dlo_eff(i+1,j,Nz)+Dlo_eff(i,j,Nz))/2. + &
                         (alpha2_o(i+1,j,Nz)+alpha2_o(i,j,Nz))/2.*(Dgo_eff(i+1,j,Nz)+Dgo_eff(i,j,Nz))/2.)/dx* &
                         (ne(i+1,j,Nz)+ne(i,j,Nz))/2./dx - &
                         (alpha1_o(i,j-1,Nz)+alpha1_o(i,j,Nz))/2.*min(v(i,j,Nz),0.)*(ne(i,j-1,Nz)+ne(i,j,Nz))/2./dy + &
                         ((alpha1_o(i,j-1,Nz)+alpha1_o(i,j,Nz))/2.*(Dlo_eff(i,j-1,Nz)+Dlo_eff(i,j,Nz))/2. + &
                         (alpha2_o(i,j-1,Nz)+alpha2_o(i,j,Nz))/2.*(Dgo_eff(i,j-1,Nz)+Dgo_eff(i,j,Nz))/2.)/dy* &
                         (ne(i,j-1,Nz)+ne(i,j,Nz))/2./dy + &
                         (alpha1_o(i,j+1,Nz)+alpha1_o(i,j,Nz))/2.*max(v(i,j+1,Nz),0.)*(ne(i,j+1,Nz)+ne(i,j,Nz))/2./dy + &
                         ((alpha1_o(i,j+1,Nz)+alpha1_o(i,j,Nz))/2.*(Dlo_eff(i,j+1,Nz)+Dlo_eff(i,j,Nz))/2. + &
                         (alpha2_o(i,j+1,Nz)+alpha2_o(i,j,Nz))/2.*(Dgo_eff(i,j+1,Nz)+Dgo_eff(i,j,Nz))/2.)/dy* &
                         (ne(i,j+1,Nz)+ne(i,j,Nz))/2./dy - &
                         (alpha1_o(i,j,Nz-1)+alpha1_o(i,j,Nz))/2.*min(w(i,j,Nz),0.)*(ne(i,j,Nz-1)+ne(i,j,Nz))/2./dz + &
                         ((alpha1_o(i,j,Nz-1)+alpha1_o(i,j,Nz))/2.*(Dlo_eff(i,j,Nz-1)+Dlo_eff(i,j,Nz))/2. + &
                         (alpha2_o(i,j,Nz-1)+alpha2_o(i,j,Nz))/2.*(Dgo_eff(i,j,Nz-1)+Dgo_eff(i,j,Nz))/2.)/dz* &
                         (ne(i,j,Nz-1)+ne(i,j,Nz))/2./dz + &
                         (alpha1_o(i,j,Nz)*Dlo_eff(i,j,Nz) + alpha2_o(i,j,Nz)*Dgo_eff(i,j,Nz))/dz*2.*ne(i,j,Nz)/dz) - &
                         dt*ne(i,j,Nz)*(-alpha1_o(i,j,Nz)*qo_0*b_o(i,j,Nz)* &
                         lc_o(i,j,Nz)/(Kc+lc_o(i,j,Nz))/(Ko+lo_o(i,j,Nz)))

                    q(i,j,Nz) = ne(i,j,Nz)*to_o(i,j,Nz) - dt * (alpha1_o(i,j,Nz)*min(w(i,j,Nz+1),0.) - &
                            2.*(alpha1_o(i,j,Nz)*Dlo_eff(i,j,Nz)+alpha2_o(i,j,Nz)*Dgo_eff(i,j,Nz))/dz)*ne(i,j,Nz)/dz*to_t(i,j)
                 endif
                 if(tobc_t.eq.2) then
                    aa(i,j,Nz) = ne(i,j,Nz) + dt*(-(alpha1_o(i-1,j,Nz)+alpha1_o(i,j,Nz))/2.*min(u(i,j,Nz),0.)* &
                         (ne(i-1,j,Nz)+ne(i,j,Nz))/2./dx + &
                         ((alpha1_o(i-1,j,Nz)+alpha1_o(i,j,Nz))/2.*(Dlo_eff(i-1,j,Nz)+Dlo_eff(i,j,Nz))/2. + &
                         (alpha2_o(i-1,j,Nz)+alpha2_o(i,j,Nz))/2.*(Dgo_eff(i-1,j,Nz)+Dgo_eff(i,j,Nz))/2.)/dx* &
                         (ne(i-1,j,Nz)+ne(i,j,Nz))/2./dx + &
                         (alpha1_o(i+1,j,Nz)+alpha1_o(i,j,Nz))/2.*max(u(i+1,j,Nz),0.)* &
                         (ne(i+1,j,Nz)+ne(i,j,Nz))/2./dx + &
                         ((alpha1_o(i+1,j,Nz)+alpha1_o(i,j,Nz))/2.*(Dlo_eff(i+1,j,Nz)+Dlo_eff(i,j,Nz))/2. + &
                         (alpha2_o(i+1,j,Nz)+alpha2_o(i,j,Nz))/2.*(Dgo_eff(i+1,j,Nz)+Dgo_eff(i,j,Nz))/2.)/dx* &
                         (ne(i+1,j,Nz)+ne(i,j,Nz))/2./dx - &
                         (alpha1_o(i,j-1,Nz)+alpha1_o(i,j,Nz))/2.*min(v(i,j,Nz),0.)*(ne(i,j-1,Nz)+ne(i,j,Nz))/2./dy + &
                         ((alpha1_o(i,j-1,Nz)+alpha1_o(i,j,Nz))/2.*(Dlo_eff(i,j-1,Nz)+Dlo_eff(i,j,Nz))/2. + &
                         (alpha2_o(i,j-1,Nz)+alpha2_o(i,j,Nz))/2.*(Dgo_eff(i,j-1,Nz)+Dgo_eff(i,j,Nz))/2.)/dy* &
                         (ne(i,j-1,Nz)+ne(i,j,Nz))/2./dy + &
                         (alpha1_o(i,j+1,Nz)+alpha1_o(i,j,Nz))/2.*max(v(i,j+1,Nz),0.)*(ne(i,j+1,Nz)+ne(i,j,Nz))/2./dy + &
                         ((alpha1_o(i,j+1,Nz)+alpha1_o(i,j,Nz))/2.*(Dlo_eff(i,j+1,Nz)+Dlo_eff(i,j,Nz))/2. + &
                         (alpha2_o(i,j+1,Nz)+alpha2_o(i,j,Nz))/2.*(Dgo_eff(i,j+1,Nz)+Dgo_eff(i,j,Nz))/2.)/dy* &
                         (ne(i,j+1,Nz)+ne(i,j,Nz))/2./dy - &
                         (alpha1_o(i,j,Nz-1)+alpha1_o(i,j,Nz))/2.*min(w(i,j,Nz),0.)*(ne(i,j,Nz-1)+ne(i,j,Nz))/2./dz + &
                         ((alpha1_o(i,j,Nz-1)+alpha1_o(i,j,Nz))/2.*(Dlo_eff(i,j,Nz-1)+Dlo_eff(i,j,Nz))/2. + &
                         (alpha2_o(i,j,Nz-1)+alpha2_o(i,j,Nz))/2.*(Dgo_eff(i,j,Nz-1)+Dgo_eff(i,j,Nz))/2.)/dz* &
                         (ne(i,j,Nz-1)+ne(i,j,Nz))/2./dz + &
                         alpha1_o(i,j,Nz)*max(w(i,j,Nz+1),0.)*ne(i,j,Nz)/dz) - &
                         dt*ne(i,j,Nz)*(-alpha1_o(i,j,Nz)*qo_0*b_o(i,j,Nz)* &
                         lc_o(i,j,Nz)/(Kc+lc_o(i,j,Nz))/(Ko+lo_o(i,j,Nz)))

                    q(i,j,Nz) = ne(i,j,Nz)*to_o(i,j,Nz)
                    endif
                 endif
              enddo
           enddo



           !! for inside

           Do i=2,Nx-1
              Do j=2,Ny-1
                 Do k=2,Nz-1
                    if(wall(i,j,k)) then
                       aa(i,j,k) = 0.
                       bb(i,j,k) = 0.
                       cc(i,j,k) = 0.
                       dd(i,j,k) = 0.
                       ee(i,j,k) = 0.
                       ff(i,j,k) = 0.
                       gg(i,j,k) = 0.
                       q(i,j,k) = 0.
                    else
                    aa(i,j,k) = ne(i,j,k) + dt*(-(alpha1_o(i-1,j,k)+alpha1_o(i,j,k))/2.*min(u(i,j,k),0.)* &
                         (ne(i-1,j,k)+ne(i,j,k))/2./dx + &
                         ((alpha1_o(i-1,j,k)+alpha1_o(i,j,k))/2.*(Dlo_eff(i-1,j,k)+Dlo_eff(i,j,k))/2. + &
                         (alpha2_o(i-1,j,k)+alpha2_o(i,j,k))/2.*(Dgo_eff(i-1,j,k)+Dgo_eff(i,j,k))/2.)/dx* &
                         (ne(i-1,j,k)+ne(i,j,k))/2./dx + &
                         (alpha1_o(i+1,j,k)+alpha1_o(i,j,k))/2.*max(u(i+1,j,k),0.)* &
                         (ne(i+1,j,k)+ne(i,j,k))/2./dx + &
                         ((alpha1_o(i+1,j,k)+alpha1_o(i,j,k))/2.*(Dlo_eff(i+1,j,k)+Dlo_eff(i,j,k))/2. + &
                         (alpha2_o(i+1,j,k)+alpha2_o(i,j,k))/2.*(Dgo_eff(i+1,j,k)+Dgo_eff(i,j,k))/2.)/dx* &
                         (ne(i+1,j,k)+ne(i,j,k))/2./dx - &
                         (alpha1_o(i,j-1,k)+alpha1_o(i,j,k))/2.*min(v(i,j,k),0.)*(ne(i,j-1,k)+ne(i,j,k))/2./dy + &
                         ((alpha1_o(i,j-1,k)+alpha1_o(i,j,k))/2.*(Dlo_eff(i,j-1,k)+Dlo_eff(i,j,k))/2. + &
                         (alpha2_o(i,j-1,k)+alpha2_o(i,j,k))/2.*(Dgo_eff(i,j-1,k)+Dgo_eff(i,j,k))/2.)/dy* &
                         (ne(i,j-1,k)+ne(i,j,k))/2./dy + &
                         (alpha1_o(i,j+1,k)+alpha1_o(i,j,k))/2.*max(v(i,j+1,k),0.)*(ne(i,j+1,k)+ne(i,j,k))/2./dy + &
                         ((alpha1_o(i,j+1,k)+alpha1_o(i,j,k))/2.*(Dlo_eff(i,j+1,k)+Dlo_eff(i,j,k))/2. + &
                         (alpha2_o(i,j+1,k)+alpha2_o(i,j,k))/2.*(Dgo_eff(i,j+1,k)+Dgo_eff(i,j,k))/2.)/dy* &
                         (ne(i,j+1,k)+ne(i,j,k))/2./dy - &
                         (alpha1_o(i,j,k-1)+alpha1_o(i,j,k))/2.*min(w(i,j,k),0.)*(ne(i,j,k-1)+ne(i,j,k))/2./dz + &
                         ((alpha1_o(i,j,k-1)+alpha1_o(i,j,k))/2.*(Dlo_eff(i,j,k-1)+Dlo_eff(i,j,k))/2. + &
                         (alpha2_o(i,j,k-1)+alpha2_o(i,j,k))/2.*(Dgo_eff(i,j,k-1)+Dgo_eff(i,j,k))/2.)/dz* &
                         (ne(i,j,k-1)+ne(i,j,k))/2./dz + &
                         (alpha1_o(i,j,k+1)+alpha1_o(i,j,k))/2.*max(w(i,j,k+1),0.)*(ne(i,j,k+1)+ne(i,j,k))/2./dz + &
                         ((alpha1_o(i,j,k+1)+alpha1_o(i,j,k))/2.*(Dlo_eff(i,j,k+1)+Dlo_eff(i,j,k))/2. + &
                         (alpha2_o(i,j,k+1)+alpha2_o(i,j,k))/2.*(Dgo_eff(i,j,k+1)+Dgo_eff(i,j,k))/2.)/dz* &
                         (ne(i,j,k+1)+ne(i,j,k))/2./dz) - &
                         dt*ne(i,j,k)*(-alpha1_o(i,j,k)*qo_0*b_o(i,j,k)* &
                         lc_o(i,j,k)/(Kc+lc_o(i,j,k))/(Ko+lo_o(i,j,k)))

                    bb(i,j,k) = -dt*((alpha1_o(i-1,j,k)+alpha1_o(i,j,k))/2.*max(u(i,j,k),0.) + &
                         ((alpha1_o(i-1,j,k)+alpha1_o(i,j,k))/2.*(Dlo_eff(i-1,j,k)+Dlo_eff(i,j,k))/2. + &
                         (alpha2_o(i-1,j,k)+alpha2_o(i,j,k))/2.*(Dgo_eff(i-1,j,k)+Dgo_eff(i,j,k))/2.)/dx)* &
                         (ne(i-1,j,k)+ne(i,j,k))/2./dx
                    cc(i,j,k) = dt*((alpha1_o(i+1,j,k)+alpha1_o(i,j,k))/2.*min(u(i+1,j,k),0.) - &
                         ((alpha1_o(i+1,j,k)+alpha1_o(i,j,k))/2.*(Dlo_eff(i+1,j,k)+Dlo_eff(i,j,k))/2. + &
                         (alpha2_o(i+1,j,k)+alpha2_o(i,j,k))/2.*(Dgo_eff(i+1,j,k)+Dgo_eff(i,j,k))/2.)/dx)* &
                         (ne(i+1,j,k)+ne(i,j,k))/2./dx
                    dd(i,j,k) = -dt*((alpha1_o(i,j-1,k)+alpha1_o(i,j,k))/2.*max(v(i,j,k),0.) + &
                         ((alpha1_o(i,j-1,k)+alpha1_o(i,j,k))/2.*(Dlo_eff(i,j-1,k)+Dlo_eff(i,j,k))/2. + &
                         (alpha2_o(i,j-1,k)+alpha2_o(i,j,k))/2.*(Dgo_eff(i,j-1,k)+Dgo_eff(i,j,k))/2.)/dy)* &
                         (ne(i,j-1,k)+ne(i,j,k))/2./dy
                    ee(i,j,k) = dt*((alpha1_o(i,j+1,k)+alpha1_o(i,j,k))/2.*min(v(i,j+1,k),0.) - &
                         ((alpha1_o(i,j+1,k)+alpha1_o(i,j,k))/2.*(Dlo_eff(i,j+1,k)+Dlo_eff(i,j,k))/2. + &
                         (alpha2_o(i,j+1,k)+alpha2_o(i,j,k))/2.*(Dgo_eff(i,j+1,k)+Dgo_eff(i,j,k))/2.)/dy)* &
                         (ne(i,j+1,k)+ne(i,j,k))/2./dy
                    ff(i,j,k) = -dt*((alpha1_o(i,j,k-1)+alpha1_o(i,j,k))/2.*max(w(i,j,k),0.) + &
                         ((alpha1_o(i,j,k-1)+alpha1_o(i,j,k))/2.*(Dlo_eff(i,j,k-1)+Dlo_eff(i,j,k))/2. + &
                         (alpha2_o(i,j,k-1)+alpha2_o(i,j,k))/2.*(Dgo_eff(i,j,k-1)+Dgo_eff(i,j,k))/2.)/dz)* &
                         (ne(i,j,k-1)+ne(i,j,k))/2./dz
                    gg(i,j,k) = dt*((alpha1_o(i,j,k+1)+alpha1_o(i,j,k))/2.*min(w(i,j,k+1),0.) - &
                         ((alpha1_o(i,j,k+1)+alpha1_o(i,j,k))/2.*(Dlo_eff(i,j,k+1)+Dlo_eff(i,j,k))/2. + &
                         (alpha2_o(i,j,k+1)+alpha2_o(i,j,k))/2.*(Dgo_eff(i,j,k+1)+Dgo_eff(i,j,k))/2.)/dz)* &
                         (ne(i,j,k+1)+ne(i,j,k))/2./dz

                    q(i,j,k) = ne(i,j,k)*to_o(i,j,k)
                    endif
                 Enddo
              Enddo
           Enddo



        else 

           Do i=2,Nx-1
              Do j=2,Ny-1
                    if(wall(i,j,1)) then
                       aa(i,j,1) = 0.
                       bb(i,j,1) = 0.
                       cc(i,j,1) = 0.
                       dd(i,j,1) = 0.
                       ee(i,j,1) = 0.
                       ff(i,j,1) = 0.
                       gg(i,j,1) = 0.
                       q(i,j,1) = 0.
                    else
                       aa(i,j,1) = ne(i,j,1) + dt*(-(alpha1_o(i-1,j,1)+alpha1_o(i,j,1))/2.*min(u(i,j,1),0.)* &
                            (ne(i-1,j,1)+ne(i,j,1))/2./dx + &
                            ((alpha1_o(i-1,j,1)+alpha1_o(i,j,1))/2.*(Dlo_eff(i-1,j,1)+Dlo_eff(i,j,1))/2. + &
                            (alpha2_o(i-1,j,1)+alpha2_o(i,j,1))/2.*(Dgo_eff(i-1,j,1)+Dgo_eff(i,j,1))/2.)/dx* &
                            (ne(i-1,j,1)+ne(i,j,1))/2./dx + &
                            (alpha1_o(i+1,j,1)+alpha1_o(i,j,1))/2.*max(u(i+1,j,1),0.)* &
                            (ne(i+1,j,1)+ne(i,j,1))/2./dx + &
                            ((alpha1_o(i+1,j,1)+alpha1_o(i,j,1))/2.*(Dlo_eff(i+1,j,1)+Dlo_eff(i,j,1))/2. + &
                            (alpha2_o(i+1,j,1)+alpha2_o(i,j,1))/2.*(Dgo_eff(i+1,j,1)+Dgo_eff(i,j,1))/2.)/dx* &
                            (ne(i+1,j,1)+ne(i,j,1))/2./dx - &
                            (alpha1_o(i,j-1,1)+alpha1_o(i,j,1))/2.*min(v(i,j,1),0.)*(ne(i,j-1,1)+ne(i,j,1))/2./dy + &
                            ((alpha1_o(i,j-1,1)+alpha1_o(i,j,1))/2.*(Dlo_eff(i,j-1,1)+Dlo_eff(i,j,1))/2. + &
                            (alpha2_o(i,j-1,1)+alpha2_o(i,j,1))/2.*(Dgo_eff(i,j-1,1)+Dgo_eff(i,j,1))/2.)/dy* &
                            (ne(i,j-1,1)+ne(i,j,1))/2./dy + &
                            (alpha1_o(i,j+1,1)+alpha1_o(i,j,1))/2.*max(v(i,j+1,1),0.)*(ne(i,j+1,1)+ne(i,j,1))/2./dy + &
                            ((alpha1_o(i,j+1,1)+alpha1_o(i,j,1))/2.*(Dlo_eff(i,j+1,1)+Dlo_eff(i,j,1))/2. + &
                            (alpha2_o(i,j+1,1)+alpha2_o(i,j,1))/2.*(Dgo_eff(i,j+1,1)+Dgo_eff(i,j,1))/2.)/dy* &
                            (ne(i,j+1,1)+ne(i,j,1))/2./dy - &
                            (tmp_alpha1_o(i,j,1)+alpha1_o(i,j,1))/2.*min(w(i,j,1),0.)*(tmp_ne(i,j,1)+ne(i,j,1))/2./dz + &
                            ((tmp_alpha1_o(i,j,1)+alpha1_o(i,j,1))/2.*(tmp_Dlo_eff(i,j,1)+Dlo_eff(i,j,1))/2. + &
                            (tmp_alpha2_o(i,j,1)+alpha2_o(i,j,1))/2.*(tmp_Dgo_eff(i,j,1)+Dgo_eff(i,j,1))/2.)/dz* &
                            (tmp_ne(i,j,1)+ne(i,j,1))/2./dz + &
                            (alpha1_o(i,j,1+1)+alpha1_o(i,j,1))/2.*max(w(i,j,1+1),0.)*(ne(i,j,1+1)+ne(i,j,1))/2./dz + &
                            ((alpha1_o(i,j,1+1)+alpha1_o(i,j,1))/2.*(Dlo_eff(i,j,1+1)+Dlo_eff(i,j,1))/2. + &
                            (alpha2_o(i,j,1+1)+alpha2_o(i,j,1))/2.*(Dgo_eff(i,j,1+1)+Dgo_eff(i,j,1))/2.)/dz* &
                            (ne(i,j,1+1)+ne(i,j,1))/2./dz) - &
                            dt*ne(i,j,1)*(-alpha1_o(i,j,1)*qo_0*b_o(i,j,1)* &
                            lc_o(i,j,1)/(Kc+lc_o(i,j,1))/(Ko+lo_o(i,j,1)))

                       bb(i,j,1) = -dt*((alpha1_o(i-1,j,1)+alpha1_o(i,j,1))/2.*max(u(i,j,1),0.) + &
                            ((alpha1_o(i-1,j,1)+alpha1_o(i,j,1))/2.*(Dlo_eff(i-1,j,1)+Dlo_eff(i,j,1))/2. + &
                            (alpha2_o(i-1,j,1)+alpha2_o(i,j,1))/2.*(Dgo_eff(i-1,j,1)+Dgo_eff(i,j,1))/2.)/dx)* &
                            (ne(i-1,j,1)+ne(i,j,1))/2./dx
                       cc(i,j,1) = dt*((alpha1_o(i+1,j,1)+alpha1_o(i,j,1))/2.*min(u(i+1,j,1),0.) - &
                            ((alpha1_o(i+1,j,1)+alpha1_o(i,j,1))/2.*(Dlo_eff(i+1,j,1)+Dlo_eff(i,j,1))/2. + &
                            (alpha2_o(i+1,j,1)+alpha2_o(i,j,1))/2.*(Dgo_eff(i+1,j,1)+Dgo_eff(i,j,1))/2.)/dx)* &
                            (ne(i+1,j,1)+ne(i,j,1))/2./dx
                       dd(i,j,1) = -dt*((alpha1_o(i,j-1,1)+alpha1_o(i,j,1))/2.*max(v(i,j,1),0.) + &
                            ((alpha1_o(i,j-1,1)+alpha1_o(i,j,1))/2.*(Dlo_eff(i,j-1,1)+Dlo_eff(i,j,1))/2. + &
                            (alpha2_o(i,j-1,1)+alpha2_o(i,j,1))/2.*(Dgo_eff(i,j-1,1)+Dgo_eff(i,j,1))/2.)/dy)* &
                            (ne(i,j-1,1)+ne(i,j,1))/2./dy
                       ee(i,j,1) = dt*((alpha1_o(i,j+1,1)+alpha1_o(i,j,1))/2.*min(v(i,j+1,1),0.) - &
                            ((alpha1_o(i,j+1,1)+alpha1_o(i,j,1))/2.*(Dlo_eff(i,j+1,1)+Dlo_eff(i,j,1))/2. + &
                            (alpha2_o(i,j+1,1)+alpha2_o(i,j,1))/2.*(Dgo_eff(i,j+1,1)+Dgo_eff(i,j,1))/2.)/dy)* &
                            (ne(i,j+1,1)+ne(i,j,1))/2./dy
                       ff(i,j,1) = 0.
                       gg(i,j,1) = dt*((alpha1_o(i,j,1+1)+alpha1_o(i,j,1))/2.*min(w(i,j,1+1),0.) - &
                            ((alpha1_o(i,j,1+1)+alpha1_o(i,j,1))/2.*(Dlo_eff(i,j,1+1)+Dlo_eff(i,j,1))/2. + &
                            (alpha2_o(i,j,1+1)+alpha2_o(i,j,1))/2.*(Dgo_eff(i,j,1+1)+Dgo_eff(i,j,1))/2.)/dz)* &
                            (ne(i,j,1+1)+ne(i,j,1))/2./dz

                       q(i,j,1) = ne(i,j,1)*to_o(i,j,1) - &
                            (-dt*((tmp_alpha1_o(i,j,1)+alpha1_o(i,j,1))/2.*max(w(i,j,1),0.) + &
                             ((tmp_alpha1_o(i,j,1)+alpha1_o(i,j,1))/2.*(tmp_Dlo_eff(i,j,1)+Dlo_eff(i,j,1))/2. + &
                            (tmp_alpha2_o(i,j,1)+alpha2_o(i,j,1))/2.*(tmp_Dgo_eff(i,j,1)+Dgo_eff(i,j,1))/2.)/dz)* &
                            (tmp_ne(i,j,1)+ne(i,j,1))/2./dz)*tmp_to_o(i,j,1)
                    endif
              Enddo
           Enddo


           Do i=2,Nx-1
              Do j=2,Ny-1
                    if(wall(i,j,Nz)) then
                       aa(i,j,Nz) = 0.
                       bb(i,j,Nz) = 0.
                       cc(i,j,Nz) = 0.
                       dd(i,j,Nz) = 0.
                       ee(i,j,Nz) = 0.
                       ff(i,j,Nz) = 0.
                       gg(i,j,Nz) = 0.
                       q(i,j,Nz) = 0.
                    else
                       aa(i,j,Nz) = ne(i,j,Nz) + dt*(-(alpha1_o(i-1,j,Nz)+alpha1_o(i,j,Nz))/2.*min(u(i,j,Nz),0.)* &
                            (ne(i-1,j,Nz)+ne(i,j,Nz))/2./dx + &
                            ((alpha1_o(i-1,j,Nz)+alpha1_o(i,j,Nz))/2.*(Dlo_eff(i-1,j,Nz)+Dlo_eff(i,j,Nz))/2. + &
                            (alpha2_o(i-1,j,Nz)+alpha2_o(i,j,Nz))/2.*(Dgo_eff(i-1,j,Nz)+Dgo_eff(i,j,Nz))/2.)/dx* &
                            (ne(i-1,j,Nz)+ne(i,j,Nz))/2./dx + &
                            (alpha1_o(i+1,j,Nz)+alpha1_o(i,j,Nz))/2.*max(u(i+1,j,Nz),0.)* &
                            (ne(i+1,j,Nz)+ne(i,j,Nz))/2./dx + &
                            ((alpha1_o(i+1,j,Nz)+alpha1_o(i,j,Nz))/2.*(Dlo_eff(i+1,j,Nz)+Dlo_eff(i,j,Nz))/2. + &
                            (alpha2_o(i+1,j,Nz)+alpha2_o(i,j,Nz))/2.*(Dgo_eff(i+1,j,Nz)+Dgo_eff(i,j,Nz))/2.)/dx* &
                            (ne(i+1,j,Nz)+ne(i,j,Nz))/2./dx - &
                            (alpha1_o(i,j-1,Nz)+alpha1_o(i,j,Nz))/2.*min(v(i,j,Nz),0.)*(ne(i,j-1,Nz)+ne(i,j,Nz))/2./dy + &
                            ((alpha1_o(i,j-1,Nz)+alpha1_o(i,j,Nz))/2.*(Dlo_eff(i,j-1,Nz)+Dlo_eff(i,j,Nz))/2. + &
                            (alpha2_o(i,j-1,Nz)+alpha2_o(i,j,Nz))/2.*(Dgo_eff(i,j-1,Nz)+Dgo_eff(i,j,Nz))/2.)/dy* &
                            (ne(i,j-1,Nz)+ne(i,j,Nz))/2./dy + &
                            (alpha1_o(i,j+1,Nz)+alpha1_o(i,j,Nz))/2.*max(v(i,j+1,Nz),0.)*(ne(i,j+1,Nz)+ne(i,j,Nz))/2./dy + &
                            ((alpha1_o(i,j+1,Nz)+alpha1_o(i,j,Nz))/2.*(Dlo_eff(i,j+1,Nz)+Dlo_eff(i,j,Nz))/2. + &
                            (alpha2_o(i,j+1,Nz)+alpha2_o(i,j,Nz))/2.*(Dgo_eff(i,j+1,Nz)+Dgo_eff(i,j,Nz))/2.)/dy* &
                            (ne(i,j+1,Nz)+ne(i,j,Nz))/2./dy - &
                            (alpha1_o(i,j,Nz-1)+alpha1_o(i,j,Nz))/2.*min(w(i,j,Nz),0.)*(ne(i,j,Nz-1)+ne(i,j,Nz))/2./dz + &
                            ((alpha1_o(i,j,Nz-1)+alpha1_o(i,j,Nz))/2.*(Dlo_eff(i,j,Nz-1)+Dlo_eff(i,j,Nz))/2. + &
                            (alpha2_o(i,j,Nz-1)+alpha2_o(i,j,Nz))/2.*(Dgo_eff(i,j,Nz-1)+Dgo_eff(i,j,Nz))/2.)/dz* &
                            (ne(i,j,Nz-1)+ne(i,j,Nz))/2./dz + &
                            (tmp_alpha1_o(i,j,2)+alpha1_o(i,j,Nz))/2.*max(w(i,j,Nz+1),0.)*(tmp_ne(i,j,2)+ne(i,j,Nz))/2./dz + &
                            ((tmp_alpha1_o(i,j,2)+alpha1_o(i,j,Nz))/2.*(tmp_Dlo_eff(i,j,2)+Dlo_eff(i,j,Nz))/2. + &
                            (tmp_alpha2_o(i,j,2)+alpha2_o(i,j,Nz))/2.*(tmp_Dgo_eff(i,j,2)+Dgo_eff(i,j,Nz))/2.)/dz* &
                            (tmp_ne(i,j,2)+ne(i,j,Nz))/2./dz) - &
                            dt*ne(i,j,Nz)* &
                            (-alpha1_o(i,j,Nz)*qo_0*b_o(i,j,Nz)* &
                            lc_o(i,j,Nz)/(Kc+lc_o(i,j,Nz))/(Ko+lo_o(i,j,Nz)))

                       bb(i,j,Nz) = -dt*((alpha1_o(i-1,j,Nz)+alpha1_o(i,j,Nz))/2.*max(u(i,j,Nz),0.) + &
                            ((alpha1_o(i-1,j,Nz)+alpha1_o(i,j,Nz))/2.*(Dlo_eff(i-1,j,Nz)+Dlo_eff(i,j,Nz))/2. + &
                            (alpha2_o(i-1,j,Nz)+alpha2_o(i,j,Nz))/2.*(Dgo_eff(i-1,j,Nz)+Dgo_eff(i,j,Nz))/2.)/dx)* &
                            (ne(i-1,j,Nz)+ne(i,j,Nz))/2./dx
                       cc(i,j,Nz) = dt*((alpha1_o(i+1,j,Nz)+alpha1_o(i,j,Nz))/2.*min(u(i+1,j,Nz),0.) - &
                            ((alpha1_o(i+1,j,Nz)+alpha1_o(i,j,Nz))/2.*(Dlo_eff(i+1,j,Nz)+Dlo_eff(i,j,Nz))/2. + &
                            (alpha2_o(i+1,j,Nz)+alpha2_o(i,j,Nz))/2.*(Dgo_eff(i+1,j,Nz)+Dgo_eff(i,j,Nz))/2.)/dx)* &
                            (ne(i+1,j,Nz)+ne(i,j,Nz))/2./dx
                       dd(i,j,Nz) = -dt*((alpha1_o(i,j-1,Nz)+alpha1_o(i,j,Nz))/2.*max(v(i,j,Nz),0.) + &
                            ((alpha1_o(i,j-1,Nz)+alpha1_o(i,j,Nz))/2.*(Dlo_eff(i,j-1,Nz)+Dlo_eff(i,j,Nz))/2. + &
                            (alpha2_o(i,j-1,Nz)+alpha2_o(i,j,Nz))/2.*(Dgo_eff(i,j-1,Nz)+Dgo_eff(i,j,Nz))/2.)/dy)* &
                            (ne(i,j-1,Nz)+ne(i,j,Nz))/2./dy
                       ee(i,j,Nz) = dt*((alpha1_o(i,j+1,Nz)+alpha1_o(i,j,Nz))/2.*min(v(i,j+1,Nz),0.) - &
                            ((alpha1_o(i,j+1,Nz)+alpha1_o(i,j,Nz))/2.*(Dlo_eff(i,j+1,Nz)+Dlo_eff(i,j,Nz))/2. + &
                            (alpha2_o(i,j+1,Nz)+alpha2_o(i,j,Nz))/2.*(Dgo_eff(i,j+1,Nz)+Dgo_eff(i,j,Nz))/2.)/dy)* &
                            (ne(i,j+1,Nz)+ne(i,j,Nz))/2./dy
                       ff(i,j,Nz) = -dt*((alpha1_o(i,j,Nz-1)+alpha1_o(i,j,Nz))/2.*max(w(i,j,Nz),0.) + &
                             ((alpha1_o(i,j,Nz-1)+alpha1_o(i,j,Nz))/2.*(Dlo_eff(i,j,Nz-1)+Dlo_eff(i,j,Nz))/2. + &
                            (alpha2_o(i,j,Nz-1)+alpha2_o(i,j,Nz))/2.*(Dgo_eff(i,j,Nz-1)+Dgo_eff(i,j,Nz))/2.)/dz)* &
                            (ne(i,j,Nz-1)+ne(i,j,Nz))/2./dz
                       gg(i,j,Nz) = 0

                       q(i,j,Nz) = ne(i,j,Nz)*to_o(i,j,Nz) - &
                            (dt*((tmp_alpha1_o(i,j,2)+alpha1_o(i,j,Nz))/2.*min(w(i,j,Nz+1),0.) - &
                            ((tmp_alpha1_o(i,j,2)+alpha1_o(i,j,Nz))/2.*(tmp_Dlo_eff(i,j,2)+Dlo_eff(i,j,Nz))/2. + &
                            (tmp_alpha2_o(i,j,2)+alpha2_o(i,j,Nz))/2.*(tmp_Dgo_eff(i,j,2)+Dgo_eff(i,j,Nz))/2.)/dz)* &
                            (tmp_ne(i,j,2)+ne(i,j,Nz))/2./dz)*tmp_to_o(i,j,2)
                    endif
              Enddo
           Enddo


           !! for inside

           Do i=2,Nx-1
              Do j=2,Ny-1
                 Do k=2,Nz-1
                    if(wall(i,j,k)) then
                       aa(i,j,k) = 0.
                       bb(i,j,k) = 0.
                       cc(i,j,k) = 0.
                       dd(i,j,k) = 0.
                       ee(i,j,k) = 0.
                       ff(i,j,k) = 0.
                       gg(i,j,k) = 0.
                       q(i,j,k) = 0.
                    else
                       aa(i,j,k) = ne(i,j,k) + dt*(-(alpha1_o(i-1,j,k)+alpha1_o(i,j,k))/2.*min(u(i,j,k),0.)* &
                            (ne(i-1,j,k)+ne(i,j,k))/2./dx + &
                            ((alpha1_o(i-1,j,k)+alpha1_o(i,j,k))/2.*(Dlo_eff(i-1,j,k)+Dlo_eff(i,j,k))/2. + &
                            (alpha2_o(i-1,j,k)+alpha2_o(i,j,k))/2.*(Dgo_eff(i-1,j,k)+Dgo_eff(i,j,k))/2.)/dx* &
                            (ne(i-1,j,k)+ne(i,j,k))/2./dx + &
                            (alpha1_o(i+1,j,k)+alpha1_o(i,j,k))/2.*max(u(i+1,j,k),0.)* &
                            (ne(i+1,j,k)+ne(i,j,k))/2./dx + &
                            ((alpha1_o(i+1,j,k)+alpha1_o(i,j,k))/2.*(Dlo_eff(i+1,j,k)+Dlo_eff(i,j,k))/2. + &
                            (alpha2_o(i+1,j,k)+alpha2_o(i,j,k))/2.*(Dgo_eff(i+1,j,k)+Dgo_eff(i,j,k))/2.)/dx* &
                            (ne(i+1,j,k)+ne(i,j,k))/2./dx - &
                            (alpha1_o(i,j-1,k)+alpha1_o(i,j,k))/2.*min(v(i,j,k),0.)*(ne(i,j-1,k)+ne(i,j,k))/2./dy + &
                            ((alpha1_o(i,j-1,k)+alpha1_o(i,j,k))/2.*(Dlo_eff(i,j-1,k)+Dlo_eff(i,j,k))/2. + &
                            (alpha2_o(i,j-1,k)+alpha2_o(i,j,k))/2.*(Dgo_eff(i,j-1,k)+Dgo_eff(i,j,k))/2.)/dy* &
                            (ne(i,j-1,k)+ne(i,j,k))/2./dy + &
                            (alpha1_o(i,j+1,k)+alpha1_o(i,j,k))/2.*max(v(i,j+1,k),0.)*(ne(i,j+1,k)+ne(i,j,k))/2./dy + &
                            ((alpha1_o(i,j+1,k)+alpha1_o(i,j,k))/2.*(Dlo_eff(i,j+1,k)+Dlo_eff(i,j,k))/2. + &
                            (alpha2_o(i,j+1,k)+alpha2_o(i,j,k))/2.*(Dgo_eff(i,j+1,k)+Dgo_eff(i,j,k))/2.)/dy* &
                            (ne(i,j+1,k)+ne(i,j,k))/2./dy - &
                            (alpha1_o(i,j,k-1)+alpha1_o(i,j,k))/2.*min(w(i,j,k),0.)*(ne(i,j,k-1)+ne(i,j,k))/2./dz + &
                            ((alpha1_o(i,j,k-1)+alpha1_o(i,j,k))/2.*(Dlo_eff(i,j,k-1)+Dlo_eff(i,j,k))/2. + &
                            (alpha2_o(i,j,k-1)+alpha2_o(i,j,k))/2.*(Dgo_eff(i,j,k-1)+Dgo_eff(i,j,k))/2.)/dz* &
                            (ne(i,j,k-1)+ne(i,j,k))/2./dz + &
                            (alpha1_o(i,j,k+1)+alpha1_o(i,j,k))/2.*max(w(i,j,k+1),0.)*(ne(i,j,k+1)+ne(i,j,k))/2./dz + &
                            ((alpha1_o(i,j,k+1)+alpha1_o(i,j,k))/2.*(Dlo_eff(i,j,k+1)+Dlo_eff(i,j,k))/2. + &
                            (alpha2_o(i,j,k+1)+alpha2_o(i,j,k))/2.*(Dgo_eff(i,j,k+1)+Dgo_eff(i,j,k))/2.)/dz* &
                            (ne(i,j,k+1)+ne(i,j,k))/2./dz) - &
                            dt*ne(i,j,k)*(-alpha1_o(i,j,k)*qo_0*b_o(i,j,k)* &
                            lc_o(i,j,k)/(Kc+lc_o(i,j,k))/(Ko+lo_o(i,j,k)))

                       bb(i,j,k) = -dt*((alpha1_o(i-1,j,k)+alpha1_o(i,j,k))/2.*max(u(i,j,k),0.) + &
                            ((alpha1_o(i-1,j,k)+alpha1_o(i,j,k))/2.*(Dlo_eff(i-1,j,k)+Dlo_eff(i,j,k))/2. + &
                            (alpha2_o(i-1,j,k)+alpha2_o(i,j,k))/2.*(Dgo_eff(i-1,j,k)+Dgo_eff(i,j,k))/2.)/dx)* &
                            (ne(i-1,j,k)+ne(i,j,k))/2./dx
                       cc(i,j,k) = dt*((alpha1_o(i+1,j,k)+alpha1_o(i,j,k))/2.*min(u(i+1,j,k),0.) - &
                            ((alpha1_o(i+1,j,k)+alpha1_o(i,j,k))/2.*(Dlo_eff(i+1,j,k)+Dlo_eff(i,j,k))/2. + &
                            (alpha2_o(i+1,j,k)+alpha2_o(i,j,k))/2.*(Dgo_eff(i+1,j,k)+Dgo_eff(i,j,k))/2.)/dx)* &
                            (ne(i+1,j,k)+ne(i,j,k))/2./dx
                       dd(i,j,k) = -dt*((alpha1_o(i,j-1,k)+alpha1_o(i,j,k))/2.*max(v(i,j,k),0.) + &
                            ((alpha1_o(i,j-1,k)+alpha1_o(i,j,k))/2.*(Dlo_eff(i,j-1,k)+Dlo_eff(i,j,k))/2. + &
                            (alpha2_o(i,j-1,k)+alpha2_o(i,j,k))/2.*(Dgo_eff(i,j-1,k)+Dgo_eff(i,j,k))/2.)/dy)* &
                            (ne(i,j-1,k)+ne(i,j,k))/2./dy
                       ee(i,j,k) = dt*((alpha1_o(i,j+1,k)+alpha1_o(i,j,k))/2.*min(v(i,j+1,k),0.) - &
                            ((alpha1_o(i,j+1,k)+alpha1_o(i,j,k))/2.*(Dlo_eff(i,j+1,k)+Dlo_eff(i,j,k))/2. + &
                            (alpha2_o(i,j+1,k)+alpha2_o(i,j,k))/2.*(Dgo_eff(i,j+1,k)+Dgo_eff(i,j,k))/2.)/dy)* &
                            (ne(i,j+1,k)+ne(i,j,k))/2./dy
                       ff(i,j,k) = -dt*((alpha1_o(i,j,k-1)+alpha1_o(i,j,k))/2.*max(w(i,j,k),0.) + &
                             ((alpha1_o(i,j,k-1)+alpha1_o(i,j,k))/2.*(Dlo_eff(i,j,k-1)+Dlo_eff(i,j,k))/2. + &
                            (alpha2_o(i,j,k-1)+alpha2_o(i,j,k))/2.*(Dgo_eff(i,j,k-1)+Dgo_eff(i,j,k))/2.)/dz)* &
                            (ne(i,j,k-1)+ne(i,j,k))/2./dz
                       gg(i,j,k) = dt*((alpha1_o(i,j,k+1)+alpha1_o(i,j,k))/2.*min(w(i,j,k+1),0.) - &
                            ((alpha1_o(i,j,k+1)+alpha1_o(i,j,k))/2.*(Dlo_eff(i,j,k+1)+Dlo_eff(i,j,k))/2. + &
                            (alpha2_o(i,j,k+1)+alpha2_o(i,j,k))/2.*(Dgo_eff(i,j,k+1)+Dgo_eff(i,j,k))/2.)/dz)* &
                            (ne(i,j,k+1)+ne(i,j,k))/2./dz

                       q(i,j,k) = ne(i,j,k)*to_o(i,j,k)
                    endif
                 Enddo
              Enddo
           Enddo


        endif
     endif







     !! cabulate seven-diagonal elements near wall


     Do i=2,Nx-1
        Do j=2,Ny-1
           Do k=1,Nz
              if(wall(i-1,j,k)) then
                 aa(i,j,k) = aa(i,j,k) - dt*((alpha1_o(i-1,j,k)+alpha1_o(i,j,k))/2.*(Dlo_eff(i-1,j,k)+Dlo_eff(i,j,k))/2. + &
                            (alpha2_o(i-1,j,k)+alpha2_o(i,j,k))/2.*(Dgo_eff(i-1,j,k)+Dgo_eff(i,j,k))/2.)/dx* &
                            (ne(i-1,j,k)+ne(i,j,k))/2./dx
                 bb(i,j,k) = 0.
              endif
              if(wall(i+1,j,k)) then
                 aa(i,j,k) = aa(i,j,k) - dt*((alpha1_o(i+1,j,k)+alpha1_o(i,j,k))/2.*(Dlo_eff(i+1,j,k)+Dlo_eff(i,j,k))/2. + &
                            (alpha2_o(i+1,j,k)+alpha2_o(i,j,k))/2.*(Dgo_eff(i+1,j,k)+Dgo_eff(i,j,k))/2.)/dx* &
                            (ne(i+1,j,k)+ne(i,j,k))/2./dx
                 cc(i,j,k) = 0.
              endif
              if(wall(i,j-1,k)) then
                 aa(i,j,k) = aa(i,j,k) - dt*((alpha1_o(i,j-1,k)+alpha1_o(i,j,k))/2.*(Dlo_eff(i,j-1,k)+Dlo_eff(i,j,k))/2. + &
                            (alpha2_o(i,j-1,k)+alpha2_o(i,j,k))/2.*(Dgo_eff(i,j-1,k)+Dgo_eff(i,j,k))/2.)/dy* &
                            (ne(i,j-1,k)+ne(i,j,k))/2./dy
                 dd(i,j,k) = 0.
              endif
              if(wall(i,j+1,k)) then
                 aa(i,j,k) = aa(i,j,k) - dt*((alpha1_o(i,j+1,k)+alpha1_o(i,j,k))/2.*(Dlo_eff(i,j+1,k)+Dlo_eff(i,j,k))/2. + &
                            (alpha2_o(i,j+1,k)+alpha2_o(i,j,k))/2.*(Dgo_eff(i,j+1,k)+Dgo_eff(i,j,k))/2.)/dy* &
                            (ne(i,j+1,k)+ne(i,j,k))/2./dy
                 ee(i,j,k) = 0.
              endif

              if(npes.eq.1) then
                 if(k.ne.1 .and. wall(i,j,k-1)) then
                    aa(i,j,k) = aa(i,j,k) - dt*((alpha1_o(i,j,k-1)+alpha1_o(i,j,k))/2.*(Dlo_eff(i,j,k-1)+Dlo_eff(i,j,k))/2. + &
                            (alpha2_o(i,j,k-1)+alpha2_o(i,j,k))/2.*(Dgo_eff(i,j,k-1)+Dgo_eff(i,j,k))/2.)/dz* &
                            (ne(i,j,k-1)+ne(i,j,k))/2./dz
                    ff(i,j,k) = 0.
                 endif
                 if(k.ne.Nz .and. wall(i,j,k+1)) then
                    aa(i,j,k) = aa(i,j,k) - dt*((alpha1_o(i,j,k+1)+alpha1_o(i,j,k))/2.*(Dlo_eff(i,j,k+1)+Dlo_eff(i,j,k))/2. + &
                            (alpha2_o(i,j,k+1)+alpha2_o(i,j,k))/2.*(Dgo_eff(i,j,k+1)+Dgo_eff(i,j,k))/2.)/dz* &
                            (ne(i,j,k+1)+ne(i,j,k))/2./dz
                    gg(i,j,k) = 0.
                 endif
              else
                 if(mype.eq.0) then
                    if(k.ne.1 .and. wall(i,j,k-1)) then
                       aa(i,j,k) = aa(i,j,k) - &
                            dt*((alpha1_o(i,j,k-1)+alpha1_o(i,j,k))/2.*(Dlo_eff(i,j,k-1)+Dlo_eff(i,j,k))/2. + &
                            (alpha2_o(i,j,k-1)+alpha2_o(i,j,k))/2.*(Dgo_eff(i,j,k-1)+Dgo_eff(i,j,k))/2.)/dz* &
                            (ne(i,j,k-1)+ne(i,j,k))/2./dz
                       ff(i,j,k) = 0.
                    endif
                    if(k.eq.Nz) then
                       if(tmp_wall(i,j,2)) then
                          aa(i,j,Nz) = aa(i,j,Nz) - &
                               dt*((tmp_alpha1_o(i,j,2)+alpha1_o(i,j,Nz))/2.*(tmp_Dlo_eff(i,j,2)+Dlo_eff(i,j,Nz))/2. + &
                               (tmp_alpha2_o(i,j,2)+alpha2_o(i,j,Nz))/2.*(tmp_Dgo_eff(i,j,2)+Dgo_eff(i,j,Nz))/2.)/dz* &
                               (tmp_ne(i,j,2)+ne(i,j,Nz))/2./dz
                          gg(i,j,Nz) = 0.
                       endif
                    else
                       if(wall(i,j,k+1)) then
                          aa(i,j,k) = aa(i,j,k) - &
                               dt*((alpha1_o(i,j,k+1)+alpha1_o(i,j,k))/2.*(Dlo_eff(i,j,k+1)+Dlo_eff(i,j,k))/2. + &
                               (alpha2_o(i,j,k+1)+alpha2_o(i,j,k))/2.*(Dgo_eff(i,j,k+1)+Dgo_eff(i,j,k))/2.)/dz* &
                               (ne(i,j,k+1)+ne(i,j,k))/2./dz
                          gg(i,j,k) = 0.
                       endif
                    endif
                 else if(mype.eq.npes-1) then
                    if(k.eq.1) then
                       if(tmp_wall(i,j,1)) then
                          aa(i,j,1) = aa(i,j,1) - &
                               dt*((tmp_alpha1_o(i,j,1)+alpha1_o(i,j,1))/2.*(tmp_Dlo_eff(i,j,1)+Dlo_eff(i,j,1))/2. + &
                               (tmp_alpha2_o(i,j,1)+alpha2_o(i,j,1))/2.*(tmp_Dgo_eff(i,j,1)+Dgo_eff(i,j,1))/2.)/dz* &
                               (tmp_ne(i,j,1)+ne(i,j,1))/2./dz
                          ff(i,j,1) = 0.
                       endif
                    else
                       if(wall(i,j,k-1)) then
                          aa(i,j,k) = aa(i,j,k) - &
                               dt*((alpha1_o(i,j,k-1)+alpha1_o(i,j,k))/2.*(Dlo_eff(i,j,k-1)+Dlo_eff(i,j,k))/2. + &
                               (alpha2_o(i,j,k-1)+alpha2_o(i,j,k))/2.*(Dgo_eff(i,j,k-1)+Dgo_eff(i,j,k))/2.)/dz* &
                               (ne(i,j,k-1)+ne(i,j,k))/2./dz
                          ff(i,j,k) = 0.
                       endif
                    endif
                    if(k.ne.Nz .and. wall(i,j,k+1)) then
                       aa(i,j,k) = aa(i,j,k) - &
                            dt*((alpha1_o(i,j,k+1)+alpha1_o(i,j,k))/2.*(Dlo_eff(i,j,k+1)+Dlo_eff(i,j,k))/2. + &
                            (alpha2_o(i,j,k+1)+alpha2_o(i,j,k))/2.*(Dgo_eff(i,j,k+1)+Dgo_eff(i,j,k))/2.)/dz* &
                            (ne(i,j,k+1)+ne(i,j,k))/2./dz
                       gg(i,j,k) = 0.
                    endif
                 else 
                    if(k.eq.1) then
                       if(tmp_wall(i,j,1)) then
                          aa(i,j,1) = aa(i,j,1) - &
                               dt*((tmp_alpha1_o(i,j,1)+alpha1_o(i,j,1))/2.*(tmp_Dlo_eff(i,j,1)+Dlo_eff(i,j,1))/2. + &
                               (tmp_alpha2_o(i,j,1)+alpha2_o(i,j,1))/2.*(tmp_Dgo_eff(i,j,1)+Dgo_eff(i,j,1))/2.)/dz* &
                               (tmp_ne(i,j,1)+ne(i,j,1))/2./dz
                          ff(i,j,1) = 0.
                       endif
                    else
                       if(wall(i,j,k-1)) then
                          aa(i,j,k) = aa(i,j,k) - &
                               dt*((alpha1_o(i,j,k-1)+alpha1_o(i,j,k))/2.*(Dlo_eff(i,j,k-1)+Dlo_eff(i,j,k))/2. + &
                               (alpha2_o(i,j,k-1)+alpha2_o(i,j,k))/2.*(Dgo_eff(i,j,k-1)+Dgo_eff(i,j,k))/2.)/dz* &
                               (ne(i,j,k-1)+ne(i,j,k))/2./dz
                          ff(i,j,k) = 0.
                       endif
                    endif
                    if(k.eq.Nz) then
                       if(tmp_wall(i,j,2)) then
                          aa(i,j,Nz) = aa(i,j,Nz) - &
                               dt*((tmp_alpha1_o(i,j,2)+alpha1_o(i,j,Nz))/2.*(tmp_Dlo_eff(i,j,2)+Dlo_eff(i,j,Nz))/2. + &
                               (tmp_alpha2_o(i,j,2)+alpha2_o(i,j,Nz))/2.*(tmp_Dgo_eff(i,j,2)+Dgo_eff(i,j,Nz))/2.)/dz* &
                               (tmp_ne(i,j,2)+ne(i,j,Nz))/2./dz
                          gg(i,j,Nz) = 0.
                       endif
                    else
                       if(wall(i,j,k+1)) then
                          aa(i,j,k) = aa(i,j,k) - &
                               dt*((alpha1_o(i,j,k+1)+alpha1_o(i,j,k))/2.*(Dlo_eff(i,j,k+1)+Dlo_eff(i,j,k))/2. + &
                               (alpha2_o(i,j,k+1)+alpha2_o(i,j,k))/2.*(Dgo_eff(i,j,k+1)+Dgo_eff(i,j,k))/2.)/dz* &
                               (ne(i,j,k+1)+ne(i,j,k))/2./dz
                          gg(i,j,k) = 0.
                       endif
                    endif
                 endif
              endif

           Enddo
        Enddo
     Enddo





     !! caloulate elements of matrix MM 


     Do k=1,Nz
        Do j=1,Ny
           Do i=1,Nx
              ma((k-1)*Nx*Ny+(j-1)*Nx+i) = aa(i,j,k)
              mb((k-1)*Nx*Ny+(j-1)*Nx+i) = bb(i,j,k)
              mc((k-1)*Nx*Ny+(j-1)*Nx+i) = cc(i,j,k)
              md((k-1)*Nx*Ny+(j-1)*Nx+i) = dd(i,j,k)
              me((k-1)*Nx*Ny+(j-1)*Nx+i) = ee(i,j,k)
              mf((k-1)*Nx*Ny+(j-1)*Nx+i) = ff(i,j,k)
              mg((k-1)*Nx*Ny+(j-1)*Nx+i) = gg(i,j,k)

              QQ((k-1)*Nx*Ny+(j-1)*Nx+i) = q(i,j,k)
              HOHO((k-1)*Nx*Ny+(j-1)*Nx+i) = to_o(i,j,k)

              WW((k-1)*Nx*Ny+(j-1)*Nx+i) = wall(i,j,k)
           Enddo
        Enddo
     Enddo








!!! solve to by successive overrelaxation (SOR)



        nn = 1

        Do while(nn.lt.1e5)

           if(WW(1)) then
              HH(1) = 0.
           else
              HH(1) = omega*(QQ(1)-mc(1)*HOHO(2)-me(1)*HOHO(1+Nx)-mg(1)*HOHO(1+Nx*Ny))/ &
                   (ma(1)+SMALL) + (1-omega)*HOHO(1)
           endif
           Do pp=2,Nx
              if(WW(pp)) then
                 HH(pp) = 0.
              else
                 HH(pp) = omega*(QQ(pp)-mb(pp)*HH(pp-1)-mc(pp)*HOHO(pp+1)-me(pp)*HOHO(pp+Nx)- &
                      mg(pp)*HOHO(pp+Nx*Ny))/(ma(pp)+SMALL) + (1-omega)*HOHO(pp)
              endif
           Enddo
           Do pp=1+Nx,Nx*Ny
              if(WW(pp)) then
                 HH(pp) = 0.
              else
                 HH(pp) = omega*(QQ(pp)-md(pp)*HH(pp-Nx)-mb(pp)*HH(pp-1)-mc(pp)*HOHO(pp+1)- &
                      me(pp)*HOHO(pp+Nx)-mg(pp)*HOHO(pp+Nx*Ny))/(ma(pp)+SMALL) + &
                      (1-omega)*HOHO(pp)
              endif
           Enddo
           Do pp=1+Nx*Ny,Nt-Nx*Ny
              if(WW(pp)) then
                 HH(pp) = 0.
              else
                 HH(pp) = omega*(QQ(pp)-mf(pp)*HH(pp-Nx*Ny)-md(pp)*HH(pp-Nx)- &
                      mb(pp)*HH(pp-1)-mc(pp)*HOHO(pp+1)-me(pp)*HOHO(pp+Nx)- &
                      mg(pp)*HOHO(pp+Nx*Ny))/(ma(pp)+SMALL) + (1-omega)*HOHO(pp)
              endif
           Enddo
           Do pp=Nt-Nx*Ny+1,Nt-Nx
              if(WW(pp)) then
                 HH(pp) = 0.
              else
                 HH(pp) = omega*(QQ(pp)-mf(pp)*HH(pp-Nx*Ny)-md(pp)*HH(pp-Nx)- &
                      mb(pp)*HH(pp-1)-mc(pp)*HOHO(pp+1)-me(pp)*HOHO(pp+Nx))/(ma(pp)+SMALL) + &
                      (1-omega)*HOHO(pp)
              endif
           Enddo
           Do pp=Nt-Nx+1,Nt-1
              if(WW(pp)) then
                 HH(pp) = 0.
              else
                 HH(pp) = omega*(QQ(pp)-mf(pp)*HH(pp-Nx*Ny)-md(pp)*HH(pp-Nx)- &
                      mb(pp)*HH(pp-1)-mc(pp)*HOHO(pp+1))/(ma(pp)+SMALL) + (1-omega)*HOHO(pp)
              endif
           Enddo
           if(WW(Nt)) then
              HH(Nt) = 0.
           else
              HH(Nt) = omega*(QQ(Nt)-mf(Nt)*HH(Nt-Nx*Ny)-md(Nt)*HH(Nt-Nx)- &
                   mb(Nt)*HH(Nt-1))/(ma(Nt)+SMALL) +(1-omega)*HOHO(Nt)
           endif


!!! calculate error
           HH_err_local = Sum(Abs(HH-HOHO))
           CALL MPI_ALLREDUCE(HH_err_local,HH_err,1,MPI_double_precision,MPI_SUM,MPI_COMM_WORLD,error)
           HH_err = HH_err/LX/LY/LZ


           If (HH_err.lt.HH_eps) exit


!!! update interface q

!!$        Do k=1,Nz
!!$           Do j=1,Ny
!!$              Do i=1,Nx
!!$                 to(i,j,k) = HH((k-1)*Nx*Ny+(j-1)*Nx+i)
!!$              Enddo
!!$           Enddo
!!$        Enddo

           Do j=1,Ny
              Do i=1,Nx
                 to(i,j,1) = HH((j-1)*Nx+i)
              Enddo
           Enddo

           Do j=1,Ny
              Do i=1,Nx
                 to(i,j,Nz) = HH((Nz-1)*Nx*Ny+(j-1)*Nx+i)
              Enddo
           Enddo




!!! exchange h during iteration

           call mpi_sendrecv(to(:,:,Nz),ilen,MPI_double_precision,inext,610, &
                tmp_to(:,:,1),ilen,MPI_double_precision,iprev,610,nallgrp,stat,error)

           call mpi_sendrecv(to(:,:,1),ilen,MPI_double_precision,iprev,620, &
                tmp_to(:,:,2),ilen,MPI_double_precision,inext,620,nallgrp,stat,error)


           CALL MPI_BARRIER(MPI_COMM_WORLD,error)


           
           if(mype.ne.0) then
              Do j=1,Ny
                 Do i=1,Nx
                    q(i,j,1) = ne(i,j,1)*to_o(i,j,1) - &
                            (-dt*((tmp_alpha1_o(i,j,1)+alpha1_o(i,j,1))/2.*max(w(i,j,1),0.) + &
                             ((tmp_alpha1_o(i,j,1)+alpha1_o(i,j,1))/2.*(tmp_Dlo_eff(i,j,1)+Dlo_eff(i,j,1))/2. + &
                            (tmp_alpha2_o(i,j,1)+alpha2_o(i,j,1))/2.*(tmp_Dgo_eff(i,j,1)+Dgo_eff(i,j,1))/2.)/dz)* &
                            (tmp_ne(i,j,1)+ne(i,j,1))/2./dz)*tmp_to(i,j,1)
                    QQ((j-1)*Nx+i) = q(i,j,1)
                 Enddo
              Enddo
           endif

           if(mype.ne.npes-1) then
              Do j=1,Ny
                 Do i=1,Nx
                    q(i,j,Nz) = ne(i,j,Nz)*to_o(i,j,Nz) - &
                            (dt*((tmp_alpha1_o(i,j,2)+alpha1_o(i,j,Nz))/2.*min(w(i,j,Nz+1),0.) - &
                            ((tmp_alpha1_o(i,j,2)+alpha1_o(i,j,Nz))/2.*(tmp_Dlo_eff(i,j,2)+Dlo_eff(i,j,Nz))/2. + &
                            (tmp_alpha2_o(i,j,2)+alpha2_o(i,j,Nz))/2.*(tmp_Dgo_eff(i,j,2)+Dgo_eff(i,j,Nz))/2.)/dz)* &
                            (tmp_ne(i,j,2)+ne(i,j,Nz))/2./dz)*tmp_to(i,j,2)
                    QQ((Nz-1)*Nx*Ny+(j-1)*Nx+i) = q(i,j,Nz)
                 Enddo
              Enddo
           endif



           HOHO = HH
           nn = nn+1

        ENDDO




        !! update h

        Do k=1,Nz
           Do j=1,Ny
              Do i=1,Nx
                 to(i,j,k) = HH((k-1)*Nx*Ny+(j-1)*Nx+i)
              Enddo
           Enddo
        Enddo





!!!! calculate lo and go 

        Do k=1,Nz
           Do j=1,Ny
              Do i=1,Nx
                 lo(i,j,k) = beta1_o(i,j,k)*to(i,j,k)
                 go(i,j,k) = beta2_o(i,j,k)*to(i,j,k)
              Enddo
           Enddo
        Enddo
        
        
        






!!!! solve td (total inorganic carbon) semi-implicitly

!!! boundary condition setting and
!!! communicate info between bottom and top planes


     call mpi_sendrecv(td_o(:,:,Nz),ilen,MPI_double_precision,inext,110, &
          tmp_td_o(:,:,1),ilen,MPI_double_precision,iprev,110,nallgrp,stat,error)
     call mpi_sendrecv(td_o(:,:,1),ilen,MPI_double_precision,iprev,120, &
          tmp_td_o(:,:,2),ilen,MPI_double_precision,inext,120,nallgrp,stat,error)


     CALL MPI_BARRIER(MPI_COMM_WORLD,error)

     

     
     if(npes.eq.1) then

        Do i=2,Nx-1
           Do j=2,Ny-1
              if(wall(i,j,1)) then
                 aa(i,j,1) = 0.
                 bb(i,j,1) = 0.
                 cc(i,j,1) = 0.
                 dd(i,j,1) = 0.
                 ee(i,j,1) = 0.
                 ff(i,j,1) = 0.
                 gg(i,j,1) = 0.
                 q(i,j,1) = 0.
              else
                 bb(i,j,1) = -dt*((alpha1_d(i-1,j,1)+alpha1_d(i,j,1))/2.*max(u(i,j,1),0.) + &
                      ((alpha1_d(i-1,j,1)+alpha1_d(i,j,1))/2.*(Dld_eff(i-1,j,1)+Dld_eff(i,j,1))/2. + &
                      (alpha2_d(i-1,j,1)+alpha2_d(i,j,1))/2.*(Dgd_eff(i-1,j,1)+Dgd_eff(i,j,1))/2.)/dx)* &
                      (ne(i-1,j,1)+ne(i,j,1))/2./dx
                 cc(i,j,1) = dt*((alpha1_d(i+1,j,1)+alpha1_d(i,j,1))/2.*min(u(i+1,j,1),0.) - &
                      ((alpha1_d(i+1,j,1)+alpha1_d(i,j,1))/2.*(Dld_eff(i+1,j,1)+Dld_eff(i,j,1))/2. + &
                      (alpha2_d(i+1,j,1)+alpha2_d(i,j,1))/2.*(Dgd_eff(i+1,j,1)+Dgd_eff(i,j,1))/2.)/dx)* &
                      (ne(i+1,j,1)+ne(i,j,1))/2./dx
                 dd(i,j,1) = -dt*((alpha1_d(i,j-1,1)+alpha1_d(i,j,1))/2.*max(v(i,j,1),0.) + &
                      ((alpha1_d(i,j-1,1)+alpha1_d(i,j,1))/2.*(Dld_eff(i,j-1,1)+Dld_eff(i,j,1))/2. + &
                      (alpha2_d(i,j-1,1)+alpha2_d(i,j,1))/2.*(Dgd_eff(i,j-1,1)+Dgd_eff(i,j,1))/2.)/dy)* &
                      (ne(i,j-1,1)+ne(i,j,1))/2./dy
                 ee(i,j,1) = dt*((alpha1_d(i,j+1,1)+alpha1_d(i,j,1))/2.*min(v(i,j+1,1),0.) - &
                      ((alpha1_d(i,j+1,1)+alpha1_d(i,j,1))/2.*(Dld_eff(i,j+1,1)+Dld_eff(i,j,1))/2. + &
                      (alpha2_d(i,j+1,1)+alpha2_d(i,j,1))/2.*(Dgd_eff(i,j+1,1)+Dgd_eff(i,j,1))/2.)/dy)* &
                      (ne(i,j+1,1)+ne(i,j,1))/2./dy
                 ff(i,j,1) = 0.
                 gg(i,j,1) = dt*((alpha1_d(i,j,1+1)+alpha1_d(i,j,1))/2.*min(w(i,j,1+1),0.) - &
                      ((alpha1_d(i,j,1+1)+alpha1_d(i,j,1))/2.*(Dld_eff(i,j,1+1)+Dld_eff(i,j,1))/2. + &
                      (alpha2_d(i,j,1+1)+alpha2_d(i,j,1))/2.*(Dgd_eff(i,j,1+1)+Dgd_eff(i,j,1))/2.)/dz)* &
                      (ne(i,j,1+1)+ne(i,j,1))/2./dz


                 if(tdbc_d.eq.1) then 
                    aa(i,j,1) = ne(i,j,1) + dt*(-(alpha1_d(i-1,j,1)+alpha1_d(i,j,1))/2.*min(u(i,j,1),0.)* &
                         (ne(i-1,j,1)+ne(i,j,1))/2./dx + &
                         ((alpha1_d(i-1,j,1)+alpha1_d(i,j,1))/2.*(Dld_eff(i-1,j,1)+Dld_eff(i,j,1))/2. + &
                         (alpha2_d(i-1,j,1)+alpha2_d(i,j,1))/2.*(Dgd_eff(i-1,j,1)+Dgd_eff(i,j,1))/2.)/dx* &
                         (ne(i-1,j,1)+ne(i,j,1))/2./dx + &
                         (alpha1_d(i+1,j,1)+alpha1_d(i,j,1))/2.*max(u(i+1,j,1),0.)* &
                         (ne(i+1,j,1)+ne(i,j,1))/2./dx + &
                         ((alpha1_d(i+1,j,1)+alpha1_d(i,j,1))/2.*(Dld_eff(i+1,j,1)+Dld_eff(i,j,1))/2. + &
                         (alpha2_d(i+1,j,1)+alpha2_d(i,j,1))/2.*(Dgd_eff(i+1,j,1)+Dgd_eff(i,j,1))/2.)/dx* &
                         (ne(i+1,j,1)+ne(i,j,1))/2./dx - &
                         (alpha1_d(i,j-1,1)+alpha1_d(i,j,1))/2.*min(v(i,j,1),0.)*(ne(i,j-1,1)+ne(i,j,1))/2./dy + &
                         ((alpha1_d(i,j-1,1)+alpha1_d(i,j,1))/2.*(Dld_eff(i,j-1,1)+Dld_eff(i,j,1))/2. + &
                         (alpha2_d(i,j-1,1)+alpha2_d(i,j,1))/2.*(Dgd_eff(i,j-1,1)+Dgd_eff(i,j,1))/2.)/dy* &
                         (ne(i,j-1,1)+ne(i,j,1))/2./dy + &
                         (alpha1_d(i,j+1,1)+alpha1_d(i,j,1))/2.*max(v(i,j+1,1),0.)*(ne(i,j+1,1)+ne(i,j,1))/2./dy + &
                         ((alpha1_d(i,j+1,1)+alpha1_d(i,j,1))/2.*(Dld_eff(i,j+1,1)+Dld_eff(i,j,1))/2. + &
                         (alpha2_d(i,j+1,1)+alpha2_d(i,j,1))/2.*(Dgd_eff(i,j+1,1)+Dgd_eff(i,j,1))/2.)/dy* &
                         (ne(i,j+1,1)+ne(i,j,1))/2./dy + &
                         (alpha1_d(i,j,1)*Dld_eff(i,j,1)+alpha2_d(i,j,1)*Dgd_eff(i,j,1))/dz*2.*ne(i,j,1)/dz + &
                         (alpha1_d(i,j,1+1)+alpha1_d(i,j,1))/2.*max(w(i,j,1+1),0.)*(ne(i,j,1+1)+ne(i,j,1))/2./dz + &
                         ((alpha1_d(i,j,1+1)+alpha1_d(i,j,1))/2.*(Dld_eff(i,j,1+1)+Dld_eff(i,j,1))/2. + &
                         (alpha2_d(i,j,1+1)+alpha2_d(i,j,1))/2.*(Dgd_eff(i,j,1+1)+Dgd_eff(i,j,1))/2.)/dz* &
                         (ne(i,j,1+1)+ne(i,j,1))/2./dz)

                    q(i,j,1) = ne(i,j,1)*td_o(i,j,1)  + &
                            dt*theta(i,j,1)* &
                            (qd_0*b_o(i,j,1)*lc_o(i,j,1)/(Kc+lc_o(i,j,1))*lo_o(i,j,1)/(Ko+lo_o(i,j,1))) - &
                            (-dt*(alpha1_d(i,j,1)*max(w(i,j,1),0.) + &
                            2.*(alpha1_d(i,j,1)*Dld_eff(i,j,1)+alpha2_d(i,j,1)*Dgd_eff(i,j,1))/dz)* &
                            ne(i,j,1)/dz)*td_d(i,j)
                 endif
                 if(tdbc_d.eq.2) then 
                    aa(i,j,1) = ne(i,j,1) + dt*(-(alpha1_d(i-1,j,1)+alpha1_d(i,j,1))/2.*min(u(i,j,1),0.)* &
                         (ne(i-1,j,1)+ne(i,j,1))/2./dx + &
                         ((alpha1_d(i-1,j,1)+alpha1_d(i,j,1))/2.*(Dld_eff(i-1,j,1)+Dld_eff(i,j,1))/2. + &
                         (alpha2_d(i-1,j,1)+alpha2_d(i,j,1))/2.*(Dgd_eff(i-1,j,1)+Dgd_eff(i,j,1))/2.)/dx* &
                         (ne(i-1,j,1)+ne(i,j,1))/2./dx + &
                         (alpha1_d(i+1,j,1)+alpha1_d(i,j,1))/2.*max(u(i+1,j,1),0.)* &
                         (ne(i+1,j,1)+ne(i,j,1))/2./dx + &
                         ((alpha1_d(i+1,j,1)+alpha1_d(i,j,1))/2.*(Dld_eff(i+1,j,1)+Dld_eff(i,j,1))/2. + &
                         (alpha2_d(i+1,j,1)+alpha2_d(i,j,1))/2.*(Dgd_eff(i+1,j,1)+Dgd_eff(i,j,1))/2.)/dx* &
                         (ne(i+1,j,1)+ne(i,j,1))/2./dx - &
                         (alpha1_d(i,j-1,1)+alpha1_d(i,j,1))/2.*min(v(i,j,1),0.)*(ne(i,j-1,1)+ne(i,j,1))/2./dy + &
                         ((alpha1_d(i,j-1,1)+alpha1_d(i,j,1))/2.*(Dld_eff(i,j-1,1)+Dld_eff(i,j,1))/2. + &
                         (alpha2_d(i,j-1,1)+alpha2_d(i,j,1))/2.*(Dgd_eff(i,j-1,1)+Dgd_eff(i,j,1))/2.)/dy* &
                         (ne(i,j-1,1)+ne(i,j,1))/2./dy + &
                         (alpha1_d(i,j+1,1)+alpha1_d(i,j,1))/2.*max(v(i,j+1,1),0.)*(ne(i,j+1,1)+ne(i,j,1))/2./dy + &
                         ((alpha1_d(i,j+1,1)+alpha1_d(i,j,1))/2.*(Dld_eff(i,j+1,1)+Dld_eff(i,j,1))/2. + &
                         (alpha2_d(i,j+1,1)+alpha2_d(i,j,1))/2.*(Dgd_eff(i,j+1,1)+Dgd_eff(i,j,1))/2.)/dy* &
                         (ne(i,j+1,1)+ne(i,j,1))/2./dy - &
                         alpha1_d(i,j,1)*min(w(i,j,1),0.)*ne(i,j,1)/dz + &
                         (alpha1_d(i,j,1+1)+alpha1_d(i,j,1))/2.*max(w(i,j,1+1),0.)*(ne(i,j,1+1)+ne(i,j,1))/2./dz + &
                         ((alpha1_d(i,j,1+1)+alpha1_d(i,j,1))/2.*(Dld_eff(i,j,1+1)+Dld_eff(i,j,1))/2. + &
                         (alpha2_d(i,j,1+1)+alpha2_d(i,j,1))/2.*(Dgd_eff(i,j,1+1)+Dgd_eff(i,j,1))/2.)/dz* &
                         (ne(i,j,1+1)+ne(i,j,1))/2./dz)
                    q(i,j,1) = ne(i,j,1)*td_o(i,j,1)  + &
                         dt*theta(i,j,1)* &
                         (qd_0*b_o(i,j,1)*lc_o(i,j,1)/(Kc+lc_o(i,j,1))*lo_o(i,j,1)/(Ko+lo_o(i,j,1)))

                 endif
              endif
           enddo
        enddo


        Do i=2,Nx-1
           Do j=2,Ny-1
              if(wall(i,j,Nz)) then
                 aa(i,j,Nz) = 0.
                 bb(i,j,Nz) = 0.
                 cc(i,j,Nz) = 0.
                 dd(i,j,Nz) = 0.
                 ee(i,j,Nz) = 0.
                 ff(i,j,Nz) = 0.
                 gg(i,j,Nz) = 0.
                 q(i,j,Nz) = 0.
              else     
                 bb(i,j,Nz) = -dt*((alpha1_d(i-1,j,Nz)+alpha1_d(i,j,Nz))/2.*max(u(i,j,Nz),0.) + &
                      ((alpha1_d(i-1,j,Nz)+alpha1_d(i,j,Nz))/2.*(Dld_eff(i-1,j,Nz)+Dld_eff(i,j,Nz))/2. + &
                      (alpha2_d(i-1,j,Nz)+alpha2_d(i,j,Nz))/2.*(Dgd_eff(i-1,j,Nz)+Dgd_eff(i,j,Nz))/2.)/dx)* &
                      (ne(i-1,j,Nz)+ne(i,j,Nz))/2./dx
                 cc(i,j,Nz) = dt*((alpha1_d(i+1,j,Nz)+alpha1_d(i,j,Nz))/2.*min(u(i+1,j,Nz),0.) - &
                      ((alpha1_d(i+1,j,Nz)+alpha1_d(i,j,Nz))/2.*(Dld_eff(i+1,j,Nz)+Dld_eff(i,j,Nz))/2. + &
                      (alpha2_d(i+1,j,Nz)+alpha2_d(i,j,Nz))/2.*(Dgd_eff(i+1,j,Nz)+Dgd_eff(i,j,Nz))/2.)/dx)* &
                      (ne(i+1,j,Nz)+ne(i,j,Nz))/2./dx
                 dd(i,j,Nz) = -dt*((alpha1_d(i,j-1,Nz)+alpha1_d(i,j,Nz))/2.*max(v(i,j,Nz),0.) + &
                      ((alpha1_d(i,j-1,Nz)+alpha1_d(i,j,Nz))/2.*(Dld_eff(i,j-1,Nz)+Dld_eff(i,j,Nz))/2. + &
                      (alpha2_d(i,j-1,Nz)+alpha2_d(i,j,Nz))/2.*(Dgd_eff(i,j-1,Nz)+Dgd_eff(i,j,Nz))/2.)/dy)* &
                      (ne(i,j-1,Nz)+ne(i,j,Nz))/2./dy
                 ee(i,j,Nz) = dt*((alpha1_d(i,j+1,Nz)+alpha1_d(i,j,Nz))/2.*min(v(i,j+1,Nz),0.) - &
                      ((alpha1_d(i,j+1,Nz)+alpha1_d(i,j,Nz))/2.*(Dld_eff(i,j+1,Nz)+Dld_eff(i,j,Nz))/2. + &
                      (alpha2_d(i,j+1,Nz)+alpha2_d(i,j,Nz))/2.*(Dgd_eff(i,j+1,Nz)+Dgd_eff(i,j,Nz))/2.)/dy)* &
                      (ne(i,j+1,Nz)+ne(i,j,Nz))/2./dy
                 ff(i,j,Nz) = -dt*((alpha1_d(i,j,Nz-1)+alpha1_d(i,j,Nz))/2.*max(w(i,j,Nz),0.) + &
                      ((alpha1_d(i,j,Nz-1)+alpha1_d(i,j,Nz))/2.*(Dld_eff(i,j,Nz-1)+Dld_eff(i,j,Nz))/2. + &
                      (alpha2_d(i,j,Nz-1)+alpha2_d(i,j,Nz))/2.*(Dgd_eff(i,j,Nz-1)+Dgd_eff(i,j,Nz))/2.)/dz)* &
                      (ne(i,j,Nz-1)+ne(i,j,Nz))/2./dz
                 gg(i,j,Nz) = 0.

                 if(tdbc_t.eq.1) then
                    aa(i,j,Nz) = ne(i,j,Nz) + dt*(-(alpha1_d(i-1,j,Nz)+alpha1_d(i,j,Nz))/2.*min(u(i,j,Nz),0.)* &
                         (ne(i-1,j,Nz)+ne(i,j,Nz))/2./dx + &
                         ((alpha1_d(i-1,j,Nz)+alpha1_d(i,j,Nz))/2.*(Dld_eff(i-1,j,Nz)+Dld_eff(i,j,Nz))/2. + &
                         (alpha2_d(i-1,j,Nz)+alpha2_d(i,j,Nz))/2.*(Dgd_eff(i-1,j,Nz)+Dgd_eff(i,j,Nz))/2.)/dx* &
                         (ne(i-1,j,Nz)+ne(i,j,Nz))/2./dx + &
                         (alpha1_d(i+1,j,Nz)+alpha1_d(i,j,Nz))/2.*max(u(i+1,j,Nz),0.)* &
                         (ne(i+1,j,Nz)+ne(i,j,Nz))/2./dx + &
                         ((alpha1_d(i+1,j,Nz)+alpha1_d(i,j,Nz))/2.*(Dld_eff(i+1,j,Nz)+Dld_eff(i,j,Nz))/2. + &
                         (alpha2_d(i+1,j,Nz)+alpha2_d(i,j,Nz))/2.*(Dgd_eff(i+1,j,Nz)+Dgd_eff(i,j,Nz))/2.)/dx* &
                         (ne(i+1,j,Nz)+ne(i,j,Nz))/2./dx - &
                         (alpha1_d(i,j-1,Nz)+alpha1_d(i,j,Nz))/2.*min(v(i,j,Nz),0.)*(ne(i,j-1,Nz)+ne(i,j,Nz))/2./dy + &
                         ((alpha1_d(i,j-1,Nz)+alpha1_d(i,j,Nz))/2.*(Dld_eff(i,j-1,Nz)+Dld_eff(i,j,Nz))/2. + &
                         (alpha2_d(i,j-1,Nz)+alpha2_d(i,j,Nz))/2.*(Dgd_eff(i,j-1,Nz)+Dgd_eff(i,j,Nz))/2.)/dy* &
                         (ne(i,j-1,Nz)+ne(i,j,Nz))/2./dy + &
                         (alpha1_d(i,j+1,Nz)+alpha1_d(i,j,Nz))/2.*max(v(i,j+1,Nz),0.)*(ne(i,j+1,Nz)+ne(i,j,Nz))/2./dy + &
                         ((alpha1_d(i,j+1,Nz)+alpha1_d(i,j,Nz))/2.*(Dld_eff(i,j+1,Nz)+Dld_eff(i,j,Nz))/2. + &
                         (alpha2_d(i,j+1,Nz)+alpha2_d(i,j,Nz))/2.*(Dgd_eff(i,j+1,Nz)+Dgd_eff(i,j,Nz))/2.)/dy* &
                         (ne(i,j+1,Nz)+ne(i,j,Nz))/2./dy - &
                         (alpha1_d(i,j,Nz-1)+alpha1_d(i,j,Nz))/2.*min(w(i,j,Nz),0.)*(ne(i,j,Nz-1)+ne(i,j,Nz))/2./dz + &
                         ((alpha1_d(i,j,Nz-1)+alpha1_d(i,j,Nz))/2.*(Dld_eff(i,j,Nz-1)+Dld_eff(i,j,Nz))/2. + &
                         (alpha2_d(i,j,Nz-1)+alpha2_d(i,j,Nz))/2.*(Dgd_eff(i,j,Nz-1)+Dgd_eff(i,j,Nz))/2.)/dz* &
                         (ne(i,j,Nz-1)+ne(i,j,Nz))/2./dz + &
                         (alpha1_d(i,j,Nz)*Dld_eff(i,j,Nz) + alpha2_d(i,j,Nz)*Dgd_eff(i,j,Nz))/dz*2.*ne(i,j,Nz)/dz)
                    q(i,j,Nz) = ne(i,j,Nz)*td_o(i,j,Nz) + &
                         dt*theta(i,j,Nz)* &
                         (qd_0*b_o(i,j,Nz)*lc_o(i,j,Nz)/(kc+lc_o(i,j,Nz))*lo_o(i,j,Nz)/(ko+lo_o(i,j,Nz))) - &
                         dt * (alpha1_d(i,j,Nz)*min(w(i,j,Nz+1),0.) - &
                         2.*(alpha1_d(i,j,Nz)*Dld_eff(i,j,Nz)+alpha2_d(i,j,Nz)*Dgd_eff(i,j,Nz))/dz)*ne(i,j,Nz)/dz*td_t(i,j)
                 endif
                 if(tdbc_t.eq.2) then
                    aa(i,j,Nz) = ne(i,j,Nz) + dt*(-(alpha1_d(i-1,j,Nz)+alpha1_d(i,j,Nz))/2.*min(u(i,j,Nz),0.)* &
                         (ne(i-1,j,Nz)+ne(i,j,Nz))/2./dx + &
                         ((alpha1_d(i-1,j,Nz)+alpha1_d(i,j,Nz))/2.*(Dld_eff(i-1,j,Nz)+Dld_eff(i,j,Nz))/2. + &
                         (alpha2_d(i-1,j,Nz)+alpha2_d(i,j,Nz))/2.*(Dgd_eff(i-1,j,Nz)+Dgd_eff(i,j,Nz))/2.)/dx* &
                         (ne(i-1,j,Nz)+ne(i,j,Nz))/2./dx + &
                         (alpha1_d(i+1,j,Nz)+alpha1_d(i,j,Nz))/2.*max(u(i+1,j,Nz),0.)* &
                         (ne(i+1,j,Nz)+ne(i,j,Nz))/2./dx + &
                         ((alpha1_d(i+1,j,Nz)+alpha1_d(i,j,Nz))/2.*(Dld_eff(i+1,j,Nz)+Dld_eff(i,j,Nz))/2. + &
                         (alpha2_d(i+1,j,Nz)+alpha2_d(i,j,Nz))/2.*(Dgd_eff(i+1,j,Nz)+Dgd_eff(i,j,Nz))/2.)/dx* &
                         (ne(i+1,j,Nz)+ne(i,j,Nz))/2./dx - &
                         (alpha1_d(i,j-1,Nz)+alpha1_d(i,j,Nz))/2.*min(v(i,j,Nz),0.)*(ne(i,j-1,Nz)+ne(i,j,Nz))/2./dy + &
                         ((alpha1_d(i,j-1,Nz)+alpha1_d(i,j,Nz))/2.*(Dld_eff(i,j-1,Nz)+Dld_eff(i,j,Nz))/2. + &
                         (alpha2_d(i,j-1,Nz)+alpha2_d(i,j,Nz))/2.*(Dgd_eff(i,j-1,Nz)+Dgd_eff(i,j,Nz))/2.)/dy* &
                         (ne(i,j-1,Nz)+ne(i,j,Nz))/2./dy + &
                         (alpha1_d(i,j+1,Nz)+alpha1_d(i,j,Nz))/2.*max(v(i,j+1,Nz),0.)*(ne(i,j+1,Nz)+ne(i,j,Nz))/2./dy + &
                         ((alpha1_d(i,j+1,Nz)+alpha1_d(i,j,Nz))/2.*(Dld_eff(i,j+1,Nz)+Dld_eff(i,j,Nz))/2. + &
                         (alpha2_d(i,j+1,Nz)+alpha2_d(i,j,Nz))/2.*(Dgd_eff(i,j+1,Nz)+Dgd_eff(i,j,Nz))/2.)/dy* &
                         (ne(i,j+1,Nz)+ne(i,j,Nz))/2./dy - &
                         (alpha1_d(i,j,Nz-1)+alpha1_d(i,j,Nz))/2.*min(w(i,j,Nz),0.)*(ne(i,j,Nz-1)+ne(i,j,Nz))/2./dz + &
                         ((alpha1_d(i,j,Nz-1)+alpha1_d(i,j,Nz))/2.*(Dld_eff(i,j,Nz-1)+Dld_eff(i,j,Nz))/2. + &
                         (alpha2_d(i,j,Nz-1)+alpha2_d(i,j,Nz))/2.*(Dgd_eff(i,j,Nz-1)+Dgd_eff(i,j,Nz))/2.)/dz* &
                         (ne(i,j,Nz-1)+ne(i,j,Nz))/2./dz + &
                         alpha1_d(i,j,Nz)*max(w(i,j,Nz+1),0.)*ne(i,j,Nz)/dz)   !!!?????
                    q(i,j,Nz) = ne(i,j,Nz)*td_o(i,j,Nz) + &
                         dt*theta(i,j,Nz)* &
                         (qd_0*b_o(i,j,Nz)*lc_o(i,j,Nz)/(kc+lc_o(i,j,Nz))*lo_o(i,j,Nz)/(ko+lo_o(i,j,Nz)))
                 endif


              endif
           enddo
        enddo


        !! for inside

        Do i=2,Nx-1
           Do j=2,Ny-1
              Do k=2,Nz-1
                 if(wall(i,j,k)) then
                    aa(i,j,k) = 0.
                    bb(i,j,k) = 0.
                    cc(i,j,k) = 0.
                    dd(i,j,k) = 0.
                    ee(i,j,k) = 0.
                    ff(i,j,k) = 0.
                    gg(i,j,k) = 0.
                    q(i,j,k) = 0.
                 else
                    aa(i,j,k) = ne(i,j,k) + dt*(-(alpha1_d(i-1,j,k)+alpha1_d(i,j,k))/2.*min(u(i,j,k),0.)* &
                         (ne(i-1,j,k)+ne(i,j,k))/2./dx + &
                         ((alpha1_d(i-1,j,k)+alpha1_d(i,j,k))/2.*(Dld_eff(i-1,j,k)+Dld_eff(i,j,k))/2. + &
                         (alpha2_d(i-1,j,k)+alpha2_d(i,j,k))/2.*(Dgd_eff(i-1,j,k)+Dgd_eff(i,j,k))/2.)/dx* &
                         (ne(i-1,j,k)+ne(i,j,k))/2./dx + &
                         (alpha1_d(i+1,j,k)+alpha1_d(i,j,k))/2.*max(u(i+1,j,k),0.)* &
                         (ne(i+1,j,k)+ne(i,j,k))/2./dx + &
                         ((alpha1_d(i+1,j,k)+alpha1_d(i,j,k))/2.*(Dld_eff(i+1,j,k)+Dld_eff(i,j,k))/2. + &
                         (alpha2_d(i+1,j,k)+alpha2_d(i,j,k))/2.*(Dgd_eff(i+1,j,k)+Dgd_eff(i,j,k))/2.)/dx* &
                         (ne(i+1,j,k)+ne(i,j,k))/2./dx - &
                         (alpha1_d(i,j-1,k)+alpha1_d(i,j,k))/2.*min(v(i,j,k),0.)*(ne(i,j-1,k)+ne(i,j,k))/2./dy + &
                         ((alpha1_d(i,j-1,k)+alpha1_d(i,j,k))/2.*(Dld_eff(i,j-1,k)+Dld_eff(i,j,k))/2. + &
                         (alpha2_d(i,j-1,k)+alpha2_d(i,j,k))/2.*(Dgd_eff(i,j-1,k)+Dgd_eff(i,j,k))/2.)/dy* &
                         (ne(i,j-1,k)+ne(i,j,k))/2./dy + &
                         (alpha1_d(i,j+1,k)+alpha1_d(i,j,k))/2.*max(v(i,j+1,k),0.)*(ne(i,j+1,k)+ne(i,j,k))/2./dy + &
                         ((alpha1_d(i,j+1,k)+alpha1_d(i,j,k))/2.*(Dld_eff(i,j+1,k)+Dld_eff(i,j,k))/2. + &
                         (alpha2_d(i,j+1,k)+alpha2_d(i,j,k))/2.*(Dgd_eff(i,j+1,k)+Dgd_eff(i,j,k))/2.)/dy* &
                         (ne(i,j+1,k)+ne(i,j,k))/2./dy - &
                         (alpha1_d(i,j,k-1)+alpha1_d(i,j,k))/2.*min(w(i,j,k),0.)*(ne(i,j,k-1)+ne(i,j,k))/2./dz + &
                         ((alpha1_d(i,j,k-1)+alpha1_d(i,j,k))/2.*(Dld_eff(i,j,k-1)+Dld_eff(i,j,k))/2. + &
                         (alpha2_d(i,j,k-1)+alpha2_d(i,j,k))/2.*(Dgd_eff(i,j,k-1)+Dgd_eff(i,j,k))/2.)/dz* &
                         (ne(i,j,k-1)+ne(i,j,k))/2./dz + &
                         (alpha1_d(i,j,k+1)+alpha1_d(i,j,k))/2.*max(w(i,j,k+1),0.)*(ne(i,j,k+1)+ne(i,j,k))/2./dz + &
                         ((alpha1_d(i,j,k+1)+alpha1_d(i,j,k))/2.*(Dld_eff(i,j,k+1)+Dld_eff(i,j,k))/2. + &
                         (alpha2_d(i,j,k+1)+alpha2_d(i,j,k))/2.*(Dgd_eff(i,j,k+1)+Dgd_eff(i,j,k))/2.)/dz* &
                         (ne(i,j,k+1)+ne(i,j,k))/2./dz)
                    bb(i,j,k) = -dt*((alpha1_d(i-1,j,k)+alpha1_d(i,j,k))/2.*max(u(i,j,k),0.) + &
                         ((alpha1_d(i-1,j,k)+alpha1_d(i,j,k))/2.*(Dld_eff(i-1,j,k)+Dld_eff(i,j,k))/2. + &
                         (alpha2_d(i-1,j,k)+alpha2_d(i,j,k))/2.*(Dgd_eff(i-1,j,k)+Dgd_eff(i,j,k))/2.)/dx)* &
                         (ne(i-1,j,k)+ne(i,j,k))/2./dx
                    cc(i,j,k) = dt*((alpha1_d(i+1,j,k)+alpha1_d(i,j,k))/2.*min(u(i+1,j,k),0.) - &
                         ((alpha1_d(i+1,j,k)+alpha1_d(i,j,k))/2.*(Dld_eff(i+1,j,k)+Dld_eff(i,j,k))/2. + &
                         (alpha2_d(i+1,j,k)+alpha2_d(i,j,k))/2.*(Dgd_eff(i+1,j,k)+Dgd_eff(i,j,k))/2.)/dx)* &
                         (ne(i+1,j,k)+ne(i,j,k))/2./dx
                    dd(i,j,k) = -dt*((alpha1_d(i,j-1,k)+alpha1_d(i,j,k))/2.*max(v(i,j,k),0.) + &
                         ((alpha1_d(i,j-1,k)+alpha1_d(i,j,k))/2.*(Dld_eff(i,j-1,k)+Dld_eff(i,j,k))/2. + &
                         (alpha2_d(i,j-1,k)+alpha2_d(i,j,k))/2.*(Dgd_eff(i,j-1,k)+Dgd_eff(i,j,k))/2.)/dy)* &
                         (ne(i,j-1,k)+ne(i,j,k))/2./dy
                    ee(i,j,k) = dt*((alpha1_d(i,j+1,k)+alpha1_d(i,j,k))/2.*min(v(i,j+1,k),0.) - &
                         ((alpha1_d(i,j+1,k)+alpha1_d(i,j,k))/2.*(Dld_eff(i,j+1,k)+Dld_eff(i,j,k))/2. + &
                         (alpha2_d(i,j+1,k)+alpha2_d(i,j,k))/2.*(Dgd_eff(i,j+1,k)+Dgd_eff(i,j,k))/2.)/dy)* &
                         (ne(i,j+1,k)+ne(i,j,k))/2./dy
                    ff(i,j,k) = -dt*((alpha1_d(i,j,k-1)+alpha1_d(i,j,k))/2.*max(w(i,j,k),0.) + &
                         ((alpha1_d(i,j,k-1)+alpha1_d(i,j,k))/2.*(Dld_eff(i,j,k-1)+Dld_eff(i,j,k))/2. + &
                         (alpha2_d(i,j,k-1)+alpha2_d(i,j,k))/2.*(Dgd_eff(i,j,k-1)+Dgd_eff(i,j,k))/2.)/dz)* &
                         (ne(i,j,k-1)+ne(i,j,k))/2./dz
                    gg(i,j,k) = dt*((alpha1_d(i,j,k+1)+alpha1_d(i,j,k))/2.*min(w(i,j,k+1),0.) - &
                         ((alpha1_d(i,j,k+1)+alpha1_d(i,j,k))/2.*(Dld_eff(i,j,k+1)+Dld_eff(i,j,k))/2. + &
                         (alpha2_d(i,j,k+1)+alpha2_d(i,j,k))/2.*(Dgd_eff(i,j,k+1)+Dgd_eff(i,j,k))/2.)/dz)* &
                         (ne(i,j,k+1)+ne(i,j,k))/2./dz

                    q(i,j,k) = ne(i,j,k)*td_o(i,j,k) + &
                         dt*theta(i,j,k)* &
                         (qd_0*b_o(i,j,k)*lc_o(i,j,k)/(Kc+lc_o(i,j,k))*lo_o(i,j,k)/(Ko+lo_o(i,j,k)))
                 endif
              Enddo
           Enddo
        Enddo


     else 
        if(mype.eq.0) then

           Do i=2,Nx-1
              Do j=2,Ny-1
                 if(wall(i,j,1)) then
                    aa(i,j,1) = 0.
                    bb(i,j,1) = 0.
                    cc(i,j,1) = 0.
                    dd(i,j,1) = 0.
                    ee(i,j,1) = 0.
                    ff(i,j,1) = 0.
                    gg(i,j,1) = 0.
                    q(i,j,1) = 0.
                 else
                 bb(i,j,1) = -dt*((alpha1_d(i-1,j,1)+alpha1_d(i,j,1))/2.*max(u(i,j,1),0.) + &
                      ((alpha1_d(i-1,j,1)+alpha1_d(i,j,1))/2.*(Dld_eff(i-1,j,1)+Dld_eff(i,j,1))/2. + &
                      (alpha2_d(i-1,j,1)+alpha2_d(i,j,1))/2.*(Dgd_eff(i-1,j,1)+Dgd_eff(i,j,1))/2.)/dx)* &
                      (ne(i-1,j,1)+ne(i,j,1))/2./dx
                 cc(i,j,1) = dt*((alpha1_d(i+1,j,1)+alpha1_d(i,j,1))/2.*min(u(i+1,j,1),0.) - &
                      ((alpha1_d(i+1,j,1)+alpha1_d(i,j,1))/2.*(Dld_eff(i+1,j,1)+Dld_eff(i,j,1))/2. + &
                      (alpha2_d(i+1,j,1)+alpha2_d(i,j,1))/2.*(Dgd_eff(i+1,j,1)+Dgd_eff(i,j,1))/2.)/dx)* &
                      (ne(i+1,j,1)+ne(i,j,1))/2./dx
                 dd(i,j,1) = -dt*((alpha1_d(i,j-1,1)+alpha1_d(i,j,1))/2.*max(v(i,j,1),0.) + &
                      ((alpha1_d(i,j-1,1)+alpha1_d(i,j,1))/2.*(Dld_eff(i,j-1,1)+Dld_eff(i,j,1))/2. + &
                      (alpha2_d(i,j-1,1)+alpha2_d(i,j,1))/2.*(Dgd_eff(i,j-1,1)+Dgd_eff(i,j,1))/2.)/dy)* &
                      (ne(i,j-1,1)+ne(i,j,1))/2./dy
                 ee(i,j,1) = dt*((alpha1_d(i,j+1,1)+alpha1_d(i,j,1))/2.*min(v(i,j+1,1),0.) - &
                      ((alpha1_d(i,j+1,1)+alpha1_d(i,j,1))/2.*(Dld_eff(i,j+1,1)+Dld_eff(i,j,1))/2. + &
                      (alpha2_d(i,j+1,1)+alpha2_d(i,j,1))/2.*(Dgd_eff(i,j+1,1)+Dgd_eff(i,j,1))/2.)/dy)* &
                      (ne(i,j+1,1)+ne(i,j,1))/2./dy
                 ff(i,j,1) = 0.
                 gg(i,j,1) = dt*((alpha1_d(i,j,1+1)+alpha1_d(i,j,1))/2.*min(w(i,j,1+1),0.) - &
                      ((alpha1_d(i,j,1+1)+alpha1_d(i,j,1))/2.*(Dld_eff(i,j,1+1)+Dld_eff(i,j,1))/2. + &
                      (alpha2_d(i,j,1+1)+alpha2_d(i,j,1))/2.*(Dgd_eff(i,j,1+1)+Dgd_eff(i,j,1))/2.)/dz)* &
                      (ne(i,j,1+1)+ne(i,j,1))/2./dz


                 if(tdbc_d.eq.1) then 
                    aa(i,j,1) = ne(i,j,1) + dt*(-(alpha1_d(i-1,j,1)+alpha1_d(i,j,1))/2.*min(u(i,j,1),0.)* &
                         (ne(i-1,j,1)+ne(i,j,1))/2./dx + &
                         ((alpha1_d(i-1,j,1)+alpha1_d(i,j,1))/2.*(Dld_eff(i-1,j,1)+Dld_eff(i,j,1))/2. + &
                         (alpha2_d(i-1,j,1)+alpha2_d(i,j,1))/2.*(Dgd_eff(i-1,j,1)+Dgd_eff(i,j,1))/2.)/dx* &
                         (ne(i-1,j,1)+ne(i,j,1))/2./dx + &
                         (alpha1_d(i+1,j,1)+alpha1_d(i,j,1))/2.*max(u(i+1,j,1),0.)* &
                         (ne(i+1,j,1)+ne(i,j,1))/2./dx + &
                         ((alpha1_d(i+1,j,1)+alpha1_d(i,j,1))/2.*(Dld_eff(i+1,j,1)+Dld_eff(i,j,1))/2. + &
                         (alpha2_d(i+1,j,1)+alpha2_d(i,j,1))/2.*(Dgd_eff(i+1,j,1)+Dgd_eff(i,j,1))/2.)/dx* &
                         (ne(i+1,j,1)+ne(i,j,1))/2./dx - &
                         (alpha1_d(i,j-1,1)+alpha1_d(i,j,1))/2.*min(v(i,j,1),0.)*(ne(i,j-1,1)+ne(i,j,1))/2./dy + &
                         ((alpha1_d(i,j-1,1)+alpha1_d(i,j,1))/2.*(Dld_eff(i,j-1,1)+Dld_eff(i,j,1))/2. + &
                         (alpha2_d(i,j-1,1)+alpha2_d(i,j,1))/2.*(Dgd_eff(i,j-1,1)+Dgd_eff(i,j,1))/2.)/dy* &
                         (ne(i,j-1,1)+ne(i,j,1))/2./dy + &
                         (alpha1_d(i,j+1,1)+alpha1_d(i,j,1))/2.*max(v(i,j+1,1),0.)*(ne(i,j+1,1)+ne(i,j,1))/2./dy + &
                         ((alpha1_d(i,j+1,1)+alpha1_d(i,j,1))/2.*(Dld_eff(i,j+1,1)+Dld_eff(i,j,1))/2. + &
                         (alpha2_d(i,j+1,1)+alpha2_d(i,j,1))/2.*(Dgd_eff(i,j+1,1)+Dgd_eff(i,j,1))/2.)/dy* &
                         (ne(i,j+1,1)+ne(i,j,1))/2./dy + &
                         (alpha1_d(i,j,1)*Dld_eff(i,j,1)+alpha2_d(i,j,1)*Dgd_eff(i,j,1))/dz*2.*ne(i,j,1)/dz + &
                         (alpha1_d(i,j,1+1)+alpha1_d(i,j,1))/2.*max(w(i,j,1+1),0.)*(ne(i,j,1+1)+ne(i,j,1))/2./dz + &
                         ((alpha1_d(i,j,1+1)+alpha1_d(i,j,1))/2.*(Dld_eff(i,j,1+1)+Dld_eff(i,j,1))/2. + &
                         (alpha2_d(i,j,1+1)+alpha2_d(i,j,1))/2.*(Dgd_eff(i,j,1+1)+Dgd_eff(i,j,1))/2.)/dz* &
                         (ne(i,j,1+1)+ne(i,j,1))/2./dz)
                    q(i,j,1) = ne(i,j,1)*td_o(i,j,1)  + &
                            dt*theta(i,j,1)* &
                            (qd_0*b_o(i,j,1)*lc_o(i,j,1)/(Kc+lc_o(i,j,1))*lo_o(i,j,1)/(Ko+lo_o(i,j,1))) - &
                            (-dt*(alpha1_d(i,j,1)*max(w(i,j,1),0.) + &
                            2.*(alpha1_d(i,j,1)*Dld_eff(i,j,1)+alpha2_d(i,j,1)*Dgd_eff(i,j,1))/dz)* &
                            ne(i,j,1)/dz)*td_d(i,j)
                 endif
                 if(tdbc_d.eq.2) then 
                    aa(i,j,1) = ne(i,j,1) + dt*(-(alpha1_d(i-1,j,1)+alpha1_d(i,j,1))/2.*min(u(i,j,1),0.)* &
                         (ne(i-1,j,1)+ne(i,j,1))/2./dx + &
                         ((alpha1_d(i-1,j,1)+alpha1_d(i,j,1))/2.*(Dld_eff(i-1,j,1)+Dld_eff(i,j,1))/2. + &
                         (alpha2_d(i-1,j,1)+alpha2_d(i,j,1))/2.*(Dgd_eff(i-1,j,1)+Dgd_eff(i,j,1))/2.)/dx* &
                         (ne(i-1,j,1)+ne(i,j,1))/2./dx + &
                         (alpha1_d(i+1,j,1)+alpha1_d(i,j,1))/2.*max(u(i+1,j,1),0.)* &
                         (ne(i+1,j,1)+ne(i,j,1))/2./dx + &
                         ((alpha1_d(i+1,j,1)+alpha1_d(i,j,1))/2.*(Dld_eff(i+1,j,1)+Dld_eff(i,j,1))/2. + &
                         (alpha2_d(i+1,j,1)+alpha2_d(i,j,1))/2.*(Dgd_eff(i+1,j,1)+Dgd_eff(i,j,1))/2.)/dx* &
                         (ne(i+1,j,1)+ne(i,j,1))/2./dx - &
                         (alpha1_d(i,j-1,1)+alpha1_d(i,j,1))/2.*min(v(i,j,1),0.)*(ne(i,j-1,1)+ne(i,j,1))/2./dy + &
                         ((alpha1_d(i,j-1,1)+alpha1_d(i,j,1))/2.*(Dld_eff(i,j-1,1)+Dld_eff(i,j,1))/2. + &
                         (alpha2_d(i,j-1,1)+alpha2_d(i,j,1))/2.*(Dgd_eff(i,j-1,1)+Dgd_eff(i,j,1))/2.)/dy* &
                         (ne(i,j-1,1)+ne(i,j,1))/2./dy + &
                         (alpha1_d(i,j+1,1)+alpha1_d(i,j,1))/2.*max(v(i,j+1,1),0.)*(ne(i,j+1,1)+ne(i,j,1))/2./dy + &
                         ((alpha1_d(i,j+1,1)+alpha1_d(i,j,1))/2.*(Dld_eff(i,j+1,1)+Dld_eff(i,j,1))/2. + &
                         (alpha2_d(i,j+1,1)+alpha2_d(i,j,1))/2.*(Dgd_eff(i,j+1,1)+Dgd_eff(i,j,1))/2.)/dy* &
                         (ne(i,j+1,1)+ne(i,j,1))/2./dy - &
                         alpha1_d(i,j,1)*min(w(i,j,1),0.)*ne(i,j,1)/dz + &
                         (alpha1_d(i,j,1+1)+alpha1_d(i,j,1))/2.*max(w(i,j,1+1),0.)*(ne(i,j,1+1)+ne(i,j,1))/2./dz + &
                         ((alpha1_d(i,j,1+1)+alpha1_d(i,j,1))/2.*(Dld_eff(i,j,1+1)+Dld_eff(i,j,1))/2. + &
                         (alpha2_d(i,j,1+1)+alpha2_d(i,j,1))/2.*(Dgd_eff(i,j,1+1)+Dgd_eff(i,j,1))/2.)/dz* &
                         (ne(i,j,1+1)+ne(i,j,1))/2./dz)
                    q(i,j,1) = ne(i,j,1)*td_o(i,j,1) + &
                         dt*theta(i,j,1)* &
                         (qd_0*b_o(i,j,1)*lc_o(i,j,1)/(Kc+lc_o(i,j,1))*lo_o(i,j,1)/(Ko+lo_o(i,j,1)))
                    endif

                 endif
              enddo
           enddo


           Do i=2,Nx-1
              Do j=2,Ny-1
                    if(wall(i,j,Nz)) then
                       aa(i,j,Nz) = 0.
                       bb(i,j,Nz) = 0.
                       cc(i,j,Nz) = 0.
                       dd(i,j,Nz) = 0.
                       ee(i,j,Nz) = 0.
                       ff(i,j,Nz) = 0.
                       gg(i,j,Nz) = 0.
                       q(i,j,Nz) = 0.
                    else
                       aa(i,j,Nz) = ne(i,j,Nz) + dt*(-(alpha1_d(i-1,j,Nz)+alpha1_d(i,j,Nz))/2.*min(u(i,j,Nz),0.)* &
                            (ne(i-1,j,Nz)+ne(i,j,Nz))/2./dx + &
                            ((alpha1_d(i-1,j,Nz)+alpha1_d(i,j,Nz))/2.*(Dld_eff(i-1,j,Nz)+Dld_eff(i,j,Nz))/2. + &
                            (alpha2_d(i-1,j,Nz)+alpha2_d(i,j,Nz))/2.*(Dgd_eff(i-1,j,Nz)+Dgd_eff(i,j,Nz))/2.)/dx* &
                            (ne(i-1,j,Nz)+ne(i,j,Nz))/2./dx + &
                            (alpha1_d(i+1,j,Nz)+alpha1_d(i,j,Nz))/2.*max(u(i+1,j,Nz),0.)* &
                            (ne(i+1,j,Nz)+ne(i,j,Nz))/2./dx + &
                            ((alpha1_d(i+1,j,Nz)+alpha1_d(i,j,Nz))/2.*(Dld_eff(i+1,j,Nz)+Dld_eff(i,j,Nz))/2. + &
                            (alpha2_d(i+1,j,Nz)+alpha2_d(i,j,Nz))/2.*(Dgd_eff(i+1,j,Nz)+Dgd_eff(i,j,Nz))/2.)/dx* &
                            (ne(i+1,j,Nz)+ne(i,j,Nz))/2./dx - &
                            (alpha1_d(i,j-1,Nz)+alpha1_d(i,j,Nz))/2.*min(v(i,j,Nz),0.)*(ne(i,j-1,Nz)+ne(i,j,Nz))/2./dy + &
                            ((alpha1_d(i,j-1,Nz)+alpha1_d(i,j,Nz))/2.*(Dld_eff(i,j-1,Nz)+Dld_eff(i,j,Nz))/2. + &
                            (alpha2_d(i,j-1,Nz)+alpha2_d(i,j,Nz))/2.*(Dgd_eff(i,j-1,Nz)+Dgd_eff(i,j,Nz))/2.)/dy* &
                            (ne(i,j-1,Nz)+ne(i,j,Nz))/2./dy + &
                            (alpha1_d(i,j+1,Nz)+alpha1_d(i,j,Nz))/2.*max(v(i,j+1,Nz),0.)*(ne(i,j+1,Nz)+ne(i,j,Nz))/2./dy + &
                            ((alpha1_d(i,j+1,Nz)+alpha1_d(i,j,Nz))/2.*(Dld_eff(i,j+1,Nz)+Dld_eff(i,j,Nz))/2. + &
                            (alpha2_d(i,j+1,Nz)+alpha2_d(i,j,Nz))/2.*(Dgd_eff(i,j+1,Nz)+Dgd_eff(i,j,Nz))/2.)/dy* &
                            (ne(i,j+1,Nz)+ne(i,j,Nz))/2./dy - &
                            (alpha1_d(i,j,Nz-1)+alpha1_d(i,j,Nz))/2.*min(w(i,j,Nz),0.)*(ne(i,j,Nz-1)+ne(i,j,Nz))/2./dz + &
                            ((alpha1_d(i,j,Nz-1)+alpha1_d(i,j,Nz))/2.*(Dld_eff(i,j,Nz-1)+Dld_eff(i,j,Nz))/2. + &
                            (alpha2_d(i,j,Nz-1)+alpha2_d(i,j,Nz))/2.*(Dgd_eff(i,j,Nz-1)+Dgd_eff(i,j,Nz))/2.)/dz* &
                            (ne(i,j,Nz-1)+ne(i,j,Nz))/2./dz + &
                            (tmp_alpha1_d(i,j,2)+alpha1_d(i,j,Nz))/2.*max(w(i,j,Nz+1),0.)*(tmp_ne(i,j,2)+ne(i,j,Nz))/2./dz + &
                            ((tmp_alpha1_d(i,j,2)+alpha1_d(i,j,Nz))/2.*(tmp_Dld_eff(i,j,2)+Dld_eff(i,j,Nz))/2. + &
                            (tmp_alpha2_d(i,j,2)+alpha2_d(i,j,Nz))/2.*(tmp_Dgd_eff(i,j,2)+Dgd_eff(i,j,Nz))/2.)/dz* &
                            (tmp_ne(i,j,2)+ne(i,j,Nz))/2./dz)
                       bb(i,j,Nz) = -dt*((alpha1_d(i-1,j,Nz)+alpha1_d(i,j,Nz))/2.*max(u(i,j,Nz),0.) + &
                            ((alpha1_d(i-1,j,Nz)+alpha1_d(i,j,Nz))/2.*(Dld_eff(i-1,j,Nz)+Dld_eff(i,j,Nz))/2. + &
                            (alpha2_d(i-1,j,Nz)+alpha2_d(i,j,Nz))/2.*(Dgd_eff(i-1,j,Nz)+Dgd_eff(i,j,Nz))/2.)/dx)* &
                            (ne(i-1,j,Nz)+ne(i,j,Nz))/2./dx
                       cc(i,j,Nz) = dt*((alpha1_d(i+1,j,Nz)+alpha1_d(i,j,Nz))/2.*min(u(i+1,j,Nz),0.) - &
                            ((alpha1_d(i+1,j,Nz)+alpha1_d(i,j,Nz))/2.*(Dld_eff(i+1,j,Nz)+Dld_eff(i,j,Nz))/2. + &
                            (alpha2_d(i+1,j,Nz)+alpha2_d(i,j,Nz))/2.*(Dgd_eff(i+1,j,Nz)+Dgd_eff(i,j,Nz))/2.)/dx)* &
                            (ne(i+1,j,Nz)+ne(i,j,Nz))/2./dx
                       dd(i,j,Nz) = -dt*((alpha1_d(i,j-1,Nz)+alpha1_d(i,j,Nz))/2.*max(v(i,j,Nz),0.) + &
                            ((alpha1_d(i,j-1,Nz)+alpha1_d(i,j,Nz))/2.*(Dld_eff(i,j-1,Nz)+Dld_eff(i,j,Nz))/2. + &
                            (alpha2_d(i,j-1,Nz)+alpha2_d(i,j,Nz))/2.*(Dgd_eff(i,j-1,Nz)+Dgd_eff(i,j,Nz))/2.)/dy)* &
                            (ne(i,j-1,Nz)+ne(i,j,Nz))/2./dy
                       ee(i,j,Nz) = dt*((alpha1_d(i,j+1,Nz)+alpha1_d(i,j,Nz))/2.*min(v(i,j+1,Nz),0.) - &
                            ((alpha1_d(i,j+1,Nz)+alpha1_d(i,j,Nz))/2.*(Dld_eff(i,j+1,Nz)+Dld_eff(i,j,Nz))/2. + &
                            (alpha2_d(i,j+1,Nz)+alpha2_d(i,j,Nz))/2.*(Dgd_eff(i,j+1,Nz)+Dgd_eff(i,j,Nz))/2.)/dy)* &
                            (ne(i,j+1,Nz)+ne(i,j,Nz))/2./dy
                       ff(i,j,Nz) = -dt*((alpha1_d(i,j,Nz-1)+alpha1_d(i,j,Nz))/2.*max(w(i,j,Nz),0.) + &
                             ((alpha1_d(i,j,Nz-1)+alpha1_d(i,j,Nz))/2.*(Dld_eff(i,j,Nz-1)+Dld_eff(i,j,Nz))/2. + &
                            (alpha2_d(i,j,Nz-1)+alpha2_d(i,j,Nz))/2.*(Dgd_eff(i,j,Nz-1)+Dgd_eff(i,j,Nz))/2.)/dz)* &
                            (ne(i,j,Nz-1)+ne(i,j,Nz))/2./dz
                       gg(i,j,Nz) = 0

                       q(i,j,Nz) = ne(i,j,Nz)*td_o(i,j,Nz) + &
                            dt*theta(i,j,Nz)* &
                            (qd_0*b_o(i,j,Nz)*lc_o(i,j,Nz)/(Kc+lc_o(i,j,Nz))*lo_o(i,j,Nz)/(Ko+lo_o(i,j,Nz))) - &
                            (dt*((tmp_alpha1_d(i,j,2)+alpha1_d(i,j,Nz))/2.*min(w(i,j,Nz+1),0.) - &
                            ((tmp_alpha1_d(i,j,2)+alpha1_d(i,j,Nz))/2.*(tmp_Dld_eff(i,j,2)+Dld_eff(i,j,Nz))/2. + &
                            (tmp_alpha2_d(i,j,2)+alpha2_d(i,j,Nz))/2.*(tmp_Dgd_eff(i,j,2)+Dgd_eff(i,j,Nz))/2.)/dz)* &
                            (tmp_ne(i,j,2)+ne(i,j,Nz))/2./dz)*tmp_td_o(i,j,2)
                    endif
                 Enddo
              Enddo


           !! for inside

           Do i=2,Nx-1
              Do j=2,Ny-1
                 Do k=2,Nz-1
                    if(wall(i,j,k)) then
                       aa(i,j,k) = 0.
                       bb(i,j,k) = 0.
                       cc(i,j,k) = 0.
                       dd(i,j,k) = 0.
                       ee(i,j,k) = 0.
                       ff(i,j,k) = 0.
                       gg(i,j,k) = 0.
                       q(i,j,k) = 0.
                    else
                       aa(i,j,k) = ne(i,j,k) + dt*(-(alpha1_d(i-1,j,k)+alpha1_d(i,j,k))/2.*min(u(i,j,k),0.)* &
                            (ne(i-1,j,k)+ne(i,j,k))/2./dx + &
                            ((alpha1_d(i-1,j,k)+alpha1_d(i,j,k))/2.*(Dld_eff(i-1,j,k)+Dld_eff(i,j,k))/2. + &
                            (alpha2_d(i-1,j,k)+alpha2_d(i,j,k))/2.*(Dgd_eff(i-1,j,k)+Dgd_eff(i,j,k))/2.)/dx* &
                            (ne(i-1,j,k)+ne(i,j,k))/2./dx + &
                            (alpha1_d(i+1,j,k)+alpha1_d(i,j,k))/2.*max(u(i+1,j,k),0.)* &
                            (ne(i+1,j,k)+ne(i,j,k))/2./dx + &
                            ((alpha1_d(i+1,j,k)+alpha1_d(i,j,k))/2.*(Dld_eff(i+1,j,k)+Dld_eff(i,j,k))/2. + &
                            (alpha2_d(i+1,j,k)+alpha2_d(i,j,k))/2.*(Dgd_eff(i+1,j,k)+Dgd_eff(i,j,k))/2.)/dx* &
                            (ne(i+1,j,k)+ne(i,j,k))/2./dx - &
                            (alpha1_d(i,j-1,k)+alpha1_d(i,j,k))/2.*min(v(i,j,k),0.)*(ne(i,j-1,k)+ne(i,j,k))/2./dy + &
                            ((alpha1_d(i,j-1,k)+alpha1_d(i,j,k))/2.*(Dld_eff(i,j-1,k)+Dld_eff(i,j,k))/2. + &
                            (alpha2_d(i,j-1,k)+alpha2_d(i,j,k))/2.*(Dgd_eff(i,j-1,k)+Dgd_eff(i,j,k))/2.)/dy* &
                            (ne(i,j-1,k)+ne(i,j,k))/2./dy + &
                            (alpha1_d(i,j+1,k)+alpha1_d(i,j,k))/2.*max(v(i,j+1,k),0.)*(ne(i,j+1,k)+ne(i,j,k))/2./dy + &
                            ((alpha1_d(i,j+1,k)+alpha1_d(i,j,k))/2.*(Dld_eff(i,j+1,k)+Dld_eff(i,j,k))/2. + &
                            (alpha2_d(i,j+1,k)+alpha2_d(i,j,k))/2.*(Dgd_eff(i,j+1,k)+Dgd_eff(i,j,k))/2.)/dy* &
                            (ne(i,j+1,k)+ne(i,j,k))/2./dy - &
                            (alpha1_d(i,j,k-1)+alpha1_d(i,j,k))/2.*min(w(i,j,k),0.)*(ne(i,j,k-1)+ne(i,j,k))/2./dz + &
                            ((alpha1_d(i,j,k-1)+alpha1_d(i,j,k))/2.*(Dld_eff(i,j,k-1)+Dld_eff(i,j,k))/2. + &
                            (alpha2_d(i,j,k-1)+alpha2_d(i,j,k))/2.*(Dgd_eff(i,j,k-1)+Dgd_eff(i,j,k))/2.)/dz* &
                            (ne(i,j,k-1)+ne(i,j,k))/2./dz + &
                            (alpha1_d(i,j,k+1)+alpha1_d(i,j,k))/2.*max(w(i,j,k+1),0.)*(ne(i,j,k+1)+ne(i,j,k))/2./dz + &
                            ((alpha1_d(i,j,k+1)+alpha1_d(i,j,k))/2.*(Dld_eff(i,j,k+1)+Dld_eff(i,j,k))/2. + &
                            (alpha2_d(i,j,k+1)+alpha2_d(i,j,k))/2.*(Dgd_eff(i,j,k+1)+Dgd_eff(i,j,k))/2.)/dz* &
                            (ne(i,j,k+1)+ne(i,j,k))/2./dz)
                       bb(i,j,k) = -dt*((alpha1_d(i-1,j,k)+alpha1_d(i,j,k))/2.*max(u(i,j,k),0.) + &
                            ((alpha1_d(i-1,j,k)+alpha1_d(i,j,k))/2.*(Dld_eff(i-1,j,k)+Dld_eff(i,j,k))/2. + &
                            (alpha2_d(i-1,j,k)+alpha2_d(i,j,k))/2.*(Dgd_eff(i-1,j,k)+Dgd_eff(i,j,k))/2.)/dx)* &
                            (ne(i-1,j,k)+ne(i,j,k))/2./dx
                       cc(i,j,k) = dt*((alpha1_d(i+1,j,k)+alpha1_d(i,j,k))/2.*min(u(i+1,j,k),0.) - &
                            ((alpha1_d(i+1,j,k)+alpha1_d(i,j,k))/2.*(Dld_eff(i+1,j,k)+Dld_eff(i,j,k))/2. + &
                            (alpha2_d(i+1,j,k)+alpha2_d(i,j,k))/2.*(Dgd_eff(i+1,j,k)+Dgd_eff(i,j,k))/2.)/dx)* &
                            (ne(i+1,j,k)+ne(i,j,k))/2./dx
                       dd(i,j,k) = -dt*((alpha1_d(i,j-1,k)+alpha1_d(i,j,k))/2.*max(v(i,j,k),0.) + &
                            ((alpha1_d(i,j-1,k)+alpha1_d(i,j,k))/2.*(Dld_eff(i,j-1,k)+Dld_eff(i,j,k))/2. + &
                            (alpha2_d(i,j-1,k)+alpha2_d(i,j,k))/2.*(Dgd_eff(i,j-1,k)+Dgd_eff(i,j,k))/2.)/dy)* &
                            (ne(i,j-1,k)+ne(i,j,k))/2./dy
                       ee(i,j,k) = dt*((alpha1_d(i,j+1,k)+alpha1_d(i,j,k))/2.*min(v(i,j+1,k),0.) - &
                            ((alpha1_d(i,j+1,k)+alpha1_d(i,j,k))/2.*(Dld_eff(i,j+1,k)+Dld_eff(i,j,k))/2. + &
                            (alpha2_d(i,j+1,k)+alpha2_d(i,j,k))/2.*(Dgd_eff(i,j+1,k)+Dgd_eff(i,j,k))/2.)/dy)* &
                            (ne(i,j+1,k)+ne(i,j,k))/2./dy
                       ff(i,j,k) = -dt*((alpha1_d(i,j,k-1)+alpha1_d(i,j,k))/2.*max(w(i,j,k),0.) + &
                             ((alpha1_d(i,j,k-1)+alpha1_d(i,j,k))/2.*(Dld_eff(i,j,k-1)+Dld_eff(i,j,k))/2. + &
                            (alpha2_d(i,j,k-1)+alpha2_d(i,j,k))/2.*(Dgd_eff(i,j,k-1)+Dgd_eff(i,j,k))/2.)/dz)* &
                            (ne(i,j,k-1)+ne(i,j,k))/2./dz
                       gg(i,j,k) = dt*((alpha1_d(i,j,k+1)+alpha1_d(i,j,k))/2.*min(w(i,j,k+1),0.) - &
                            ((alpha1_d(i,j,k+1)+alpha1_d(i,j,k))/2.*(Dld_eff(i,j,k+1)+Dld_eff(i,j,k))/2. + &
                            (alpha2_d(i,j,k+1)+alpha2_d(i,j,k))/2.*(Dgd_eff(i,j,k+1)+Dgd_eff(i,j,k))/2.)/dz)* &
                            (ne(i,j,k+1)+ne(i,j,k))/2./dz

                       q(i,j,k) = ne(i,j,k)*td_o(i,j,k) + &
                            dt*theta(i,j,k)* &
                            (qd_0*b_o(i,j,k)*lc_o(i,j,k)/(Kc+lc_o(i,j,k))*lo_o(i,j,k)/(Ko+lo_o(i,j,k)))
                    endif
                 Enddo
              Enddo
           Enddo



        else if(mype.eq.npes-1) then


           Do i=2,Nx-1
              Do j=2,Ny-1
                    if(wall(i,j,1)) then
                       aa(i,j,1) = 0.
                       bb(i,j,1) = 0.
                       cc(i,j,1) = 0.
                       dd(i,j,1) = 0.
                       ee(i,j,1) = 0.
                       ff(i,j,1) = 0.
                       gg(i,j,1) = 0.
                       q(i,j,1) = 0.
                    else
                       aa(i,j,1) = ne(i,j,1) + dt*(-(alpha1_d(i-1,j,1)+alpha1_d(i,j,1))/2.*min(u(i,j,1),0.)* &
                            (ne(i-1,j,1)+ne(i,j,1))/2./dx + &
                            ((alpha1_d(i-1,j,1)+alpha1_d(i,j,1))/2.*(Dld_eff(i-1,j,1)+Dld_eff(i,j,1))/2. + &
                            (alpha2_d(i-1,j,1)+alpha2_d(i,j,1))/2.*(Dgd_eff(i-1,j,1)+Dgd_eff(i,j,1))/2.)/dx* &
                            (ne(i-1,j,1)+ne(i,j,1))/2./dx + &
                            (alpha1_d(i+1,j,1)+alpha1_d(i,j,1))/2.*max(u(i+1,j,1),0.)* &
                            (ne(i+1,j,1)+ne(i,j,1))/2./dx + &
                            ((alpha1_d(i+1,j,1)+alpha1_d(i,j,1))/2.*(Dld_eff(i+1,j,1)+Dld_eff(i,j,1))/2. + &
                            (alpha2_d(i+1,j,1)+alpha2_d(i,j,1))/2.*(Dgd_eff(i+1,j,1)+Dgd_eff(i,j,1))/2.)/dx* &
                            (ne(i+1,j,1)+ne(i,j,1))/2./dx - &
                            (alpha1_d(i,j-1,1)+alpha1_d(i,j,1))/2.*min(v(i,j,1),0.)*(ne(i,j-1,1)+ne(i,j,1))/2./dy + &
                            ((alpha1_d(i,j-1,1)+alpha1_d(i,j,1))/2.*(Dld_eff(i,j-1,1)+Dld_eff(i,j,1))/2. + &
                            (alpha2_d(i,j-1,1)+alpha2_d(i,j,1))/2.*(Dgd_eff(i,j-1,1)+Dgd_eff(i,j,1))/2.)/dy* &
                            (ne(i,j-1,1)+ne(i,j,1))/2./dy + &
                            (alpha1_d(i,j+1,1)+alpha1_d(i,j,1))/2.*max(v(i,j+1,1),0.)*(ne(i,j+1,1)+ne(i,j,1))/2./dy + &
                            ((alpha1_d(i,j+1,1)+alpha1_d(i,j,1))/2.*(Dld_eff(i,j+1,1)+Dld_eff(i,j,1))/2. + &
                            (alpha2_d(i,j+1,1)+alpha2_d(i,j,1))/2.*(Dgd_eff(i,j+1,1)+Dgd_eff(i,j,1))/2.)/dy* &
                            (ne(i,j+1,1)+ne(i,j,1))/2./dy - &
                            (tmp_alpha1_d(i,j,1)+alpha1_d(i,j,1))/2.*min(w(i,j,1),0.)*(tmp_ne(i,j,1)+ne(i,j,1))/2./dz + &
                            ((tmp_alpha1_d(i,j,1)+alpha1_d(i,j,1))/2.*(tmp_Dld_eff(i,j,1)+Dld_eff(i,j,1))/2. + &
                            (tmp_alpha2_d(i,j,1)+alpha2_d(i,j,1))/2.*(tmp_Dgd_eff(i,j,1)+Dgd_eff(i,j,1))/2.)/dz* &
                            (tmp_ne(i,j,1)+ne(i,j,1))/2./dz + &
                            (alpha1_d(i,j,1+1)+alpha1_d(i,j,1))/2.*max(w(i,j,1+1),0.)*(ne(i,j,1+1)+ne(i,j,1))/2./dz + &
                            ((alpha1_d(i,j,1+1)+alpha1_d(i,j,1))/2.*(Dld_eff(i,j,1+1)+Dld_eff(i,j,1))/2. + &
                            (alpha2_d(i,j,1+1)+alpha2_d(i,j,1))/2.*(Dgd_eff(i,j,1+1)+Dgd_eff(i,j,1))/2.)/dz* &
                            (ne(i,j,1+1)+ne(i,j,1))/2./dz)
                       bb(i,j,1) = -dt*((alpha1_d(i-1,j,1)+alpha1_d(i,j,1))/2.*max(u(i,j,1),0.) + &
                            ((alpha1_d(i-1,j,1)+alpha1_d(i,j,1))/2.*(Dld_eff(i-1,j,1)+Dld_eff(i,j,1))/2. + &
                            (alpha2_d(i-1,j,1)+alpha2_d(i,j,1))/2.*(Dgd_eff(i-1,j,1)+Dgd_eff(i,j,1))/2.)/dx)* &
                            (ne(i-1,j,1)+ne(i,j,1))/2./dx
                       cc(i,j,1) = dt*((alpha1_d(i+1,j,1)+alpha1_d(i,j,1))/2.*min(u(i+1,j,1),0.) - &
                            ((alpha1_d(i+1,j,1)+alpha1_d(i,j,1))/2.*(Dld_eff(i+1,j,1)+Dld_eff(i,j,1))/2. + &
                            (alpha2_d(i+1,j,1)+alpha2_d(i,j,1))/2.*(Dgd_eff(i+1,j,1)+Dgd_eff(i,j,1))/2.)/dx)* &
                            (ne(i+1,j,1)+ne(i,j,1))/2./dx
                       dd(i,j,1) = -dt*((alpha1_d(i,j-1,1)+alpha1_d(i,j,1))/2.*max(v(i,j,1),0.) + &
                            ((alpha1_d(i,j-1,1)+alpha1_d(i,j,1))/2.*(Dld_eff(i,j-1,1)+Dld_eff(i,j,1))/2. + &
                            (alpha2_d(i,j-1,1)+alpha2_d(i,j,1))/2.*(Dgd_eff(i,j-1,1)+Dgd_eff(i,j,1))/2.)/dy)* &
                            (ne(i,j-1,1)+ne(i,j,1))/2./dy
                       ee(i,j,1) = dt*((alpha1_d(i,j+1,1)+alpha1_d(i,j,1))/2.*min(v(i,j+1,1),0.) - &
                            ((alpha1_d(i,j+1,1)+alpha1_d(i,j,1))/2.*(Dld_eff(i,j+1,1)+Dld_eff(i,j,1))/2. + &
                            (alpha2_d(i,j+1,1)+alpha2_d(i,j,1))/2.*(Dgd_eff(i,j+1,1)+Dgd_eff(i,j,1))/2.)/dy)* &
                            (ne(i,j+1,1)+ne(i,j,1))/2./dy
                       ff(i,j,1) = 0.
                       gg(i,j,1) = dt*((alpha1_d(i,j,1+1)+alpha1_d(i,j,1))/2.*min(w(i,j,1+1),0.) - &
                            ((alpha1_d(i,j,1+1)+alpha1_d(i,j,1))/2.*(Dld_eff(i,j,1+1)+Dld_eff(i,j,1))/2. + &
                            (alpha2_d(i,j,1+1)+alpha2_d(i,j,1))/2.*(Dgd_eff(i,j,1+1)+Dgd_eff(i,j,1))/2.)/dz)* &
                            (ne(i,j,1+1)+ne(i,j,1))/2./dz

                       q(i,j,1) = ne(i,j,1)*td_o(i,j,1) + &
                            dt*theta(i,j,1)* &
                            (qd_0*b_o(i,j,1)*lc_o(i,j,1)/(Kc+lc_o(i,j,1))*lo_o(i,j,1)/(Ko+lo_o(i,j,1))) - &
                            (-dt*((tmp_alpha1_d(i,j,1)+alpha1_d(i,j,1))/2.*max(w(i,j,1),0.) + &
                             ((tmp_alpha1_d(i,j,1)+alpha1_d(i,j,1))/2.*(tmp_Dld_eff(i,j,1)+Dld_eff(i,j,1))/2. + &
                            (tmp_alpha2_d(i,j,1)+alpha2_d(i,j,1))/2.*(tmp_Dgd_eff(i,j,1)+Dgd_eff(i,j,1))/2.)/dz)* &
                            (tmp_ne(i,j,1)+ne(i,j,1))/2./dz)*tmp_td_o(i,j,1)
                    endif
              Enddo
           Enddo
           
           
           
           Do i=2,Nx-1
              Do j=2,Ny-1
                 if(wall(i,j,Nz)) then
                    aa(i,j,Nz) = 0.
                    bb(i,j,Nz) = 0.
                    cc(i,j,Nz) = 0.
                    dd(i,j,Nz) = 0.
                    ee(i,j,Nz) = 0.
                    ff(i,j,Nz) = 0.
                    gg(i,j,Nz) = 0.
                    q(i,j,Nz) = 0.
                 else     
                 bb(i,j,Nz) = -dt*((alpha1_d(i-1,j,Nz)+alpha1_d(i,j,Nz))/2.*max(u(i,j,Nz),0.) + &
                      ((alpha1_d(i-1,j,Nz)+alpha1_d(i,j,Nz))/2.*(Dld_eff(i-1,j,Nz)+Dld_eff(i,j,Nz))/2. + &
                      (alpha2_d(i-1,j,Nz)+alpha2_d(i,j,Nz))/2.*(Dgd_eff(i-1,j,Nz)+Dgd_eff(i,j,Nz))/2.)/dx)* &
                      (ne(i-1,j,Nz)+ne(i,j,Nz))/2./dx
                 cc(i,j,Nz) = dt*((alpha1_d(i+1,j,Nz)+alpha1_d(i,j,Nz))/2.*min(u(i+1,j,Nz),0.) - &
                      ((alpha1_d(i+1,j,Nz)+alpha1_d(i,j,Nz))/2.*(Dld_eff(i+1,j,Nz)+Dld_eff(i,j,Nz))/2. + &
                      (alpha2_d(i+1,j,Nz)+alpha2_d(i,j,Nz))/2.*(Dgd_eff(i+1,j,Nz)+Dgd_eff(i,j,Nz))/2.)/dx)* &
                      (ne(i+1,j,Nz)+ne(i,j,Nz))/2./dx
                 dd(i,j,Nz) = -dt*((alpha1_d(i,j-1,Nz)+alpha1_d(i,j,Nz))/2.*max(v(i,j,Nz),0.) + &
                      ((alpha1_d(i,j-1,Nz)+alpha1_d(i,j,Nz))/2.*(Dld_eff(i,j-1,Nz)+Dld_eff(i,j,Nz))/2. + &
                      (alpha2_d(i,j-1,Nz)+alpha2_d(i,j,Nz))/2.*(Dgd_eff(i,j-1,Nz)+Dgd_eff(i,j,Nz))/2.)/dy)* &
                      (ne(i,j-1,Nz)+ne(i,j,Nz))/2./dy
                 ee(i,j,Nz) = dt*((alpha1_d(i,j+1,Nz)+alpha1_d(i,j,Nz))/2.*min(v(i,j+1,Nz),0.) - &
                      ((alpha1_d(i,j+1,Nz)+alpha1_d(i,j,Nz))/2.*(Dld_eff(i,j+1,Nz)+Dld_eff(i,j,Nz))/2. + &
                      (alpha2_d(i,j+1,Nz)+alpha2_d(i,j,Nz))/2.*(Dgd_eff(i,j+1,Nz)+Dgd_eff(i,j,Nz))/2.)/dy)* &
                      (ne(i,j+1,Nz)+ne(i,j,Nz))/2./dy
                 ff(i,j,Nz) = -dt*((alpha1_d(i,j,Nz-1)+alpha1_d(i,j,Nz))/2.*max(w(i,j,Nz),0.) + &
                      ((alpha1_d(i,j,Nz-1)+alpha1_d(i,j,Nz))/2.*(Dld_eff(i,j,Nz-1)+Dld_eff(i,j,Nz))/2. + &
                      (alpha2_d(i,j,Nz-1)+alpha2_d(i,j,Nz))/2.*(Dgd_eff(i,j,Nz-1)+Dgd_eff(i,j,Nz))/2.)/dz)* &
                      (ne(i,j,Nz-1)+ne(i,j,Nz))/2./dz
                 gg(i,j,Nz) = 0.

                 if(tdbc_t.eq.1) then
                    aa(i,j,Nz) = ne(i,j,Nz) + dt*(-(alpha1_d(i-1,j,Nz)+alpha1_d(i,j,Nz))/2.*min(u(i,j,Nz),0.)* &
                         (ne(i-1,j,Nz)+ne(i,j,Nz))/2./dx + &
                         ((alpha1_d(i-1,j,Nz)+alpha1_d(i,j,Nz))/2.*(Dld_eff(i-1,j,Nz)+Dld_eff(i,j,Nz))/2. + &
                         (alpha2_d(i-1,j,Nz)+alpha2_d(i,j,Nz))/2.*(Dgd_eff(i-1,j,Nz)+Dgd_eff(i,j,Nz))/2.)/dx* &
                         (ne(i-1,j,Nz)+ne(i,j,Nz))/2./dx + &
                         (alpha1_d(i+1,j,Nz)+alpha1_d(i,j,Nz))/2.*max(u(i+1,j,Nz),0.)* &
                         (ne(i+1,j,Nz)+ne(i,j,Nz))/2./dx + &
                         ((alpha1_d(i+1,j,Nz)+alpha1_d(i,j,Nz))/2.*(Dld_eff(i+1,j,Nz)+Dld_eff(i,j,Nz))/2. + &
                         (alpha2_d(i+1,j,Nz)+alpha2_d(i,j,Nz))/2.*(Dgd_eff(i+1,j,Nz)+Dgd_eff(i,j,Nz))/2.)/dx* &
                         (ne(i+1,j,Nz)+ne(i,j,Nz))/2./dx - &
                         (alpha1_d(i,j-1,Nz)+alpha1_d(i,j,Nz))/2.*min(v(i,j,Nz),0.)*(ne(i,j-1,Nz)+ne(i,j,Nz))/2./dy + &
                         ((alpha1_d(i,j-1,Nz)+alpha1_d(i,j,Nz))/2.*(Dld_eff(i,j-1,Nz)+Dld_eff(i,j,Nz))/2. + &
                         (alpha2_d(i,j-1,Nz)+alpha2_d(i,j,Nz))/2.*(Dgd_eff(i,j-1,Nz)+Dgd_eff(i,j,Nz))/2.)/dy* &
                         (ne(i,j-1,Nz)+ne(i,j,Nz))/2./dy + &
                         (alpha1_d(i,j+1,Nz)+alpha1_d(i,j,Nz))/2.*max(v(i,j+1,Nz),0.)*(ne(i,j+1,Nz)+ne(i,j,Nz))/2./dy + &
                         ((alpha1_d(i,j+1,Nz)+alpha1_d(i,j,Nz))/2.*(Dld_eff(i,j+1,Nz)+Dld_eff(i,j,Nz))/2. + &
                         (alpha2_d(i,j+1,Nz)+alpha2_d(i,j,Nz))/2.*(Dgd_eff(i,j+1,Nz)+Dgd_eff(i,j,Nz))/2.)/dy* &
                         (ne(i,j+1,Nz)+ne(i,j,Nz))/2./dy - &
                         (alpha1_d(i,j,Nz-1)+alpha1_d(i,j,Nz))/2.*min(w(i,j,Nz),0.)*(ne(i,j,Nz-1)+ne(i,j,Nz))/2./dz + &
                         ((alpha1_d(i,j,Nz-1)+alpha1_d(i,j,Nz))/2.*(Dld_eff(i,j,Nz-1)+Dld_eff(i,j,Nz))/2. + &
                         (alpha2_d(i,j,Nz-1)+alpha2_d(i,j,Nz))/2.*(Dgd_eff(i,j,Nz-1)+Dgd_eff(i,j,Nz))/2.)/dz* &
                         (ne(i,j,Nz-1)+ne(i,j,Nz))/2./dz + &
                         (alpha1_d(i,j,Nz)*Dld_eff(i,j,Nz) + alpha2_d(i,j,Nz)*Dgd_eff(i,j,Nz))/dz*2.*ne(i,j,Nz)/dz)
                    q(i,j,Nz) = ne(i,j,Nz)*td_o(i,j,Nz) + &
                            dt*theta(i,j,Nz)* &
                            (qd_0*b_o(i,j,Nz)*lc_o(i,j,Nz)/(kc+lc_o(i,j,Nz))*lo_o(i,j,Nz)/(ko+lo_o(i,j,Nz))) - &
                            dt * (alpha1_d(i,j,Nz)*min(w(i,j,Nz+1),0.) - &
                            2.*(alpha1_d(i,j,Nz)*Dld_eff(i,j,Nz)+alpha2_d(i,j,Nz)*Dgd_eff(i,j,Nz))/dz)*ne(i,j,Nz)/dz*td_t(i,j)
                 endif
                 if(tdbc_t.eq.2) then
                    aa(i,j,Nz) = ne(i,j,Nz) + dt*(-(alpha1_d(i-1,j,Nz)+alpha1_d(i,j,Nz))/2.*min(u(i,j,Nz),0.)* &
                         (ne(i-1,j,Nz)+ne(i,j,Nz))/2./dx + &
                         ((alpha1_d(i-1,j,Nz)+alpha1_d(i,j,Nz))/2.*(Dld_eff(i-1,j,Nz)+Dld_eff(i,j,Nz))/2. + &
                         (alpha2_d(i-1,j,Nz)+alpha2_d(i,j,Nz))/2.*(Dgd_eff(i-1,j,Nz)+Dgd_eff(i,j,Nz))/2.)/dx* &
                         (ne(i-1,j,Nz)+ne(i,j,Nz))/2./dx + &
                         (alpha1_d(i+1,j,Nz)+alpha1_d(i,j,Nz))/2.*max(u(i+1,j,Nz),0.)* &
                         (ne(i+1,j,Nz)+ne(i,j,Nz))/2./dx + &
                         ((alpha1_d(i+1,j,Nz)+alpha1_d(i,j,Nz))/2.*(Dld_eff(i+1,j,Nz)+Dld_eff(i,j,Nz))/2. + &
                         (alpha2_d(i+1,j,Nz)+alpha2_d(i,j,Nz))/2.*(Dgd_eff(i+1,j,Nz)+Dgd_eff(i,j,Nz))/2.)/dx* &
                         (ne(i+1,j,Nz)+ne(i,j,Nz))/2./dx - &
                         (alpha1_d(i,j-1,Nz)+alpha1_d(i,j,Nz))/2.*min(v(i,j,Nz),0.)*(ne(i,j-1,Nz)+ne(i,j,Nz))/2./dy + &
                         ((alpha1_d(i,j-1,Nz)+alpha1_d(i,j,Nz))/2.*(Dld_eff(i,j-1,Nz)+Dld_eff(i,j,Nz))/2. + &
                         (alpha2_d(i,j-1,Nz)+alpha2_d(i,j,Nz))/2.*(Dgd_eff(i,j-1,Nz)+Dgd_eff(i,j,Nz))/2.)/dy* &
                         (ne(i,j-1,Nz)+ne(i,j,Nz))/2./dy + &
                         (alpha1_d(i,j+1,Nz)+alpha1_d(i,j,Nz))/2.*max(v(i,j+1,Nz),0.)*(ne(i,j+1,Nz)+ne(i,j,Nz))/2./dy + &
                         ((alpha1_d(i,j+1,Nz)+alpha1_d(i,j,Nz))/2.*(Dld_eff(i,j+1,Nz)+Dld_eff(i,j,Nz))/2. + &
                         (alpha2_d(i,j+1,Nz)+alpha2_d(i,j,Nz))/2.*(Dgd_eff(i,j+1,Nz)+Dgd_eff(i,j,Nz))/2.)/dy* &
                         (ne(i,j+1,Nz)+ne(i,j,Nz))/2./dy - &
                         (alpha1_d(i,j,Nz-1)+alpha1_d(i,j,Nz))/2.*min(w(i,j,Nz),0.)*(ne(i,j,Nz-1)+ne(i,j,Nz))/2./dz + &
                         ((alpha1_d(i,j,Nz-1)+alpha1_d(i,j,Nz))/2.*(Dld_eff(i,j,Nz-1)+Dld_eff(i,j,Nz))/2. + &
                         (alpha2_d(i,j,Nz-1)+alpha2_d(i,j,Nz))/2.*(Dgd_eff(i,j,Nz-1)+Dgd_eff(i,j,Nz))/2.)/dz* &
                         (ne(i,j,Nz-1)+ne(i,j,Nz))/2./dz + &
                         alpha1_d(i,j,Nz)*max(w(i,j,Nz+1),0.)*ne(i,j,Nz)/dz)
                    q(i,j,Nz) = ne(i,j,Nz)*td_o(i,j,Nz) + &
                            dt*theta(i,j,Nz)* &
                            (qd_0*b_o(i,j,Nz)*lc_o(i,j,Nz)/(kc+lc_o(i,j,Nz))*lo_o(i,j,Nz)/(ko+lo_o(i,j,Nz)))
                    endif
                 endif
              enddo
           enddo



           !! for inside

           Do i=2,Nx-1
              Do j=2,Ny-1
                 Do k=2,Nz-1
                    if(wall(i,j,k)) then
                       aa(i,j,k) = 0.
                       bb(i,j,k) = 0.
                       cc(i,j,k) = 0.
                       dd(i,j,k) = 0.
                       ee(i,j,k) = 0.
                       ff(i,j,k) = 0.
                       gg(i,j,k) = 0.
                       q(i,j,k) = 0.
                    else
                    aa(i,j,k) = ne(i,j,k) + dt*(-(alpha1_d(i-1,j,k)+alpha1_d(i,j,k))/2.*min(u(i,j,k),0.)* &
                         (ne(i-1,j,k)+ne(i,j,k))/2./dx + &
                         ((alpha1_d(i-1,j,k)+alpha1_d(i,j,k))/2.*(Dld_eff(i-1,j,k)+Dld_eff(i,j,k))/2. + &
                         (alpha2_d(i-1,j,k)+alpha2_d(i,j,k))/2.*(Dgd_eff(i-1,j,k)+Dgd_eff(i,j,k))/2.)/dx* &
                         (ne(i-1,j,k)+ne(i,j,k))/2./dx + &
                         (alpha1_d(i+1,j,k)+alpha1_d(i,j,k))/2.*max(u(i+1,j,k),0.)* &
                         (ne(i+1,j,k)+ne(i,j,k))/2./dx + &
                         ((alpha1_d(i+1,j,k)+alpha1_d(i,j,k))/2.*(Dld_eff(i+1,j,k)+Dld_eff(i,j,k))/2. + &
                         (alpha2_d(i+1,j,k)+alpha2_d(i,j,k))/2.*(Dgd_eff(i+1,j,k)+Dgd_eff(i,j,k))/2.)/dx* &
                         (ne(i+1,j,k)+ne(i,j,k))/2./dx - &
                         (alpha1_d(i,j-1,k)+alpha1_d(i,j,k))/2.*min(v(i,j,k),0.)*(ne(i,j-1,k)+ne(i,j,k))/2./dy + &
                         ((alpha1_d(i,j-1,k)+alpha1_d(i,j,k))/2.*(Dld_eff(i,j-1,k)+Dld_eff(i,j,k))/2. + &
                         (alpha2_d(i,j-1,k)+alpha2_d(i,j,k))/2.*(Dgd_eff(i,j-1,k)+Dgd_eff(i,j,k))/2.)/dy* &
                         (ne(i,j-1,k)+ne(i,j,k))/2./dy + &
                         (alpha1_d(i,j+1,k)+alpha1_d(i,j,k))/2.*max(v(i,j+1,k),0.)*(ne(i,j+1,k)+ne(i,j,k))/2./dy + &
                         ((alpha1_d(i,j+1,k)+alpha1_d(i,j,k))/2.*(Dld_eff(i,j+1,k)+Dld_eff(i,j,k))/2. + &
                         (alpha2_d(i,j+1,k)+alpha2_d(i,j,k))/2.*(Dgd_eff(i,j+1,k)+Dgd_eff(i,j,k))/2.)/dy* &
                         (ne(i,j+1,k)+ne(i,j,k))/2./dy - &
                         (alpha1_d(i,j,k-1)+alpha1_d(i,j,k))/2.*min(w(i,j,k),0.)*(ne(i,j,k-1)+ne(i,j,k))/2./dz + &
                         ((alpha1_d(i,j,k-1)+alpha1_d(i,j,k))/2.*(Dld_eff(i,j,k-1)+Dld_eff(i,j,k))/2. + &
                         (alpha2_d(i,j,k-1)+alpha2_d(i,j,k))/2.*(Dgd_eff(i,j,k-1)+Dgd_eff(i,j,k))/2.)/dz* &
                         (ne(i,j,k-1)+ne(i,j,k))/2./dz + &
                         (alpha1_d(i,j,k+1)+alpha1_d(i,j,k))/2.*max(w(i,j,k+1),0.)*(ne(i,j,k+1)+ne(i,j,k))/2./dz + &
                         ((alpha1_d(i,j,k+1)+alpha1_d(i,j,k))/2.*(Dld_eff(i,j,k+1)+Dld_eff(i,j,k))/2. + &
                         (alpha2_d(i,j,k+1)+alpha2_d(i,j,k))/2.*(Dgd_eff(i,j,k+1)+Dgd_eff(i,j,k))/2.)/dz* &
                         (ne(i,j,k+1)+ne(i,j,k))/2./dz)
                    bb(i,j,k) = -dt*((alpha1_d(i-1,j,k)+alpha1_d(i,j,k))/2.*max(u(i,j,k),0.) + &
                         ((alpha1_d(i-1,j,k)+alpha1_d(i,j,k))/2.*(Dld_eff(i-1,j,k)+Dld_eff(i,j,k))/2. + &
                         (alpha2_d(i-1,j,k)+alpha2_d(i,j,k))/2.*(Dgd_eff(i-1,j,k)+Dgd_eff(i,j,k))/2.)/dx)* &
                         (ne(i-1,j,k)+ne(i,j,k))/2./dx
                    cc(i,j,k) = dt*((alpha1_d(i+1,j,k)+alpha1_d(i,j,k))/2.*min(u(i+1,j,k),0.) - &
                         ((alpha1_d(i+1,j,k)+alpha1_d(i,j,k))/2.*(Dld_eff(i+1,j,k)+Dld_eff(i,j,k))/2. + &
                         (alpha2_d(i+1,j,k)+alpha2_d(i,j,k))/2.*(Dgd_eff(i+1,j,k)+Dgd_eff(i,j,k))/2.)/dx)* &
                         (ne(i+1,j,k)+ne(i,j,k))/2./dx
                    dd(i,j,k) = -dt*((alpha1_d(i,j-1,k)+alpha1_d(i,j,k))/2.*max(v(i,j,k),0.) + &
                         ((alpha1_d(i,j-1,k)+alpha1_d(i,j,k))/2.*(Dld_eff(i,j-1,k)+Dld_eff(i,j,k))/2. + &
                         (alpha2_d(i,j-1,k)+alpha2_d(i,j,k))/2.*(Dgd_eff(i,j-1,k)+Dgd_eff(i,j,k))/2.)/dy)* &
                         (ne(i,j-1,k)+ne(i,j,k))/2./dy
                    ee(i,j,k) = dt*((alpha1_d(i,j+1,k)+alpha1_d(i,j,k))/2.*min(v(i,j+1,k),0.) - &
                         ((alpha1_d(i,j+1,k)+alpha1_d(i,j,k))/2.*(Dld_eff(i,j+1,k)+Dld_eff(i,j,k))/2. + &
                         (alpha2_d(i,j+1,k)+alpha2_d(i,j,k))/2.*(Dgd_eff(i,j+1,k)+Dgd_eff(i,j,k))/2.)/dy)* &
                         (ne(i,j+1,k)+ne(i,j,k))/2./dy
                    ff(i,j,k) = -dt*((alpha1_d(i,j,k-1)+alpha1_d(i,j,k))/2.*max(w(i,j,k),0.) + &
                         ((alpha1_d(i,j,k-1)+alpha1_d(i,j,k))/2.*(Dld_eff(i,j,k-1)+Dld_eff(i,j,k))/2. + &
                         (alpha2_d(i,j,k-1)+alpha2_d(i,j,k))/2.*(Dgd_eff(i,j,k-1)+Dgd_eff(i,j,k))/2.)/dz)* &
                         (ne(i,j,k-1)+ne(i,j,k))/2./dz
                    gg(i,j,k) = dt*((alpha1_d(i,j,k+1)+alpha1_d(i,j,k))/2.*min(w(i,j,k+1),0.) - &
                         ((alpha1_d(i,j,k+1)+alpha1_d(i,j,k))/2.*(Dld_eff(i,j,k+1)+Dld_eff(i,j,k))/2. + &
                         (alpha2_d(i,j,k+1)+alpha2_d(i,j,k))/2.*(Dgd_eff(i,j,k+1)+Dgd_eff(i,j,k))/2.)/dz)* &
                         (ne(i,j,k+1)+ne(i,j,k))/2./dz

                    q(i,j,k) = ne(i,j,k)*td_o(i,j,k) + &
                         dt*theta(i,j,k)* &
                         (qd_0*b_o(i,j,k)*lc_o(i,j,k)/(Kc+lc_o(i,j,k))*lo_o(i,j,k)/(Ko+lo_o(i,j,k)))
                    endif
                 Enddo
              Enddo
           Enddo



        else 

           Do i=2,Nx-1
              Do j=2,Ny-1
                    if(wall(i,j,1)) then
                       aa(i,j,1) = 0.
                       bb(i,j,1) = 0.
                       cc(i,j,1) = 0.
                       dd(i,j,1) = 0.
                       ee(i,j,1) = 0.
                       ff(i,j,1) = 0.
                       gg(i,j,1) = 0.
                       q(i,j,1) = 0.
                    else
                       aa(i,j,1) = ne(i,j,1) + dt*(-(alpha1_d(i-1,j,1)+alpha1_d(i,j,1))/2.*min(u(i,j,1),0.)* &
                            (ne(i-1,j,1)+ne(i,j,1))/2./dx + &
                            ((alpha1_d(i-1,j,1)+alpha1_d(i,j,1))/2.*(Dld_eff(i-1,j,1)+Dld_eff(i,j,1))/2. + &
                            (alpha2_d(i-1,j,1)+alpha2_d(i,j,1))/2.*(Dgd_eff(i-1,j,1)+Dgd_eff(i,j,1))/2.)/dx* &
                            (ne(i-1,j,1)+ne(i,j,1))/2./dx + &
                            (alpha1_d(i+1,j,1)+alpha1_d(i,j,1))/2.*max(u(i+1,j,1),0.)* &
                            (ne(i+1,j,1)+ne(i,j,1))/2./dx + &
                            ((alpha1_d(i+1,j,1)+alpha1_d(i,j,1))/2.*(Dld_eff(i+1,j,1)+Dld_eff(i,j,1))/2. + &
                            (alpha2_d(i+1,j,1)+alpha2_d(i,j,1))/2.*(Dgd_eff(i+1,j,1)+Dgd_eff(i,j,1))/2.)/dx* &
                            (ne(i+1,j,1)+ne(i,j,1))/2./dx - &
                            (alpha1_d(i,j-1,1)+alpha1_d(i,j,1))/2.*min(v(i,j,1),0.)*(ne(i,j-1,1)+ne(i,j,1))/2./dy + &
                            ((alpha1_d(i,j-1,1)+alpha1_d(i,j,1))/2.*(Dld_eff(i,j-1,1)+Dld_eff(i,j,1))/2. + &
                            (alpha2_d(i,j-1,1)+alpha2_d(i,j,1))/2.*(Dgd_eff(i,j-1,1)+Dgd_eff(i,j,1))/2.)/dy* &
                            (ne(i,j-1,1)+ne(i,j,1))/2./dy + &
                            (alpha1_d(i,j+1,1)+alpha1_d(i,j,1))/2.*max(v(i,j+1,1),0.)*(ne(i,j+1,1)+ne(i,j,1))/2./dy + &
                            ((alpha1_d(i,j+1,1)+alpha1_d(i,j,1))/2.*(Dld_eff(i,j+1,1)+Dld_eff(i,j,1))/2. + &
                            (alpha2_d(i,j+1,1)+alpha2_d(i,j,1))/2.*(Dgd_eff(i,j+1,1)+Dgd_eff(i,j,1))/2.)/dy* &
                            (ne(i,j+1,1)+ne(i,j,1))/2./dy - &
                            (tmp_alpha1_d(i,j,1)+alpha1_d(i,j,1))/2.*min(w(i,j,1),0.)*(tmp_ne(i,j,1)+ne(i,j,1))/2./dz + &
                            ((tmp_alpha1_d(i,j,1)+alpha1_d(i,j,1))/2.*(tmp_Dld_eff(i,j,1)+Dld_eff(i,j,1))/2. + &
                            (tmp_alpha2_d(i,j,1)+alpha2_d(i,j,1))/2.*(tmp_Dgd_eff(i,j,1)+Dgd_eff(i,j,1))/2.)/dz* &
                            (tmp_ne(i,j,1)+ne(i,j,1))/2./dz + &
                            (alpha1_d(i,j,1+1)+alpha1_d(i,j,1))/2.*max(w(i,j,1+1),0.)*(ne(i,j,1+1)+ne(i,j,1))/2./dz + &
                            ((alpha1_d(i,j,1+1)+alpha1_d(i,j,1))/2.*(Dld_eff(i,j,1+1)+Dld_eff(i,j,1))/2. + &
                            (alpha2_d(i,j,1+1)+alpha2_d(i,j,1))/2.*(Dgd_eff(i,j,1+1)+Dgd_eff(i,j,1))/2.)/dz* &
                            (ne(i,j,1+1)+ne(i,j,1))/2./dz)
                       bb(i,j,1) = -dt*((alpha1_d(i-1,j,1)+alpha1_d(i,j,1))/2.*max(u(i,j,1),0.) + &
                            ((alpha1_d(i-1,j,1)+alpha1_d(i,j,1))/2.*(Dld_eff(i-1,j,1)+Dld_eff(i,j,1))/2. + &
                            (alpha2_d(i-1,j,1)+alpha2_d(i,j,1))/2.*(Dgd_eff(i-1,j,1)+Dgd_eff(i,j,1))/2.)/dx)* &
                            (ne(i-1,j,1)+ne(i,j,1))/2./dx
                       cc(i,j,1) = dt*((alpha1_d(i+1,j,1)+alpha1_d(i,j,1))/2.*min(u(i+1,j,1),0.) - &
                            ((alpha1_d(i+1,j,1)+alpha1_d(i,j,1))/2.*(Dld_eff(i+1,j,1)+Dld_eff(i,j,1))/2. + &
                            (alpha2_d(i+1,j,1)+alpha2_d(i,j,1))/2.*(Dgd_eff(i+1,j,1)+Dgd_eff(i,j,1))/2.)/dx)* &
                            (ne(i+1,j,1)+ne(i,j,1))/2./dx
                       dd(i,j,1) = -dt*((alpha1_d(i,j-1,1)+alpha1_d(i,j,1))/2.*max(v(i,j,1),0.) + &
                            ((alpha1_d(i,j-1,1)+alpha1_d(i,j,1))/2.*(Dld_eff(i,j-1,1)+Dld_eff(i,j,1))/2. + &
                            (alpha2_d(i,j-1,1)+alpha2_d(i,j,1))/2.*(Dgd_eff(i,j-1,1)+Dgd_eff(i,j,1))/2.)/dy)* &
                            (ne(i,j-1,1)+ne(i,j,1))/2./dy
                       ee(i,j,1) = dt*((alpha1_d(i,j+1,1)+alpha1_d(i,j,1))/2.*min(v(i,j+1,1),0.) - &
                            ((alpha1_d(i,j+1,1)+alpha1_d(i,j,1))/2.*(Dld_eff(i,j+1,1)+Dld_eff(i,j,1))/2. + &
                            (alpha2_d(i,j+1,1)+alpha2_d(i,j,1))/2.*(Dgd_eff(i,j+1,1)+Dgd_eff(i,j,1))/2.)/dy)* &
                            (ne(i,j+1,1)+ne(i,j,1))/2./dy
                       ff(i,j,1) = 0.
                       gg(i,j,1) = dt*((alpha1_d(i,j,1+1)+alpha1_d(i,j,1))/2.*min(w(i,j,1+1),0.) - &
                            ((alpha1_d(i,j,1+1)+alpha1_d(i,j,1))/2.*(Dld_eff(i,j,1+1)+Dld_eff(i,j,1))/2. + &
                            (alpha2_d(i,j,1+1)+alpha2_d(i,j,1))/2.*(Dgd_eff(i,j,1+1)+Dgd_eff(i,j,1))/2.)/dz)* &
                            (ne(i,j,1+1)+ne(i,j,1))/2./dz

                       q(i,j,1) = ne(i,j,1)*td_o(i,j,1) + &
                            dt*theta(i,j,1)* &
                            (qd_0*b_o(i,j,1)*lc_o(i,j,1)/(Kc+lc_o(i,j,1))*lo_o(i,j,1)/(Ko+lo_o(i,j,1))) - &
                            (-dt*((tmp_alpha1_d(i,j,1)+alpha1_d(i,j,1))/2.*max(w(i,j,1),0.) + &
                             ((tmp_alpha1_d(i,j,1)+alpha1_d(i,j,1))/2.*(tmp_Dld_eff(i,j,1)+Dld_eff(i,j,1))/2. + &
                            (tmp_alpha2_d(i,j,1)+alpha2_d(i,j,1))/2.*(tmp_Dgd_eff(i,j,1)+Dgd_eff(i,j,1))/2.)/dz)* &
                            (tmp_ne(i,j,1)+ne(i,j,1))/2./dz)*tmp_td_o(i,j,1)
                    endif
              Enddo
           Enddo


           Do i=2,Nx-1
              Do j=2,Ny-1
                    if(wall(i,j,Nz)) then
                       aa(i,j,Nz) = 0.
                       bb(i,j,Nz) = 0.
                       cc(i,j,Nz) = 0.
                       dd(i,j,Nz) = 0.
                       ee(i,j,Nz) = 0.
                       ff(i,j,Nz) = 0.
                       gg(i,j,Nz) = 0.
                       q(i,j,Nz) = 0.
                    else
                       aa(i,j,Nz) = ne(i,j,Nz) + dt*(-(alpha1_d(i-1,j,Nz)+alpha1_d(i,j,Nz))/2.*min(u(i,j,Nz),0.)* &
                            (ne(i-1,j,Nz)+ne(i,j,Nz))/2./dx + &
                            ((alpha1_d(i-1,j,Nz)+alpha1_d(i,j,Nz))/2.*(Dld_eff(i-1,j,Nz)+Dld_eff(i,j,Nz))/2. + &
                            (alpha2_d(i-1,j,Nz)+alpha2_d(i,j,Nz))/2.*(Dgd_eff(i-1,j,Nz)+Dgd_eff(i,j,Nz))/2.)/dx* &
                            (ne(i-1,j,Nz)+ne(i,j,Nz))/2./dx + &
                            (alpha1_d(i+1,j,Nz)+alpha1_d(i,j,Nz))/2.*max(u(i+1,j,Nz),0.)* &
                            (ne(i+1,j,Nz)+ne(i,j,Nz))/2./dx + &
                            ((alpha1_d(i+1,j,Nz)+alpha1_d(i,j,Nz))/2.*(Dld_eff(i+1,j,Nz)+Dld_eff(i,j,Nz))/2. + &
                            (alpha2_d(i+1,j,Nz)+alpha2_d(i,j,Nz))/2.*(Dgd_eff(i+1,j,Nz)+Dgd_eff(i,j,Nz))/2.)/dx* &
                            (ne(i+1,j,Nz)+ne(i,j,Nz))/2./dx - &
                            (alpha1_d(i,j-1,Nz)+alpha1_d(i,j,Nz))/2.*min(v(i,j,Nz),0.)*(ne(i,j-1,Nz)+ne(i,j,Nz))/2./dy + &
                            ((alpha1_d(i,j-1,Nz)+alpha1_d(i,j,Nz))/2.*(Dld_eff(i,j-1,Nz)+Dld_eff(i,j,Nz))/2. + &
                            (alpha2_d(i,j-1,Nz)+alpha2_d(i,j,Nz))/2.*(Dgd_eff(i,j-1,Nz)+Dgd_eff(i,j,Nz))/2.)/dy* &
                            (ne(i,j-1,Nz)+ne(i,j,Nz))/2./dy + &
                            (alpha1_d(i,j+1,Nz)+alpha1_d(i,j,Nz))/2.*max(v(i,j+1,Nz),0.)*(ne(i,j+1,Nz)+ne(i,j,Nz))/2./dy + &
                            ((alpha1_d(i,j+1,Nz)+alpha1_d(i,j,Nz))/2.*(Dld_eff(i,j+1,Nz)+Dld_eff(i,j,Nz))/2. + &
                            (alpha2_d(i,j+1,Nz)+alpha2_d(i,j,Nz))/2.*(Dgd_eff(i,j+1,Nz)+Dgd_eff(i,j,Nz))/2.)/dy* &
                            (ne(i,j+1,Nz)+ne(i,j,Nz))/2./dy - &
                            (alpha1_d(i,j,Nz-1)+alpha1_d(i,j,Nz))/2.*min(w(i,j,Nz),0.)*(ne(i,j,Nz-1)+ne(i,j,Nz))/2./dz + &
                            ((alpha1_d(i,j,Nz-1)+alpha1_d(i,j,Nz))/2.*(Dld_eff(i,j,Nz-1)+Dld_eff(i,j,Nz))/2. + &
                            (alpha2_d(i,j,Nz-1)+alpha2_d(i,j,Nz))/2.*(Dgd_eff(i,j,Nz-1)+Dgd_eff(i,j,Nz))/2.)/dz* &
                            (ne(i,j,Nz-1)+ne(i,j,Nz))/2./dz + &
                            (tmp_alpha1_d(i,j,2)+alpha1_d(i,j,Nz))/2.*max(w(i,j,Nz+1),0.)*(tmp_ne(i,j,2)+ne(i,j,Nz))/2./dz + &
                            ((tmp_alpha1_d(i,j,2)+alpha1_d(i,j,Nz))/2.*(tmp_Dld_eff(i,j,2)+Dld_eff(i,j,Nz))/2. + &
                            (tmp_alpha2_d(i,j,2)+alpha2_d(i,j,Nz))/2.*(tmp_Dgd_eff(i,j,2)+Dgd_eff(i,j,Nz))/2.)/dz* &
                            (tmp_ne(i,j,2)+ne(i,j,Nz))/2./dz)
                       bb(i,j,Nz) = -dt*((alpha1_d(i-1,j,Nz)+alpha1_d(i,j,Nz))/2.*max(u(i,j,Nz),0.) + &
                            ((alpha1_d(i-1,j,Nz)+alpha1_d(i,j,Nz))/2.*(Dld_eff(i-1,j,Nz)+Dld_eff(i,j,Nz))/2. + &
                            (alpha2_d(i-1,j,Nz)+alpha2_d(i,j,Nz))/2.*(Dgd_eff(i-1,j,Nz)+Dgd_eff(i,j,Nz))/2.)/dx)* &
                            (ne(i-1,j,Nz)+ne(i,j,Nz))/2./dx
                       cc(i,j,Nz) = dt*((alpha1_d(i+1,j,Nz)+alpha1_d(i,j,Nz))/2.*min(u(i+1,j,Nz),0.) - &
                            ((alpha1_d(i+1,j,Nz)+alpha1_d(i,j,Nz))/2.*(Dld_eff(i+1,j,Nz)+Dld_eff(i,j,Nz))/2. + &
                            (alpha2_d(i+1,j,Nz)+alpha2_d(i,j,Nz))/2.*(Dgd_eff(i+1,j,Nz)+Dgd_eff(i,j,Nz))/2.)/dx)* &
                            (ne(i+1,j,Nz)+ne(i,j,Nz))/2./dx
                       dd(i,j,Nz) = -dt*((alpha1_d(i,j-1,Nz)+alpha1_d(i,j,Nz))/2.*max(v(i,j,Nz),0.) + &
                            ((alpha1_d(i,j-1,Nz)+alpha1_d(i,j,Nz))/2.*(Dld_eff(i,j-1,Nz)+Dld_eff(i,j,Nz))/2. + &
                            (alpha2_d(i,j-1,Nz)+alpha2_d(i,j,Nz))/2.*(Dgd_eff(i,j-1,Nz)+Dgd_eff(i,j,Nz))/2.)/dy)* &
                            (ne(i,j-1,Nz)+ne(i,j,Nz))/2./dy
                       ee(i,j,Nz) = dt*((alpha1_d(i,j+1,Nz)+alpha1_d(i,j,Nz))/2.*min(v(i,j+1,Nz),0.) - &
                            ((alpha1_d(i,j+1,Nz)+alpha1_d(i,j,Nz))/2.*(Dld_eff(i,j+1,Nz)+Dld_eff(i,j,Nz))/2. + &
                            (alpha2_d(i,j+1,Nz)+alpha2_d(i,j,Nz))/2.*(Dgd_eff(i,j+1,Nz)+Dgd_eff(i,j,Nz))/2.)/dy)* &
                            (ne(i,j+1,Nz)+ne(i,j,Nz))/2./dy
                       ff(i,j,Nz) = -dt*((alpha1_d(i,j,Nz-1)+alpha1_d(i,j,Nz))/2.*max(w(i,j,Nz),0.) + &
                             ((alpha1_d(i,j,Nz-1)+alpha1_d(i,j,Nz))/2.*(Dld_eff(i,j,Nz-1)+Dld_eff(i,j,Nz))/2. + &
                            (alpha2_d(i,j,Nz-1)+alpha2_d(i,j,Nz))/2.*(Dgd_eff(i,j,Nz-1)+Dgd_eff(i,j,Nz))/2.)/dz)* &
                            (ne(i,j,Nz-1)+ne(i,j,Nz))/2./dz
                       gg(i,j,Nz) = 0

                       q(i,j,Nz) = ne(i,j,Nz)*td_o(i,j,Nz) + &
                            dt*theta(i,j,Nz)* &
                            (qd_0*b_o(i,j,Nz)*lc_o(i,j,Nz)/(Kc+lc_o(i,j,Nz))*lo_o(i,j,Nz)/(Ko+lo_o(i,j,Nz))) - &
                            (dt*((tmp_alpha1_d(i,j,2)+alpha1_d(i,j,Nz))/2.*min(w(i,j,Nz+1),0.) - &
                            ((tmp_alpha1_d(i,j,2)+alpha1_d(i,j,Nz))/2.*(tmp_Dld_eff(i,j,2)+Dld_eff(i,j,Nz))/2. + &
                            (tmp_alpha2_d(i,j,2)+alpha2_d(i,j,Nz))/2.*(tmp_Dgd_eff(i,j,2)+Dgd_eff(i,j,Nz))/2.)/dz)* &
                            (tmp_ne(i,j,2)+ne(i,j,Nz))/2./dz)*tmp_td_o(i,j,2)
                    endif
              Enddo
           Enddo


           !! for inside

           Do i=2,Nx-1
              Do j=2,Ny-1
                 Do k=2,Nz-1
                    if(wall(i,j,k)) then
                       aa(i,j,k) = 0.
                       bb(i,j,k) = 0.
                       cc(i,j,k) = 0.
                       dd(i,j,k) = 0.
                       ee(i,j,k) = 0.
                       ff(i,j,k) = 0.
                       gg(i,j,k) = 0.
                       q(i,j,k) = 0.
                    else
                       aa(i,j,k) = ne(i,j,k) + dt*(-(alpha1_d(i-1,j,k)+alpha1_d(i,j,k))/2.*min(u(i,j,k),0.)* &
                            (ne(i-1,j,k)+ne(i,j,k))/2./dx + &
                            ((alpha1_d(i-1,j,k)+alpha1_d(i,j,k))/2.*(Dld_eff(i-1,j,k)+Dld_eff(i,j,k))/2. + &
                            (alpha2_d(i-1,j,k)+alpha2_d(i,j,k))/2.*(Dgd_eff(i-1,j,k)+Dgd_eff(i,j,k))/2.)/dx* &
                            (ne(i-1,j,k)+ne(i,j,k))/2./dx + &
                            (alpha1_d(i+1,j,k)+alpha1_d(i,j,k))/2.*max(u(i+1,j,k),0.)* &
                            (ne(i+1,j,k)+ne(i,j,k))/2./dx + &
                            ((alpha1_d(i+1,j,k)+alpha1_d(i,j,k))/2.*(Dld_eff(i+1,j,k)+Dld_eff(i,j,k))/2. + &
                            (alpha2_d(i+1,j,k)+alpha2_d(i,j,k))/2.*(Dgd_eff(i+1,j,k)+Dgd_eff(i,j,k))/2.)/dx* &
                            (ne(i+1,j,k)+ne(i,j,k))/2./dx - &
                            (alpha1_d(i,j-1,k)+alpha1_d(i,j,k))/2.*min(v(i,j,k),0.)*(ne(i,j-1,k)+ne(i,j,k))/2./dy + &
                            ((alpha1_d(i,j-1,k)+alpha1_d(i,j,k))/2.*(Dld_eff(i,j-1,k)+Dld_eff(i,j,k))/2. + &
                            (alpha2_d(i,j-1,k)+alpha2_d(i,j,k))/2.*(Dgd_eff(i,j-1,k)+Dgd_eff(i,j,k))/2.)/dy* &
                            (ne(i,j-1,k)+ne(i,j,k))/2./dy + &
                            (alpha1_d(i,j+1,k)+alpha1_d(i,j,k))/2.*max(v(i,j+1,k),0.)*(ne(i,j+1,k)+ne(i,j,k))/2./dy + &
                            ((alpha1_d(i,j+1,k)+alpha1_d(i,j,k))/2.*(Dld_eff(i,j+1,k)+Dld_eff(i,j,k))/2. + &
                            (alpha2_d(i,j+1,k)+alpha2_d(i,j,k))/2.*(Dgd_eff(i,j+1,k)+Dgd_eff(i,j,k))/2.)/dy* &
                            (ne(i,j+1,k)+ne(i,j,k))/2./dy - &
                            (alpha1_d(i,j,k-1)+alpha1_d(i,j,k))/2.*min(w(i,j,k),0.)*(ne(i,j,k-1)+ne(i,j,k))/2./dz + &
                            ((alpha1_d(i,j,k-1)+alpha1_d(i,j,k))/2.*(Dld_eff(i,j,k-1)+Dld_eff(i,j,k))/2. + &
                            (alpha2_d(i,j,k-1)+alpha2_d(i,j,k))/2.*(Dgd_eff(i,j,k-1)+Dgd_eff(i,j,k))/2.)/dz* &
                            (ne(i,j,k-1)+ne(i,j,k))/2./dz + &
                            (alpha1_d(i,j,k+1)+alpha1_d(i,j,k))/2.*max(w(i,j,k+1),0.)*(ne(i,j,k+1)+ne(i,j,k))/2./dz + &
                            ((alpha1_d(i,j,k+1)+alpha1_d(i,j,k))/2.*(Dld_eff(i,j,k+1)+Dld_eff(i,j,k))/2. + &
                            (alpha2_d(i,j,k+1)+alpha2_d(i,j,k))/2.*(Dgd_eff(i,j,k+1)+Dgd_eff(i,j,k))/2.)/dz* &
                            (ne(i,j,k+1)+ne(i,j,k))/2./dz)
                       bb(i,j,k) = -dt*((alpha1_d(i-1,j,k)+alpha1_d(i,j,k))/2.*max(u(i,j,k),0.) + &
                            ((alpha1_d(i-1,j,k)+alpha1_d(i,j,k))/2.*(Dld_eff(i-1,j,k)+Dld_eff(i,j,k))/2. + &
                            (alpha2_d(i-1,j,k)+alpha2_d(i,j,k))/2.*(Dgd_eff(i-1,j,k)+Dgd_eff(i,j,k))/2.)/dx)* &
                            (ne(i-1,j,k)+ne(i,j,k))/2./dx
                       cc(i,j,k) = dt*((alpha1_d(i+1,j,k)+alpha1_d(i,j,k))/2.*min(u(i+1,j,k),0.) - &
                            ((alpha1_d(i+1,j,k)+alpha1_d(i,j,k))/2.*(Dld_eff(i+1,j,k)+Dld_eff(i,j,k))/2. + &
                            (alpha2_d(i+1,j,k)+alpha2_d(i,j,k))/2.*(Dgd_eff(i+1,j,k)+Dgd_eff(i,j,k))/2.)/dx)* &
                            (ne(i+1,j,k)+ne(i,j,k))/2./dx
                       dd(i,j,k) = -dt*((alpha1_d(i,j-1,k)+alpha1_d(i,j,k))/2.*max(v(i,j,k),0.) + &
                            ((alpha1_d(i,j-1,k)+alpha1_d(i,j,k))/2.*(Dld_eff(i,j-1,k)+Dld_eff(i,j,k))/2. + &
                            (alpha2_d(i,j-1,k)+alpha2_d(i,j,k))/2.*(Dgd_eff(i,j-1,k)+Dgd_eff(i,j,k))/2.)/dy)* &
                            (ne(i,j-1,k)+ne(i,j,k))/2./dy
                       ee(i,j,k) = dt*((alpha1_d(i,j+1,k)+alpha1_d(i,j,k))/2.*min(v(i,j+1,k),0.) - &
                            ((alpha1_d(i,j+1,k)+alpha1_d(i,j,k))/2.*(Dld_eff(i,j+1,k)+Dld_eff(i,j,k))/2. + &
                            (alpha2_d(i,j+1,k)+alpha2_d(i,j,k))/2.*(Dgd_eff(i,j+1,k)+Dgd_eff(i,j,k))/2.)/dy)* &
                            (ne(i,j+1,k)+ne(i,j,k))/2./dy
                       ff(i,j,k) = -dt*((alpha1_d(i,j,k-1)+alpha1_d(i,j,k))/2.*max(w(i,j,k),0.) + &
                             ((alpha1_d(i,j,k-1)+alpha1_d(i,j,k))/2.*(Dld_eff(i,j,k-1)+Dld_eff(i,j,k))/2. + &
                            (alpha2_d(i,j,k-1)+alpha2_d(i,j,k))/2.*(Dgd_eff(i,j,k-1)+Dgd_eff(i,j,k))/2.)/dz)* &
                            (ne(i,j,k-1)+ne(i,j,k))/2./dz
                       gg(i,j,k) = dt*((alpha1_d(i,j,k+1)+alpha1_d(i,j,k))/2.*min(w(i,j,k+1),0.) - &
                            ((alpha1_d(i,j,k+1)+alpha1_d(i,j,k))/2.*(Dld_eff(i,j,k+1)+Dld_eff(i,j,k))/2. + &
                            (alpha2_d(i,j,k+1)+alpha2_d(i,j,k))/2.*(Dgd_eff(i,j,k+1)+Dgd_eff(i,j,k))/2.)/dz)* &
                            (ne(i,j,k+1)+ne(i,j,k))/2./dz

                       q(i,j,k) = ne(i,j,k)*td_o(i,j,k) + &
                            dt*theta(i,j,k)* &
                            (qd_0*b_o(i,j,k)*lc_o(i,j,k)/(Kc+lc_o(i,j,k))*lo_o(i,j,k)/(Ko+lo_o(i,j,k)))
                    endif
                 Enddo
              Enddo
           Enddo


        endif
     endif







     !! cabulate seven-diagdnal elements near wall


     Do i=2,Nx-1
        Do j=2,Ny-1
           Do k=1,Nz
              if(wall(i-1,j,k)) then
                 aa(i,j,k) = aa(i,j,k) - dt*((alpha1_d(i-1,j,k)+alpha1_d(i,j,k))/2.*(Dld_eff(i-1,j,k)+Dld_eff(i,j,k))/2. + &
                            (alpha2_d(i-1,j,k)+alpha2_d(i,j,k))/2.*(Dgd_eff(i-1,j,k)+Dgd_eff(i,j,k))/2.)/dx* &
                            (ne(i-1,j,k)+ne(i,j,k))/2./dx
                 bb(i,j,k) = 0.
              endif
              if(wall(i+1,j,k)) then
                 aa(i,j,k) = aa(i,j,k) - dt*((alpha1_d(i+1,j,k)+alpha1_d(i,j,k))/2.*(Dld_eff(i+1,j,k)+Dld_eff(i,j,k))/2. + &
                            (alpha2_d(i+1,j,k)+alpha2_d(i,j,k))/2.*(Dgd_eff(i+1,j,k)+Dgd_eff(i,j,k))/2.)/dx* &
                            (ne(i+1,j,k)+ne(i,j,k))/2./dx
                 cc(i,j,k) = 0.
              endif
              if(wall(i,j-1,k)) then
                 aa(i,j,k) = aa(i,j,k) - dt*((alpha1_d(i,j-1,k)+alpha1_d(i,j,k))/2.*(Dld_eff(i,j-1,k)+Dld_eff(i,j,k))/2. + &
                            (alpha2_d(i,j-1,k)+alpha2_d(i,j,k))/2.*(Dgd_eff(i,j-1,k)+Dgd_eff(i,j,k))/2.)/dy* &
                            (ne(i,j-1,k)+ne(i,j,k))/2./dy
                 dd(i,j,k) = 0.
              endif
              if(wall(i,j+1,k)) then
                 aa(i,j,k) = aa(i,j,k) - dt*((alpha1_d(i,j+1,k)+alpha1_d(i,j,k))/2.*(Dld_eff(i,j+1,k)+Dld_eff(i,j,k))/2. + &
                            (alpha2_d(i,j+1,k)+alpha2_d(i,j,k))/2.*(Dgd_eff(i,j+1,k)+Dgd_eff(i,j,k))/2.)/dy* &
                            (ne(i,j+1,k)+ne(i,j,k))/2./dy
                 ee(i,j,k) = 0.
              endif

              if(npes.eq.1) then
                 if(k.ne.1 .and. wall(i,j,k-1)) then
                    aa(i,j,k) = aa(i,j,k) - dt*((alpha1_d(i,j,k-1)+alpha1_d(i,j,k))/2.*(Dld_eff(i,j,k-1)+Dld_eff(i,j,k))/2. + &
                            (alpha2_d(i,j,k-1)+alpha2_d(i,j,k))/2.*(Dgd_eff(i,j,k-1)+Dgd_eff(i,j,k))/2.)/dz* &
                            (ne(i,j,k-1)+ne(i,j,k))/2./dz
                    ff(i,j,k) = 0.
                 endif
                 if(k.ne.Nz .and. wall(i,j,k+1)) then
                    aa(i,j,k) = aa(i,j,k) - dt*((alpha1_d(i,j,k+1)+alpha1_d(i,j,k))/2.*(Dld_eff(i,j,k+1)+Dld_eff(i,j,k))/2. + &
                            (alpha2_d(i,j,k+1)+alpha2_d(i,j,k))/2.*(Dgd_eff(i,j,k+1)+Dgd_eff(i,j,k))/2.)/dz* &
                            (ne(i,j,k+1)+ne(i,j,k))/2./dz
                    gg(i,j,k) = 0.
                 endif
              else
                 if(mype.eq.0) then
                    if(k.ne.1 .and. wall(i,j,k-1)) then
                       aa(i,j,k) = aa(i,j,k) - &
                            dt*((alpha1_d(i,j,k-1)+alpha1_d(i,j,k))/2.*(Dld_eff(i,j,k-1)+Dld_eff(i,j,k))/2. + &
                            (alpha2_d(i,j,k-1)+alpha2_d(i,j,k))/2.*(Dgd_eff(i,j,k-1)+Dgd_eff(i,j,k))/2.)/dz* &
                            (ne(i,j,k-1)+ne(i,j,k))/2./dz
                       ff(i,j,k) = 0.
                    endif
                    if(k.eq.Nz) then
                       if(tmp_wall(i,j,2)) then
                          aa(i,j,Nz) = aa(i,j,Nz) - &
                               dt*((tmp_alpha1_d(i,j,2)+alpha1_d(i,j,Nz))/2.*(tmp_Dld_eff(i,j,2)+Dld_eff(i,j,Nz))/2. + &
                               (tmp_alpha2_d(i,j,2)+alpha2_d(i,j,Nz))/2.*(tmp_Dgd_eff(i,j,2)+Dgd_eff(i,j,Nz))/2.)/dz* &
                               (tmp_ne(i,j,2)+ne(i,j,Nz))/2./dz
                          gg(i,j,Nz) = 0.
                       endif
                    else
                       if(wall(i,j,k+1)) then
                          aa(i,j,k) = aa(i,j,k) - &
                               dt*((alpha1_d(i,j,k+1)+alpha1_d(i,j,k))/2.*(Dld_eff(i,j,k+1)+Dld_eff(i,j,k))/2. + &
                               (alpha2_d(i,j,k+1)+alpha2_d(i,j,k))/2.*(Dgd_eff(i,j,k+1)+Dgd_eff(i,j,k))/2.)/dz* &
                               (ne(i,j,k+1)+ne(i,j,k))/2./dz
                          gg(i,j,k) = 0.
                       endif
                    endif
                 else if(mype.eq.npes-1) then
                    if(k.eq.1) then
                       if(tmp_wall(i,j,1)) then
                          aa(i,j,1) = aa(i,j,1) - &
                               dt*((tmp_alpha1_d(i,j,1)+alpha1_d(i,j,1))/2.*(tmp_Dld_eff(i,j,1)+Dld_eff(i,j,1))/2. + &
                               (tmp_alpha2_d(i,j,1)+alpha2_d(i,j,1))/2.*(tmp_Dgd_eff(i,j,1)+Dgd_eff(i,j,1))/2.)/dz* &
                               (tmp_ne(i,j,1)+ne(i,j,1))/2./dz
                          ff(i,j,1) = 0.
                       endif
                    else
                       if(wall(i,j,k-1)) then
                          aa(i,j,k) = aa(i,j,k) - &
                               dt*((alpha1_d(i,j,k-1)+alpha1_d(i,j,k))/2.*(Dld_eff(i,j,k-1)+Dld_eff(i,j,k))/2. + &
                               (alpha2_d(i,j,k-1)+alpha2_d(i,j,k))/2.*(Dgd_eff(i,j,k-1)+Dgd_eff(i,j,k))/2.)/dz* &
                               (ne(i,j,k-1)+ne(i,j,k))/2./dz
                          ff(i,j,k) = 0.
                       endif
                    endif
                    if(k.ne.Nz .and. wall(i,j,k+1)) then
                       aa(i,j,k) = aa(i,j,k) - &
                            dt*((alpha1_d(i,j,k+1)+alpha1_d(i,j,k))/2.*(Dld_eff(i,j,k+1)+Dld_eff(i,j,k))/2. + &
                            (alpha2_d(i,j,k+1)+alpha2_d(i,j,k))/2.*(Dgd_eff(i,j,k+1)+Dgd_eff(i,j,k))/2.)/dz* &
                            (ne(i,j,k+1)+ne(i,j,k))/2./dz
                       gg(i,j,k) = 0.
                    endif
                 else 
                    if(k.eq.1) then
                       if(tmp_wall(i,j,1)) then
                          aa(i,j,1) = aa(i,j,1) - &
                               dt*((tmp_alpha1_d(i,j,1)+alpha1_d(i,j,1))/2.*(tmp_Dld_eff(i,j,1)+Dld_eff(i,j,1))/2. + &
                               (tmp_alpha2_d(i,j,1)+alpha2_d(i,j,1))/2.*(tmp_Dgd_eff(i,j,1)+Dgd_eff(i,j,1))/2.)/dz* &
                               (tmp_ne(i,j,1)+ne(i,j,1))/2./dz
                          ff(i,j,1) = 0.
                       endif
                    else
                       if(wall(i,j,k-1)) then
                          aa(i,j,k) = aa(i,j,k) - &
                               dt*((alpha1_d(i,j,k-1)+alpha1_d(i,j,k))/2.*(Dld_eff(i,j,k-1)+Dld_eff(i,j,k))/2. + &
                               (alpha2_d(i,j,k-1)+alpha2_d(i,j,k))/2.*(Dgd_eff(i,j,k-1)+Dgd_eff(i,j,k))/2.)/dz* &
                               (ne(i,j,k-1)+ne(i,j,k))/2./dz
                          ff(i,j,k) = 0.
                       endif
                    endif
                    if(k.eq.Nz) then
                       if(tmp_wall(i,j,2)) then
                          aa(i,j,Nz) = aa(i,j,Nz) - &
                               dt*((tmp_alpha1_d(i,j,2)+alpha1_d(i,j,Nz))/2.*(tmp_Dld_eff(i,j,2)+Dld_eff(i,j,Nz))/2. + &
                               (tmp_alpha2_d(i,j,2)+alpha2_d(i,j,Nz))/2.*(tmp_Dgd_eff(i,j,2)+Dgd_eff(i,j,Nz))/2.)/dz* &
                               (tmp_ne(i,j,2)+ne(i,j,Nz))/2./dz
                          gg(i,j,Nz) = 0.
                       endif
                    else
                       if(wall(i,j,k+1)) then
                          aa(i,j,k) = aa(i,j,k) - &
                               dt*((alpha1_d(i,j,k+1)+alpha1_d(i,j,k))/2.*(Dld_eff(i,j,k+1)+Dld_eff(i,j,k))/2. + &
                               (alpha2_d(i,j,k+1)+alpha2_d(i,j,k))/2.*(Dgd_eff(i,j,k+1)+Dgd_eff(i,j,k))/2.)/dz* &
                               (ne(i,j,k+1)+ne(i,j,k))/2./dz
                          gg(i,j,k) = 0.
                       endif
                    endif
                 endif
              endif

           Enddo
        Enddo
     Enddo





     !! caloulate elements of matrix MM 


     Do k=1,Nz
        Do j=1,Ny
           Do i=1,Nx
              ma((k-1)*Nx*Ny+(j-1)*Nx+i) = aa(i,j,k)
              mb((k-1)*Nx*Ny+(j-1)*Nx+i) = bb(i,j,k)
              mc((k-1)*Nx*Ny+(j-1)*Nx+i) = cc(i,j,k)
              md((k-1)*Nx*Ny+(j-1)*Nx+i) = dd(i,j,k)
              me((k-1)*Nx*Ny+(j-1)*Nx+i) = ee(i,j,k)
              mf((k-1)*Nx*Ny+(j-1)*Nx+i) = ff(i,j,k)
              mg((k-1)*Nx*Ny+(j-1)*Nx+i) = gg(i,j,k)

              QQ((k-1)*Nx*Ny+(j-1)*Nx+i) = q(i,j,k)
              HOHO((k-1)*Nx*Ny+(j-1)*Nx+i) = td_o(i,j,k)

              WW((k-1)*Nx*Ny+(j-1)*Nx+i) = wall(i,j,k)
           Enddo
        Enddo
     Enddo










!!! solve td by successive overrelaxation (SOR)



        nn = 1

        Do while(nn.lt.1e5)

           if(WW(1)) then
              HH(1) = 0.
           else
              HH(1) = omega*(QQ(1)-mc(1)*HOHO(2)-me(1)*HOHO(1+Nx)-mg(1)*HOHO(1+Nx*Ny))/ &
                   (ma(1)+SMALL) + (1-omega)*HOHO(1)
           endif
           Do pp=2,Nx
              if(WW(pp)) then
                 HH(pp) = 0.
              else
                 HH(pp) = omega*(QQ(pp)-mb(pp)*HH(pp-1)-mc(pp)*HOHO(pp+1)-me(pp)*HOHO(pp+Nx)- &
                      mg(pp)*HOHO(pp+Nx*Ny))/(ma(pp)+SMALL) + (1-omega)*HOHO(pp)
              endif
           Enddo
           Do pp=1+Nx,Nx*Ny
              if(WW(pp)) then
                 HH(pp) = 0.
              else
                 HH(pp) = omega*(QQ(pp)-md(pp)*HH(pp-Nx)-mb(pp)*HH(pp-1)-mc(pp)*HOHO(pp+1)- &
                      me(pp)*HOHO(pp+Nx)-mg(pp)*HOHO(pp+Nx*Ny))/(ma(pp)+SMALL) + &
                      (1-omega)*HOHO(pp)
              endif
           Enddo
           Do pp=1+Nx*Ny,Nt-Nx*Ny
              if(WW(pp)) then
                 HH(pp) = 0.
              else
                 HH(pp) = omega*(QQ(pp)-mf(pp)*HH(pp-Nx*Ny)-md(pp)*HH(pp-Nx)- &
                      mb(pp)*HH(pp-1)-mc(pp)*HOHO(pp+1)-me(pp)*HOHO(pp+Nx)- &
                      mg(pp)*HOHO(pp+Nx*Ny))/(ma(pp)+SMALL) + (1-omega)*HOHO(pp)
              endif
           Enddo
           Do pp=Nt-Nx*Ny+1,Nt-Nx
              if(WW(pp)) then
                 HH(pp) = 0.
              else
                 HH(pp) = omega*(QQ(pp)-mf(pp)*HH(pp-Nx*Ny)-md(pp)*HH(pp-Nx)- &
                      mb(pp)*HH(pp-1)-mc(pp)*HOHO(pp+1)-me(pp)*HOHO(pp+Nx))/(ma(pp)+SMALL) + &
                      (1-omega)*HOHO(pp)
              endif
           Enddo
           Do pp=Nt-Nx+1,Nt-1
              if(WW(pp)) then
                 HH(pp) = 0.
              else
                 HH(pp) = omega*(QQ(pp)-mf(pp)*HH(pp-Nx*Ny)-md(pp)*HH(pp-Nx)- &
                      mb(pp)*HH(pp-1)-mc(pp)*HOHO(pp+1))/(ma(pp)+SMALL) + (1-omega)*HOHO(pp)
              endif
           Enddo
           if(WW(Nt)) then
              HH(Nt) = 0.
           else
              HH(Nt) = omega*(QQ(Nt)-mf(Nt)*HH(Nt-Nx*Ny)-md(Nt)*HH(Nt-Nx)- &
                   mb(Nt)*HH(Nt-1))/(ma(Nt)+SMALL) +(1-omega)*HOHO(Nt)
           endif


!!! calculate error
           HH_err_local = Sum(Abs(HH-HOHO))
           CALL MPI_ALLREDUCE(HH_err_local,HH_err,1,MPI_double_precision,MPI_SUM,MPI_COMM_WORLD,error)
           HH_err = HH_err/LX/LY/LZ


           If (HH_err.lt.HH_eps) exit


!!! update interface q

!!$        Do k=1,Nz
!!$           Do j=1,Ny
!!$              Do i=1,Nx
!!$                 td(i,j,k) = HH((k-1)*Nx*Ny+(j-1)*Nx+i)
!!$              Enddo
!!$           Enddo
!!$        Enddo

           Do j=1,Ny
              Do i=1,Nx
                 td(i,j,1) = HH((j-1)*Nx+i)
              Enddo
           Enddo

           Do j=1,Ny
              Do i=1,Nx
                 td(i,j,Nz) = HH((Nz-1)*Nx*Ny+(j-1)*Nx+i)
              Enddo
           Enddo




!!! exchange h during iteration

           call mpi_sendrecv(td(:,:,Nz),ilen,MPI_double_precision,inext,610, &
                tmp_td(:,:,1),ilen,MPI_double_precision,iprev,610,nallgrp,stat,error)

           call mpi_sendrecv(td(:,:,1),ilen,MPI_double_precision,iprev,620, &
                tmp_td(:,:,2),ilen,MPI_double_precision,inext,620,nallgrp,stat,error)


           CALL MPI_BARRIER(MPI_COMM_WORLD,error)


           
           if(mype.ne.0) then
              Do j=1,Ny
                 Do i=1,Nx
                    q(i,j,1) = ne(i,j,1)*td_o(i,j,1) + &
                         dt*theta(i,j,1)* &
                         (qd_0*b_o(i,j,1)*lc_o(i,j,1)/(Kc+lc_o(i,j,1))*lo_o(i,j,1)/(Ko+lo_o(i,j,1))) - &
                            (-dt*((tmp_alpha1_d(i,j,1)+alpha1_d(i,j,1))/2.*max(w(i,j,1),0.) + &
                             ((tmp_alpha1_d(i,j,1)+alpha1_d(i,j,1))/2.*(tmp_Dld_eff(i,j,1)+Dld_eff(i,j,1))/2. + &
                            (tmp_alpha2_d(i,j,1)+alpha2_d(i,j,1))/2.*(tmp_Dgd_eff(i,j,1)+Dgd_eff(i,j,1))/2.)/dz)* &
                            (tmp_ne(i,j,1)+ne(i,j,1))/2./dz)*tmp_td(i,j,1)
                    QQ((j-1)*Nx+i) = q(i,j,1)
                 Enddo
              Enddo
           endif

           if(mype.ne.npes-1) then
              Do j=1,Ny
                 Do i=1,Nx
                    q(i,j,Nz) = ne(i,j,Nz)*td_o(i,j,Nz) + &
                         dt*theta(i,j,Nz)* &
                         (qd_0*b_o(i,j,Nz)*lc_o(i,j,Nz)/(Kc+lc_o(i,j,Nz))*lo_o(i,j,Nz)/(Ko+lo_o(i,j,Nz))) - &
                            (dt*((tmp_alpha1_d(i,j,2)+alpha1_d(i,j,Nz))/2.*min(w(i,j,Nz+1),0.) - &
                            ((tmp_alpha1_d(i,j,2)+alpha1_d(i,j,Nz))/2.*(tmp_Dld_eff(i,j,2)+Dld_eff(i,j,Nz))/2. + &
                            (tmp_alpha2_d(i,j,2)+alpha2_d(i,j,Nz))/2.*(tmp_Dgd_eff(i,j,2)+Dgd_eff(i,j,Nz))/2.)/dz)* &
                            (tmp_ne(i,j,2)+ne(i,j,Nz))/2./dz)*tmp_td(i,j,2)
                    QQ((Nz-1)*Nx*Ny+(j-1)*Nx+i) = q(i,j,Nz)
                 Enddo
              Enddo
           endif



           HOHO = HH
           nn = nn+1

        ENDDO




        !! update h

        Do k=1,Nz
           Do j=1,Ny
              Do i=1,Nx
                 td(i,j,k) = HH((k-1)*Nx*Ny+(j-1)*Nx+i)
              Enddo
           Enddo
        Enddo





!!!! calculate ld and gd 

        Do k=1,Nz
           Do j=1,Ny
              Do i=1,Nx
                 ld(i,j,k) = beta1_d(i,j,k)*td(i,j,k)
                 gd(i,j,k) = beta2_d(i,j,k)*td(i,j,k)
              Enddo
           Enddo
        Enddo
        
        
        










!!! calculate carbon degradation rate 


!!$     Do k=1,Nz
!!$        Do j=1,Ny
!!$           Do i=1,Nx
!!$              if (wall(i,j,k)) then
!!$                 cr(i,j,k) = 0.0
!!$              else
!!$                 cr(i,j,k) = (c(i,j,k)-c_o(i,j,k))/dt
!!$              endif
!!$           Enddo
!!$        Enddo
!!$     Enddo









!!!!! calclate flux and mass of o2 through top surface

!!$     Do i=1,Nx
!!$        Do j=1,Ny
!!$              if(lobc_t.eq.1) then
!!$                 flo(i,j) = w_t*lo_t + 2.*Dlo_eff(i,j,Nz)*(lo(i,j,Nz)-lo_t)/dz 
!!$                 mlo(i,j) = flo(i,j)*dt*theta(i,j,Nz)*dx*dy
!!$              endif
!!$              if(gobc_t.eq.1) then
!!$                 fgo(i,j) = 2.*Dgo_eff(i,j,Nz)*(go(i,j,Nz)-go_t)/dz 
!!$                 mgo(i,j) = fgo(i,j)*dt*epsi(i,j,Nz)*dx*dy
!!$              endif
!!$        Enddo
!!$     Enddo



!!!!! calclate total mass transfer of o2 through top surface

!!$     Do i=1,Nx
!!$        Do j=1,Ny
!!$           tmo = tmo + mlo(i,j) + mgo(i,j)
!!$        Enddo
!!$     Enddo



!notde
!print*,'//',Dld_eff(nxx,nyy,Nz),ld(nxx,nyy,Nz),ld_t




!!!!! calclate flux and mass of co2-c through top surface
         
        if(mype.eq.npes-1) then
           Do i=1,Nx
              Do j=1,Ny
                 if(tdbc_t.eq.1) then
                    !              fld(i,j) = (w(i,j,Nz)*ld_t + 2.*Dld_eff(i,j,Nz)*(ld(i,j,Nz)-ld_t)/dz)/KPH*44./12.
                    fld(i,j) = w(i,j,Nz+1)*ld_t + 2.*Dld_eff(i,j,Nz)*(ld(i,j,Nz)-ld_t)/dz/KPH*theta(i,j,Nz)*dx*dy
                    mld(i,j) = fld(i,j)*dt
                 endif
                 if(tdbc_t.eq.1) then
                    fgd(i,j) = 2.*Dgd_eff(i,j,Nz)*(gd(i,j,Nz)-gd_t)/dz*epsi(i,j,Nz)*dx*dy
                    mgd(i,j) = fgd(i,j)*dt
                 endif
              Enddo
           Enddo
           tfld = sum(fld)
           tfgd = sum(fgd)
           tmd = tmd + sum(mld) + sum(mgd)
        endif


  CALL MPI_BCAST(tfld,1,MPI_double_precision,npes-1,MPI_COMM_WORLD,error)
  CALL MPI_BCAST(tfgd,1,MPI_double_precision,npes-1,MPI_COMM_WORLD,error)
  CALL MPI_BCAST(tmd,1,MPI_double_precision,npes-1,MPI_COMM_WORLD,error)
  CALL MPI_BARRIER(MPI_COMM_WORLD,error)


!!     tmd = (sum(c0-c) - sum(d-d0))/KPH*44./12.*dx*dy*dz
!!     tmd = (sum(c0-c) - sum(d-d0))*44./12.*dx*dy*dz





!!$!!! assume aqueous co2 and aqueous co2 reach equilibrium immediately and
!!$!!! satisfy Henry's law
!!$
!!$
!!$     Do i=1,Nx
!!$        Do j=1,Ny
!!$           Do k=1,Nz
!!$              if(.not.wall(i,j,k)) then
!!$                 d(i,j,k) = ld(i,j,k)*theta(i,j,k) + 12./44.*gd(i,j,k)*epsi(i,j,k)
!!$                 ld(i,j,k) = KPH*KH*R_igl*T_abs*d(i,j,k)/(KPH*KH*R_igl*T_abs*theta(i,j,k)+100.*epsi(i,j,k))
!!$                 gd(i,j,k) = 4400.*d(i,j,k)/12./(KPH*KH*R_igl*T_abs*theta(i,j,k)+100.*epsi(i,j,k))
!!$              endif
!!$           Enddo
!!$        Enddo
!!$     Enddo






     if(mod(istep,Print_step).eq.0) then
        Print*,'mype=',mype,'c---------','nxx=',nxx,'nyy=',nyy
        Print*,c(nxx,nyy,:)
        Print*,'mype=',mype,'lc---------','nxx=',nxx,'nyy=',nyy
        Print*,lc(nxx,nyy,:)
        Print*,'mype=',mype,'sc---------','nxx=',nxx,'nyy=',nyy
        Print*,sc(nxx,nyy,:)
        Print*,'mype=',mype,'b---------','nxx=',nxx,'nyy=',nyy
        Print*,b(nxx,nyy,:)
        Print*,'mype=',mype,'to---------','nxx=',nxx,'nyy=',nyy
        Print*,to(nxx,nyy,:)
        Print*,'mype=',mype,'lo---------','nxx=',nxx,'nyy=',nyy
        Print*,lo(nxx,nyy,:)
        Print*,'mype=',mype,'go---------','nxx=',nxx,'nyy=',nyy
        Print*,go(nxx,nyy,:)
        Print*,'mype=',mype,'td---------','nxx=',nxx,'nyy=',nyy
        Print*,td(nxx,nyy,:)
        Print*,'mype=',mype,'ld---------','nxx=',nxx,'nyy=',nyy
        Print*,ld(nxx,nyy,:)
        Print*,'mype=',mype,'gd---------','nxx=',nxx,'nyy=',nyy
        Print*,gd(nxx,nyy,:)
     endif





!!!*********** output data ***********

     IF(mod(istep,Print_step)==0) THEN

        clocalmax = maxval(c)
        CALL MPI_REDUCE(clocalmax,cmax,1,MPI_double_precision,MPI_MAX,0,MPI_COMM_WORLD,error)
        lclocalmax = maxval(lc)
        CALL MPI_REDUCE(lclocalmax,lcmax,1,MPI_double_precision,MPI_MAX,0,MPI_COMM_WORLD,error)
        sclocalmax = maxval(sc)
        CALL MPI_REDUCE(sclocalmax,scmax,1,MPI_double_precision,MPI_MAX,0,MPI_COMM_WORLD,error)
        blocalmax = maxval(b)
        CALL MPI_REDUCE(blocalmax,bmax,1,MPI_double_precision,MPI_MAX,0,MPI_COMM_WORLD,error)
        tolocalmax = maxval(to)
        CALL MPI_REDUCE(tolocalmax,tomax,1,MPI_double_precision,MPI_MAX,0,MPI_COMM_WORLD,error)
        lolocalmax = maxval(lo)
        CALL MPI_REDUCE(lolocalmax,lomax,1,MPI_double_precision,MPI_MAX,0,MPI_COMM_WORLD,error)
        golocalmax = maxval(go)
        CALL MPI_REDUCE(golocalmax,gomax,1,MPI_double_precision,MPI_MAX,0,MPI_COMM_WORLD,error)
        tdlocalmax = maxval(td)
        CALL MPI_REDUCE(tdlocalmax,tdmax,1,MPI_double_precision,MPI_MAX,0,MPI_COMM_WORLD,error)
        ldlocalmax = maxval(ld)
        CALL MPI_REDUCE(ldlocalmax,ldmax,1,MPI_double_precision,MPI_MAX,0,MPI_COMM_WORLD,error)
        gdlocalmax = maxval(gd)
        CALL MPI_REDUCE(gdlocalmax,gdmax,1,MPI_double_precision,MPI_MAX,0,MPI_COMM_WORLD,error)

!!$        colocalsum = sum(c_o)
!!$        CALL MPI_REDUCE(colocalsum,cosum,1,MPI_double_precision,MPI_SUM,0,MPI_COMM_WORLD,error)
        clocalsum = sum(c)
        CALL MPI_REDUCE(clocalsum,csum,1,MPI_double_precision,MPI_SUM,0,MPI_COMM_WORLD,error)
        lclocalsum = sum(lc)
        CALL MPI_REDUCE(lclocalsum,lcsum,1,MPI_double_precision,MPI_SUM,0,MPI_COMM_WORLD,error)
        sclocalsum = sum(sc)
        CALL MPI_REDUCE(sclocalsum,scsum,1,MPI_double_precision,MPI_SUM,0,MPI_COMM_WORLD,error)
        blocalsum = sum(b)
        CALL MPI_REDUCE(blocalsum,bsum,1,MPI_double_precision,MPI_SUM,0,MPI_COMM_WORLD,error)
        tolocalsum = sum(to)
        CALL MPI_REDUCE(tolocalsum,tosum,1,MPI_double_precision,MPI_SUM,0,MPI_COMM_WORLD,error)
        lolocalsum = sum(lo)
        CALL MPI_REDUCE(lolocalsum,losum,1,MPI_double_precision,MPI_SUM,0,MPI_COMM_WORLD,error)
        golocalsum = sum(go)
        CALL MPI_REDUCE(golocalsum,gosum,1,MPI_double_precision,MPI_SUM,0,MPI_COMM_WORLD,error)
        tdlocalsum = sum(td)
        CALL MPI_REDUCE(tdlocalsum,tdsum,1,MPI_double_precision,MPI_SUM,0,MPI_COMM_WORLD,error)
        ldlocalsum = sum(ld)
        CALL MPI_REDUCE(ldlocalsum,ldsum,1,MPI_double_precision,MPI_SUM,0,MPI_COMM_WORLD,error)
        gdlocalsum = sum(gd)
        CALL MPI_REDUCE(gdlocalsum,gdsum,1,MPI_double_precision,MPI_SUM,0,MPI_COMM_WORLD,error)

        CALL MPI_BARRIER(MPI_COMM_WORLD,error)



        IF(mype.eq.0) THEN
        print*,'maximum concentration of total organic carbon: ',cmax
        print*,'cumulative concentration of total organic carbon: ',csum
        print*,'maximum concentration of aqueous organic carbon: ',lcmax
        print*,'cumulative concentration of aqueous organic carbon: ',lcsum
        print*,'maximum concentration of soil organic carbon: ',scmax
        print*,'cumulative concentration of soil organic carbon: ',scsum
        print*,'maximum concentration of bacteria: ',bmax
        print*,'cumulative concentration of bacteria: ',bsum
        print*,'maximum concentration of o2: ',tomax
        print*,'cumulative concentration of o2: ',tosum
        print*,'maximum concentration of aqueous o2: ',lomax
        print*,'cumulative concentration of aqueous o2: ',losum
        print*,'maximum concentration of air o2: ',gomax
        print*,'cumulative concentration of air o2: ',gosum
        print*,'maximum concentration of total inorganic carbon: ',tdmax
        print*,'cumulative concentration of total inorganic carbon: ',tdsum
        print*,'maximum concentration of aqueous inorganic carbon: ',ldmax
        print*,'cumulative concentration of aqueous inorganic carbon: ',ldsum
        print*,'maximum concentration of air co2-c: ',gdmax
        print*,'cumulative concentration of air co2-c: ',gdsum
!        print*,'cumulative flux of aqueous o2 through surface: ',sum(flo)
!        print*,'cumulative flux of air o2 through surface: ',sum(fgo)
!        print*,'total mass of o2 through surface: ',tmo
        print*,'cumulative flux of aqueous co2-c through surface: ',tfld
        print*,'cumulative flux of air co2-c through surface: ',tfgd
        print*,'total mass of co2-c through surface: ',tmd
        endif

!!$        IF(mype.eq.npes-1) THEN
!        print*,'cumulative flux of aqueous o2 on surface: ',sum(flo)
!        print*,'total mass of o2 through top surface: ',tmo
!!$           print*,'cumulative flux of aqueous co2 on surface: ',tfd
!!$           print*,'total mass of co2 through top surface: ',tmd
!!$        endif


     endif


!!$!!!*********** write output files ***********

     IF(mod(istep,File_step)==0) THEN
        call OutputData !(c,lc,sc,b,o,lo,go,d,ld,gd,flo,fgo,fld,fgd,Nz,Nx,Ny,Nstart,MEMO,OUT_DIR,PMname,istep)
     endif




!!!!! update 

     lc_o = lc
     sc_o = sc
     b_o = b
     lo_o = lo
     to_o = to
     td_o = td
!



  EndDo



  endtime = MPI_WTIME()


  print*,'End computation'
  print '("Computation time = ",f10.4," seconds.")',endtime-starttime


  CALL MPI_FINALIZE (error)


END PROGRAM USMSSolver



!--------------------------------------------------------------
SUBROUTINE input
  !--------------------------------------------------------------    

  Use GLOBAL
  Implicit None
  include 'mpif.h'

  call MPI_INIT(error)
  if (error.ne.0) stop "ERROR: can't initialize mpi" 
  call MPI_COMM_SIZE(MPI_COMM_WORLD,npes,error)
  call MPI_COMM_RANK(MPI_COMM_WORLD,mype,error)
  nallgrp=MPI_COMM_WORLD
  
  
  if (mype.eq.0) then

     Read(*,*),LX,LY,LZ
     Write(*,*)'## INPUT'
     Write(*,*)'## Spatial resolution LX, LY, LZ = ',LX,LY,LZ
     Read(*,*),dx,dy,dz
     Write(*,*)'## Spatial length dx,dy,dz = ',dx,dy,dz
     Read(*,*),dt
     Write(*,*)'## Time length dt = ',dt

     Read(*,*),lcbc_d,lcbc_t
     Write(*,*)'## Carbon conc boundary type on down and top planes = ',lcbc_d,lcbc_t
     Read(*,*),lc_d,lc_t
     Write(*,*)'## Conc on down and top planes lc_d,lc_t = ',lc_d,lc_t
!!$  Read(*,*),bbc_d,bbc_t
!!$  Write(*,*)'## Bacteria conc boundary type on down and top planes = ',bbc_d,bbc_t
!!$  Read(*,*),b_d,b_t
!!$  Write(*,*)'## Conc on down and top planes b_d,b_t = ',b_d,b_t
     Read(*,*),tobc_d,tobc_t
     Write(*,*)'## total oxygen conc boundary type on down and top planes = ',tobc_d,tobc_t
     Read(*,*),lo_d,lo_t
     Write(*,*)'## aqueous oxygen conc on down and top planes lo_d,lo_t = ',lo_d,lo_t
!     Read(*,*),gobc_d,gobc_t
!     Write(*,*)'## Gas oxygen conc boundary type on down and top planes = ',gobc_d,gobc_t
     Read(*,*),go_d,go_t
     Write(*,*)'## Gas oxygen conc on down and top planes go_d,go_t = ',go_d,go_t
     Read(*,*),tdbc_d,tdbc_t
     Write(*,*)'## TD conc boundary type on down and top planes = ',tdbc_d,tdbc_t
     Read(*,*),ld_d,ld_t
     Write(*,*)'## gas co2 conc on down and top planes lo_d,lo_t = ',ld_d,ld_t
!     Read(*,*),gdbc_d,gdbc_t
!     Write(*,*)'## Gas co2 conc boundary type on down and top planes = ',gdbc_d,gdbc_t
     Read(*,*),gd_d,gd_t
     Write(*,*)'## Gas co2 conc on down and top planes go_d,go_t = ',gd_d,gd_t

     Read(*,*),Dlc
     Write(*,*)'## Melocular diffusion coefficient of carbon = ',Dlc
!!$  Read(*,*),Db
!!$  Write(*,*)'## Melocular diffusion coefficient of bacteria = ',Db
     Read(*,*),Dlo
     Write(*,*)'## Melocular diffusion coefficient of oxygen in aqueous = ',Dlo
     Read(*,*),Dgo
     Write(*,*)'## Melocular diffusion coefficient of oxygen in air = ',Dgo
     Read(*,*),Dld
     Write(*,*)'## Melocular diffusion coefficient of co2 in aqueous = ',Dld
     Read(*,*),Dgd
     Write(*,*)'## Melocular diffusion coefficient of co2 in air = ',Dgd
     Read(*,*),mmm_l,nnn_l
     Write(*,*)'## exponent for solute =',mmm_l,nnn_l
     Read(*,*),mmm_g,nnn_g
     Write(*,*)'## exponent for gas =',mmm_g,nnn_g
     Read(*,*),Kcc
     Write(*,*)'## adsorption/desorption equilibrium constant Kcc =',Kcc
     Read(*,*),SA
     Write(*,*)'## specific surface area SA =',SA
     Read(*,*),a1,b1
     Write(*,*)'## coefficient for transfer rate =',a1,b1

     Read(*,*),Kho
     Write(*,*)'## Henrys law constant of oxygen Kho =',Kho
     Read(*,*),Khd
     Write(*,*)'## Henrys law constant of co2 Khd =',Khd
     Read(*,*),Kph
     Write(*,*)'## equilibrium coefficient related to pH =',Kph
     Read(*,*),rhos
     Write(*,*)'## soil density rhos =',rhos

     Read(*,*),qc_0
     Write(*,*)'## Maximum rate of contaminant consumption = ',qc_0
     Read(*,*),qo_0
     Write(*,*)'## Maximum rate of oxygen consumption = ',qo_0
     Read(*,*),qd_0
     Write(*,*)'## Maximum rate of co2 production = ',qd_0
     Read(*,*),Kc
     Write(*,*)'## Half-saturation coefficient of carbon = ',Kc
     Read(*,*),Ko
     Write(*,*)'## Half-saturation coefficient of oxygen = ',Ko
     Read(*,*),Yd
     Write(*,*)'## Yield coefficient = ',Yd
     Read(*,*),Dk
     Write(*,*)'## Bacterial decay rate = ',Dk


  Read(*,*),OM
  Write(*,*)'## Organic matter = ',OM
     Read(*,*),omega
     Write(*,*)'## Over-relaxation factor = ',omega
     Read(*,*),HH_eps
     Write(*,*)'## inner iteration criteria = ',HH_eps

     Read(*,*), Max_step,Print_step,File_step
     Write(*,*)'## Maximum # of time steps =',Max_step
     Write(*,*)'## # of time steps between print =',Print_step
     Write(*,*)'## # of time steps between file output =',File_step
     Read(*,*), PMname
     Write(*,*)'## Porous medium file: ',Trim(PMname)
     Read(*,*), IN_DIR
     Write(*,*)'## Input directory:  ',Trim(IN_DIR)
     Read(*,*), OUT_DIR
     Write(*,*)'## Output directory: ',Trim(OUT_DIR)
     Read(*,*),t0
     Write(*,*)'## initial time steps = ',t0



  endif

  inext=mod(mype+npes+1,npes)
  iprev=mod(mype+npes-1,npes)


  CALL MPI_BCAST(LX,1,MPI_INTEGER,0,MPI_COMM_WORLD,error)
  CALL MPI_BCAST(LY,1,MPI_INTEGER,0,MPI_COMM_WORLD,error)
  CALL MPI_BCAST(LZ,1,MPI_INTEGER,0,MPI_COMM_WORLD,error)

  CALL MPI_BCAST(dx,1,MPI_double_precision,0,MPI_COMM_WORLD,error)
  CALL MPI_BCAST(dy,1,MPI_double_precision,0,MPI_COMM_WORLD,error)
  CALL MPI_BCAST(dz,1,MPI_double_precision,0,MPI_COMM_WORLD,error)
  CALL MPI_BCAST(dt,1,MPI_double_precision,0,MPI_COMM_WORLD,error)

  CALL MPI_BCAST(lcbc_d,1,MPI_INTEGER,0,MPI_COMM_WORLD,error)
  CALL MPI_BCAST(lcbc_t,1,MPI_INTEGER,0,MPI_COMM_WORLD,error)
  CALL MPI_BCAST(lc_d,1,MPI_double_precision,0,MPI_COMM_WORLD,error)
  CALL MPI_BCAST(lc_t,1,MPI_double_precision,0,MPI_COMM_WORLD,error)

  CALL MPI_BCAST(tobc_d,1,MPI_INTEGER,0,MPI_COMM_WORLD,error)
  CALL MPI_BCAST(tobc_t,1,MPI_INTEGER,0,MPI_COMM_WORLD,error)
  CALL MPI_BCAST(lo_d,1,MPI_double_precision,0,MPI_COMM_WORLD,error)
  CALL MPI_BCAST(lo_t,1,MPI_double_precision,0,MPI_COMM_WORLD,error)
  CALL MPI_BCAST(go_d,1,MPI_double_precision,0,MPI_COMM_WORLD,error)
  CALL MPI_BCAST(go_t,1,MPI_double_precision,0,MPI_COMM_WORLD,error)

  CALL MPI_BCAST(tdbc_d,1,MPI_INTEGER,0,MPI_COMM_WORLD,error)
  CALL MPI_BCAST(tdbc_t,1,MPI_INTEGER,0,MPI_COMM_WORLD,error)
  CALL MPI_BCAST(ld_d,1,MPI_double_precision,0,MPI_COMM_WORLD,error)
  CALL MPI_BCAST(ld_t,1,MPI_double_precision,0,MPI_COMM_WORLD,error)
!  CALL MPI_BCAST(gdbc_d,1,MPI_INTEGER,0,MPI_COMM_WORLD,error)
!  CALL MPI_BCAST(gdbc_t,1,MPI_INTEGER,0,MPI_COMM_WORLD,error)
  CALL MPI_BCAST(gd_d,1,MPI_double_precision,0,MPI_COMM_WORLD,error)
  CALL MPI_BCAST(gd_t,1,MPI_double_precision,0,MPI_COMM_WORLD,error)

  CALL MPI_BCAST(Dlc,1,MPI_double_precision,0,MPI_COMM_WORLD,error)
  CALL MPI_BCAST(Dlo,1,MPI_double_precision,0,MPI_COMM_WORLD,error)
  CALL MPI_BCAST(Dgo,1,MPI_double_precision,0,MPI_COMM_WORLD,error)
  CALL MPI_BCAST(Dld,1,MPI_double_precision,0,MPI_COMM_WORLD,error)
  CALL MPI_BCAST(Dgd,1,MPI_double_precision,0,MPI_COMM_WORLD,error)

  CALL MPI_BCAST(mmm_l,1,MPI_double_precision,0,MPI_COMM_WORLD,error)
  CALL MPI_BCAST(nnn_l,1,MPI_double_precision,0,MPI_COMM_WORLD,error)
  CALL MPI_BCAST(mmm_g,1,MPI_double_precision,0,MPI_COMM_WORLD,error)
  CALL MPI_BCAST(nnn_g,1,MPI_double_precision,0,MPI_COMM_WORLD,error)
  CALL MPI_BCAST(Kcc,1,MPI_double_precision,0,MPI_COMM_WORLD,error)
  CALL MPI_BCAST(SA,1,MPI_double_precision,0,MPI_COMM_WORLD,error)
  CALL MPI_BCAST(a1,1,MPI_double_precision,0,MPI_COMM_WORLD,error)
  CALL MPI_BCAST(b1,1,MPI_double_precision,0,MPI_COMM_WORLD,error)
  CALL MPI_BCAST(Kho,1,MPI_double_precision,0,MPI_COMM_WORLD,error)
  CALL MPI_BCAST(Khd,1,MPI_double_precision,0,MPI_COMM_WORLD,error)
  CALL MPI_BCAST(KPH,1,MPI_double_precision,0,MPI_COMM_WORLD,error)
  CALL MPI_BCAST(rhos,1,MPI_double_precision,0,MPI_COMM_WORLD,error)

  CALL MPI_BCAST(qc_0,1,MPI_double_precision,0,MPI_COMM_WORLD,error)
  CALL MPI_BCAST(qo_0,1,MPI_double_precision,0,MPI_COMM_WORLD,error)
  CALL MPI_BCAST(qd_0,1,MPI_double_precision,0,MPI_COMM_WORLD,error)
  CALL MPI_BCAST(Kc,1,MPI_double_precision,0,MPI_COMM_WORLD,error)
  CALL MPI_BCAST(Ko,1,MPI_double_precision,0,MPI_COMM_WORLD,error)
  CALL MPI_BCAST(Yd,1,MPI_double_precision,0,MPI_COMM_WORLD,error)
  CALL MPI_BCAST(Dk,1,MPI_double_precision,0,MPI_COMM_WORLD,error)
!!$
!!$  CALL MPI_BCAST(tau1,1,MPI_double_precision,0,MPI_COMM_WORLD,error)
!!$  CALL MPI_BCAST(lambda1,1,MPI_double_precision,0,MPI_COMM_WORLD,error)
!!$  CALL MPI_BCAST(Sm,1,MPI_double_precision,0,MPI_COMM_WORLD,error)

  CALL MPI_BCAST(OM,1,MPI_double_precision,0,MPI_COMM_WORLD,error)

  CALL MPI_BCAST(omega,1,MPI_double_precision,0,MPI_COMM_WORLD,error)
  CALL MPI_BCAST(HH_eps,1,MPI_double_precision,0,MPI_COMM_WORLD,error)
!!$  CALL MPI_BCAST(cc_eps,1,MPI_double_precision,0,MPI_COMM_WORLD,error)
!!$  CALL MPI_BCAST(lo_eps,1,MPI_double_precision,0,MPI_COMM_WORLD,error)
!!$  CALL MPI_BCAST(ld_eps,1,MPI_double_precision,0,MPI_COMM_WORLD,error)


  CALL MPI_BCAST(Max_step,1,MPI_INTEGER,0,MPI_COMM_WORLD,error)
  CALL MPI_BCAST(Print_step,1,MPI_INTEGER,0,MPI_COMM_WORLD,error)
  CALL MPI_BCAST(File_step,1,MPI_INTEGER,0,MPI_COMM_WORLD,error)

  CALL MPI_BCAST(PMname,64,MPI_CHARACTER,0,MPI_COMM_WORLD,error)
  CALL MPI_BCAST(IN_DIR,64,MPI_CHARACTER,0,MPI_COMM_WORLD,error)
  CALL MPI_BCAST(OUT_DIR,64,MPI_CHARACTER,0,MPI_COMM_WORLD,error)
  CALL MPI_BCAST(t0,1,MPI_INTEGER,0,MPI_COMM_WORLD,error)


  CALL MPI_BARRIER(MPI_COMM_WORLD,error)

  Nstart = mype * LZ / npes + 1
  Nx = LX
  Ny = LY
  Nz = (mype+1) * LZ / npes - mype * LZ / npes

  print *,'Processor =',mype,'Domain start =',Nstart,'Domain size =',Nz



END SUBROUTINE input





!--------------------------------------------------------------
Subroutine ReadDensity(den,myblock,LX,LY,sta,MEMO,FileName)
  !--------------------------------------------------------------      
  Implicit None
  Integer:: LX,LY,myblock,MEMO,RecLength,IO,idz,sta
  Double precision,Intent(OUT):: den(LX,LY,myblock)
  Character*(*),Intent(IN):: FileName

  RecLength=MEMO*LX*LY
  Open(30,ACCESS='DIRECT',STATUS='OLD',RECL=RecLength, &
       FILE=FileName,IOSTAT=IO)
  If (IO.Ne.0) Then
     Write(*,*)'Error opening density file, IO = ',IO
     Stop
  End If
  Do idz = 1,myblock
     Read(30,REC=sta+idz-1)den(:,:,idz)
  End Do
  Close(30)
End Subroutine ReadDensity






!--------------------------------------------------------------
Subroutine OutputData 
  !--------------------------------------------------------------   

  use GLOBAL
  implicit none
  include 'mpif.h'


  Rootname = trim(OUT_DIR)//'/'//trim(PMname)//'_c_'// &
       't'// &   
       CHAR(mod(istep,10000000)/1000000+48)//&
       CHAR(mod(istep,1000000)/100000+48)//&
       CHAR(mod(istep,100000)/10000+48)//&
       CHAR(mod(istep,10000)/1000+48)//&
       CHAR(mod(istep,1000)/100+48)//&
       CHAR(mod(istep,100)/10+48)//&
       CHAR(mod(istep,10)+48)
  CALL WriteDEN(c(:,:,:),Nz,Nx,Ny,Nstart,MEMO,trim(Rootname)//'_'//&
       CHAR(mype/1000+48)//CHAR(mod(mype,1000)/100+48)//&
       CHAR(mod(mype,100)/10+48)//CHAR(mod(mype,10)+48))

  CALL MPI_BARRIER(MPI_COMM_WORLD,error)

  if(mype.eq.0) then
     CALL SYSTEM('cat '//trim(Rootname)//'_* > ' &
          //trim(Rootname))
     CALL SYSTEM('rm '//trim(Rootname)//'_*')
  endif


  Rootname = trim(OUT_DIR)//'/'//trim(PMname)//'_lc_'// &
       't'// &   
       CHAR(mod(istep,10000000)/1000000+48)//&
       CHAR(mod(istep,1000000)/100000+48)//&
       CHAR(mod(istep,100000)/10000+48)//&
       CHAR(mod(istep,10000)/1000+48)//&
       CHAR(mod(istep,1000)/100+48)//&
       CHAR(mod(istep,100)/10+48)//&
       CHAR(mod(istep,10)+48)
  CALL WriteDEN(lc(:,:,:),Nz,Nx,Ny,Nstart,MEMO,trim(Rootname)//'_'//&
       CHAR(mype/1000+48)//CHAR(mod(mype,1000)/100+48)//&
       CHAR(mod(mype,100)/10+48)//CHAR(mod(mype,10)+48))

  CALL MPI_BARRIER(MPI_COMM_WORLD,error)

  if(mype.eq.0) then
     CALL SYSTEM('cat '//trim(Rootname)//'_* > ' &
          //trim(Rootname))
     CALL SYSTEM('rm '//trim(Rootname)//'_*')
  endif

  Rootname = trim(OUT_DIR)//'/'//trim(PMname)//'_sc_'// &
       't'// &   
       CHAR(mod(istep,10000000)/1000000+48)//&
       CHAR(mod(istep,1000000)/100000+48)//&
       CHAR(mod(istep,100000)/10000+48)//&
       CHAR(mod(istep,10000)/1000+48)//&
       CHAR(mod(istep,1000)/100+48)//&
       CHAR(mod(istep,100)/10+48)//&
       CHAR(mod(istep,10)+48)
  CALL WriteDEN(sc(:,:,:),Nz,Nx,Ny,Nstart,MEMO,trim(Rootname)//'_'//&
       CHAR(mype/1000+48)//CHAR(mod(mype,1000)/100+48)//&
       CHAR(mod(mype,100)/10+48)//CHAR(mod(mype,10)+48))

  CALL MPI_BARRIER(MPI_COMM_WORLD,error)

  if(mype.eq.0) then
     CALL SYSTEM('cat '//trim(Rootname)//'_* > ' &
          //trim(Rootname))
     CALL SYSTEM('rm '//trim(Rootname)//'_*')
  endif

  Rootname = trim(OUT_DIR)//'/'//trim(PMname)//'_b_'// &
       't'// &   
       CHAR(mod(istep,10000000)/1000000+48)//&
       CHAR(mod(istep,1000000)/100000+48)//&
       CHAR(mod(istep,100000)/10000+48)//&
       CHAR(mod(istep,10000)/1000+48)//&
       CHAR(mod(istep,1000)/100+48)//&
       CHAR(mod(istep,100)/10+48)//&
       CHAR(mod(istep,10)+48)
  CALL WriteDEN(b(:,:,:),Nz,Nx,Ny,Nstart,MEMO,trim(Rootname)//'_'//&
       CHAR(mype/1000+48)//CHAR(mod(mype,1000)/100+48)//&
       CHAR(mod(mype,100)/10+48)//CHAR(mod(mype,10)+48))

  CALL MPI_BARRIER(MPI_COMM_WORLD,error)

  if(mype.eq.0) then
     CALL SYSTEM('cat '//trim(Rootname)//'_* > ' &
          //trim(Rootname))
     CALL SYSTEM('rm '//trim(Rootname)//'_*')
  endif

  Rootname = trim(OUT_DIR)//'/'//trim(PMname)//'_to_'// &
       't'// &   
       CHAR(mod(istep,10000000)/1000000+48)//&
       CHAR(mod(istep,1000000)/100000+48)//&
       CHAR(mod(istep,100000)/10000+48)//&
       CHAR(mod(istep,10000)/1000+48)//&
       CHAR(mod(istep,1000)/100+48)//&
       CHAR(mod(istep,100)/10+48)//&
       CHAR(mod(istep,10)+48)
  CALL WriteDEN(to(:,:,:),Nz,Nx,Ny,Nstart,MEMO,trim(Rootname)//'_'//&
       CHAR(mype/1000+48)//CHAR(mod(mype,1000)/100+48)//&
       CHAR(mod(mype,100)/10+48)//CHAR(mod(mype,10)+48))

  CALL MPI_BARRIER(MPI_COMM_WORLD,error)

  if(mype.eq.0) then
     CALL SYSTEM('cat '//trim(Rootname)//'_* > ' &
          //trim(Rootname))
     CALL SYSTEM('rm '//trim(Rootname)//'_*')
  endif

  Rootname = trim(OUT_DIR)//'/'//trim(PMname)//'_lo_'// &
       't'// &   
       CHAR(mod(istep,10000000)/1000000+48)//&
       CHAR(mod(istep,1000000)/100000+48)//&
       CHAR(mod(istep,100000)/10000+48)//&
       CHAR(mod(istep,10000)/1000+48)//&
       CHAR(mod(istep,1000)/100+48)//&
       CHAR(mod(istep,100)/10+48)//&
       CHAR(mod(istep,10)+48)
  CALL WriteDEN(lo(:,:,:),Nz,Nx,Ny,Nstart,MEMO,trim(Rootname)//'_'//&
       CHAR(mype/1000+48)//CHAR(mod(mype,1000)/100+48)//&
       CHAR(mod(mype,100)/10+48)//CHAR(mod(mype,10)+48))

  CALL MPI_BARRIER(MPI_COMM_WORLD,error)

  if(mype.eq.0) then
     CALL SYSTEM('cat '//trim(Rootname)//'_* > ' &
          //trim(Rootname))
     CALL SYSTEM('rm '//trim(Rootname)//'_*')
  endif


  Rootname = trim(OUT_DIR)//'/'//trim(PMname)//'_go_'// &
       't'// &   
       CHAR(mod(istep,10000000)/1000000+48)//&
       CHAR(mod(istep,1000000)/100000+48)//&
       CHAR(mod(istep,100000)/10000+48)//&
       CHAR(mod(istep,10000)/1000+48)//&
       CHAR(mod(istep,1000)/100+48)//&
       CHAR(mod(istep,100)/10+48)//&
       CHAR(mod(istep,10)+48)
  CALL WriteDEN(go(:,:,:),Nz,Nx,Ny,Nstart,MEMO,trim(Rootname)//'_'//&
       CHAR(mype/1000+48)//CHAR(mod(mype,1000)/100+48)//&
       CHAR(mod(mype,100)/10+48)//CHAR(mod(mype,10)+48))

  CALL MPI_BARRIER(MPI_COMM_WORLD,error)

  if(mype.eq.0) then
     CALL SYSTEM('cat '//trim(Rootname)//'_* > ' &
          //trim(Rootname))
     CALL SYSTEM('rm '//trim(Rootname)//'_*')
  endif


  Rootname = trim(OUT_DIR)//'/'//trim(PMname)//'_td_'// &
       't'// &   
       CHAR(mod(istep,10000000)/1000000+48)//&
       CHAR(mod(istep,1000000)/100000+48)//&
       CHAR(mod(istep,100000)/10000+48)//&
       CHAR(mod(istep,10000)/1000+48)//&
       CHAR(mod(istep,1000)/100+48)//&
       CHAR(mod(istep,100)/10+48)//&
       CHAR(mod(istep,10)+48)
  CALL WriteDEN(td(:,:,:),Nz,Nx,Ny,Nstart,MEMO,trim(Rootname)//'_'//&
       CHAR(mype/1000+48)//CHAR(mod(mype,1000)/100+48)//&
       CHAR(mod(mype,100)/10+48)//CHAR(mod(mype,10)+48))

  CALL MPI_BARRIER(MPI_COMM_WORLD,error)

  if(mype.eq.0) then
     CALL SYSTEM('cat '//trim(Rootname)//'_* > ' &
          //trim(Rootname))
     CALL SYSTEM('rm '//trim(Rootname)//'_*')
  endif


  Rootname = trim(OUT_DIR)//'/'//trim(PMname)//'_ld_'// &
       't'// &   
       CHAR(mod(istep,10000000)/1000000+48)//&
       CHAR(mod(istep,1000000)/100000+48)//&
       CHAR(mod(istep,100000)/10000+48)//&
       CHAR(mod(istep,10000)/1000+48)//&
       CHAR(mod(istep,1000)/100+48)//&
       CHAR(mod(istep,100)/10+48)//&
       CHAR(mod(istep,10)+48)
  CALL WriteDEN(ld(:,:,:),Nz,Nx,Ny,Nstart,MEMO,trim(Rootname)//'_'//&
       CHAR(mype/1000+48)//CHAR(mod(mype,1000)/100+48)//&
       CHAR(mod(mype,100)/10+48)//CHAR(mod(mype,10)+48))

  CALL MPI_BARRIER(MPI_COMM_WORLD,error)

  if(mype.eq.0) then
     CALL SYSTEM('cat '//trim(Rootname)//'_* > ' &
          //trim(Rootname))
     CALL SYSTEM('rm '//trim(Rootname)//'_*')
  endif


  Rootname = trim(OUT_DIR)//'/'//trim(PMname)//'_gd_'// &
       't'// &   
       CHAR(mod(istep,10000000)/1000000+48)//&
       CHAR(mod(istep,1000000)/100000+48)//&
       CHAR(mod(istep,100000)/10000+48)//&
       CHAR(mod(istep,10000)/1000+48)//&
       CHAR(mod(istep,1000)/100+48)//&
       CHAR(mod(istep,100)/10+48)//&
       CHAR(mod(istep,10)+48)
  CALL WriteDEN(gd(:,:,:),Nz,Nx,Ny,Nstart,MEMO,trim(Rootname)//'_'//&
       CHAR(mype/1000+48)//CHAR(mod(mype,1000)/100+48)//&
       CHAR(mod(mype,100)/10+48)//CHAR(mod(mype,10)+48))

  CALL MPI_BARRIER(MPI_COMM_WORLD,error)

  if(mype.eq.0) then
     CALL SYSTEM('cat '//trim(Rootname)//'_* > ' &
          //trim(Rootname))
     CALL SYSTEM('rm '//trim(Rootname)//'_*')
  endif

!!$  Rootname = trim(OUT_DIR)//'/'//trim(PMname)//'_cr_'// &
!!$       't'// &   
!!$       CHAR(mod(istep,10000000)/1000000+48)//&
!!$       CHAR(mod(istep,1000000)/100000+48)//&
!!$       CHAR(mod(istep,100000)/10000+48)//&
!!$       CHAR(mod(istep,10000)/1000+48)//&
!!$       CHAR(mod(istep,1000)/100+48)//&
!!$       CHAR(mod(istep,100)/10+48)//&
!!$       CHAR(mod(istep,10)+48)
!!$  CALL WriteDEN(cr(:,:,:),Nz,Nx,Ny,Nstart,MEMO,trim(Rootname)//'_'//&
!!$       CHAR(mype/1000+48)//CHAR(mod(mype,1000)/100+48)//&
!!$       CHAR(mod(mype,100)/10+48)//CHAR(mod(mype,10)+48))
!!$
!!$  CALL MPI_BARRIER(MPI_COMM_WORLD,error)
!!$
!!$  if(mype.eq.0) then
!!$     CALL SYSTEM('cat '//trim(Rootname)//'_* > ' &
!!$          //trim(Rootname))
!!$     CALL SYSTEM('rm '//trim(Rootname)//'_*')
!!$  endif


!!$  Rootname = trim(OUT_DIR)//'/'//trim(FileName)//'_flo_'// &
!!$       't'// &   
!!$       CHAR(mod(istep,10000000)/1000000+48)//&
!!$       CHAR(mod(istep,1000000)/100000+48)//&
!!$       CHAR(mod(istep,100000)/10000+48)//&
!!$       CHAR(mod(istep,10000)/1000+48)//&
!!$       CHAR(mod(istep,1000)/100+48)//&
!!$       CHAR(mod(istep,100)/10+48)//&
!!$       CHAR(mod(istep,10)+48)
!!$  CALL WriteFlux(flo(:,:),Nx,Ny,Nstart,MEMO,trim(Rootname))
!!$
!!$  Rootname = trim(OUT_DIR)//'/'//trim(FileName)//'_fgo_'// &
!!$       't'// &   
!!$       CHAR(mod(istep,10000000)/1000000+48)//&
!!$       CHAR(mod(istep,1000000)/100000+48)//&
!!$       CHAR(mod(istep,100000)/10000+48)//&
!!$       CHAR(mod(istep,10000)/1000+48)//&
!!$       CHAR(mod(istep,1000)/100+48)//&
!!$       CHAR(mod(istep,100)/10+48)//&
!!$       CHAR(mod(istep,10)+48)
!!$  CALL WriteFlux(fgo(:,:),Nx,Ny,Nstart,MEMO,trim(Rootname))
!!$
!!$  Rootname = trim(OUT_DIR)//'/'//trim(FileName)//'_fld_'// &
!!$       't'// &   
!!$       CHAR(mod(istep,10000000)/1000000+48)//&
!!$       CHAR(mod(istep,1000000)/100000+48)//&
!!$       CHAR(mod(istep,100000)/10000+48)//&
!!$       CHAR(mod(istep,10000)/1000+48)//&
!!$       CHAR(mod(istep,1000)/100+48)//&
!!$       CHAR(mod(istep,100)/10+48)//&
!!$       CHAR(mod(istep,10)+48)
!!$  CALL WriteFlux(fld(:,:),Nx,Ny,Nstart,MEMO,trim(Rootname))
!!$
!!$  Rootname = trim(OUT_DIR)//'/'//trim(FileName)//'_fgd_'// &
!!$       't'// &   
!!$       CHAR(mod(istep,10000000)/1000000+48)//&
!!$       CHAR(mod(istep,1000000)/100000+48)//&
!!$       CHAR(mod(istep,100000)/10000+48)//&
!!$       CHAR(mod(istep,10000)/1000+48)//&
!!$       CHAR(mod(istep,1000)/100+48)//&
!!$       CHAR(mod(istep,100)/10+48)//&
!!$       CHAR(mod(istep,10)+48)
!!$  CALL WriteFlux(fgd(:,:),Nx,Ny,Nstart,MEMO,trim(Rootname))

  





end Subroutine OutputData




!--------------------------------------------------------------
SUBROUTINE WriteDEN(DEN,myblock,LX,LY,Nstart,MEMO,FileName)
  !--------------------------------------------------------------

  IMPLICIT none
  INTEGER:: LX,LY,MEMO,RecLength,IO,idz,myblock,Nstart
  Double precision:: DEN(LX,LY,myblock)
  CHARACTER*(*),INTENT(IN):: FileName

  RecLength=MEMO*LX*LY
  OPEN(20,ACCESS='DIRECT',STATUS='UNKNOWN',RECL=RecLength, &
       FILE=FileName,IOSTAT=IO)
  DO idz=1,myblock
     WRITE(20,REC=idz)DEN(:,:,idz)
  ENDDO
  CLOSE(20)
  RETURN
END SUBROUTINE WriteDEN




!--------------------------------------------------------------
SUBROUTINE WriteFlux(DEN,LX,LY,Nstart,MEMO,FileName)
  !--------------------------------------------------------------

  IMPLICIT none
  INTEGER:: LX,LY,MEMO,RecLength,IO,idz,Nstart
  Double precision:: DEN(LX,LY)
  CHARACTER*(*),INTENT(IN):: FileName

  RecLength=MEMO*LX*LY
  OPEN(20,ACCESS='DIRECT',STATUS='UNKNOWN',RECL=RecLength, &
       FILE=FileName,IOSTAT=IO)
  WRITE(20,REC=1)DEN(:,:)
  CLOSE(20)
  RETURN
END SUBROUTINE WriteFlux






!--------------------------------------------------------------
SUBROUTINE WriteMatlab()
  !--------------------------------------------------------------
  USE GLOBAL
  IMPLICIT none
  INTEGER :: idx,time

  OPEN(20,FILE=trim(OUT_DIR)//'/'//'LoadGeo.m')
  WRITE(20,'(A,I4,A,I4,A,I4,A,A,A,A,A)')'por = ReadReal(' ,LX,',', LY, ',' ,LZ,&
       ',''../',trim(IN_DIR),'/',trim(PMname),''',''float64'');'
  CLOSE(20)



  OPEN(20,FILE=trim(OUT_DIR)//'/'//'LoadCon.m')

  WRITE(20,'(A,I8,A)')'DT = ', File_step,';'

  WRITE(20,'(A,I8,A,I8,A,I8,A,I8,A,A,A)')'c = zeros( ',Max_step/File_step+1,&
       ',', LX,',', LY,',', LZ,');'
  WRITE(20,'(A,I8,A,I8,A,I8,A,I8,A,A,A)')'lc = zeros( ',Max_step/File_step+1,&
       ',', LX,',', LY,',', LZ,');'
  WRITE(20,'(A,I8,A,I8,A,I8,A,I8,A,A,A)')'sc = zeros( ',Max_step/File_step+1,&
       ',', LX,',', LY,',', LZ,');'
  WRITE(20,'(A,I8,A,I8,A,I8,A,I8,A,A,A)')'b = zeros( ',Max_step/File_step+1,&
       ',', LX,',', LY,',', LZ,');'
  WRITE(20,'(A,I8,A,I8,A,I8,A,I8,A,A,A)')'to = zeros( ',Max_step/File_step+1,&
       ',', LX,',', LY,',', LZ,');'
  WRITE(20,'(A,I8,A,I8,A,I8,A,I8,A,A,A)')'lo = zeros( ',Max_step/File_step+1,&
       ',', LX,',', LY,',', LZ,');'
  WRITE(20,'(A,I8,A,I8,A,I8,A,I8,A,A,A)')'go = zeros( ',Max_step/File_step+1,&
       ',', LX,',', LY,',', LZ,');'
  WRITE(20,'(A,I8,A,I8,A,I8,A,I8,A,A,A)')'td = zeros( ',Max_step/File_step+1,&
       ',', LX,',', LY,',', LZ,');'
  WRITE(20,'(A,I8,A,I8,A,I8,A,I8,A,A,A)')'ld = zeros( ',Max_step/File_step+1,&
       ',', LX,',', LY,',', LZ,');'
  WRITE(20,'(A,I8,A,I8,A,I8,A,I8,A,A,A)')'gd = zeros( ',Max_step/File_step+1,&
       ',', LX,',', LY,',', LZ,');'


  Rootname = '../'//trim(IN_DIR)//'/'//trim(PMname)//'_c'
  WRITE(20,'(A,I8,A,I8,A,I8,A,I8,A,A,A)')'c(', 1,',:,:,:) = ReadReal( ',&
       LX,',', LY,',', LZ,',''', trim(Rootname),''',''float64'');'
  Rootname = '../'//trim(IN_DIR)//'/'//trim(PMname)//'_lc'
  WRITE(20,'(A,I8,A,I8,A,I8,A,I8,A,A,A)')'lc(', 1,',:,:,:) = ReadReal( ',&
       LX,',', LY,',', LZ,',''', trim(Rootname),''',''float64'');'
  Rootname = '../'//trim(IN_DIR)//'/'//trim(PMname)//'_sc'
  WRITE(20,'(A,I8,A,I8,A,I8,A,I8,A,A,A)')'sc(', 1,',:,:,:) = ReadReal( ',&
       LX,',', LY,',', LZ,',''', trim(Rootname),''',''float64'');'
  Rootname = '../'//trim(IN_DIR)//'/'//trim(PMname)//'_b'
  WRITE(20,'(A,I8,A,I8,A,I8,A,I8,A,A,A)')'b(', 1,',:,:,:) = ReadReal( ',&
       LX,',', LY,',', LZ,',''', trim(Rootname),''',''float64'');'
  Rootname = '../'//trim(IN_DIR)//'/'//trim(PMname)//'_to'
  WRITE(20,'(A,I8,A,I8,A,I8,A,I8,A,A,A)')'o(', 1,',:,:,:) = ReadReal( ',&
       LX,',', LY,',', LZ,',''', trim(Rootname),''',''float64'');'
  Rootname = '../'//trim(IN_DIR)//'/'//trim(PMname)//'_lo'
  WRITE(20,'(A,I8,A,I8,A,I8,A,I8,A,A,A)')'lo(', 1,',:,:,:) = ReadReal( ',&
       LX,',', LY,',', LZ,',''', trim(Rootname),''',''float64'');'
  Rootname = '../'//trim(IN_DIR)//'/'//trim(PMname)//'_go'
  WRITE(20,'(A,I8,A,I8,A,I8,A,I8,A,A,A)')'go(', 1,',:,:,:) = ReadReal( ',&
       LX,',', LY,',', LZ,',''', trim(Rootname),''',''float64'');'
  Rootname = '../'//trim(IN_DIR)//'/'//trim(PMname)//'_td'
  WRITE(20,'(A,I8,A,I8,A,I8,A,I8,A,A,A)')'d(', 1,',:,:,:) = ReadReal( ',&
       LX,',', LY,',', LZ,',''', trim(Rootname),''',''float64'');'
  Rootname = '../'//trim(IN_DIR)//'/'//trim(PMname)//'_ld'
  WRITE(20,'(A,I8,A,I8,A,I8,A,I8,A,A,A)')'ld(', 1,',:,:,:) = ReadReal( ',&
       LX,',', LY,',', LZ,',''', trim(Rootname),''',''float64'');'
  Rootname = '../'//trim(IN_DIR)//'/'//trim(PMname)//'_gd'
  WRITE(20,'(A,I8,A,I8,A,I8,A,I8,A,A,A)')'gd(', 1,',:,:,:) = ReadReal( ',&
       LX,',', LY,',', LZ,',''', trim(Rootname),''',''float64'');'





  DO idx = 1, Max_step/File_step

     time=idx*File_step
     Rootname = trim(PMname)//'_c_'// &
          't'// &   
          CHAR(mod(time,10000000)/1000000+48)//&
          CHAR(mod(time,1000000)/100000+48)//&
          CHAR(mod(time,100000)/10000+48)//&
          CHAR(mod(time,10000)/1000+48)//&
          CHAR(mod(time,1000)/100+48)//&
          CHAR(mod(time,100)/10+48)//&
          CHAR(mod(time,10)+48)
     WRITE(20,'(A,I8,A,I8,A,I8,A,I8,A,A,A)')'c(', idx+1,',:,:,:) = ReadReal( ',&
          LX,',', LY,',', LZ,',''', trim(Rootname),''',''float64'');'

     time=idx*File_step
     Rootname = trim(PMname)//'_lc_'// &
          't'// &   
          CHAR(mod(time,10000000)/1000000+48)//&
          CHAR(mod(time,1000000)/100000+48)//&
          CHAR(mod(time,100000)/10000+48)//&
          CHAR(mod(time,10000)/1000+48)//&
          CHAR(mod(time,1000)/100+48)//&
          CHAR(mod(time,100)/10+48)//&
          CHAR(mod(time,10)+48)
     WRITE(20,'(A,I8,A,I8,A,I8,A,I8,A,A,A)')'lc(', idx+1,',:,:,:) = ReadReal( ',&
          LX,',', LY,',', LZ,',''', trim(Rootname),''',''float64'');'

     time=idx*File_step
     Rootname = trim(PMname)//'_sc_'// &
          't'// &   
          CHAR(mod(time,10000000)/1000000+48)//&
          CHAR(mod(time,1000000)/100000+48)//&
          CHAR(mod(time,100000)/10000+48)//&
          CHAR(mod(time,10000)/1000+48)//&
          CHAR(mod(time,1000)/100+48)//&
          CHAR(mod(time,100)/10+48)//&
          CHAR(mod(time,10)+48)
     WRITE(20,'(A,I8,A,I8,A,I8,A,I8,A,A,A)')'sc(', idx+1,',:,:,:) = ReadReal( ',&
          LX,',', LY,',', LZ,',''', trim(Rootname),''',''float64'');'


     time=idx*File_step
     Rootname = trim(PMname)//'_b_'// &
          't'// &   
          CHAR(mod(time,10000000)/1000000+48)//&
          CHAR(mod(time,1000000)/100000+48)//&
          CHAR(mod(time,100000)/10000+48)//&
          CHAR(mod(time,10000)/1000+48)//&
          CHAR(mod(time,1000)/100+48)//&
          CHAR(mod(time,100)/10+48)//&
          CHAR(mod(time,10)+48)
     WRITE(20,'(A,I8,A,I8,A,I8,A,I8,A,A,A)')'b(', idx+1,',:,:,:) = ReadReal( ',&
          LX,',', LY,',', LZ,',''', trim(Rootname),''',''float64'');'


     time=idx*File_step
     Rootname = trim(PMname)//'_to_'// &
          't'// &   
          CHAR(mod(time,10000000)/1000000+48)//&
          CHAR(mod(time,1000000)/100000+48)//&
          CHAR(mod(time,100000)/10000+48)//&
          CHAR(mod(time,10000)/1000+48)//&
          CHAR(mod(time,1000)/100+48)//&
          CHAR(mod(time,100)/10+48)//&
          CHAR(mod(time,10)+48)
     WRITE(20,'(A,I8,A,I8,A,I8,A,I8,A,A,A)')'o(', idx+1,',:,:,:) = ReadReal( ',&
          LX,',', LY,',', LZ,',''', trim(Rootname),''',''float64'');'



     time=idx*File_step
     Rootname = trim(PMname)//'_lo_'// &
          't'// &   
          CHAR(mod(time,10000000)/1000000+48)//&
          CHAR(mod(time,1000000)/100000+48)//&
          CHAR(mod(time,100000)/10000+48)//&
          CHAR(mod(time,10000)/1000+48)//&
          CHAR(mod(time,1000)/100+48)//&
          CHAR(mod(time,100)/10+48)//&
          CHAR(mod(time,10)+48)
     WRITE(20,'(A,I8,A,I8,A,I8,A,I8,A,A,A)')'lo(', idx+1,',:,:,:) = ReadReal( ',&
          LX,',', LY,',', LZ,',''', trim(Rootname),''',''float64'');'


     time=idx*File_step
     Rootname = trim(PMname)//'_go_'// &
          't'// &   
          CHAR(mod(time,10000000)/1000000+48)//&
          CHAR(mod(time,1000000)/100000+48)//&
          CHAR(mod(time,100000)/10000+48)//&
          CHAR(mod(time,10000)/1000+48)//&
          CHAR(mod(time,1000)/100+48)//&
          CHAR(mod(time,100)/10+48)//&
          CHAR(mod(time,10)+48)
     WRITE(20,'(A,I8,A,I8,A,I8,A,I8,A,A,A)')'go(', idx+1,',:,:,:) = ReadReal( ',&
          LX,',', LY,',', LZ,',''', trim(Rootname),''',''float64'');'



     time=idx*File_step
     Rootname = trim(PMname)//'_td_'// &
          't'// &   
          CHAR(mod(time,10000000)/1000000+48)//&
          CHAR(mod(time,1000000)/100000+48)//&
          CHAR(mod(time,100000)/10000+48)//&
          CHAR(mod(time,10000)/1000+48)//&
          CHAR(mod(time,1000)/100+48)//&
          CHAR(mod(time,100)/10+48)//&
          CHAR(mod(time,10)+48)
     WRITE(20,'(A,I8,A,I8,A,I8,A,I8,A,A,A)')'d(', idx+1,',:,:,:) = ReadReal( ',&
          LX,',', LY,',', LZ,',''', trim(Rootname),''',''float64'');'



     time=idx*File_step
     Rootname = trim(PMname)//'_ld_'// &
          't'// &   
          CHAR(mod(time,10000000)/1000000+48)//&
          CHAR(mod(time,1000000)/100000+48)//&
          CHAR(mod(time,100000)/10000+48)//&
          CHAR(mod(time,10000)/1000+48)//&
          CHAR(mod(time,1000)/100+48)//&
          CHAR(mod(time,100)/10+48)//&
          CHAR(mod(time,10)+48)
     WRITE(20,'(A,I8,A,I8,A,I8,A,I8,A,A,A)')'ld(', idx+1,',:,:,:) = ReadReal( ',&
          LX,',', LY,',', LZ,',''', trim(Rootname),''',''float64'');'


     time=idx*File_step
     Rootname = trim(PMname)//'_gd_'// &
          't'// &   
          CHAR(mod(time,10000000)/1000000+48)//&
          CHAR(mod(time,1000000)/100000+48)//&
          CHAR(mod(time,100000)/10000+48)//&
          CHAR(mod(time,10000)/1000+48)//&
          CHAR(mod(time,1000)/100+48)//&
          CHAR(mod(time,100)/10+48)//&
          CHAR(mod(time,10)+48)
     WRITE(20,'(A,I8,A,I8,A,I8,A,I8,A,A,A)')'gd(', idx+1,',:,:,:) = ReadReal( ',&
          LX,',', LY,',', LZ,',''', trim(Rootname),''',''float64'');'




  enddo


  close(20)






  OPEN(20,FILE=trim(OUT_DIR)//'/'//'LoadFlux.m')

  WRITE(20,'(A,I8,A)')'DT = ', File_step,';'

  WRITE(20,'(A,I8,A,I8,A,I8,A,I8,A)')'flo = zeros( ',Max_step/File_step+1,&
       ',', LX,',', LY,');'
  WRITE(20,'(A,I8,A,I8,A,I8,A,I8,A)')'fgo = zeros( ',Max_step/File_step+1,&
       ',', LX,',', LY,');'
  WRITE(20,'(A,I8,A,I8,A,I8,A,I8,A)')'fld = zeros( ',Max_step/File_step+1,&
       ',', LX,',', LY,');'
  WRITE(20,'(A,I8,A,I8,A,I8,A,I8,A)')'fgd = zeros( ',Max_step/File_step+1,&
       ',', LX,',', LY,');'



!!$  Rootname = '../'//trim(IN_DIR)//'/'//trim(PMname)//'_flo'
!!$  WRITE(20,'(A,I8,A,I8,A,I8,A,A,A)')'flo(', 1,',:,:) = ReadFlux( ',&
!!$       LX,',', LY,',''', trim(Rootname),''',''float64'');'
!!$  Rootname = '../'//trim(IN_DIR)//'/'//trim(PMname)//'_fgo'
!!$  WRITE(20,'(A,I8,A,I8,A,I8,A,A,A)')'fgo(', 1,',:,:) = ReadFlux( ',&
!!$       LX,',', LY,',''', trim(Rootname),''',''float64'');'
!!$  Rootname = '../'//trim(IN_DIR)//'/'//trim(PMname)//'_fld'
!!$  WRITE(20,'(A,I8,A,I8,A,I8,A,A,A)')'fld(', 1,',:,:) = ReadFlux( ',&
!!$       LX,',', LY,',''', trim(Rootname),''',''float64'');'
!!$  Rootname = '../'//trim(IN_DIR)//'/'//trim(PMname)//'_fgd'
!!$  WRITE(20,'(A,I8,A,I8,A,I8,A,A,A)')'fgd(', 1,',:,:) = ReadFlux( ',&
!!$       LX,',', LY,',''', trim(Rootname),''',''float64'');'





  DO idx = 1, Max_step/File_step


     time=idx*File_step
     Rootname = trim(PMname)//'_flo_'// &
          't'// &   
          CHAR(mod(time,10000000)/1000000+48)//&
          CHAR(mod(time,1000000)/100000+48)//&
          CHAR(mod(time,100000)/10000+48)//&
          CHAR(mod(time,10000)/1000+48)//&
          CHAR(mod(time,1000)/100+48)//&
          CHAR(mod(time,100)/10+48)//&
          CHAR(mod(time,10)+48)
     WRITE(20,'(A,I8,A,I8,A,I8,A,A,A)')'flo(', idx+1,',:,:) = ReadFlux( ',&
          LX,',', LY,',''', trim(Rootname),''',''float64'');'

     time=idx*File_step
     Rootname = trim(PMname)//'_fgo_'// &
          't'// &   
          CHAR(mod(time,10000000)/1000000+48)//&
          CHAR(mod(time,1000000)/100000+48)//&
          CHAR(mod(time,100000)/10000+48)//&
          CHAR(mod(time,10000)/1000+48)//&
          CHAR(mod(time,1000)/100+48)//&
          CHAR(mod(time,100)/10+48)//&
          CHAR(mod(time,10)+48)
     WRITE(20,'(A,I8,A,I8,A,I8,A,A,A)')'fgo(', idx+1,',:,:) = ReadFlux( ',&
          LX,',', LY,',''', trim(Rootname),''',''float64'');'

     time=idx*File_step
     Rootname = trim(PMname)//'_fld_'// &
          't'// &   
          CHAR(mod(time,10000000)/1000000+48)//&
          CHAR(mod(time,1000000)/100000+48)//&
          CHAR(mod(time,100000)/10000+48)//&
          CHAR(mod(time,10000)/1000+48)//&
          CHAR(mod(time,1000)/100+48)//&
          CHAR(mod(time,100)/10+48)//&
          CHAR(mod(time,10)+48)
     WRITE(20,'(A,I8,A,I8,A,I8,A,A,A)')'fld(', idx+1,',:,:) = ReadFlux( ',&
          LX,',', LY,',''', trim(Rootname),''',''float64'');'

     time=idx*File_step
     Rootname = trim(PMname)//'_fgd_'// &
          't'// &   
          CHAR(mod(time,10000000)/1000000+48)//&
          CHAR(mod(time,1000000)/100000+48)//&
          CHAR(mod(time,100000)/10000+48)//&
          CHAR(mod(time,10000)/1000+48)//&
          CHAR(mod(time,1000)/100+48)//&
          CHAR(mod(time,100)/10+48)//&
          CHAR(mod(time,10)+48)
     WRITE(20,'(A,I8,A,I8,A,I8,A,A,A)')'fgd(', idx+1,',:,:) = ReadFlux( ',&
          LX,',', LY,',''', trim(Rootname),''',''float64'');'




  enddo


  close(20)







  RETURN
END SUBROUTINE WriteMatlab








!--------------------------------------------------------------
SUBROUTINE SOR(ma,mb,mc,md,me,mf,mg,HH,QQ,omega,HH_star,Nx,Ny,Nz,Nt,HH_eps,WW)
  !--------------------------------------------------------------    

  !  Use GLOBAL
  Implicit None
  Integer::pp,nn,Nx,Ny,Nz,Nt
  Double precision::omega,HH_err,HH_eps
  Double precision::SMALL_NUM
  Double precision,Dimension(Nt)::ma,mb,mc,md,me,mf,mg,HH,QQ,HH_star
  Logical,Dimension(Nt)::WW


  nn = 1
  SMALL_NUM = 1e-15

  Do While(nn.lt.1e5)



     if(WW(1)) then
        HH(1) = 0.
     else
        HH(1) = omega*(QQ(1)-mc(1)*HH_star(2)-me(1)*HH_star(1+Nx)-mg(1)*HH_star(1+Nx*Ny))/ &
             (ma(1)+SMALL_NUM) + (1-omega)*HH_star(1)
     endif
     Do pp=2,Nx
        if(WW(pp)) then
           HH(pp) = 0.
        else
           HH(pp) = omega*(QQ(pp)-mb(pp)*HH(pp-1)-mc(pp)*HH_star(pp+1)-me(pp)*HH_star(pp+Nx)- &
                mg(pp)*HH_star(pp+Nx*Ny))/(ma(pp)+SMALL_NUM) + (1-omega)*HH_star(pp)
        endif
     Enddo
     Do pp=1+Nx,Nx*Ny
        if(WW(pp)) then
           HH(pp) = 0.
        else
           HH(pp) = omega*(QQ(pp)-md(pp)*HH(pp-Nx)-mb(pp)*HH(pp-1)-mc(pp)*HH_star(pp+1)- &
                me(pp)*HH_star(pp+Nx)-mg(pp)*HH_star(pp+Nx*Ny))/(ma(pp)+SMALL_NUM) + &
                (1-omega)*HH_star(pp)
        endif
     Enddo
     Do pp=1+Nx*Ny,Nt-Nx*Ny
        if(WW(pp)) then
           HH(pp) = 0.
        else
           HH(pp) = omega*(QQ(pp)-mf(pp)*HH(pp-Nx*Ny)-md(pp)*HH(pp-Nx)- &
                mb(pp)*HH(pp-1)-mc(pp)*HH_star(pp+1)-me(pp)*HH_star(pp+Nx)- &
                mg(pp)*HH_star(pp+Nx*Ny))/(ma(pp)+SMALL_NUM) + (1-omega)*HH_star(pp)
        endif
     Enddo
     Do pp=Nt-Nx*Ny+1,Nt-Nx
        if(WW(pp)) then
           HH(pp) = 0.
        else
           HH(pp) = omega*(QQ(pp)-mf(pp)*HH(pp-Nx*Ny)-md(pp)*HH(pp-Nx)- &
                mb(pp)*HH(pp-1)-mc(pp)*HH_star(pp+1)-me(pp)*HH_star(pp+Nx))/(ma(pp)+SMALL_NUM) + &
                (1-omega)*HH_star(pp)
        endif
     Enddo
     Do pp=Nt-Nx+1,Nt-1
        if(WW(pp)) then
           HH(pp) = 0.
        else
           HH(pp) = omega*(QQ(pp)-mf(pp)*HH(pp-Nx*Ny)-md(pp)*HH(pp-Nx)- &
                mb(pp)*HH(pp-1)-mc(pp)*HH_star(pp+1))/(ma(pp)+SMALL_NUM) + (1-omega)*HH_star(pp)
        endif
     Enddo
     if(WW(Nt)) then
        HH(Nt) = 0.
     else
        HH(Nt) = omega*(QQ(Nt)-mf(Nt)*HH(Nt-Nx*Ny)-md(Nt)*HH(Nt-Nx)- &
             mb(Nt)*HH(Nt-1))/(ma(Nt)+SMALL_NUM) +(1-omega)*HH_star(Nt)
     endif




     HH_err = Sum(Abs(HH-HH_star))/Nt

!!$
!!$     print*,'inner iteration, nn=',nn,'HH_err=',HH_err
!!$     print*,'HH_star=',maxval(HH_star),minval(HH_star)
!!$     print*,'HH=',maxval(HH),minval(HH)

!!$     do pp=1,Nt
!!$        print*,pp+121,HH(pp)
!!$     enddo
!!$ 
     If (HH_err.lt.HH_eps) exit

     HH_star = HH
     nn = nn + 1

  Enddo
  RETURN
END SUBROUTINE SOR




