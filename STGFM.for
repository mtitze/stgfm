c--------------------------------------------------------------------------------------
c
c     Symplectic tracking through an extended multipole
c     the mathematics have been derived by Malte Titze
c       Note of January 9th, 2014
c
c     The equations are different as compared the earlier versions 
c       where the formalism did not describe the displacement correctly.
c
c     ATTENTION: Due to increased complexity of the equations the key equations 
c       are programmed in complex notation, now, in contrast to the earlier code
c
c     input data in cartesian coordinates
c        x0, y0     positions of particle in mm
c        xp0, yp0   angles of particle in mrad
c
c--------------------------------------------------------------------------------------

      include 'STGFM.cmn'
      real*8 x0ll,y0ll,xp0ll,yp0ll
      real*8 xfll,yfll,xpfll,ypfll
      real*8 bx,by,bz,r_loc,phi_loc,z_loc,xx,xlint0
      integer isteps,i_new_fields,ispecial
      
      dimension z0s(10000),xpfs(10000),xpfs0(10000)              
      dimension x0l(10001),xp0l(10001),y0l(10001),yp0l(10001)
      dimension xfl(0:10001),xpfl(0:10001),yfl(0:10001),ypfl(0:10001)

c------------- initialize

      write(6,*) '                                '
      write(6,*) ' *****************************  '
      write(6,*) '  STGFM                         '
      write(6,*) ' *****************************  '
      write(6,*) ' Experimental Code - 06.04.2015 '
      write(6,*) ' (c) by M. Titze and J. Bahrdt  '
      write(6,*) '                                '

      Pi=4.d0*datan(1.d0)
      small=1.d0-10
      iverdat=0    ! set this to one just before any dndnh001-call to
                   ! write down the corresponding data at this point
      ireverse=1   ! 1: change sign of all fields
      ilinear=0     ! 1: X002=0, Y002=0
                    ! 10: X002=0, Y002=0 & X011=X101=Y011=Y101=0
      i_only_c0=0   ! 1: use only c0 Fourier coefficent 
     
      ifldzero=0    ! 1: c0 is turned off 
      isinoff=1     ! 1: no sin-terms in longitudinal expansion
      i_old=0       ! 1: old code
      iverbose=1  
      bz_zf_off=0    ! 1: force B_z=0 at zf, 
                     ! 2: correct potentials by constant A(x0,y0,zf)
                     !    this corrects observed
                     !    asymmetry problem for quadrupoles and is
                     !    correct (see theory)
                     ! else: do nothing.

      iwbnbnh=1      ! Write bn and bnh to file.
      i_new_fields=1 ! The relations between the fields hold only
                     ! for infinite order in p. If the summation
                     ! is stopped at some point, we would introduce
                     ! an error. Set this to 1 to avoid this error.

      ispecial=0     ! 1: larger lattice

      call get_input
      if(ireverse.eq.1)scal=-scal

      call get_xkxn
      if(i_readjust.eq.1)call readjust_iinterval 
c------------- z0 and zf are set after this point

      call get_xkxn    
      call cn_from_cntsnt_new   
      call scal_cn      ! include stiffness of beam
      call get_Gm0  
      call apnapnh            
      call facexp  
       
c      call get_psi
      
c------------- putting numbers      
c     imode   = 0: single point
c     imode   = 1: x, xp
c     imode   = 2: y, yp
c     imode   = 3: x, y
c     imode   = 4: xp, yp
c     imode   = 5: starting with x0, y0, tracking many turns

1000  imode=5
c      imode=1234

      write(6,*)' imode, igf mode  = ',imode,igf
      write(6,*)'-------------------------'      

c      if(imode.eq.1234)then
c      do imode=1,4
          
      if(imode.eq.0)then        
      x00=0.d-10       ! mm     avoid division by zero
      xp00=0.d0      ! mrad
      y00=0.d-10       ! mm
      yp00=0.d0      ! mrad        
      write(6,*)'   '
c      write(6,*)' xp0,yp0 (mrad) '
c      read(5,*)xp00,yp00      
c      write(6,*)' x00,y00 (mm) '
c      read(5,*)x00,y00
      write(6,*)' x0,xp0 (mm) '
      read(5,*)x00,xp00      
c      write(6,*)' y0,yp0 (mm) '
c      read(5,*)y00,yp00      
      iloop = 1
      dphi=0.d0
      goto 8888
      endif     
      
      if(imode.eq.1)then        
      x00=9.d0       ! mm
      xp00=10.d0     ! mrad
      y00=0.d0       ! mm
      yp00=0.d0      ! mrad      
      endif
      
      if(imode.eq.2)then        
      x00=0.d0       ! mm
      xp00=0.d0      ! mrad
      y00=9.d0       ! mm
      yp00=10.d0     ! mrad      
      endif
      
      if(imode.eq.3)then        
      x00=9.d0       ! mm
      xp00=0.d0      ! mrad
      y00=9.d0       ! mm
      yp00=0.d0      ! mrad      
      endif
             
      if(imode.eq.4)then        
      x00=1.d-5       ! mm
      xp00=10.d0     ! mrad
      y00=1.d-5       ! mm
      yp00=10.d0     ! mrad      
      endif
      
      if(imode.eq.5)then        
      x00=30.0d0       ! mm
      xp00=0.d0      ! mrad
      y00=0.d0      ! mm
      yp00=0.d0     ! mrad
      endif
          
      iloop=63
      dphi=(2.d0*Pi)/dfloat(iloop-1)
      if(imode.eq.5)iloop=1000
      
8888  continue
    
      if(iloop.gt.1)iverbose=0

      if(imode.eq.5)then
      xfl(0)=x00
      xpfl(0)=xp00
      yfl(0)=y00
      ypfl(0)=yp00      
      endif
      
      z00=z0
      zf0=zf
      xlint0=xlint
      
c---------------------------------------------------------------------
      Do ii=1,iloop
c---------------------------------------------------------------------
      phi=(ii-1)*dphi

      if(imode.eq .0)then
      x0l(ii)=x00
      xp0l(ii)=xp00
      y0l(ii)=y00
      yp0l(ii)=yp00
      endif
      
      if((imode.eq.1).or.(imode.eq.2))then
      x0ll=x00*dcos(phi)
      xp0ll=xp00*dsin(phi)
      y0ll=y00*dcos(phi)
      yp0ll=yp00*dsin(phi)     
      endif
      
      if((imode.eq.3).or.(imode.eq.4))then
      x0ll=x00*dcos(phi)
      xp0ll=xp00*dcos(phi)
      y0ll=y00*dsin(phi)
      yp0ll=yp00*dsin(phi)      
      endif 

      if(imode.eq.5)then
      x0a=xfl(ii-1)
      y0a=yfl(ii-1)
      xp0a=xpfl(ii-1)
      yp0a=ypfl(ii-1)  
    
c----- first half of SR
      call track_ring(x0a,xp0a,y0a,yp0a,
     &                x0ll,xp0ll,y0ll,yp0ll,2)

        if(ispecial.eq.1)then
        x0a=x0ll
        xp0a=xp0ll
        y0a=y0ll
        yp0a=yp0ll      
        call track_ring(x0a,xp0a,y0a,yp0a,
     &                x0ll,xp0ll,y0ll,yp0ll,1)   
        x0a=x0ll
        xp0a=xp0ll
        y0a=y0ll
        yp0a=yp0ll           
        call track_ring(x0a,xp0a,y0a,yp0a,
     &                x0ll,xp0ll,y0ll,yp0ll,2)  
        endif

      endif
      
      if(imode.ne.0)then
       x0l(ii)=x0ll
       y0l(ii)=y0ll
       xp0l(ii)=xp0ll
       yp0l(ii)=yp0ll
      endif
      
       x0=x0l(ii)/1000.d0     ! get m
       xp0=xp0l(ii)/1000.d0   ! get m
       y0=y0l(ii)/1000.d0     ! get rad
       yp0=yp0l(ii)/1000.d0   ! get rad
       
c======================= start 3D-multipole ========================
      xf=x0
      xpf=xp0      
      yf=y0
      ypf=yp0
      
      do isteps=1,itwo_steps
      x0=xf
      y0=yf
      xp0=xpf
      yp0=ypf
          
      z0=z00
      zf=zf0   
      
      if(itwo_steps.eq.2)then  
      xlint=xlint0/2.d0
      if(isteps.eq.1)then
       z0=z00
       zf=0.d0
      else
       z0=0.d0
       zf=zf0
      endif
      endif

      call facexp

      write(6,*)'  '
      write(6,*)' start multipole segment ',isteps
      write(6,*)'-------------------------'     
      write(6,*)' z0, zf ',z0,zf
      write(6,*)' xlint ',xlint
      write(6,*)' x0, y0 (mm) ',(1000.d0)*x0,(1000.d0)*y0
      write(6,*)' xp0, yp0 (mrad) ',(1000.d0)*xp0,(1000.d0)*yp0
     
c ---- initialisiere mit erweiterter Ordnung

      call zerobsdbs               
      call zerodsdds        
      call zeroesdes
           
      call cartesian2polar(x0,y0,r0,phi0) 
      if (iverbose.eq.1) then
        write(6,*)' r0 (mm)', (1000.d0)*r0
        write(6,*)' phi0', phi0
      endif
c-------- Felder und 1., 2., 3. Ableitung 
c--------- erste Ordnung

      call bnbnh(r0,phi0)

      if(i_new_fields.eq.1)then
        call drrrbnbnh_JB(r0,phi0)
        call dbnbnh_JB(r0,phi0)
        call ddbnbnh_JB(r0,phi0)
        call dddbnbnh_JB(r0,phi0)
      else
        call dbnbnh(r0,phi0)
        call ddbnbnh(r0,phi0)
        call dddbnbnh(r0,phi0)
      endif 

      call dndnh001(r0,phi0)

      if (bz_zf_off.eq.2) then
        call get_AxAy(phi0,zf,axrp0zf,ayrp0zf)
          if (iverbose.eq.1) write(6,*)' axrp0zf, ayrp0zf', 
     & axrp0zf, ayrp0zf
        dn001(0)=dn001(0) - dsin(phi0)*axrp0zf + dcos(phi0)*ayrp0zf
        dnh001(0)=dnh001(0) + dcos(phi0)*axrp0zf + dsin(phi0)*ayrp0zf
      endif

      if((iverbose.eq.1).and.(iverdat.eq.1))then
        iverdat=0
        write(6,*)' writing dn001(n) and dnh001(n) to dndnh.dat'
        open(unit=10,file='dndnh.dat')
        do i=-iord1,iord1
        write(10,*)dn001(i),dnh001(i)   
        enddo 
        close(10)
      endif

      call ddndnh001  
      call dddndnh001
      call ddddndnh001
      
c--------- zweite Ordnung    

      call eta2(r0,phi0)    ! Verwendung von dn001 und dnh001
      call deta2(r0,phi0)   ! Verwendung von ddn001 und ddnh001
      call ddeta2(r0,phi0)  ! Verwendung von dddn001 und dddnh001
      call dddeta2(r0,phi0) ! Verwendung von ddddn001 und ddddnh001
      
      call dndnh2(r0,phi0)  ! 1. Ableitungen von eta011 eta101 eta002
      call ddndnh2(r0,phi0) ! 1.-3. Ableitung von eta011 eta101 eta002
      
c--------- dritte Ordnung

      call eta3(r0,phi0)    ! Verwendung von dn101, dn011 dnh101 dnh011
      call deta3(r0,phi0)   ! Verwendung von ddn101 ddn011 ddn002 dddn002

      call dndnh3(r0,phi0)  ! Verwendung von deta201 deta021 deta111 
                            ! deta102 deta012
      call ddndnh3(r0,phi0)  

c--------- vierte Ordnung

      call eta4(r0,phi0)
      call deta4(r0,phi0)  
      
c-----------------------------------------
      call get_fqpk

      call gx(r0,phi0)  ! gross x
      call gy(r0,phi0)  ! gross y
      call get_AxAy(phi0,z0,Ax0,Ay0)

c      write(6,*) ' Ax0 ==== ', Ax0          
      call get_px0py0(xp0,yp0,Ax0,Ay0,px0,py0)  
      if(i_old.ne.1)then
        call get_pxfpyf_newton(px0,py0,pxf,pyf)  
        else
        call get_pxfpyf_old(Ax0,Ay0,xp0,yp0,px0,py0,pxf,pyf)  
      endif
        write(6,*)' px0, py0 ',px0,py0  
        write(6,*)' pxf, pyf ',pxf,pyf       
      if(i_old.ne.1)then
        call get_xfyf(x0,y0,pxf,pyf,xf,yf) 
        else
        call get_xfyf_old(x0,y0,xp0,yp0,Ax0,Ay0,pxf,pyf,xf,yf)
      endif
      call cartesian2polar(xf,yf,rf,phif)
      if (iverbose.eq.1) then
        write(6,*)' rf (mm)', (1000.d0)*rf
        write(6,*)' phif', phif
      endif

      call zerobsdbs
      call zerodsdds

      call bnbnh(rf,phif)

      call dndnh001(rf,phif)
      if (bz_zf_off.eq.2) then
        dn001(0)=dn001(0) - dsin(phi0)*axrp0zf + dcos(phi0)*ayrp0zf
        dnh001(0)=dnh001(0) + dcos(phi0)*axrp0zf + dsin(phi0)*ayrp0zf
      endif

      call get_AxAy(phif,zf,Axf,Ayf)
      
      xpf = pxf - Axf
      ypf = pyf - Ayf
      
      write(6,*)'  '
      write(6,*)' end multipole segment ',isteps
      write(6,*)'-------------------------'     
      write(6,*)' x0, y0 (mm) ',(1000.d0)*x0,(1000.d0)*y0
      write(6,*)' xf, yf (mm) ',(1000.d0)*xf,(1000.d0)*yf
      write(6,*)' xp0, yp0 (mrad) ',(1000.d0)*xp0,(1000.d0)*yp0
      write(6,*)' xpf, ypf (mrad) ',(1000.d0)*xpf,(1000.d0)*ypf  

      write(6,*)' Ax0, Ay0 ',Ax0,Ay0      
      write(6,*)' Axf, Ayf ',Axf,Ayf   

      enddo ! itwo_steps
      
      write(6,*)'  '  
      write(6,*)' end multipole '
      write(6,*)'-------------------------'     
      write(6,*)' x0, y0 (mm) ',(1000.d0)*x0,(1000.d0)*y0
      write(6,*)' xf, yf (mm) ',(1000.d0)*xf,(1000.d0)*yf
      write(6,*)' xp0, yp0 (mrad) ',(1000.d0)*xp0,(1000.d0)*yp0
      write(6,*)' xpf, ypf (mrad) ',(1000.d0)*xpf,(1000.d0)*ypf

      write(6,*)' Ax0, Ay0 ',Ax0,Ay0      
      write(6,*)' Axf, Ayf ',Axf,Ayf
     
c================== end 3D-multipole =====================================

c----- convert to mm and mrad
      xfa=1000.d0*xf
      yfa=1000.d0*yf
      xpfa=1000.d0*xpf
      ypfa=1000.d0*ypf
     
      if(imode.eq.5)then
c----- second half of SR      
      call track_ring(xfa,xpfa,yfa,ypfa,
     &                xfll,xpfll,yfll,ypfll,2)

        if(ispecial.eq.1)then
        x0a=xfll
        xp0a=xpfll
        y0a=yfll
        yp0a=ypfll      
        call track_ring(x0a,xp0a,y0a,yp0a,
     &                xfll,xpfll,yfll,ypfll,1)   
        x0a=xfll
        xp0a=xpfll
        y0a=yfll
        yp0a=ypfll           
        call track_ring(x0a,xp0a,y0a,yp0a,
     &                xfll,xpfll,yfll,ypfll,2)
        endif 
       
      else
      xfll=xfa
      yfll=yfa
      xpfll=xpfa
      ypfll=ypfa
      endif
      
c----- save data in mm mrad
      xfl(ii)=xfll
      yfl(ii)=yfll
      xpfl(ii)=xpfll
      ypfl(ii)=ypfll
      
      enddo       ! end loop ===========================================
      
c     imode   = 0: single point
c     imode   = 1: x, xp
c     imode   = 2: y, yp
c     imode   = 3: x, y
c     imode   = 4: xp, yp
1     format(4e15.6)
2     format(6e15.6)      
      
      if(imode.eq.0)goto 900
      if(imode.eq.1)goto 100
      if(imode.eq.2)goto 200     
      if(imode.eq.3)goto 300
      if(imode.eq.4)goto 400   
      if(imode.eq.5)goto 500        
      
100   open(unit=10,file='PSXGF.dat')
      do ii=1,iloop
      write(10,1)x0l(ii),xp0l(ii),xfl(ii),xpfl(ii)    
      enddo    
      close(10)
      goto 900
      
200   open(unit=10,file='PSYGF.dat')
      do ii=1,iloop
      write(10,1)y0l(ii),yp0l(ii),yfl(ii),ypfl(ii)    
      enddo    
      close(10)
      goto 900
      
300   open(unit=10,file='PSXYGF.dat')
      do ii=1,iloop
      write(10,2)x0l(ii),y0l(ii),xfl(ii),xpfl(ii),yfl(ii),ypfl(ii)     
      enddo    
      close(10)
      goto 900
      
400   open(unit=10,file='PSXPYPGF.dat')
      do ii=1,iloop
      write(10,2)xp0l(ii),yp0l(ii),xfl(ii),xpfl(ii),yfl(ii),ypfl(ii)
      enddo    
      close(10)
      goto 900      

500   open(unit=10,file='symplecticity.dat')
      do ii=1,iloop
      write(10,2)y0l(ii),yp0l(ii),xfl(ii),xpfl(ii),yfl(ii),ypfl(ii)
      enddo    
      close(10)
      
      open(unit=10,file='tune-x-STGFM.dat')
      do ii=1,iloop
      write(10,111)dfloat(ii),xfl(ii)
      enddo    
      close(10) 
      
      open(unit=10,file='tune-y-STGFM.dat')
      do ii=1,iloop
      write(10,111)dfloat(ii),yfl(ii)
      enddo    
      close(10)   
            
111   format(2d20.12)      
      open(unit=10,file='tune-ref.dat')
      do ii=1,iloop
      xx=100.d0*((2.d0*pi)/dfloat(iloop))
      write(10,111)dfloat(ii),dsin(xx*dfloat(ii))
      enddo    
      close(10)     
      
      goto 900            
            
900         continue
      
c      r_loc=rf
c      phi_loc=phif
c      z_loc=zf
c      write(6,*)' r_loc, z_loc ',r_loc,z_loc
c      call bnbnh(r_loc,phi_loc)
c      call dndnh001(r_loc,phi_loc)
c      call ddndnh001 
c      call get_b(bx,by,bz,r_loc,z_loc)
c      call get_AxAy(rf,phif,zf,Axf,Ayf)     
      
c      write(6,*)' bz = ',bz
c      write(6,*)' zf-z0 = ',zf-z0   
c      write(6,*)' Ax0, Ay0 ',Ax0,Ay0      
c      write(6,*)' Axf, Ayf ',Axf,Ayf
      
      stop
      end

c--------------------------------
      subroutine track_ring(x0a,xp0a,y0a,yp0a,x0b,xp0b,y0b,yp0b,iii)
c--------------------------------        
c
c     Uebergaben in mm und mrad
c     innerhalb der Routine m und rad (und Tesla)
c
      include 'STGFM.cmn' 
      integer iii
      real*8 matGesF(2,2),matGesD(2,2),matFk(2,2),matDk(2,2),
     & matL1(2,2),matL2(2,2),
     & qbb,qdx,br,ql,dl,wk,xk,cc,ss,ch,sh,ffact

      dist1=100.d0
      dist2=100.d0

c umrechnen in m und rad 
      x0a=x0a/1000.d0
      xp0a=xp0a/1000.d0
      y0a=y0a/1000.d0
      yp0a=yp0a/1000.d0      
      dist1=dist1/1000.d0
      dist2=dist2/1000.d0

      xk=(2.2307d0/2.d0)*atrack
      if(ireverse.eq.1)xk=-xk   
      
      matFk(1,1)=1.d0
      matFk(1,2)=0.d0
      matFk(2,1)=-xk
      matFk(2,2)=1.d0
      
      matDk(1,1)=1.d0
      matDk(1,2)=0.d0
      matDk(2,1)=xk
      matDk(2,2)=1.d0
      
      matL1(1,1)=1.d0
      matL1(1,2)=dist1
      matL1(2,1)=0.d0
      matL1(2,2)=1.d0
      
      matL2(1,1)=1.d0
      matL2(1,2)=dist2
      matL2(2,1)=0.d0
      matL2(2,2)=1.d0
      
      if(iii.eq.1)then
      x0a1=matL1(1,1)*x0a+matL1(1,2)*xp0a
      xp0a1=matL1(2,1)*x0a+matL1(2,2)*xp0a
      y0a1=matL1(1,1)*y0a+matL1(1,2)*yp0a
      yp0a1=matL1(2,1)*y0a+matL1(2,2)*yp0a
      
      x0a2=matDk(1,1)*x0a1+matDk(1,2)*xp0a1
      xp0a2=matDk(2,1)*x0a1+matDk(2,2)*xp0a1
      y0a2=matFk(1,1)*y0a1+matFk(1,2)*yp0a1
      yp0a2=matFk(2,1)*y0a1+matFk(2,2)*yp0a1 
   
      x0b=matL2(1,1)*x0a2+matL2(1,2)*xp0a2
      xp0b=matL2(2,1)*x0a2+matL2(2,2)*xp0a2
      y0b=matL2(1,1)*y0a2+matL2(1,2)*yp0a2
      yp0b=matL2(2,1)*y0a2+matL2(2,2)*yp0a2   
      endif
 
      if(iii.eq.2)then
      x0a1=matL2(1,1)*x0a+matL2(1,2)*xp0a
      xp0a1=matL2(2,1)*x0a+matL2(2,2)*xp0a
      y0a1=matL2(1,1)*y0a+matL2(1,2)*yp0a
      yp0a1=matL2(2,1)*y0a+matL2(2,2)*yp0a
      
      x0a2=matFk(1,1)*x0a1+matFk(1,2)*xp0a1
      xp0a2=matFk(2,1)*x0a1+matFk(2,2)*xp0a1
      y0a2=matDk(1,1)*y0a1+matDk(1,2)*yp0a1
      yp0a2=matDk(2,1)*y0a1+matDk(2,2)*yp0a1 
      
      x0b=matL1(1,1)*x0a2+matL1(1,2)*xp0a2
      xp0b=matL1(2,1)*x0a2+matL1(2,2)*xp0a2
      y0b=matL1(1,1)*y0a2+matL1(1,2)*yp0a2
      yp0b=matL1(2,1)*y0a2+matL1(2,2)*yp0a2        
      endif
           
      x0a=x0a*1000.d0
      xp0a=xp0a*1000.d0
      y0a=y0a*1000.d0
      yp0a=yp0a*1000.d0      
      
      x0b=x0b*1000.d0
      xp0b=xp0b*1000.d0
      y0b=y0b*1000.d0
      yp0b=yp0b*1000.d0           

      dist1=dist1*1000.d0
      dist2=dist2*1000.d0
      
      return
      end
      
c--------------------------------
      subroutine get_psi
c--------------------------------        
      include 'STGFM.cmn'
      complex*16 psi(0:10), psi_zero(0:10)
      integer j,k,nmax
      
c------ j=1
      j=1      ! j =/= 0. if j=0 we have psi_zero
      nmax=1
      
      psi(0)=(1.d0/(xi*xkxn(j)))*
     &  (cdexp(-xi*xkxn(j)*s1)-cdexp(-xi*xkxn(j)*s2))
      psi_zero(0)=s2-s1
      do k=1,nmax
        psi(k)=xi*(cdexp(-(xi*xkxn(j)*zf))/xkxn(j))*
     &  ((zf-s2)**k*cdexp(xi*xkxn(j)*(zf-s2))-
     &  (zf-s1)**k*cdexp(xi*xkxn(j)*(zf-s1)))+
     &  xi*(dfloat(k)/xkxn(j))*psi(k-1)

	psi_zero(k)=1.d0/(dfloat(k+1))*((zf-s1)**(k+1)-
     &   (zf-s2)**(k+1))
      enddo
      
1     format(' k,j,psi(k,j) ',2I4,2d20.10)
2     format(' k,psi_zero(k) ',I4,2d20.10)     
c      do k=1,nmax
c          write(6,1)k,j,dreal(psi(k)),dimag(psi(k))
c      enddo
c      do k=1,nmax
c          write(6,2)k,dreal(psi_zero(k)),dimag(psi_zero(k))
c      enddo
c      
c      write(6,*)'-------------------------'
      
c-------- j=-1
      j=-1      ! j =/= 0. if j=0 we have psi_zero
      
      psi(0)=(1.d0/(xi*xkxn(j)))*
     &  (cdexp(-xi*xkxn(j)*s1)-cdexp(-xi*xkxn(j)*s2))
      psi_zero(0)=s2-s1
      do k=1,nmax
        psi(k)=xi*(cdexp(-(xi*xkxn(j)*zf))/xkxn(j))*
     &  ((zf-s2)**k*cdexp(xi*xkxn(j)*(zf-s2))-
     &  (zf-s1)**k*cdexp(xi*xkxn(j)*(zf-s1)))+
     &  xi*(dfloat(k)/xkxn(j))*psi(k-1)

	psi_zero(k)=1.d0/(dfloat(k+1))*((zf-s1)**(k+1)-
     &   (zf-s2)**(k+1))
      enddo
      
c      do k=1,nmax
c          write(6,1)k,j,dreal(psi(k)),dimag(psi(k))
c      enddo
c      
c      write(6,*)'-------------------------'
c      
c      write(6,*)' dppeta102(0) ',dppeta102(0)
c      write(6,*)' dppeta102(0) ',dppeta102(0)
c      write(6,*)' dppeta102(1) ',dppeta102(1)
c      write(6,*)'j, drdnh002(j) ',j,drdnh002(j)
c      write(6,*)'j, drdn002(j) ',j,drdn002(j)    

c      write(6,*)' eta102(1) ',eta102(1) 
      
c      write(6,*)' dnh002(1) ',dnh002(1)
c      write(6,*)' dpdn002(1) ',dpdn002(1)      
c      write(6,*)' dppdnh002(1) ',dppdnh002(1)     
      
c      write(6,*)' eta002(1) ',eta002(1)
      
c      write(6,*)' dnh101(1) ',dnh101(1)     
c      write(6,*)' dphieta101(1) ',dphieta101(1)
      
c      write(6,*)' dphidnh001(1) ',dphidnh001(1)
      
c      write(6,*)' dppdnh101(1) ',dppdnh101(1)
c      write(6,*)' dppreta101(1) ',dppreta101(1)      

c      write(6,*)' eta101(0) ',eta101(0)
c      write(6,*)' dreta101(0) ',dreta101(0)
c      write(6,*)' dppeta101(0) ',dppeta101(0)
c      write(6,*)' dppreta101(0) ',dppreta101(0)
      
c      write(6,*)' dpprdnh001(0) ',dpprdnh001(0)

c      write(6,*)' drdnh102(0) ',drdnh102(0)
c      write(6,*)' drreta102(0) ',drreta102(0)
c      write(6,*)' drrreta002(0) ',drrreta002(0)

c      write(6,*)' drrrdnh001(0) ',drrrdnh001(0)
c      write(6,*)' drrrdn001(0) ',drrrdn001(0)

c      write(6,*)' drrrbnh(0) ',drrrbnh(0)

      return
      end
      
c--------------------------------
      subroutine get_b(bx,by,bz,r,z)
c--------------------------------        
      include 'STGFM.cmn'
      complex*16 psi(0:10), psi_zero(0:10)
      real*8 bx,by,bz
      integer j,k,nmax

      bx=0.d0
      by=0.d0
      bz=0.d0
      
      do nfou=-nfour,nfour
      bz=bz+
     & (-drdn001(nfou)+(1.d0/r)*dphidnh001(nfou)-(1.d0/r)*dn001(nfou))*
     & cdexp(xi*xkxn(nfou)*z)
      
c      write(6,*)' get_b ',drdn001(nfou)+(1.d0/r)*dphidnh001(nfou)
c      write(6,*)' get_b1 ',dn001(nfou)     
      
      enddo
  
      return
      end
      
c--------------------------------
      subroutine drift(x1,y1,xp1,yp1)
c--------------------------------        
c     in dieser Routine mm und mrad
c
      include 'STGFM.cmn'
      real*8 x1,y1,xp1,yp1,drft,xx,yy,xxp,yyp
      
      drft=-0.5d0*(zf-z0)   ! in m
      write(6,*)' drift = ',drft
      write(6,*)' xx0 ',x1
      xx=x1+drft*xp1    
      yy=y1+drft*yp1      
      xxp=xp1
      yyp=yp1

      x1=xx
      y1=yy
      xp1=xxp
      yp1=yyp
      
      write(6,*)' xx1 ',x1      
      
      return
      end

c--------------------------------
      subroutine get_input     
c-------------------------------- 
      include 'STGFM.cmn'
      integer clen1      

      xm1=-1.d0
      xi=cdsqrt(xm1)

      write(6,*)' reading input data from STGFM.par'
1     format(a80)       
      open(unit=10,file='STGFM.par')
      read(10,*)energy        ! energy in GeV
        emass=9.10938291d-31
        clight=2.99792458d8  
        echarge=1.60217657d-19      
        gamma=(1.0d9*energy*echarge)/(emass*clight**2) 
        beta=dsqrt(1.d0-1.d0/gamma*2)
        speed=beta*clight
        scal=(gamma*emass*speed)/echarge    
        write(6,*)' Brho = ',scal
        
      read(10,*)iord1    ! expansion order of exponentials
      read(10,*)mo      ! multipole order (1=dipole, 2=quadrupole etc)
        xmo=mo
      read(10,*)z0,zf   ! integration interval
      read(10,*)i_readjust
      read(10,*)s1,s2   ! expansion interval for Fourier coefficients
        xlint=zf-z0     ! length of multipole including fringe fields
        xl=s2-s1        ! 1st harmonic period for 
                        ! Fourier expansion (can be different from xl (FF)) 
      read(10,*)mp0     ! maximum p-value for transverse field expansion 
      read(10,1)dfile   ! file name of Fourier coefficients
        ifile=clen1(dfile)
        fcfile(1:ifile)=dfile(1:ifile)
        write(6,*)' reading Fourier coefficients from'
        write(6,1)fcfile  
      read(10,*)r00     ! reference radius for Fourier decomposition
      read(10,*)nfour   ! number of Fourier coef. for longitudinal exp.
                        ! not including the zeroth term
      read(10,*)igf	! determine which generating functions are taken
c      write(6,*)' nfour input = ',nfour
      read(10,*)atrack  ! additional factor for all elements in imode = 5   
      read(10,*)itwo_steps
      close(10)

c------- the energy order is max 3 in this program.
      iord_max=iord1    ! 4*iord1  order for all 
                        ! terms who do not require z0, zf
      iord_newton=iord1  ! order used for Newton-Routine
      iord2=iord1   ! 2*iord1  order used in computation of 
                    !  terms wrt. energy order 2
      iord3=iord1   ! 3*iord1    "             3

      do nfou=0,iord_max
        cnt(nfou)=0.d0
        snt(nfou)=0.d0
      enddo
      
      open(unit=10,file=fcfile,status='old')
      do nfou=0,nfour
      read(10,*)nn,cnt(nfou),snt(nfou)
      enddo 
      close(10)
      
c------- the following operation compnesates an error in FOUR_NEW      
      cnt(0)=cnt(0)/2.d0     ! snt(0)=0
      if(ifldzero.eq.1)cnt(0)=0.d0     ! snt(0)=0            

c------- atrack-scaling 
      do nfou=0,nfour
      cnt(nfou)=cnt(nfou)*atrack
      snt(nfou)=snt(nfou)*atrack
      enddo

c------------ del sin-terms for test purposes
      if(isinoff.eq.1)then
      do nfou=0,nfour
      snt(nfou)=0.d0
      enddo       
      endif

      return      
      end   

c------------------------------------        
      subroutine get_xkxn     
c------------------------------------
      include 'STGFM.cmn'

      xl=s2-s1
      xpi=4.d0*datan(1.d0)
      xkx0=(2.d0*xpi)/xl

      xkxn(0)=0.d0
      xkxn2(0)=0.d0
      do nfou=1,iord_max
        xkxn(nfou)=xkx0*dfloat(nfou)
        xkxn2(nfou)=xkxn(nfou)*xkxn(nfou)
        xkxn(-nfou)=-xkxn(nfou)
        xkxn2(-nfou)=xkxn2(nfou)
      enddo 

      return
      end

c------------------------------------        
      subroutine readjust_iinterval     
c------------------------------------      
      include 'STGFM.cmn'

      small_loc=0.05d0*((s2-s1)/dfloat(iord_newton))
      dd=(s2-s1)/dfloat(4*iord_newton)
      
c----- ajust left boundary z0
      z00=z0
      za=z0    
100   call fdf(za,func,dfunc)
      if(dabs(dfunc).lt.small_loc)then
        za=za-dd
        else
        delta=-func/dfunc
        za=za+delta
      endif
c     write(6,*)' left ',delta
      if(dabs(delta).gt.1.d-15)goto 100
      z0=za
      write(6,*)'-------------------------'
      write(6,*)' new integration boundaries'
      write(6,*)' old / new z0 = ',z00,z0
      
c----- ajust right boundary zf
      zf0=zf
      za=zf
      
200   call fdf(za,func,dfunc)
      if(dabs(dfunc).lt.small_loc)then
        za=za+dd
        else
        delta=-func/dfunc
        za=za+delta
      endif
c     write(6,*)' right ',delta
      if(dabs(delta).gt.1.d-15)goto 200
      zf=za
      write(6,*)' old / new zf = ',zf0,zf
      if(iverbose.eq.1)then
      write(6,*)' readjust_iinterval: z represented up to order', 
     & iord_newton
      endif
      write(6,*)'-------------------------'      
      xlint=zf-z0


      return
      end

c------------------------------------        
      subroutine fdf(za,func,dfunc)     
c------------------------------------      
      include 'STGFM.cmn'      
      
      func=-za
      dfunc=-1.d0

      do nfou=1,iord_newton
      func=func+(2.d0/xkxn(nfou))*(
     &  dsin(xkxn(nfou)*s1)*dcos(xkxn(nfou)*za)-
     &  dcos(xkxn(nfou)*s1)*dsin(xkxn(nfou)*za))
      dfunc=dfunc-2.d0*(
     &  dsin(xkxn(nfou)*s1)*dsin(xkxn(nfou)*za)+
     &  dcos(xkxn(nfou)*s1)*dcos(xkxn(nfou)*za))
      enddo
 
      return
      end

c------------------------------------        
      subroutine cn_from_cntsnt     
c------------------------------------      
      include 'STGFM.cmn'      
      dimension deno(0:20)
      r02=r00*r00

      deno(0)=2.d0*(r00**(mo-1)) ! Achtung: m! kürzt sich weg  
      cnr(0)=cnt(0)/(deno(0)/2.d0)
      cni(0)=0.d0
      cn(0)=cnr(0)+xi*cni(0)
      
      do ip=1,mp0
        deno(ip)=deno(ip-1)*(r02/dfloat(4*(mo+ip)*ip))*xkxn2(ip)      
      enddo
      
      denom=0.d0
      do ip=0,mp0
        denom=denom+deno(ip)
      enddo

      do nfou=1,nfour        
        cnr(nfou)=cnt(nfou)/denom
        cni(nfou)=-snt(nfou)/denom
        cn(nfou)=cnr(nfou)+xi*cni(nfou)
        cn(-nfou)=dconjg(cn(nfou))               
      enddo
      
      do nfou=-nfour,nfour
      if((i_only_c0.eq.1).and.(nfou.ne.0))cn(nfou)=0.d0
      enddo      

      if(iverbose.eq.1)write(6,*)'writing CN*r00 to CN.DAT ',nfour  
      open(unit=10,file='CN.DAT',status='unknown')    
      do nfou=-nfour,nfour
      write(10,*)nfou,r00*cn(nfou)
      enddo      
      close(10)
      
      return      
      end         

c------------------------------------        
      subroutine cn_from_cntsnt_new     
c------------------------------------      
      include 'STGFM.cmn'  
      dimension deno(0:20)

      r02=r00*r00
      deno(0)=xmo*(r00**(mo-1))

      cnr(0)=cnt(0)/(deno(0)/dfloat(ifacult(mo)))
      cni(0)=0.d0

      cn(0)=cnr(0)+xi*cni(0)
c initialize cn's with zeros for n!=0, 300 is hard-coded in .cmn
      do nfou=1,iord_max
        cn(-nfou)=0.d0
        cn(nfou)=0.d0
      enddo

      do ip=1,mp0
        deno(ip)=deno(ip-1)*(r02/dfloat(4*(mo+ip)*ip))      
      enddo
      if(i_only_c0.ne.1)then
c overwrite cn's from -nfour to nfour
        do nfou=1,nfour
        denom=0.d0
        do ip=0,mp0
          denom=denom+deno(ip)*xkxn2(nfou)**ip
        enddo
          cnr(nfou)=cnt(nfou)/denom
          cni(nfou)=-snt(nfou)/denom
          cn(nfou)=cnr(nfou)+xi*cni(nfou)
          cn(-nfou)=dconjg(cn(nfou))          
        enddo
      endif  

      if(iverbose.eq.1)write(6,*)'writing CN*r00 to CN.DAT ',nfour
      open(unit=10,file='CN.DAT',status='unknown')    
      do nfou=-nfour,nfour
      write(10,*)nfou,r00*cn(nfou)
      enddo      
      close(10)
      
      return      
      end       

c------------------------------------        
      subroutine cn_from_cntsnt_new_JB     
c------------------------------------      
      include 'STGFM.cmn'  
      dimension deno(0:20)
      r02=r00*r00
      
      deno(0)=xmo*(r00**(mo-1))
      cnr(0)=cnt(0)/(deno(0)/dfloat(ifacult(mo)))

      cni(0)=0.d0
      cn(0)=cnr(0)+xi*cni(0)

      do ip=1,mp0
        deno(ip)=deno(ip-1)*(r02/dfloat(4*(mo+ip)*ip))      
      enddo

      do nfou=1,nfour 
      denom=0.d0
      do ip=0,mp0
        denom=denom+deno(ip)*xkxn2(nfou)**ip
      enddo
        cnr(nfou)=cnt(nfou)/denom
        cni(nfou)=-snt(nfou)/denom
        cn(nfou)=cnr(nfou)+xi*cni(nfou)
        cn(-nfou)=dconjg(cn(nfou))          
      enddo
      
      do nfou=-nfour,nfour
      if((i_only_c0.eq.1).and.(nfou.ne.0))cn(nfou)=0.d0
      enddo      
  

      if(iverbose.eq.1)write(6,*)'writing CN*r00 to CN.DAT ',nfour  
      open(unit=10,file='CN.DAT',status='unknown')    
      do nfou=-nfour,nfour
      write(10,*)nfou,r00*cn(nfou)
      enddo      
      close(10)
      
      return      
      end    

c--------------------------------        
      subroutine get_Gm0     
c--------------------------------      
      include 'STGFM.cmn'    
     
      ianzGm0=401
      zmin=z0
      dz=xlint/dfloat(ianzGm0-1)
      
      do i=1,ianzGm0
      zz(i)=zmin+dfloat(i-1)*dz
      do nfou=-nfour,nfour
      Gm0(i)=Gm0(i)+cn(nfou)*cdexp(xi*xkxn(nfou)*zz(i))              
      enddo
      
      enddo
      
      if(iverbose.eq.1)then
      open(unit=10,file='GM0.DAT',status='replace')
      write(6,*)'writing Gm0*r00 to GM0.DAT'
      write(6,*)'-------------------------'
      do i=1,ianzGm0
      write(10,*)zz(i),r00*dreal(Gm0(i)),r00*dimag(Gm0(i))
      enddo
      close(10)   
      endif
      
      return
      end

c------------------------------------------        
      subroutine scal_cn      
c------------------------------------------      

      include 'STGFM.cmn'   

      do nfou=-iord_max,iord_max
        cn(nfou)=cn(nfou)/scal     
      enddo

      return 
      end   
          
c------------------------------------------        
      subroutine only_c0      
c------------------------------------------      
      include 'STGFM.cmn'   
      
      do nfou=1,iord_max
        cn(nfou)=0.d0
        cn(-nfou)=0.d0
      enddo

      return 
      end   
          
c------------------------------------------        
      subroutine apnapnh      
c------------------------------------------      
c
c     define apn(ip,nfou) & apnh(ip,nfou)
c           for nfou >/= 1
c
c------------------------------------------
      include 'STGFM.cmn'      

c---- define ap & aph; indices from 0 to mp0
      ap(0)=1.0d0/dfloat(ifacult(mo-1))
      aph(0)=ap(0)   
      do ip=1,mp0          
        ap(ip)=1.0d0/dfloat(4**ip*ifacult(mo+ip)*ifacult(ip))          
        aph(ip)=ap(ip)
      enddo
 
      do ip=1,mp0
        ap(ip)=ap(ip)*(mo+2.d0*ip)
        aph(ip)=aph(ip)*mo
      enddo

c---- define apn & apnh           
      do nfou=1,iord_max
        fac=1.d0/xkxn2(nfou)
        do ip=0,mp0
          fac=fac*xkxn2(nfou)       
          apn(ip,nfou)=ap(ip)*fac
          apnh(ip,nfou)=aph(ip)*fac      
        enddo     
      
      enddo

1     format(' ip, ap, aph ',i4,2f15.6)
2     format(' ip, apn, apnh ',i4,2f15.6)
      
c      write(6,*)'-------- nfou = 1 --------------------'
c      do ip=0,mp0
c      write(6,1)ip,ap(ip),aph(ip)
c      enddo
c     
c      do ip=0,mp0
c      write(6,2)ip,apn(ip,1),apnh(ip,1)
c      enddo
c
c      write(6,*)'--------------------------------------'     
           
      return      
      end  
      
c--------------------------------        
      subroutine zerobsdbs      
c--------------------------------      
      include 'STGFM.cmn'      

      do nfou=-iord_max,iord_max
      bn(nfou)=0.d0
      bnh(nfou)=0.d0
      
      drbn(nfou)=0.d0
      drbnh(nfou)=0.d0
      dphibn(nfou)=0.d0
      dphibnh(nfou)=0.d0
      
      dppbn(nfou)=0.d0
      dprbn(nfou)=0.d0
      drrbn(nfou)=0.d0     
      dppbnh(nfou)=0.d0
      dprbnh(nfou)=0.d0
      drrbnh(nfou)=0.d0 
      
      dpppbn(nfou)=0.d0
      dpppbnh(nfou)=0.d0
      dpprbn(nfou)=0.d0
      dpprbnh(nfou)=0.d0
      dprrbn(nfou)=0.d0
      dprrbnh(nfou)=0.d0
      drrrbn(nfou)=0.d0
      drrrbnh(nfou)=0.d0  
      enddo
              
      return
      end

c--------------------------------        
      subroutine bnbnh(r,phi)      
c------------------------------------------------------------------------------   
c
c     Attention: sin(m phi), cos(m phi) terms are skipped at this point 
c     to avoid checking for division by zero
c     the terms will be regarded later
c
c------------------------------------------------------------------------------   
      include 'STGFM.cmn'
      complex*16 facc,si    
     
      r2=r*r
      
      bn(0) =(1.d0/dfloat(ifacult(mo-1)))*cn(0)*r**(mo-1)
      bnh(0)=(1.d0/dfloat(ifacult(mo-1)))*cn(0)*r**(mo-1)
         
      do nfou=1,nfour
      bn(nfou)=0.d0
      bnh(nfou)=0.d0
      facc=cn(nfou)
      re=r**(mo-3) 
      
      do ip=0,mp0
      re=re*r2   
      bn(nfou)=bn(nfou)+apn(ip,nfou)*re
      bnh(nfou)=bnh(nfou)+apnh(ip,nfou)*re 
      enddo
                  
      bn(nfou)=bn(nfou)*facc    
      bnh(nfou)=bnh(nfou)*facc     
      
      bn(-nfou)=dconjg(bn(nfou))    
      bnh(-nfou)=dconjg(bnh(nfou))           

      enddo
     
      if((iverbose.eq.1).and.(iwbnbnh.eq.1))then
      iwbnbnh=0
      write(6,*)'writing BN to BN.DAT and BNH to BNH.DAT',nfour
      open(unit=10,file='BN.DAT',status='replace')
      do nfou=-nfour,nfour
      write(10,*)nfou,bn(nfou)
      enddo
      close(10)   
      
      open(unit=10,file='BNH.DAT',status='replace')
      do nfou=-nfour,nfour
      write(10,*)nfou,bnh(nfou)
      enddo
      close(10)    
      endif
      
      return      
      end      

c--------------------------------        
      subroutine dbnbnh(r,phi)      
c-------------------------------- 
c
c     missing sin(m*phi), cos(m*phi) are implemented 
c
c---------------------------------     
      include 'STGFM.cmn'      

      do nfou=-nfour,nfour
      drbn(nfou)=(xmo/r)*dsin(xmo*phi)*
     &  (1.d0+(xkxn(nfou)*r)**2/xmo**2)*bnh(nfou)-
     &  (bn(nfou)*dsin(xmo*phi))/r
      dphibn(nfou)=xmo*dcos(xmo*phi)*bn(nfou)
      drbnh(nfou)=-(bnh(nfou)/r)*dcos(xmo*phi)+
     &  xmo*dcos(xmo*phi)*(bn(nfou)/r)
      dphibnh(nfou)=-xmo*dsin(xmo*phi)*bnh(nfou)               
      enddo
      
      return 
      end
c--------------------------------        
      subroutine ddbnbnh(r,phi)      
c-------------------------------- 
c
c     missing sin(m*phi), cos(m*phi) are implemented
c     avoid numeric simulations with big numbers
c
c---------------------------------     
      include 'STGFM.cmn'      
      complex*16 tmp
      real*8 small_loc1
 
      small_loc1=1.d-12
      
      do nfou=-nfour,nfour
      dppbn(nfou)=-dsin(xmo*phi)*xmo**2*bn(nfou)
      dprbn(nfou)=dcos(xmo*phi)*(1.d0+((xkxn(nfou)*r)/xmo)**2)*
     &     (xmo**2/r)*bnh(nfou)-((xmo*dcos(xmo*phi))/r)*bn(nfou)
      tmp=-((3.d0*xmo**2))*(dsin(xmo*phi)/xmo)*bnh(nfou)+
     &     dsin(xmo*phi)*((xmo**2+2.d0))*bn(nfou)   
      if(cdabs(tmp).lt.small_loc1)tmp=0.d0
      drrbn(nfou)=tmp/r**2-
     &    (xkxn(nfou)**2)*(dsin(xmo*phi)/xmo)*bnh(nfou)+
     &    dsin(xmo*phi)*(xkxn(nfou)**2)*bn(nfou)   
      
      dppbnh(nfou)=-dcos(xmo*phi)*xmo**2*bnh(nfou)
      dprbnh(nfou)=(xmo/r)*dsin(xmo*phi)*bnh(nfou)-
     &    dsin(xmo*phi)*(xmo**2/r)*bn(nfou)

      tmp=dcos(xmo*phi)*((xmo**2+2.d0))*bnh(nfou)-
     &    (3.d0*xmo*dcos(xmo*phi))*bn(nfou)
      if(cdabs(tmp).lt.small_loc1)tmp=0.d0
      drrbnh(nfou)=tmp/r**2+xkxn(nfou)**2*bnh(nfou)*dcos(xmo*phi)
      enddo
      
c      if(mp0.eq.0)then
c      do nfou=-nfour,nfour
c      drrbn(nfou)=0.d0
c      drrbnh(nfou)=0.d0
c      enddo
c      endif
      
      if(iverbose.eq.1)then
      write(6,*)'-----------------------------'      
      write(6,*)' bn(0) ',bn(0)
      write(6,*)' bnh(0) ',bnh(0)    
      write(6,*)'---'
      write(6,*)' bn(1) ',bn(1)
      write(6,*)' bnh(1) ',bnh(1)
      write(6,*)'---'
      write(6,*)' drrbn(1) ',drrbn(1)
      write(6,*)' drrbnh(1) ',drrbnh(1)       
      write(6,*)'-----------------------------'
      endif
      
      return 
      end

c--------------------------------        
      subroutine dddbnbnh(r,phi)      
c-------------------------------- 
c
c     missing sin(m*phi), cos(m*phi) are implenmeted 
c
c---------------------------------     
      include 'STGFM.cmn'      
      complex*16 tmp
      real*8 small_loc1

      small_loc1=1.d-12
      
      do nfou=-nfour,nfour
      dpppbnh(nfou)=-xmo**2*dphibnh(nfou)
      dpprbnh(nfou)=-xmo**2*drbnh(nfou)
      
      tmp=((xmo**2+2.d0))*dphibnh(nfou)+
     &  ((3.d0*xmo**2))*dsin(xmo*phi)*bn(nfou)
      if(cdabs(tmp).lt.small_loc1)tmp=0.d0     
      dprrbnh(nfou)=tmp/r**2+(xkxn2(nfou))*dphibnh(nfou)
      
      tmp=-(5.d0*xmo**2+4.d0)*bnh(nfou)*dcos(xmo*phi)+
     &     9.d0*xmo*dcos(xmo*phi)*bn(nfou)+
     &     (r*(xmo**2+2.d0))*drbnh(nfou)

      if(cdabs(tmp).lt.small_loc1)tmp=0.d0       
      drrrbnh(nfou)=tmp/r**3-(3.d0/r)*xkxn2(nfou)*bnh(nfou)*
     & dcos(xmo*phi)+(xkxn2(nfou))*drbnh(nfou)

      dpppbn(nfou)=-xmo**2*dphibn(nfou)
      dpprbn(nfou)=-xmo**2*drbn(nfou)
      
      tmp=-(3.d0*xmo**2)*dcos(xmo*phi)*bnh(nfou)+
     &  (xmo**2+2.d0)*dphibn(nfou)
      if(cdabs(tmp).lt.small_loc1)tmp=0.d0      
      dprrbn(nfou)=tmp/r**2-(xkxn2(nfou))*dcos(xmo*phi)*
     &  bnh(nfou)+(xkxn2(nfou))*dphibn(nfou)
      
      tmp=dsin(xmo*phi)*bnh(nfou)*
     &  ((9.d0*xmo))+dsin(xmo*phi)*bn(nfou)*
     &  (-(5.d0*xmo**2+4.d0))+      
     &  r*(xmo**2+2.d0)*drbn(nfou)
      
      if(cdabs(tmp).lt.small_loc1)tmp=0.d0  
      drrrbn(nfou)=tmp/r**3+dsin(xmo*phi)*bnh(nfou)*
     &  (xkxn2(nfou)/(xmo*r))+
     &  dsin(xmo*phi)*bn(nfou)*
     &  (-xkxn2(nfou)/r)+(xkxn2(nfou))*drbn(nfou)
      enddo
      
c      if(mp0.eq.0)then
c      do nfou=-nfour,nfour
c      drrrbn(nfou)=0.d0
c      drrrbnh(nfou)=0.d0
c      dprrbn(nfou)=0.d0
c      dprrbnh(nfou)=0.d0
c      enddo
c      endif
      
c      write(6,*)' kn**2(0) ',xkxn2(0)
c      write(6,*)' drrrbn(0),  drrrbnh(0)  ',drrrbn(0), drrrbnh(0) 
      
      return
      end
      
c------------------------------------------------------------------------------        
      subroutine drrrbnbnh_JB(r,phi)
c------------------------------------------------------------------------------   
c
c     sin(m phi), cos(m phi) are impemented at this point
c
c------------------------------------------------------------------------------   
      include 'STGFM.cmn'
      complex*16 facc,si,sx,cx,ff1,ffh1,ff2,ffh2,ff3,ffh3    
      real*8 rr
 
      sx=dsin(xmo*phi)     
      cx=dcos(xmo*phi) 
      
      if((mo-2).ge.0)then
        ff1=dfloat(mo-1)*sx
        ffh1=dfloat(mo-1)*cx      
      else
        ff1=0.d0
        ffh1=0.d0       
      endif
      
      drbn(0) =(ff1/dfloat(ifacult(mo-1)))*cn(0)*r**(mo-2)
      drbnh(0)=(ffh1/dfloat(ifacult(mo-1)))*cn(0)*r**(mo-2)
      
      if((mo-3).ge.0)then
        ff2=dfloat(mo-1)*dfloat(mo-2)*sx
        ffh2=dfloat(mo-1)*dfloat(mo-2)*cx        
      else
        ff2=0.d0
        ffh2=0.d0        
      endif
      
      drrbn(0) =(ff2/dfloat(ifacult(mo-1)))*cn(0)*r**(mo-3)
      drrbnh(0)=(ffh2/dfloat(ifacult(mo-1)))*cn(0)*r**(mo-3)
      
      if((mo-4).ge.0)then
        ff3=dfloat(mo-1)*dfloat(mo-2)*dfloat(mo-3)*sx
        ffh3=dfloat(mo-1)*dfloat(mo-2)*dfloat(mo-3)*cx       
      else
        ff3=0.d0
        ffh3=0.d0        
      endif
      
      drrrbn(0) =(ff3/dfloat(ifacult(mo-1)))*cn(0)*r**(mo-4)
      drrrbnh(0)=(ffh3/dfloat(ifacult(mo-1)))*cn(0)*r**(mo-4)
      
      do nfou=1,nfour
      drbn(nfou)=0.d0
      drbnh(nfou)=0.d0

      drrbn(nfou)=0.d0
      drrbnh(nfou)=0.d0
      
      drrrbn(nfou)=0.d0
      drrrbnh(nfou)=0.d0

      facc=cn(nfou)  
      
      do ip=0,mp0
      
      if((2*ip+mo-2).ge.0)then
        ff1=dfloat(2*ip+mo-1)*sx
        ffh1=dfloat(2*ip+mo-1)*cx      
      else
        ff1=0.d0
        ffh1=0.d0       
      endif     
      
      rr=r**(2*ip+mo-2)
      drbn(nfou)=drbn(nfou)+ff1*apn(ip,nfou)*rr
      drbnh(nfou)=drbnh(nfou)+ffh1*apnh(ip,nfou)*rr
      
      if((2*ip+mo-3).ge.0)then
        ff2=dfloat(2*ip+mo-1)*dfloat(2*ip+mo-2)*sx
        ffh2=dfloat(2*ip+mo-1)*dfloat(2*ip+mo-2)*cx      
      else
        ff2=0.d0
        ffh2=0.d0       
      endif           
          
      rr=rr/r
      drrbn(nfou)=drrbn(nfou)+ff2*apn(ip,nfou)*rr
      drrbnh(nfou)=drrbnh(nfou)+ffh2*apnh(ip,nfou)*rr       
 
      if((2*ip+mo-4).ge.0)then
        ff3=dfloat(2*ip+mo-1)*dfloat(2*ip+mo-2)*dfloat(2*ip+mo-3)*sx
        ffh3=dfloat(2*ip+mo-1)*dfloat(2*ip+mo-2)*dfloat(2*ip+mo-3)*cx
      else
        ff3=0.d0
        ffh3=0.d0       
      endif                
       
      rr=rr/r
      drrrbn(nfou)=drrrbn(nfou)+ff3*apn(ip,nfou)*rr
      drrrbnh(nfou)=drrrbnh(nfou)+ffh3*apnh(ip,nfou)*rr

      enddo
                  
      drbn(nfou)=drbn(nfou)*facc        
      drbnh(nfou)=drbnh(nfou)*facc          

      drrbn(nfou)=drrbn(nfou)*facc        
      drrbnh(nfou)=drrbnh(nfou)*facc 
      
      drrrbn(nfou)=drrrbn(nfou)*facc        
      drrrbnh(nfou)=drrrbnh(nfou)*facc                      

      drbn(-nfou)=dconjg(drbn(nfou))    
      drbnh(-nfou)=dconjg(drbnh(nfou))      
      
      drrbn(-nfou)=dconjg(drrbn(nfou))    
      drrbnh(-nfou)=dconjg(drrbnh(nfou))  

      drrrbn(-nfou)=dconjg(drrrbn(nfou))    
      drrrbnh(-nfou)=dconjg(drrrbnh(nfou))
      enddo
      
      return      
      end            

c--------------------------------        
      subroutine dbnbnh_JB(r,phi)      
c-------------------------------- 
c
c     missing sin(m*phi), cos(m*phi) are implenmeted 
c
c---------------------------------     
      include 'STGFM.cmn'      
 
      do nfou=-nfour,nfour
       
c      write(6,*)' dbn nfou ',nfou
c      write(6,*)' neu ',drbn(nfou)
c      drbn(nfou)=(xmo/r)*dsin(xmo*phi)*
c     &  (1.d0+(xkxn(nfou)*r)**2/xmo**2)*bnh(nfou)-
c     &  (bn(nfou)*dsin(xmo*phi))/r          
c      write(6,*)' alt ',drbn(nfou)
c      write(6,*)'-----'
      
      dphibn(nfou)=xmo*dcos(xmo*phi)*bn(nfou)
      
c      write(6,*)' dbnh nfou ',nfou
c      write(6,*)' neu ',drbnh(nfou)      

c      drbnh(nfou)=-(bnh(nfou)/r)*dcos(xmo*phi)+
c     &  xmo*dcos(xmo*phi)*(bn(nfou)/r)      
      
c      write(6,*)' alt ',drbnh(nfou)
c      write(6,*)'-----'
      
      dphibnh(nfou)=-xmo*dsin(xmo*phi)*bnh(nfou) 
      
      enddo
      
      return 
      end      

c--------------------------------        
      subroutine ddbnbnh_JB(r,phi)      
c-------------------------------- 
c
c     missing sin(m*phi), cos(m*phi) are implemented
c     avoid numeric simulations with big numbers
c
c---------------------------------     
      include 'STGFM.cmn'      
      complex*16 tmp,sx,cx
      real*8 small_loc1

      small_loc1=1.d-12
      
      sx=dsin(xmo*phi)     
      cx=dcos(xmo*phi)       
      
      do nfou=-nfour,nfour
      dppbn(nfou)=-dsin(xmo*phi)*xmo**2*bn(nfou)
      dppbnh(nfou)=-dcos(xmo*phi)*xmo**2*bnh(nfou)      
      
      if(abs(sx).gt.small_loc1)then
          dprbn(nfou)=drbn(nfou)*(cx/sx)*xmo
      else
          dprbn(nfou)=0.d0
      endif

      if(abs(cx).gt.small_loc1)then
          dprbnh(nfou)=drbnh(nfou)*(-sx/cx)*xmo
      else
          dprbnh(nfou)=0.d0
      endif      
      
      enddo
      
      return 
      end      

c--------------------------------        
      subroutine dddbnbnh_JB(r,phi)      
c-------------------------------- 
c
c     missing sin(m*phi), cos(m*phi) are implenmeted 
c
c---------------------------------     
      include 'STGFM.cmn'      
      complex*16 tmp,sx,cx
      real*8 small_loc1

      small_loc1=1.d-12
      
      sx=dsin(xmo*phi)     
      cx=dcos(xmo*phi)       
      
      do nfou=-nfour,nfour
      dpppbnh(nfou)=-xmo**2*dphibnh(nfou)
      dpppbn(nfou)=-xmo**2*dphibn(nfou)      
      
      dpprbn(nfou)=-xmo**2*drbn(nfou)
      dpprbnh(nfou)=-xmo**2*drbnh(nfou)
  
      if(abs(sx).gt.small_loc1)then
        dprrbn(nfou)=xmo*drrbn(nfou)*(cx/sx)
      else
        dprrbn(nfou)=0.d0
      endif
      
      if(abs(cx).gt.small_loc1)then      
        dprrbnh(nfou)=xmo*drrbnh(nfou)*(-sx/cx)      
      else
        dprrbnh(nfou)=0.d0
      endif
      
      enddo
      
      return
      end

c--------------------------------        
      subroutine zerodsdds      
c--------------------------------      
      include 'STGFM.cmn'      

      do nfou=-iord_max,iord_max
      dn001(nfou)=0.d0
      dnh001(nfou)=0.d0
      dn011(nfou)=0.d0
      dnh011(nfou)=0.d0
      dn101(nfou)=0.d0
      dnh101(nfou)=0.d0
      
      drdn001(nfou)=0.d0
      drdnh001(nfou)=0.d0
      dphidn001(nfou)=0.d0
      dphidnh001(nfou)=0.d0  
      drdn101(nfou)=0.d0
      drdnh101(nfou)=0.d0 
      drdn011(nfou)=0.d0
      drdnh011(nfou)=0.d0      
      
      dpdn101(nfou)=0.d0
      dpdnh101(nfou)=0.d0   
      dpdnh011(nfou)=0.d0
      dpdn011(nfou)=0.d0     
      
      drrdn001(nfou)=0.d0
      drrdnh001(nfou)=0.d0
      dprdn001(nfou)=0.d0
      dprdnh001(nfou)=0.d0  
      dppdn001(nfou)=0.d0
      dppdnh001(nfou)=0.d0 

      dn002(nfou)=0.d0
      dnh002(nfou)=0.d0     
      drdn002(nfou)=0.d0
      drdnh002(nfou)=0.d0    
      dpdn002(nfou)=0.d0
      dpdnh002(nfou)=0.d0         

      dppdnh011(nfou)=0.d0
      dppdnh101(nfou)=0.d0
      dppdnh002(nfou)=0.d0
      dppdn011(nfou)=0.d0
      dppdn101(nfou)=0.d0
      dppdn002(nfou)=0.d0
      dprdnh011(nfou)=0.d0
      dprdnh101(nfou)=0.d0
      dprdnh002(nfou)=0.d0
      dprdn011(nfou)=0.d0
      dprdn101(nfou)=0.d0
      dprdn002(nfou)=0.d0
      drrdnh011(nfou)=0.d0
      drrdnh101(nfou)=0.d0
      drrdnh002(nfou)=0.d0
      drrdn011(nfou)=0.d0
      drrdn101(nfou)=0.d0
      drrdn002(nfou)=0.d0
      dpppdn001(nfou)=0.d0
      dpppdnh001(nfou)=0.d0
      dpprdn001(nfou)=0.d0
      dpprdnh001(nfou)=0.d0
      dprrdn001(nfou)=0.d0
      dprrdnh001(nfou)=0.d0
      drrrdn001(nfou)=0.d0
      drrrdnh001(nfou)=0.d0    

      dnh201(nfou)=0.d0
      dnh021(nfou)=0.d0
      dnh111(nfou)=0.d0
      dnh102(nfou)=0.d0
      dnh012(nfou)=0.d0
      dn201(nfou)=0.d0
      dn021(nfou)=0.d0
      dn111(nfou)=0.d0
      dn102(nfou)=0.d0
      dn012(nfou)=0.d0
      drdnh201(nfou)=0.d0
      drdnh021(nfou)=0.d0
      drdnh111(nfou)=0.d0
      drdnh102(nfou)=0.d0
      drdnh012(nfou)=0.d0
      drdn201(nfou)=0.d0
      drdn021(nfou)=0.d0
      drdn111(nfou)=0.d0
      drdn102(nfou)=0.d0
      drdn012(nfou)=0.d0
      dpdnh201(nfou)=0.d0
      dpdnh021(nfou)=0.d0
      dpdnh111(nfou)=0.d0
      dpdnh102(nfou)=0.d0
      dpdnh012(nfou)=0.d0
      dpdn201(nfou)=0.d0
      dpdn021(nfou)=0.d0
      dpdn111(nfou)=0.d0
      dpdn102(nfou)=0.d0
      dpdn012(nfou)=0.d0
      enddo
      
      return
      end
      
c--------------------------------        
      subroutine zeroesdes      
c--------------------------------      
      include 'STGFM.cmn'      

      do nfou=-iord_max,iord_max
      eta101(nfou)=0.d0
      eta011(nfou)=0.d0
      eta111(nfou)=0.d0
      eta201(nfou)=0.d0
      eta021(nfou)=0.d0
      dreta101(nfou)=0.d0
      dreta011(nfou)=0.d0
      dphieta101(nfou)=0.d0
      dphieta011(nfou)=0.d0     
      drreta101(nfou)=0.d0
      drreta011(nfou)=0.d0     
      dpreta101(nfou)=0.d0
      dpreta011(nfou)=0.d0         
      dppeta101(nfou)=0.d0
      dppeta011(nfou)=0.d0
      drreta002(nfou)=0.d0
      dpreta002(nfou)=0.d0
      dppeta002(nfou)=0.d0     
      dreta111(nfou)=0.d0
      dreta201(nfou)=0.d0
      dreta021(nfou)=0.d0     
      dpeta111(nfou)=0.d0
      dpeta201(nfou)=0.d0
      dpeta021(nfou)=0.d0
      drreta102(nfou)=0.d0
      drreta012(nfou)=0.d0
      dpreta102(nfou)=0.d0
      dpreta012(nfou)=0.d0
      dppeta102(nfou)=0.d0
      dppeta012(nfou)=0.d0      
      eta002(nfou)=0.d0
      eta012(nfou)=0.d0
      eta102(nfou)=0.d0
      dreta002(nfou)=0.d0
      dphieta002(nfou)=0.d0   
      dreta102(nfou)=0.d0
      dpeta102(nfou)=0.d0
      dreta012(nfou)=0.d0
      dpeta012(nfou)=0.d0     
      dpppeta101(nfou)=0.d0
      dpppeta011(nfou)=0.d0
      dpppeta002(nfou)=0.d0
      dppreta101(nfou)=0.d0
      dppreta011(nfou)=0.d0
      dppreta002(nfou)=0.d0
      dprreta101(nfou)=0.d0
      dprreta011(nfou)=0.d0
      dprreta002(nfou)=0.d0
      drrreta101(nfou)=0.d0
      drrreta011(nfou)=0.d0
      drrreta002(nfou)=0.d0     
      eta202(nfou)=0.d0
      eta112(nfou)=0.d0
      eta022(nfou)=0.d0   
      eta003(nfou)=0.d0
      dreta202(nfou)=0.d0
      dreta112(nfou)=0.d0
      dreta022(nfou)=0.d0
      dreta003(nfou)=0.d0
      dpeta202(nfou)=0.d0
      dpeta112(nfou)=0.d0
      dpeta022(nfou)=0.d0
      dpeta003(nfou)=0.d0
      drreta201(nfou)=0.d0
      drreta021(nfou)=0.d0
      drreta111(nfou)=0.d0
      dpreta201(nfou)=0.d0
      dpreta021(nfou)=0.d0
      dpreta111(nfou)=0.d0
      dppeta201(nfou)=0.d0
      dppeta021(nfou)=0.d0
      dppeta111(nfou)=0.d0
      enddo  

      return
      end 

c--------------------------------        
      subroutine dndnh001(r,phi)      
c--------------------------------      
c
c     missing sin(m*phi), cos(m*phi) are implemented 
c
c---------------------------------   
      include 'STGFM.cmn' 

      dn001(0)=0.5d0*(s1+s2)*bn(0)
      dnh001(0)=-0.5d0*(s1+s2)*bnh(0)

      do nfou=-iord1,-1
      dn001(nfou)=(xi/xkxn(nfou))*dsin(xmo*phi)*
     &   (cdexp(-xi*xkxn(nfou)*s1)*bn(0)-bn(nfou))     
      dnh001(nfou)=-(xi/xkxn(nfou))*dcos(xmo*phi)*
     &   (cdexp(-xi*xkxn(nfou)*s1)*bnh(0)-bnh(nfou))
      enddo
 
      do nfou=1,iord1
      dn001(nfou)=(xi/xkxn(nfou))*dsin(xmo*phi)*
     &   (cdexp(-xi*xkxn(nfou)*s1)*bn(0)-bn(nfou))   
      dnh001(nfou)=-(xi/xkxn(nfou))*dcos(xmo*phi)*
     &   (cdexp(-xi*xkxn(nfou)*s1)*bnh(0)-bnh(nfou))  
      enddo            
  
      if(bz_zf_off==1)then
        do nfou=-iord1,-1
        dn001(0)=dn001(0)-dn001(nfou)*cdexp(xi*xkxn(nfou)*zf)
        dnh001(0)=dnh001(0)-dnh001(nfou)*cdexp(xi*xkxn(nfou)*zf)
        enddo
        do nfou=1,iord1
        dn001(0)=dn001(0)-dn001(nfou)*cdexp(xi*xkxn(nfou)*zf)
        dnh001(0)=dnh001(0)-dnh001(nfou)*cdexp(xi*xkxn(nfou)*zf)
        enddo
      endif
    
      return
      end    

c--------------------------------        
      subroutine ddndnh001      
c--------------------------------      
      include 'STGFM.cmn'
      drdn001(0)=0.5d0*(s1+s2)*drbn(0)
      drdnh001(0)=-0.5d0*(s1+s2)*drbnh(0)
         
      dphidn001(0)=0.5d0*(s1+s2)*dphibn(0)
      dphidnh001(0)=-0.5d0*(s1+s2)*dphibnh(0)

      do nfou=-iord1,-1
      drdn001(nfou)=(xi/xkxn(nfou))*
     &   (cdexp(-xi*xkxn(nfou)*s1)*drbn(0)-drbn(nfou))      
      dphidn001(nfou)=(xi/xkxn(nfou))*
     &   (cdexp(-xi*xkxn(nfou)*s1)*dphibn(0)-dphibn(nfou))
      drdnh001(nfou)=-(xi/xkxn(nfou))*
     &   (cdexp(-xi*xkxn(nfou)*s1)*drbnh(0)-drbnh(nfou))
      dphidnh001(nfou)=-(xi/xkxn(nfou))*
     &   (cdexp(-xi*xkxn(nfou)*s1)*dphibnh(0)-dphibnh(nfou))               
      enddo
      
      do nfou=1,iord1
      drdn001(nfou)=(xi/xkxn(nfou))*
     &   (cdexp(-xi*xkxn(nfou)*s1)*drbn(0)-drbn(nfou))       
      dphidn001(nfou)=(xi/xkxn(nfou))*
     &   (cdexp(-xi*xkxn(nfou)*s1)*dphibn(0)-dphibn(nfou)) 
      drdnh001(nfou)=-(xi/xkxn(nfou))*
     &   (cdexp(-xi*xkxn(nfou)*s1)*drbnh(0)-drbnh(nfou))
      dphidnh001(nfou)=-(xi/xkxn(nfou))*
     &   (cdexp(-xi*xkxn(nfou)*s1)*dphibnh(0)-dphibnh(nfou))                        
      enddo
     
      if(bz_zf_off==1)then
        do nfou=-iord1,-1
        drdn001(0)=drdn001(0)-drdn001(nfou)*cdexp(xi*xkxn(nfou)*zf)  
        dphidn001(0)=dphidn001(0)-dphidn001(nfou)*
     & cdexp(xi*xkxn(nfou)*zf)
 
        drdnh001(0)=drdnh001(0)-drdnh001(nfou)*cdexp(xi*xkxn(nfou)*zf)         
        dphidnh001(0)=dphidnh001(0)-dphidnh001(nfou)*
     & cdexp(xi*xkxn(nfou)*zf)      
        enddo
        do nfou=1,iord1
        drdn001(0)=drdn001(0)-drdn001(nfou)*cdexp(xi*xkxn(nfou)*zf)  
        dphidn001(0)=dphidn001(0)-dphidn001(nfou)*
     & cdexp(xi*xkxn(nfou)*zf)
 
        drdnh001(0)=drdnh001(0)-drdnh001(nfou)*cdexp(xi*xkxn(nfou)*zf)         
        dphidnh001(0)=dphidnh001(0)-dphidnh001(nfou)*
     & cdexp(xi*xkxn(nfou)*zf)      
        enddo
      endif
      
      return
      end
              
c--------------------------------        
      subroutine dddndnh001      
c--------------------------------      
      include 'STGFM.cmn' 

      dppdn001(0)=0.5d0*(s1+s2)*dppbn(0)
      dppdnh001(0)=-0.5d0*(s1+s2)*dppbnh(0)
      
      dprdn001(0)=0.5d0*(s1+s2)*dprbn(0)
      dprdnh001(0)=-0.5d0*(s1+s2)*dprbnh(0)
      
      drrdn001(0)=0.5d0*(s1+s2)*drrbn(0)
      drrdnh001(0)=-0.5d0*(s1+s2)*drrbnh(0)

      do nfou=-iord1,-1
      dppdn001(nfou)=(xi/xkxn(nfou))*
     &   (cdexp(-xi*xkxn(nfou)*s1)*dppbn(0)-dppbn(nfou))       
      dprdn001(nfou)=(xi/xkxn(nfou))*
     &   (cdexp(-xi*xkxn(nfou)*s1)*dprbn(0)-dprbn(nfou)) 
      drrdn001(nfou)=(xi/xkxn(nfou))*
     &   (cdexp(-xi*xkxn(nfou)*s1)*drrbn(0)-drrbn(nfou))       
      
      dppdnh001(nfou)=-(xi/xkxn(nfou))*
     &   (cdexp(-xi*xkxn(nfou)*s1)*dppbnh(0)-dppbnh(nfou))
      dprdnh001(nfou)=-(xi/xkxn(nfou))*
     &   (cdexp(-xi*xkxn(nfou)*s1)*dprbnh(0)-dprbnh(nfou))                
      drrdnh001(nfou)=-(xi/xkxn(nfou))*
     &   (cdexp(-xi*xkxn(nfou)*s1)*drrbnh(0)-drrbnh(nfou))
      enddo
      
      do nfou=1,iord1
      dppdn001(nfou)=(xi/xkxn(nfou))*
     &   (cdexp(-xi*xkxn(nfou)*s1)*dppbn(0)-dppbn(nfou))       
      dprdn001(nfou)=(xi/xkxn(nfou))*
     &   (cdexp(-xi*xkxn(nfou)*s1)*dprbn(0)-dprbn(nfou)) 
      drrdn001(nfou)=(xi/xkxn(nfou))*
     &   (cdexp(-xi*xkxn(nfou)*s1)*drrbn(0)-drrbn(nfou))       
      
      dppdnh001(nfou)=-(xi/xkxn(nfou))*
     &   (cdexp(-xi*xkxn(nfou)*s1)*dppbnh(0)-dppbnh(nfou)) 
      dprdnh001(nfou)=-(xi/xkxn(nfou))*
     &   (cdexp(-xi*xkxn(nfou)*s1)*dprbnh(0)-dprbnh(nfou))                 
      drrdnh001(nfou)=-(xi/xkxn(nfou))*
     &   (cdexp(-xi*xkxn(nfou)*s1)*drrbnh(0)-drrbnh(nfou))
      enddo
    
      if(bz_zf_off==1)then
        do nfou=-iord1,-1     
        dppdn001(0)=dppdn001(0)-dppdn001(nfou)*cdexp(xi*xkxn(nfou)*zf)  
        dprdn001(0)=dprdn001(0)-dprdn001(nfou)*cdexp(xi*xkxn(nfou)*zf)   
        drrdn001(0)=drrdn001(0)-drrdn001(nfou)*cdexp(xi*xkxn(nfou)*zf)
  
        dppdnh001(0)=dppdnh001(0)-dppdnh001(nfou)*
     & cdexp(xi*xkxn(nfou)*zf)           
        dprdnh001(0)=dprdnh001(0)-dprdnh001(nfou)*
     & cdexp(xi*xkxn(nfou)*zf)      
        drrdnh001(0)=drrdnh001(0)-drrdnh001(nfou)*
     & cdexp(xi*xkxn(nfou)*zf) 
        enddo
        do nfou=1,iord1     
        dppdn001(0)=dppdn001(0)-dppdn001(nfou)*
     & cdexp(xi*xkxn(nfou)*zf)  
        dprdn001(0)=dprdn001(0)-dprdn001(nfou)*
     & cdexp(xi*xkxn(nfou)*zf)   
        drrdn001(0)=drrdn001(0)-drrdn001(nfou)*
     & cdexp(xi*xkxn(nfou)*zf)
  
        dppdnh001(0)=dppdnh001(0)-dppdnh001(nfou)*
     & cdexp(xi*xkxn(nfou)*zf)           
        dprdnh001(0)=dprdnh001(0)-dprdnh001(nfou)*
     & cdexp(xi*xkxn(nfou)*zf)      
        drrdnh001(0)=drrdnh001(0)-drrdnh001(nfou)*
     & cdexp(xi*xkxn(nfou)*zf) 
        enddo
      endif
      
      return
      end
       
      
c--------------------------------        
      subroutine ddddndnh001      
c--------------------------------      
      include 'STGFM.cmn' 

      dpppdn001(0)=0.5d0*(s1+s2)*dpppbn(0)
      dpppdnh001(0)=-0.5d0*(s1+s2)*dpppbnh(0)
 
      dpprdn001(0)=0.5d0*(s1+s2)*dpprbn(0)
      dpprdnh001(0)=-0.5d0*(s1+s2)*dpprbnh(0)

      dprrdn001(0)=0.5d0*(s1+s2)*dprrbn(0)
      dprrdnh001(0)=-0.5d0*(s1+s2)*dprrbnh(0)   
   
      drrrdn001(0)=0.5d0*(s1+s2)*drrrbn(0)
      drrrdnh001(0)=-0.5d0*(s1+s2)*drrrbnh(0)
   
c------ dritte Ableitungen nach phi und r

      do nfou=-iord1,-1
      dpppdn001(nfou)=(xi/xkxn(nfou))*
     &   (cdexp(-xi*xkxn(nfou)*s1)*dpppbn(0)-dpppbn(nfou))     
      dpprdn001(nfou)=(xi/xkxn(nfou))*
     &   (cdexp(-xi*xkxn(nfou)*s1)*dpprbn(0)-dpprbn(nfou))     
      dprrdn001(nfou)=(xi/xkxn(nfou))*
     &   (cdexp(-xi*xkxn(nfou)*s1)*dprrbn(0)-dprrbn(nfou)) 
      drrrdn001(nfou)=(xi/xkxn(nfou))*
     &   (cdexp(-xi*xkxn(nfou)*s1)*drrrbn(0)-drrrbn(nfou))     
 
      dpppdnh001(nfou)=-(xi/xkxn(nfou))*
     &   (cdexp(-xi*xkxn(nfou)*s1)*dpppbnh(0)-dpppbnh(nfou))            
      dpprdnh001(nfou)=-(xi/xkxn(nfou))*
     &   (cdexp(-xi*xkxn(nfou)*s1)*dpprbnh(0)-dpprbnh(nfou)) 
      dprrdnh001(nfou)=-(xi/xkxn(nfou))*
     &   (cdexp(-xi*xkxn(nfou)*s1)*dprrbnh(0)-dprrbnh(nfou))                
      drrrdnh001(nfou)=-(xi/xkxn(nfou))*
     &   (cdexp(-xi*xkxn(nfou)*s1)*drrrbnh(0)-drrrbnh(nfou))
      enddo

      do nfou=1,iord1
      dpppdn001(nfou)=(xi/xkxn(nfou))*
     &   (cdexp(-xi*xkxn(nfou)*s1)*dpppbn(0)-dpppbn(nfou))     
      dpprdn001(nfou)=(xi/xkxn(nfou))*
     &   (cdexp(-xi*xkxn(nfou)*s1)*dpprbn(0)-dpprbn(nfou))     
      dprrdn001(nfou)=(xi/xkxn(nfou))*
     &   (cdexp(-xi*xkxn(nfou)*s1)*dprrbn(0)-dprrbn(nfou)) 
      drrrdn001(nfou)=(xi/xkxn(nfou))*
     &   (cdexp(-xi*xkxn(nfou)*s1)*drrrbn(0)-drrrbn(nfou))     
 
      dpppdnh001(nfou)=-(xi/xkxn(nfou))*
     &   (cdexp(-xi*xkxn(nfou)*s1)*dpppbnh(0)-dpppbnh(nfou))            
      dpprdnh001(nfou)=-(xi/xkxn(nfou))*
     &   (cdexp(-xi*xkxn(nfou)*s1)*dpprbnh(0)-dpprbnh(nfou)) 
      dprrdnh001(nfou)=-(xi/xkxn(nfou))*
     &   (cdexp(-xi*xkxn(nfou)*s1)*dprrbnh(0)-dprrbnh(nfou))                 
      drrrdnh001(nfou)=-(xi/xkxn(nfou))*
     &   (cdexp(-xi*xkxn(nfou)*s1)*drrrbnh(0)-drrrbnh(nfou)) 
      enddo     

      if(bz_zf_off==1)then
        do nfou=-iord1,-1    
        dpppdn001(0)=dpppdn001(0)-dpppdn001(nfou)*
     & cdexp(xi*xkxn(nfou)*zf)     
        dpprdn001(0)=dpprdn001(0)-dpprdn001(nfou)*
     & cdexp(xi*xkxn(nfou)*zf)
        dprrdn001(0)=dprrdn001(0)-dprrdn001(nfou)*
     & cdexp(xi*xkxn(nfou)*zf)      
        drrrdn001(0)=drrrdn001(0)-drrrdn001(nfou)*
     & cdexp(xi*xkxn(nfou)*zf)
 
        dpppdnh001(0)=dpppdnh001(0)-dpppdnh001(nfou)*
     & cdexp(xi*xkxn(nfou)*zf)            
        dpprdnh001(0)=dpprdnh001(0)-dpprdnh001(nfou)*
     & cdexp(xi*xkxn(nfou)*zf)           
        dprrdnh001(0)=dprrdnh001(0)-dprrdnh001(nfou)*
     & cdexp(xi*xkxn(nfou)*zf)      
        drrrdnh001(0)=drrrdnh001(0)-drrrdnh001(nfou)*
     & cdexp(xi*xkxn(nfou)*zf) 
        enddo
        do nfou=1,iord1    
        dpppdn001(0)=dpppdn001(0)-dpppdn001(nfou)*
     & cdexp(xi*xkxn(nfou)*zf)     
        dpprdn001(0)=dpprdn001(0)-dpprdn001(nfou)*
     & cdexp(xi*xkxn(nfou)*zf)
        dprrdn001(0)=dprrdn001(0)-dprrdn001(nfou)*
     & cdexp(xi*xkxn(nfou)*zf)      
        drrrdn001(0)=drrrdn001(0)-drrrdn001(nfou)*
     & cdexp(xi*xkxn(nfou)*zf)
 
        dpppdnh001(0)=dpppdnh001(0)-dpppdnh001(nfou)*
     & cdexp(xi*xkxn(nfou)*zf)            
        dpprdnh001(0)=dpprdnh001(0)-dpprdnh001(nfou)*
     & cdexp(xi*xkxn(nfou)*zf)           
        dprrdnh001(0)=dprrdnh001(0)-dprrdnh001(nfou)*
     & cdexp(xi*xkxn(nfou)*zf)      
        drrrdnh001(0)=drrrdnh001(0)-drrrdnh001(nfou)*
     & cdexp(xi*xkxn(nfou)*zf) 
        enddo
      endif
            
      return
      end

c--------------------------------        
      subroutine eta2(r,phi)      
c--------------------------------      
      include 'STGFM.cmn'           

      ds=dsin(phi)
      dc=dcos(phi)
  
      do nfou=-iord1,iord1
      eta101(nfou)=-dc*dnh001(nfou)+ds*dn001(nfou)
      eta011(nfou)=-ds*dnh001(nfou)-dc*dn001(nfou)
      enddo
      
      do nfou=-iord2,iord2      
      do m=-iord1,iord1
        if(iabs(nfou-m).le.iord1)then
       eta002(nfou)=eta002(nfou)-0.5d0*(
     & dnh001(m)*dnh001(nfou-m)+dn001(m)*dn001(nfou-m))
       endif
      enddo
      enddo
      
      if(iverbose.eq.1)then
      write(6,*)' eta101(0) ',eta101(0)
      write(6,*)' dnh001(0) ',dnh001(0)
      write(6,*)' dn001(0) ',dn001(0)
      endif
      
      return
      end      
      
c--------------------------------        
      subroutine deta2(r,phi)      
c--------------------------------      
      include 'STGFM.cmn'   

      do nfou=-iord1,iord1
      dreta101(nfou)=-dcos(phi)*drdnh001(nfou)+dsin(phi)*drdn001(nfou)
c      if(iverbose.eq.1)write(6,*)' dreta101 ',nfou,dreta101(nfou)       
      dreta011(nfou)=-dsin(phi)*drdnh001(nfou)-dcos(phi)*drdn001(nfou)
c      if(iverbose.eq.1)write(6,*)' dreta011 ',nfou,dreta011(nfou)
      enddo 

c     write(6,*)'ini dreta002(0), dreta002(1) ',dreta002(0), dreta002(1)
      do nfou=-iord2,iord2      
      do m=-iord1,iord1
        if(iabs(nfou-m).le.iord1)then
        dreta002(nfou)=dreta002(nfou)-0.5d0*(
     &  drdnh001(m)*dnh001(nfou-m)+dnh001(m)*drdnh001(nfou-m)+
     &  drdn001(m)*dn001(nfou-m)+dn001(m)*drdn001(nfou-m))     
        endif
      enddo    
      enddo

      do nfou=-iord1,iord1
      dphieta101(nfou)=dsin(phi)*dnh001(nfou)+dcos(phi)*dn001(nfou)-
     &   dcos(phi)*dphidnh001(nfou)+dsin(phi)*dphidn001(nfou)     
c      if(iverbose.eq.1)write(6,*)' dphieta101 ',nfou,dphieta101(nfou)      
      dphieta011(nfou)=-dcos(phi)*dnh001(nfou)+dsin(phi)*dn001(nfou)-
     &   dsin(phi)*dphidnh001(nfou)-dcos(phi)*dphidn001(nfou)
      enddo
      
      do nfou=-iord2,iord2   
      do m=-iord1,iord1
        if(iabs(nfou-m).le.iord1)then
        dphieta002(nfou)=dphieta002(nfou)-0.5d0*(
     &  dphidnh001(m)*dnh001(nfou-m)+dnh001(m)*dphidnh001(nfou-m)+
     &  dphidn001(m)*dn001(nfou-m)+dn001(m)*dphidn001(nfou-m)) 
        endif    
      enddo
      enddo
      
      return
      end

c--------------------------------        
      subroutine ddeta2(r,phi)      
c--------------------------------      
      include 'STGFM.cmn'           

      do nfou=-iord1,iord1
      drreta101(nfou)=-dcos(phi)*drrdnh001(nfou)+
     &     dsin(phi)*drrdn001(nfou)    
      dpreta101(nfou)=dsin(phi)*drdnh001(nfou)+dcos(phi)*drdn001(nfou)-
     &     dcos(phi)*dprdnh001(nfou)+dsin(phi)*dprdn001(nfou)     
      dppeta101(nfou)=dcos(phi)*dnh001(nfou)-dsin(phi)*dn001(nfou)+
     &     dsin(phi)*dphidnh001(nfou)+dcos(phi)*dphidn001(nfou)+
     &     dsin(phi)*dphidnh001(nfou)+dcos(phi)*dphidn001(nfou)-
     &     dcos(phi)*dppdnh001(nfou)+dsin(phi)*dppdn001(nfou)
 
      drreta011(nfou)=-dsin(phi)*drrdnh001(nfou)-
     &     dcos(phi)*drrdn001(nfou)    
      dpreta011(nfou)=-dcos(phi)*drdnh001(nfou)+dsin(phi)*drdn001(nfou)-
     &     dsin(phi)*dprdnh001(nfou)-dcos(phi)*dprdn001(nfou)
      dppeta011(nfou)=dsin(phi)*dnh001(nfou)+dcos(phi)*dn001(nfou)-
     &     dcos(phi)*dphidnh001(nfou)+dsin(phi)*dphidn001(nfou)-
     &     dcos(phi)*dphidnh001(nfou)+dsin(phi)*dphidn001(nfou)-
     &     dsin(phi)*dppdnh001(nfou)-dcos(phi)*dppdn001(nfou)      
      enddo
         
      do nfou=-iord2,iord2      
      do m=-iord1,iord1
        if(iabs(nfou-m).le.iord1)then
        drreta002(nfou)=drreta002(nfou)-0.5d0*(
     &  drrdnh001(m)*dnh001(nfou-m)+drdnh001(m)*drdnh001(nfou-m)+
     &  drdnh001(m)*drdnh001(nfou-m)+dnh001(m)*drrdnh001(nfou-m)+
     &  drrdn001(m)*dn001(nfou-m)+drdn001(m)*drdn001(nfou-m)+   
     &  drdn001(m)*drdn001(nfou-m)+dn001(m)*drrdn001(nfou-m))
        
        dpreta002(nfou)=dpreta002(nfou)-0.5d0*(
     &  dprdnh001(m)*dnh001(nfou-m)+drdnh001(m)*dphidnh001(nfou-m)+
     &  dphidnh001(m)*drdnh001(nfou-m)+dnh001(m)*dprdnh001(nfou-m)+
     &  dprdn001(m)*dn001(nfou-m)+drdn001(m)*dphidn001(nfou-m)+   
     &  dphidn001(m)*drdn001(nfou-m)+dn001(m)*dprdn001(nfou-m))       
        
        dppeta002(nfou)=dppeta002(nfou)-0.5d0*(
     &  dppdnh001(m)*dnh001(nfou-m)+dphidnh001(m)*dphidnh001(nfou-m)+
     &  dphidnh001(m)*dphidnh001(nfou-m)+dnh001(m)*dppdnh001(nfou-m)+
     &  dppdn001(m)*dn001(nfou-m)+dphidn001(m)*dphidn001(nfou-m)+   
     &  dphidn001(m)*dphidn001(nfou-m)+dn001(m)*dppdn001(nfou-m)) 
        endif
      enddo    
      enddo
       
ctesttest
      if(iverbose.eq.1)then
      write(6,*)'-------------------------------'
      write(6,*)' drreta101(1) ',drreta101(1)
      write(6,*)' drreta011(1) ',drreta011(1)
      write(6,*)'-------------------------------'
      
      write(6,*)' drrdnh101(1) ',drrdnh101(1)
      write(6,*)' drrdnh011(1) ',drrdnh011(1)
      write(6,*)'-------------------------------'    
      write(6,*)' drrdn101(1) ',drrdn101(1)
      write(6,*)' drrdn011(1) ',drrdn011(1)
      write(6,*)'-------------------------------'    
           
      write(6,*)' drrdnh001(1) ',drrdnh001(1)
      write(6,*)' drrdn001(1) ',drrdn001(1)
      write(6,*)'-------------------------------'     
      endif
      
      return
      end

c--------------------------------        
      subroutine dddeta2(r,phi)      
c--------------------------------      
      include 'STGFM.cmn'           

c--------- Ableitungen nach r      
      
      do nfou=-iord1,iord1
      drrreta101(nfou)=-dcos(phi)*drrrdnh001(nfou)+
     &     dsin(phi)*drrrdn001(nfou)    
      dprreta101(nfou)=dsin(phi)*drrdnh001(nfou)+
     &     dcos(phi)*drrdn001(nfou)-
     &     dcos(phi)*dprrdnh001(nfou)+dsin(phi)*dprrdn001(nfou)     
      dppreta101(nfou)=dcos(phi)*drdnh001(nfou)-dsin(phi)*drdn001(nfou)+
     &     dsin(phi)*dprdnh001(nfou)+dcos(phi)*dprdn001(nfou)+
     &     dsin(phi)*dprdnh001(nfou)+dcos(phi)*dprdn001(nfou)-
     &     dcos(phi)*dpprdnh001(nfou)+dsin(phi)*dpprdn001(nfou) 

c      if(nfou.eq.0)then
c       write(6,*)' dppreta101(0) kontrolle ',dppreta101(0)
c       write(6,*)' '
c       write(6,*)' non-trivial summands in dppreta101(0) at phi=0: '
c       write(6,*)' drdnh001(nfou) ',drdnh001(nfou)
c       write(6,*)' dprdn001(nfou) ',dprdn001(nfou)
c       write(6,*)' dpprdnh001(nfou) ',dpprdnh001(nfou)
c       write(6,*)' '
c      endif
      drrreta011(nfou)=-dsin(phi)*drrrdnh001(nfou)-
     &     dcos(phi)*drrrdn001(nfou)    
      dprreta011(nfou)=-dcos(phi)*drrdnh001(nfou)+
     &     dsin(phi)*drrdn001(nfou)-
     &     dsin(phi)*dprrdnh001(nfou)-dcos(phi)*dprrdn001(nfou)
      dppreta011(nfou)=dsin(phi)*drdnh001(nfou)+dcos(phi)*drdn001(nfou)-
     &     dcos(phi)*dprdnh001(nfou)+dsin(phi)*dprdn001(nfou)-
     &     dcos(phi)*dprdnh001(nfou)+dsin(phi)*dprdn001(nfou)-
     &     dsin(phi)*dpprdnh001(nfou)-dcos(phi)*dpprdn001(nfou)      
      enddo
         
      do nfou=-iord2,iord2      
      do m=-iord1,iord1
        if(iabs(nfou-m).le.iord1)then
        drrreta002(nfou)=drrreta002(nfou)-0.5d0*(
     &  drrrdnh001(m)*dnh001(nfou-m)+drrdnh001(m)*drdnh001(nfou-m)+
     &  drrdnh001(m)*drdnh001(nfou-m)+drdnh001(m)*drrdnh001(nfou-m)+
     &  drrrdn001(m)*dn001(nfou-m)+drrdn001(m)*drdn001(nfou-m)+   
     &  drrdn001(m)*drdn001(nfou-m)+drdn001(m)*drrdn001(nfou-m)+
     &  drrdnh001(m)*drdnh001(nfou-m)+drdnh001(m)*drrdnh001(nfou-m)+
     &  drdnh001(m)*drrdnh001(nfou-m)+dnh001(m)*drrrdnh001(nfou-m)+
     &  drrdn001(m)*drdn001(nfou-m)+drdn001(m)*drrdn001(nfou-m)+   
     &  drdn001(m)*drrdn001(nfou-m)+dn001(m)*drrrdn001(nfou-m))   
 
        dprreta002(nfou)=dprreta002(nfou)-0.5d0*(
     &  dprrdnh001(m)*dnh001(nfou-m)+drrdnh001(m)*dphidnh001(nfou-m)+
     &  dprdnh001(m)*drdnh001(nfou-m)+drdnh001(m)*dprdnh001(nfou-m)+
     &  dprrdn001(m)*dn001(nfou-m)+drrdn001(m)*dphidn001(nfou-m)+   
     &  dprdn001(m)*drdn001(nfou-m)+drdn001(m)*dprdn001(nfou-m)+
     &  dprdnh001(m)*drdnh001(nfou-m)+drdnh001(m)*dprdnh001(nfou-m)+
     &  dphidnh001(m)*drrdnh001(nfou-m)+dnh001(m)*dprrdnh001(nfou-m)+
     &  dprdn001(m)*drdn001(nfou-m)+drdn001(m)*dprdn001(nfou-m)+   
     &  dphidn001(m)*drrdn001(nfou-m)+dn001(m)*dprrdn001(nfou-m))      
        
        dppreta002(nfou)=dppreta002(nfou)-0.5d0*(
     &  dpprdnh001(m)*dnh001(nfou-m)+dprdnh001(m)*dphidnh001(nfou-m)+
     &  dprdnh001(m)*dphidnh001(nfou-m)+drdnh001(m)*dppdnh001(nfou-m)+
     &  dpprdn001(m)*dn001(nfou-m)+dprdn001(m)*dphidn001(nfou-m)+   
     &  dprdn001(m)*dphidn001(nfou-m)+drdn001(m)*dppdn001(nfou-m)+
     &  dppdnh001(m)*drdnh001(nfou-m)+dphidnh001(m)*dprdnh001(nfou-m)+
     &  dphidnh001(m)*dprdnh001(nfou-m)+dnh001(m)*dpprdnh001(nfou-m)+
     &  dppdn001(m)*drdn001(nfou-m)+dphidn001(m)*dprdn001(nfou-m)+   
     &  dphidn001(m)*dprdn001(nfou-m)+dn001(m)*dpprdn001(nfou-m)) 
        endif
      enddo    
      enddo

c---------- Ableitungen nach phi

      do nfou=-iord1,iord1
      dpppeta101(nfou)=-dsin(phi)*dnh001(nfou)-dcos(phi)*dn001(nfou)+
     &     dcos(phi)*dphidnh001(nfou)-dsin(phi)*dphidn001(nfou)+
     &     dcos(phi)*dphidnh001(nfou)-dsin(phi)*dphidn001(nfou)+
     &     dsin(phi)*dppdnh001(nfou)+dcos(phi)*dppdn001(nfou)+      
     &     dcos(phi)*dphidnh001(nfou)-dsin(phi)*dphidn001(nfou)+
     &     dsin(phi)*dppdnh001(nfou)+dcos(phi)*dppdn001(nfou)+
     &     dsin(phi)*dppdnh001(nfou)+dcos(phi)*dppdn001(nfou)-
     &     dcos(phi)*dpppdnh001(nfou)+dsin(phi)*dpppdn001(nfou)      
      
      dpppeta011(nfou)=dcos(phi)*dnh001(nfou)-dsin(phi)*dn001(nfou)+
     &     dsin(phi)*dphidnh001(nfou)+dcos(phi)*dphidn001(nfou)+
     &     dsin(phi)*dphidnh001(nfou)+dcos(phi)*dphidn001(nfou)-
     &     dcos(phi)*dppdnh001(nfou)+dsin(phi)*dppdn001(nfou)+      
     &     dsin(phi)*dphidnh001(nfou)+dcos(phi)*dphidn001(nfou)-
     &     dcos(phi)*dppdnh001(nfou)+dsin(phi)*dppdn001(nfou)-
     &     dcos(phi)*dppdnh001(nfou)+dsin(phi)*dppdn001(nfou)-
     &     dsin(phi)*dpppdnh001(nfou)-dcos(phi)*dpppdn001(nfou)    
      enddo
         
      do nfou=-iord2,iord2      
      do m=-iord1,iord1
        if(iabs(nfou-m).le.iord1)then
        dpppeta002(nfou)=dpppeta002(nfou)-0.5d0*(
     &  dpppdnh001(m)*dnh001(nfou-m)+dppdnh001(m)*dphidnh001(nfou-m)+
     &  dppdnh001(m)*dphidnh001(nfou-m)+dphidnh001(m)*dppdnh001(nfou-m)+
     &  dpppdn001(m)*dn001(nfou-m)+dppdn001(m)*dphidn001(nfou-m)+   
     &  dppdn001(m)*dphidn001(nfou-m)+dphidn001(m)*dppdn001(nfou-m)+
     &  dppdnh001(m)*dphidnh001(nfou-m)+dphidnh001(m)*dppdnh001(nfou-m)+
     &  dphidnh001(m)*dppdnh001(nfou-m)+dnh001(m)*dpppdnh001(nfou-m)+
     &  dppdn001(m)*dphidn001(nfou-m)+dphidn001(m)*dppdn001(nfou-m)+   
     &  dphidn001(m)*dppdn001(nfou-m)+dn001(m)*dpppdn001(nfou-m)) 
        endif
      enddo    
      enddo
             
      
      return
      end
      
c-------------------------------------------
      subroutine dndnh2(r,phi)
c-------------------------------------------     
      include 'STGFM.cmn'           

      do nfou=-iord1,-1
      dnh011(nfou)=(xi/xkxn(nfou))*
     &  (cdexp(-xi*xkxn(nfou)*s1)*dreta011(0)-dreta011(nfou))
      dnh011(0)=dnh011(0)-dnh011(nfou)*cdexp(xi*xkxn(nfou)*zf)
      enddo
      
      do nfou=1,iord1
      dnh011(nfou)=(xi/xkxn(nfou))*
     &  (cdexp(-xi*xkxn(nfou)*s1)*dreta011(0)-dreta011(nfou))
      dnh011(0)=dnh011(0)-dnh011(nfou)*cdexp(xi*xkxn(nfou)*zf)
      enddo
  
      do nfou=-iord1,-1
      dnh101(nfou)=(xi/xkxn(nfou))*
     &  (cdexp(-xi*xkxn(nfou)*s1)*dreta101(0)-dreta101(nfou))
      dnh101(0)=dnh101(0)-dnh101(nfou)*cdexp(xi*xkxn(nfou)*zf)
      enddo
      
      do nfou=1,iord1
      dnh101(nfou)=(xi/xkxn(nfou))*
     &  (cdexp(-xi*xkxn(nfou)*s1)*dreta101(0)-dreta101(nfou))
      dnh101(0)=dnh101(0)-dnh101(nfou)*cdexp(xi*xkxn(nfou)*zf)
      enddo

c---------------------------------

      do nfou=-iord1,-1
      dn011(nfou)=(1.d0/r)*(xi/xkxn(nfou))*
     &  (cdexp(-xi*xkxn(nfou)*s1)*dphieta011(0)-dphieta011(nfou))
      dn011(0)=dn011(0)-dn011(nfou)*cdexp(xi*xkxn(nfou)*zf)
      enddo
      
      do nfou=1,iord1
      dn011(nfou)=(1.d0/r)*(xi/xkxn(nfou))*
     &  (cdexp(-xi*xkxn(nfou)*s1)*dphieta011(0)-dphieta011(nfou))
      dn011(0)=dn011(0)-dn011(nfou)*cdexp(xi*xkxn(nfou)*zf)
      enddo
  
      do nfou=-iord1,-1
      dn101(nfou)=(1.d0/r)*(xi/xkxn(nfou))*
     &  (cdexp(-xi*xkxn(nfou)*s1)*dphieta101(0)-dphieta101(nfou))
      dn101(0)=dn101(0)-dn101(nfou)*cdexp(xi*xkxn(nfou)*zf)
      enddo
      
      do nfou=1,iord1
      dn101(nfou)=(1.d0/r)*(xi/xkxn(nfou))*
     &  (cdexp(-xi*xkxn(nfou)*s1)*dphieta101(0)-dphieta101(nfou))
      dn101(0)=dn101(0)-dn101(nfou)*cdexp(xi*xkxn(nfou)*zf)
      enddo
      
      do nfou=-iord2,-1
      dnh002(nfou)=(xi/xkxn(nfou))*
     &  (cdexp(-xi*xkxn(nfou)*s1)*dreta002(0)-dreta002(nfou))
      dnh002(0)=dnh002(0)-dnh002(nfou)*cdexp(xi*xkxn(nfou)*zf)
      enddo

      do nfou=1,iord2
      dnh002(nfou)=(xi/xkxn(nfou))*
     &  (cdexp(-xi*xkxn(nfou)*s1)*dreta002(0)-dreta002(nfou))
      dnh002(0)=dnh002(0)-dnh002(nfou)*cdexp(xi*xkxn(nfou)*zf)
      enddo
      
c--------------------------------------------

      do nfou=-iord2,-1
      dn002(nfou)=(1.d0/r)*(xi/xkxn(nfou))*
     &  (cdexp(-xi*xkxn(nfou)*s1)*dphieta002(0)-dphieta002(nfou))
      dn002(0)=dn002(0)-dn002(nfou)*cdexp(xi*xkxn(nfou)*zf)
      enddo

      do nfou=1,iord2
      dn002(nfou)=(1.d0/r)*(xi/xkxn(nfou))*
     &  (cdexp(-xi*xkxn(nfou)*s1)*dphieta002(0)-dphieta002(nfou))
      dn002(0)=dn002(0)-dn002(nfou)*cdexp(xi*xkxn(nfou)*zf)
      enddo
 
      return
      end

c--------------------------------        
      subroutine eta3(r,phi)      
c--------------------------------      
      include 'STGFM.cmn'           

      ds=dsin(phi)
      dc=dcos(phi)

      do nfou=-iord1,iord1
      eta111(nfou)=-dc*dnh011(nfou)+ds*dn011(nfou)-
     &              ds*dnh101(nfou)-dc*dn101(nfou)
      eta201(nfou)=-dc*dnh101(nfou)+ds*dn101(nfou)
      eta021(nfou)=-ds*dnh011(nfou)-dc*dn011(nfou)
      enddo    
      
      do nfou=-iord2,iord2
      eta102(nfou)=-dcos(phi)*dnh002(nfou)+dsin(phi)*dn002(nfou)
      do m=-iord1,iord1
       if(iabs(nfou-m).le.iord1)then
         eta102(nfou)=eta102(nfou)-0.5d0*(
     &   dnh001(m)*dnh101(nfou-m)+dn001(m)*dn101(nfou-m)+
     &   dnh001(nfou-m)*dnh101(m)+dn001(nfou-m)*dn101(m))
       endif
      enddo
      enddo    
      
c      write(6,*)' dn002(1) ',dn002(1)
c      write(6,*)' dnh002(1) ',dnh002(1)  
c      write(6,*)'-------------------------'
      do nfou=-iord2,iord2
      eta012(nfou)=-dsin(phi)*dnh002(nfou)-dcos(phi)*dn002(nfou)
      do m=-iord1,iord1
       if(iabs(nfou-m).le.iord1)then
         eta012(nfou)=eta012(nfou)-0.5d0*(
     &   dnh001(m)*dnh011(nfou-m)+dn001(m)*dn011(nfou-m)+
     &   dnh001(nfou-m)*dnh011(m)+dn001(nfou-m)*dn011(m))
       endif
      enddo
c      if(abs(nfou).lt.3.d0)then
c       write(6,*)' nfou, eta012'
c       write(6,*)'   ',nfou,dreal(eta012(nfou)),dimag(eta012(nfou))
c      endif
      enddo 
 
      do nfou=-iord3,iord3
      do m=-iord2,iord2
       if(iabs(nfou-m).le.iord2)then
         eta003(nfou)=eta003(nfou)-0.5d0*(
     &   dnh001(m)*dnh002(nfou-m)+dnh002(m)*dnh001(nfou-m)+
     &   dn001(m)*dn002(nfou-m)+dn002(m)*dn001(nfou-m))
       endif
      enddo
      enddo 
      
      return
      end      
  
c--------------------------------        
      subroutine deta3(r,phi)      
c--------------------------------      
      include 'STGFM.cmn'           
      
      ds=dsin(phi)
      dc=dcos(phi)

      do nfou=-iord1,iord1
      dreta111(nfou)=-dc*drdnh011(nfou)+ds*drdn011(nfou)-
     &                ds*drdnh101(nfou)-dc*drdn101(nfou)
      dreta201(nfou)=-dc*drdnh101(nfou)+ds*drdn101(nfou)
      dreta021(nfou)=-ds*drdnh011(nfou)-dc*drdn011(nfou)
      
      dpeta111(nfou)=-dc*dpdnh011(nfou)+ds*dpdn011(nfou)-
     &                ds*dpdnh101(nfou)-dc*dpdn101(nfou)+
     &                ds*dnh011(nfou)+dc*dn011(nfou)-
     &                dc*dnh101(nfou)+ds*dn101(nfou)      
      dpeta201(nfou)=-dc*dpdnh101(nfou)+ds*dpdn101(nfou)+     
     &                ds*dnh101(nfou)+dc*dn101(nfou)
      dpeta021(nfou)=-ds*dpdnh011(nfou)-dc*dpdn011(nfou)-     
     &                dc*dnh011(nfou)+ds*dn011(nfou)         
      enddo
      
      do nfou=-iord2,iord2      
      dreta102(nfou)=-dc*drdnh002(nfou)+ds*drdn002(nfou)
      dreta012(nfou)=-ds*drdnh002(nfou)-dc*drdn002(nfou)
      dpeta102(nfou)=-dc*dpdnh002(nfou)+ds*dpdn002(nfou)+
     &   ds*dnh002(nfou)+dc*dn002(nfou) 
      dpeta012(nfou)=-ds*dpdnh002(nfou)-dc*dpdn002(nfou)- 
     &   dc*dnh002(nfou)+ds*dn002(nfou)   
      
      drreta102(nfou)=-dc*drrdnh002(nfou)+ds*drrdn002(nfou)
      drreta012(nfou)=-ds*drrdnh002(nfou)-dc*drrdn002(nfou)
      dpreta102(nfou)=-dc*dprdnh002(nfou)+ds*dprdn002(nfou)+
     &   ds*drdnh002(nfou)+dc*drdn002(nfou) 
      dpreta012(nfou)=-ds*dprdnh002(nfou)-dc*dprdn002(nfou)- 
     &   dc*drdnh002(nfou)+ds*drdn002(nfou)     
      dppeta102(nfou)=-dc*dppdnh002(nfou)+ds*dppdn002(nfou)+
     &   ds*dpdnh002(nfou)+dc*dpdn002(nfou) 
     &               +ds*dpdnh002(nfou)+dc*dpdn002(nfou)+
     &   dc*dnh002(nfou)-ds*dn002(nfou)     
      dppeta012(nfou)=-ds*dppdnh002(nfou)-dc*dppdn002(nfou)- 
     &   dc*dpdnh002(nfou)+ds*dpdn002(nfou)        
     &                -dc*dpdnh002(nfou)+ds*dpdn002(nfou)+ 
     &   ds*dnh002(nfou)+dc*dn002(nfou)
      enddo        
c----- 2. Ableitungen von eta201, eta111, eta021
      do nfou=-iord1,iord1             
c----- Ableitungen nach r
      drreta111(nfou)=-dc*drrdnh011(nfou)+ds*drrdn011(nfou)-
     &                ds*drrdnh101(nfou)-dc*drrdn101(nfou)
      drreta201(nfou)=-dc*drrdnh101(nfou)+ds*drrdn101(nfou)
      drreta021(nfou)=-ds*drrdnh011(nfou)-dc*drrdn011(nfou)
      
      dpreta111(nfou)=-dc*dprdnh011(nfou)+ds*dprdn011(nfou)-
     &                ds*dprdnh101(nfou)-dc*dprdn101(nfou)+
     &                ds*drdnh011(nfou)+dc*drdn011(nfou)-
     &                dc*drdnh101(nfou)+ds*drdn101(nfou)      
      dpreta201(nfou)=-dc*dprdnh101(nfou)+ds*dprdn101(nfou)+     
     &                ds*drdnh101(nfou)+dc*drdn101(nfou)
      dpreta021(nfou)=-ds*dprdnh011(nfou)-dc*dprdn011(nfou)-     
     &                dc*drdnh011(nfou)+ds*drdn011(nfou)              
      
c------ Ableitungen nach phi           
      dppeta111(nfou)=ds*dpdnh011(nfou)+dc*dpdn011(nfou)-
     &                dc*dpdnh101(nfou)+ds*dpdn101(nfou)+
     &                dc*dnh011(nfou)-ds*dn011(nfou)+
     &                ds*dnh101(nfou)+dc*dn101(nfou)-     
     &                dc*dppdnh011(nfou)+ds*dppdn011(nfou)-
     &                ds*dppdnh101(nfou)-dc*dppdn101(nfou)+
     &                ds*dpdnh011(nfou)+dc*dpdn011(nfou)-
     &                dc*dpdnh101(nfou)+ds*dpdn101(nfou)        
      dppeta201(nfou)=ds*dpdnh101(nfou)+dc*dpdn101(nfou)+     
     &                dc*dnh101(nfou)-ds*dn101(nfou)-
     &                dc*dppdnh101(nfou)+ds*dppdn101(nfou)+     
     &                ds*dpdnh101(nfou)+dc*dpdn101(nfou)
      dppeta021(nfou)=-dc*dpdnh011(nfou)+ds*dpdn011(nfou)+     
     &                ds*dnh011(nfou)+dc*dn011(nfou)-               
     &                ds*dppdnh011(nfou)-dc*dppdn011(nfou)-     
     &                dc*dpdnh011(nfou)+ds*dpdn011(nfou)   
      enddo

      do nfou=-iord2,iord2
      do m=-iord1,iord1
       if(iabs(nfou-m).le.iord1)then
        dreta102(nfou)=dreta102(nfou)-0.5d0*(
     &   drdnh001(m)*dnh101(nfou-m)+dnh001(m)*drdnh101(nfou-m)+
     &   drdnh101(m)*dnh001(nfou-m)+dnh101(m)*drdnh001(nfou-m)+
     &   drdn001(m)*dn101(nfou-m)+dn001(m)*drdn101(nfou-m)+
     &   drdn101(m)*dn001(nfou-m)+dn101(m)*drdn001(nfou-m))
       
        dreta012(nfou)=dreta012(nfou)-0.5d0*(
     &   drdnh001(m)*dnh011(nfou-m)+dnh001(m)*drdnh011(nfou-m)+
     &   drdnh011(m)*dnh001(nfou-m)+dnh011(m)*drdnh001(nfou-m)+
     &   drdn001(m)*dn011(nfou-m)+dn001(m)*drdn011(nfou-m)+
     &   drdn011(m)*dn001(nfou-m)+dn011(m)*drdn001(nfou-m))
 
        dpeta102(nfou)=dpeta102(nfou)-0.5d0*(
     &   dphidnh001(m)*dnh101(nfou-m)+dnh001(m)*dpdnh101(nfou-m)+
     &   dpdnh101(m)*dnh001(nfou-m)+dnh101(m)*dphidnh001(nfou-m)+
     &   dphidn001(m)*dn101(nfou-m)+dn001(m)*dpdn101(nfou-m)+
     &   dpdn101(m)*dn001(nfou-m)+dn101(m)*dphidn001(nfou-m))
       
        dpeta012(nfou)=dpeta012(nfou)-0.5d0*(
     &   dphidnh001(m)*dnh011(nfou-m)+dnh001(m)*dpdnh011(nfou-m)+
     &   dpdnh011(m)*dnh001(nfou-m)+dnh011(m)*dphidnh001(nfou-m)+
     &   dphidn001(m)*dn011(nfou-m)+dn001(m)*dpdn011(nfou-m)+
     &   dpdn011(m)*dn001(nfou-m)+dn011(m)*dphidn001(nfou-m))

c  doppelte Ableitungen nach r und phi      

        drreta102(nfou)=drreta102(nfou)-0.5d0*(
     &   drrdnh001(m)*dnh101(nfou-m)+dnh001(m)*drrdnh101(nfou-m)+
     &   drrdnh101(m)*dnh001(nfou-m)+dnh101(m)*drrdnh001(nfou-m)+
     &   drrdn001(m)*dn101(nfou-m)+dn001(m)*drrdn101(nfou-m)+
     &   drrdn101(m)*dn001(nfou-m)+dn101(m)*drrdn001(nfou-m)+        
     &   drdnh001(m)*drdnh101(nfou-m)+drdnh001(m)*drdnh101(nfou-m)+
     &   drdnh101(m)*drdnh001(nfou-m)+drdnh101(m)*drdnh001(nfou-m)+
     &   drdn001(m)*drdn101(nfou-m)+drdn001(m)*drdn101(nfou-m)+
     &   drdn101(m)*drdn001(nfou-m)+drdn101(m)*drdn001(nfou-m))
       
        drreta012(nfou)=drreta012(nfou)-0.5d0*(
     &   drrdnh001(m)*dnh011(nfou-m)+dnh001(m)*drrdnh011(nfou-m)+
     &   drrdnh011(m)*dnh001(nfou-m)+dnh011(m)*drrdnh001(nfou-m)+
     &   drrdn001(m)*dn011(nfou-m)+dn001(m)*drrdn011(nfou-m)+
     &   drrdn011(m)*dn001(nfou-m)+dn011(m)*drrdn001(nfou-m)+
     &   drdnh001(m)*drdnh011(nfou-m)+drdnh001(m)*drdnh011(nfou-m)+
     &   drdnh011(m)*drdnh001(nfou-m)+drdnh011(m)*drdnh001(nfou-m)+
     &   drdn001(m)*drdn011(nfou-m)+drdn001(m)*drdn011(nfou-m)+
     &   drdn011(m)*drdn001(nfou-m)+drdn011(m)*drdn001(nfou-m))
 
        dpreta102(nfou)=dpreta102(nfou)-0.5d0*(
     &   dprdnh001(m)*dnh101(nfou-m)+dnh001(m)*dprdnh101(nfou-m)+
     &   dprdnh101(m)*dnh001(nfou-m)+dnh101(m)*dprdnh001(nfou-m)+
     &   dprdn001(m)*dn101(nfou-m)+dn001(m)*dprdn101(nfou-m)+
     &   dprdn101(m)*dn001(nfou-m)+dn101(m)*dprdn001(nfou-m)+
     &   dphidnh001(m)*drdnh101(nfou-m)+drdnh001(m)*dpdnh101(nfou-m)+
     &   dpdnh101(m)*drdnh001(nfou-m)+drdnh101(m)*dphidnh001(nfou-m)+
     &   dphidn001(m)*drdn101(nfou-m)+drdn001(m)*dpdn101(nfou-m)+
     &   dpdn101(m)*drdn001(nfou-m)+drdn101(m)*dphidn001(nfou-m))
       
        dpreta012(nfou)=dpreta012(nfou)-0.5d0*(
     &   dprdnh001(m)*dnh011(nfou-m)+dnh001(m)*dprdnh011(nfou-m)+
     &   dprdnh011(m)*dnh001(nfou-m)+dnh011(m)*dprdnh001(nfou-m)+
     &   dprdn001(m)*dn011(nfou-m)+dn001(m)*dprdn011(nfou-m)+
     &   dprdn011(m)*dn001(nfou-m)+dn011(m)*dprdn001(nfou-m)+
     &   dphidnh001(m)*drdnh011(nfou-m)+drdnh001(m)*dpdnh011(nfou-m)+
     &   dpdnh011(m)*drdnh001(nfou-m)+drdnh011(m)*dphidnh001(nfou-m)+
     &   dphidn001(m)*drdn011(nfou-m)+drdn001(m)*dpdn011(nfou-m)+
     &   dpdn011(m)*drdn001(nfou-m)+drdn011(m)*dphidn001(nfou-m))
        
c        if(nfou.eq.1)then
c         write(6,*)' nfou, dppeta102(nfou) ',nfou, dppeta102(nfou)
c        endif
         
        dppeta102(nfou)=dppeta102(nfou)-0.5d0*(
     &   dppdnh001(m)*dnh101(nfou-m)+dnh001(m)*dppdnh101(nfou-m)+
     &   dppdnh101(m)*dnh001(nfou-m)+dnh101(m)*dppdnh001(nfou-m)+
     &   dppdn001(m)*dn101(nfou-m)+dn001(m)*dppdn101(nfou-m)+
     &   dppdn101(m)*dn001(nfou-m)+dn101(m)*dppdn001(nfou-m)+
     &   dphidnh001(m)*dpdnh101(nfou-m)+dphidnh001(m)*dpdnh101(nfou-m)+
     &   dpdnh101(m)*dphidnh001(nfou-m)+dpdnh101(m)*dphidnh001(nfou-m)+
     &   dphidn001(m)*dpdn101(nfou-m)+dphidn001(m)*dpdn101(nfou-m)+
     &   dpdn101(m)*dphidn001(nfou-m)+dpdn101(m)*dphidn001(nfou-m))
       
        dppeta012(nfou)=dppeta012(nfou)-0.5d0*(
     &   dppdnh001(m)*dnh011(nfou-m)+dnh001(m)*dppdnh011(nfou-m)+
     &   dppdnh011(m)*dnh001(nfou-m)+dnh011(m)*dppdnh001(nfou-m)+
     &   dppdn001(m)*dn011(nfou-m)+dn001(m)*dppdn011(nfou-m)+
     &   dppdn011(m)*dn001(nfou-m)+dn011(m)*dppdn001(nfou-m)+
     &   dphidnh001(m)*dpdnh011(nfou-m)+dphidnh001(m)*dpdnh011(nfou-m)+
     &   dpdnh011(m)*dphidnh001(nfou-m)+dpdnh011(m)*dphidnh001(nfou-m)+
     &   dphidn001(m)*dpdn011(nfou-m)+dphidn001(m)*dpdn011(nfou-m)+
     &   dpdn011(m)*dphidn001(nfou-m)+dpdn011(m)*dphidn001(nfou-m))
       endif
          
      enddo
      enddo    
      
 
      do nfou=-iord3,iord3
      do m=-iord2,iord2
       if(iabs(nfou-m).le.iord2)then
         dreta003(nfou)=dreta003(nfou)-0.5d0*(
     &   drdnh001(m)*dnh002(nfou-m)+drdnh002(m)*dnh001(nfou-m)+
     &   drdn001(m)*dn002(nfou-m)+drdn002(m)*dn001(nfou-m)+
     &   dnh001(m)*drdnh002(nfou-m)+dnh002(m)*drdnh001(nfou-m)+
     &   dn001(m)*drdn002(nfou-m)+dn002(m)*drdn001(nfou-m))
         
         dpeta003(nfou)=dpeta003(nfou)-0.5d0*(
     &   dphidnh001(m)*dnh002(nfou-m)+dpdnh002(m)*dnh001(nfou-m)+
     &   dphidn001(m)*dn002(nfou-m)+dpdn002(m)*dn001(nfou-m)+
     &   dnh001(m)*dpdnh002(nfou-m)+dnh002(m)*dphidnh001(nfou-m)+
     &   dn001(m)*dpdn002(nfou-m)+dn002(m)*dphidn001(nfou-m))   
       endif
      enddo
      enddo 
           
ctesttesttest
      if(iverbose.eq.1)then
      write(6,*)'------------------------------'
      write(6,*)' drdnh101(1) ',drdnh101(1)
      write(6,*)' drdnh011(1) ',drdnh011(1)
      write(6,*)'------------------------------'      
      write(6,*)' drdn101(1) ',drdn101(1)
      write(6,*)' drdn011(1) ',drdn011(1)           
      write(6,*)'------------------------------'
      endif
      
      return
      end      

c--------------------------------        
      subroutine eta4(r,phi)      
c--------------------------------      
      include 'STGFM.cmn'           

      ds=dsin(phi) 
      dc=dcos(phi)     

      do nfou=-iord2,iord2         
      eta202(nfou)=-dc*dnh102(nfou)+ds*dn102(nfou)         
      eta112(nfou)=-dc*dnh012(nfou)+ds*dn012(nfou)-
     &              ds*dnh102(nfou)-dc*dn102(nfou)  
      eta022(nfou)=-ds*dnh012(nfou)-dc*dn012(nfou)
    
      do m=-iord1,iord1
       if(iabs(nfou-m).le.iord1)then
         eta202(nfou)=eta202(nfou)-0.5d0*(
     &   dnh001(m)*dnh201(nfou-m)+dnh101(m)*dnh101(nfou-m)+
     &   dnh201(m)*dnh001(nfou-m)+
     &   dn001(m)*dn201(nfou-m)+dn101(m)*dn101(nfou-m)+
     &   dn201(m)*dn001(nfou-m))         
         
         eta022(nfou)=eta022(nfou)-0.5d0*(
     &   dnh001(m)*dnh021(nfou-m)+dnh011(m)*dnh011(nfou-m)+
     &   dnh021(m)*dnh001(nfou-m)+
     &   dn001(m)*dn021(nfou-m)+dn011(m)*dn011(nfou-m)+
     &   dn021(m)*dn001(nfou-m))         
                  
         eta112(nfou)=eta112(nfou)-0.5d0*(
     &   dnh001(m)*dnh111(nfou-m)+dnh011(m)*dnh101(nfou-m)+
     &   dnh101(m)*dnh011(nfou-m)+
     &   dnh111(m)*dnh001(nfou-m)+
     &   dn001(m)*dn111(nfou-m)+dn011(m)*dn101(nfou-m)+
     &   dn101(m)*dn011(nfou-m)+
     &   dn111(m)*dn001(nfou-m))
       endif
      enddo
      enddo
      
      return
      end
      
c--------------------------------        
      subroutine deta4(r,phi)      
c--------------------------------      
      include 'STGFM.cmn'           

      ds=dsin(phi) 
      dc=dcos(phi)     

c Ableitungen nach r      
      do nfou=-iord2,iord2          
      dreta202(nfou)=-dc*drdnh102(nfou)+ds*drdn102(nfou)         
      dreta112(nfou)=-dc*drdnh012(nfou)+ds*drdn012(nfou)-
     &              ds*drdnh102(nfou)-dc*drdn102(nfou)  
      dreta022(nfou)=-ds*drdnh012(nfou)-dc*drdn012(nfou)
      
c      if(nfou.eq.0)then
c          write(6,*)' dreta202(0) ini ',dreta202(0)  
c          write(6,*)' term1 ',-dc*drdnh102(nfou)
c          write(6,*)' term2',ds*drdn102(nfou)
c      endif
    
      do m=-iord1,iord1
       if(iabs(nfou-m).le.iord1)then
         dreta202(nfou)=dreta202(nfou)-0.5d0*(
     &   drdnh001(m)*dnh201(nfou-m)+drdnh101(m)*dnh101(nfou-m)+
     &   drdnh201(m)*dnh001(nfou-m)+
     &   drdn001(m)*dn201(nfou-m)+drdn101(m)*dn101(nfou-m)+
     &   drdn201(m)*dn001(nfou-m)+
     &   dnh001(m)*drdnh201(nfou-m)+dnh101(m)*drdnh101(nfou-m)+
     &   dnh201(m)*drdnh001(nfou-m)+
     &   dn001(m)*drdn201(nfou-m)+dn101(m)*drdn101(nfou-m)+
     &   dn201(m)*drdn001(nfou-m))
 
         dreta022(nfou)=dreta022(nfou)-0.5d0*(
     &   drdnh001(m)*dnh021(nfou-m)+drdnh011(m)*dnh011(nfou-m)+
     &   drdnh021(m)*dnh001(nfou-m)+
     &   drdn001(m)*dn021(nfou-m)+drdn011(m)*dn011(nfou-m)+
     &   drdn021(m)*dn001(nfou-m)+   
     &   dnh001(m)*drdnh021(nfou-m)+dnh011(m)*drdnh011(nfou-m)+
     &   dnh021(m)*drdnh001(nfou-m)+
     &   dn001(m)*drdn021(nfou-m)+dn011(m)*drdn011(nfou-m)+
     &   dn021(m)*drdn001(nfou-m))                
         
         dreta112(nfou)=dreta112(nfou)-0.5d0*(
     &   drdnh001(m)*dnh111(nfou-m)+drdnh011(m)*dnh101(nfou-m)+
     &   drdnh101(m)*dnh011(nfou-m)+drdnh111(m)*dnh001(nfou-m)+
     &   drdn001(m)*dn111(nfou-m)+drdn011(m)*dn101(nfou-m)+
     &   drdn101(m)*dn011(nfou-m)+drdn111(m)*dn001(nfou-m)+
     &   dnh001(m)*drdnh111(nfou-m)+dnh011(m)*drdnh101(nfou-m)+
     &   dnh101(m)*drdnh011(nfou-m)+dnh111(m)*drdnh001(nfou-m)+
     &   dn001(m)*drdn111(nfou-m)+dn011(m)*drdn101(nfou-m)+
     &   dn101(m)*drdn011(nfou-m)+dn111(m)*drdn001(nfou-m))         
       endif
      enddo
      enddo
    
c      write(6,*)' dreta202(0) ',dreta202(0)

c------ Ableitungen nach phi
      do nfou=-iord2,iord2          
      dpeta202(nfou)=ds*dnh102(nfou)+dc*dn102(nfou)-
     & dc*dpdnh102(nfou)+ds*dpdn102(nfou)
      
c      if(nfou.eq.0)then
c          write(6,*)' dpeta202(0)/r ini ',dpeta202(0)/r     
c          write(6,*)' term1/r ',ds*dnh102(nfou)/r
c          write(6,*)' term2/r',dc*dn102(nfou)/r
c          write(6,*)' term3/r ',-dc*dpdnh102(nfou)/r
c          write(6,*)' term4/r ',ds*dpdn102(nfou)/r
c      endif
      
      dpeta112(nfou)=ds*dnh012(nfou)+dc*dn012(nfou)-
     &               dc*dnh102(nfou)+ds*dn102(nfou)-
     &               dc*dpdnh012(nfou)+ds*dpdn012(nfou)-
     &               ds*dpdnh102(nfou)-dc*dpdn102(nfou)
      dpeta022(nfou)=-dc*dnh012(nfou)+ds*dn012(nfou)-
     &                ds*dpdnh012(nfou)-dc*dpdn012(nfou) 
    
      do m=-iord1,iord1
       if(iabs(nfou-m).le.iord1)then
         dpeta202(nfou)=dpeta202(nfou)-0.5d0*(
     &   dphidnh001(m)*dnh201(nfou-m)+dpdnh101(m)*dnh101(nfou-m)+
     &   dpdnh201(m)*dnh001(nfou-m)+
     &   dphidn001(m)*dn201(nfou-m)+dpdn101(m)*dn101(nfou-m)+
     &   dpdn201(m)*dn001(nfou-m)+
     &   dnh001(m)*dpdnh201(nfou-m)+dnh101(m)*dpdnh101(nfou-m)+
     &   dnh201(m)*dphidnh001(nfou-m)+
     &   dn001(m)*dpdn201(nfou-m)+dn101(m)*dpdn101(nfou-m)+
     &   dn201(m)*dphidn001(nfou-m))
 
         dpeta022(nfou)=dpeta022(nfou)-0.5d0*(
     &   dphidnh001(m)*dnh021(nfou-m)+dpdnh011(m)*dnh011(nfou-m)+
     &   dpdnh021(m)*dnh001(nfou-m)+
     &   dphidn001(m)*dn021(nfou-m)+dpdn011(m)*dn011(nfou-m)+
     &   dpdn021(m)*dn001(nfou-m)+   
     &   dnh001(m)*dpdnh021(nfou-m)+dnh011(m)*dpdnh011(nfou-m)+
     &   dnh021(m)*dphidnh001(nfou-m)+
     &   dn001(m)*dpdn021(nfou-m)+dn011(m)*dpdn011(nfou-m)+
     &   dn021(m)*dphidn001(nfou-m))                
         
         dpeta112(nfou)=dpeta112(nfou)-0.5d0*(
     &   dphidnh001(m)*dnh111(nfou-m)+dpdnh011(m)*dnh101(nfou-m)+
     &   dpdnh101(m)*dnh011(nfou-m)+dpdnh111(m)*dnh001(nfou-m)+
     &   dphidn001(m)*dn111(nfou-m)+dpdn011(m)*dn101(nfou-m)+
     &   dpdn101(m)*dn011(nfou-m)+dpdn111(m)*dn001(nfou-m)+
     &   dnh001(m)*dpdnh111(nfou-m)+dnh011(m)*dpdnh101(nfou-m)+
     &   dnh101(m)*dpdnh011(nfou-m)+dnh111(m)*dphidnh001(nfou-m)+
     &   dn001(m)*dpdn111(nfou-m)+dn011(m)*dpdn101(nfou-m)+
     &   dn101(m)*dpdn011(nfou-m)+dn111(m)*dphidn001(nfou-m))  
       endif
      enddo
      enddo
      
c      write(6,*)' dpeta112(1)/r komplett ',dpeta112(1)/r      

      return
      end
            
c-------------------------------------------
      subroutine ddndnh2(r,phi)
c-------------------------------------------     
      include 'STGFM.cmn'      
      complex*16 tmp
      real*8 small_loc1

      small_loc1=1.d-12
      
      r3i=1.d0/r**3
      r2i=1.d0/r**2
      ri=1.d0/r

c--------- einfache Ableitungen nach phi 

      do nfou=-iord1,-1
      dpdnh011(nfou)=(xi/xkxn(nfou))*
     &  (cdexp(-xi*xkxn(nfou)*s1)*dpreta011(0)-dpreta011(nfou))
      dpdnh011(0)=dpdnh011(0)-dpdnh011(nfou)*cdexp(xi*xkxn(nfou)*zf)
      enddo
      
      do nfou=1,iord1
      dpdnh011(nfou)=(xi/xkxn(nfou))*
     &  (cdexp(-xi*xkxn(nfou)*s1)*dpreta011(0)-dpreta011(nfou))
      dpdnh011(0)=dpdnh011(0)-dpdnh011(nfou)*cdexp(xi*xkxn(nfou)*zf)
      enddo

      do nfou=-iord1,-1
      dpdnh101(nfou)=(xi/xkxn(nfou))*
     &  (cdexp(-xi*xkxn(nfou)*s1)*dpreta101(0)-dpreta101(nfou))
      dpdnh101(0)=dpdnh101(0)-dpdnh101(nfou)*cdexp(xi*xkxn(nfou)*zf)
      enddo
      
      do nfou=1,iord1
      dpdnh101(nfou)=(xi/xkxn(nfou))*
     &  (cdexp(-xi*xkxn(nfou)*s1)*dpreta101(0)-dpreta101(nfou))
      dpdnh101(0)=dpdnh101(0)-dpdnh101(nfou)*cdexp(xi*xkxn(nfou)*zf)
      enddo
  
      do nfou=-iord2,-1
      dpdnh002(nfou)=(xi/xkxn(nfou))*
     &  (cdexp(-xi*xkxn(nfou)*s1)*dpreta002(0)-dpreta002(nfou))
      dpdnh002(0)=dpdnh002(0)-dpdnh002(nfou)*cdexp(xi*xkxn(nfou)*zf)
      enddo
      
      do nfou=1,iord2
      dpdnh002(nfou)=(xi/xkxn(nfou))*
     &  (cdexp(-xi*xkxn(nfou)*s1)*dpreta002(0)-dpreta002(nfou))
      dpdnh002(0)=dpdnh002(0)-dpdnh002(nfou)*cdexp(xi*xkxn(nfou)*zf)
      enddo     
c---------------------------------

      do nfou=-iord1,-1
      dpdn011(nfou)=(1.d0/r)*(xi/xkxn(nfou))*
     &  (cdexp(-xi*xkxn(nfou)*s1)*dppeta011(0)-dppeta011(nfou))
      dpdn011(0)=dpdn011(0)-dpdn011(nfou)*cdexp(xi*xkxn(nfou)*zf)
      enddo
      
      do nfou=1,iord1
      dpdn011(nfou)=(1.d0/r)*(xi/xkxn(nfou))*
     &  (cdexp(-xi*xkxn(nfou)*s1)*dppeta011(0)-dppeta011(nfou))
      dpdn011(0)=dpdn011(0)-dpdn011(nfou)*cdexp(xi*xkxn(nfou)*zf)
      enddo
  
      do nfou=-iord1,-1
      dpdn101(nfou)=(1.d0/r)*(xi/xkxn(nfou))*
     &  (cdexp(-xi*xkxn(nfou)*s1)*dppeta101(0)-dppeta101(nfou))
      dpdn101(0)=dpdn101(0)-dpdn101(nfou)*cdexp(xi*xkxn(nfou)*zf)
      enddo
      
      do nfou=1,iord1
      dpdn101(nfou)=(1.d0/r)*(xi/xkxn(nfou))*
     &  (cdexp(-xi*xkxn(nfou)*s1)*dppeta101(0)-dppeta101(nfou))
      dpdn101(0)=dpdn101(0)-dpdn101(nfou)*cdexp(xi*xkxn(nfou)*zf)
      enddo

      do nfou=-iord2,-1
      dpdn002(nfou)=(1.d0/r)*(xi/xkxn(nfou))*
     &  (cdexp(-xi*xkxn(nfou)*s1)*dppeta002(0)-dppeta002(nfou))
      dpdn002(0)=dpdn002(0)-dpdn002(nfou)*cdexp(xi*xkxn(nfou)*zf)
      enddo
      
      do nfou=1,iord2
      dpdn002(nfou)=(1.d0/r)*(xi/xkxn(nfou))*
     &  (cdexp(-xi*xkxn(nfou)*s1)*dppeta002(0)-dppeta002(nfou))
      dpdn002(0)=dpdn002(0)-dpdn002(nfou)*cdexp(xi*xkxn(nfou)*zf)
      enddo
        
c------ einfache Ableitungen nach r
      do nfou=-iord1,-1
      drdnh011(nfou)=(xi/xkxn(nfou))*
     &  (cdexp(-xi*xkxn(nfou)*s1)*drreta011(0)-drreta011(nfou))     
      drdnh011(0)=drdnh011(0)-drdnh011(nfou)*cdexp(xi*xkxn(nfou)*zf)
      enddo
      
      do nfou=1,iord1
      drdnh011(nfou)=(xi/xkxn(nfou))*
     &  (cdexp(-xi*xkxn(nfou)*s1)*drreta011(0)-drreta011(nfou))    
      drdnh011(0)=drdnh011(0)-drdnh011(nfou)*cdexp(xi*xkxn(nfou)*zf)
      enddo

      do nfou=-iord1,-1
      drdnh101(nfou)=(xi/xkxn(nfou))*
     &  (cdexp(-xi*xkxn(nfou)*s1)*drreta101(0)-drreta101(nfou))
      drdnh101(0)=drdnh101(0)-drdnh101(nfou)*cdexp(xi*xkxn(nfou)*zf)
      enddo
      
      do nfou=1,iord1
      drdnh101(nfou)=(xi/xkxn(nfou))*
     &  (cdexp(-xi*xkxn(nfou)*s1)*drreta101(0)-drreta101(nfou))
      drdnh101(0)=drdnh101(0)-drdnh101(nfou)*cdexp(xi*xkxn(nfou)*zf)
      enddo
      
      do nfou=-iord2,-1
      drdnh002(nfou)=(xi/xkxn(nfou))*
     &  (cdexp(-xi*xkxn(nfou)*s1)*drreta002(0)-drreta002(nfou))     
      drdnh002(0)=drdnh002(0)-drdnh002(nfou)*cdexp(xi*xkxn(nfou)*zf)
      enddo
      
      do nfou=1,iord2
      drdnh002(nfou)=(xi/xkxn(nfou))*
     &  (cdexp(-xi*xkxn(nfou)*s1)*drreta002(0)-drreta002(nfou))    
      drdnh002(0)=drdnh002(0)-drdnh002(nfou)*cdexp(xi*xkxn(nfou)*zf)
      enddo      
c---------------------------------

      do nfou=-iord1,-1
      drdn011(nfou)=(xi/xkxn(nfou))*(
     & -r2i*(cdexp(-xi*xkxn(nfou)*s1)*dphieta011(0)-dphieta011(nfou))+
     &  ri*(cdexp(-xi*xkxn(nfou)*s1)*dpreta011(0)-dpreta011(nfou)))     
      drdn011(0)=drdn011(0)-drdn011(nfou)*cdexp(xi*xkxn(nfou)*zf)
      enddo
       
      do nfou=1,iord1
      drdn011(nfou)=(xi/xkxn(nfou))*(
     & -r2i*(cdexp(-xi*xkxn(nfou)*s1)*dphieta011(0)-dphieta011(nfou))+
     &  ri*(cdexp(-xi*xkxn(nfou)*s1)*dpreta011(0)-dpreta011(nfou)))     
      drdn011(0)=drdn011(0)-drdn011(nfou)*cdexp(xi*xkxn(nfou)*zf)
      enddo

      do nfou=-iord1,-1
      drdn101(nfou)=(xi/xkxn(nfou))*(
     & -r2i*(cdexp(-xi*xkxn(nfou)*s1)*dphieta101(0)-dphieta101(nfou))+
     &  ri*(cdexp(-xi*xkxn(nfou)*s1)*dpreta101(0)-dpreta101(nfou))) 
      drdn101(0)=drdn101(0)-drdn101(nfou)*cdexp(xi*xkxn(nfou)*zf)
      enddo
      
      do nfou=1,iord1
      drdn101(nfou)=(xi/xkxn(nfou))*(
     & -r2i*(cdexp(-xi*xkxn(nfou)*s1)*dphieta101(0)-dphieta101(nfou))+
     &  ri*(cdexp(-xi*xkxn(nfou)*s1)*dpreta101(0)-dpreta101(nfou))) 
      drdn101(0)=drdn101(0)-drdn101(nfou)*cdexp(xi*xkxn(nfou)*zf)
      enddo

      do nfou=-iord2,-1
      drdn002(nfou)=(xi/xkxn(nfou))*(
     & -r2i*(cdexp(-xi*xkxn(nfou)*s1)*dphieta002(0)-dphieta002(nfou))+
     &  ri*(cdexp(-xi*xkxn(nfou)*s1)*dpreta002(0)-dpreta002(nfou)))     
      drdn002(0)=drdn002(0)-drdn002(nfou)*cdexp(xi*xkxn(nfou)*zf)
      enddo
      
      do nfou=1,iord2
      drdn002(nfou)=(xi/xkxn(nfou))*(
     & -r2i*(cdexp(-xi*xkxn(nfou)*s1)*dphieta002(0)-dphieta002(nfou))+
     &  ri*(cdexp(-xi*xkxn(nfou)*s1)*dpreta002(0)-dpreta002(nfou)))     
      drdn002(0)=drdn002(0)-drdn002(nfou)*cdexp(xi*xkxn(nfou)*zf)
      enddo
      
c--------- Ableitung nach r     
c--------- einfache Ableitungen nach phi 

      do nfou=-iord1,-1
      dprdnh011(nfou)=(xi/xkxn(nfou))*
     &  (cdexp(-xi*xkxn(nfou)*s1)*dprreta011(0)-dprreta011(nfou))
      dprdnh011(0)=dprdnh011(0)-dprdnh011(nfou)*cdexp(xi*xkxn(nfou)*zf)
      enddo
      
      do nfou=1,iord1
      dprdnh011(nfou)=(xi/xkxn(nfou))*
     &  (cdexp(-xi*xkxn(nfou)*s1)*dprreta011(0)-dprreta011(nfou))
      dprdnh011(0)=dprdnh011(0)-dprdnh011(nfou)*cdexp(xi*xkxn(nfou)*zf)
      enddo

      do nfou=-iord1,-1
      dprdnh101(nfou)=(xi/xkxn(nfou))*
     &  (cdexp(-xi*xkxn(nfou)*s1)*dprreta101(0)-dprreta101(nfou))
      dprdnh101(0)=dprdnh101(0)-dprdnh101(nfou)*cdexp(xi*xkxn(nfou)*zf)
      enddo
      
      do nfou=1,iord1
      dprdnh101(nfou)=(xi/xkxn(nfou))*
     &  (cdexp(-xi*xkxn(nfou)*s1)*dprreta101(0)-dprreta101(nfou))
      dprdnh101(0)=dprdnh101(0)-dprdnh101(nfou)*cdexp(xi*xkxn(nfou)*zf)
      enddo
  
      do nfou=-iord2,-1
      dprdnh002(nfou)=(xi/xkxn(nfou))*
     &  (cdexp(-xi*xkxn(nfou)*s1)*dprreta002(0)-dprreta002(nfou))
      dprdnh002(0)=dprdnh002(0)-dprdnh002(nfou)*cdexp(xi*xkxn(nfou)*zf)
      enddo
      
      do nfou=1,iord2
      dprdnh002(nfou)=(xi/xkxn(nfou))*
     &  (cdexp(-xi*xkxn(nfou)*s1)*dprreta002(0)-dprreta002(nfou))
      dprdnh002(0)=dprdnh002(0)-dprdnh002(nfou)*cdexp(xi*xkxn(nfou)*zf)
      enddo     
c---------------------------------

      do nfou=-iord1,-1
      dprdn011(nfou)=(-r2i)*(xi/xkxn(nfou))*
     &  (cdexp(-xi*xkxn(nfou)*s1)*dppeta011(0)-dppeta011(nfou))+
     &              (1.d0/r)*(xi/xkxn(nfou))*
     &  (cdexp(-xi*xkxn(nfou)*s1)*dppreta011(0)-dppreta011(nfou))      
      dprdn011(0)=dprdn011(0)-dprdn011(nfou)*cdexp(xi*xkxn(nfou)*zf)
      enddo
      
      do nfou=1,iord1
      dprdn011(nfou)=(-r2i)*(xi/xkxn(nfou))*
     &  (cdexp(-xi*xkxn(nfou)*s1)*dppeta011(0)-dppeta011(nfou))+
     &              (1.d0/r)*(xi/xkxn(nfou))*
     &  (cdexp(-xi*xkxn(nfou)*s1)*dppreta011(0)-dppreta011(nfou))      
      dprdn011(0)=dprdn011(0)-dprdn011(nfou)*cdexp(xi*xkxn(nfou)*zf)
      enddo
      
      do nfou=-iord1,-1
      dprdn101(nfou)=(-r2i)*(xi/xkxn(nfou))*
     &  (cdexp(-xi*xkxn(nfou)*s1)*dppeta101(0)-dppeta101(nfou))+
     &              (1.d0/r)*(xi/xkxn(nfou))*
     &  (cdexp(-xi*xkxn(nfou)*s1)*dppreta101(0)-dppreta101(nfou))      
      dprdn101(0)=dprdn101(0)-dprdn101(nfou)*cdexp(xi*xkxn(nfou)*zf)
      enddo
      
      do nfou=1,iord1
      dprdn101(nfou)=(-r2i)*(xi/xkxn(nfou))*
     &  (cdexp(-xi*xkxn(nfou)*s1)*dppeta101(0)-dppeta101(nfou))+
     &              (1.d0/r)*(xi/xkxn(nfou))*
     &  (cdexp(-xi*xkxn(nfou)*s1)*dppreta101(0)-dppreta101(nfou))      
      dprdn101(0)=dprdn101(0)-dprdn101(nfou)*cdexp(xi*xkxn(nfou)*zf)
      enddo
           
      do nfou=-iord2,-1
      dprdn002(nfou)=(-r2i)*(xi/xkxn(nfou))*
     &  (cdexp(-xi*xkxn(nfou)*s1)*dppeta002(0)-dppeta002(nfou))+
     &              (1.d0/r)*(xi/xkxn(nfou))*
     &  (cdexp(-xi*xkxn(nfou)*s1)*dppreta002(0)-dppreta002(nfou))      
      dprdn002(0)=dprdn002(0)-dprdn002(nfou)*cdexp(xi*xkxn(nfou)*zf)
      enddo
 
      do nfou=1,iord2
      dprdn002(nfou)=(-r2i)*(xi/xkxn(nfou))*
     &  (cdexp(-xi*xkxn(nfou)*s1)*dppeta002(0)-dppeta002(nfou))+
     &              (1.d0/r)*(xi/xkxn(nfou))*
     &  (cdexp(-xi*xkxn(nfou)*s1)*dppreta002(0)-dppreta002(nfou))      
      dprdn002(0)=dprdn002(0)-dprdn002(nfou)*cdexp(xi*xkxn(nfou)*zf)
      enddo           
             
c------ einfache Ableitungen nach r

      do nfou=-iord1,-1
      drrdnh011(nfou)=(xi/xkxn(nfou))*
     &  (cdexp(-xi*xkxn(nfou)*s1)*drrreta011(0)-drrreta011(nfou))     
      drrdnh011(0)=drrdnh011(0)-drrdnh011(nfou)*cdexp(xi*xkxn(nfou)*zf)
      enddo
      
      do nfou=1,iord1
      drrdnh011(nfou)=(xi/xkxn(nfou))*
     &  (cdexp(-xi*xkxn(nfou)*s1)*drrreta011(0)-drrreta011(nfou))     
      drrdnh011(0)=drrdnh011(0)-drrdnh011(nfou)*cdexp(xi*xkxn(nfou)*zf)
      enddo

      do nfou=-iord1,-1
      drrdnh101(nfou)=(xi/xkxn(nfou))*
     &  (cdexp(-xi*xkxn(nfou)*s1)*drrreta101(0)-drrreta101(nfou))
      drrdnh101(0)=drrdnh101(0)-drrdnh101(nfou)*cdexp(xi*xkxn(nfou)*zf)
      enddo
      
      do nfou=1,iord1
      drrdnh101(nfou)=(xi/xkxn(nfou))*
     &  (cdexp(-xi*xkxn(nfou)*s1)*drrreta101(0)-drrreta101(nfou))
      drrdnh101(0)=drrdnh101(0)-drrdnh101(nfou)*cdexp(xi*xkxn(nfou)*zf)
      enddo
      
      do nfou=-iord2,-1
      drrdnh002(nfou)=(xi/xkxn(nfou))*
     &  (cdexp(-xi*xkxn(nfou)*s1)*drrreta002(0)-drrreta002(nfou))     
      drrdnh002(0)=drrdnh002(0)-drrdnh002(nfou)*cdexp(xi*xkxn(nfou)*zf)
      enddo
      
      do nfou=1,iord2
      drrdnh002(nfou)=(xi/xkxn(nfou))*
     &  (cdexp(-xi*xkxn(nfou)*s1)*drrreta002(0)-drrreta002(nfou))     
      drrdnh002(0)=drrdnh002(0)-drrdnh002(nfou)*cdexp(xi*xkxn(nfou)*zf)
      enddo      
c---------------------------------

      do nfou=-iord1,-1
      tmp=-ri*(cdexp(-xi*xkxn(nfou)*s1)*dpreta011(0)-dpreta011(nfou))+
     &  (cdexp(-xi*xkxn(nfou)*s1)*dprreta011(0)-dprreta011(nfou))+
     &  2.d0*r2i*(cdexp(-xi*xkxn(nfou)*s1)*
     &  dphieta011(0)-dphieta011(nfou))-
     &  ri*(cdexp(-xi*xkxn(nfou)*s1)*dpreta011(0)-dpreta011(nfou))
       if(cdabs(tmp).lt.small_loc1)tmp=0.d0     
      drrdn011(nfou)=(xi/xkxn(nfou))*tmp*ri
      enddo

      do nfou=1,iord1
       tmp=-ri*(cdexp(-xi*xkxn(nfou)*s1)*dpreta011(0)-dpreta011(nfou))+
     &  (cdexp(-xi*xkxn(nfou)*s1)*dprreta011(0)-dprreta011(nfou))+
     &  2.d0*r2i*(cdexp(-xi*xkxn(nfou)*s1)*
     &  dphieta011(0)-dphieta011(nfou))-
     &  ri*(cdexp(-xi*xkxn(nfou)*s1)*dpreta011(0)-dpreta011(nfou))
       if(cdabs(tmp).lt.small_loc1)tmp=0.d0     
      drrdn011(nfou)=(xi/xkxn(nfou))*tmp*ri
      enddo
  
      do nfou=-iord1,-1
      tmp=-ri*(cdexp(-xi*xkxn(nfou)*s1)*dpreta101(0)-dpreta101(nfou))+
     &  (cdexp(-xi*xkxn(nfou)*s1)*dprreta101(0)-dprreta101(nfou))+
     &  2.d0*r2i*(cdexp(-xi*xkxn(nfou)*s1)*
     &  dphieta101(0)-dphieta101(nfou))-
     &  ri*(cdexp(-xi*xkxn(nfou)*s1)*dpreta101(0)-dpreta101(nfou))
      drrdn101(nfou)=(xi/xkxn(nfou))*tmp*ri
      drrdn101(0)=drrdn101(0)-drrdn101(nfou)*cdexp(xi*xkxn(nfou)*zf)
      enddo      
      
      do nfou=1,iord1
      tmp=-ri*(cdexp(-xi*xkxn(nfou)*s1)*dpreta101(0)-dpreta101(nfou))+
     &  (cdexp(-xi*xkxn(nfou)*s1)*dprreta101(0)-dprreta101(nfou))+
     &  2.d0*r2i*(cdexp(-xi*xkxn(nfou)*s1)*
     &  dphieta101(0)-dphieta101(nfou))-
     &  ri*(cdexp(-xi*xkxn(nfou)*s1)*dpreta101(0)-dpreta101(nfou))
      drrdn101(nfou)=(xi/xkxn(nfou))*tmp*ri
      drrdn101(0)=drrdn101(0)-drrdn101(nfou)*cdexp(xi*xkxn(nfou)*zf)
      enddo
      
      do nfou=-iord2,-1
      tmp=-ri*(cdexp(-xi*xkxn(nfou)*s1)*dpreta002(0)-dpreta002(nfou))+
     &  (cdexp(-xi*xkxn(nfou)*s1)*dprreta002(0)-dprreta002(nfou))+
     &  2.d0*r2i*(cdexp(-xi*xkxn(nfou)*s1)*
     &  dphieta002(0)-dphieta002(nfou))-
     &  ri*(cdexp(-xi*xkxn(nfou)*s1)*dpreta002(0)-dpreta002(nfou))
      drrdn002(nfou)=(xi/xkxn(nfou))*tmp*ri
      drrdn002(0)=drrdn002(0)-drrdn002(nfou)*cdexp(xi*xkxn(nfou)*zf)
      enddo
      
      do nfou=1,iord2
      tmp=-ri*(cdexp(-xi*xkxn(nfou)*s1)*dpreta002(0)-dpreta002(nfou))+
     &  (cdexp(-xi*xkxn(nfou)*s1)*dprreta002(0)-dprreta002(nfou))+
     &  2.d0*r2i*(cdexp(-xi*xkxn(nfou)*s1)*
     &  dphieta002(0)-dphieta002(nfou))-
     &  ri*(cdexp(-xi*xkxn(nfou)*s1)*dpreta002(0)-dpreta002(nfou))
      drrdn002(nfou)=(xi/xkxn(nfou))*tmp*ri
      drrdn002(0)=drrdn002(0)-drrdn002(nfou)*cdexp(xi*xkxn(nfou)*zf)
      enddo     
   
c--------- Ableitung nach phi
c--------- einfache Ableitungen nach phi 

      do nfou=-iord1,-1
      dppdnh011(nfou)=(xi/xkxn(nfou))*
     &  (cdexp(-xi*xkxn(nfou)*s1)*dppreta011(0)-dppreta011(nfou))
      dppdnh011(0)=dppdnh011(0)-dppdnh011(nfou)*cdexp(xi*xkxn(nfou)*zf)
      enddo
      
      do nfou=1,iord1
      dppdnh011(nfou)=(xi/xkxn(nfou))*
     &  (cdexp(-xi*xkxn(nfou)*s1)*dppreta011(0)-dppreta011(nfou))
      dppdnh011(0)=dppdnh011(0)-dppdnh011(nfou)*cdexp(xi*xkxn(nfou)*zf)
      enddo

      do nfou=-iord1,-1
      dppdnh101(nfou)=(xi/xkxn(nfou))*
     &  (cdexp(-xi*xkxn(nfou)*s1)*dppreta101(0)-dppreta101(nfou))
      dppdnh101(0)=dppdnh101(0)-dppdnh101(nfou)*cdexp(xi*xkxn(nfou)*zf)
      enddo
      
      do nfou=1,iord1
       dppdnh101(nfou)=(xi/xkxn(nfou))*
     &  (cdexp(-xi*xkxn(nfou)*s1)*dppreta101(0)-dppreta101(nfou))
      dppdnh101(0)=dppdnh101(0)-dppdnh101(nfou)*cdexp(xi*xkxn(nfou)*zf)
      enddo
  
      do nfou=-iord2,-1
      dppdnh002(nfou)=(xi/xkxn(nfou))*
     &  (cdexp(-xi*xkxn(nfou)*s1)*dppreta002(0)-dppreta002(nfou))
      dppdnh002(0)=dppdnh002(0)-dppdnh002(nfou)*cdexp(xi*xkxn(nfou)*zf)
      enddo
      
      do nfou=1,iord2
      dppdnh002(nfou)=(xi/xkxn(nfou))*
     &  (cdexp(-xi*xkxn(nfou)*s1)*dppreta002(0)-dppreta002(nfou))
      dppdnh002(0)=dppdnh002(0)-dppdnh002(nfou)*cdexp(xi*xkxn(nfou)*zf)
      enddo     
c---------------------------------

      do nfou=-iord1,-1
      dppdn011(nfou)=(1.d0/r)*(xi/xkxn(nfou))*
     &  (cdexp(-xi*xkxn(nfou)*s1)*dpppeta011(0)-dpppeta011(nfou))
      dppdn011(0)=dppdn011(0)-dppdn011(nfou)*cdexp(xi*xkxn(nfou)*zf)
      enddo
      
      do nfou=1,iord1
      dppdn011(nfou)=(1.d0/r)*(xi/xkxn(nfou))*
     &  (cdexp(-xi*xkxn(nfou)*s1)*dpppeta011(0)-dpppeta011(nfou))
      dppdn011(0)=dppdn011(0)-dppdn011(nfou)*cdexp(xi*xkxn(nfou)*zf)
      enddo
  
      do nfou=-iord1,-1
      dppdn101(nfou)=(1.d0/r)*(xi/xkxn(nfou))*
     &  (cdexp(-xi*xkxn(nfou)*s1)*dpppeta101(0)-dpppeta101(nfou))
      dppdn101(0)=dppdn101(0)-dppdn101(nfou)*cdexp(xi*xkxn(nfou)*zf)
      enddo
      
      do nfou=1,iord1
      dppdn101(nfou)=(1.d0/r)*(xi/xkxn(nfou))*
     &  (cdexp(-xi*xkxn(nfou)*s1)*dpppeta101(0)-dpppeta101(nfou))
      dppdn101(0)=dppdn101(0)-dppdn101(nfou)*cdexp(xi*xkxn(nfou)*zf)
      enddo

      do nfou=-iord2,-1
      dppdn002(nfou)=(1.d0/r)*(xi/xkxn(nfou))*
     &  (cdexp(-xi*xkxn(nfou)*s1)*dpppeta002(0)-dpppeta002(nfou))
      dppdn002(0)=dppdn002(0)-dppdn002(nfou)*cdexp(xi*xkxn(nfou)*zf)
      enddo
      
      do nfou=1,iord2
      dppdn002(nfou)=(1.d0/r)*(xi/xkxn(nfou))*
     &  (cdexp(-xi*xkxn(nfou)*s1)*dpppeta002(0)-dpppeta002(nfou))
      dppdn002(0)=dppdn002(0)-dppdn002(nfou)*cdexp(xi*xkxn(nfou)*zf)
      enddo 
      
      return
      end

c-------------------------------------------
      subroutine dndnh3(r,phi)
c-------------------------------------------     
      include 'STGFM.cmn'      

      ri=1.d0/r

      do nfou=-iord1,-1
      dnh201(nfou)=(xi/xkxn(nfou))*
     &  (cdexp(-xi*xkxn(nfou)*s1)*dreta201(0)-dreta201(nfou))
      dnh201(0)=dnh201(0)-dnh201(nfou)*cdexp(xi*xkxn(nfou)*zf)
      enddo
      
      do nfou=1,iord1
      dnh201(nfou)=(xi/xkxn(nfou))*
     &  (cdexp(-xi*xkxn(nfou)*s1)*dreta201(0)-dreta201(nfou))
      dnh201(0)=dnh201(0)-dnh201(nfou)*cdexp(xi*xkxn(nfou)*zf)
      enddo   
c--
      do nfou=-iord1,-1
      dnh021(nfou)=(xi/xkxn(nfou))*
     &  (cdexp(-xi*xkxn(nfou)*s1)*dreta021(0)-dreta021(nfou))
      dnh021(0)=dnh021(0)-dnh021(nfou)*cdexp(xi*xkxn(nfou)*zf)
      enddo
      
      do nfou=1,iord1
      dnh021(nfou)=(xi/xkxn(nfou))*
     &  (cdexp(-xi*xkxn(nfou)*s1)*dreta021(0)-dreta021(nfou))
      dnh021(0)=dnh021(0)-dnh021(nfou)*cdexp(xi*xkxn(nfou)*zf)
      enddo         
      
c--
      do nfou=-iord1,-1
      dnh111(nfou)=(xi/xkxn(nfou))*
     &  (cdexp(-xi*xkxn(nfou)*s1)*dreta111(0)-dreta111(nfou))
      dnh111(0)=dnh111(0)-dnh111(nfou)*cdexp(xi*xkxn(nfou)*zf)
      enddo
      
      do nfou=1,iord1
      dnh111(nfou)=(xi/xkxn(nfou))*
     &  (cdexp(-xi*xkxn(nfou)*s1)*dreta111(0)-dreta111(nfou))
      dnh111(0)=dnh111(0)-dnh111(nfou)*cdexp(xi*xkxn(nfou)*zf)
      enddo
   
c--
      do nfou=-iord2,-1
      dnh102(nfou)=(xi/xkxn(nfou))*
     &  (cdexp(-xi*xkxn(nfou)*s1)*dreta102(0)-dreta102(nfou))
      dnh102(0)=dnh102(0)-dnh102(nfou)*cdexp(xi*xkxn(nfou)*zf)
      enddo
      
      do nfou=1,iord2
      dnh102(nfou)=(xi/xkxn(nfou))*
     &  (cdexp(-xi*xkxn(nfou)*s1)*dreta102(0)-dreta102(nfou))
      dnh102(0)=dnh102(0)-dnh102(nfou)*cdexp(xi*xkxn(nfou)*zf)
      enddo     
c--
      do nfou=-iord2,-1
      dnh012(nfou)=(xi/xkxn(nfou))*
     &  (cdexp(-xi*xkxn(nfou)*s1)*dreta012(0)-dreta012(nfou))
      dnh012(0)=dnh012(0)-dnh012(nfou)*cdexp(xi*xkxn(nfou)*zf)
      enddo
      
      do nfou=1,iord2
      dnh012(nfou)=(xi/xkxn(nfou))*
     &  (cdexp(-xi*xkxn(nfou)*s1)*dreta012(0)-dreta012(nfou))
      dnh012(0)=dnh012(0)-dnh012(nfou)*cdexp(xi*xkxn(nfou)*zf)
      enddo           
c----------------------------------------------------------------

      do nfou=-iord1,-1
      dn201(nfou)=ri*(xi/xkxn(nfou))*
     &  (cdexp(-xi*xkxn(nfou)*s1)*dpeta201(0)-dpeta201(nfou))
      dn201(0)=dn201(0)-dn201(nfou)*cdexp(xi*xkxn(nfou)*zf)
      enddo
      
      do nfou=1,iord1
      dn201(nfou)=ri*(xi/xkxn(nfou))*
     &  (cdexp(-xi*xkxn(nfou)*s1)*dpeta201(0)-dpeta201(nfou))
      dn201(0)=dn201(0)-dn201(nfou)*cdexp(xi*xkxn(nfou)*zf)
      enddo   
c--
      do nfou=-iord1,-1
      dn021(nfou)=ri*(xi/xkxn(nfou))*
     &  (cdexp(-xi*xkxn(nfou)*s1)*dpeta021(0)-dpeta021(nfou))
      dn021(0)=dn021(0)-dn021(nfou)*cdexp(xi*xkxn(nfou)*zf)
      enddo
      
      do nfou=1,iord1
      dn021(nfou)=ri*(xi/xkxn(nfou))*
     &  (cdexp(-xi*xkxn(nfou)*s1)*dpeta021(0)-dpeta021(nfou))
      dn021(0)=dn021(0)-dn021(nfou)*cdexp(xi*xkxn(nfou)*zf)
      enddo         
      
c--
      do nfou=-iord1,-1
      dn111(nfou)=ri*(xi/xkxn(nfou))*
     &  (cdexp(-xi*xkxn(nfou)*s1)*dpeta111(0)-dpeta111(nfou))
      dn111(0)=dn111(0)-dn111(nfou)*cdexp(xi*xkxn(nfou)*zf)
      enddo
      
      do nfou=1,iord1
      dn111(nfou)=ri*(xi/xkxn(nfou))*
     &  (cdexp(-xi*xkxn(nfou)*s1)*dpeta111(0)-dpeta111(nfou))
      dn111(0)=dn111(0)-dn111(nfou)*cdexp(xi*xkxn(nfou)*zf)
      enddo
 
c--
      do nfou=-iord2,-1
      dn102(nfou)=ri*(xi/xkxn(nfou))*
     &  (cdexp(-xi*xkxn(nfou)*s1)*dpeta102(0)-dpeta102(nfou))
      dn102(0)=dn102(0)-dn102(nfou)*cdexp(xi*xkxn(nfou)*zf)
      enddo
      
      do nfou=1,iord2
      dn102(nfou)=ri*(xi/xkxn(nfou))*
     &  (cdexp(-xi*xkxn(nfou)*s1)*dpeta102(0)-dpeta102(nfou))
      dn102(0)=dn102(0)-dn102(nfou)*cdexp(xi*xkxn(nfou)*zf)
      enddo
                
c--
      do nfou=-iord2,-1
      dn012(nfou)=ri*(xi/xkxn(nfou))*
     &  (cdexp(-xi*xkxn(nfou)*s1)*dpeta012(0)-dpeta012(nfou))
      dn012(0)=dn012(0)-dn012(nfou)*cdexp(xi*xkxn(nfou)*zf)
      enddo
      
      do nfou=1,iord2
      dn012(nfou)=ri*(xi/xkxn(nfou))*
     &  (cdexp(-xi*xkxn(nfou)*s1)*dpeta012(0)-dpeta012(nfou))
      dn012(0)=dn012(0)-dn012(nfou)*cdexp(xi*xkxn(nfou)*zf)
      enddo     
      
      return
      end
  
c-------------------------------------------
      subroutine ddndnh3(r,phi)
c-------------------------------------------     
      include 'STGFM.cmn'      

      ri=1.d0/r
      r2i=1.d0/r**2

c-------- Ableitungen nach phi
      do nfou=-iord1,-1
      dpdnh201(nfou)=(xi/xkxn(nfou))*
     &  (cdexp(-xi*xkxn(nfou)*s1)*dpreta201(0)-dpreta201(nfou))
      dpdnh201(0)=dpdnh201(0)-dpdnh201(nfou)*cdexp(xi*xkxn(nfou)*zf)
      enddo
      
      do nfou=1,iord1
      dpdnh201(nfou)=(xi/xkxn(nfou))*
     &  (cdexp(-xi*xkxn(nfou)*s1)*dpreta201(0)-dpreta201(nfou))
      dpdnh201(0)=dpdnh201(0)-dpdnh201(nfou)*cdexp(xi*xkxn(nfou)*zf)
      enddo   
c--
      do nfou=-iord1,-1
      dpdnh021(nfou)=(xi/xkxn(nfou))*
     &  (cdexp(-xi*xkxn(nfou)*s1)*dpreta021(0)-dpreta021(nfou))
      dpdnh021(0)=dpdnh021(0)-dpdnh021(nfou)*cdexp(xi*xkxn(nfou)*zf)
      enddo
      
      do nfou=1,iord1
      dpdnh021(nfou)=(xi/xkxn(nfou))*
     &  (cdexp(-xi*xkxn(nfou)*s1)*dpreta021(0)-dpreta021(nfou))
      dpdnh021(0)=dpdnh021(0)-dpdnh021(nfou)*cdexp(xi*xkxn(nfou)*zf)
      enddo         
      
c--
      do nfou=-iord1,-1
      dpdnh111(nfou)=(xi/xkxn(nfou))*
     &  (cdexp(-xi*xkxn(nfou)*s1)*dpreta111(0)-dpreta111(nfou))
      dpdnh111(0)=dpdnh111(0)-dpdnh111(nfou)*cdexp(xi*xkxn(nfou)*zf)
      enddo
      
      do nfou=1,iord1
      dpdnh111(nfou)=(xi/xkxn(nfou))*
     &  (cdexp(-xi*xkxn(nfou)*s1)*dpreta111(0)-dpreta111(nfou))
      dpdnh111(0)=dpdnh111(0)-dpdnh111(nfou)*cdexp(xi*xkxn(nfou)*zf)
      enddo
   
c--
      do nfou=-iord2,-1
      dpdnh102(nfou)=(xi/xkxn(nfou))*
     &  (cdexp(-xi*xkxn(nfou)*s1)*dpreta102(0)-dpreta102(nfou))
      dpdnh102(0)=dpdnh102(0)-dpdnh102(nfou)*cdexp(xi*xkxn(nfou)*zf)
      enddo
      
      do nfou=1,iord2
      dpdnh102(nfou)=(xi/xkxn(nfou))*
     &  (cdexp(-xi*xkxn(nfou)*s1)*dpreta102(0)-dpreta102(nfou))
      dpdnh102(0)=dpdnh102(0)-dpdnh102(nfou)*cdexp(xi*xkxn(nfou)*zf)
      enddo     
c--
      do nfou=-iord2,-1
      dpdnh012(nfou)=(xi/xkxn(nfou))*
     &  (cdexp(-xi*xkxn(nfou)*s1)*dpreta012(0)-dpreta012(nfou))
      dpdnh012(0)=dpdnh012(0)-dpdnh012(nfou)*cdexp(xi*xkxn(nfou)*zf)
      enddo
      
      do nfou=1,iord2
      dpdnh012(nfou)=(xi/xkxn(nfou))*
     &  (cdexp(-xi*xkxn(nfou)*s1)*dpreta012(0)-dpreta012(nfou))
      dpdnh012(0)=dpdnh012(0)-dpdnh012(nfou)*cdexp(xi*xkxn(nfou)*zf)
      enddo           
c----------------------------------------------------------------

      do nfou=-iord1,-1
      dpdn201(nfou)=ri*(xi/xkxn(nfou))*
     &  (cdexp(-xi*xkxn(nfou)*s1)*dppeta201(0)-dppeta201(nfou))
      dpdn201(0)=dpdn201(0)-dpdn201(nfou)*cdexp(xi*xkxn(nfou)*zf)
      enddo
      
      do nfou=1,iord1
      dpdn201(nfou)=ri*(xi/xkxn(nfou))*
     &  (cdexp(-xi*xkxn(nfou)*s1)*dppeta201(0)-dppeta201(nfou))
      dpdn201(0)=dpdn201(0)-dpdn201(nfou)*cdexp(xi*xkxn(nfou)*zf)
      enddo   
c--
      do nfou=-iord1,-1
      dpdn021(nfou)=ri*(xi/xkxn(nfou))*
     &  (cdexp(-xi*xkxn(nfou)*s1)*dppeta021(0)-dppeta021(nfou))
      dpdn021(0)=dpdn021(0)-dpdn021(nfou)*cdexp(xi*xkxn(nfou)*zf)
      enddo
      
      do nfou=1,iord1
      dpdn021(nfou)=ri*(xi/xkxn(nfou))*
     &  (cdexp(-xi*xkxn(nfou)*s1)*dppeta021(0)-dppeta021(nfou))
      dpdn021(0)=dpdn021(0)-dpdn021(nfou)*cdexp(xi*xkxn(nfou)*zf)
      enddo         
      
c--
      do nfou=-iord1,-1
      dpdn111(nfou)=ri*(xi/xkxn(nfou))*
     &  (cdexp(-xi*xkxn(nfou)*s1)*dppeta111(0)-dppeta111(nfou))
      dpdn111(0)=dpdn111(0)-dpdn111(nfou)*cdexp(xi*xkxn(nfou)*zf)
      enddo
      
      do nfou=1,iord1
      dpdn111(nfou)=ri*(xi/xkxn(nfou))*
     &  (cdexp(-xi*xkxn(nfou)*s1)*dppeta111(0)-dppeta111(nfou))
      dpdn111(0)=dpdn111(0)-dpdn111(nfou)*cdexp(xi*xkxn(nfou)*zf)
      enddo
 
c--
      do nfou=-iord2,-1
      dpdn102(nfou)=ri*(xi/xkxn(nfou))*
     &  (cdexp(-xi*xkxn(nfou)*s1)*dppeta102(0)-dppeta102(nfou))
      dpdn102(0)=dpdn102(0)-dpdn102(nfou)*cdexp(xi*xkxn(nfou)*zf)    
      enddo
      
      do nfou=1,iord2
      dpdn102(nfou)=ri*(xi/xkxn(nfou))*
     &  (cdexp(-xi*xkxn(nfou)*s1)*dppeta102(0)-dppeta102(nfou))
      dpdn102(0)=dpdn102(0)-dpdn102(nfou)*cdexp(xi*xkxn(nfou)*zf)       
      enddo    
c--
      do nfou=-iord2,-1
      dpdn012(nfou)=ri*(xi/xkxn(nfou))*
     &  (cdexp(-xi*xkxn(nfou)*s1)*dppeta012(0)-dppeta012(nfou))
      dpdn012(0)=dpdn012(0)-dpdn012(nfou)*cdexp(xi*xkxn(nfou)*zf)
      enddo
      
      do nfou=1,iord2
      dpdn012(nfou)=ri*(xi/xkxn(nfou))*
     &  (cdexp(-xi*xkxn(nfou)*s1)*dppeta012(0)-dppeta012(nfou))
      dpdn012(0)=dpdn012(0)-dpdn012(nfou)*cdexp(xi*xkxn(nfou)*zf)
      enddo     
          
c------ Ableitungen nach r
      do nfou=-iord1,-1
      drdnh201(nfou)=(xi/xkxn(nfou))*
     &  (cdexp(-xi*xkxn(nfou)*s1)*drreta201(0)-drreta201(nfou))
      drdnh201(0)=drdnh201(0)-drdnh201(nfou)*cdexp(xi*xkxn(nfou)*zf)
      enddo
      
      do nfou=1,iord1
      drdnh201(nfou)=(xi/xkxn(nfou))*
     &  (cdexp(-xi*xkxn(nfou)*s1)*drreta201(0)-drreta201(nfou))
      drdnh201(0)=drdnh201(0)-drdnh201(nfou)*cdexp(xi*xkxn(nfou)*zf)
      enddo   
c--
      do nfou=-iord1,-1
      drdnh021(nfou)=(xi/xkxn(nfou))*
     &  (cdexp(-xi*xkxn(nfou)*s1)*drreta021(0)-drreta021(nfou))
      drdnh021(0)=drdnh021(0)-drdnh021(nfou)*cdexp(xi*xkxn(nfou)*zf)
      enddo
      
      do nfou=1,iord1
      drdnh021(nfou)=(xi/xkxn(nfou))*
     &  (cdexp(-xi*xkxn(nfou)*s1)*drreta021(0)-drreta021(nfou))
      drdnh021(0)=drdnh021(0)-drdnh021(nfou)*cdexp(xi*xkxn(nfou)*zf)
      enddo         
      
c--
      do nfou=-iord1,-1
      drdnh111(nfou)=(xi/xkxn(nfou))*
     &  (cdexp(-xi*xkxn(nfou)*s1)*drreta111(0)-drreta111(nfou))
      drdnh111(0)=drdnh111(0)-drdnh111(nfou)*cdexp(xi*xkxn(nfou)*zf)
      enddo
      
      do nfou=1,iord1
      drdnh111(nfou)=(xi/xkxn(nfou))*
     &  (cdexp(-xi*xkxn(nfou)*s1)*drreta111(0)-drreta111(nfou))
      drdnh111(0)=drdnh111(0)-drdnh111(nfou)*cdexp(xi*xkxn(nfou)*zf)
      enddo
   
c--
      do nfou=-iord2,-1
      drdnh102(nfou)=(xi/xkxn(nfou))*
     &  (cdexp(-xi*xkxn(nfou)*s1)*drreta102(0)-drreta102(nfou))
      drdnh102(0)=drdnh102(0)-drdnh102(nfou)*cdexp(xi*xkxn(nfou)*zf)
      enddo
      
      do nfou=1,iord2
      drdnh102(nfou)=(xi/xkxn(nfou))*
     &  (cdexp(-xi*xkxn(nfou)*s1)*drreta102(0)-drreta102(nfou))
      drdnh102(0)=drdnh102(0)-drdnh102(nfou)*cdexp(xi*xkxn(nfou)*zf)
      enddo     
c--
      do nfou=-iord2,-1
      drdnh012(nfou)=(xi/xkxn(nfou))*
     &  (cdexp(-xi*xkxn(nfou)*s1)*drreta012(0)-drreta012(nfou))
      drdnh012(0)=drdnh012(0)-drdnh012(nfou)*cdexp(xi*xkxn(nfou)*zf)
      enddo
      
      do nfou=1,iord2
      drdnh012(nfou)=(xi/xkxn(nfou))*
     &  (cdexp(-xi*xkxn(nfou)*s1)*drreta012(0)-drreta012(nfou))
      drdnh012(0)=drdnh012(0)-drdnh012(nfou)*cdexp(xi*xkxn(nfou)*zf)
      enddo           
c----------------------------------------------------------------

      do nfou=-iord1,-1
      drdn201(nfou)=-r2i*(xi/xkxn(nfou))*
     &  (cdexp(-xi*xkxn(nfou)*s1)*dpeta201(0)-dpeta201(nfou))+
     &             ri*(xi/xkxn(nfou))*
     &  (cdexp(-xi*xkxn(nfou)*s1)*dpreta201(0)-dpreta201(nfou))   
      drdn201(0)=drdn201(0)-drdn201(nfou)*cdexp(xi*xkxn(nfou)*zf)
      enddo
      
      do nfou=1,iord1
      drdn201(nfou)=-r2i*(xi/xkxn(nfou))*
     &  (cdexp(-xi*xkxn(nfou)*s1)*dpeta201(0)-dpeta201(nfou))+
     &             ri*(xi/xkxn(nfou))*
     &  (cdexp(-xi*xkxn(nfou)*s1)*dpreta201(0)-dpreta201(nfou))   
      drdn201(0)=drdn201(0)-drdn201(nfou)*cdexp(xi*xkxn(nfou)*zf)
      enddo   
c--
      do nfou=-iord1,-1
       drdn021(nfou)=-r2i*(xi/xkxn(nfou))*
     &  (cdexp(-xi*xkxn(nfou)*s1)*dpeta021(0)-dpeta021(nfou))+
     &             ri*(xi/xkxn(nfou))*
     &  (cdexp(-xi*xkxn(nfou)*s1)*dpreta021(0)-dpreta021(nfou))   
      drdn021(0)=drdn021(0)-drdn021(nfou)*cdexp(xi*xkxn(nfou)*zf)
      enddo
      
      do nfou=1,iord1
      drdn021(nfou)=-r2i*(xi/xkxn(nfou))*
     &  (cdexp(-xi*xkxn(nfou)*s1)*dpeta021(0)-dpeta021(nfou))+
     &             ri*(xi/xkxn(nfou))*
     &  (cdexp(-xi*xkxn(nfou)*s1)*dpreta021(0)-dpreta021(nfou))   
      drdn021(0)=drdn021(0)-drdn021(nfou)*cdexp(xi*xkxn(nfou)*zf)
      enddo         
      
c--
      do nfou=-iord1,-1
       drdn111(nfou)=-r2i*(xi/xkxn(nfou))*
     &  (cdexp(-xi*xkxn(nfou)*s1)*dpeta111(0)-dpeta111(nfou))+
     &             ri*(xi/xkxn(nfou))*
     &  (cdexp(-xi*xkxn(nfou)*s1)*dpreta111(0)-dpreta111(nfou))   
      drdn111(0)=drdn111(0)-drdn111(nfou)*cdexp(xi*xkxn(nfou)*zf)
      enddo
      
      do nfou=1,iord1
      drdn111(nfou)=-r2i*(xi/xkxn(nfou))*
     &  (cdexp(-xi*xkxn(nfou)*s1)*dpeta111(0)-dpeta111(nfou))+
     &             ri*(xi/xkxn(nfou))*
     &  (cdexp(-xi*xkxn(nfou)*s1)*dpreta111(0)-dpreta111(nfou))   
      drdn111(0)=drdn111(0)-drdn111(nfou)*cdexp(xi*xkxn(nfou)*zf)
      enddo
 
c--
      do nfou=-iord2,-1
      drdn102(nfou)=-r2i*(xi/xkxn(nfou))*
     &  (cdexp(-xi*xkxn(nfou)*s1)*dpeta102(0)-dpeta102(nfou))+
     &             ri*(xi/xkxn(nfou))*
     &  (cdexp(-xi*xkxn(nfou)*s1)*dpreta102(0)-dpreta102(nfou))   
      drdn102(0)=drdn102(0)-drdn102(nfou)*cdexp(xi*xkxn(nfou)*zf)
            
      enddo
      
      do nfou=1,iord2
      drdn102(nfou)=-r2i*(xi/xkxn(nfou))*
     &  (cdexp(-xi*xkxn(nfou)*s1)*dpeta102(0)-dpeta102(nfou))+
     &             ri*(xi/xkxn(nfou))*
     &  (cdexp(-xi*xkxn(nfou)*s1)*dpreta102(0)-dpreta102(nfou))   
      drdn102(0)=drdn102(0)-drdn102(nfou)*cdexp(xi*xkxn(nfou)*zf)
      
      enddo
                
c--
      do nfou=-iord2,-1
      drdn012(nfou)=-r2i*(xi/xkxn(nfou))*
     &  (cdexp(-xi*xkxn(nfou)*s1)*dpeta012(0)-dpeta012(nfou))+
     &             ri*(xi/xkxn(nfou))*
     &  (cdexp(-xi*xkxn(nfou)*s1)*dpreta012(0)-dpreta012(nfou))   
      drdn012(0)=drdn012(0)-drdn012(nfou)*cdexp(xi*xkxn(nfou)*zf)
      enddo
      
      do nfou=1,iord2
      drdn012(nfou)=-r2i*(xi/xkxn(nfou))*
     &  (cdexp(-xi*xkxn(nfou)*s1)*dpeta012(0)-dpeta012(nfou))+
     &             ri*(xi/xkxn(nfou))*
     &  (cdexp(-xi*xkxn(nfou)*s1)*dpreta012(0)-dpreta012(nfou))   
      drdn012(0)=drdn012(0)-drdn012(nfou)*cdexp(xi*xkxn(nfou)*zf)
      enddo     
      
      return
      end
            
c-------------------------------------------
      subroutine facexp 
c-------------------------------------------
      include 'STGFM.cmn'      

      facex(0)=0.d0
      
      do nfou=-iord_max,-1
      facex(nfou)=
     &  (cdexp(xi*xkxn(nfou)*z0)-
     &  cdexp(xi*xkxn(nfou)*zf))/(xi*xkxn(nfou))
      facex(-nfou)=
     &  (cdexp(xi*xkxn(-nfou)*z0)-
     &   cdexp(xi*xkxn(-nfou)*zf))/(xi*xkxn(-nfou))
      enddo

      return
      end

c-------------------------------------------
      subroutine gx(r,phi) 
c-------------------------------------------
      include 'STGFM.cmn'

      X011=(dcos(phi)*dreta011(0)-(dsin(phi)/r)*dphieta011(0))*(z0-zf)
      X101=(dcos(phi)*dreta101(0)-(dsin(phi)/r)*dphieta101(0))*(z0-zf)  
      X002=(dcos(phi)*dreta002(0)-(dsin(phi)/r)*dphieta002(0))*(z0-zf)
      
      X102=(dcos(phi)*dreta102(0)-(dsin(phi)/r)*dpeta102(0))*(z0-zf)
      X012=(dcos(phi)*dreta012(0)-(dsin(phi)/r)*dpeta012(0))*(z0-zf)
      X201=(dcos(phi)*dreta201(0)-(dsin(phi)/r)*dpeta201(0))*(z0-zf)  
      X021=(dcos(phi)*dreta021(0)-(dsin(phi)/r)*dpeta021(0))*(z0-zf)
      X003=(dcos(phi)*dreta003(0)-(dsin(phi)/r)*dpeta003(0))*(z0-zf)
      
c      write(6,*)'-----------------'      
c      write(6,*)' x021 ini ',x021
      
      X111=(dcos(phi)*dreta111(0)-(dsin(phi)/r)*dpeta111(0))*(z0-zf)  

      x202=(dcos(phi)*dreta202(0)-(dsin(phi)/r)*dpeta202(0))*(z0-zf)
      x112=(dcos(phi)*dreta112(0)-(dsin(phi)/r)*dpeta112(0))*(z0-zf)
      x022=(dcos(phi)*dreta022(0)-(dsin(phi)/r)*dpeta022(0))*(z0-zf)
      
c      write(6,*)' x202 ini ',x202

      do nfou=-iord1,iord1
      X011=X011+facex(nfou)*
     &  (dcos(phi)*dreta011(nfou)-(dsin(phi)/r)*dphieta011(nfou))
      X101=X101+facex(nfou)*
     &  (dcos(phi)*dreta101(nfou)-(dsin(phi)/r)*dphieta101(nfou))
      X201=X201+facex(nfou)*
     &  (dcos(phi)*dreta201(nfou)-(dsin(phi)/r)*dpeta201(nfou)) 
      X021=X021+facex(nfou)*
     &  (dcos(phi)*dreta021(nfou)-(dsin(phi)/r)*dpeta021(nfou))   
      X111=X111+facex(nfou)*
     &  (dcos(phi)*dreta111(nfou)-(dsin(phi)/r)*dpeta111(nfou)) 
      enddo

      do nfou=-iord2,iord2
      X002=X002+facex(nfou)*
     &  (dcos(phi)*dreta002(nfou)-(dsin(phi)/r)*dphieta002(nfou))      
      X102=X102+facex(nfou)*
     &  (dcos(phi)*dreta102(nfou)-(dsin(phi)/r)*dpeta102(nfou))
      X012=X012+facex(nfou)*
     &  (dcos(phi)*dreta012(nfou)-(dsin(phi)/r)*dpeta012(nfou))        
      
      X202=X202+facex(nfou)*
     &  (dcos(phi)*dreta202(nfou)-(dsin(phi)/r)*dpeta202(nfou))
      X112=X202+facex(nfou)*
     &  (dcos(phi)*dreta112(nfou)-(dsin(phi)/r)*dpeta112(nfou))       
      X022=X202+facex(nfou)*
     &  (dcos(phi)*dreta022(nfou)-(dsin(phi)/r)*dpeta022(nfou))       
      enddo

      do nfou=-iord3,iord3
      X003=X003+facex(nfou)*
     &  (dcos(phi)*dreta003(nfou)-(dsin(phi)/r)*dpeta003(nfou))  
      enddo 

      if(iverbose.eq.1)then
      write(6,*)'-----------------'
      write(6,*)' x002 ',x002
      write(6,*)' x101 ',x101
      write(6,*)' x011 ',x011   
      write(6,*)' x201 ',x201
      write(6,*)' x111 ',x111
      write(6,*)' x021 ',x021
      write(6,*)' x102 ',x102
      write(6,*)' x012 ',x012
      write(6,*)' x003 ',x003
      write(6,*)' x202 ',x202
      write(6,*)' x112 ',x112
      write(6,*)' x022 ',x022
      endif
      
      return
      end
                                                                
c-------------------------------------------
      subroutine gy(r,phi) 
c-------------------------------------------
      include 'STGFM.cmn'      

      Y011=(dsin(phi)*dreta011(0)+(dcos(phi)/r)*dphieta011(0))*(z0-zf)
      Y101=(dsin(phi)*dreta101(0)+(dcos(phi)/r)*dphieta101(0))*(z0-zf)
      Y002=(dsin(phi)*dreta002(0)+(dcos(phi)/r)*dphieta002(0))*(z0-zf)
      
      Y102=(dsin(phi)*dreta102(0)+(dcos(phi)/r)*dpeta102(0))*(z0-zf)
      Y012=(dsin(phi)*dreta012(0)+(dcos(phi)/r)*dpeta012(0))*(z0-zf)
      Y201=(dsin(phi)*dreta201(0)+(dcos(phi)/r)*dpeta201(0))*(z0-zf)
      Y021=(dsin(phi)*dreta021(0)+(dcos(phi)/r)*dpeta021(0))*(z0-zf)  
      Y111=(dsin(phi)*dreta111(0)+(dcos(phi)/r)*dpeta111(0))*(z0-zf)  
      Y003=(dsin(phi)*dreta003(0)+(dcos(phi)/r)*dpeta003(0))*(z0-zf)  

      y202=(dsin(phi)*dreta202(0)+(dcos(phi)/r)*dpeta202(0))*(z0-zf) 
      y112=(dsin(phi)*dreta112(0)+(dcos(phi)/r)*dpeta112(0))*(z0-zf)
      y022=(dsin(phi)*dreta022(0)+(dcos(phi)/r)*dpeta022(0))*(z0-zf)
      
c      write(6,*)' dpeta112 (0) ',dpeta112(0)
      
      do nfou=-iord1,iord1
      Y011=Y011+facex(nfou)*
     &  (dsin(phi)*dreta011(nfou)+(dcos(phi)/r)*dphieta011(nfou))
      Y101=Y101+facex(nfou)*
     &  (dsin(phi)*dreta101(nfou)+(dcos(phi)/r)*dphieta101(nfou))
      Y201=Y201+facex(nfou)*
     &  (dsin(phi)*dreta201(nfou)+(dcos(phi)/r)*dpeta201(nfou))    
      Y021=Y021+facex(nfou)*
     &  (dsin(phi)*dreta021(nfou)+(dcos(phi)/r)*dpeta021(nfou))   
      Y111=Y111+facex(nfou)*
     &  (dsin(phi)*dreta111(nfou)+(dcos(phi)/r)*dpeta111(nfou))
      enddo

      do nfou=-iord2,iord2 
      Y002=Y002+facex(nfou)*
     &  (dsin(phi)*dreta002(nfou)+(dcos(phi)/r)*dphieta002(nfou))      
      Y102=Y102+facex(nfou)*
     &  (dsin(phi)*dreta102(nfou)+(dcos(phi)/r)*dpeta102(nfou))
      Y012=Y012+facex(nfou)*
     &  (dsin(phi)*dreta012(nfou)+(dcos(phi)/r)*dpeta012(nfou))
      Y202=Y202+facex(nfou)*
     &  (dsin(phi)*dreta202(nfou)+(dcos(phi)/r)*dpeta202(nfou))  
      Y112=Y112+facex(nfou)*
     &  (dsin(phi)*dreta112(nfou)+(dcos(phi)/r)*dpeta112(nfou))
      Y022=Y022+facex(nfou)*
     &  (dsin(phi)*dreta022(nfou)+(dcos(phi)/r)*dpeta022(nfou))
      enddo

      do nfou=-iord3,iord3 
      Y003=Y003+facex(nfou)*
     &  (dsin(phi)*dreta003(nfou)+(dcos(phi)/r)*dpeta003(nfou))    
      enddo   

      if(iverbose.eq.1)then
      write(6,*)'-----------------' 
      write(6,*)' y002 ',y002
      write(6,*)' y101 ',y101
      write(6,*)' y011 ',y011   
      write(6,*)' y201 ',y201
      write(6,*)' y111 ',y111
      write(6,*)' y102 ',y102
      write(6,*)' y012 ',y012    
      write(6,*)' y021 ',y021    
      write(6,*)' y003 ',y003     
      write(6,*)' y202 ',y202
      write(6,*)' y112 ',y112
      write(6,*)' y022 ',y022
      write(6,*)'-----------------'  
      
      write(6,*)' dreta012(1) ',dreta012(1) 
      write(6,*)' dpeta012(1) ',dpeta012(1)   
      write(6,*)'-----------------'     
      write(6,*)' dreta102(1) ',dreta102(1) 
      write(6,*)' dpeta102(1) ',dpeta102(1)    
      write(6,*)'-----------------'     

      endif
      
ctesttesttest
c      y002=0.d0
c      y102=0.d0
c      y012=0.d0
      
      
      
      return
      end
  
c-------------------------------------------
      subroutine get_AxAy(phi,z,Ax,Ay)      
c-------------------------------------------
      include 'STGFM.cmn'      
      complex*16 axc,ayc,draxc,dpaxc,dparc
      real*8 drax,dpax,dpar

      Axc=0.d0
      Ayc=0.d0
      
      do nfou=-iord1,iord1
      Axc=Axc-(dcos(phi)*dnh001(nfou)-dsin(phi)*dn001(nfou))*
     &  cdexp(xi*xkxn(nfou)*z)
      Ayc=Ayc-(sin(phi)*dnh001(nfou)+dcos(phi)*dn001(nfou))*
     &  cdexp(xi*xkxn(nfou)*z)    
      
      enddo
      
      Ax=dreal(Axc)
      Ay=dreal(Ayc)
      
      dparc=0.d0
      do nfou=-iord1,iord1
      dpArc=dpArc-drdn001(nfou)*cdexp(xi*xkxn(nfou)*z)
      
      enddo
      
      dpAr=dreal(dpArc)
 
c       write(6,*)' dpAr ',dpAr
               
      return
      end

c-------------------------------------------
      subroutine get_px0py0(xp0,yp0,Ax0,Ay0,px0,py0)      
c-------------------------------------------
      include 'STGFM.cmn'

      px0 = xp0 + Ax0
      py0 = yp0 + Ay0
      
      return
      end
      
c-------------------------------------------
      subroutine get_pxfpyf_newton(px0,py0,pxf,pyf)      
c-------------------------------------------
      include 'STGFM.cmn'      
      
      complex*16 x00n,x10,x01,x20,x11,x02,
     & y00n,y10,y01,y20,y11,y02,
     & g11,g12,g21,g22,hx,hy,
     & f11,f12,f21,f22,ac,cc,detie
      
      real*8 small_loc1,dpxf,dpyf,check_x,check_y
      
      integer imsteps,isteps

      small_loc1=1.d-12
      imsteps=12  ! maximal number of iteration steps until termination
      isteps=0
      
      if(igf.eq.2)then
      x00n=x002
      x10=x101
      x01=x011
      x20=0.d0
      x11=0.d0
      x02=0.d0
      
      y00n=y002    
      y10=y101
      y01=y011
      y20=0.d0
      y11=0.d0
      y02=0.d0
      endif

      if(igf.eq.21)then
      x00n=0.d0
      x10=x101
      x01=x011   
      x20=x201
      x11=x111
      x02=x021
      
      y00n=0.d0
      y10=y101
      y01=y011    
      y20=y201
      y11=y111
      y02=y021  
      endif
      
      if(igf.eq.11)then          
      x00n=0d0
      x10=x101
      x01=x011 
      x20=0.d0
      x11=0.d0
      x02=0.d0
      
      y00n=0.d0
      y10=y101
      y01=y011 
      y20=0.d0
      y11=0.d0
      y02=0.d0  
      endif
      
      if(igf.eq.12)then          
      x00n=x002
      x10=x101+x102
      x01=x011+x012  
      x20=0.d0
      x11=0.d0
      x02=0.d0
      
      y00n=y002
      y10=y101+y102
      y01=y011+y012  
      y20=0.d0
      y11=0.d0
      y02=0.d0  
      endif

      if(igf.eq.22)then
      x00n=x002
      x10=x101+x102
      x01=x011+x012    
      x20=x201+x202
      x11=x111+x112
      x02=x021+x022
      
      y00n=y002
      y10=y101+y102
      y01=y011+y012    
      y20=y201+y202
      y11=y111+y112
      y02=y021+y022
      endif
      
      if(igf.eq.0)then    ! alles mitnehmen (corresponds to igf=3)
      x00n=x002+x003
      x10=x101+x102
      x01=x011+x012    
      x20=x201+x202
      x11=x111+x112
      x02=x021+x022
      
      y00n=y002+y003
      y10=y101+y102
      y01=y011+y012    
      y20=y201+y202
      y11=y111+y112
      y02=y021+y022
      endif

      pxf=px0
      pyf=py0

100   continue
      isteps=isteps+1

      g11=1.d0+x10+2.d0*x20*pxf+x11*pyf
      g12=x01+2.d0*x02*pyf+x11*pxf
      g21=y10+2.d0*y20*pxf+y11*pyf
      g22=1.d0+y01+2.d0*y02*pyf+y11*pxf

      if(iverbose.eq.1)then
      write(6,*)' g12 = ',g12
      write(6,*)' g21 = ',g21
      write(6,*)' inverse Determinante', 1.d0/(g11*g22-g12*g21) 
      endif
      
      hx=x00n+(1.d0+x10)*pxf+x01*pyf+
     &  x02*pyf**2+x20*pxf**2+x11*pxf*pyf-px0
      hy=y00n+(1.d0+y01)*pyf+y10*pxf+
     &  y20*pxf**2+y02*pyf**2+y11*pxf*pyf-py0
      
      dpxf=dreal(-(1.d0/(g11*g22-g12*g21))*(g22*hx-g12*hy))
      dpyf=dreal(-(1.d0/(g11*g22-g12*g21))*(-g21*hx+g11*hy))  
      
c      detie=1.d0-x101-y011-x102-y012+x101**2+y011**2+x101*y011+x011*y101
c      write(6,*) 'inverse Determinante entwickelt',detie       

c      write(6,*)' 1.d0+x10 ',1.d0+x10
c      write(6,*)' 2.d0*x20*pxf+x11*pyf ',2.d0*x20*pxf+x11*pyf
c      write(6,*)' g11, g12 ',g11,g12
c      write(6,*)' g21, g22 ',g21,g22         
c      write(6,*)' g22*hx-g12*hy ',g22*hx-g12*hy
c      write(6,*)' -g21*hx+g11*hy ',-g21*hx+g11*hy
      
      pxf=pxf+dpxf
      pyf=pyf+dpyf
      
c      write(6,*)' small_loc1 ',small_loc1      
      write(6,*)' dpxf, dpyf ',dpxf,dpyf
      write(6,*)'------------------------------------'    
      
      if(isteps.gt.imsteps)then
          write(6,*)'too many iterations. program terminated. '
          stop
      endif
      
      if((dabs(dpxf).gt.small_loc1).or.
     &   (dabs(dpyf).gt.small_loc1))goto 100
      
c---------- check
      check_x=-px0+pxf+x00n+x10*pxf+x01*pyf+x02*pyf**2+
     &  x20*pxf**2+x11*pxf*pyf
      
      check_y=-py0+pyf+y00n+y10*pxf+y01*pyf+y02*pyf**2+
     &  y20*pxf**2+y11*pxf*pyf
 
      if(iverbose.eq.1)then
c      write(6,*)'------------------------------------'      
c      write(6,*)' check_px,check_py new ',check_x,check_y
c      write(6,*)'------------------------------------'         
      endif
      
      return 
      end
     
c-------------------------------------------
      subroutine get_pxfpyf_old(Ax0,Ay0,xp0,yp0,px0,py0,
     &           pxf,pyf)      
c-------------------------------------------
c     Imaginärteile von pxf, pyf sind null wie es sein soll
c-------------------------------------------
      include 'STGFM.cmn'      
      complex*16 f11,f12,f21,f22,
     & x00n,x10,x01,x20,x11,x02,
     & y00n,y10,y01,y20,y11,y02
      complex*16 ac0,cc0,ac1,cc1,ac2,cc2,check_x,check_y      
      
      if(iverbose.eq.1)then
      write(6,*)'------------------------------------'    
      write(6,*)' x011 ',x011
      write(6,*)' x101 ',x101
      write(6,*)' x012 ',x012
      write(6,*)' x102 ',x102
      write(6,*)' x002 ',x002           
      write(6,*)'------------------------------------'      
      write(6,*)' y011 ',y011
      write(6,*)' y101 ',y101
      write(6,*)' y012 ',y012
      write(6,*)' y102 ',y102
      write(6,*)' y002 ',y002           
      write(6,*)'------------------------------------'      
      endif
      
      if(ilinear.eq.1)then
      x002=0.d0
      y002=0.d0
      endif
 
      if(ilinear.eq.10)then
      x002=0.d0
      y002=0.d0
      x011=0.d0
      y011=0.d0 
      x101=0.d0
      y101=0.d0          
      endif
      
      f11=1.d0-X101-x102+x011*y101+x101**2
      f12=-X011-x012+x011*(x101+y011)
      f21=-Y101-y102+y101*(x101+y011)
      f22=1.d0-Y011-y012+x011*y101+y011**2

      f11lin=1.d0-X101
      f12lin=-X011
      f21lin=-Y101
      f22lin=1.d0-Y011 
      
c      ac=f11*xp0+f21*yp0-x002+ax0*(1.d0-x101)-x011*ay0
c      cc=f21*xp0+f22*yp0-y002+ay0*(1.d0-y011)-y101*ax0

c lin bedeutet totale Ordnung </= 2
c      aclin=f11lin*xp0+f21lin*yp0-x002+ax0*(1.d0-x101)-x011*ay0     
c      cclin=f21lin*xp0+f22lin*yp0-y002+ay0*(1.d0-y011)-y101*ax0     
      
c      pxf=dreal(aclin)
c      pyf=dreal(cclin)
      
c      pxflin=dreal(aclin)
c      pyflin=dreal(cclin)         

      ac0=xp0
      cc0=yp0
      
      ac1=-x101*xp0-x011*yp0+Ax0
      cc1=-y101*xp0-y011*yp0+Ay0
      
      ac2=(-x102+x101**2+x011*y101)*xp0+(-x012+x011*(x101+y011))*yp0-
     & x002-x011*ay0-x101*ax0
      cc2=(-y102+y101*(x101+y011))*xp0+(-y012+y011**2+x011*y101)*yp0-
     & y002-y101*ax0-y011*ay0 

      pxf0=dreal(ac0)
      pyf0=dreal(cc0)
      
      pxf1=dreal(ac1)
      pyf1=dreal(cc1)
           
      pxf2=dreal(ac2)
      pyf2=dreal(cc2)
      
      pxf=pxf0+pxf1+pxf2
      pyf=pyf0+pyf1+pyf2
     

c---------- check

      x00n=x002
      x10=x101+x102
      x01=x011+x012  
      x20=0.d0
      x11=0.d0
      x02=0.d0
      
      y00n=y002
      y10=y101+y102
      y01=y011+y012  
      y20=0.d0
      y11=0.d0
      y02=0.d0      

      check_x=-px0+pxf+x00n+x10*pxf+x01*pyf+x02*pyf**2+
     &  x20*pxf**2+x11*pxf*pyf
      
      check_y=-py0+pyf+y00n+y10*pxf+y01*pyf+y02*pyf**2+
     &  y20*pxf**2+y11*pxf*pyf
      
      write(6,*)'------------------------------------'      
      write(6,*)' check_px,check_py old ',check_x,check_y
      write(6,*)'------------------------------------'        
      
      return
      end
 
c-------------------------------------------
      subroutine get_fqpk    
c-------------------------------------------
c     Imaginärteile der fqpk sind alle null wie es sein soll
c-------------------------------------------
      include 'STGFM.cmn'   
    
      f101=eta101(0)*(z0-zf)
      f011=eta011(0)*(z0-zf)      
      f111=eta111(0)*(z0-zf)
      f102=eta102(0)*(z0-zf)
      f012=eta012(0)*(z0-zf)
      f201=eta201(0)*(z0-zf)      
      f021=eta021(0)*(z0-zf)    
      f202=eta202(0)*(z0-zf)    
      f112=eta112(0)*(z0-zf)         
      f022=eta022(0)*(z0-zf)      
c      write(6,*) 'start0', eta101(0)*(z0-zf)

      do nfou=-iord1,iord1   ! facex(0) = 0
        f101=f101+facex(nfou)*eta101(nfou)
        f011=f011+facex(nfou)*eta011(nfou)  
        f111=f111+facex(nfou)*eta111(nfou)  
        f201=f201+facex(nfou)*eta201(nfou)       
        f021=f021+facex(nfou)*eta021(nfou)
c        write(6,*) ' +++ '
c        write(6,*) 'nfou', nfou
c        write(6,*) 'facex(nfou)', facex(nfou)
c        write(6,*) 'eta101(nfou)', eta101(nfou)
c        write(6,*) 'f101', f101
      enddo
      
      do nfou=-iord2,iord2              
        f102=f102+facex(nfou)*eta102(nfou)       
        f012=f012+facex(nfou)*eta012(nfou)
        f202=f202+facex(nfou)*eta202(nfou)     
        f112=f112+facex(nfou)*eta112(nfou)
        f022=f022+facex(nfou)*eta022(nfou) 
      enddo

      if(iverbose.eq.1)then
      write(6,*)'----------------------------------'
      write(6,*)' f101, f102 = ',f101,f102
      write(6,*)' f011, f012 = ',f011,f012
      write(6,*)' f201, f202 = ',f201,f202
      write(6,*)' f111, f112 = ',f111,f112
      write(6,*)' f021, f022 = ',f021,f022
      write(6,*)'----------------------------------'     
      endif
      
ctesttesttest  
c      f102=0.d0      
c      f012=0.d0     ! der Term scheint falsch zu sein
      
      return
      end
      
c-------------------------------------------
      subroutine get_xfyf(x0,y0,pxf,pyf,xf,yf)      
c-------------------------------------------
      include 'STGFM.cmn'      

      complex*16 f10,f01,f20,f11,f02
      
      if(igf.eq.2)then    
      f10=f101
      f01=f011
      f20=0.d0
      f11=0.d0
      f02=0.d0
      endif
      
      if(igf.eq.11)then          
      f10=f101
      f01=f011
      f20=0.d0
      f11=0.d0
      f02=0.d0     
      endif   

      if(igf.eq.12)then          
      f10=f101+f102
      f01=f011+f012  
      f20=0.d0
      f11=0.d0
      f02=0.d0     
      endif
      
      if(igf.eq.21)then
      f10=f101
      f01=f011    
      f20=f201
      f11=f111
      f02=f021
      endif
      
      if(igf.eq.22)then
      f10=f101+f102
      f01=f011+f012    
      f20=f201+f202
      f11=f111+f112
      f02=f021+f022
      endif

      if(igf.eq.0)then
      f10=f101+f102
      f01=f011+f012    
      f20=f201+f202
      f11=f111+f112
      f02=f021+f022
      endif

c      write(6,*) 'x0, pxf, f101', x0, pxf, f101

      xf=x0+pxf*(zf-z0)+f10+f11*pyf+2.d0*f20*pxf
      yf=y0+pyf*(zf-z0)+f01+f11*pxf+2.d0*f02*pyf
      return
      end 

c-------------------------------------------
      subroutine get_xfyf_old(x0,y0,xp0,yp0,ax0,ay0,pxf,pyf,xf,yf)      
c-------------------------------------------
      include 'STGFM.cmn'      
      
      xf=x0+pxf*(zf-z0)+f101+f102+f111*(pyf0+pyf1)+f112*pyf0+
     & 2.d0*f201*(pxf0+pxf1)+2.d0*f202*pxf0
      yf=y0+pyf*(zf-z0)+f011+f012+f111*(pxf0+pxf1)+f112*pxf0+
     & 2.d0*f021*(pyf0+pyf1)+2.d0*f022*pyf0

c      xf=x0+pxf*(zf-z0)+f101+yp0*f111+2.d0*xp0*f201+
c     &   f102+(-y101*xp0-y011*yp0)*f111+
c     &   2.d0*(-x101*xp0-x011*yp0)*f201
c
c      yf=y0+pyf*(zf-z0)+f011+xp0*f111+2.d0*yp0*f021+
c     &   f012+(-x101*xp0-x011*yp0)*f111+
c     &   2.d0*(-y101*xp0-y011*yp0)*f021
      
      return
      end 
      
c-------------------------------------------
      subroutine get_xpfypf(pxf,pyf,Axf,Ayf,xpf,ypf)   
c-------------------------------------------
      include 'STGFM.cmn'      

      xpf=pxf-Axf
      ypf=pyf-Ayf
      
      return
      end                           

c-------------------------------------------
      subroutine cartesian2polar(x,y,r,phi)      
c-------------------------------------------
      include 'STGFM.cmn'      

      r=dsqrt(x*x+y*y)
      phi=datan2(y,x)
 
      return      
      end      
       
c-------------------------------------------
      subroutine polar2cartesian(r,phi,x,y)      
c-------------------------------------------      
      include 'STGFM.cmn'      
      
      x=r*dcos(phi)
      y=r*dsin(phi)
  
      return      
      end            

c----------------------------------------------------------
      integer function clen1(str)
c----------------------------------------------------------
      character*(*) str
      do 23000 i=1,80
      if(.not.(ichar(str(i:i)) .eq. 32))goto 23002
      clen1 = i-1
      goto 23001
23002 continue
      clen1 = i
23003 continue
23000 continue
23001 continue
      return
      end

c----------------------------------------------------------
      integer function ifacult(n)
c----------------------------------------------------------
      if(n.eq.0)then
      ifacult=1
      else
      ifacult=1
      do i=1,n
      ifacult=ifacult*i
      enddo
      endif
      return
      end


