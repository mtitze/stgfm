*CMZ :  3.02/04 13/03/2015  10.38.25  by  Michael Scheer
*CMZ :  2.69/02 02/11/2012  16.40.18  by  Michael Scheer
*CMZ :  2.68/05 04/09/2012  13.30.58  by  Michael Scheer
*CMZ :  2.68/04 03/09/2012  11.52.24  by  Michael Scheer
*CMZ :  2.68/03 31/08/2012  09.45.28  by  Michael Scheer
*-- Author : Michael Scheer
      subroutine urad(
     &  gammai,dgamtot,
     &  xelec,yelec,zelec,vxelec,vyelec,vzelec,
     &  xf,yf,zf,efx,efy,efz,
     &  xexit,yexit,zexit,vnxex,vnyex,vnzex,texit,ds,
     &  nstep,ndim,traxyz,
     &  xobsv,yobsv,zobsv,phelow,phehig,
     &  nphener,phener,aradx,arady,aradz,stokes,powden,
     &  ieneloss,ivelofield
     &  ,istatus)

c Author: Michael Scheer, Michael.Scheer@Helmholtz-Berlin.de

c NO WARRANTY

*KEEP,gplhint.
!******************************************************************************
!
!      Copyright 2013 Helmholtz-Zentrum Berlin (HZB)
!      Hahn-Meitner-Platz 1
!      D-14109 Berlin
!      Germany
!
!      Author Michael Scheer, Michael.Scheer@Helmholtz-Berlin.de
!
! -----------------------------------------------------------------------
!
!    This program is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy (wave_gpl.txt) of the GNU General Public
!    License along with this program.
!    If not, see <http://www.gnu.org/licenses/>.
!
!    Dieses Programm ist Freie Software: Sie koennen es unter den Bedingungen
!    der GNU General Public License, wie von der Free Software Foundation,
!    Version 3 der Lizenz oder (nach Ihrer Option) jeder spaeteren
!    veroeffentlichten Version, weiterverbreiten und/oder modifizieren.
!
!    Dieses Programm wird in der Hoffnung, dass es nuetzlich sein wird, aber
!    OHNE JEDE GEWAEHRLEISTUNG, bereitgestellt; sogar ohne die implizite
!    Gewaehrleistung der MARKTFAEHIGKEIT oder EIGNUNG FueR EINEN BESTIMMTEN ZWECK.
!    Siehe die GNU General Public License fuer weitere Details.
!
!    Sie sollten eine Kopie (wave_gpl.txt) der GNU General Public License
!    zusammen mit diesem Programm erhalten haben. Wenn nicht,
!    siehe <http://www.gnu.org/licenses/>.
!
!******************************************************************************
*KEND.

c This subroutine calculates the trajectory and the synchrotron radiation
c of an electron passing a magnetic and electric field. The electric field
c part is very preliminary and not yet fully tested. The FORTRAN random
c generator is used in uradrndm. It should be replaced by a better one
c provided by the user.

c The fields B=(bx,by,bz) and E=(ex,ey,ez) are calculated in the routine
c uradfield(x,y,z,bx,by,bz,ex,ey,ez,istatus) provided by the user.
c As an example uradbmap.f may be used, which reads a 3D map of the mag. field

c Coordinate system (right handed):
c -----------------------------------

c      x: longitudinal direction
c      y: transversal vertical direction
c      z: transversal horizontal direction

c Units:
c ------

c SI-units: m, Tesla, sec, V etc., but eV for the photon energy
c The flux-density is given in Nph/s/mm**2/0.1%BW, the power-density in W/mm**2

c Input:
c --------

c real*8 gammai: Gamma factor of the e-

c real*8 xelec:  Initial x of e-
c real*8 yelec:  Initial y of e-
c real*8 zelec:  Initial z of e-

c real*8 vxelec:  Initial velocity in x of e-
c real*8 vyelec:  Initial velocity in y of e-
c real*8 vzelec:  Initial velocity in z of e-
c The velocity is internally normalized, so the input norm does not matter

c real*8 xf: x of point in exit plane
c real*8 yf: y of point in exit plane
c real*8 zf: z of point in exit plane

c real* vnxex: x component of normal vector of exit plane
c real* vnyex: y component of normal vector of exit plane
c real* vnzex: z component of normal vector of exit plane

c The tracking stops, when the electron hits the exit plane. The size of the last
c step is corrected such that the plane is hit.

c integer nstep: If nstep > 0, the trajectory array traxyz is filled, see below
c integer ndim: Dimension of traxyz, see below

c real*8 phelow: Lowest photon energy / eV
c real*8 phehig: Higest photon energy / eV

c integer nphener: Number of equidistant photon energies

c integer ieneloss:  0: no energy loss due to synchotron radiation
c                    1: continous energy loss due to synchotron radiation
c integer ieneloss: -1: no energy loss due to synchotron radiation with quantum
c                       fluctuations

c integer ivelofield: Contral flag for the calculation of the velocity field:
c                    0: the spectrum includes the velocity field
c                    1: the specrum does not include the velocity field
c                    2: the spectrum includes only the velocity field

c Output:
c -------

c integer: istatus: Status flag:
c  0: no error found
c -1: initial gamma or velocity zero
c -2: dimension ndim of traxyz exceeded
c -3: bad value of ivelofield
c  else: status of uradfield

c real*8 xexit: x of last point of the trajectory
c real*8 yexit: y of last point of the trajectory
c real*8 zexit: z of last point of the trajectory
c real*8 texit: t of last point of the trajectory

c real*8 vnxex: x component of norm. velocity vector of last point
c real*8 vnyex: y component of norm. velocity vector of last point
c real*8 vnzex: z component of norm. velocity vector of last point

c real*8 phener(nphener): Array of equidistant photon energies

c integer nstep: Number of tracking steps in array traxyz:
c               At input time, nstep is the nth step to be stored in traxyz.
c               At output time, it is the number of stored items.

c real*8 traxyz(1:14,i): Array:
c        traxyz(1,i):  x
c        traxyz(2,i):  y
c        traxyz(3,i):  z
c        traxyz(4,i):  tracking time
c        traxyz(5,i):  x-comp. of norm. velocity vector
c        traxyz(6,i):  y-comp. of norm. velocity vector
c        traxyz(7,i):  z-comp. of norm. velocity vector
c        traxyz(8,i):  x-comp. of mag. field in the center of the step
c        traxyz(9,i):  y-comp. of mag. field in the center of the step
c        traxyz(10,i): z-comp. of mag. field in the center of the step
c        traxyz(11,i): gamma
c        traxyz(12,i): x-comp. of elec. field in the center of the step
c        traxyz(13,i): y-comp. of elec. field in the center of the step
c        traxyz(14,i): z-comp. of elec. field in the center of the step

c The phase is calculated by phase=phase0+n*dt*dphase,
c where dphase is the phase difference of the nth step dt. Phase0=
c (xobsv-xelec)/clight. The phase factor of the integrand is
c exp(i*omega*phase),where omega referes to the considered photon energy.
c complex*16 aradx(i): x-comp. of amplitude of radiation field of phener(i)
c complex*16 arady(i): y-comp. of amplitude of radiation field of phener(i)
c complex*16 aradz(i): z-comp. of amplitude of radiation field of phener(i)

c real*8 array of Stokes parameters of ith photon energy:
c        stokes(1,i): S0, i.e. total intensity
c        stokes(2,i): S1, Stokes parameter of linear +/- 90 degree polarisation
c        stokes(3,i): S1, Stokes parameter of linear +/- 45 degree polarisation
c        stokes(4,i): S1, Stokes parameter of circular polarisation
c        S0 = sqrt(S1**2+S2**2+S3**2)

c Compilation:
c ------------

c For uradbmap at least F90 is required.
c The line length exceeds 72 characters, please use an appropriate
c compiler option. It is recommended to use compiler options to initialize all
c variables to zero and to treat them as ,,saved''

      implicit none

      complex*16 aradx(nphener),arady(nphener),aradz(nphener),
     &  ziom,zi,zidom,zone,ziomr1,zicr1,zic,
     &  expom1,expom,dexpomph1,dexpomph,ddexpomph,dexpom,
     &  expomv2,vstokes(4,3),
     &  apolh,apolr,apoll,apol45,dum3

      double precision
     &  gammai,dgamtot,dt2,powden,t,phase,
     &  xelec,yelec,zelec,vxelec,vyelec,vzelec,
     &  xexit,yexit,zexit,vnxex,vnyex,vnzex,texit,
     &  xobsv,yobsv,zobsv,phelow,phehig,
     &  phener(nphener),dom2,c,rspn
     &  ,traxyz(14,ndim),stokes(4,nphener),x1,y1,z1,vx1,vy1,
     &  vz1,x2,y2,z2,vx2,vy2,vz2
     &  ,ds,dtim,bshift,x2b,y2b,z2b,bx1,by1,bz1,bx2,by2,bz2
     &  ,dgamsum,gamma,dt
     &  ,x2int,y2int,z2int,ddt,dddt,ddt2,dddt2
     &  ,vx2int,vy2int,vz2int,vxpint,vypint,vzpint
     &  ,vxp,vyp,vzp
     &  ,x3int,y3int,z3int,vx3int,vy3int,vz3int,ddddt,ddddt2
     &  ,efx,efy,efz,xf,yf,zf,dist1,dist2,disti
     &  ,dtim0,beta,vn,efx2,efy2,efz2,t1,t2,clight,c1,
     &  dgamma,vxsign,bx,by,bz,bpx,bpy,bpz,rarg(5),px,py,pz,
     &  dphase,r,rx,ry,rz
     &  ,dom1,rnbx,rnby,rnbz,dum11,rnr4,br4,b3,rnr2,br2,bet1n,
     &  rnx,rny,rnz,r1,banwid,specnor,pownor,current,
     &  stok1,stok2,stok3,stok4,om,dom,hbarev,echarge,eps0,pi,vsto,dph,
     &  r0

      integer ieneloss,istatus,icharge,nphener,ivelofield,nthstep,
     &  nstep,ndim,kstep,lstep,ifreq,isto

      save

      data bshift/0.5d0/
      data clight/2.99792458d8/
      data hbarev/6.58211889D-16/
      data banwid/1.0d-3/
      data current/0.10d0/
      data eps0/8.854187817D-12/
      data echarge/1.602176462D-19/
      data pi/3.14159265358979d0/

      data zi/(0.0d0,1.0d0)/
      data zone/(1.0d0,0.0d0)/

      dph=0.0d0
      phener(1)=phelow
      if (nphener.gt.1) dph=(phehig-phelow)/(nphener-1)
      do ifreq=2,nphener
        phener(ifreq)=phener(ifreq-1)+dph
      enddo

      istatus=0
      icharge=-1

      x1=xelec
      y1=yelec
      z1=zelec
      vx1=vxelec
      vy1=vyelec
      vz1=vzelec
      t1=0.0d0

      gamma=gammai
      beta=dsqrt((1.d0-1.d0/gamma)*(1.d0+1.d0/gamma))
      vn=sqrt(vx1*vx1+vy1*vy1+vz1*vz1)
      vx1=vx1/vn*clight*beta
      vy1=vy1/vn*clight*beta
      vz1=vz1/vn*clight*beta

c vxsign takes care for the direction of flight, since particle must gain
c energy if tracked back

      if (vx1.lt.0) then
        vxsign=-1.0d0
      else
        vxsign=1.0d0
      endif

      dgamsum=0.0d0
      dgamtot=0.0d0
      powden=0.0d0
      aradx=(0.0d0,0.0d0)
      arady=(0.0d0,0.0d0)
      aradz=(0.0d0,0.0d0)

      dtim=ds/clight
      dt=dtim
      dt2=dtim*bshift
      dtim0=dtim

      x2=x1
      y2=y1
      z2=z1
      t2=t1

      vx2=vx1
      vy2=vy1
      vz2=vz1

      x2b=x1+vx1*dt2
      y2b=y1+vy1*dt2
      z2b=z1+vz1*dt2

      call uradfield(x2b,y2b,z2b,bx2,by2,bz2,efx2,efy2,efz2,istatus)
      if (istatus.ne.0) goto 9000

      nthstep=nstep
      nstep=0
      vn=sqrt(vx2*vx2+vy2*vy2+vz2*vz2)

      if (gamma.le.0.0d0.or.vn.le.0.0d0) then
        istatus=-1
        return
      endif

      kstep=-1
      nstep=0

      if(nthstep.gt.0) then
        nstep=nstep+1
        if (nstep.gt.ndim) then
          istatus=-2
          goto 9000
        endif
        kstep=kstep+1
        if (kstep.eq.nthstep) kstep=0
        traxyz(1,nstep)=x2
        traxyz(2,nstep)=y2
        traxyz(3,nstep)=z2
        traxyz(4,nstep)=t2
        traxyz(5,nstep)=vx2/vn
        traxyz(6,nstep)=vy2/vn
        traxyz(7,nstep)=vz2/vn
        traxyz(8,nstep)=bx2
        traxyz(9,nstep)=by2
        traxyz(10,nstep)=bz2
        traxyz(11,nstep)=gamma
        traxyz(12,nstep)=efx2
        traxyz(13,nstep)=efy2
        traxyz(14,nstep)=efz2
      endif

      dom=0.0d0
      om=0.0d0
      if (nphener.gt.1) then
        om=phener(1)/hbarev
        dom=(phener(2)-phener(1))/hbarev
      else if (nphener.eq.1) then
        om=phener(1)/hbarev
      endif

      zidom=zi*dom
      ziom=zi*om
      zic=zi*clight

      lstep=0
      t=-dt
      r0=xobsv-xelec
      r=sqrt((xobsv-x1)**2+((yobsv-y1)**2+(zobsv-z1)**2))
      PHASE=(r-r0)*c1
      expom1=zone
      dexpomph1=zone
      c=clight
      c1=1.0d0/clight

c--- Loop der Trajektorie

1000  continue

      x1=x2
      y1=y2
      z1=z2

      t1=t2

      vx1=vx2
      vy1=vy2
      vz1=vz2

      bx1=bx2
      by1=by2
      bz1=bz2

      x2b=x1+vx1*dt2
      y2b=y1+vy1*dt2
      z2b=z1+vz1*dt2

      call uradfield(x2b,y2b,z2b,bx2,by2,bz2,efx2,efy2,efz2,istatus)
      if (istatus.ne.0) goto 9000

      call uradstep(x1,y1,z1,vx1,vy1,vz1,bx2,by2,bz2,efx2,efy2,efz2,
     &  dtim,
     &  x2,y2,z2,vx2,vy2,vz2,vxp,vyp,vzp,gamma,icharge,ieneloss,dgamma)

      if (ieneloss.ne.0) then
        dgamsum=dgamsum+vxsign*dgamma
        if (abs(dgamsum).gt.gamma*1.0d-8) then
          gamma=gamma+dgamsum
          dgamtot=dgamtot+dgamsum
          dgamsum=0.0d0
        endif
        beta=dsqrt((1.d0-1.d0/gamma)*(1.d0+1.d0/gamma))
        vn=sqrt(vx2*vx2+vy2*vy2+vz2*vz2)
        vx2=vx2/vn*clight*beta
        vy2=vy2/vn*clight*beta
        vz2=vz2/vn*clight*beta
      endif

      BX=VX2*C1
      BY=VY2*C1
      BZ=VZ2*C1

      BPX=VXP*C1
      BPY=VYP*C1
      BPZ=VZP*C1

C CONTRIBUTION OF TIME STEP TO SYNCHROTRON RADIATION {

C REAL PART OF INTEGRAND {

      RX=XOBSV-X2
      RY=YOBSV-Y2
      RZ=ZOBSV-Z2

      R=SQRT(RX*RX+RY*RY+RZ*RZ)
      R1=1.D0/R
      ZICR1=ZIC*R1

      RNX=RX*R1
      RNY=RY*R1
      RNZ=RZ*R1

C--- THE DISTANCE R IS INTRODUCED HERE EXPLICITLY (S. PROGRAM OF CHAOEN WANG

      BET1N=(1.D0-BX*RNX)-BY*RNY-BZ*RNZ

      br2=by**2+bz**2
      rnr2=rny**2+rnz**2
      b3=beta**3
      br4=br2**2
      rnr4=rnr2**2

      if(br2.lt.1.0d-4.and.rnr2.lt.1.0d-4) then
        bet1n=
     &    1.0d0/(1+beta)/gamma**2
     &    +beta*(rnr2/2.0d0
     &    +rnr4/8.0d0)
     &    +(br2/2.0d0
     &    -br2*rnr2/4.0d0
     &    -br2*rnr4/16.0d0)/beta
     &    +b3*br4*(1.0d0/8.0d0
     &    -rnr2/16.0d0
     &    -rnr4/64.0d0)
     &    -by*rny
     &    -bz*rnz
      endif

      DUM11=1.D0/BET1N
      DOM1=1.D0/(R*BET1N*BET1N)

      RNBX=RNX-BX
      RNBY=RNY-BY
      RNBZ=RNZ-BZ

      PX=(RNBY*BPZ-RNBZ*BPY)
      PY=(RNBZ*BPX-RNBX*BPZ)
      PZ=(RNBX*BPY-RNBY*BPX)

      IF (IVELOFIELD.EQ.0.OR.IVELOFIELD.EQ.2) THEN !2 WEGEN POWER
        DOM2=C*DOM1*R1/GAMMA**2
        RARG(1)=(RNY*PZ-RNZ*PY)*DOM1+(RNX-BX)*DOM2
        RARG(2)=(RNZ*PX-RNX*PZ)*DOM1+(RNY-BY)*DOM2
        RARG(3)=(RNX*PY-RNY*PX)*DOM1+(RNZ-BZ)*DOM2
      ELSE IF (IVELOFIELD.EQ.1) THEN
        RARG(1)=(RNY*PZ-RNZ*PY)*DOM1
        RARG(2)=(RNZ*PX-RNX*PZ)*DOM1
        RARG(3)=(RNX*PY-RNY*PX)*DOM1
      ELSE IF (IVELOFIELD.LT.0) THEN
        DOM2=C*DOM1*R1/GAMMA**2
        RARG(1)=(RNX-BX)*DOM2
        RARG(2)=(RNY-BY)*DOM2
        RARG(3)=(RNZ-BZ)*DOM2
      ELSE	!IVELOFIELD
        istatus=-3
        return
      ENDIF	!IVELOFIELD

C DO NOT USE, RESULTS IN NUMERICAL PROBLEMS	    RARG(4)=T+R*C1

      DPHASE=BET1N*DT
      PHASE=PHASE+DPHASE

      RARG(5)=(RARG(1)*RARG(1)+RARG(2)*RARG(2)+RARG(3)*RARG(3))*DUM11

C REAL PART OF INTEGRAND }

C COMPLEX PART OF INTEGRAND {

C    ASSUMES phener(I+1)=2*phener(I)   FOR IFREQ2P=2
C    OR phener(I+1)=phener(I)+DELTA    FOR IFREQ2P>2

C--- LOOP OVER ALL FREQUENCES

      IFREQ=1

      OM=phener(IFREQ)/hbarev
      ZIOM=ZI*OM

      EXPOM=EXPOM1
      DEXPOMPH1=EXP(ZIOM*DPHASE)
      DEXPOMPH=DEXPOMPH1

      IF(nphener.GT.1) THEN
        DEXPOM=EXP(ZIDOM*PHASE)
        DDEXPOMPH=EXP(ZIDOM*DPHASE)
      ENDIF	!IFREQ2P

      powden=powden+rarg(5)*dt

      IF (IVELOFIELD.NE.2) THEN

        dum3=expom*(zone-dexpomph)/om/bet1n
        aradx(ifreq)=aradx(ifreq)+rarg(1)*dum3
        arady(ifreq)=arady(ifreq)+rarg(2)*dum3
        aradz(ifreq)=aradz(ifreq)+rarg(3)*dum3

      ELSE !IVELOFIELD

        EXPOMV2=R1/BET1N*EXPOM*(ZONE-DEXPOMPH)
        ZIOMR1=ZONE+ZICR1/OM

        aradx(ifreq)=aradx(ifreq)+EXPOMV2*(BX-RNX*ZIOMR1)
        arady(ifreq)=arady(ifreq)+EXPOMV2*(BX-RNX*ZIOMR1)
        aradz(ifreq)=aradz(ifreq)+EXPOMV2*(BX-RNX*ZIOMR1)

      ENDIF !IVELOFIELD

      IF (IVELOFIELD.NE.2) THEN

        DO IFREQ=2,nphener

          OM=OM+DOM
          EXPOM=EXPOM*DEXPOM
          DEXPOMPH=DEXPOMPH*DDEXPOMPH

          EXPOMV2=1.0D0/BET1N/OM*EXPOM*(ZONE-DEXPOMPH)
          aradx(ifreq)=aradx(ifreq)+RARG(1)*EXPOMV2
          arady(ifreq)=arady(ifreq)+RARG(2)*EXPOMV2
          aradz(ifreq)=aradz(ifreq)+RARG(3)*EXPOMV2

        ENDDO   !LOOP OVER ALL FREQUENCES

      else

        DO IFREQ=2,nphener

          OM=OM+DOM
          EXPOM=EXPOM*DEXPOM
          DEXPOMPH=DEXPOMPH*DDEXPOMPH

          EXPOMV2=R1/BET1N*EXPOM*(ZONE-DEXPOMPH)
          ZIOMR1=ZONE+ZICR1/OM

          aradx(ifreq)=aradx(ifreq)+EXPOMV2*(BX-RNX*ZIOMR1)
          arady(ifreq)=aradx(ifreq)+EXPOMV2*(BX-RNX*ZIOMR1)
          aradz(ifreq)=aradx(ifreq)+EXPOMV2*(BX-RNX*ZIOMR1)

        ENDDO   !LOOP OVER ALL FREQUENCES

      ENDIF !IVELOFIELD


C COMPLEX PART OF INTEGRAND }

C CONTRIBUTION OF TIME STEP TO SYNCHROTRON RADIATION }

      EXPOM1=EXPOM1*DEXPOMPH1

      t2=t1+dtim

c ef is normal vector of perpendiculare plane at the end of the reference orbi
c dist is distance of electron to this plane
c tracking stops if trajectorie hits this plane

      dist2=(x2-xf)*efx+(y2-yf)*efy+(z2-zf)*efz

      if ( dist2.lt.0.0d0)  then
        kstep=kstep+1
        if (kstep.eq.nthstep) kstep=0
        if(kstep.eq.0) then
          nstep=nstep+1
          if (nstep.gt.ndim) then
            istatus=-2
            goto 9000
          endif
          traxyz(1,nstep)=x2
          traxyz(2,nstep)=y2
          traxyz(3,nstep)=z2
          traxyz(4,nstep)=t2
          traxyz(5,nstep)=vx2/vn
          traxyz(6,nstep)=vy2/vn
          traxyz(7,nstep)=vz2/vn
          traxyz(8,nstep)=bx2
          traxyz(9,nstep)=by2
          traxyz(10,nstep)=bz2
          traxyz(11,nstep)=gamma
          traxyz(12,nstep)=efx2
          traxyz(13,nstep)=efy2
          traxyz(14,nstep)=efz2
        endif
        goto 1000
      endif

c--- ende of trajectory, dist2 not exactly zero, correct x2

      dist1=(x1-xf)*efx+(y1-yf)*efy+(z1-zf)*efz

      ddt=dtim*dabs(dist1)/(dabs(dist1)+dabs(dist2))
      ddt2=ddt*bshift

      x2b=x1+vx1*ddt2
      y2b=y1+vy1*ddt2
      z2b=z1+vz1*ddt2

      call uradfield(x2b,y2b,z2b,bx2,by2,bz2,efx2,efy2,efz2,istatus)
      if (istatus.ne.0) goto 9000

      call uradstep(x1,y1,z1,vx1,vy1,vz1,bx2,by2,bz2,efx2,efy2,efz2,ddt,
     &  x2int,y2int,z2int,vx2int,vy2int,vz2int,
     &  vxpint,vypint,vzpint,gamma,icharge,
     &  ieneloss,dgamma)

      disti=(x2int-xf)*efx+(y2int-yf)*efy+(z2int-zf)*efz
      dddt=ddt

      if (dist1.ne.0.) dddt=ddt-ddt*disti/dabs(dist1)

      dddt2=dddt*bshift

      x2b=x1+vx1*dddt2
      y2b=y1+vy1*dddt2
      z2b=z1+vz1*dddt2

      call uradfield(x2b,y2b,z2b,bx2,by2,bz2,efx2,efy2,efz2,istatus)
      if (istatus.ne.0) goto 9000

      call uradstep(x1,y1,z1,vx1,vy1,vz1,bx2,by2,bz2,efx2,efy2,efz2,
     &  dddt,
     &  x3int,y3int,z3int,vx3int,vy3int,vz3int,
     &  vxpint,vypint,vzpint,gamma,icharge,ieneloss,dgamma)

      disti=(x3int-xf)*efx+(y3int-yf)*efy+(z3int-zf)*efz
      ddddt=dddt
      if (dist1.ne.0.) ddddt=dddt-dddt*disti/dabs(dist1)
      ddddt2=bshift*ddddt

      x2b=x1+vx1*ddddt2
      y2b=y1+vy1*ddddt2
      z2b=z1+vz1*ddddt2

      call uradfield(x2b,y2b,z2b,bx2,by2,bz2,efx2,efy2,efz2,istatus)
      if (istatus.ne.0) goto 9000

      call uradstep(x1,y1,z1,vx1,vy1,vz1,bx2,by2,bz2,efx2,efy2,efz2,
     &  ddddt,
     &  x2,y2,z2,vx2,vy2,vz2,
     &  vxp,vyp,vzp,gamma,icharge,
     &  ieneloss,dgamma)

      t2=t1+ddddt

c ieneloss already applied for last step, ignore error due to
c corrected step-length here

      if(nthstep.ge.0) then
        nstep=nstep+1
        if (nstep.gt.ndim) then
          istatus=-2
          goto 9000
        endif
        nthstep=nstep
        traxyz(1,nstep)=x2
        traxyz(2,nstep)=y2
        traxyz(3,nstep)=z2
        traxyz(4,nstep)=t2
        traxyz(5,nstep)=vx2/vn
        traxyz(6,nstep)=vy2/vn
        traxyz(7,nstep)=vz2/vn
        traxyz(8,nstep)=bx2
        traxyz(9,nstep)=by2
        traxyz(10,nstep)=bz2
        traxyz(11,nstep)=gamma
        traxyz(12,nstep)=efx2
        traxyz(13,nstep)=efy2
        traxyz(14,nstep)=efz2
      endif

      t2=t1+ddddt

      xexit=x2
      yexit=y2
      zexit=z2

      vn=sqrt(vx2*vx2+vy2*vy2+vz2*vz2)
      vnxex=vx2/vn
      vnyex=vy2/vn
      vnzex=vz2/vn

      texit=t2

9000  continue

      specnor=
     &  banwid
     & /(4.0d0*pi**2*clight*hbarev)
     & /(4.0d0*pi*eps0)
     &  *current/1.0d6 !per mm**2

      pownor=echarge/16.0d0/pi/pi/eps0/clight*current/1.0d6 !W/mm**2

      rspn=sqrt(specnor)

      vstokes(1,1)=( 0.0d0,        0.0d0)      !horizontal polarization
      vstokes(1,2)=( 0.0d0,        0.0d0)
      vstokes(1,3)=(-1.0d0,       -1.0d0)

      vstokes(2,1)=( 0.0d0,        0.0d0)      !right handed polarization
      vstokes(2,2)=( 0.0d0,       -1.0d0)
      vstokes(2,3)=(+1.0d0,        0.0d0)

      vstokes(3,1)=( 0.0d0,        0.0d0)      !left handed polarization
      vstokes(3,2)=( 0.0d0,       -1.0d0)
      vstokes(3,3)=(-1.0d0,        0.0d0)

      vstokes(4,1)=( 0.0d0,        0.0d0)      !45 degree linear polarization
      vstokes(4,2)=( 1.0d0,        0.0d0)
      vstokes(4,3)=(-1.0d0,        0.0d0)

      do isto=1,4
        vsto=dsqrt
     &    (cdabs(vstokes(isto,1))**2
     &    +cdabs(vstokes(isto,2))**2
     &    +cdabs(vstokes(isto,3))**2)
        vstokes(isto,1)=vstokes(isto,1)/vsto
        vstokes(isto,2)=vstokes(isto,2)/vsto
        vstokes(isto,3)=vstokes(isto,3)/vsto

      enddo

      do ifreq=1,nphener

        aradx(ifreq)=aradx(ifreq)*rspn
        arady(ifreq)=arady(ifreq)*rspn
        aradz(ifreq)=aradz(ifreq)*rspn

        apolh=
     &    aradx(ifreq)*conjg(vstokes(1,1))
     &    +arady(ifreq)*conjg(vstokes(1,2))
     &    +aradz(ifreq)*conjg(vstokes(1,3))

        apolr=
     &    aradx(ifreq)*conjg(vstokes(2,1))
     &    +arady(ifreq)*conjg(vstokes(2,2))
     &    +aradz(ifreq)*conjg(vstokes(2,3))

        apoll=
     &    aradx(ifreq)*conjg(vstokes(3,1))
     &    +arady(ifreq)*conjg(vstokes(3,2))
     &    +aradz(ifreq)*conjg(vstokes(3,3))

        apol45=
     &    aradx(ifreq)*conjg(vstokes(4,1))
     &    +arady(ifreq)*conjg(vstokes(4,2))
     &    +aradz(ifreq)*conjg(vstokes(4,3))

        stok1=
     &    apolr*conjg(apolr)+
     &    apoll*conjg(apoll)

        stok2=-stok1+
     &    2.*apolh*conjg(apolh)

        stok3=
     &    2.*apol45*conjg(apol45)-
     &    stok1

        stok4=
     &    apolr*conjg(apolr)-
     &    apoll*conjg(apoll)

        stokes(1,ifreq)=stok1
        stokes(2,ifreq)=stok2
        stokes(3,ifreq)=stok3
        stokes(4,ifreq)=stok4

      enddo !nphener

      powden=powden*pownor

      return
      end

c      include 'uradfield.f'
c      include 'uradbmap.f'
      include 'uradrndm.f'
*CMZ :  3.02/04 11/03/2015  13.19.06  by  Michael Scheer
*CMZ :  3.01/04 20/05/2014  12.49.15  by  Michael Scheer
*CMZ :  2.70/00 12/11/2012  11.53.09  by  Michael Scheer
*CMZ :  2.68/05 06/09/2012  15.53.22  by  Michael Scheer
*CMZ :  2.68/04 03/09/2012  09.07.11  by  Michael Scheer
*CMZ :  2.68/03 27/08/2012  09.35.31  by  Michael Scheer
*-- Author :    Michael Scheer
c*******************************************************************************
      subroutine uradstep(x1,y1,z1,vx1,vy1,vz1,bxin,byin,bzin,
     &  efieldx,efieldy,efieldz,
     &  dtim,
     &  x2,y2,z2,vx2,vy2,vz2,vxp,vyp,vzp,gamma,icharge,
     &  ieneloss,dgamma)
c*******************************************************************************

c Author: Michael Scheer, Michael.Scheer@Helmholtz-Berlin.de

c NO WARRANTY

*KEEP,gplhint.
!******************************************************************************
!
!      Copyright 2013 Helmholtz-Zentrum Berlin (HZB)
!      Hahn-Meitner-Platz 1
!      D-14109 Berlin
!      Germany
!
!      Author Michael Scheer, Michael.Scheer@Helmholtz-Berlin.de
!
! -----------------------------------------------------------------------
!
!    This program is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy (wave_gpl.txt) of the GNU General Public
!    License along with this program.
!    If not, see <http://www.gnu.org/licenses/>.
!
!    Dieses Programm ist Freie Software: Sie koennen es unter den Bedingungen
!    der GNU General Public License, wie von der Free Software Foundation,
!    Version 3 der Lizenz oder (nach Ihrer Option) jeder spaeteren
!    veroeffentlichten Version, weiterverbreiten und/oder modifizieren.
!
!    Dieses Programm wird in der Hoffnung, dass es nuetzlich sein wird, aber
!    OHNE JEDE GEWAEHRLEISTUNG, bereitgestellt; sogar ohne die implizite
!    Gewaehrleistung der MARKTFAEHIGKEIT oder EIGNUNG FueR EINEN BESTIMMTEN ZWECK.
!    Siehe die GNU General Public License fuer weitere Details.
!
!    Sie sollten eine Kopie (wave_gpl.txt) der GNU General Public License
!    zusammen mit diesem Programm erhalten haben. Wenn nicht,
!    siehe <http://www.gnu.org/licenses/>.
!
!******************************************************************************
*KEND.

      implicit none

      double precision bbet1,vsenk1,vsenk2,tz,tz2,tz3,tz4,tz5,
     &  vsenkzyk,
     &  efieldx,efieldy,efieldz

      double precision bx2,by2,bz2,bsq,bbet,bux,buy,buz,v0sq,
     &  vx1,vy1,vz1,
     &  v0bet,vpar,vparx,vpary,vparz,vparsq,vsenk,dtim,bxin,byin,bzin
      double precision x1n,y1n,z1n,x2n,y2n,z2n

      double precision x2,y2,z2,vx2,vy2,vz2,x1,y1,z1,vxp,vyp,vzp,zyk,
     &  gamma,sz,cz
     &  ,dgamma

      double precision bmovecut,dgammao,v0,v12n(3),emom,pel(3),
     &  dpphoton(3)

      double precision erad1,cgam1,pdum,echarge1,emasskg1,emassg1,
     &  clight1,pi1

      save

      integer icharge,ieneloss,ical

      data ical/0/

      data pi1/3.14159265358979d0/
      data emasskg1/9.10938188D-31/
      data emassg1/0.510998902D-3/
      data clight1/2.99792458d8/
      data echarge1/1.602176462D-19/

      data bmovecut/1.0d-7/

      if (ical.eq.0) then
        erad1=2.8179380D-15
        cgam1=4.0d0/3.0d0*pi1*erad1/emassg1**3
        pdum=cgam1/2.0d0/pi1*clight1*(clight1/1.0d9)**2*emassg1
        ical=1
      endif

      dgamma=0.0d0

      if (icharge.gt.0) then
        bx2=-bxin
        by2=-byin
        bz2=-bzin
      else
        bx2=bxin
        by2=byin
        bz2=bzin
      endif

      if(dabs(bx2).lt.bmovecut.and.dabs(by2).lt.bmovecut
     &    .and.dabs(bz2).lt.bmovecut) then

        vxp=0.d0
        vyp=0.d0
        vzp=0.d0

        x2=x1+vx1*dtim
        y2=y1+vy1*dtim
        z2=z1+vz1*dtim

        vx2=vx1
        vy2=vy1
        vz2=vz1

        goto 999

      endif	!b-cut

      bsq=bx2*bx2+by2*by2+bz2*bz2
      bbet=dsqrt(bsq)
      bbet1=1.0d0/bbet

      bux=bx2*bbet1
      buy=by2*bbet1
      buz=bz2*bbet1

c
c   Betrag von v0 paralell und senkrecht
c
      v0sq=vx1*vx1+vy1*vy1+vz1*vz1
      v0bet=dsqrt(v0sq)

c
c  vpar
c
      vpar=vx1*bux+vy1*buy+vz1*buz
      vparsq=vpar*vpar
      vparx=vpar*bux
      vpary=vpar*buy
      vparz=vpar*buz

c      vsenk2=(v0sq-vparsq)
      vsenk2=(v0bet-vpar)*(v0bet+vpar)

      if (vsenk2.le.0.0d0) then

c
c  Zeitableitung der Geschwindigkeit
c
        vxp=0.d0
        vyp=0.d0
        vzp=0.d0

c
c v(dtim),x(dtim) berechnen
c

        x2=x1+vx1*dtim
        y2=y1+vy1*dtim
        z2=z1+vz1*dtim

        vx2=vx1
        vy2=vy1
        vz2=vz1

        goto 999

      else    !(vsenk2.lt.0.0)

        vsenk=dsqrt(vsenk2)
        vsenk1=1.0d0/vsenk

c
c   vektor n1 berechnen
c

        x1n=(vx1-vpar*bux)*vsenk1
        y1n=(vy1-vpar*buy)*vsenk1
        z1n=(vz1-vpar*buz)*vsenk1

c
c  Vektor n2=(bux,buy,buz) kreuz n1
c

        x2n = buy*z1n - buz*y1n
        y2n = buz*x1n - bux*z1n
        z2n = bux*y1n - buy*x1n

c
c Zyklotronfrequenz
c
        zyk=(echarge1/(gamma*emasskg1))*bbet
c
c

        tz=zyk*dtim

        if (tz.le.0.03d0) then
          tz2=tz*tz
          tz3=tz2*tz
          tz4=tz3*tz
          tz5=tz4*tz
          cz=1.0d0-tz2/2.0d0+tz4/24.0d0
          sz=tz-tz3/6.0d0+tz5/120.0d0
        else
          cz=cos(tz)
          sz=sin(tz)
        endif

c
c  Zeitableitung der Geschwindigkeit
c

        vxp=vsenk*zyk*x2n
        vyp=vsenk*zyk*y2n
        vzp=vsenk*zyk*z2n

c
c v(dtim),x(dtim) berechnen
c

        vx2=vparx + vsenk*(x1n*cz+x2n*sz)
        vy2=vpary + vsenk*(y1n*cz+y2n*sz)
        vz2=vparz + vsenk*(z1n*cz+z2n*sz)

c
c x(dtim) berechnen
c

        vsenkzyk=vsenk/zyk

        x2=x1+vsenkzyk*(x2n+x1n*sz-x2n*cz)+vparx*dtim
        y2=y1+vsenkzyk*(y2n+y1n*sz-y2n*cz)+vpary*dtim
        z2=z1+vsenkzyk*(z2n+z1n*sz-z2n*cz)+vparz*dtim

        if (ieneloss.ne.0) then
          dgamma=-pdum*bsq*vsenk2/v0sq*gamma**2*dtim
          if (ieneloss.eq.-1) then
            v0=sqrt(v0sq)
            v12n(1)=vx2/v0
            v12n(2)=vy2/v0
            v12n(3)=vz2/v0
            call uradphoton(v12n,gamma,bx2,by2,bz2,
     &        dgamma,dtim,dpphoton) !dgamma will be overwritten
            if (dgamma.ne.0.0d0) then
              emom=emassg1*dsqrt((gamma-1)*(gamma+1)) !gev
              pel=emom*v12n+dpphoton
              dgamma=
     &          sqrt(1.0d0+(pel(1)**2+pel(2)**2+pel(3)**2)/emassg1**2)-
     &          gamma
c              emom=emassg1*dsqrt((gamma+dgamma-1)*(gamma+dgamma+1)) !Gev
              emom=sqrt(pel(1)**2+pel(2)**2+pel(3)**2)
              vx2=pel(1)/((gamma+dgamma)*emassg1)*clight1
              vy2=pel(2)/((gamma+dgamma)*emassg1)*clight1
              vz2=pel(3)/((gamma+dgamma)*emassg1)*clight1
            endif !dgamma
          endif !ieneloss .eq. -1
        endif !ieneloss

      endif   !(vsenk.lt.0.0)

999   continue

      dgammao=dgamma
      call uradestep(icharge,x2,y2,z2,vx2,vy2,vz2,
     &  efieldx,efieldy,efieldz,dtim,gamma,dgamma)
      dgamma=dgammao+dgamma

      vxp=(vx2-vx1)/dtim
      vyp=(vy2-vy1)/dtim
      vzp=(vz2-vz1)/dtim

      return
      end
*CMZ :  3.02/04 11/03/2015  13.19.06  by  Michael Scheer
*CMZ :  2.70/00 12/11/2012  11.53.09  by  Michael Scheer
*CMZ :  2.68/04 03/09/2012  13.24.02  by  Michael Scheer
*CMZ :  2.68/03 25/08/2012  14.52.34  by  Michael Scheer
*CMZ :  2.66/21 22/11/2011  13.51.05  by  Michael Scheer
*-- Author :    Michael Scheer
      subroutine uradestep(icharge,x,y,z,vx,vy,vz,
     &  efieldx,efieldy,efieldz,
     &  dt,gamma,dgamma)

c Author: Michael Scheer, Michael.Scheer@Helmholtz-Berlin.de

c NO WARRANTY

*KEEP,gplhint.
!******************************************************************************
!
!      Copyright 2013 Helmholtz-Zentrum Berlin (HZB)
!      Hahn-Meitner-Platz 1
!      D-14109 Berlin
!      Germany
!
!      Author Michael Scheer, Michael.Scheer@Helmholtz-Berlin.de
!
! -----------------------------------------------------------------------
!
!    This program is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy (wave_gpl.txt) of the GNU General Public
!    License along with this program.
!    If not, see <http://www.gnu.org/licenses/>.
!
!    Dieses Programm ist Freie Software: Sie koennen es unter den Bedingungen
!    der GNU General Public License, wie von der Free Software Foundation,
!    Version 3 der Lizenz oder (nach Ihrer Option) jeder spaeteren
!    veroeffentlichten Version, weiterverbreiten und/oder modifizieren.
!
!    Dieses Programm wird in der Hoffnung, dass es nuetzlich sein wird, aber
!    OHNE JEDE GEWAEHRLEISTUNG, bereitgestellt; sogar ohne die implizite
!    Gewaehrleistung der MARKTFAEHIGKEIT oder EIGNUNG FueR EINEN BESTIMMTEN ZWECK.
!    Siehe die GNU General Public License fuer weitere Details.
!
!    Sie sollten eine Kopie (wave_gpl.txt) der GNU General Public License
!    zusammen mit diesem Programm erhalten haben. Wenn nicht,
!    siehe <http://www.gnu.org/licenses/>.
!
!******************************************************************************
*KEND.

c Simple approach to electrical fields into account
c The energy and gamma of the particle is not changed, but returned in dgamma

      double precision x,y,z,vx,vy,vz,dt,gamma,dgamma
      double precision px,py,pz,qdt,eg,vn,pn,
     &  px0,py0,pz0,pn0,pn2,pn02,de,dpx,dpy,dpz,dtm,
     &  efieldx,efieldy,efieldz,
     &  echarge1,emasskg1,emassg1

      integer icharge

      data echarge1/1.602176462D-19/
      data emasskg1/9.10938188D-31/
      data emassg1/0.510998902D-3/

      if (efieldx.eq.0.0d0.and.efieldy.eq.0.0d0.
     &    and.efieldz.eq.0.0d0) then
        dgamma=0.0d0
        return
      endif

      qdt=icharge*echarge1*dt
      eg=emasskg1*gamma

      px0=vx*eg !SI-units
      py0=vy*eg
      pz0=vz*eg
      pn02=px0*px0+py0*py0+pz0*pz0
      pn0=sqrt(pn02)

      dpx=efieldx*qdt !SI-units
      dpy=efieldy*qdt
      dpz=efieldz*qdt

      px=px0+dpx !SI-units
      py=py0+dpy
      pz=pz0+dpz

      vn=sqrt(vx*vx+vy*vy+vz*vz)
      pn2=px*px+py*py+pz*pz
      pn=sqrt(pn2)

c total momentum and energy are kept!!
      vx=px/pn*vn
      vy=py/pn*vn
      vz=pz/pn*vn

      dtm=0.5d0*dt/(emasskg1*gamma) ! F=dp/dt, m=m_e*gamma, a=F/m*dp/dt/m_e/gamma

      x=x+dtm*dpx
      y=y+dtm*dpy
      z=z+dtm*dpz

      de=(vx*dpx+dpy*vy+dpz*vz)/echarge1/1.0d9 !GeV
      dgamma=de/emassg1

      return
      end
*CMZ :          25/03/2015  09.51.10  by  Michael Scheer
*CMZ :  3.02/00 28/08/2014  15.45.26  by  Michael Scheer
*CMZ :  2.70/00 12/11/2012  11.53.09  by  Michael Scheer
*CMZ :  2.68/03 01/09/2012  16.13.42  by  Michael Scheer
*CMZ :  2.68/02 08/06/2012  09.54.11  by  Michael Scheer
*-- Author : Michael Scheer
      subroutine uradphoton(veln,gamma,bx,by,bz,dgamma,dtim,dpphoton)

*KEEP,gplhint.
!******************************************************************************
!
!      Copyright 2013 Helmholtz-Zentrum Berlin (HZB)
!      Hahn-Meitner-Platz 1
!      D-14109 Berlin
!      Germany
!
!      Author Michael Scheer, Michael.Scheer@Helmholtz-Berlin.de
!
! -----------------------------------------------------------------------
!
!    This program is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy (wave_gpl.txt) of the GNU General Public
!    License along with this program.
!    If not, see <http://www.gnu.org/licenses/>.
!
!    Dieses Programm ist Freie Software: Sie koennen es unter den Bedingungen
!    der GNU General Public License, wie von der Free Software Foundation,
!    Version 3 der Lizenz oder (nach Ihrer Option) jeder spaeteren
!    veroeffentlichten Version, weiterverbreiten und/oder modifizieren.
!
!    Dieses Programm wird in der Hoffnung, dass es nuetzlich sein wird, aber
!    OHNE JEDE GEWAEHRLEISTUNG, bereitgestellt; sogar ohne die implizite
!    Gewaehrleistung der MARKTFAEHIGKEIT oder EIGNUNG FueR EINEN BESTIMMTEN ZWECK.
!    Siehe die GNU General Public License fuer weitere Details.
!
!    Sie sollten eine Kopie (wave_gpl.txt) der GNU General Public License
!    zusammen mit diesem Programm erhalten haben. Wenn nicht,
!    siehe <http://www.gnu.org/licenses/>.
!
!******************************************************************************
*KEND.

c NO WARRANTY

      implicit none

      integer ical,i,nbing1

      double precision veln(3),bmag(3),bx,by,bz,ebeam,elmom,gamma,
     &  bparn,bper(3),bpern,epho,eec,bpervn(3),
     &  dgamma,b2per,dpphoton(3),erad1,
     &  dtim,ec,photons,de,deecg1,eecg1,g1,yrnint10

      real rnrn
      double precision, dimension (:), allocatable ::
     &  xrn,yrn,yrnint,coef,work1,work2,work3,work4

      save

      double precision cgam1,pdum,emassg1,emasse1,clight1,pi1,
     &  eecmaxg1,ecdipkev,hbarev1

      data eecmaxg1/5.0d0/
      data nbing1/1000/
      data hbarev1/6.58211889D-16/
      data pi1/3.14159265358979d0/
      data emassg1/0.510998902D-3/
      data emasse1/0.510998902D6/
      data clight1/2.99792458d8/

      if (ical.eq.0) then

        ecdipkev=3.0d0/2.0d0*hbarev1*(clight1/emasse1)**2/emasse1*1.0d15

        allocate(xrn(nbing1))
        allocate(yrn(nbing1))
        allocate(yrnint(nbing1))
        allocate(coef(nbing1))
        allocate(work1(nbing1))
        allocate(work2(nbing1))
        allocate(work3(nbing1))
        allocate(work4(nbing1))

        deecg1=eecmaxg1/(nbing1-1)

        eecg1=0.0d0
        do i=1,10
          eecg1=eecg1+deecg1/10.0d0
          call util_g1_static(eecg1,g1)
          xrn(i)=eecg1
          yrn(i)=g1/eecg1
        enddo

        yrnint(1)=0.0d0
        do i=2,10
          yrnint(i)=yrnint(i-1)+(yrn(i)+yrn(i-1))
     &    /2.0d0*(xrn(i)-xrn(i-1))
        enddo
        yrnint10=yrnint(10)

        eecg1=0.0d0
        do i=10,nbing1
          eecg1=eecg1+deecg1
          call util_g1_static(eecg1,g1)
          xrn(i)=eecg1
          yrn(i)=g1/eecg1
        enddo

        call util_spline_running_integral(
     &    xrn(10:nbing1),yrn(10:nbing1),nbing1-10+1,yrnint(10:nbing1),
     &    coef,work1,work2,work3,work4)

        yrnint(10)=yrnint10
        yrnint(11:nbing1)=yrnint(11:nbing1)+yrnint10
        yrnint=yrnint/yrnint(nbing1)

        do i=2,nbing1
          if (yrnint(i).le.yrnint(i-1)) then
            stop '*** Error in photon: Bad integration of G1 ***'
          endif
        enddo

        call util_spline_coef(
     &    yrnint,xrn,nbing1,0.0d0,0.0d0,
     &    coef,work1,work2,work3,work4)

        erad1=2.8179380D-15
        cgam1=4.0d0/3.0d0*pi1*erad1/emassg1**3
        !dgamma=-pdum*gamma**2*b**2*dt
        pdum=cgam1/2.0d0/pi1*clight1*(clight1/1.0d9)**2*emassg1

        ical=1

      endif !ical

      bmag(1)=bx
      bmag(2)=by
      bmag(3)=bz

      elmom=emassg1*dsqrt((gamma-1)*(gamma+1)) !GeV
      ebeam=emassg1*gamma !GeV

      bparn=(bmag(1)*veln(1)+bmag(2)*veln(2)+bmag(3)*veln(3))
      bper=bmag-bparn*veln
      b2per=bper(1)**2+bper(2)**2+bper(3)**2
      bpern=sqrt(b2per)
      bpervn=bper/bpern

      ec=ecdipkev*bpern*ebeam**2*1.0d-6 !GeV

      call uradrndm(rnrn)  !S. 39

      !dgamma = pdum * gamma**2 * b2per * dtim
      !dN = 15*sqrt(3)/8 * dE/Ec = 3.2476 * de/ec

      de=pdum*b2per*gamma*ebeam*dtim !GeV

      if (ec.ne.0.0d0) then
        photons=3.2476d0*de/ec !number of photons
      else
        photons=0.0d0
      endif

      call uradrndm(rnrn)

      if(rnrn.le.photons) then

        call uradrndm(rnrn)  !s. 39
        call util_spline_inter(yrnint,xrn,coef,nbing1,
     &    dble(rnrn),eec,-1)
        if (eec.lt.0.0d0) then
          print*,'*** Warning in PHOTON: Negative 
     &    photon energy occured ***'
          print*,'rnrn:',rnrn
          print*,'setting Epho/Ec = 1.e-6'
          eec=1.0d-6
        endif

        epho=eec*ec
        dgamma=-epho/ebeam*gamma
        dpphoton=-veln*epho

      else

        dpphoton=0.0d0
        dgamma=0.0d0

      endif !(rnrn.le.wrad)

      return
      end
*CMZ :          27/03/2015  15.13.24  by  Michael Scheer
*CMZ : 00.00/15 03/09/2012  09.26.37  by  Michael Scheer
*-- Author :    Michael Scheer   10/05/2012
      subroutine util_g1_static(y,g1)

c calculates G1(y) with an estimated precision of about 1.5e-3 for y<=30,
c and 1.5e-2 for y>30.

      implicit none

      integer ical,npoi,ipoi,npoilow,npoihigh,ndatp
      parameter(ndatp=29)

      double precision y,g1,g1_30,c_30,g1_5em5,c_5em5,y_30,
     &  y_5em5,ydum,g1dum,r1

      double precision
     &  ywlow(ndatp),ywhigh(ndatp),coeflow(ndatp),coefhigh(ndatp),
     &  g1wlow(ndatp),g1whigh(ndatp),g1walow(ndatp),
     &  g1wahigh(ndatp),r1low(ndatp),r1high(ndatp),
     &  w1(ndatp),w2(ndatp),w3(ndatp),w4(ndatp),
     &  ystat(ndatp),g1stat(ndatp)

      save

      data ical/0/
      data y_30/30.0d0/
      data y_5em5/5.0d-5/
      data g1_30/6.580794488121591d-013/ !WAVE
      data g1_5em5/7.909860755922665E-002/ !WAVE

* Numerisch mit WAVE berechnet 11.5.2012 (ISPECDIP=2)
        data ystat/
     &    5.0D-005, 7.0D-005, 2.0D-004, 5.0D-004, 1.0D-003,
     &    2.0D-003, 5.0D-003, 1.0D-002, 2.0D-002, 5.0D-002,
     &    0.10D0, 0.20D0, 0.50D0, 1.0D0, 2.0D0,
     &    3.0D0, 4.0D0, 5.0D0, 6.0D0, 7.0D0,
     &    8.0D0, 9.0D0, 10.00D0, 20.00D0, 30.00D0,
     &    40.0D0, 50.0D0, 60.0D0, 70.0D0
     &    /

        data g1stat/
     &    7.909860755922665D-002, 8.846122555950733D-002,
     &    0.125342415210665d0, 0.169701277907238d0, 0.2131391d0,
     &    0.2671962d0, 0.3584969d0, 0.4449725d0, 0.5472394d0,
     &    0.7015719d0, 0.8181855d0, 0.9033860d0, 0.8708191d0,
     &    0.65142282d0, 0.30163590d0, 0.128565710002655d0,
     &    5.282739666852105D-002,
     &    2.12481297D-002, 8.426079715722744D-003,
     &    3.307610970763407D-003, 1.288451614441198D-003,
     &    4.988932935072772D-004, 1.92238264D-004,
     &    1.19686345D-008, 6.58079455D-013, 3.42988745D-017,
     &    1.73478519D-021, 8.60693915D-026, 4.21333348D-030
     &    /

      if (ical.eq.0) then

        c_5em5=g1_5em5/y_5em5**(1.0d0/3.0d0)
        c_30=g1_30/(sqrt(y_30)*exp(-y_30))

        npoilow=0
        npoihigh=0

        do npoi=1,ndatp
          ydum=ystat(npoi)
          if (ydum.ge.y_5em5.and.ydum.le.4.0d0) then !zwei Abfragen wegen 4.0
            npoilow=npoilow+1
          endif
          if (ydum.ge.4.0d0.and.ydum.le.y_30) then
            npoihigh=npoihigh+1
          endif
        enddo

        npoi=ndatp

        npoilow=0
        npoihigh=0
        do ipoi=1,npoi
          ydum=ystat(ipoi)
          g1dum=g1stat(ipoi)
          if (ydum.ge.y_5em5.and.ydum.le.4.0d0) then
            npoilow=npoilow+1
            ywlow(npoilow)=ydum
            g1wlow(npoilow)=g1dum
            g1walow(npoilow)=
     &        391.8d0 * ydum**(1.0d0/3.0d0) * exp(-ydum*0.8307d0)
     &        -192.0d0 * sqrt(ydum) * exp(-ydum*0.7880d0)
          endif
          if (ydum.ge.4.0d0.and.ydum.le.y_30) then
            npoihigh=npoihigh+1
            ywhigh(npoihigh)=ydum
            g1whigh(npoihigh)=g1dum
            g1wahigh(npoihigh)=164.0d0*sqrt(ydum)*EXP(-ydum)
          endif
        enddo

        r1low(1:npoilow)=g1wlow(1:npoilow)/g1walow(1:npoilow)
        r1high(1:npoihigh)=g1whigh(1:npoihigh)/g1wahigh(1:npoihigh)

        call util_spline_coef(ywlow,r1low,npoilow,0.0d0,0.0d0,coeflow,
     &    w1,w2,w3,w4)
        call util_spline_coef(ywhigh,r1high,npoihigh,0.0d0
     &    ,0.0d0,coefhigh,
     &    w1,w2,w3,w4)

        ical=1
      endif

      if (y.le.5.0d-5) then
        g1=c_5em5*y**(1.0d0/3.0d0)
      else if (y.ge.30.0d0) then
        g1=c_30*sqrt(y)*exp(-y)
      else

        if (y.ge.y_5em5.and.y.lt.4.0d0) then
          call util_spline_inter(ywlow,r1low,coeflow,npoilow,y,r1,-1)
          g1=r1*(
     &      391.8d0 * y**(1.0d0/3.0d0) * exp(-y*0.8307d0)
     &      -192.0d0 * sqrt(y) * exp(-y*0.7880d0))
        else if (y.ge.4.0d0.and.y.le.y_30) then
          call util_spline_inter(ywhigh,r1high,
     &         coefhigh,npoihigh,y,r1,-1)
          g1=r1*(164.0d0*sqrt(y)* EXP(-y))
        endif

      endif

      return
9999  stop '*** File wave-g1.dat not found ***'
      end
*CMZ :          19/03/2014  12.14.18  by  Michael Scheer
*CMZ : 00.00/15 03/09/2012  09.26.58  by  Michael Scheer
*CMZ : 00.00/07 05/03/2008  15.43.44  by  Michael Scheer
*CMZ : 00.00/02 14/08/2006  13.22.55  by  Michael Scheer
*-- Author :    Michael Scheer   23/01/2004
      subroutine util_skip_comment_end(lun,ieof)

      implicit none

      integer lun,ieof
      character com

      ieof=0

1     read(lun,'(a)',end=99) com

      if (com.ne.'!'.and.com.ne.'*'.and.com.ne.'#'
     &    .and.com.ne.'%'.and.com.ne.'@') then
        backspace(lun)
      else
        goto 1
      endif

      return

99    ieof=1

      return
      end
*CMZ :          19/03/2014  12.30.26  by  Michael Scheer
*CMZ : 00.00/15 03/09/2012  09.27.13  by  Michael Scheer
*CMZ : 00.00/07 22/03/2010  15.28.00  by  Michael Scheer
*CMZ : 00.00/02 26/03/97  10.23.11  by  Michael Scheer
*CMZ : 00.00/00 10/01/95  15.27.40  by  Michael Scheer
*-- Author :
      SUBROUTINE UTIL_PARABEL(Xin,Yin,A,YP,XOPT,yopt,IFAIL)

C--- CALCULATES A(1),A(2),A(3), THE DERIVATIVES YP(X(1)),YP(X(2)),YP(X(3)),
C    AND THE EXTREMUM (XOPT,A(XOPT)) OF PARABOLA Y=A1+A2*X+A3*X**2
C    FROM COORDINATES OF THE THREE POINTS (X(1),Y(1)),(X(2),Y(2)),(X(3),Y(3))
C

      IMPLICIT NONE

      INTEGER IFAIL

      REAL*8 A(3),X(3),Y(3),DXM,DXP,x0,a1,a2,dxm2,dxp2,dxmax,dymax,
     &  xin(3),yin(3)
      REAL*8 DET,YP(3),XOPT,yopt,a22,fm,fp,f0

      IFAIL=0

c calculate f=a0+a1*(x-x0)+a2*(x-x0)**2
c  = a0 + a1*x - a0*x0 + a2*x**2 - 2*a2*x*x0 + a2*x0**2
c  = a0 + (a2*x0 -a0)*x0 + (a1 - 2*a2*x0 )*x+ a2*x**2

c change system: (x0,s0)->(0,0), i.e.
c calculate f=a1*dx+a2*dx**2
c  df/dx=a1+2*a2*dx_max =! 0, dx_max=-a1/2/a2

      x=xin
      y=yin

      call util_sort_func(3,x,y)

      x0=x(2)
      f0=y(2)

      fm=y(1)-f0
      fp=y(3)-f0

      dxm=x(1)-x0
      dxp=x(3)-x0

c fm=a1*dxm+a2*dxm**2
c fp=a1*dxp+a2*dxp**2

c (dxm dxm2) (a1) = (y(1))
c (dxp dxp2) (a2) = (y(3))

      dxm2=dxm*dxm
      dxp2=dxp*dxp

      det=dxm*dxp2-dxp*dxm2

      if (det.ne.0.0d0) then
        a1=(fm*dxp2-fp*dxm2)/det
        a2=(fp*dxm-fm*dxp)/det
      else
        ifail=1
        return
      endif

      if (a2.ne.0.0d0) then
        dxmax=-a1/(2.0d0*a2)
        dymax=(a1+a2*dxmax)*dxmax
        xopt=x0+dxmax
        yopt=f0+dymax
      endif

c calculate f=f0+a1*dx+a2*dx**2
c = a1*x - a1*x0 + a2*x**2 + a2*x0**2 - 2*a2*x*x0
c  f = f0 + (a2*x0 -a1)*x0 + (a1 - 2*a2*x0 )*x+ a2*x**2

      a22=2.0d0*a2

      a(1)=f0 + (a2*x0 -a1)*x0
      a(2)=a1 - a22*x0
      a(3)=a2

c calculate yp=a1+2*a2*dx

      yp(1)=a1+a22*dxm
      yp(2)=a1
      yp(3)=a1+a22*dxp

      RETURN
      END
*CMZ : 00.00/15 05/01/2012  13.52.39  by  Michael Scheer
*CMZ : 00.00/00 11/01/95  11.41.04  by  Michael Scheer
*-- Author :
      SUBROUTINE UTIL_SORT_FUNC(N,RA,YA)

C--- HEAPSORT ROUTINE; SEE NUMERICAL RECIPES 8.2 S 231
C--- ARRAY YA IS FUNCTION OF RA AND SORTED ACCORDINGLY

      IMPLICIT NONE

      INTEGER N,L,IR,I,J

      REAL*8 RA(N),RRA
      REAL*8 YA(N),YYA

      if (n.lt.2) return

      L=N/2+1
      IR=N
10    CONTINUE
        IF(L.GT.1)THEN
          L=L-1
          RRA=RA(L)
          YYA=YA(L)
        ELSE
          RRA=RA(IR)
          YYA=YA(IR)
          RA(IR)=RA(1)
          YA(IR)=YA(1)
          IR=IR-1
          IF(IR.EQ.1)THEN
            RA(1)=RRA
            YA(1)=YYA
            RETURN
          ENDIF
        ENDIF
        I=L
        J=L+L
20      IF(J.LE.IR)THEN
          IF(J.LT.IR)THEN
            IF(RA(J).LT.RA(J+1))J=J+1
          ENDIF
          IF(RRA.LT.RA(J))THEN
            RA(I)=RA(J)
            YA(I)=YA(J)
            I=J
            J=J+J
          ELSE
            J=IR+1
          ENDIF
        GO TO 20
        ENDIF
        RA(I)=RRA
        YA(I)=YYA
      GO TO 10
      END
*CMZ :          19/03/2014  12.16.04  by  Michael Scheer
*CMZ : 00.00/15 12/10/2013  12.11.47  by  Michael Scheer
*CMZ : 00.00/11 11/02/2011  15.34.09  by  Michael Scheer
*-- Author :    Michael Scheer   11/02/2011
      SUBROUTINE util_spline_running_integral(X,Y,N,resultat
     &                                 ,COEF,WORK1,WORK2,WORK3,WORK4)

C---  CALCULATES RUNNING INTERGRAL OF Y(X) VIA SPLINES

      IMPLICIT NONE

      INTEGER I,N
      REAL*8 X(N),Y(N),resultat(n)
      REAL*8 COEF(N),WORK1(N),WORK2(N),WORK3(N),WORK4(N)

C---  SPLINE-COEFFICIENTS

      CALL UTIL_SPLINE_COEF(X,Y,N,9999.0d0,9999.0d0,COEF,
     &  WORK1,WORK2,WORK3,WORK4)

C--- INTEGRATION

      resultat(1)=0.0D0
      DO I=1,N-1

        resultat(i+1)=resultat(i)
     &          +(X(I+1)-X(I))*0.5D0
     &          *(Y(I)+Y(I+1))
     &          -(X(I+1)-X(I))**3/24.D0
     &          *(COEF(I)+COEF(I+1))

      ENDDO

      RETURN
      END
*CMZ :          19/03/2014  12.06.16  by  Michael Scheer
*CMZ : 00.00/15 03/09/2012  09.27.48  by  Michael Scheer
*CMZ : 00.00/07 12/10/2009  12.17.45  by  Michael Scheer
*CMZ : 00.00/02 14/04/2003  12.46.09  by  Michael Scheer
*CMZ : 00.00/00 10/01/95  15.27.48  by  Michael Scheer
*-- Author :
      SUBROUTINE UTIL_SPLINE_COEF(X,Y,N,YP1,YPN,Y2,AA,BB,CC,C)

C--- CALCULATES SPLINE COEFFICIENTS

C--   INPUT:

C-       N: NUMBER OF X,Y-VALUES
C-       X: ARRAY OF X-VALUES
C-       Y: ARRAY OF Y-VALUES
C-       YP1:  SECOND DERIVATIVE AT FIRST X-VALUE
C-       YPN:  SECOND DERIVATIVE AT LAST X-VALUE

C--   OUPUT:

C-       Y2:   SPLINE-COEFFICIENTS

C--   WORKINGSPACE: AA(N),BB(N),CC(N),C(N)


      IMPLICIT NONE

      INTEGER N,J
      REAL*8  X(N),Y(N),Y2(N),AA(N),BB(N),CC(N),C(N)

      REAL*8 YP1,YPN

      double precision xx(3),yy(3),a(3),yp(3),xopt,yopt
      INTEGER ifail

      IF (N.LT.3) then
        if (abs(yp1).eq.9999.0d0) then
          y2(1)=0.0d0
        else
          y2(1)=yp1
        endif
        if (abs(ypn).eq.9999.0d0) then
          y2(n)=0.0d0
        else
          y2(n)=ypn
        endif
        RETURN
      endif

      if (abs(yp1).eq.9999.0d0) then
        xx=x(1:3)
        yy=y(1:3)
        call UTIL_PARABEL(xx,yy,A,YP,XOPT,yopt,IFAIL)
        if (ifail.eq.0) then
          y2(1)=2.0d0*a(3)
        else
          y2(1)=0.0d0
        endif
      else
        Y2(1)=YP1
      endif

      if (abs(ypn).eq.9999.0d0) then
        xx=x(n-2:n)
        yy=y(n-2:n)
        call UTIL_PARABEL(xx,yy,A,YP,XOPT,yopt,IFAIL)
        if (ifail.eq.0) then
          y2(n)=2.0d0*a(3)
        else
          y2(N)=0.0d0
        endif
      else
        Y2(N)=YPN
      endif

      C(1)=Y2(1)
      C(N)=y2(n)

      BB(1)=1.D0
      CC(1)=0.D0
      CC(N)=1.D0

      DO J=2,N-1
        if(x(j+1).eq.x(j)) then
          write(6,*)
          write(6,*)'*** Error in util_spline_coef: 
     & Intervall of zero length'
          write(6,*)'j, x(j), x(j+1):',j,x(j),x(j+1)
          write(6,*)
          stop
        endif
          AA(J)=(X(J  )-X(J-1))/6.D0
          BB(J)=(X(J+1)-X(J-1))/3.D0
          CC(J)=(X(J+1)-X(J  ))/6.D0
          C(J)=(Y(J+1)-Y(J  ))/(X(J+1)-X(J  ))
     &          -(Y(J  )-Y(J-1))/(X(J  )-X(J-1))
      ENDDO !J

      DO J=2,N-1

          BB(J)=BB(J)-AA(J)*CC(J-1)
           C(J)= C(J)-AA(J)* C(J-1)
C          AA(J)=AA(J)-AA(J)*BB(J-1)

          CC(J)=CC(J)/BB(J)
           C(J)= C(J)/BB(J)
          BB(J)=1.D0

      ENDDO !J

      DO J=N-1,2,-1
         Y2(J)=C(J)-CC(J)*Y2(J+1)
      ENDDO

      RETURN
      END
*CMZ :          19/03/2014  12.14.02  by  Michael Scheer
*CMZ : 00.00/15 03/09/2012  09.28.03  by  Michael Scheer
*CMZ : 00.00/11 07/06/2011  14.38.25  by  Michael Scheer
*CMZ : 00.00/08 15/12/2010  14.05.16  by  Michael Scheer
*CMZ : 00.00/07 07/05/2008  14.28.10  by  Michael Scheer
*CMZ : 00.00/02 25/08/2006  15.27.06  by  Michael Scheer
*CMZ : 00.00/01 23/02/96  14.56.50  by  Michael Scheer
*CMZ : 00.00/00 10/01/95  15.27.54  by  Michael Scheer
*-- Author :
      SUBROUTINE UTIL_SPLINE_INTER(XA,YA,Y2A,N,X,Y,MODE)

C---  INTERPOLATES Y(X) VIA SPLINE

C--   INPUT:

C-       N: NUMBER OF X,Y-VALUES
C-       XA:   ARRAY OF X-VALUES
C-       YA:   ARRAY OF Y-VALUES
C-       YA2:  ARRAY SPLINE COEFFICIENTS
C-       X: Y(X) IS CALCULATED
C-       MODE: CONTROL FLAG:
C-             MODE.GE.0: USE VALUES OF LAST CALL TO START WITH
C-             MODE.LT.0: NEW INITIALIZATION

C--   OUTPUT:

C-       Y: Y(X) IS CALCULATED

      IMPLICIT NONE

      INTEGER NOLD,N,KLO,KHI,KLOLD,K,MODE,NORDER

      REAL*8 Y,X,XA1OLD,XANOLD,H,A,B

      REAL*8 XA(N),YA(N),Y2A(N),EPS,XX

      save klold,nold,xa1old,xanold

      DATA KLOLD/1/,NOLD/-99/
      DATA XA1OLD/-9999.D0/,XANOLD/-9999./

      EPS=ABS(XA(N)-XA(1))/1.0D10
      XX=X

      IF(XA(1).LT.XA(N)) THEN

        IF(XX.LT.XA(1).AND.XX.GT.XA(1)-EPS) THEN
          XX=XA(1)
        ELSE IF(XX.GT.XA(N).AND.XX.LT.XA(N)+EPS) THEN
          XX=XA(N)
        ENDIF

        IF(XX.LT.XA(1).OR.XX.GT.XA(N)) THEN
          WRITE(6,*)'XA(1), XA(N):',XA(1), XA(N)
          WRITE(6,*)'X:',X
          WRITE(6 ,*)'***ERROR IN UTIL_SPLINE_INTER: X OUT OF RANGE ***'
          STOP
        ENDIF

      ELSE

        IF(XX.LT.XA(N).AND.XX.GT.XA(N)-EPS) THEN
          XX=XA(N)
        ELSE IF(XX.GT.XA(1).AND.XX.LT.XA(N)+EPS) THEN
          XX=XA(1)
        ENDIF

        IF(XX.LT.XA(N).OR.XX.GT.XA(1)) THEN
          WRITE(6,*)'XA(1), XA(N):',XA(1), XA(N)
          WRITE(6,*)'X:',X
          WRITE(6 ,*)'***ERROR IN UTIL_SPLINE_INTER: X OUT OF RANGE ***'
          STOP
        ENDIF

      ENDIF

      norder=1
      if (xa(n).lt.xa(1)) then
        norder=-1
      endif

      if (norder.eq.1) then

        IF (MODE.LT.0.OR.KLOLD.GE.N) THEN
          KLO=1
        ELSE IF(NOLD.EQ.N
     &      .AND. XA(1).EQ.XA1OLD
     &      .AND. XA(N).EQ.XANOLD
     &      .AND. XX.GT.XA(KLOLD)
     &      ) THEN
          KLO=KLOLD
        ELSE
          KLO=1
        ENDIF

        IF (XX.LT.XA(KLO+1)) THEN
          KHI=KLO+1
          GOTO 2
        ENDIF

        KHI=N
1       IF (KHI-KLO.GT.1) THEN
          K=(KHI+KLO)/2
          IF(XA(K).GT.XX)THEN
            KHI=K
          ELSE
            KLO=K
          ENDIF
          GOTO 1
        ENDIF

2       H=XA(KHI)-XA(KLO)

        IF (H.le.0.0D0) THEN
          WRITE(6 ,*)'*** ERROR IN UTIL_SPLINE_INTER: BAD INPUT ***'
          STOP
        ENDIF

        A=(XA(KHI)-XX)/H
        B=(XX-XA(KLO))/H
        Y=A*YA(KLO)+B*YA(KHI)+
     &    (A*(A+1.D0)*(A-1.D0)*Y2A(KLO)+B*(B+1.D0)*
     &    (B-1.D0)*Y2A(KHI))*(H**2)/6.D0

      KLOLD=KLO
      NOLD=N
      XA1OLD=XA(1)
      XANOLD=XA(N)

      else !(norder.eq.1) then

        IF (MODE.LT.0.or.nold.ne.n) THEN
          KLO=1
        ELSE IF(
     &      XA(1).EQ.XA1OLD
     &      .AND. XA(N).EQ.XANOLD
     &      .AND. XX.lt.XA(KLOLD)
     &      ) THEN
          KLO=KLOLD
        ELSE
          KLO=1
        ENDIF

        IF (XX.gt.XA(KLO+1)) THEN
          KHI=KLO+1
          GOTO 21
        ENDIF

        KHI=N
11      IF (KHI-KLO.GT.1) THEN
          K=(KHI+KLO)/2
          IF(XA(K).LT.XX)THEN
            KHI=K
          ELSE
            KLO=K
          ENDIF
          GOTO 11
        ENDIF

21      H=XA(KHI)-XA(KLO)

        IF (H.ge.0.0D0) THEN
          WRITE(6 ,*)'*** ERROR IN UTIL_SPLINE_INTER: BAD INPUT ***'
          STOP
        ENDIF

        A=(XA(KHI)-XX)/H
        B=(XX-XA(KLO))/H
        Y=A*YA(KLO)+B*YA(KHI)+
     &    (A*(A+1.D0)*(A-1.D0)*Y2A(KLO)+B*(B+1.D0)*
     &    (B-1.D0)*Y2A(KHI))*(H**2)/6.D0

        KLOLD=KLO
        NOLD=N
        XA1OLD=XA(1)
        XANOLD=XA(N)

      endif !(norder.eq.1) then

      RETURN
      END
*CMZ :  2.68/02 02/07/2012  11.19.55  by  Michael Scheer
*CMZ :  2.66/09 22/03/2010  15.23.05  by  Michael Scheer
*CMZ : 00.00/02 26/03/97  10.23.11  by  Michael Scheer
*CMZ : 00.00/00 10/01/95  15.27.40  by  Michael Scheer
*-- Author : Michael Scheer
      SUBROUTINE parabel_short(x,y,a)
*KEEP,gplhint.
!******************************************************************************
!
!      Copyright 2013 Helmholtz-Zentrum Berlin (HZB)
!      Hahn-Meitner-Platz 1
!      D-14109 Berlin
!      Germany
!
!      Author Michael Scheer, Michael.Scheer@Helmholtz-Berlin.de
!
! -----------------------------------------------------------------------
!
!    This program is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy (wave_gpl.txt) of the GNU General Public
!    License along with this program.
!    If not, see <http://www.gnu.org/licenses/>.
!
!    Dieses Programm ist Freie Software: Sie koennen es unter den Bedingungen
!    der GNU General Public License, wie von der Free Software Foundation,
!    Version 3 der Lizenz oder (nach Ihrer Option) jeder spaeteren
!    veroeffentlichten Version, weiterverbreiten und/oder modifizieren.
!
!    Dieses Programm wird in der Hoffnung, dass es nuetzlich sein wird, aber
!    OHNE JEDE GEWAEHRLEISTUNG, bereitgestellt; sogar ohne die implizite
!    Gewaehrleistung der MARKTFAEHIGKEIT oder EIGNUNG FueR EINEN BESTIMMTEN ZWECK.
!    Siehe die GNU General Public License fuer weitere Details.
!
!    Sie sollten eine Kopie (wave_gpl.txt) der GNU General Public License
!    zusammen mit diesem Programm erhalten haben. Wenn nicht,
!    siehe <http://www.gnu.org/licenses/>.
!
!******************************************************************************
*KEND.

C--- CALCULATES A(1),A(2),A(3)

      IMPLICIT NONE

      REAL*8 A(3),X(3),Y(3),DXM,DXP,x0,a1,a2,dxm2,dxp2
      REAL*8 DET,a22,fm,fp,f0

      x0=x(2)
      f0=y(2)

      fm=y(1)-f0
      fp=y(3)-f0

      dxm=x(1)-x0
      dxp=x(3)-x0

c fm=a1*dxm+a2*dxm**2
c fp=a1*dxp+a2*dxp**2

c (dxm dxm2) (a1) = (y(1))
c (dxp dxp2) (a2) = (y(3))

      dxm2=dxm*dxm
      dxp2=dxp*dxp

      det=dxm*dxp2-dxp*dxm2

      a1=(fm*dxp2-fp*dxm2)/det
      a2=(fp*dxm-fm*dxp)/det
      a22=2.0d0*a2

      a(1)=f0 + (a2*x0 -a1)*x0
      a(2)=a1 - a22*x0
      a(3)=a2

      RETURN
      END
*CMZ :  3.01/02 25/02/2014  14.52.02  by  Michael Scheer
*CMZ :  2.70/06 07/01/2013  13.42.44  by  Michael Scheer
*CMZ :  2.70/05 02/01/2013  13.41.59  by  Michael Scheer
*CMZ :  2.66/20 18/10/2011  08.06.47  by  Michael Scheer
*CMZ :  2.66/18 01/12/2010  10.11.23  by  Michael Scheer
*CMZ :  2.66/13 23/11/2010  09.59.43  by  Michael Scheer
*CMZ :  2.63/00 10/01/2008  12.40.35  by  Michael Scheer
*CMZ :  2.58/01 17/01/2007  10.59.07  by  Michael Scheer
*CMZ :  2.57/04 13/01/2006  11.10.13  by  Michael Scheer
*CMZ :  2.57/03 09/12/2005  11.19.05  by  Michael Scheer
*CMZ :  2.52/15 03/01/2005  15.51.05  by  Michael Scheer
*CMZ :  2.47/22 03/12/2003  10.29.33  by  Michael Scheer
*CMZ :  2.46/02 21/01/2003  16.21.35  by  Michael Scheer
*CMZ :  2.41/09 14/08/2002  17.26.05  by  Michael Scheer
*CMZ :  2.39/00 03/01/2002  12.30.43  by  Michael Scheer
*CMZ :  2.20/01 02/01/2001  11.38.55  by  Michael Scheer
*CMZ :  2.15/00 01/05/2000  11.48.08  by  Michael Scheer
*CMZ :  2.13/04 24/01/2000  17.57.30  by  Michael Scheer
*CMZ :  2.13/00 06/10/99  16.32.34  by  Michael Scheer
*CMZ :  1.03/06 01/07/98  10.29.18  by  Michael Scheer
*CMZ : 00.00/04 29/04/94  17.57.32  by  Michael Scheer
*CMZ : 00.00/00 28/04/94  16.14.41  by  Michael Scheer
*-- Author : Michael Scheer
      SUBROUTINE ZEIT(LUN)
*KEEP,gplhint.
!******************************************************************************
!
!      Copyright 2013 Helmholtz-Zentrum Berlin (HZB)
!      Hahn-Meitner-Platz 1
!      D-14109 Berlin
!      Germany
!
!      Author Michael Scheer, Michael.Scheer@Helmholtz-Berlin.de
!
! -----------------------------------------------------------------------
!
!    This program is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy (wave_gpl.txt) of the GNU General Public
!    License along with this program.
!    If not, see <http://www.gnu.org/licenses/>.
!
!    Dieses Programm ist Freie Software: Sie koennen es unter den Bedingungen
!    der GNU General Public License, wie von der Free Software Foundation,
!    Version 3 der Lizenz oder (nach Ihrer Option) jeder spaeteren
!    veroeffentlichten Version, weiterverbreiten und/oder modifizieren.
!
!    Dieses Programm wird in der Hoffnung, dass es nuetzlich sein wird, aber
!    OHNE JEDE GEWAEHRLEISTUNG, bereitgestellt; sogar ohne die implizite
!    Gewaehrleistung der MARKTFAEHIGKEIT oder EIGNUNG FueR EINEN BESTIMMTEN ZWECK.
!    Siehe die GNU General Public License fuer weitere Details.
!
!    Sie sollten eine Kopie (wave_gpl.txt) der GNU General Public License
!    zusammen mit diesem Programm erhalten haben. Wenn nicht,
!    siehe <http://www.gnu.org/licenses/>.
!
!******************************************************************************
*KEND.

C To determine date and time and write it to logical unit LUN

      IMPLICIT NONE
      INTEGER LUN

      CHARACTER SPACER(50)
*KEEP,datetime.
      character(10) dtday,dttime,dtzone
      integer idatetime(8)
      common/datec/idatetime,dtday,dttime,dtzone
*KEND.

      DATA SPACER/50*' '/

      CALL date_and_time(dtday,dttime,dtzone,idatetime)

      WRITE(LUN,*)
      WRITE(LUN,*)SPACER,dttime(1:2),':',dttime(3:4),':',dttime(5:6),' '
     &,dtday(7:8),'.',dtday(5:6),'.',dtday(3:4)
      WRITE(LUN,*)


      RETURN
      END
