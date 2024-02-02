      program dsmc
c     Integrate the langevin equation equivalent to the FP equation of a molecular gas with nonlineal drag
c     Study the velocity distribution function after a quench from very high temperature.

      implicit none

      integer id,d,iseed
      integer N,i,j,jmax,l,k,jmaxx

      real ran3

      real(8) vmod,varg,pi,rN,x,y,chi2,z0,zeff
      real*8 phi,deltac,deltacx
      real*8 time,a2,a3,cj,c,v2th,vth,cx,cy
      real*8 tboltz2d,val
      real*8 v2kk,v2,v4,temp,temp0,g,v6,tempf,thetai
      real*8 Deltat,tau,tmed,smax,s,sctrl,deltas

      parameter (g=0.1d0,jmax=200,deltac=0.01d0,deltat=1.d-5)
      parameter (N=1*10**6,d=2,temp0=10000.d0,tempf=1.d0)
      parameter (z0=1.772453851d0*dsqrt(tempf),thetai=temp0/tempf)
      parameter (smax=12.5d0,deltas=0.1d0,deltacx=0.05d0,jmaxx=25)

      double precision v(N,d),fdvx(-jmaxx:jmaxx,-jmaxx:jmaxx)
      double precision vcm(d),sigma(d),fdv(jmax),cmoments(10)

      open(1,file='kk.dat',status="unknown")
      open(2,file='fdv_01g_10000T0_1Tf_1e6N_12s_2D.dat'
     *,status="unknown")
      open(3,file='Moments_01g_10000T0_1Tf_1e6N_12s_2D.dat'
     *,status="unknown")
      open(4,file='fdvx_01g_10000T0_1Tf_1e6N_12s_2D.dat'
     *,status="unknown")

      pi=4.d0*datan(1.d0)
      rN=dble(N)
      iseed=-123456

!     Initialize some parameters
      time=0.d0
      sctrl=deltas

!     initially: gaussian velocity distribution with T=temp0
      do i=1,N
13      continue
        x=ran3(iseed)
        if (x.eq.0) goto 13
        y=ran3(iseed)
        vmod=dsqrt(-2.d0*temp0*dlog(x))
        varg=2.d0*pi*y
        v(i,1)=vmod*dcos(varg)
        v(i,2)=vmod*dsin(varg)
c14      x=ran3(iseed)
c        if (x.eq.0) goto 14
c        y=ran3(iseed)
c        vmod=dsqrt(-2.d0*temp0*dlog(x))
c        varg=2.d0*pi*y
c        v(i,3)=vmod*dcos(varg)
      enddo

!     Velocity of the Center of mass and change to the CM frame
      do id=1,d
        vcm(id)=0.d0
      enddo

      do id=1,d
        do i=1,N
            vcm(id)=vcm(id)+v(i,id)
        enddo

        vcm(id)=vcm(id)/rN
      enddo

      do i=1,N
        do id=1,d
            v(i,id)=v(i,id)-vcm(id)
        enddo
      enddo

c     Compute Temperature and cumulants
      v2=0.d0
      v4=0.d0
      v6=0.d0
      do i=1,N
            v2kk=0.d0
            do id=1,d
              v(i,id)=v(i,id)-vcm(id)
              v2kk=v2kk+v(i,id)**2.d0
            enddo
            v2=v2+v2kk
            v4=v4+v2kk**2
            v6=v6+v2kk**3
      enddo
      v2=v2/rN
      temp=v2/d
      v4=v4/rN
      v6=v6/rN/(2.d0*temp)**3
      a2=v4/temp**2.d0/dble(d*(d+2))-1.d0
      a3=1.d0+3.d0*a2-8.d0/d*v6/(d+2.d0)/(d+4.d0)
      tau=z0*time
      write(1,199) tau,g*thetai*tau,temp,temp/temp0,a2,a3

!----------------------------------------------------
!     Relaxation
!----------------------------------------------------
      do while (s.lt.smax)

       time=time+deltat
       tau=z0*time
       s=tau*g*thetai

!     Apply the nonlinear drag and kick all the particles at random
         do i=1,N
           v2=v(i,1)**2+v(i,2)**2
           zeff=z0*(1.d0-2.d0*g+g*v2/tempf)
11         continue
           v(i,1)=v(i,1)*(1.d0-zeff*deltat)
           v(i,2)=v(i,2)*(1.d0-zeff*deltat)
c           v(i,3)=v(i,3)*(1.d0-zeff*deltat)
           chi2=2*z0*tempf*(1.d0+g*v2/tempf)
           val=dsqrt(chi2*deltat)
           x=ran3(iseed)
           if (x.eq.0.d0) goto 11
           y=ran3(iseed)
           v(i,1)=v(i,1)+val*dsqrt(-2.d0*dlog(x))*dcos(2.d0*pi*y)
           v(i,2)=v(i,2)+val*dsqrt(-2.d0*dlog(x))*dsin(2.d0*pi*y)
c12         x=ran3(iseed)
c           if (x.eq.0) goto 12
c           y=ran3(iseed)
c           v(i,3)=v(i,3)+val*dsqrt(-2.d0*dlog(x))*dcos(2.d0*pi*y)
         enddo

c     Compute Temperature every deltas
        if (s.ge.sctrl) then
         do id=1,d
          vcm(id)=0.d0
          do i=1,N
            vcm(id)=vcm(id)+v(i,id)
          enddo
          vcm(id)=vcm(id)/rN
         enddo

         v2=0.d0
         v4=0.d0
         v6=0.d0
         do i=1,N
            v2kk=0.d0
            do id=1,d
              v(i,id)=v(i,id)-vcm(id)
              v2kk=v2kk+v(i,id)**2.d0
            enddo
            v2=v2+v2kk
            v4=v4+v2kk**2
            v6=v6+v2kk**3
         enddo
         v2=v2/rN
         temp=v2/d
         v4=v4/rN
         v6=v6/rN/(2.d0*temp)**3
         a2=v4/temp**2.d0/dble(d*(d+2))-1.d0
         a3=1.d0+3.d0*a2-8.d0/d*v6/(d+2.d0)/(d+4.d0)
         write(1,199) tau,s,temp,temp0/temp,a2,a3
         sctrl=sctrl+deltas
        endif

      enddo
C     End of  relaxation


c     Compute fdv
      fdv=0.d0
      cmoments=0.d0
      v2th=2.d0*v2/d
      do i=1,N
            v2kk=0.d0
            do id=1,d
              v2kk=v2kk+v(i,id)**2.d0
            enddo
            c=dsqrt(v2kk/v2th)
            j=int(c/deltac)+1
            if (j.le.jmax) fdv(j)=fdv(j)+1.d0

c           Moments
            do 22 l=1,10
            cmoments(l)=c**l+cmoments(l)
22          continue

      enddo
      fdv=fdv/rN/deltac
      do j=1,jmax
        cj=deltac*(j-0.5d0)
      write(2,199)  cj, fdv(j)
      enddo

      cmoments=cmoments/rN
      do l=1,10
       write(3,*)  l,cmoments(l)
      enddo

c     Probability Density of cx and cy
      fdvx=0.d0
      vth=dsqrt(v2th)
      do i=1,N
        cx=v(i,1)/vth
        j=int(cx/deltacx)
        if(cx.gt.0) j=j+1
        cy=v(i,2)/vth
        k=int(cy/deltacx)
        if(cy.gt.0) k=k+1
        if ((abs(j).le.jmaxx).and.(abs(k).le.jmaxx)) then
         fdvx(j,k)=fdvx(j,k)+1.d0
        endif

      enddo
      fdvx=fdvx/rN/deltacx
      do j=-jmaxx,jmaxx
        cx=deltacx*(j-0.5d0)
        do k=-jmaxx,jmaxx
         cy=deltacx*(k-0.5d0)
         write(4,199)  cx,cy,fdvx(j,k)
        enddo
      enddo

      close(30)
199   format(10(e14.7,2x))
20    format(5(e14.7,2x),2x,i8)

      stop
      end


c====================================================================
!     Subrutina Ran3 numeros aleatorios
C====================================================================
C       PROGRAM: ran3.f
C       TYPE   : function
C       PURPOSE: generate random numbers
C       VERSION: 17 June 94
C       COMMENT: Initialize idum with negative integer
C======================================================================
        real function ran3(idum)
        integer mbig,mseed,mz,ma,mj,mk,i,ii,k,inext,inextp,iff,idum
        double precision fac
        Parameter (mbig=1000000000,Mseed=161803398,Mz=0,fac=1./Mbig)
        Dimension MA(55)
        save
        if (idum.lt.0.or.iff.eq.0) then
       iff=1
       mj=mseed-iabs(idum)
       mj=mod(mj,mbig)
       ma(55)=mj
       mk=1
       do 11 i=1,54
          ii=mod(21*i,55)
          ma(ii)=mk
          mk=mj-mk
          if (mk.lt.mz) mk=mk+mbig
          mj=ma(ii)
 11    continue
       do 13 k=1,4
          do 12 i=1,55
         ma(i)=ma(i)-ma(1+mod(i+30,55))
         if (ma(i).lt.mz) ma(i)=ma(i)+mbig
 12       continue
 13    continue
       inext=0
       inextp=31
       idum=1
    end if
    inext=inext+1
    if (inext.eq.56) inext=1
    inextp=inextp+1
    if (inextp.eq.56) inextp=1
    mj=ma(inext)-ma(inextp)
    if (mj.lt.mz) mj=mj+mbig
    ma(inext)=mj
    ran3=mj*fac
    return
    end
