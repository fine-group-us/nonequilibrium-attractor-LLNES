      program   brownian3D_LLNES
c     Integrate the langevin equation of a brownian particle in a 3D nonlinear potential
c     Compute the probability density  after a quench from very high temperature.

      implicit none

      integer d,iseed,k,j
      integer N,i,jmax,l

      real ran3

      real(8) xmod,xarg,pi,rN,y1,y2,chi,cx,cy,cz
      real*8 xmed,ymed,x0,y0,z0,zmed
      real*8 g,tempf,deltacx,rmed,temp0
      real*8 time,Deltat,tllnes

      parameter (g=1.d0,deltat=1.d-6)
      parameter (N=10**5,d=3,tempf=1.d-3,temp0=1.d4)
      parameter (tllnes=1.d-5)
      parameter (deltacx=0.05d0,jmax=50)

      double precision x(N,d),y(d)
      double precision fdxy(-jmax:jmax,-jmax:jmax,-jmax:jmax)


      open(4,file='fdxyz_1d-3Tf_1d5N_1d-5t.dat'
     *,status="unknown")

c      open(4,file="test1d-2.dat",status="unknown")

!     Initialize some parameters
      pi=4.d0*datan(1.d0)
      rN=dble(N)
      iseed=-123456
      time=0.d0
      chi=dsqrt(2.d0*tempf*deltat/g)

!     Initial condition: gaussian distribution at temp0
      do i=1,N
13      continue
        y1=ran3(iseed)
        if (y1.eq.0) goto 13
        y2=ran3(iseed)
        xmod=dsqrt(-2.d0*temp0*dlog(y1))
        xarg=2.d0*pi*y2
        x(i,1)=xmod*dcos(xarg)
        x(i,2)=xmod*dsin(xarg)

12      y1=ran3(iseed)
        if (y1.eq.0) goto 12
        y2=ran3(iseed)
        xmod=dsqrt(-2.d0*temp0*dlog(y1))
        xarg=2.d0*pi*y2
        x(i,3)=xmod*dcos(xarg)
      enddo

!----------------------------------------------------
!     Relaxation until LLNES
!----------------------------------------------------
      time=time+deltat
      do while (time.lt.tllnes)
         do i=1,N
14         continue
           y1=ran3(iseed)
           if (y1.eq.0) goto 14
           y2=ran3(iseed)
           xmod=dsqrt(-2.d0*dlog(y1))
           xarg=2.d0*pi*y2
           y(1)=xmod*dcos(xarg)
           y(2)=xmod*dsin(xarg)
15         y1=ran3(iseed)
           if (y1.eq.0) goto 15
           y2=ran3(iseed)
           xmod=dsqrt(-2.d0*dlog(y1))
           xarg=2.d0*pi*y2
           y(3)=xmod*dsin(xarg)
           x0=x(i,1)
           y0=x(i,2)
           z0=x(i,3)
           x(i,1)=x0-(x0**3+x0*y0**2+3.d0*x0*z0**2)*deltat
           x(i,1)=x(i,1)+y(1)*chi
           x(i,2)=y0-(y0**3+x0**2*y0+3.d0*y0*z0**2)*deltat
           x(i,2)=x(i,2)+y(2)*chi
           x(i,3)=z0-(9.d0*z0**3+3.d0*x0**2*z0+3.d0*y0**2*z0)*deltat
           x(i,3)=x(i,3)+y(3)*chi
         enddo
        time=time+deltat
      enddo
C     End of  relaxation

c     Compute <r(t)>
      rmed=0.d0
      xmed=0.d0
      ymed=0.d0
      zmed=0.d0
      do i=1,N
         xmed=xmed+dabs(x(i,1))
         ymed=ymed+dabs(x(i,2))
         zmed=zmed+dabs(x(i,3))
         rmed=rmed+dsqrt(x(i,1)**2+x(i,2)**2+x(i,3)**2)
      enddo
      rmed=rmed/rN
      xmed=xmed/rN
      ymed=ymed/rN
      zmed=zmed/rN


c     Probability Density of cx, cy and cz
      fdxy=0.d0
      do i=1,N
        cx=x(i,1)/xmed
        j=int(cx/deltacx)
        if(cx.gt.0) j=j+1
        cy=x(i,2)/ymed
        k=int(cy/deltacx)
        if(cy.gt.0) k=k+1
        cz=x(i,3)/zmed
        l=int(cz/deltacx)
        if(cz.gt.0) l=l+1

      if((abs(j).le.jmax).and.(abs(k).le.jmax).and.(abs(l).le.jmax))then
          fdxy(j,k,l)=fdxy(j,k,l)+1.d0
        endif
      enddo

      fdxy=fdxy/rN/deltacx
      do j=-jmax,jmax
        cx=deltacx*(j-0.5d0)
        do k=-jmax,jmax
         cy=deltacx*(k-0.5d0)
         do l=-jmax,jmax
           cz=deltacx*(l-0.5d0)
           write(4,199)  cx,cy,cz,fdxy(j,k,l)
         enddo
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
