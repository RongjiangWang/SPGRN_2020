      subroutine getgrn(depth,lats,lons,latr,lonr,
     &            is,m0,strike,dip,rake,mtt,mpp,mrr,mtp,mpr,mrt,
     &            iout,t0)
      use ggalloc
      implicit none
c
      integer*4 is,iout
      real*8 depth,lats,lons,latr,lonr,t0
      real*8 m0,strike,dip,rake,mtt,mpp,mrr,mtp,mpr,mrt
c
      integer*4 i,j,ir,irsel,ig,igsel,offset
      integer*4 iz1,iz2,ir1,ir2,i1,i2
      real*4 t0r4
      real*8 rmin,dmin,rn,re,z,z1,z2,d,d1,d2
      real*8 expl,clvd,ss12,ss11,ds31,ds23
      real*8 dis,azi,csa,ssa,cs2a,ss2a,t
      real*4 a,b,tp0,tpg,dtp
c
      integer*4 fseek
c
      igsel=1
      dmin=dabs(depth-deps(1))
      do ig=2,ng
        if(dmin.gt.dabs(depth-deps(ig)))then
          igsel=ig
          dmin=dabs(depth-deps(ig))
        endif
      enddo
c
      call disazi(6371.d0,lats,lons,latr,lonr,rn,re)
c
c     dis = epicentral distance
c     azi = azimuth (from south to east) of receiver relative to source
c
      dis=dsqrt(rn*rn+re*re)
      azi=datan2(re,-rn)
      ssa=dsin(azi)
      csa=dcos(azi)
      ss2a=dsin(2.d0*azi)
      cs2a=dcos(2.d0*azi)
c
      if(dis.lt.r(1).or.dis.gt.r(nr))then
        stop ' Error in getgrn: distance exceeds the GF range!'
      else if(depth.lt.ztp1.or.depth.gt.ztp2)then
        stop ' Error in getgrn: depth exceeds the tp_table range!'
      else if(dis.lt.dtp1.or.dis.gt.dtp2)then
        stop ' Error in getgrn: distance exceeds the tp_table range!'
      endif
      irsel=1
      rmin=dabs(dis-r(1))
      do ir=2,nr
        if(rmin.gt.dabs(dis-r(ir)))then
          irsel=ir
          rmin=dabs(dis-r(ir))
        endif
      enddo
c
      open(20,file=grnfile(igsel),form='unformatted',status='old')
      offset=(irsel-1)*(12+(nt*4+8)*10)
      if(offset.gt.0)i=fseek(20,offset,0)
      read(20)t0r4
      t0=dble(t0r4)
      read(20)(rex(i),i=1,nt)
      read(20)(tex(i),i=1,nt)
      read(20)(rss(i),i=1,nt)
      read(20)(tss(i),i=1,nt)
      read(20)(pss(i),i=1,nt)
      read(20)(rds(i),i=1,nt)
      read(20)(tds(i),i=1,nt)
      read(20)(pds(i),i=1,nt)
      read(20)(rcl(i),i=1,nt)
      read(20)(tcl(i),i=1,nt)
      close(20)
c
      if(is.eq.2)then
        call moments(m0,strike,dip,rake,
     &               mtt,mpp,mrr,mtp,mpr,mrt)
      endif
      expl=(mtt+mpp+mrr)/3.d0
      clvd=mrr-expl
      ss12=mtp
      ss11=(mtt-mpp)/2.d0
      ds31=mrt
      ds23=mpr
c
      do i=1,nt
        vrtp(i,1)=expl*dble(rex(i))+clvd*dble(rcl(i))
     &           +(ss12*ss2a+ss11*cs2a)*dble(rss(i))
     &           +(ds31*csa+ds23*ssa)*dble(rds(i))
        vrtp(i,2)=expl*dble(tex(i))+clvd*dble(tcl(i))
     &           +(ss12*ss2a+ss11*cs2a)*dble(tss(i))
     &           +(ds31*csa+ds23*ssa)*dble(tds(i))
        vrtp(i,3)=(ss12*cs2a-ss11*ss2a)*dble(pss(i))
     &           +(ds31*ssa-ds23*csa)*dble(pds(i))
      enddo
      if(iout.eq.0)then
        do j=1,3
          vrtp(1,j)=0.d0
          do i=2,nt
            vrtp(i,j)=vrtp(i-1,j)+vrtp(i,j)*dt
          enddo
        enddo
      else if(iout.eq.2)then
        do j=1,3
          do i=nt,2,-1
            vrtp(i,j)=(vrtp(i,j)-vrtp(i-1,j))/dt
          enddo
          vrtp(1,j)=vrtp(2,j)
        enddo
      endif
c
c     Green function travel times
c
      z=deps(igsel)
      d=r(irsel)
c
      iz1=1+idint((z-ztp1)/dztp)
      iz2=iz1+1
      ir1=1+idint((d-dtp1)/ddtp)
      ir2=ir1+1
c
      z1=ztp1+dble(iz1-1)*dztp
      z2=z1+dztp
      d1=dtp1+dble(ir1-1)*ddtp
      d2=d1+ddtp
c
      tpg=((d2-d)*(z-z1)*tp(ir1,iz2)
     &    +(d-d1)*(z-z1)*tp(ir2,iz2)
     &    +(d2-d)*(z2-z)*tp(ir1,iz1)
     &    +(d-d1)*(z2-z)*tp(ir2,iz1))/(dztp*ddtp)
c
      z=depth
      d=dis
c
      iz1=1+idint((z-ztp1)/dztp)
      iz2=iz1+1
      ir1=1+idint((d-dtp1)/ddtp)
      ir2=ir1+1
c
      z1=ztp1+dble(iz1-1)*dztp
      z2=z1+dztp
      d1=dtp1+dble(ir1-1)*ddtp
      d2=d1+ddtp
c
      tp0=((d2-d)*(z-z1)*tp(ir1,iz2)
     &    +(d-d1)*(z-z1)*tp(ir2,iz2)
     &    +(d2-d)*(z2-z)*tp(ir1,iz1)
     &    +(d-d1)*(z2-z)*tp(ir2,iz1))/(dztp*ddtp)
c
      dtp=tp0-tpg
c
      do i=1,nt
        t=dble(i-1)*dt
        t=dmin1(dmax1(0.d0,t-dtp),twin)
        i1=min0(1+idint(t/dt),nt)
        i2=min0(i1+1,nt)
        b=dmod(t/dt,1.d0)
        a=1.0-b
        do j=1,3
          dvaseis(i,j)=a*vrtp(i1,j)+b*vrtp(i2,j)
        enddo
      enddo
c
      return
      end