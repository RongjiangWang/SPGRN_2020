      program ggmain
      use ggalloc
      implicit none
c
      integer*4 i,j,m,ia,ib,ig,ir,flen,is,irec,iout,ntout,ierr
      real*8 latr1,lonr1,latr2,lonr2,dlatr,dlonr,latr,lonr
      real*8 depth,lats,lons,m0,strike,dip,rake
      real*8 t,t0,t1,t2,a,b
      real*8 mtt,mpp,mrr,mtp,mpr,mrt
      character*1 dva
      character*80 outfile,outcmp(3),stdgrndir,tptable,tstable
c
      write(*,'(a,$)')' space-time domain Green function database: '
      read(*,'(a)')stdgrndir
      write(*,'(a,$)')' travel time table for the first P arrival: '
      read(*,'(a)')tptable
      write(*,'(a)')' select input format of focal mechanism'
      write(*,'(a,$)')' (1 = moment tensor, 2 = double couple): '
      read(*,*)is
      if(is.eq.1)then
        write(*,'(a,$)')' Mtt,Mpp,Mrr,Mtp,Mpr,Mrt[Nm]: '
        read(*,*)mtt,mpp,mrr,mtp,mpr,mrt
      else
        write(*,'(a,$)')' moment[Nm], strike, dip, rake[deg]: '
        read(*,*)m0,strike,dip,rake
      endif
      write(*,'(a,$)')' source depth [km]: '
      read(*,*)depth
      write(*,'(a,$)')' source location (lat, lon)[deg]: '
      read(*,*)lats,lons
      write(*,'(a,$)')' number of samples of receiver profile: '
      read(*,*)nrec
      write(*,'(a,$)')' start receiver location (lat, lon)[deg]: '
      read(*,*)latr1,lonr1
      write(*,'(a,$)')' end receiver location (lat, lon)[deg]: '
      read(*,*)latr2,lonr2
      write(*,'(a,$)')' output in dis(0), vel(1) or acc(2): '
      read(*,*)iout
      write(*,'(a,$)')' output file (without extension): '
      read(*,'(a)')outfile
c
      allocate(stdtxt(nrec),stat=ierr)
      if(ierr.ne.0)stop 'error in ggmain: stdtxt not allocated!'
      allocate(tred(nrec),stat=ierr)
      if(ierr.ne.0)stop 'error in ggmain: tred not allocated!'
      allocate(dvaout(nrec),stat=ierr)
      if(ierr.ne.0)stop 'error in ggmain: dvaout not allocated!'
c
      do m=1,80
        if(outfile(m:m).eq.' ')goto 10
      enddo
10    m=m-1
      outcmp(1)=outfile(1:m)//'_z.dat'
      outcmp(2)=outfile(1:m)//'_r.dat'
      outcmp(3)=outfile(1:m)//'_t.dat'
c
      if(iout.eq.0)then
        dva='d'
      else if(iout.eq.1)then
        dva='v'
      else if(iout.eq.2)then
        dva='a'
      endif
      do irec=1,nrec
        stdtxt(irec)(1:1)=dva
        i=irec/1000
        stdtxt(irec)(2:2)=char(ichar('0')+i)
        i=mod(irec,1000)/100
        stdtxt(irec)(3:3)=char(ichar('0')+i)
        i=mod(irec,100)/10
        stdtxt(irec)(4:4)=char(ichar('0')+i)
        i=mod(irec,10)
        stdtxt(irec)(5:5)=char(ichar('0')+i)
	enddo
c
      do flen=80,1,-1
        if(stdgrndir(flen:flen).ne.' ')goto 100
      enddo
100   if(stdgrndir(flen:flen).ne.'/'.or.
     &   stdgrndir(flen:flen).ne.'\')then
        stdgrndir=stdgrndir(1:flen)//'/'
        flen=flen+1
        if(flen.gt.80)then
          stop 'too long name of Green function database!'
        endif 
      endif
c
      open(10,file=stdgrndir(1:flen)//'GreenInfo.dat',status='old')
      call skipdoc(10)
      read(10,*)twin,dt,nt
c
      allocate(dvaseis(nt,3),stat=ierr)
      if(ierr.ne.0)stop 'error in ggmain: dvaseis not allocated!'
      allocate(dvaswap(nt,3,nrec),stat=ierr)
      if(ierr.ne.0)stop 'error in ggmain: dvaswap not allocated!'
      allocate(vrtp(nt,3),stat=ierr)
      if(ierr.ne.0)stop 'error in ggmain: vrtp not allocated!'
      allocate(rex(nt),stat=ierr)
      if(ierr.ne.0)stop 'error in ggmain: rex not allocated!'
      allocate(tex(nt),stat=ierr)
      if(ierr.ne.0)stop 'error in ggmain: tex not allocated!'
      allocate(rss(nt),stat=ierr)
      if(ierr.ne.0)stop 'error in ggmain: rss not allocated!'
      allocate(tss(nt),stat=ierr)
      if(ierr.ne.0)stop 'error in ggmain: tss not allocated!'
      allocate(pss(nt),stat=ierr)
      if(ierr.ne.0)stop 'error in ggmain: pss not allocated!'
      allocate(rds(nt),stat=ierr)
      if(ierr.ne.0)stop 'error in ggmain: rds not allocated!'
      allocate(tds(nt),stat=ierr)
      if(ierr.ne.0)stop 'error in ggmain: tds not allocated!'
      allocate(pds(nt),stat=ierr)
      if(ierr.ne.0)stop 'error in ggmain: pds not allocated!'
      allocate(rcl(nt),stat=ierr)
      if(ierr.ne.0)stop 'error in ggmain: rcl not allocated!'
      allocate(tcl(nt),stat=ierr)
      if(ierr.ne.0)stop 'error in ggmain: tcl not allocated!'
c
      call skipdoc(10)
      read(10,*)nr
c
      allocate(r(nr),stat=ierr)
c
      read(10,*)(r(i),i=1,nr)
      call skipdoc(10)
      read(10,*)ng
c
      allocate(deps(ng),stat=ierr)
      allocate(grnfile(ng),stat=ierr)
c
      read(10,*)((deps(ig),grnfile(ig)),ig=1,ng)
      close(10)
c
      open(21,file=tptable,status='old')
      call skipdoc(21)
      read(21,*)ztp1,ztp2,nztp
      dztp=(ztp2-ztp1)/dble(nztp-1)
      call skipdoc(21)
      read(21,*)dtp1,dtp2,ndtp
      ddtp=(dtp2-dtp1)/dble(ndtp-1)
c
      allocate(tp(ndtp,nztp),stat=ierr)
      if(ierr.ne.0)stop 'error in ggmain: tp not allocated!'
c
      call skipdoc(21)
      read(21,*)((tp(i,j),j=1,nztp),i=1,ndtp)
      close(21)
c
      if(nrec.eq.1)then
        dlatr=0.d0
        dlonr=0.d0
      else
        dlatr=(latr2-latr1)/dble(nrec-1)
        dlonr=(lonr2-lonr1)/dble(nrec-1)
      endif
      do irec=1,nrec
        latr=latr1+dble(irec-1)*dlatr
        lonr=lonr1+dble(irec-1)*dlonr
        call getgrn(depth,lats,lons,latr,lonr,
     &              is,m0,strike,dip,rake,mtt,mpp,mrr,mtp,mpr,mrt,
     &              iout,t0)
        tred(irec)=t0
        do j=1,3
          do i=1,nt
            dvaswap(i,j,irec)=dvaseis(i,j)
          enddo
        enddo
      enddo
c
      twin=dble(nt-1)*dt
      t1=tred(1)
      t2=tred(1)
      do irec=2,nrec
        t1=dmin1(t1,tred(irec))
        t2=dmax1(t2,tred(irec))
      enddo
      t2=t2+twin
      ntout=1+idint((t2-t1)/dt)
c
      do j=1,3
        open(30,file=outcmp(j),status='unknown')
        write(30,'(a,$)')'      time'
        do irec=1,nrec-1
          write(30,'(a14,$)')stdtxt(irec)
        enddo
        write(30,'(a14)')stdtxt(nrec)
        do i=1,ntout
          t=t1+dble(i-1)*dt
          write(30,'(f10.2,$)')t
          do irec=1,nrec
            if(t.lt.tred(irec).or.t.gt.tred(irec)+twin)then
              dvaout(irec)=0.d0
            else
              ia=1+idint((t-tred(irec))/dt)
              ib=min0(ia+1,nt)
              b=dmod((t-tred(irec))/dt,1.d0)
              a=1.d0-b
              dvaout(irec)=a*dvaswap(ia,j,irec)+b*dvaswap(ib,j,irec)
            endif
          enddo
          do irec=1,nrec-1
            write(30,'(E14.6,$)')dvaout(irec)
          enddo
          write(30,'(E14.6)')dvaout(nrec)
        enddo
        close(30)
      enddo
c
      stop
      end
