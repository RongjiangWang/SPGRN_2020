      module ggalloc
      integer*4 nt,nr,nrec,ng,nztp,ndtp
      real*8 twin,dt
      real*8 ztp1,ztp2,dztp,dtp1,dtp2,ddtp
c===================================================================
c     allocatable variables
c===================================================================
      real*8, allocatable:: r(:)
      real*8, allocatable:: deps(:)
      real*8, allocatable:: tp(:,:),tred(:)
      character*80, allocatable:: grnfile(:)*80
      character*5, allocatable:: stdtxt(:)
      real*8, allocatable:: dvaseis(:,:),dvaswap(:,:,:),dvaout(:)
      real*8, allocatable:: vrtp(:,:)
      real*4, allocatable:: rex(:)
      real*4, allocatable:: tex(:)
      real*4, allocatable:: rss(:)
      real*4, allocatable:: tss(:)
      real*4, allocatable:: pss(:)
      real*4, allocatable:: rcl(:)
      real*4, allocatable:: tcl(:)
      real*4, allocatable:: rds(:)
      real*4, allocatable:: tds(:)
      real*4, allocatable:: pds(:)
c
      end module
