program main
 implicit none

 ! parameters 
 integer, parameter :: nrep=10

 ! other
 integer, dimension(1:nrep) :: temperature
 integer, dimension(1:5) :: dum
 integer, dimension(:,:), allocatable :: repwalk
 integer, dimension(:), allocatable :: nlist, temp
 real(8) :: sumt, sumt2
 character(len=10) :: dum4
 character(len=7) :: title
 logical :: eof
 integer :: ierr, cnt, nfr, ifr, rep, i, j, k, l

 ! initialize
 nfr=0
 allocate(nlist(0:nfr))
 cnt=0

 ! read in temperatures
 open(24,file='RepTemp.dat')
 do rep=1,nrep,1
  read(24,*) temperature(rep)
 end do  
 close(24)

 ! get number of exchanges per frame
 open(24,file='RunStatus')
 do
  read(24,'(1a7)',iostat=ierr) title
  eof=.false.
  if(ierr.ne.0) eof=.true.
  if((title.eq.'Replica').or.(eof.eqv..true.)) then
   nlist(nfr)=cnt
   if(eof) exit
   allocate(temp(0:nfr))
   temp=nlist
   deallocate(nlist)
   nfr=nfr+1
   allocate(nlist(0:nfr))
   nlist(0:nfr-1)=temp
   deallocate(temp)
   cnt=0
  else
   cnt=cnt+1
  end if
 end do
 close(24) 

 ! initialize repwalk
 allocate(repwalk(0:nfr,1:nrep))
 do rep=1,nrep,1
  repwalk(0,rep)=rep
 end do

 ! get repwalk
 l=0
 open(24,file='RunStatus')
 open(25,file='repwalk.dat')
 write(25,'(40i3)') repwalk(0,:)
 do ifr=1,nfr,1
  if((nlist(ifr).ne.(nrep/2)).and.(nlist(ifr).ne.((nrep/2)-1))) then
   print*, 'Error: nlist incorrect.', nlist(ifr); stop
  end if
  !print*, ifr
  k=0
  read(24,'()')
  repwalk(ifr,:)=repwalk(ifr-1,:)
  do i=1,nlist(ifr),1
   read(24,*) dum(1:3), dum4, dum(5)
   !print*, dum(1:3), dum4, dum(5)
   k=k+dum(5)
   if(dum(5).eq.1) then
    dum(1)=repwalk(ifr,dum(2))
    repwalk(ifr,dum(2))=repwalk(ifr,dum(3))
    repwalk(ifr,dum(3))=dum(1)
   end if
  end do 
  write(25,'(40i3)') repwalk(ifr,:)
  !if(k.eq.0) then
  ! print*, ifr, k
  !end if
 end do
 close(24)
 close(25)
 !print*, l, real(nfr,8) 

 ! get hansmann
 open(24,file='hansmann.dat')
 do i=1,nrep,1
  sumt=0.d0
  sumt2=0.d0
  do j=1,nrep,1
   cnt=0
   do ifr=0,nfr,1
   !do ifr=4001,nfr,1
    if(repwalk(ifr,i).eq.j) cnt=cnt+1
   end do
   sumt=sumt+cnt
   sumt2=sumt2+cnt*cnt
  end do
  write(24,'(2x,1i3,2x,1f15.7)') temperature(i), 1.d0-sqrt(sumt2)/sumt
 end do
 close(24)
end program main
