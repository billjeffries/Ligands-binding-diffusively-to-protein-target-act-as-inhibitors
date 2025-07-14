program compute_cluster
  implicit none

  integer, parameter :: long = selected_real_kind(15, 307)
  ! real(long), parameter :: rc = 1.d0
  integer, parameter :: nclmax = 10000, nstr = 10000

  type cluster
    integer, dimension(nstr) :: list
    integer :: centroid
    real :: sumsq
    real :: sumnat
    real :: sumnonnat
  end type cluster
  type(cluster), dimension(1:nclmax) :: cl

  logical, dimension(nstr, nstr) :: cxn
  !real, dimension(nstr, nstr) :: rmsd
  integer, dimension(nstr) :: tr, fr, pep, tempcl, maxcl, used, nat, nonnat
  real, dimension(nstr,3) :: xyz !used for making output xyz
  character(len=10) :: argv
  real(long) :: rc
  real :: r
  integer :: i, j, k, n, dum, ncl
  
  ! read rc as argument
  call getarg(1, argv)
  read(argv, *) rc
  print*, 'rc=', rc
  !read in native/nonnative counts
  open(16, file = "kkpkke_contacts.dat")
  do i = 1, nstr
    read(16,*) nat(i), nonnat(i)
  enddo
  close(16)
  ! read in rmsd
  cxn = .false.
  i = 1
  j = 2
  n = 0
  open(24, file='rmsd_pept_list.dat')
  do
    read(24, *, iostat=k) r
    if (k .ne. 0) exit
    if (r .lt. rc) then
      cxn(i, j) = .true.
      cxn(j, i) = .true.
    end if

    !rmsd(i, j) = r
    !rmsd(j, i) = r

    j = j + 1
    if (j .gt. nstr) then
      i = i + 1
      j = i + 1
    end if
    n = n + 1
  end do
  if (n .ne. nstr * (nstr - 1) / 2) then
    print*, 'n incorrect!', n
    print*, 'nstr is:', nstr
  end if

  !do i = 1, nstr - 1
  !  do j = i + 1, nstr 
  !    read(24, *) r
  !    if (r .lt. rc) then
  !      cxn(i, j) = .true.
  !      cxn(j, i) = .true.
  !    end if
  !    rmsd(i, j) = r
  !    rmsd(j, i) = r
  !  end do
  !end do

  close(24)

  ! Reading in bound status for combined trajectories
  used = 0
  !used(1) = 1
  !open(24, file = "bound_N3_combined.dat")
  !do i = 2, nstr, 1
  !  read(24, *) j
!	if (j .eq. 0) used(i) = 1
 ! end do
  ! perform clustering
  ncl = 1
  do
    cl(ncl)%list = 0
    maxcl = 0
    k = 0
    do i = 1, nstr, 1
      if (used(i) .ne. 0) cycle
      tempcl = 0
      tempcl(i) = 1
      do j = 1, nstr, 1
        if (i .eq. j) cycle
        if (used(j) .ne. 0) cycle
        if (cxn(i, j) .eqv. .false.) cycle
        tempcl(j) = 1
      end do
      if (sum(tempcl) .lt. sum(maxcl)) cycle
      maxcl = tempcl
      k = i  ! k is the centroid
    end do
    if (sum(maxcl) .lt. int(0.01 * nstr)) then
      ncl = ncl - 1  ! erase this cluster from existance 
      exit
    end if
    cl(ncl)%list = maxcl
    cl(ncl)%centroid = k
    cl(ncl)%sumsq = 0.
    cl(ncl)%sumnat = 0.
    cl(ncl)%sumnonnat = 0.
    do i = 1, nstr, 1
      if (maxcl(i) .ne. 1) cycle
      if (used(i) .ge. 1) then
        print*, 'Error: structure already used.'
        stop
      end if
      used(i) = ncl
      cl(ncl)%sumnat = cl(ncl)%sumnat + nat(i)
      cl(ncl)%sumnonnat = cl(ncl)%sumnonnat + nonnat(i)
      ! cl(ncl)%sumsq = cl(ncl)%sumsq + rmsd(i, k)**2.d0
    end do
    if (sum(used) .eq. nstr) exit
    if (ncl .eq. nclmax) exit
    ncl = ncl + 1
  end do

  ! output clusters
  do i = 1, ncl, 1
    if (sum(cl(i)%list) .eq. 1) cycle
    print*, cl(i)%centroid, 100.d0 * sum(cl(i)%list) / real(nstr, long), &
            cl(i)%sumnat / sum(cl(i)%list), cl(i)%sumnonnat / sum(cl(i)%list)
  end do
  ! output list of cluster ids per frame
  open(24, file = "output/cluster_ids.dat")
  do i = 1, nstr
    write(24,*) used(i)
  enddo
  close(24)
  ! output xyz file by cluster id
  ! readin xyz coords
  !open(25, file = "kkpkke_i1.xyz")
  !read(25,*) dum
  !read(25,*) dum
  !do i = 1, nstr
  !  read(25,*) dum, xyz(i,:)
  !enddo 
  !close(25)
  !write to file
  !open (26, file = "output/kkpkke_i1_cluster.xyz")
  !write(26,*) nstr
  !write(26,*) " "
  !do i = 1, nstr
   ! write(26,*) used(i), xyz(i,:)
   ! enddo

end program compute_cluster
