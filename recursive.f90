implicit none
integer :: animal  (0:100), padre  (0:100), madre  (0:100)
integer :: ganimal (0:200), gpadre (0:200), gmadre (0:200)
integer :: n, nn, i, j, k, llamadas, est(20), s(20), igam(4), m(20)
real*8  :: par, sum
integer :: prof  !recursion depth
integer :: hay, k1, i1, i2
logical :: checklista
integer :: caso(15, 4)
real*8  :: phi(15), delta(15)  ! phi are multiple relationships
                               ! delta are identity coefficients
real*8  :: suma, suma2, temp
integer :: cuantos(15)
character :: arg*50

! calculate the 15 partitions and store in "caso"
print *, 'Calculating partitions of four elements'
do i = 1, 4
  s(i) = 1
  m(i) = 1
enddo
do i=1,4
  caso(1,i) = s(i)
enddo
cuantos(1) = 1
write(*,'(i10,a,i2,a,17i4)')1,' -- ',1,' -- ',(caso(1,j),j=1,4)
do k = 1, 15
  phi(k)   = 0.0d0
  delta(k) = 0.0d0
enddo
do k = 2, 15
  do i = 1, 4
    est(i) = s(i)
  enddo 
  call next(s, m, 4, hay)
  if (hay.eq.0) then
    print *,'>>> ',k  !combination not found
  else
    do i=1,4
      caso(k,i) = s(i)
    enddo
    cuantos(k) = 0
    do i = 1, 4
      do j = 1, 4
        if (caso(k,j).eq.i)then
          cuantos(k) = cuantos(k) + 1
          exit
        end if
      end do
    enddo
    write(*,'(i10,a,i2,a,4i4)')k,' -- ',cuantos(k),' -- ', (caso(k,j),j=1,4)
  end if
enddo

! padre and madre means sire and dam (individuals)
! gpadre and gmadre the same for gametes (double size)

animal  (0) = 0
padre   (0) = 0
madre   (0) = 0
ganimal (0) = 0
gpadre  (0) = 0
gmadre  (0) = 0
call getarg(1, arg)
open (5, file=arg)
call getarg(2, arg)
read(arg,*)i1 
call getarg(3, arg)
read(arg,*)i2 
k = 0
read(5,*) n
do i = 1, n
  read(5,*)animal(i), padre(i), madre(i)
  k = k + 1
  gpadre  (k) = 2*padre(i)-1
  if (gpadre(k).lt.0)gpadre(k)=0
  gmadre  (k) = 2*padre(i)
  if (gmadre(k).lt.0)gmadre(k)=0
  k = k + 1
  ganimal (k) = k
  gpadre  (k) = 2*madre(i)-1
  if (gpadre(k).lt.0)gpadre(k)=0
  gmadre  (k) = 2*madre(i)
  if (gmadre(k).lt.0)gmadre(k)=0
end do
 
do i = 1, 4
  s(i) = 1
  m(i) = 1
enddo

print *, 'Calculating multiple relationships'
print *,'Partition - number of recursive calls - multiple relationship - partition'
do k = 1, 15 
  llamadas = 0
  do i = 1, 4
    est(i) = caso(k,i)
  enddo
  !the four gametes
  igam(1) = 2*i1-1
  igam(2) = 2*i1
  igam(3) = 2*i2-1
  igam(4) = 2*i2
  nn = 4
  call calc_iden(nn,est,igam,par,gpadre,gmadre,llamadas,1)
  write(*,'(2i10,a,f20.16,17i4)')k,llamadas,' -- ', par,(caso(k,j),j=1,4)
   if (k.le. 10) then
  write(72,'(i10,a,f20.16,17i4)')llamadas,' -- ', par,(caso(k,j),j=1,4)
   end if
   if (k.ge. 991 .and. k .le. 1000) then
  write(73,'(i10,a,f20.16,17i4)')llamadas,' -- ', par,(caso(k,j),j=1,4)
   end if
   if (k.ge. 4131) then
  write(74,'(i10,a,f20.16,17i4)')llamadas,' -- ', par,(caso(k,j),j=1,4)
   end if
  delta(k) = par - phi(k)
  phi(k) = par

  do j = k+1, 15
    if (checklista(4,caso(j,:),caso(k,:))) then
      phi(j) = phi(j) + delta(k)
    end if
  enddo
enddo

print '(/,a,/)','delta - multiple relationship - partition'

do i = 1, 15
  print '(3h== ,i5,2f12.5,20i4)',i,delta(i),phi(i),(caso(i,j),j=1,4)
  write(55,'(i5,2f12.6,20i4)')i,delta(i),phi(i),(caso(i,j),j=1,4)
enddo

stop
end

!-------------------------------------------------------------------------------
function unif (x1)
!generacion de un numero uniforme u[0, 1] x1 es la semilla
real*8 x1, unif
x1 = dmod(16807.0d0 * x1, 0.2147483647d+10)
unif = x1 / 0.2147483647d+10
return
end

function checklista(n,a,b)
!busca si la lista a es caso particular de la b
!
! por ejemplo (112232332) es caso particular de (1122222222)
! porque todo lo que es equivalente en la primera es equivalente en la segunda
integer :: n,a(20),b(20)
logical :: checklista

do i = 1, n
  do j = i+1, n
    if (a(i).eq.a(j)) then
      if (b(i) .ne. b(j)) then
        checklista = .false.
        return
      end if  
    end if
  enddo
enddo
checklista = .true.
return
end

recursive subroutine calc_iden(nn,est,igam,par,gpadre,gmadre,llamadas,prof)
implicit none
integer nn
integer :: gpadre (0:200), gmadre (0:200)
integer :: nn1,nn2, nn_viejo, i, j, k, k1, k2, k3, llamadas, est(20), igam(4)
real*8  :: par, par1, par2
integer :: igam1(20), igam2(20), est1(20), est2(20),im(20)
integer :: prof, prof1, prof2

llamadas = llamadas + 1
nn_viejo = 1000
do while (nn_viejo .ne. nn)
  nn_viejo = nn
  do i = 1, nn
    do j = i+1, nn
      if (igam(i).eq.igam(j)) then
        k1 = est(i)
        k2 = est(j)
        k3 = j
        do k = 1, nn
          if (est(k) .eq. k2) est(k) = k1
        enddo
        do k = k3+1, nn
          igam(k-1) = igam(k)
          est(k-1) = est(k)
        enddo
        nn = nn - 1
        goto 7
      end if
    enddo
  enddo 
7 continue
  ! buscar el primer grupo con un solo gameto y eliminarlo  
  if (nn.ne.1) then 
    do i = 1, nn
      k = est(i)
      k1 = 0
      do j = 1, nn
        if (est(j).eq.k) k1 = k1 + 1
      enddo
      if (k1 .eq. 1) then
        !quitar el bicho "i"
        do k3 = i+1, nn
          igam(k3-1) = igam(k3)
          est (k3-1) = est (k3)
        enddo
        nn = nn - 1
        goto 6
      end if
    enddo
  end if
6 continue  
enddo
if (nn.eq.1) then
  par = 1.0d0
  return
end if

i=-9999
do j=1,nn
  if (igam(j).gt.i)then
    i = igam(j) 
    k = j       
  end if
enddo

do i = 1, nn
  est1(i) = est(i)
  est2(i) = est(i)
  igam1(i) = igam(i)
  igam2(i) = igam(i)
enddo
nn1 = nn
nn2 = nn
prof1 = prof + 1
prof2 = prof + 1
if (gpadre(igam(k)).gt.0) then
  igam1(k) = gpadre(igam(k))
  call calc_iden(nn1,est1,igam1,par1,gpadre,gmadre,llamadas,prof1)
else
  par1=0.0d0
end if
if (gmadre(igam(k)).gt.0) then
  igam2(k) = gmadre(igam(k))
  call calc_iden(nn2,est2,igam2,par2,gpadre,gmadre,llamadas,prof2)
else
  par2=0.0d0
end if
par = 0.5d0*(par1+par2)
return
end

subroutine next(s,m,n,hay)
! look for next partition
integer s(20), m(20), n, hay
integer i
i = 0
s(i+1) = s(i+1) + 1
do while (i < n-1 .and. s(i+1) > m(i+2)+1)
  s(i+1) = 1
  i = i + 1
  s(i+1) = s(i+1) + 1
enddo
if (i .eq. n - 1) then
   hay = 0
   return
end if
if ( s(i+1) .gt. m(i+1) ) then
    m(i+1) = s(i+1)
end if
do j=i-1, 0, -1
  m(j+1) = m(i+1)
enddo
hay = 1
! hay = 0 --> partition not found
! hay = 1 --> partition found
return
end


