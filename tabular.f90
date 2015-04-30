      !'nnn' is the the maximum number of individuals
      !'nnnn' is the maximum number of gametes
      !'n' is the number of individuals on data file
      !'nn' is the number of gametes on data file
      implicit none
      integer n, nn, nnn, nnnn
      parameter (nnn = 50, nnnn = 2 * nnn)
      integer i, j, ii, jj, k, l 
      integer i1, i2, i3, i4
      character*1 ijk
      integer ia (0:nnn)  , ip (0:nnn)  , im (0:nnn)
      integer iag(0:nnnn) , ipg(0:nnnn) , img(0:nnnn)
      integer copies(0:nnnn,2), gcopies(0:nnnn), id, nd, dico
      real*8 mat(0:nnnn, 0:nnnn), mat3(0:nnnn, 0:nnnn, 0:nnnn), x1
      real*8 mat4 (0:nnnn, 0:nnnn, 0:nnnn, 0:nnnn)
      real*8 mat22(0:nnnn, 0:nnnn, 0:nnnn, 0:nnnn)
      real*8 matg(0:nnnn, 0:nnnn, 0:nnnn, 0:nnnn), delta(15)
      real t1, t2
      character arg*50
      call cpu_time(t1)
      x1 = 0.23847362d+10
      nd = 10**8
      do i = 0, nnn
        do j = 0, nnn
          mat (i,j) = 0.0d0
        enddo
      enddo
      call getarg(1, arg)
      open (5,  file=arg)
      k = 0
      read(5,*)n
      nn = 2 * n
      do i = 1, n
        read(5,*) ia(i), ip(i), im(i)
        k = k + 1
        iag (k) = k
        ipg (k) = 2*ip(i) - 1
        img (k) = 2*ip(i)
        k = k + 1
        iag (k) = k
        ipg (k) = 2*im(i) - 1
        img (k) = 2*im(i) 
      enddo
      do i=1,nn
        if (ipg(i) .lt. 0) ipg(i) = 0
        if (img(i) .lt. 0) img(i) = 0
      enddo
      close(5)
      ! tabular conventional algorithm for gametes
      do i = 1, nn
        mat(i,i) = 1.0d0 
        do j = i+1, nn
          mat(i,j) = 0.5d0 * (mat(i,ipg(j)) + mat(i,img(j)))
          mat(j,i) = mat(i,j)
        enddo
      enddo

      !3D tabular for gametes
      do i = 1, nn
        do j = 1, nn
          do k = 1, nn
              mat3(i,j,k) = 0.0d0
          enddo
        enddo
      enddo
      do i = 1, nn
        do j = i, nn
          do k = j, nn
            if (i.eq.j .and. j.eq.k) then
              !3d corner
              mat3(i,i,i) = 1.0d0 
            else if (i.ne.j .and. j.eq.k) then
              !edge
              mat3(i,j,j) = mat(i,j)
              mat3(j,i,j) = mat3(i,j,j)
              mat3(j,j,i) = mat3(i,j,j)
            else
              !2d column
              mat3(i,j,k) = 0.5d0 * (mat3(i,j,ipg(k))+mat3(i,j,img(k)))
              mat3(i,k,j) = mat3(i,j,k)
              mat3(j,i,k) = mat3(i,j,k)
              mat3(j,k,i) = mat3(i,j,k)
              mat3(k,i,j) = mat3(i,j,k)
              mat3(k,j,i) = mat3(i,j,k)
            end if
          enddo
        enddo
      enddo

      !4D tabular for gametes
      do i = 0, nn
        do j = 0, nn
          do k = 0, nn
            do l = 0, nn
              mat4(i,j,k,l) = 0.0d0
            enddo
          enddo
        enddo
      enddo
      do i = 1, nn
        do j = i, nn
          do k = j, nn
            do l = k, nn
              if (i.eq.j .and. j.eq.k .and. k.eq.l) then
                !4d corner
                mat4(i,i,i,i) = 1.0d0 
              else if (i.ne.j .and. j.eq.k .and. k.eq.l) then
                !3d edge
                mat4(i,j,j,j) = mat(i,j)
                mat4(j,i,j,j) = mat4(i,j,j,j)
                mat4(j,j,i,j) = mat4(i,j,j,j)
                mat4(j,j,j,i) = mat4(i,j,j,j)
              else if (j.ne.k .and. k.eq.l) then
                !2d edge (plane)
                mat4(i,j,k,k) = mat3(i,j,k)
                mat4(j,i,k,k) = mat4(i,j,k,k)
                mat4(k,k,i,j) = mat4(i,j,k,k)
                mat4(k,k,j,i) = mat4(i,j,k,k)
                mat4(k,i,k,j) = mat4(i,j,k,k)
                mat4(k,j,k,i) = mat4(i,j,k,k)
                mat4(k,i,j,k) = mat4(i,j,k,k)
                mat4(k,j,i,k) = mat4(i,j,k,k)
                mat4(i,k,k,j) = mat4(i,j,k,k)
                mat4(j,k,k,i) = mat4(i,j,k,k)
                mat4(i,k,j,k) = mat4(i,j,k,k)
                mat4(j,k,i,k) = mat4(i,j,k,k)
              else
                mat4(i,j,k,l) = 0.5d0*(mat4(i,j,k,ipg(l))+mat4(i,j,k,img(l)))
                mat4(i,j,l,k) = mat4(i,j,k,l)
                mat4(i,l,j,k) = mat4(i,j,k,l)
                mat4(i,k,j,l) = mat4(i,j,k,l)
                mat4(i,k,l,j) = mat4(i,j,k,l)
                mat4(i,l,k,j) = mat4(i,j,k,l)
                mat4(j,i,k,l) = mat4(i,j,k,l)
                mat4(j,i,l,k) = mat4(i,j,k,l)
                mat4(k,i,j,l) = mat4(i,j,k,l)
                mat4(l,i,j,k) = mat4(i,j,k,l)
                mat4(k,i,l,j) = mat4(i,j,k,l)
                mat4(l,i,k,j) = mat4(i,j,k,l)
                mat4(j,k,i,l) = mat4(i,j,k,l)
                mat4(j,l,i,k) = mat4(i,j,k,l)
                mat4(k,j,i,l) = mat4(i,j,k,l)
                mat4(l,j,i,k) = mat4(i,j,k,l)
                mat4(k,l,i,j) = mat4(i,j,k,l)
                mat4(l,k,i,j) = mat4(i,j,k,l)
                mat4(j,k,l,i) = mat4(i,j,k,l)
                mat4(j,l,k,i) = mat4(i,j,k,l)
                mat4(k,j,l,i) = mat4(i,j,k,l)
                mat4(l,j,k,i) = mat4(i,j,k,l)
                mat4(k,l,j,i) = mat4(i,j,k,l)
                mat4(l,k,j,i) = mat4(i,j,k,l)
              end if
            enddo
          enddo
        enddo
      enddo
      
      ! 2-2D tabular for gametes
      do i = 0, nn
        do j = 0, nn
          do k = 0, nn
            do l = 0, nn
              mat22(i,j,k,l) = 0.0d0
            enddo
          enddo
        enddo
      enddo
      do i = 1, nn
      do k = 1, nn     
      do j = i, nn
      do l = k, nn
        if (i.eq.j .and. j.eq.k .and. k.eq.l) then
		  ijk = 'E'
        else if (i.lt.j .and. j.eq.k .and. k.eq.l) then
		  ijk = 'D'
        else if (k.lt.l .and. l.eq.i .and. i.eq.j) then
		  ijk = 'D'
        else if (j.eq.l .and. k.ne.l .and. i.ne.l) then
		  ijk = 'C'
        else if (i.lt.k .and. j.lt.k .and. k.eq.l) then
		  ijk = 'B'
        else if (l.lt.j .and. k.lt.j .and. i.eq.j) then
		  ijk = 'B'
		else if (l .gt. j) then
		  ijk = 'A'
		else 
		  ijk = 'F'
		endif

        if (ijk .eq. 'E') then
          !4D corner
          mat22(i,i,i,i) = 1.0d0 
	      !print *,'-------   ',i,j,k,l,'E'
        else if (ijk .eq. 'D') then
          !3D edge
          mat22(i,j,j,j) = mat(i,j)
          mat22(j,i,j,j) = mat22(i,j,j,j)
          mat22(j,j,i,j) = mat22(i,j,j,j)
          mat22(j,j,j,i) = mat22(i,j,j,j)
	      else if (ijk .eq. 'B') then
          !2D edge (plane)
          mat22(i,j,k,k) = mat(i,j)
          mat22(j,i,k,k) = mat22(i,j,k,k)
          mat22(k,k,i,j) = mat22(i,j,k,k)
          mat22(k,k,j,i) = mat22(i,j,k,k)
		else if (ijk .eq. 'C') then 
		  mat22(i,j,k,j) = mat3(i,j,k)
          mat22(k,j,i,j) = mat22(i,j,k,j)
          mat22(k,j,j,i) = mat22(i,j,k,j)
          mat22(i,j,j,k) = mat22(i,j,k,j)
          mat22(j,i,j,k) = mat22(i,j,k,j)
          mat22(j,k,j,i) = mat22(i,j,k,j)
          mat22(j,k,i,j) = mat22(i,j,k,j)
          mat22(j,i,k,j) = mat22(i,j,k,j)
		else
          if (ijk .eq. 'A') then
            mat22(i,j,k,l)=0.5d0*(mat22(i,j,k,ipg(l))+mat22(i,j,k,img(l)))
          else
            mat22(i,j,k,l)=0.5d0*(mat22(i,ipg(j),k,l)+mat22(i,img(j),k,l))
          end if
          mat22(i,j,l,k) = mat22(i,j,k,l)
          mat22(j,i,k,l) = mat22(i,j,k,l)
          mat22(j,i,l,k) = mat22(i,j,k,l)
          mat22(k,l,i,j) = mat22(i,j,k,l)
          mat22(l,k,i,j) = mat22(i,j,k,l)
          mat22(k,l,j,i) = mat22(i,j,k,l)
          mat22(l,k,j,i) = mat22(i,j,k,l)
		end if

      enddo
      enddo
      enddo
      enddo

      call getarg(2,arg)
      read(arg,*)i
      call getarg(3,arg)
      read(arg,*)j
      i1 = i * 2 - 1
      i2 = i * 2
      i3 = j * 2 - 1
      i4 = j * 2
      delta(1)  = mat4 (i1,i2,i3,i4)
      delta(2)  = mat4 (i1,i2,i3,i3) - delta(1)
      delta(3)  = mat4 (i1,i2,i2,i4) - delta(1)
      delta(4)  = mat4 (i1,i1,i3,i4) - delta(1)
      delta(5)  = mat4 (i2,i2,i3,i4) - delta(1)
      delta(6)  = mat22(i1,i2,i3,i4) - delta(1)
      delta(7)  = mat22(i1,i2,i3,i3) - delta(6) - delta(1) - delta(2) - delta(3)
      delta(8)  = mat22(i1,i1,i3,i4) - delta(6) - delta(1) - delta(4) - delta(5)
      delta(9)  = mat22(i1,i3,i2,i4) - delta(1)
      delta(10) = mat22(i1,i3,i2,i2) - delta(9) - delta(1) - delta(2) - delta(4)
      delta(11) = mat22(i1,i1,i2,i4) - delta(9) - delta(1) - delta(3) - delta(5)
      delta(12) = mat22(i1,i4,i2,i3) - delta(1)
      delta(13) = mat22(i1,i4,i2,i2) - delta(12) - delta(1) - delta(3) - delta(4)
      delta(14) = mat22(i1,i1,i2,i3) - delta(12) - delta(1) - delta(2) - delta(5)
      delta(15) = mat22(i1,i1,i1,i1) - delta(1) - delta(2) - delta(3) - delta(4) &
                                   & - delta(5) - delta(6) - delta(7) - delta(8) &
                                   & - delta(9) - delta(10) - delta(11) - delta(12) &
                                   & - delta(13) - delta(14)
       
      print '(11hdelta 1  = ,f7.5, 3x, f7.5)', delta(1)
      print '(11hdelta 2  = ,f7.5, 3x, f7.5)', delta(2)
      print '(11hdelta 3  = ,f7.5, 3x, f7.5)', delta(3)
      print '(11hdelta 4  = ,f7.5, 3x, f7.5)', delta(4)
      print '(11hdelta 5  = ,f7.5, 3x, f7.5)', delta(5)
      print '(11hdelta 6  = ,f7.5, 3x, f7.5)', delta(6)
      print '(11hdelta 7  = ,f7.5, 3x, f7.5)', delta(7)
      print '(11hdelta 8  = ,f7.5, 3x, f7.5)', delta(8)
      print '(11hdelta 9  = ,f7.5, 3x, f7.5)', delta(9)
      print '(11hdelta 10 = ,f7.5, 3x, f7.5)', delta(10)
      print '(11hdelta 11 = ,f7.5, 3x, f7.5)', delta(11)
      print '(11hdelta 12 = ,f7.5, 3x, f7.5)', delta(12)
      print '(11hdelta 13 = ,f7.5, 3x, f7.5)', delta(13)
      print '(11hdelta 14 = ,f7.5, 3x, f7.5)', delta(14)
      print '(11hdelta 15 = ,f7.5, 3x, f7.5)', delta(15)

      stop
      end

