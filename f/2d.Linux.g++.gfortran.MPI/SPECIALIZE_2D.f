

c
c     This function copies floating point numbers from 1.0D0 array to another.
c
      subroutine fastcopy (dest,dest_l1, dest_l2, dest_h1, dest_h2,imin,
     & jmin, imax, jmax,src,src_l1, src_l2, src_h1, src_h2,imn, jmn,nc
     &omp)

      implicit none
c
c     Bounds to fill in dest
c
      integer imin, jmin, imax, jmax
      integer dest_l1, dest_l2, dest_h1, dest_h2
c
c     Bounds to fill from src
c
      integer imn,  jmn
      integer src_l1, src_l2, src_h1, src_h2
      integer ncomp

      DOUBLE PRECISION  dest(dest_l1:dest_h1, dest_l2:dest_h2,ncomp)
      DOUBLE PRECISION  src(src_l1:src_h1, src_l2:src_h2,ncomp)
c
c     Local variables
c
      integer i,j,k,ioff,joff

      ioff=imn-imin
      joff=jmn-jmin

      do k = 1, ncomp
         do j = jmin,jmax
            do i = imin,imax
               dest(i,j,k) = src(i+ioff,j+joff,k)
            end do
         end do
      end do

      end
c
c     This function copies from a 2D array to a 1D 1.0D0.
c
      subroutine fastcopytomem (lo,hi,data,data_l1, data_l2, data_h1, da
     &ta_h2,ncomp,dst)

      implicit none

      integer ncomp
      integer lo(2), hi(2)
      integer data_l1, data_l2, data_h1, data_h2
      DOUBLE PRECISION data(data_l1:data_h1, data_l2:data_h2,ncomp), dst
     &(*)
c
c     Local variables.
c
      integer i,j,k,cnt

      cnt = 1

      do k = 1, ncomp
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               dst(cnt) = data(i,j,k)
               cnt = cnt + 1
            end do
         end do
      end do

      end
c
c     This function copies from a 1D array into a 2D array.
c
      subroutine fastcopyfrommem (lo,hi,data,data_l1, data_l2, data_h1, 
     &data_h2,ncomp,src)

      implicit none

      integer ncomp
      integer lo(2), hi(2)
      integer data_l1, data_l2, data_h1, data_h2
      DOUBLE PRECISION data(data_l1:data_h1, data_l2:data_h2,ncomp), src
     &(*)
c
c     Local variables.
c
      integer i,j,k,cnt

      cnt = 1

      do k = 1, ncomp
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               data(i,j,k) = src(cnt)
               cnt = cnt + 1
            end do
         end do
      end do

      end
c
c     This function sets a section of an array to a value.
c
      subroutine fastsetval (val,lo,hi,dest,dest_l1, dest_l2, dest_h1, d
     &est_h2,ncomp)

      implicit none

      integer ncomp
      integer lo(2), hi(2)
      integer dest_l1, dest_l2, dest_h1, dest_h2
      DOUBLE PRECISION  val
      DOUBLE PRECISION  dest(dest_l1:dest_h1, dest_l2:dest_h2,ncomp)
c
c     Local variables
c 
      integer i,j,k

      do k = 1, ncomp
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               dest(i,j,k) = val
            end do
         end do
      end do

      end
c
c     Calculate max(abs) over section of an array.
c
      subroutine fastzeronorm (src,src_l1, src_l2, src_h1, src_h2,lo,hi,
     &ncomp,nrm)

      implicit none

      integer ncomp
      integer lo(3), hi(3)
      integer src_l1, src_l2, src_h1, src_h2
      DOUBLE PRECISION  src(src_l1:src_h1, src_l2:src_h2,ncomp), nrm
c
c     Local variables.
c
      integer i,j,k

      nrm = 0.0d0

      do k = 1,ncomp
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               nrm = max(nrm, abs(src(i,j,k)))
            end do
         end do
      end do

      end
c
c     Calculate sum(abs) over section of an array.
c
      subroutine fastonenorm (src,src_l1, src_l2, src_h1, src_h2,lo,hi,n
     &comp,nrm)

      implicit none

      integer ncomp
      integer lo(3), hi(3)
      integer src_l1, src_l2, src_h1, src_h2
      DOUBLE PRECISION  src(src_l1:src_h1, src_l2:src_h2,ncomp), nrm
c
c     Local variables.
c
      integer i,j,k

      nrm = 0.0d0

      do k = 1, ncomp
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               nrm = nrm + abs(src(i,j,k))
            end do
         end do
      end do

      end
c
c     Calculate sum() over section of an array.
c
      subroutine fastsum (src,src_l1, src_l2, src_h1, src_h2,lo,hi,ncomp
     &,sm)

      implicit none

      integer ncomp
      integer lo(3), hi(3)
      integer src_l1, src_l2, src_h1, src_h2
      DOUBLE PRECISION  src(src_l1:src_h1, src_l2:src_h2,ncomp), sm
c
c     Local variables.
c
      integer i,j,k

      sm = 0.0d0

      do k = 1, ncomp
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               sm = sm + src(i,j,k)
            end do
         end do
      end do

      end
c
c     This function adds floating point numbers from 1.0D0 array to another.
c
      subroutine fastplus (dest,dest_l1, dest_l2, dest_h1, dest_h2,imin,
     & jmin, imax, jmax,src,src_l1, src_l2, src_h1, src_h2,imn, jmn,nc
     &omp)

      implicit none
c
c     Bounds to fill in dest
c
      integer imin, jmin, imax, jmax
      integer dest_l1, dest_l2, dest_h1, dest_h2
c
c     Bounds to fill from src
c
      integer imn,  jmn
      integer src_l1, src_l2, src_h1, src_h2
      integer ncomp

      DOUBLE PRECISION  dest(dest_l1:dest_h1, dest_l2:dest_h2,ncomp)
      DOUBLE PRECISION  src(src_l1:src_h1, src_l2:src_h2,ncomp)
c
c     Local variables
c
      integer i,j,k,ioff,joff

      ioff=imn-imin
      joff=jmn-jmin

      do k = 1, ncomp
         do j = jmin,jmax
            do i = imin,imax
               dest(i,j,k) = dest(i,j,k) + src(i+ioff,j+joff,k)
            end do
         end do
      end do

      end
c
c     This function multiplys floating point numbers from 1.0D0 array to another.
c
      subroutine fastmult (dest,dest_l1, dest_l2, dest_h1, dest_h2,imin,
     & jmin, imax, jmax,src,src_l1, src_l2, src_h1, src_h2,imn, jmn,nc
     &omp)

      implicit none
c
c     Bounds to fill in dest
c
      integer imin, jmin, imax, jmax
      integer dest_l1, dest_l2, dest_h1, dest_h2
c
c     Bounds to fill from src
c
      integer imn,  jmn
      integer src_l1, src_l2, src_h1, src_h2
      integer ncomp

      DOUBLE PRECISION  dest(dest_l1:dest_h1, dest_l2:dest_h2,ncomp)
      DOUBLE PRECISION  src(src_l1:src_h1, src_l2:src_h2,ncomp)
c
c     Local variables
c
      integer i,j,k,ioff,joff

      ioff=imn-imin
      joff=jmn-jmin

      do k = 1, ncomp
         do j = jmin,jmax
            do i = imin,imax
               dest(i,j,k) = dest(i,j,k) * src(i+ioff,j+joff,k)
            end do
         end do
      end do

      end
c
c     This function adds scaled floating point numbers from 1.0D0 array to another.
c
      subroutine fastsaxpy (dest,dest_l1, dest_l2, dest_h1, dest_h2,imin
     &, jmin, imax, jmax,a, src,src_l1, src_l2, src_h1, src_h2,imn, jm
     &n,ncomp)

      implicit none
c
c     Bounds to fill in dest
c
      integer imin, jmin, imax, jmax
      integer dest_l1, dest_l2, dest_h1, dest_h2
c
c     Bounds to fill from src
c
      integer imn,  jmn
      integer src_l1, src_l2, src_h1, src_h2
      integer ncomp

      DOUBLE PRECISION  dest(dest_l1:dest_h1, dest_l2:dest_h2,ncomp)
      DOUBLE PRECISION  a
      DOUBLE PRECISION  src(src_l1:src_h1, src_l2:src_h2,ncomp)
c
c     Local variables
c
      integer i,j,k,ioff,joff

      ioff=imn-imin
      joff=jmn-jmin

      do k = 1, ncomp
         do j = jmin, jmax
            do i = imin, imax
               dest(i,j,k) = dest(i,j,k) + a * src(i+ioff,j+joff,k)
            end do
         end do
      end do

      end
c
c     This function subtracts floating point numbers from 1.0D0 array to another.
c
      subroutine fastminus (dest,dest_l1, dest_l2, dest_h1, dest_h2,imin
     &, jmin, imax, jmax,src,src_l1, src_l2, src_h1, src_h2,imn, jmn,n
     &comp)

      implicit none
c
c     Bounds to fill in dest
c
      integer imin, jmin, imax, jmax
      integer dest_l1, dest_l2, dest_h1, dest_h2
c
c     Bounds to fill from src
c
      integer imn,  jmn
      integer src_l1, src_l2, src_h1, src_h2
      integer ncomp

      DOUBLE PRECISION  dest(dest_l1:dest_h1, dest_l2:dest_h2,ncomp)
      DOUBLE PRECISION  src(src_l1:src_h1, src_l2:src_h2,ncomp)
c
c     Local variables
c
      integer i,j,k,ioff,joff

      ioff=imn-imin
      joff=jmn-jmin

      do k = 1, ncomp
         do j = jmin,jmax
            do i = imin,imax
               dest(i,j,k) = dest(i,j,k) - src(i+ioff,j+joff,k)
            end do
         end do
      end do

      end
c
c     This function divides floating point numbers from 1.0D0 array by another
c
      subroutine fastdivide (dest,dest_l1, dest_l2, dest_h1, dest_h2,imi
     &n, jmin, imax, jmax,src,src_l1, src_l2, src_h1, src_h2,imn, jmn,
     &ncomp)

      implicit none
c
c     Bounds to fill in dest
c
      integer imin, jmin, imax, jmax
      integer dest_l1, dest_l2, dest_h1, dest_h2
c
c     Bounds to fill from src
c
      integer imn,  jmn
      integer src_l1, src_l2, src_h1, src_h2
      integer ncomp

      DOUBLE PRECISION  dest(dest_l1:dest_h1, dest_l2:dest_h2,ncomp)
      DOUBLE PRECISION  src(src_l1:src_h1, src_l2:src_h2,ncomp)
c
c     Local variables
c
      integer i,j,k,ioff,joff

      ioff=imn-imin
      joff=jmn-jmin

      do k = 1, ncomp
         do j = jmin,jmax
            do i = imin,imax
               dest(i,j,k) = dest(i,j,k) / src(i+ioff,j+joff,k)
            end do
         end do
      end do

      end
c
c     This function divides floating point numbers from 1.0D0 array by another
c     wherever the denominator is non-0.0D0.
c
      subroutine fastprotdivide (dest,dest_l1, dest_l2, dest_h1, dest_h2
     &,imin, jmin, imax, jmax,src,src_l1, src_l2, src_h1, src_h2,imn, 
     &jmn,ncomp)

      implicit none
c
c     Bounds to fill in dest
c
      integer imin, jmin, imax, jmax
      integer dest_l1, dest_l2, dest_h1, dest_h2
c
c     Bounds to fill from src
c
      integer imn,  jmn
      integer src_l1, src_l2, src_h1, src_h2
      integer ncomp

      DOUBLE PRECISION  dest(dest_l1:dest_h1, dest_l2:dest_h2,ncomp)
      DOUBLE PRECISION  src(src_l1:src_h1, src_l2:src_h2,ncomp)
c
c     Local variables
c
      integer i,j,k,ioff,joff

      ioff=imn-imin
      joff=jmn-jmin

      do k = 1, ncomp
         do j = jmin,jmax
            do i = imin,imax
               if (src(i+ioff,j+joff,k) .ne. 0.0D0) then
                  dest(i,j,k) = dest(i,j,k) / src(i+ioff,j+joff,k)
               end if
            end do
         end do
      end do

      end

      subroutine fastinvert (src,src_l1, src_l2, src_h1, src_h2,lo,hi,va
     &l,ncomp)

      implicit none

      integer ncomp
      integer lo(2), hi(2)
      integer src_l1, src_l2, src_h1, src_h2
      DOUBLE PRECISION  src(src_l1:src_h1, src_l2:src_h2,ncomp), val
c
c     Local variables.
c
      integer i,j,k

      do k = 1,ncomp
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               src(i,j,k) = val / src(i,j,k)
            end do
         end do
      end do

      end

