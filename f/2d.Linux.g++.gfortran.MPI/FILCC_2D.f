

c ::: -----------------------------------------------------------
c ::: This routine is intended to be a generic fill function
c ::: for cell-centered data.  It knows how to extrapolate
c ::: and reflect data and is used to supplement the problem-specific
c ::: fill functions which call it.
c ::: 
c ::: INPUTS/OUTPUTS:
c ::: q           <=  array to fill
c ::: lo,hi        => index extent of q array
c ::: domlo,domhi  => index extent of problem domain
c ::: dx           => cell spacing
c ::: xlo          => physical location of lower left hand
c :::	              corner of q array
c ::: bc	   => array of boundary flags bc(SPACEDIM,lo:hi)
c ::: 
c ::: NOTE: all corner as well as edge data is filled if not 3
c ::: -----------------------------------------------------------
      subroutine filcc(q,q_l1, q_l2, q_h1, q_h2,domlo,domhi,dx,xlo,bc)

      implicit none

      integer    q_l1, q_l2, q_h1, q_h2
      integer    domlo(2), domhi(2)
      integer    bc(2,2)
      DOUBLE PRECISION     xlo(2), dx(2)
      DOUBLE PRECISION     q(q_l1:q_h1, q_l2:q_h2)

      integer    nlft, nrgt, nbot, ntop
      integer    ilo, ihi, jlo, jhi
      integer    i, j
      integer    is, ie, js, je

      nlft = max(0,domlo(1)-q_l1)
      nrgt = max(0,q_h1-domhi(1))
      nbot = max(0,domlo(2)-q_l2)
      ntop = max(0,q_h2-domhi(2))

      is = max(q_l1,domlo(1))
      ie = min(q_h1,domhi(1))
      js = max(q_l2,domlo(2))
      je = min(q_h2,domhi(2))

c     ::::: first fill sides
      if (nlft .gt. 0) then
         ilo = domlo(1)

         if (bc(1,1) .eq. 2) then
            do i = 1, nlft
            do j = q_l2, q_h2
               q(ilo-i,j) = q(ilo,j)
            end do
            end do
         else if (bc(1,1) .eq. 4) then
            do i = 2, nlft
            do j = q_l2, q_h2
               q(ilo-i,j) = q(ilo,j) 
            end do 
            end do 
            if (ilo+2 .le. ie) then 
             do j = q_l2, q_h2
                q(ilo-1,j) = (15*q(ilo,j) - 10*q(ilo+1,j) + 3*q(ilo+2,j)
     &) * 0.125D0
             end do 
            else 
             do j = q_l2, q_h2
               q(ilo-1,j) = 0.5D0*(3*q(ilo,j) - q(ilo+1,j))
             end do
            end if
         else if (bc(1,1) .eq. 1) then
            do i = 1, nlft
             do j = q_l2, q_h2
               q(ilo-i,j) = q(ilo+i-1,j)
            end do
            end do
         else if (bc(1,1) .eq. -1) then
            do i = 1, nlft
            do j = q_l2, q_h2
               q(ilo-i,j) = -q(ilo+i-1,j)
            end do
            end do
         end if
      end if

      if (nrgt .gt. 0) then
         ihi = domhi(1)

         if (bc(1,2) .eq. 2) then
            do i = 1, nrgt
            do j = q_l2, q_h2
               q(ihi+i,j) = q(ihi,j)
            end do
            end do
         else if (bc(1,2) .eq. 4) then
            do i = 2, nrgt
            do j = q_l2, q_h2
               q(ihi+i,j) = q(ihi,j)
            end do
            end do
            if (ihi-2 .ge. is) then
             do j = q_l2, q_h2
               q(ihi+1,j) = (15*q(ihi,j) - 10*q(ihi-1,j) + 3*q(ihi-2,j))
     & * 0.125D0
             end do
            else
             do j = q_l2, q_h2
               q(ihi+1,j) = 0.5D0*(3*q(ihi,j) - q(ihi-1,j))
             end do
            end if
         else if (bc(1,2) .eq. 1) then
            do i = 1, nrgt
            do j = q_l2, q_h2
               q(ihi+i,j) = q(ihi-i+1,j)
            end do
            end do
         else if (bc(1,2) .eq. -1) then
            do i = 1, nrgt
            do j = q_l2, q_h2
               q(ihi+i,j) = -q(ihi-i+1,j)
            end do
            end do
         end if
      end if

      if (nbot .gt. 0) then
         jlo = domlo(2)

         if (bc(2,1) .eq. 2) then
            do j = 1, nbot
            do i = q_l1, q_h1
               q(i,jlo-j) = q(i,jlo)
            end do
            end do
         else if (bc(2,1) .eq. 4) then
            do j = 2, nbot
            do i = q_l1, q_h1
               q(i,jlo-j) = q(i,jlo)
            end do
            end do
            if (jlo+2 .le. je) then
             do i = q_l1, q_h1
               q(i,jlo-1) = (15*q(i,jlo) - 10*q(i,jlo+1) + 3*q(i,jlo+2))
     & * 0.125D0
             end do
            else
             do i = q_l1, q_h1
               q(i,jlo-1) = 0.5D0*(3*q(i,jlo) - q(i,jlo+1))
             end do
            end if
         else if (bc(2,1) .eq. 1) then
            do j = 1, nbot
            do i = q_l1, q_h1
               q(i,jlo-j) = q(i,jlo+j-1)
            end do
            end do
         else if (bc(2,1) .eq. -1) then
            do j = 1, nbot
            do i = q_l1, q_h1
               q(i,jlo-j) = -q(i,jlo+j-1)
            end do
            end do
         end if
      end if

      if (ntop .gt. 0) then
         jhi = domhi(2)

         if (bc(2,2) .eq. 2) then
            do j = 1, ntop
            do i = q_l1, q_h1
               q(i,jhi+j) = q(i,jhi)
            end do
            end do
         else if (bc(2,2) .eq. 4) then
            do j = 2, ntop
            do i = q_l1, q_h1
               q(i,jhi+j) = q(i,jhi)
            end do
            end do
            if (jhi-2 .ge. js) then
            do i = q_l1, q_h1
               q(i,jhi+1) = (15*q(i,jhi) - 10*q(i,jhi-1) + 3*q(i,jhi-2))
     & * 0.125D0
             end do
            else
             do i = q_l1, q_h1
               q(i,jhi+1) = 0.5D0*(3*q(i,jhi) - q(i,jhi-1))
             end do
            end if
         else if (bc(2,2) .eq. 1) then
            do j = 1, ntop
            do i = q_l1, q_h1
               q(i,jhi+j) = q(i,jhi-j+1)
            end do
            end do
         else if (bc(2,2) .eq. -1) then
            do j = 1, ntop
            do i = q_l1, q_h1
               q(i,jhi+j) = -q(i,jhi-j+1)
            end do
            end do
         end if
      end if

      if ((nlft .gt. 0 .and. bc(1,1) .eq. 4) .and.(nbot .gt. 0 .and. bc(
     &2,1) .eq. 4) ) then

        if (jlo+2 .le. je) then 
          q(ilo-1,jlo-1) = 0.5D0 * 0.125D0 * (15*q(ilo-1,jlo) - 10*q(ilo
     &-1,jlo+1) + 3*q(ilo-1,jlo+2))
        else
          q(ilo-1,jlo-1) = 0.5D0 * 0.5D0 * (3*q(ilo-1,jlo) - q(ilo-1,jlo
     &+1))
        end if

        if (ilo+2 .le. ie) then 
          q(ilo-1,jlo-1) = q(ilo-1,jlo-1) + 0.5D0 * 0.125D0 * (15*q(ilo,
     &jlo-1) - 10*q(ilo+1,jlo-1) + 3*q(ilo+2,jlo-1)) 
        else
          q(ilo-1,jlo-1) = q(ilo-1,jlo-1) + 0.5D0 * 0.5D0 * (3*q(ilo,jlo
     &-1) - q(ilo+1,jlo-1))
        end if

      end if

      if ((nlft .gt. 0 .and. bc(1,1) .eq. 4) .and.(ntop .gt. 0 .and. bc(
     &2,2) .eq. 4) ) then

        if (jhi-2 .ge. js) then 
          q(ilo-1,jhi+1) = 0.5D0 * 0.125D0 * (15*q(ilo-1,jhi) - 10*q(ilo
     &-1,jhi-1) + 3*q(ilo-1,jhi-2))
        else
          q(ilo-1,jhi+1) = 0.5D0 * 0.5D0 * (3*q(ilo-1,jhi) - q(ilo-1,jhi
     &-1))
        end if

        if (ilo+2 .le. ie) then 
          q(ilo-1,jhi+1) = q(ilo-1,jhi+1) + 0.5D0 * 0.125D0 * (15*q(ilo,
     &jhi+1) - 10*q(ilo+1,jhi+1) + 3*q(ilo+2,jhi+1))
        else
          q(ilo-1,jhi+1) = q(ilo-1,jhi+1) + 0.5D0 * 0.5D0 * (3*q(ilo,jhi
     &+1) - q(ilo+1,jhi+1))
        end if
      end if

      if ((nrgt .gt. 0 .and. bc(1,2) .eq. 4) .and.(nbot .gt. 0 .and. bc(
     &2,1) .eq. 4) ) then
        if (jlo+2 .le. je) then 
          q(ihi+1,jlo-1) = 0.5D0 * 0.125D0 * (15*q(ihi+1,jlo) - 10*q(ihi
     &+1,jlo+1) + 3*q(ihi+1,jlo+2))
        else
          q(ihi+1,jlo-1) = 0.5D0 * 0.5D0 * (3*q(ihi+1,jlo) - q(ihi+1,jlo
     &+1))
        end if

        if (ihi-2 .ge. is) then 
          q(ihi+1,jlo-1) = q(ihi+1,jlo-1) + 0.5D0 * 0.125D0 * (15*q(ihi,
     &jlo-1) - 10*q(ihi-1,jlo-1) + 3*q(ihi-2,jlo-1))
        else
          q(ihi+1,jlo-1) = q(ihi+1,jlo-1) + 0.5D0 * 0.5D0 * (3*q(ihi,jlo
     &-1) - q(ihi-1,jlo-1))
        end if
      end if

      if ((nrgt .gt. 0 .and. bc(1,2) .eq. 4) .and.(ntop .gt. 0 .and. bc(
     &2,2) .eq. 4) ) then

        if (jhi-2 .ge. js) then 
          q(ihi+1,jhi+1) = 0.5D0 * 0.125D0 * (15*q(ihi+1,jhi) - 10*q(ihi
     &+1,jhi-1) + 3*q(ihi+1,jhi-2))
        else
          q(ihi+1,jhi+1) = 0.5D0 * 0.5D0 * (3*q(ihi+1,jhi) - q(ihi+1,jhi
     &-1))
        end if

        if (ihi-2 .ge. is) then 
          q(ihi+1,jhi+1) = q(ihi+1,jhi+1) + 0.5D0 * 0.125D0 * (15*q(ihi,
     &jhi+1) - 10*q(ihi-1,jhi+1) + 3*q(ihi-2,jhi+1))
        else
          q(ihi+1,jhi+1) = q(ihi+1,jhi+1) + 0.5D0 * 0.5D0 * (3*q(ihi,jhi
     &+1) - q(ihi-1,jhi+1))
        end if

      end if

      return
      end

      subroutine hoextraptocc(q,q_l1, q_l2, q_h1, q_h2,domlo,domhi,dx,xl
     &o)

      implicit none

      integer    q_l1, q_l2, q_h1, q_h2
      integer    domlo(2), domhi(2)
      DOUBLE PRECISION     xlo(2), dx(2)
      DOUBLE PRECISION     q(q_l1:q_h1, q_l2:q_h2)

      integer    nlft, nrgt, nbot, ntop
      integer    ilo, ihi, jlo, jhi
      integer    i, j
      integer    is, ie, js, je

      nlft = max(0,domlo(1)-q_l1)
      nrgt = max(0,q_h1-domhi(1))
      nbot = max(0,domlo(2)-q_l2)
      ntop = max(0,q_h2-domhi(2))

      is = max(q_l1,domlo(1))
      ie = min(q_h1,domhi(1))
      js = max(q_l2,domlo(2))
      je = min(q_h2,domhi(2))

c
c     Set these to invalid values, they shouldn't be used if not reset
c
      ilo = -10
      jlo = -10
      ihi = 100000000
      jhi = 100000000

c
c     First fill sides.
c
      if (nlft .gt. 0) then
         ilo = domlo(1)
         do i = 2, nlft
            do j = q_l2, q_h2
               q(ilo-i,j) = q(ilo,j) 
            end do 
         end do 
         if (ilo+2 .le. ie) then 
            do j = q_l2, q_h2
               q(ilo-1,j) = 3*q(ilo,j) - 3*q(ilo+1,j) + q(ilo+2,j)
            end do 
         else 
            do j = q_l2, q_h2
               q(ilo-1,j) = 2*q(ilo,j) - q(ilo+1,j)
            end do
         end if
      end if

      if (nrgt .gt. 0) then
         ihi = domhi(1)
         do i = 2, nrgt
            do j = q_l2, q_h2
               q(ihi+i,j) = q(ihi,j)
            end do
         end do
         if (ihi-2 .ge. is) then
            do j = q_l2, q_h2
               q(ihi+1,j) = 3*q(ihi,j) - 3*q(ihi-1,j) + q(ihi-2,j)
            end do
         else
            do j = q_l2, q_h2
               q(ihi+1,j) = 2*q(ihi,j) - q(ihi-1,j)
            end do
         end if
      end if

      if (nbot .gt. 0) then
         jlo = domlo(2)
         do j = 2, nbot
            do i = q_l1, q_h1
               q(i,jlo-j) = q(i,jlo)
            end do
         end do
         if (jlo+2 .le. je) then
            do i = q_l1, q_h1
               q(i,jlo-1) = 3*q(i,jlo) - 3*q(i,jlo+1) + q(i,jlo+2)
            end do
         else
            do i = q_l1, q_h1
               q(i,jlo-1) = 2*q(i,jlo) - q(i,jlo+1)
            end do
         end if
      end if

      if (ntop .gt. 0) then
         jhi = domhi(2)
         do j = 2, ntop
            do i = q_l1, q_h1
               q(i,jhi+j) = q(i,jhi)
            end do
         end do
         if (jhi-2 .ge. js) then
            do i = q_l1, q_h1
               q(i,jhi+1) = 3*q(i,jhi) - 3*q(i,jhi-1) + q(i,jhi-2)
            end do
         else
            do i = q_l1, q_h1
               q(i,jhi+1) = 2*q(i,jhi) - q(i,jhi-1)
            end do
         end if
      end if

      if ((nlft .gt. 0) .and. (nbot .gt. 0)) then
          if (jlo+2 .le. je) then
             q(ilo-1,jlo-1) = 0.5D0 *(3*q(ilo-1,jlo) - 3*q(ilo-1,jlo+1) 
     &+ q(ilo-1,jlo+2))
          else
             q(ilo-1,jlo-1) = 0.5D0 * (2*q(ilo-1,jlo) - q(ilo-1,jlo+1))
          end if

          if (ilo+2 .le. ie) then 
             q(ilo-1,jlo-1) = q(ilo-1,jlo-1) + 0.5D0 *(3*q(ilo,jlo-1) - 
     &3*q(ilo+1,jlo-1) + q(ilo+2,jlo-1)) 
          else
             q(ilo-1,jlo-1) = q(ilo-1,jlo-1) + 0.5D0 *(2*q(ilo,jlo-1) - 
     &q(ilo+1,jlo-1))
          end if
      end if

      if ((nlft .gt. 0) .and. (ntop .gt. 0)) then 
          if (jhi-2 .ge. js) then 
             q(ilo-1,jhi+1) = 0.5D0 *(3*q(ilo-1,jhi) - 3*q(ilo-1,jhi-1) 
     &+ q(ilo-1,jhi-2))
          else
             q(ilo-1,jhi+1) = 0.5D0 * (2*q(ilo-1,jhi) - q(ilo-1,jhi-1))
          end if

          if (ilo+2 .le. ie) then 
             q(ilo-1,jhi+1) = q(ilo-1,jhi+1) + 0.5D0 *(3*q(ilo,jhi+1) - 
     &3*q(ilo+1,jhi+1) + q(ilo+2,jhi+1))
          else
             q(ilo-1,jhi+1) = q(ilo-1,jhi+1) + 0.5D0 *(2*q(ilo,jhi+1) - 
     &q(ilo+1,jhi+1))
          end if
      end if

      if ((nrgt .gt. 0) .and. (nbot .gt. 0)) then 
          if (jlo+2 .le. je) then 
             q(ihi+1,jlo-1) = 0.5D0 *(3*q(ihi+1,jlo) - 3*q(ihi+1,jlo+1) 
     &+ q(ihi+1,jlo+2))
          else
             q(ihi+1,jlo-1) = 0.5D0 * (2*q(ihi+1,jlo) - q(ihi+1,jlo+1))
          end if

          if (ihi-2 .ge. is) then 
             q(ihi+1,jlo-1) = q(ihi+1,jlo-1) + 0.5D0 *(3*q(ihi,jlo-1) - 
     &3*q(ihi-1,jlo-1) + q(ihi-2,jlo-1))
          else
             q(ihi+1,jlo-1) = q(ihi+1,jlo-1) + 0.5D0 *(2*q(ihi,jlo-1) - 
     &q(ihi-1,jlo-1))
          end if
      end if

      if ((nrgt .gt. 0) .and. (ntop .gt. 0)) then 
          if (jhi-2 .ge. js) then 
             q(ihi+1,jhi+1) = 0.5D0 *(3*q(ihi+1,jhi) - 3*q(ihi+1,jhi-1) 
     &+ q(ihi+1,jhi-2))
          else
             q(ihi+1,jhi+1) = 0.5D0 * (2*q(ihi+1,jhi) - q(ihi+1,jhi-1))
          end if

          if (ihi-2 .ge. is) then 
             q(ihi+1,jhi+1) = q(ihi+1,jhi+1) + 0.5D0 *(3*q(ihi,jhi+1) - 
     &3*q(ihi-1,jhi+1) + q(ihi-2,jhi+1))
          else
             q(ihi+1,jhi+1) = q(ihi+1,jhi+1) + 0.5D0 *(2*q(ihi,jhi+1) - 
     &q(ihi-1,jhi+1))
          end if
      end if

      end
