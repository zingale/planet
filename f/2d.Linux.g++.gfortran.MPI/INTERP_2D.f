

c ::: --------------------------------------------------------------
c ::: nbinterp:  node based bilinear interpolation
c :::
c ::: INPUTS/OUTPUTS
c ::: fine        <=>  (modify) fine grid array
c ::: fine_l1, fine_l2, fine_h1, fine_h2   =>  (const)  index limits of fine grid
c ::: fblo,fbhi    =>  (const)  subregion of fine grid to get values
c :::
c ::: crse         =>  (const)  coarse grid data widened by 1 zone
c ::: crse_l1, crse_l2, crse_h1, crse_h2   =>  (const)  index limits of coarse grid
c :::
c ::: lratio(3)    =>  (const)  refinement ratio between levels
c ::: nvar         =>  (const)  number of components in array
c ::: num_slp      =>  (const)  number of types of slopes
c :::
c ::: TEMPORARY ARRAYS
c ::: sl           =>  num_slp 1-D slope arrays
c ::: --------------------------------------------------------------
c :::
      subroutine nbinterp (crse, crse_l1, crse_l2, crse_h1, crse_h2, cb_
     &l1, cb_l2, cb_h1, cb_h2,fine, fine_l1, fine_l2, fine_h1, fine_h2
     &, fb_l1, fb_l2, fb_h1, fb_h2,lratiox, lratioy, nvar,sl, num_slp,
     &actual_comp,actual_state)

      implicit none

      integer crse_l1, crse_l2, crse_h1, crse_h2
      integer cb_l1, cb_l2, cb_h1, cb_h2
      integer fine_l1, fine_l2, fine_h1, fine_h2
      integer fb_l1, fb_l2, fb_h1, fb_h2
      integer lratiox, lratioy, nvar
      integer num_slp
      integer actual_comp,actual_state
      DOUBLE PRECISION  fine(fine_l1:fine_h1, fine_l2:fine_h2,nvar)
      DOUBLE PRECISION  crse(crse_l1:crse_h1, crse_l2:crse_h2,nvar)
      DOUBLE PRECISION  sl(cb_l1:cb_h1,num_slp)

      integer lx, ly
      integer i, j, ifn, jfn, n
      integer ilo, ihi, jlo, jhi
      integer jstrtFine, jstopFine, istrtFine, istopFine

      DOUBLE PRECISION fx, fy
      DOUBLE PRECISION RX, RY, RXY
      DOUBLE PRECISION dx0, d0x, dx1
      DOUBLE PRECISION slope

      slope(i,j,n,fx,fy) = crse(i,j,n) +fx*sl(i,1) + fy*sl(i,2) + fx*fy*
     &sl(i,3)

      RX = 1.0D0/dble(lratiox)
      RY = 1.0D0/dble(lratioy)
      RXY = RX*RY
c
c     NOTES:
c         1) (i, j) loop over the coarse cells
c         2) ?strtFine and ?stopFine are the beginning and ending fine cell
c            indices corresponding to the current coarse cell.  ?stopFine
c            is restricted for the last coarse cell in each direction since
c            for this cell we only need to do the face and not the fine nodes
c            inside this cell.
c         3) (lx, ly) as well as ?lo and ?hi refer to the fine node indices
c            as an offset from ?strtFine.
c
      do 100 n = 1, nvar
        do 120 j = cb_l2, cb_h2
          jstrtFine = j * lratioy
          jstopFine = jstrtFine + lratioy - 1
          if (j .eq. cb_h2) jstopFine = jstrtFine

          jlo = max(fb_l2,jstrtFine) - jstrtFine
          jhi = min(fb_h2,jstopFine) - jstrtFine
c
c         ::::: compute slopes :::::
c
c         NOTE: The IF logic in the calculation of the slopes is to
c               prevent stepping out of bounds on the coarse data when
c               computing the slopes on the ARG_H?(cb) cells.  These
c               slopes actually are not used since they are multiplied by
c               0.0D0.
c
          do i = cb_l1, cb_h1
            dx0 = 0.0D0
            if (i .NE. cb_h1) dx0 = crse(i+1,j,n) - crse(i,j,n)

            d0x = 0.0D0
            if (j .NE. cb_h2) d0x = crse(i,j+1,n) - crse(i,j,n)

            dx1 = 0.0D0
            if (i .NE. cb_h1 .and. j .NE. cb_h2)dx1 = crse(i+1,j+1,n) - 
     &crse(i,j+1,n)

            sl(i,1) = RX*dx0
            sl(i,2) = RY*d0x
            sl(i,3) = RXY*(dx1 - dx0)
          end do
c
c         ::::: compute fine strip of interpolated data
c
          do ly = jlo, jhi
            jfn = lratioy * j + ly
            fy = dble(ly)

            do i = cb_l1, cb_h1
              istrtFine = i * lratiox
              istopFine = istrtFine + lratiox - 1
              if (i .eq. cb_h1) istopFine = istrtFine

              ilo = max(fb_l1,istrtFine) - istrtFine
              ihi = min(fb_h1,istopFine) - istrtFine

              do lx = ilo, ihi
                ifn = lratiox * i + lx
                fx = dble(lx)

                fine(ifn,jfn,n) = slope(i,j,n,fx,fy)
              end do
            end do
          end do
120     continue
100   continue

      end
c ::: 
c ::: --------------------------------------------------------------
c ::: cbinterp:  cell centered bilinear interpolation
c ::: 
c ::: NOTE: it is assumed that the coarse grid array is
c ::: large enough to define interpolated values
c ::: in the region fblo:fbhi on the fine grid
c ::: 
c ::: Inputs/Outputs
c ::: fine        <=>  (modify) fine grid array
c ::: fine_l1, fine_l2, fine_h1, fine_h2   =>  (const)  index limits of fine grid
c ::: fb_l1, fb_l2, fb_h1, fb_h2     =>  (const)  subregion of fine grid to get values
c ::: 
c ::: crse         =>  (const)  coarse grid data 
c ::: crse_l1, crse_l2, crse_h1, crse_h2   =>  (const)  index limits of coarse grid
c ::: 
c ::: lratio(2)    =>  (const)  refinement ratio between levels
c ::: nvar         =>  (const)  number of components in array
c ::: 
c ::: TEMPORARY ARRAYS
c ::: slx,sly,slxy =>  1-D slope arrays
c ::: strip        =>  1-D temp array
c ::: --------------------------------------------------------------
c ::: 
      subroutine cbinterp (crse, crse_l1, crse_l2, crse_h1, crse_h2, cb_
     &l1, cb_l2, cb_h1, cb_h2,fine, fine_l1, fine_l2, fine_h1, fine_h2
     &, fb_l1, fb_l2, fb_h1, fb_h2,lratiox, lratioy, nvar,sl, num_slp,
     & strip, strip_lo, strip_hi,actual_comp,actual_state)

      implicit none

      integer crse_l1, crse_l2, crse_h1, crse_h2
      integer cb_l1, cb_l2, cb_h1, cb_h2
      integer fine_l1, fine_l2, fine_h1, fine_h2
      integer fb_l1, fb_l2, fb_h1, fb_h2
      integer lratiox, lratioy, nvar
      integer num_slp
      integer actual_comp,actual_state
      integer strip_lo, strip_hi
      DOUBLE PRECISION  fine(fine_l1:fine_h1, fine_l2:fine_h2, nvar)
      DOUBLE PRECISION  crse(crse_l1:crse_h1, crse_l2:crse_h2, nvar)
      DOUBLE PRECISION  sl(cb_l1:cb_h1,num_slp)
      DOUBLE PRECISION  strip(strip_lo:strip_hi)

      integer lx, ly, hratx, hraty, ic, jc, jfn, jfc, i, n
      DOUBLE PRECISION x, y, denomx, denomy

      denomx = 1.0D0/dble(2*lratiox)
      denomy = 1.0D0/dble(2*lratioy)

      hratx = lratiox/2
      hraty = lratioy/2

      do n = 1, nvar 
         do jc = cb_l2, cb_h2-1 
            do ic = cb_l1, cb_h1-1
               sl(ic,1) = crse(ic+1,jc,n)-crse(ic,jc,n)
               sl(ic,2) = crse(ic,jc+1,n)-crse(ic,jc,n)
               sl(ic,3) = crse(ic+1,jc+1,n)-crse(ic+1,jc,n)- crse(ic ,jc
     &+1,n)+crse(ic ,jc,n)
            end do
            do ly = 0, lratioy-1 
               jfn = jc*lratioy + ly
               jfc = jfn + hraty
               if (jfc .ge. fb_l2  .and.  jfc .le. fb_h2) then
                  y = denomy*(2.0D0*ly + 1.0D0)
                  do lx = 0, lratiox-1
                     do ic = cb_l1, cb_h1-1
                        i = ic*lratiox + lx
                        x = denomx*(2.0D0*lx + 1.0D0)
                        strip(i) = crse(ic,jc,n) + x*sl(ic,1) +y*sl(ic,2
     &) + x*y*sl(ic,3)
                     end do
                  end do
                  do i = fb_l1, fb_h1 
                     fine(i,jfc,n) = strip(i-hratx)
                  end do
               end if
            end do
         end do
      end do

      end

c ::: 
c ::: --------------------------------------------------------------
c ::: linccinterp:   linear conservative interpolation from coarse grid to
c ::: subregion of fine grid defined by (fblo,fbhi)
c ::: 
c ::: The interpolation is linear in that it uses a
c ::: a limiting scheme that preserves the value of 
c ::: any linear combination of the
c ::: coarse grid data components--e.g.,
c ::: if sum_ivar a(ic,jc,ivar)*fab(ic,jc,ivar) = 0, then
c ::: sum_ivar a(ic,jc,ivar)*fab(if,jf,ivar) = 0 is satisfied
c ::: in all fine cells if,jf covering coarse cell ic,jc.
c ::: 
c ::: If lin_limit = 0, the interpolation scheme is identical to
c ::: the used in ccinterp for limslope=1; the results should
c ::: be exactly the same -- difference = hard 0.
c ::: 
c ::: Inputs/Outputs
c ::: fine        <=>  (modify) fine grid array
c ::: flo,fhi      =>  (const)  index limits of fine grid
c ::: fblo,fbhi    =>  (const)  subregion of fine grid to get values
c ::: nvar         =>  (const)  number of variables in state vector
c ::: lratio(2)    =>  (const)  refinement ratio between levels
c ::: 
c ::: crse         =>  (const)  coarse grid data widended by 1 zone
c ::: clo,chi      =>  (const)  index limits of crse grid
c ::: cslo,cshi    =>  (const)  coarse grid index limits where
c :::				slopes are to be defined. This is
c :::				the projection of (fblo,fbhi) down
c :::				to the coarse level 
c ::: ucslope      =>  (modify) temp array of unlimited coarse grid slopes
c ::: lcslope      =>  (modify) temp array of limited coarse grid slopes
c ::: slope_factor =>  (modify) temp array of slope limiting factors
c ::: lin_limit    =>  (const)  != 0 => do linear slope limiting scheme
c :::
c ::: --------------------------------------------------------------
c ::: 
       subroutine linccinterp (fine, fine_l1, fine_l2, fine_h1, fine_h2,
     & fblo, fbhi,fvcb_l1, fvcb_l2, fvcb_h1, fvcb_h2,crse, crse_l1, c
     &rse_l2, crse_h1, crse_h2, cvcb_l1, cvcb_l2, cvcb_h1, cvcb_h2,uc
     &_xslope, lc_xslope, xslope_factor,uc_yslope, lc_yslope, yslope_
     &factor,cslope_l1, cslope_l2, cslope_h1, cslope_h2,cslopelo, csl
     &opehi,nvar, lratiox, lratioy,bc, lim_slope, lin_limit,fvcx, fvc
     &y, cvcx, cvcy,voffx, voffy, alpha, cmax, cmin,actual_comp,actua
     &l_state)

       implicit none

       integer fine_l1, fine_l2, fine_h1, fine_h2
       integer crse_l1, crse_l2, crse_h1, crse_h2
       integer fvcb_l1, fvcb_l2, fvcb_h1, fvcb_h2
       integer cvcb_l1, cvcb_l2, cvcb_h1, cvcb_h2
       integer cslope_l1, cslope_l2, cslope_h1, cslope_h2
       integer fblo(2), fbhi(2)
       integer cslopelo(2), cslopehi(2)
       integer lratiox, lratioy, nvar
       integer lim_slope, lin_limit
       integer bc(2,2,nvar)
       integer actual_comp,actual_state
       DOUBLE PRECISION fine(fine_l1:fine_h1, fine_l2:fine_h2,nvar)
       DOUBLE PRECISION crse(crse_l1:crse_h1, crse_l2:crse_h2, nvar)
       DOUBLE PRECISION uc_xslope(cslope_l1:cslope_h1, cslope_l2:cslope_
     &h2,nvar)
       DOUBLE PRECISION lc_xslope(cslope_l1:cslope_h1, cslope_l2:cslope_
     &h2,nvar)
       DOUBLE PRECISION xslope_factor(cslope_l1:cslope_h1, cslope_l2:csl
     &ope_h2)
       DOUBLE PRECISION uc_yslope(cslope_l1:cslope_h1, cslope_l2:cslope_
     &h2,nvar)
       DOUBLE PRECISION lc_yslope(cslope_l1:cslope_h1, cslope_l2:cslope_
     &h2,nvar)
       DOUBLE PRECISION yslope_factor(cslope_l1:cslope_h1, cslope_l2:csl
     &ope_h2)
       DOUBLE PRECISION alpha(cslope_l1:cslope_h1, cslope_l2:cslope_h2,n
     &var)
       DOUBLE PRECISION cmax(cslope_l1:cslope_h1, cslope_l2:cslope_h2,nv
     &ar)
       DOUBLE PRECISION cmin(cslope_l1:cslope_h1, cslope_l2:cslope_h2,nv
     &ar)
       DOUBLE PRECISION fvcx(fvcb_l1:fvcb_h1)
       DOUBLE PRECISION fvcy(fvcb_l2:fvcb_h2)
       DOUBLE PRECISION voffx(fvcb_l1:fvcb_h1)
       DOUBLE PRECISION voffy(fvcb_l2:fvcb_h2)
       DOUBLE PRECISION cvcx(cvcb_l1:cvcb_h1)
       DOUBLE PRECISION cvcy(cvcb_l2:cvcb_h2)

       integer n
       integer i, ic
       integer j, jc
       DOUBLE PRECISION cen, forw, back, slp
       DOUBLE PRECISION factorn, denom
       DOUBLE PRECISION fxcen, cxcen, fycen, cycen
       DOUBLE PRECISION orig_corr_fact,corr_fact
       DOUBLE PRECISION dummy_fine
       logical xok, yok
       integer ncbx, ncby
       integer ioff,joff
       integer voff_lo(2),voff_hi(2)

       ncbx = cslopehi(1)-cslopelo(1)+1
       ncby = cslopehi(2)-cslopelo(2)+1

       voff_lo(1) = cslopelo(1) * lratiox
       voff_lo(2) = cslopelo(2) * lratioy
       voff_hi(1) = (cslopehi(1)+1) * lratiox - 1
       voff_hi(2) = (cslopehi(2)+1) * lratioy - 1

       xok = (ncbx .ge. 2)
       yok = (ncby .ge. 2)

       do j = voff_lo(2),voff_hi(2)
         jc = (j+lratioy*iabs(j))/lratioy-iabs(j)
         fycen = 0.5D0*(fvcy(j)+fvcy(j+1))
         cycen = 0.5D0*(cvcy(jc)+cvcy(jc+1))
         voffy(j) = (fycen-cycen)/(cvcy(jc+1)-cvcy(jc))
       end do
       do i = voff_lo(1),voff_hi(1)
          ic = (i+lratiox*iabs(i))/lratiox-iabs(i)
          fxcen = 0.5D0*(fvcx(i)+fvcx(i+1))
          cxcen = 0.5D0*(cvcx(ic)+cvcx(ic+1))
          voffx(i) = (fxcen-cxcen)/(cvcx(ic+1)-cvcx(ic))
       end do

       do n = 1, nvar
c
c ...     Prevent underflow for small crse values.
c
          do j = cslopelo(2)-1,cslopehi(2)+1
             do i = cslopelo(1)-1, cslopehi(1)+1
                crse(i,j,n) = merge(crse(i,j,n),0.0D0,abs(crse(i,j,n)).g
     &t.1.0d-50)
             end do
          end do
c
c ...     Initialize alpha = 1 and define cmax and cmin as neighborhood max/mins.
c
          do j = cslopelo(2),cslopehi(2)
             do i = cslopelo(1), cslopehi(1)
                alpha(i,j,n) = 1.d0
                cmax(i,j,n) = crse(i,j,n)
                cmin(i,j,n) = crse(i,j,n)
                do joff = -1,1
                do ioff = -1,1
                  cmax(i,j,n) = max(cmax(i,j,n),crse(i+ioff,j+joff,n))
                  cmin(i,j,n) = min(cmin(i,j,n),crse(i+ioff,j+joff,n))
                end do
                end do
             end do
          end do

       end do
c
c ...  Compute unlimited and limited slopes
c
       do n = 1, nvar

          do j=cslopelo(2), cslopehi(2)
             do i=cslopelo(1), cslopehi(1)
                uc_xslope(i,j,n) = 0.5D0*(crse(i+1,j,n)-crse(i-1,j,n))
                cen  = uc_xslope(i,j,n)
                forw = 2.0D0*(crse(i+1,j,n)-crse(i,j,n))
                back = 2.0D0*(crse(i,j,n)-crse(i-1,j,n))
                slp  = min(abs(forw),abs(back))
                slp  = merge(slp,0.0D0,forw*back>=0.0D0)
                lc_xslope(i,j,n)=sign(1.0D0,cen)*min(slp,abs(cen))
             end do
          end do

          if (bc(1,1,n) .eq. 3 .or. bc(1,1,n).eq.4) then
            i = cslopelo(1)
            if (xok) then
                do j=cslopelo(2), cslopehi(2)
                   uc_xslope(i,j,n) = -16.0D0/15.0D0*crse(i-1,j,n)+ 0.5D
     &0*crse(i,j,n)+ 0.66666666666666667D0*crse(i+1,j,n) 
     &- 0.1D0*crse(i+2,j,n)
                end do
            else
                do j=cslopelo(2), cslopehi(2)
                   uc_xslope(i,j,n) = 0.25D0 * (crse(i+1,j,n) + 5.0D0*cr
     &se(i,j,n) - 6.0D0*crse(i-1,j,n) )
                end do
            endif
            do j=cslopelo(2), cslopehi(2)
               cen  = uc_xslope(i,j,n)
               forw = 2.0D0*(crse(i+1,j,n)-crse(i,j,n))
               back = 2.0D0*(crse(i,j,n)-crse(i-1,j,n))
               slp  = min(abs(forw),abs(back))
               slp  = merge(slp,0.0D0,forw*back>=0.0D0)
               lc_xslope(i,j,n)=sign(1.0D0,cen)*min(slp,abs(cen))
            end do
          end if

          if (bc(1,2,n) .eq. 3 .or. bc(1,2,n).eq.4) then
            i = cslopehi(1)
            if (xok) then
                do j=cslopelo(2), cslopehi(2)
                   uc_xslope(i,j,n) = 16.0D0/15.0D0*crse(i+1,j,n)- 0.5D0
     &*crse(i,j,n)- 0.66666666666666667D0*crse(i-1,j,n) +
     & 0.1D0*crse(i-2,j,n)
                end do
            else
                do j=cslopelo(2), cslopehi(2)
                   uc_xslope(i,j,n) = -0.25D0 * (crse(i-1,j,n) + 5.0D0*c
     &rse(i,j,n) - 6.0D0*crse(i+1,j,n) )
                end do
            endif
            do j=cslopelo(2), cslopehi(2)
               cen  = uc_xslope(i,j,n)
               forw = 2.0D0*(crse(i+1,j,n)-crse(i,j,n))
               back = 2.0D0*(crse(i,j,n)-crse(i-1,j,n))
               slp  = min(abs(forw),abs(back))
               slp  = merge(slp,0.0D0,forw*back>=0.0D0)
               lc_xslope(i,j,n)=sign(1.0D0,cen)*min(slp,abs(cen))
            end do
          end if

          do j=cslopelo(2), cslopehi(2)
             do i=cslopelo(1), cslopehi(1)
                uc_yslope(i,j,n) = 0.5D0*(crse(i,j+1,n)-crse(i,j-1,n))
                cen  = uc_yslope(i,j,n)
                forw = 2.0D0*(crse(i,j+1,n)-crse(i,j,n))
                back = 2.0D0*(crse(i,j,n)-crse(i,j-1,n))
                slp  = min(abs(forw),abs(back))
                slp  = merge(slp,0.0D0,forw*back>=0.0D0)
                lc_yslope(i,j,n)=sign(1.0D0,cen)*min(slp,abs(cen))
             end do
          end do

          if (bc(2,1,n) .eq. 3 .or. bc(2,1,n).eq.4) then
             j = cslopelo(2)
             if (yok) then
                do i=cslopelo(1), cslopehi(1)
                   uc_yslope(i,j,n) = -16.0D0/15.0D0*crse(i,j-1,n)+ 0.5D
     &0*crse(i,j,n)+ 0.66666666666666667D0*crse(i,j+1,n) 
     &- 0.1D0*crse(i,j+2,n)
                end do
             else
                do i=cslopelo(1), cslopehi(1)
                   uc_yslope(i,j,n) = 0.25D0 * (crse(i,j+1,n) + 5.0D0*cr
     &se(i,j,n) - 6.0D0*crse(i,j-1,n) )
                end do
             endif
             do i=cslopelo(1), cslopehi(1)
                cen  = uc_yslope(i,j,n)
                forw = 2.0D0*(crse(i,j+1,n)-crse(i,j,n))
                back = 2.0D0*(crse(i,j,n)-crse(i,j-1,n))
                slp  = min(abs(forw),abs(back))
                slp  = merge(slp,0.0D0,forw*back>=0.0D0)
                lc_yslope(i,j,n)=sign(1.0D0,cen)*min(slp,abs(cen))
             end do
          end if

          if (bc(2,2,n) .eq. 3 .or. bc(2,2,n).eq.4) then
             j = cslopehi(2)
             if (yok) then
                do i=cslopelo(1), cslopehi(1)
                   uc_yslope(i,j,n) = 16.0D0/15.0D0*crse(i,j+1,n)- 0.5D0
     &*crse(i,j,n)- 0.66666666666666667D0*crse(i,j-1,n) +
     & 0.1D0*crse(i,j-2,n)
                end do
             else
                do i=cslopelo(1), cslopehi(1)
                   uc_yslope(i,j,n) = -0.25D0 * (crse(i,j-1,n) + 5.0D0*c
     &rse(i,j,n) - 6.0D0*crse(i,j+1,n) )
                end do
             endif
             do i=cslopelo(1), cslopehi(1)
                cen  = uc_yslope(i,j,n)
                forw = 2.0D0*(crse(i,j+1,n)-crse(i,j,n))
                back = 2.0D0*(crse(i,j,n)-crse(i,j-1,n))
                slp  = min(abs(forw),abs(back))
                slp  = merge(slp,0.0D0,forw*back>=0.0D0)
                lc_yslope(i,j,n)=sign(1.0D0,cen)*min(slp,abs(cen))
             end do
          end if

       end do

       if (lim_slope.eq.0) then
c
c ...    Do the interpolation using unlimited slopes.
c
          do n = 1, nvar
             do j = fblo(2), fbhi(2)
                jc = (j+lratioy*iabs(j))/lratioy-iabs(j)
                do i = fblo(1), fbhi(1)
                   ic = (i+lratiox*iabs(i))/lratiox-iabs(i)
                   fine(i,j,n) = crse(ic,jc,n) + voffx(i)*uc_xslope(ic,j
     &c,n)+ voffy(j)*uc_yslope(ic,jc,n)
                end do
             end do
          end do

       else 

         if (lin_limit.eq.1) then
c
c ...      compute linear limited slopes
c          Note that the limited and the unlimited slopes
c          have the same sign, and it is assumed that they do.
c
c ... --> compute slope factors
c
           do j=cslopelo(2), cslopehi(2)
             do i=cslopelo(1), cslopehi(1)
                xslope_factor(i,j) = 1.0D0
                yslope_factor(i,j) = 1.0D0
             end do
           end do

           do n = 1, nvar
             do j=cslopelo(2), cslopehi(2)
                do i=cslopelo(1), cslopehi(1)
                   denom = uc_xslope(i,j,n)
                   denom = merge(denom,1.0D0,denom.ne.0.0D0)
                   factorn = lc_xslope(i,j,n)/denom
                   factorn = merge(1.0D0,factorn,denom.eq.0.0D0)
                   xslope_factor(i,j) = min(xslope_factor(i,j),factorn)
                   denom = uc_yslope(i,j,n)
                   denom = merge(denom,1.0D0,denom.ne.0.0D0)
                   factorn = lc_yslope(i,j,n)/denom
                   factorn = merge(1.0D0,factorn,denom.eq.0.0D0)
                   yslope_factor(i,j) = min(yslope_factor(i,j),factorn)
                end do
             end do
           end do
c
c ... -->  compute linear limited slopes
c
           do n = 1, nvar
             do j=cslopelo(2), cslopehi(2)
                do i=cslopelo(1), cslopehi(1)
                   lc_xslope(i,j,n) = xslope_factor(i,j)*uc_xslope(i,j,n
     &)
                   lc_yslope(i,j,n) = yslope_factor(i,j)*uc_yslope(i,j,n
     &)
                end do
             end do
           end do

         else
c
c          Limit slopes so as to not introduce new maxs or mins.
c

            do n = 1, nvar
               do j = voff_lo(2),voff_hi(2)
                  jc = (j+lratioy*iabs(j))/lratioy-iabs(j)

                  do i = voff_lo(1),voff_hi(1)
                     ic = (i+lratiox*iabs(i))/lratiox-iabs(i)

                     orig_corr_fact = voffx(i)*lc_xslope(ic,jc,n)+ voffy
     &(j)*lc_yslope(ic,jc,n) 
                     dummy_fine = crse(ic,jc,n) + orig_corr_fact
                     if ( (dummy_fine .gt. cmax(ic,jc,n)) .and.(abs(orig
     &_corr_fact) .gt. 1.e-10*abs(crse(ic,jc,n)))) then
                        corr_fact = (cmax(ic,jc,n) - crse(ic,jc,n)) / or
     &ig_corr_fact
                        alpha(ic,jc,n) = min(alpha(ic,jc,n),corr_fact)
                     endif
                     if ( (dummy_fine .lt. cmin(ic,jc,n)) .and.(abs(orig
     &_corr_fact) .gt. 1.e-10*abs(crse(ic,jc,n)))) then
                        corr_fact = (cmin(ic,jc,n) - crse(ic,jc,n)) / or
     &ig_corr_fact
                        alpha(ic,jc,n) = min(alpha(ic,jc,n),corr_fact)
                     endif

                  end do
               end do
            end do

         end if
c
c ...    Do the interpolation with limited slopes.
c
          do n = 1, nvar
            do j = fblo(2), fbhi(2)
               jc = (j+lratioy*iabs(j))/lratioy-iabs(j)
               do i = fblo(1), fbhi(1)
                  ic = (i+lratiox*iabs(i))/lratiox-iabs(i)
                  fine(i,j,n) = crse(ic,jc,n) + alpha(ic,jc,n)*( voffx(i
     &)*lc_xslope(ic,jc,n)+voffy(j)*lc_yslope(ic,jc,n) )
               end do
            end do
          end do

       end if

       end

      subroutine cqinterp (fine, fine_l1, fine_l2, fine_h1, fine_h2, fb_
     &l1, fb_l2, fb_h1, fb_h2,nvar, lratiox, lratioy, crse, clo, chi, 
     &cb_l1, cb_l2, cb_h1, cb_h2,fslo, fshi, cslope, clen, fslope, fda
     &t,flen, voff, bc, limslope,fvcx, fvcy, cvcx, cvcy,actual_comp,ac
     &tual_state)

      implicit none

      integer fine_l1, fine_l2, fine_h1, fine_h2
      integer fslo(2), fshi(2)
      integer fb_l1, fb_l2, fb_h1, fb_h2
      integer cb_l1, cb_l2, cb_h1, cb_h2
      integer clo, chi
      integer lratiox, lratioy, nvar, clen, flen, limslope
      integer bc(2,2,nvar)
      integer actual_comp,actual_state
      DOUBLE PRECISION fine(fine_l1:fine_h1, fine_l2:fine_h2,nvar)
      DOUBLE PRECISION crse(clo:chi, nvar)
      DOUBLE PRECISION cslope(clo:chi, 5)
      DOUBLE PRECISION fslope(flen, 5)
      DOUBLE PRECISION fdat(flen)
      DOUBLE PRECISION voff(flen)
      DOUBLE PRECISION fvcx(fb_l1:fb_h1+1)
      DOUBLE PRECISION fvcy(fb_l2:fb_h2+1)
      DOUBLE PRECISION cvcx(cb_l1:cb_h1+1)
      DOUBLE PRECISION cvcy(cb_l2:cb_h2+1)

      integer n, fn
      integer i, ic, ioff
      integer j, jc, joff
      integer ist, jst
      DOUBLE PRECISION cen
      DOUBLE PRECISION fcen, ccen
      DOUBLE PRECISION diffxy,diffxx,diffyy
      DOUBLE PRECISION yoff
      integer ncbx, ncby
      integer ncsx
      integer jslo
      integer icc, istart, iend
      logical xok, yok

      ncbx = cb_h1-cb_l1+1
      ncby = cb_h2-cb_l2+1
      xok = (ncbx .ge. 2)
      yok = (ncby .ge. 2)
      ncsx = ncbx+2
      ist = 1
      jst = ncsx
      jslo = cb_l2-1

      do i = fb_l1, fb_h1 
         fn = i-fslo(1)+1
         ic = (i+lratiox*iabs(i))/lratiox-iabs(i)
         fcen = 0.5D0*(fvcx(i)+fvcx(i+1))
         ccen = 0.5D0*(cvcx(ic)+cvcx(ic+1))
         voff(fn) = (fcen-ccen)/(cvcx(ic+1)-cvcx(ic))
      end do   

      do n = 1, nvar 
        do i = clo, chi 
          crse(i,n) = merge(crse(i,n),0.0D0,abs(crse(i,n)).gt.1.0d-50)
        end do
      end do

      do 290 n = 1, nvar 

            do i = 1, clen 
               cen = 0.5D0*(crse(i+ist,n)-crse(i-ist,n))
               diffxy = 0.25D0*(crse(i+ist+jst,n)+crse(i-ist-jst,n)-crse
     &(i-ist+jst,n)-crse(i+ist-jst,n))
               diffxx = crse(i+ist,n)-2.0D0*crse(i,n)+crse(i-ist,n)
               cslope(i,1)=cen
               cslope(i,3)=diffxx
               cslope(i,5)=diffxy
            end do
            if (xok) then
               if (bc(1,1,n) .eq. 3 .or. bc(1,1,n).eq.4) then
                  do i = 1, clen, jst 
                     cen = -16.0D0/15.0D0*crse(i-ist,n) + 0.5D0*crse(i,n
     &)+ 0.66666666666666667D0*crse(i+ist,n) - 0.1D0*cr
     &se(i+2*ist,n)
                     cslope(i,1)=cen
                     cslope(i,3)=0.0D0
                     cslope(i,5)=0.0D0
                  end do
               end if
               if (bc(1,2,n) .eq. 3 .or. bc(1,2,n).eq.4) then
                  do i = ncbx, clen, jst 
                     cen = 16.0D0/15.0D0*crse(i+ist,n) - 0.5D0*crse(i,n)
     &- 0.66666666666666667D0*crse(i-ist,n) + 0.1D0*crs
     &e(i-2*ist,n)
                     cslope(i,1)=cen
                     cslope(i,3)=0.0D0
                     cslope(i,5)=0.0D0
                  end do
               end if
            end if

            do i = 1, clen 
               cen  = 0.5D0*(crse(i+jst,n)-crse(i-jst,n))
               diffyy = crse(i+jst,n)-2.0D0*crse(i,n)+crse(i-jst,n)
               cslope(i,2)=cen
               cslope(i,4)=diffyy
            end do
            if (yok) then
               if (bc(2,1,n) .eq. 3 .or. bc(2,1,n).eq.4) then
                  do i = 1, ncbx 
                     cen = -16.0D0/15.0D0*crse(i-jst,n) + 0.5D0*crse(i,n
     &)+ 0.66666666666666667D0*crse(i+jst,n) - 0.1D0*cr
     &se(i+2*jst,n)
                     cslope(i,2)=cen
                     cslope(i,4)=0.0D0
                     cslope(i,5)=0.0D0
                  end do
               end if
               if (bc(2,2,n) .eq. 3 .or. bc(2,2,n).eq.4) then
                  do i = clen-ncbx,clen 
                     cen = 16.0D0/15.0D0*crse(i+jst,n) - 0.5D0*crse(i,n)
     &- 0.66666666666666667D0*crse(i-jst,n) + 0.1D0*crs
     &e(i-2*jst,n)
                     cslope(i,2)=cen
                     cslope(i,4)=0.0D0
                     cslope(i,5)=0.0D0
                  end do
               end if
            end if

            do 360 jc = cb_l2, cb_h2 
               do 370 ioff = 1, lratiox 
                  icc = clo + ist + jst*(jc-jslo)
                  istart = ioff
                  iend = ioff + (ncbx-1)*lratiox
                  do 380 fn = istart, iend, lratiox 
                     fslope(fn,1) = cslope(icc,1)
                     fslope(fn,2) = cslope(icc,2)
                     fslope(fn,3) = cslope(icc,3)
                     fslope(fn,4) = cslope(icc,4)
                     fslope(fn,5) = cslope(icc,5)
                     fdat(fn) = crse(icc,n)
                     icc = icc + ist
380               continue
370            continue

               do 390 joff = 0, lratioy-1 
                  j = lratioy*jc + joff
                  if ((j.lt.fb_l2).or.(j.gt.fb_h2)) goto 390
                  fcen = 0.5D0*(fvcy(j)+fvcy(j+1))
                  ccen = 0.5D0*(cvcy(jc)+cvcy(jc+1))
                  yoff = (fcen-ccen)/(cvcy(jc+1)-cvcy(jc))

                  do 400 i = fb_l1, fb_h1 
                     fn = i-fslo(1)+1
                     fine(i,j,n) = fdat(fn) + voff(fn) *fslope(fn,1)+ yo
     &ff *fslope(fn,2)+ 0.5D0*voff(fn)*voff(fn)*fslope(
     &fn,3)+ 0.5D0*yoff *yoff *fslope(fn,4)+ voff(fn)*y
     &off *fslope(fn,5)
400               continue
390            continue
360         continue

290   continue

      end
c ::: 
c ::: --------------------------------------------------------------
c ::: pcinterp:  cell centered piecewise constant interpolation
c ::: 
c ::: Inputs/Outputs
c ::: fine        <=>  (modify) fine grid array
c ::: flo,fhi      =>  (const)  index limits of fine grid
c ::: fblo,fbhi    =>  (const)  subregion of fine grid to get values
c ::: 
c ::: crse         =>  (const)  coarse grid data 
c ::: clo,chi      =>  (const)  index limits of coarse grid
c ::: cblo,cbhi    =>  (const) coarse grid region containing fblo,fbhi
c ::: 
c ::: longdir      =>  (const)  which index direction is longest (1 or 2)
c ::: lratio(2)    =>  (const)  refinement ratio between levels
c ::: nvar         =>  (const)  number of components in array
c ::: 
c ::: TEMPORARY ARRAYS
c ::: ftmp         =>  1-D temp array
c ::: --------------------------------------------------------------
c ::: 
      subroutine pcinterp (crse,crse_l1, crse_l2, crse_h1, crse_h2,cblo,
     &cbhi,fine,fine_l1, fine_l2, fine_h1, fine_h2,fblo,fbhi,longdir,l
     &ratiox,lratioy,nvar,ftmp,ftmp_lo,ftmp_hi,actual_comp,actual_stat
     &e)

      implicit none

      integer crse_l1, crse_l2, crse_h1, crse_h2
      integer cblo(2), cbhi(2)
      integer fine_l1, fine_l2, fine_h1, fine_h2
      integer fblo(2), fbhi(2)
      integer ftmp_lo, ftmp_hi
      integer nvar, lratiox, lratioy, longdir
      integer actual_comp,actual_state
      DOUBLE PRECISION  crse(crse_l1:crse_h1, crse_l2:crse_h2, nvar)
      DOUBLE PRECISION  fine(fine_l1:fine_h1, fine_l2:fine_h2, nvar)
      DOUBLE PRECISION  ftmp(ftmp_lo:ftmp_hi)

      integer i, j, ic, jc, ioff, joff, n

      if (longdir .eq. 1) then
         do n = 1, nvar
         do jc = cblo(2), cbhi(2)
            j = jc*lratioy
            do ioff = 0, lratiox-1
               do ic = cblo(1), cbhi(1)
                  i = lratiox*ic + ioff
                  ftmp(i) = crse(ic,jc,n)
               end do
            end do
            do joff = 0, lratioy-1
               j = lratioy*jc + joff
               if (j.ge.fblo(2).and.j.le.fbhi(2)) then
                  do i = fblo(1), fbhi(1)
                     fine(i,j,n) = ftmp(i)
                  end do
               end if
            end do
         end do
         end do
      else
         do n = 1, nvar
         do ic = cblo(1), cbhi(1)
            i = ic*lratiox
            do joff = 0, lratioy-1
               do jc = cblo(2), cbhi(2)
                  j = lratioy*jc + joff
                  ftmp(j) = crse(ic,jc,n)
               end do
            end do
            do ioff = 0, lratiox-1
               i = lratiox*ic + ioff
               if (i.ge.fblo(1).and.i.le.fbhi(1)) then
                  do j = fblo(2), fbhi(2)
                     fine(i,j,n) = ftmp(j)
                  end do
               end if
            end do
         end do
         end do
      end if
      end

c ::: 
c ::: --------------------------------------------------------------
c ::: protect_interp:   redo interpolation if the result of linccinterp
c ::: generates under- or overshoots.
c ::: 
c ::: 
c ::: Inputs/Outputs
c ::: fine        <=>  (modify) fine grid array
c ::: flo,fhi      =>  (const)  index limits of fine grid
c ::: fblo,fbhi    =>  (const)  subregion of fine grid to get values
c ::: cblo,cbhi    =>  (const)  coarse equivalent of fblo,fbhi
c ::: nvar         =>  (const)  number of variables in state vector
c ::: lratio(3)    =>  (const)  refinement ratio between levels
c ::: 
c ::: crse         =>  (const)  coarse grid data widended by 1 zone
c ::: clo,chi      =>  (const)  index limits of crse grid
c :::
c ::: --------------------------------------------------------------
c ::: 
      subroutine printerp (fine, fine_l1, fine_l2, fine_h1, fine_h2, fbl
     &o, fbhi, crse, crse_l1, crse_l2, crse_h1, crse_h2, cblo, cbhi,fv
     &cx, fvcy, fb_l1, fb_l2, fb_h1, fb_h2,cvcx, cvcy, cb_l1, cb_l2, c
     &b_h1, cb_h2,fine_state, state_l1, state_l2, state_h1, state_h2, 
     &nvar, lratiox, lratioy, bc)

      implicit none

      integer fine_l1, fine_l2, fine_h1, fine_h2
      integer crse_l1, crse_l2, crse_h1, crse_h2
      integer state_l1, state_l2, state_h1, state_h2
      integer fblo(2), fbhi(2)
      integer cblo(2), cbhi(2)
      integer fb_l1, fb_l2, fb_h1, fb_h2
      integer cb_l1, cb_l2, cb_h1, cb_h2
      integer lratiox, lratioy, nvar
      integer bc(2,2,nvar)
      DOUBLE PRECISION fine(fine_l1:fine_h1, fine_l2:fine_h2,nvar)
      DOUBLE PRECISION crse(crse_l1:crse_h1, crse_l2:crse_h2, nvar)
      DOUBLE PRECISION fine_state(state_l1:state_h1, state_l2:state_h2, 
     &nvar)
      DOUBLE PRECISION fvcx(fb_l1:fb_h1)
      DOUBLE PRECISION fvcy(fb_l2:fb_h2)
      DOUBLE PRECISION cvcx(cb_l1:cb_h1)
      DOUBLE PRECISION cvcy(cb_l2:cb_h2)

      integer rMAX
      parameter (rMAX = 32)
      DOUBLE PRECISION alpha, sumN, sumP, negVal, posVal
      DOUBLE PRECISION crseTot, crseTotnew
      DOUBLE PRECISION orig_fine(0:rMAX-1,0:rMAX-1)
      DOUBLE PRECISION fvol,cvol
      integer redo_me
      integer ilo,ihi,jlo,jhi
      integer i,j,ic,jc,n
      integer icase

      if (MAX(lratiox,lratioy).gt.rMAX) then
         print *,'rMAX in INTERP_2D::FORT_PROTECT_INTERP must be >= ',MA
     &X(lratiox,lratioy)
         call bl_abort(" ")
      endif

      do jc = cblo(2), cbhi(2)
      do ic = cblo(1), cbhi(1)

         ilo = max(lratiox*ic            ,fine_l1)
         ihi = min(lratiox*ic+(lratiox-1),fine_h1)
         jlo = max(lratioy*jc            ,fine_l2)
         jhi = min(lratioy*jc+(lratioy-1),fine_h2)

         do n = 2, nvar-1

            redo_me = 0
            do j = jlo,jhi
            do i = ilo,ihi
               if ((fine_state(i,j,n)+fine(i,j,n)) .lt. 0.d0) redo_me = 
     &1
            enddo
            enddo
c
c ****************************************************************************************
c
c           If all the fine values are non-negative after the original interpolated 
c            correction, then we do nothing here.
c
c           If any of the fine values are negative after the original interpolated
c            correction, then we do our best.
c
c           Special cases:
c
c             1) Coarse correction > 0, and fine_state has some cells with 
c                negative values which will be filled before adding to the other cells.
c                Use the correction to bring negative cells to 0.0D0, then
c                distribute the remaining positive proportionally.
c
c             2) Coarse correction > 0, and correction can not make them all
c                positive.  Add correction only to the negative cells, in proportion
c                to their magnitude.
c
c             3) Coarse correction < 0, and fine_state DOES NOT have enough
c                  have enough positive state to absorb it.  Here we bring
c                  all the positive fine cells to 0.0D0 then distribute the remaining
c                  negative amount in such a way as to make them all as close to the
c                  same negative value as possible.
c
c             4) Coarse correction < 0, fine_state has enough
c                  positive state to absorb it without making any fine 
c                  cells negative, BUT fine_state+fine is currently negative
c                  in at least 1.0D0 fine cell.  Here just take a constant percentage
c                  away from each positive and don't touch the negatives.
c
c             crseTot = volume-weighted sum of all interpolated values of the correction,
c                       which is equivalent to the total volume-weighted coarse correction
c             SumN = volume-weighted sum of all negative values of fine_state
c             SumP = volume-weighted sum of all positive values of fine_state
c
c ****************************************************************************************
c

            if (redo_me .eq. 1) then

               icase = 0

               do j = jlo,jhi
               do i = ilo,ihi
                  orig_fine(i-ilo,j-jlo) = fine(i,j,n)
               enddo
               enddo

               crseTot = 0.d0
               do j = jlo,jhi
               do i = ilo,ihi
                  fvol = (fvcx(i+1)-fvcx(i)) * (fvcy(j+1)-fvcy(j))
                  crseTot = crseTot + fvol * fine(i,j,n)
               enddo
               enddo

               cvol = (cvcx(ic+1)-cvcx(ic)) * (cvcy(jc+1)-cvcy(jc))

               sumN = 0.0D0
               sumP = 0.0D0
               do j = jlo,jhi
               do i = ilo,ihi
                  fvol = (fvcx(i+1)-fvcx(i)) * (fvcy(j+1)-fvcy(j))
                  if (fine_state(i,j,n) .le. 0.d0) then
                    sumN = SumN + fvol * fine_state(i,j,n)
                  else
                    sumP = sumP + fvol * fine_state(i,j,n)
                  endif
               enddo
               enddo

               if (crseTot .gt. 0.d0 .and. crseTot .ge. abs(sumN)) then
c              Here we want to fill in the negative values first, then add
c                the remaining positive proportionally.

                   icase = 1
                   do j = jlo,jhi
                   do i = ilo,ihi
                      if (fine_state(i,j,n) .le. 0.d0) then
                        fine(i,j,n) = -fine_state(i,j,n)
                      endif
                   enddo
                   enddo

                   if (sumP > 0.d0) then

                     alpha = (crseTot - abs(sumN)) / sumP

                     do j = jlo,jhi
                     do i = ilo,ihi
                       if (fine_state(i,j,n) .ge. 0.d0) then
                         fine(i,j,n) = alpha * fine_state(i,j,n)
                       endif
                     enddo
                     enddo

                   else

                     posVal = (crseTot - abs(sumN)) / cvol

                     do j = jlo,jhi
                     do i = ilo,ihi
                       fine(i,j,n) = fine(i,j,n) + posVal
                     enddo
                     enddo

                   endif

                 endif

               if (crseTot .gt. 0.d0. and. crseTot .lt. abs(sumN)) then
c              Here we don't have enough positive correction to fill all the
c                negative values of state, so we just try to fill them proportionally
c                and don't add any correction to the states already positive.

                   icase = 2
                   alpha = crseTot / abs(sumN)

                   do j = jlo,jhi
                   do i = ilo,ihi
                     if (fine_state(i,j,n) .lt. 0.d0) then
                       fine(i,j,n) = alpha * abs(fine_state(i,j,n))
                     else 
                       fine(i,j,n) = 0.d0
                     endif
                   enddo
                   enddo

               endif

               if (crseTot .lt. 0.d0. and. abs(crseTot) .gt. sumP) then
c              Here we don't have enough positive states to absorb all the
c                negative correction, so we want to end up with all the fine
c                cells having the same negative value.

                   icase = 3
                   negVal = (sumP + sumN + crseTot)/cvol

                   do j = jlo,jhi
                   do i = ilo,ihi
                      fine(i,j,n) = negVal - fine_state(i,j,n)
                   enddo
                   enddo

               endif

               if (crseTot .lt. 0.d0 .and. abs(crseTot) .lt. sumP.and. (
     &sumP+sumN+crseTot) .gt. 0.d0) then
c              Here we have enough positive states to absorb all the
c                negative correction *and* redistribute to make negative cells
c                positive. 

                   icase = 4
                   alpha = (crseTot + sumN) / sumP

                   do j = jlo,jhi
                   do i = ilo,ihi
                      if (fine_state(i,j,n) .lt. 0.d0) then
                        fine(i,j,n) = -fine_state(i,j,n)
                      else
                        fine(i,j,n) = alpha * fine_state(i,j,n)
                      endif  
                   enddo
                   enddo

               endif

               if (crseTot .lt. 0.d0. and. abs(crseTot) .lt. sumP.and. (
     &sumP+sumN+crseTot) .le. 0.d0) then
c              Here we have enough positive states to absorb all the
c                negative correction, but not to fix the states already negative. 
c                We bring all the positive states to 0.0D0, and use whatever 
c                remaining positiveness from the states to help the negative states.

                   icase = 5
                   alpha = (crseTot + sumP) / sumN

                   do j = jlo,jhi
                   do i = ilo,ihi
                      if (fine_state(i,j,n) .gt. 0.d0) then
                        fine(i,j,n) = -fine_state(i,j,n)
                      else 
                        fine(i,j,n) = alpha * fine_state(i,j,n)
                      endif
                   enddo
                   enddo

               endif

               crseTotnew   = 0.d0
               do j = jlo,jhi
               do i = ilo,ihi
                  fvol = (fvcx(i+1)-fvcx(i)) * (fvcy(j+1)-fvcy(j))
                  crseTotnew   = crseTotnew   + fvol * fine(i,j,n)
               enddo
               enddo

               if (abs(crseTotnew - crseTot)/cvol .gt. 1.e-8) then
                  print *,' '
                  print *,'BLEW CONSERVATION with ICASE = ',icase
                  print *,'AT COARSE CELL ',ic,jc,' AND COMPONENT ',n
                  print *,'CRSETOT NEW OLD ',crseTotnew, crseTot
                  print *,'CVOL ',cvol
                  print *,'SUMP SUMN ',sumP,sumN
                  do j = jlo,jhi
                  do i = ilo,ihi
                     fvol = (fvcx(i+1)-fvcx(i)) * (fvcy(j+1)-fvcy(j))
                     print *,'FINE OLD NEW ',i,j,orig_fine(i-ilo,j-jlo),
     &fine(i,j,n), fine_state(i,j,n),fvol
                     if (abs(fvol) .lt. 1.d-50) then
                       print *,'MAKING FVOL ',fvcx(i+1),fvcx(i),fvcy(j+1
     &),fvcy(j)
                     endif
                  enddo
                  enddo
               endif

c              do j = jlo,jhi
c              do i = ilo,ihi
c                 if ((fine_state(i,j,n) + fine(i,j,n)) .lt. 0.d0) then
c                    print *,'STILL NEGATIVE AT ',i,j,n
c                    print *,'AT COARSE CELL ',ic,jc
c                    print *,'FINE STATE ',fine_state(i,j,n)
c                    print *,'FINE CORRECTION ',fine(i,j,n)
c                    print *,'CRSETOT ',crseTot
c                    print *,'SUMN / SUMP ',sumN, sumP
c                    print *,' '
c                 endif
c              enddo
c              enddo
c              enddo
c           End (if redo .eq. 1)
            endif

         enddo

c     Set sync for density (n=1) to sum of spec sync (2:nvar-1)
         do j = jlo,jhi
         do i = ilo,ihi
            fine(i,j,1) = 0.d0
            do n = 2,nvar-1
               fine(i,j,1) = fine(i,j,1) + fine(i,j,n)
            enddo
         enddo
         enddo

c     End of coarse index loops
      enddo
      enddo
      end

c ::: 
c ::: --------------------------------------------------------------
c ::: quartinterp: quartic conservative interpolation from coarse grid to
c ::: subregion of fine grid defined by (fblo,fbhi)
c ::: 
c ::: Inputs/Outputs
c ::: fine        <=>  (modify) fine grid array
c ::: flo,fhi      =>  (const)  index limits of fine grid
c ::: fblo,fbhi    =>  (const)  subregion of fine grid to get values
c ::: nvar         =>  (const)  number of variables in state vector
c ::: lratio[xy]   =>  (const)  refinement ratio between levels
c ::: 
c ::: crse         =>  (const)  coarse grid data
c ::: clo,chi      =>  (const)  index limits of crse grid
c ::: cblo,cbhi    =>  (const)  coarse grid region containing fblo,fbhi and widen by 2 or 4 cells
c :::
c ::: cb2lo,cb2hi  =>  (const)  coarse grid region containing fblo,fbhi
c ::: fb2lo,fb2hi  =>  (const)  fine version of cb2. It could be wider than fb
c ::: 
c ::: TEMPORARY ARRAYS
c ::: ftmp         =>  1-D temp array
c ::: ctmp         =>  2-D temp array
c ::: --------------------------------------------------------------
c ::: 
       subroutine quartinterp (fine, fine_l1, fine_l2, fine_h1, fine_h2,
     & fblo, fbhi, fb2lo, fb2hi,crse, crse_l1, crse_l2, crse_h1, crse
     &_h2, cblo, cbhi, cb2lo, cb2hi,nvar, lratiox, lratioy,ftmp, ctmp
     &,bc,actual_comp,actual_state)

       implicit none

       integer fine_l1, fine_l2, fine_h1, fine_h2
       integer crse_l1, crse_l2, crse_h1, crse_h2
       integer fblo(2), fbhi(2), fb2lo(2), fb2hi(2)
       integer cblo(2), cbhi(2), cb2lo(2), cb2hi(2)
       integer nvar,lratiox,lratioy
       integer bc(2,2,nvar)
       integer actual_comp,actual_state
       DOUBLE PRECISION fine(fine_l1:fine_h1, fine_l2:fine_h2,nvar)
       DOUBLE PRECISION crse(crse_l1:crse_h1, crse_l2:crse_h2,nvar)
       DOUBLE PRECISION ftmp(fb2lo(1):fb2hi(1))
       DOUBLE PRECISION ctmp(cblo(1):cbhi(1),0:lratioy-1)

c      Local variables
       integer i,j,ii,jj,n,iry
       DOUBLE PRECISION cL(-2:2)
c       DOUBLE PRECISION cR(-2:2)
       data cL/ -0.01171875D0, 0.0859375D0, 0.5d0, -0.0859375D0,0.011718
     &75D0 /
c       data cR/  0.01171875D0, -0.0859375D0, 0.5d0,  0.0859375D0,
c     $          -0.01171875D0 /

       if (lratiox.eq.2 .and. lratioy.eq.2) then
          do n = 1, nvar
             do j = cb2lo(2), cb2hi(2)
                do i = cblo(1), cbhi(1)
                   ctmp(i,0) = 2.d0*(cL(-2)*crse(i,j-2,n) + cL(-1)*crse(
     &i,j-1,n)+ cL( 0)*crse(i,j ,n)+ cL( 1)*crse(i,j+1,n)
     &+ cL( 2)*crse(i,j+2,n))
                   ctmp(i,1) = 2.d0*crse(i,j,n) - ctmp(i,0)
c$$$                   ctmp(i,1) = 2.d0*(cR(-2)*crse(i,j-2,n) 
c$$$     $                  +            cR(-1)*crse(i,j-1,n)
c$$$     $                  +            cR( 0)*crse(i,j  ,n)
c$$$     $                  +            cR( 1)*crse(i,j+1,n)
c$$$     $                  +            cR( 2)*crse(i,j+2,n))
                enddo
                do iry = 0, 1
                   jj = j*2+iry
                   if (jj.ge.fblo(2).and.jj.le.fbhi(2)) then
                      do i = cb2lo(1), cb2hi(1)
                         ii = 2*i
                         ftmp(ii ) = 2.d0*(cL(-2)*ctmp(i-2,iry) + cL(-1)
     &*ctmp(i-1,iry)+ cL( 0)*ctmp(i ,iry)+ cL( 1)*c
     &tmp(i+1,iry)+ cL( 2)*ctmp(i+2,iry))
                         ftmp(ii+1) = 2.d0*ctmp(i,iry) - ftmp(ii)
c$$$                         ftmp(ii+1) = 2.d0*(cR(-2)*ctmp(i-2,iry) 
c$$$     $                        +             cR(-1)*ctmp(i-1,iry)
c$$$     $                        +             cR( 0)*ctmp(i  ,iry)
c$$$     $                        +             cR( 1)*ctmp(i+1,iry)
c$$$     $                        +             cR( 2)*ctmp(i+2,iry))
                      enddo
                      do ii = fblo(1), fbhi(1)
                         fine(ii,jj,n) = ftmp(ii)
                      enddo
                   endif
                enddo
             enddo
          enddo
       else if (lratiox.eq.4 .and. lratioy.eq.4) then
c      todo
          write(6,*) 'FORT_QUARTINTERP: refinement ratio = 4 TODO'
          stop
       else
          write(6,*) 'FORT_QUARTINTERP: unsupported refinement ratio'
          stop
       endif

       end
