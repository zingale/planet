

c ---------------------------------------------------------------
c ::  bdintrpxlo : Interpolation on Xlo Face
c ::       Quadratic Interpolation from crse data
c ::       in directions transverse to face of grid
c ::
c ::  Inputs/Outputs:
c ::  bdry       <=  fine grid bndry data strip
c ::  bdry_l1, bdry_l2, bdry_h1, bdry_h2  => index limits of bdry
c ::  lo,hi       => index limits of grd interior
c ::  cb_l1, cb_l2, cb_h1, cb_h2    => index limits of coarsened grid interior
c ::  nvar        => number of variables to interpolate
c ::  ratios(2)   => refinement ratios
c ::  not_covered => mask is set to this value if cell is not
c ::                 covered by another fine grid and not outside the domain.
c ::  mask        => fine grid mask bndry strip
c ::  mask_l1, mask_l2, mask_h1, mask_h2  => index limits of mask array
c ::  crse        => crse grid bndry data strip
c ::  crse_l1, crse_l2, crse_h1, crse_h2  => index limits of crse array
c ::  derives     => crse grid tmp array
c ---------------------------------------------------------------
      subroutine bdintrpxlo (bdry,bdry_l1, bdry_l2, bdry_h1, bdry_h2,lo,
     &hi,cb_l1, cb_l2, cb_h1, cb_h2,nvar,ratios,not_covered,mask,mask_
     &l1, mask_l2, mask_h1, mask_h2,crse,crse_l1, crse_l2, crse_h1, cr
     &se_h2,derives,max_order)

      implicit none
      integer  nvar, ratios(2), not_covered,max_order
      integer  lo(2), hi(2)
      integer  bdry_l1, bdry_l2, bdry_h1, bdry_h2
      integer  mask_l1, mask_l2, mask_h1, mask_h2
      integer  crse_l1, crse_l2, crse_h1, crse_h2
      integer  cb_l1, cb_l2, cb_h1, cb_h2
      DOUBLE PRECISION   bdry(bdry_l1:bdry_h1, bdry_l2:bdry_h2,nvar)
      DOUBLE PRECISION   derives(cb_l2:cb_h2,2)      
      integer  mask(mask_l1:mask_h1, mask_l2:mask_h2)
      DOUBLE PRECISION   crse(crse_l1:crse_h1, crse_l2:crse_h2,nvar)

      integer  i, j, ic, jc, off, n
      integer  jclo, jchi, ratioy

      integer Norder, NN, m
      parameter (Norder = 3)
      DOUBLE PRECISION x(Norder), y(Norder), c(Norder), xInt

      ratioy = ratios(2)

      jclo = cb_l2
      jchi = cb_h2
      ic   = cb_l1-1
      i    = lo(1)-1

      if (max_order.eq.1) then
         do n = 1, nvar
            do off = 0, ratioy - 1
               do jc = jclo, jchi
                  j = ratioy*jc + off
                  bdry(i,j,n) = crse(ic,jc,n)
               end do
            end do
         end do
      else

      do n=1,nvar
         do jc=jclo,jchi
            j = ratioy*jc

            NN = 1
            y(NN) = crse(ic,jc,n)
            x(NN) = 0.0D0

            if (mask(i,j-1).eq.not_covered) then
               NN=NN+1
               y(NN) = crse(ic,jc-1,n)
               x(NN) = -1.0D0
            else if (mask(i,j+2*ratioy).eq.not_covered) then
               NN=NN+1
               y(NN) = crse(ic,jc+2,n)
               x(NN) = 2.0D0
            endif

            if (mask(i,j+ratioy).eq.not_covered) then
               NN=NN+1
               y(NN) = crse(ic,jc+1,n)
               x(NN) = 1.0D0
            else if (mask(i,j-ratioy-1).eq.not_covered) then
               NN=NN+1
               y(NN) = crse(ic,jc-2,n)
               x(NN) = -2.0D0
            endif

            if ( (mask(i,j-1).ne.not_covered).and.(mask(i,j+ratioy).ne.n
     &ot_covered) ) NN = 1

            do off = 0,ratioy-1
               xInt = (dble(off - ratioy/2) + 0.5D0)/ratioy
               call polyInterpCoeff(xInt, x, NN, c)
               bdry(i,j+off,n) = 0.0D0
               do m=1,NN
                  bdry(i,j+off,n) = bdry(i,j+off,n) + c(m)*y(m)
               end do
            end do
         end do
      end do

      endif
      end

c ---------------------------------------------------------------
c ::  bdintrpxhi : Interpolation on Xhi Face
c ::       Quadratic Interpolation from crse data
c ::       in directions transverse to face of grid
c ::
c ::  Inputs/Outputs:
c ::  bdry       <=  fine grid bndry data strip
c ::  bdry_l1, bdry_l2, bdry_h1, bdry_h2  => index limits of bdry
c ::  lo,hi       => index limits of grd interior
c ::  cb_l1, cb_l2, cb_h1, cb_h2    => index limits of coarsened grid interior
c ::  nvar        => number of variables to interpolate
c ::  ratios(2)   => refinement ratios
c ::  not_covered => mask is set to this value if cell is not
c ::                 covered by another fine grid and not outside the domain.
c ::  mask        => fine grid mask bndry strip
c ::  mask_l1, mask_l2, mask_h1, mask_h2  => index limits of mask array
c ::  crse        => crse grid bndry data strip
c ::  crse_l1, crse_l2, crse_h1, crse_h2  => index limits of crse array
c ::  derives     => crse grid tmp array
c ---------------------------------------------------------------
      subroutine bdintrpxhi (bdry,bdry_l1, bdry_l2, bdry_h1, bdry_h2,lo,
     &hi,cb_l1, cb_l2, cb_h1, cb_h2,nvar,ratios,not_covered,mask,mask_
     &l1, mask_l2, mask_h1, mask_h2,crse,crse_l1, crse_l2, crse_h1, cr
     &se_h2,derives,max_order)
      implicit none
      integer  nvar, ratios(2), not_covered,max_order
      integer  lo(2), hi(2)
      integer  bdry_l1, bdry_l2, bdry_h1, bdry_h2
      integer  mask_l1, mask_l2, mask_h1, mask_h2
      integer  cb_l1, cb_l2, cb_h1, cb_h2
      integer  crse_l1, crse_l2, crse_h1, crse_h2
      DOUBLE PRECISION   bdry(bdry_l1:bdry_h1, bdry_l2:bdry_h2,nvar)
      DOUBLE PRECISION   derives(cb_l2:cb_h2,2)      
      integer  mask(mask_l1:mask_h1, mask_l2:mask_h2)
      DOUBLE PRECISION   crse(crse_l1:crse_h1, crse_l2:crse_h2,nvar)

      integer  i, j, ic, jc, off, n
      integer  jclo, jchi, ratioy

      integer Norder, NN, m
      parameter (Norder = 3)
      DOUBLE PRECISION x(Norder), y(Norder), c(Norder), xInt

      ratioy = ratios(2)

      jclo = cb_l2
      jchi = cb_h2
      ic   = cb_h1+1
      i    = hi(1)+1

      if (max_order.eq.1) then
         do n = 1, nvar
            do off = 0, ratioy - 1
               do jc = jclo, jchi
                  j = ratioy*jc + off
                  bdry(i,j,n) = crse(ic,jc,n)
               end do
            end do
         end do
      else

      do n=1,nvar
         do jc=jclo,jchi
            j = ratioy*jc

            NN = 1
            y(NN) = crse(ic,jc,n)
            x(NN) = 0.0D0

            if (mask(i,j-1).eq.not_covered) then
               NN=NN+1
               y(NN) = crse(ic,jc-1,n)
               x(NN) = -1.0D0
            else if (mask(i,j+2*ratioy).eq.not_covered) then
               NN=NN+1
               y(NN) = crse(ic,jc+2,n)
               x(NN) = 2.0D0
            endif

            if (mask(i,j+ratioy).eq.not_covered) then
               NN=NN+1
               y(NN) = crse(ic,jc+1,n)
               x(NN) = 1.0D0
            else if (mask(i,j-ratioy-1).eq.not_covered) then
               NN=NN+1
               y(NN) = crse(ic,jc-2,n)
               x(NN) = -2.0D0
            endif

            if ( (mask(i,j-1).ne.not_covered).and.(mask(i,j+ratioy).ne.n
     &ot_covered) ) NN = 1

            do off = 0,ratioy-1
               xInt = (dble(off - ratioy/2) + 0.5D0)/ratioy
               call polyInterpCoeff(xInt, x, NN, c)
               bdry(i,j+off,n) = 0.0D0
               do m=1,NN
                  bdry(i,j+off,n) = bdry(i,j+off,n) + c(m)*y(m)
               end do
            end do
         end do
      end do

      endif
      end

c ---------------------------------------------------------------
c ::  bdintrpylo : Interpolation on Ylo Face
c ::       Quadratic Interpolation from crse data
c ::       in directions transverse to face of grid
c ::
c ::  Inputs/Outputs:
c ::  bdry       <=  fine grid bndry data strip
c ::  bdry_l1, bdry_l2, bdry_h1, bdry_h2  => index limits of bdry
c ::  lo,hi       => index limits of grd interior
c ::  cb_l1, cb_l2, cb_h1, cb_h2    => index limits of coarsened grid interior
c ::  nvar        => number of variables to interpolate
c ::  ratios(2)   => refinement ratios
c ::  not_covered => mask is set to this value if cell is not
c ::                 covered by another fine grid and not outside the domain.
c ::  mask        => fine grid mask bndry strip
c ::  mask_l1, mask_l2, mask_h1, mask_h2  => index limits of mask array
c ::  crse        => crse grid bndry data strip
c ::  crse_l1, crse_l2, crse_h1, crse_h2  => index limits of crse array
c ::  derives     => crse grid tmp array
c ---------------------------------------------------------------
      subroutine bdintrpylo (bdry,bdry_l1, bdry_l2, bdry_h1, bdry_h2,lo,
     &hi,cb_l1, cb_l2, cb_h1, cb_h2,nvar,ratios,not_covered,mask,mask_
     &l1, mask_l2, mask_h1, mask_h2,crse,crse_l1, crse_l2, crse_h1, cr
     &se_h2,derives,max_order)
      implicit none
      integer  nvar, ratios(2), not_covered,max_order
      integer  lo(2), hi(2)
      integer  bdry_l1, bdry_l2, bdry_h1, bdry_h2
      integer  mask_l1, mask_l2, mask_h1, mask_h2
      integer  cb_l1, cb_l2, cb_h1, cb_h2
      integer  crse_l1, crse_l2, crse_h1, crse_h2
      DOUBLE PRECISION   bdry(bdry_l1:bdry_h1, bdry_l2:bdry_h2,nvar)
      DOUBLE PRECISION   derives(cb_l1:cb_h1,2)
      integer  mask(mask_l1:mask_h1, mask_l2:mask_h2)
      DOUBLE PRECISION   crse(crse_l1:crse_h1, crse_l2:crse_h2,nvar)

      integer  i, j, ic, jc, off, n
      integer  iclo, ichi, ratiox

      integer Norder, NN, m
      parameter (Norder = 3)
      DOUBLE PRECISION x(Norder), y(Norder), c(Norder), xInt

      ratiox = ratios(1)

      iclo = cb_l1
      ichi = cb_h1
      jc   = cb_l2-1
      j    = lo(2)-1

      if (max_order.eq.1) then
      do n = 1, nvar
         do off = 0, ratiox - 1
            do ic = iclo, ichi
               i = ratiox*ic + off
               bdry(i,j,n) = crse(ic,jc,n)
            end do
         end do
      end do
      else

      do n=1,nvar
         do ic=iclo,ichi
            i = ratiox*ic

            NN = 1
            y(NN) = crse(ic,jc,n)
            x(NN) = 0.0D0

            if (mask(i-1,j).eq.not_covered) then
               NN=NN+1
               y(NN) = crse(ic-1,jc,n)
               x(NN) = -1.0D0
            else if (mask(i+2*ratiox,j).eq.not_covered) then
               NN=NN+1
               y(NN) = crse(ic+2,jc,n)
               x(NN) = 2.0D0
            endif

            if (mask(i+ratiox,j).eq.not_covered) then
               NN=NN+1
               y(NN) = crse(ic+1,jc,n)
               x(NN) = 1.0D0
            else if (mask(i-ratiox-1,j).eq.not_covered) then
               NN=NN+1
               y(NN) = crse(ic-2,jc,n)
               x(NN) = -2.0D0
            endif

            if ( (mask(i-1,j).ne.not_covered).and.(mask(i+ratiox,j).ne.n
     &ot_covered) ) NN = 1

            do off = 0,ratiox-1
               xInt = (dble(off - ratiox/2) + 0.5D0)/ratiox
               call polyInterpCoeff(xInt, x, NN, c)
               bdry(i+off,j,n) = 0.0D0
               do m=1,NN
                  bdry(i+off,j,n) = bdry(i+off,j,n) + c(m)*y(m)
               end do
            end do
         end do
      end do

      endif
      end

c ---------------------------------------------------------------
c ::  bdintrpyhi : Interpolation on Yhi Face
c ::       Quadratic Interpolation from crse data
c ::       in directions transverse to face of grid
c ::
c ::  Inputs/Outputs:
c ::  bdry       <=  fine grid bndry data strip
c ::  bdry_l1, bdry_l2, bdry_h1, bdry_h2  => index limits of bdry
c ::  lo,hi       => index limits of grd interior
c ::  cb_l1, cb_l2, cb_h1, cb_h2    => index limits of coarsened grid interior
c ::  nvar        => number of variables to interpolate
c ::  ratios(2)   => refinement ratios
c ::  not_covered => mask is set to this value if cell is not
c ::                 covered by another fine grid and not outside the domain.
c ::  mask        => fine grid mask bndry strip
c ::  mask_l1, mask_l2, mask_h1, mask_h2  => index limits of mask array
c ::  crse        => crse grid bndry data strip
c ::  crse_l1, crse_l2, crse_h1, crse_h2  => index limits of crse array
c ::  derives     => crse grid tmp array
c ---------------------------------------------------------------
      subroutine bdintrpyhi (bdry,bdry_l1, bdry_l2, bdry_h1, bdry_h2,lo,
     &hi,cb_l1, cb_l2, cb_h1, cb_h2,nvar,ratios,not_covered,mask,mask_
     &l1, mask_l2, mask_h1, mask_h2,crse,crse_l1, crse_l2, crse_h1, cr
     &se_h2,derives,max_order)
      implicit none
      integer  nvar, ratios(2), not_covered,max_order
      integer  lo(2), hi(2)
      integer  bdry_l1, bdry_l2, bdry_h1, bdry_h2
      integer  mask_l1, mask_l2, mask_h1, mask_h2
      integer  cb_l1, cb_l2, cb_h1, cb_h2
      integer  crse_l1, crse_l2, crse_h1, crse_h2
      DOUBLE PRECISION   bdry(bdry_l1:bdry_h1, bdry_l2:bdry_h2,nvar)
      DOUBLE PRECISION   derives(cb_l1:cb_h1,2)
      integer  mask(mask_l1:mask_h1, mask_l2:mask_h2)
      DOUBLE PRECISION   crse(crse_l1:crse_h1, crse_l2:crse_h2,nvar)

      integer  i, j, ic, jc, off, n
      integer  iclo, ichi, ratiox

      integer Norder, NN, m
      parameter (Norder = 3)
      DOUBLE PRECISION x(Norder), y(Norder), c(Norder), xInt

      ratiox = ratios(1)

      iclo = cb_l1
      ichi = cb_h1
      jc   = cb_h2+1
      j    = hi(2)+1

      if (max_order.eq.1) then
         do n = 1, nvar
            do off = 0, ratiox - 1
               do ic = iclo, ichi
                  i = ratiox*ic + off
                  bdry(i,j,n) = crse(ic,jc,n)
               end do
            end do
         end do
      else

      do n=1,nvar
         do ic=iclo,ichi
            i = ratiox*ic

            NN = 1
            y(NN) = crse(ic,jc,n)
            x(NN) = 0.0D0

            if (mask(i-1,j).eq.not_covered) then
               NN=NN+1
               y(NN) = crse(ic-1,jc,n)
               x(NN) = -1.0D0
            else if (mask(i+2*ratiox,j).eq.not_covered) then
               NN=NN+1
               y(NN) = crse(ic+2,jc,n)
               x(NN) = 2.0D0
            endif

            if (mask(i+ratiox,j).eq.not_covered) then
               NN=NN+1
               y(NN) = crse(ic+1,jc,n)
               x(NN) = 1.0D0
            else if (mask(i-ratiox-1,j).eq.not_covered) then
               NN=NN+1
               y(NN) = crse(ic-2,jc,n)
               x(NN) = -2.0D0
            endif

            if ( (mask(i-1,j).ne.not_covered).and.(mask(i+ratiox,j).ne.n
     &ot_covered) ) NN = 1

            do off = 0,ratiox-1
               xInt = (dble(off - ratiox/2) + 0.5D0)/ratiox
               call polyInterpCoeff(xInt, x, NN, c)
               bdry(i+off,j,n) = 0.0D0
               do m=1,NN
                  bdry(i+off,j,n) = bdry(i+off,j,n) + c(m)*y(m)
               end do
            end do
         end do
      end do

      endif
      end

