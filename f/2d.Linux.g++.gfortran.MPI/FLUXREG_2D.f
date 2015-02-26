

c ::: -----------------------------------------------------------
c ::: Add fine grid flux to flux register.  Flux array is a fine grid
c ::: edge based object, Register is a coarse grid edge based object.
c ::: It is assumed that the coarsened flux region contains the register
c ::: region.
c :::
c ::: INPUTS/OUTPUTS:
c ::: reg       <=> edge centered coarse grid flux register
c ::: reg_l1, reg_l2, reg_h1, reg_h2  => index limits for reg
c ::: flx        => edge centered fine grid flux array
c ::: flx_l1, flx_l2, flx_h1, flx_h2  => index limits for flx
c ::: numcomp    => number of components to update
c ::: dir        => direction normal to flux register
c ::: ratio(2)   => refinement ratios between coarse and fine
c ::: mult       => scalar multiplicative factor
c ::: -----------------------------------------------------------
      subroutine frfineadd(reg,reg_l1, reg_l2, reg_h1, reg_h2,flx,flx_l1
     &, flx_l2, flx_h1, flx_h2,numcomp,dir,ratio,mult)

      implicit none

      integer    reg_l1, reg_l2, reg_h1, reg_h2
      integer    flx_l1, flx_l2, flx_h1, flx_h2
      integer    ratio(2), dir, numcomp
      DOUBLE PRECISION     mult
      DOUBLE PRECISION     reg(reg_l1:reg_h1, reg_l2:reg_h2,numcomp)
      DOUBLE PRECISION     flx(flx_l1:flx_h1, flx_l2:flx_h2,numcomp)

      integer    n, i, j, ic, jc, off
      integer    ratiox, ratioy

      ratiox = ratio(1)
      ratioy = ratio(2)

      if (dir .eq. 0) then
c        ::::: flux normal to X direction
         ic = reg_l1
         i = ic*ratiox
         if (reg_l1 .ne. reg_h1) then
            call bl_abort("FORT_FRFINEADD: bad register direction")
         end if
         if (i .lt. flx_l1 .or. i .gt. flx_h1) then
            call bl_abort("FORT_FRFINEADD: index outside flux range")
         end if
         do n = 1, numcomp
            do off = 0, ratioy-1
               do jc = reg_l2, reg_h2
                  j = ratioy*jc + off
                  reg(ic,jc,n) = reg(ic,jc,n) + mult*flx(i,j,n)
               end do
            end do
         end do
      else
c        ::::: flux normal to Y direction
         jc = reg_l2
         j = jc*ratioy
         if (reg_l2 .ne. reg_h2) then
            call bl_abort("FORT_FRFINEADD: bad register direction")
         end if
         if (j .lt. flx_l2 .or. j .gt. flx_h2) then
            call bl_abort("FORT_FRFINEADD: index outside flux range")
         end if
         do n = 1, numcomp
            do off = 0, ratiox-1
               do ic = reg_l1, reg_h1
                  i = ratiox*ic + off
                  reg(ic,jc,n) = reg(ic,jc,n) + mult*flx(i,j,n)
               end do
            end do
         end do
      end if

      end

c ::: -----------------------------------------------------------
c ::: Add fine grid flux times area to flux register.
c ::: Flux array is a fine grid edge based object, Register is a
c ::: coarse grid edge based object.
c ::: It is assumed that the coarsened flux region contains the register
c ::: region.
c :::
c ::: INPUTS/OUTPUTS:
c ::: reg       <=> edge centered coarse grid flux register
c ::: rlo,rhi    => index limits for reg
c ::: flx        => edge centered fine grid flux array
c ::: flx_l1, flx_l2, flx_h1, flx_h2  => index limits for flx
c ::: area       => edge centered area array
c ::: area_l1, area_l2, area_h1, area_h2 => index limits for area
c ::: numcomp    => number of components to update
c ::: dir        => direction normal to flux register
c ::: ratio(2)   => refinements ratio between coarse and fine
c ::: mult       => scalar multiplicative factor
c ::: -----------------------------------------------------------
      subroutine frfaadd(reg,reg_l1, reg_l2, reg_h1, reg_h2,flx,flx_l1, 
     &flx_l2, flx_h1, flx_h2,area,area_l1, area_l2, area_h1, area_h2,n
     &umcomp,dir,ratio,mult)

      implicit none

      integer    reg_l1, reg_l2, reg_h1, reg_h2
      integer    flx_l1, flx_l2, flx_h1, flx_h2
      integer    area_l1, area_l2, area_h1, area_h2
      integer    ratio(2), dir, numcomp
      DOUBLE PRECISION     mult
      DOUBLE PRECISION     reg(reg_l1:reg_h1, reg_l2:reg_h2,numcomp)
      DOUBLE PRECISION     flx(flx_l1:flx_h1, flx_l2:flx_h2,numcomp)
      DOUBLE PRECISION     area(area_l1:area_h1, area_l2:area_h2)

      integer    n, i, j, ic, jc, off
      integer    ratiox, ratioy

      ratiox = ratio(1)
      ratioy = ratio(2)

      if (dir .eq. 0) then
c        ::::: flux normal to X direction
         ic = reg_l1
         i = ic*ratiox
         if (reg_l1 .ne. reg_h1) then
            call bl_abort("FORT_FRFAADD: bad register direction")
         end if
         if (i .lt. flx_l1 .or. i .gt. flx_h1) then
            call bl_abort("FORT_FRFAADD: index outside flux range")
         end if
         do n = 1, numcomp
            do off = 0, ratioy-1
               do jc = reg_l2, reg_h2
                  j = ratioy*jc + off
                  reg(ic,jc,n) = reg(ic,jc,n) + mult*area(i,j)*flx(i,j,n
     &)
               end do
            end do
         end do
      else
c        ::::: flux normal to Y direction
         jc = reg_l2
         j = jc*ratioy
         if (reg_l2 .ne. reg_h2) then
            call bl_abort("FORT_FRFAADD: bad register direction")
         end if
         if (j .lt. flx_l2 .or. j .gt. flx_h2) then
            call bl_abort("FORT_FRFAADD: index outside flux range")
         end if
         do n = 1, numcomp
            do off = 0, ratiox-1
               do ic = reg_l1, reg_h1
                  i = ratiox*ic + off
                  reg(ic,jc,n) = reg(ic,jc,n) + mult*area(i,j)*flx(i,j,n
     &)
               end do
            end do
         end do
      end if

      end

c ::
c :: --------------------------------------------------------------
c :: reflux:   reflux the data on the outer boundary of
c ::           a fine grid.
c ::
c :: Inputs/Outputs
c :: s           <=>  state data array
c :: slo,shi      =>  index limits of s array
c :: vol          =>  volume array
c :: vlo,vhi      =>  index limits of vol array
c :: reg          =>  flux register
c :: rlo,rhi      =>  index limits of reg array
c :: lo,hi        =>  subregion of s array to be updated
c :: numcomp      =>  number of components to update
c :: mult         =>  multiplative factor (+1 or -1 depending on nomal)
c :: --------------------------------------------------------------
c ::
      subroutine frreflux (s,s_l1, s_l2, s_h1, s_h2,vol,vol_l1, vol_l2, 
     &vol_h1, vol_h2,reg,reg_l1, reg_l2, reg_h1, reg_h2,lo,hi,shft,num
     &comp,mult)

      implicit none

      integer    s_l1, s_l2, s_h1, s_h2
      integer    vol_l1, vol_l2, vol_h1, vol_h2
      integer    reg_l1, reg_l2, reg_h1, reg_h2
      integer    lo(2), hi(2), shft(2)
      integer    numcomp
      DOUBLE PRECISION     mult
      DOUBLE PRECISION     reg(reg_l1:reg_h1, reg_l2:reg_h2,numcomp)
      DOUBLE PRECISION       s(s_l1:s_h1, s_l2:s_h2,numcomp)
      DOUBLE PRECISION     vol(vol_l1:vol_h1, vol_l2:vol_h2)

      integer n, i, j, ilo, jlo

      ilo = shft(1)
      jlo = shft(2)

      do n = 1, numcomp
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               s(i-ilo,j-jlo,n) = s(i-ilo,j-jlo,n) + mult*reg(i,j,n)/vol
     &(i-ilo,j-jlo)
            end do
         end do
      end do

      end
