

c :: ----------------------------------------------------------
c :: SETVOL
c ::             Compute the volume of each cell
c ::
c :: INPUTS / OUTPUTS:
c ::  vol         <=  volume array
c ::  vlo,vhi      => index limits of vol array
c ::  offset       => shift to origin of computational domain
c ::  dx           => cell size
c ::  coord        => coordinate flag (0 = cartesian, 1 = RZ, 2 = RTHETA)
c :: ----------------------------------------------------------
c ::
       subroutine setvol(vol,vol_l1, vol_l2, vol_h1, vol_h2,offset,dx,co
     &ord)
       implicit none
       integer    vol_l1, vol_l2, vol_h1, vol_h2
       integer    coord
       DOUBLE PRECISION     dx(2), offset(2)
       DOUBLE PRECISION     vol(vol_l1:vol_h1, vol_l2:vol_h2)

       integer    i, j
       DOUBLE PRECISION     ri, ro, pi, po, v
       DOUBLE PRECISION     RZFACTOR
       parameter (RZFACTOR = 2.d0*3.14159265358979323846d0)

       if (coord .eq. 0) then
c
c         ::::: cartesian
c
          v = dx(1)*dx(2)
          do j = vol_l2, vol_h2
             do i = vol_l1, vol_h1
                vol(i,j) = v
             end do
          end do
       elseif(coord .eq. 1) then
c
c         ::::: R-Z
c
          do i = vol_l1, vol_h1
             ri = offset(1) + dx(1)*i
             ro = ri + dx(1)
             v = (0.5D0*RZFACTOR)*dx(2)*dx(1)*(ro + ri)
             do j = vol_l2, vol_h2
                vol(i,j) = abs(v)
             end do
          end do
       elseif(coord .eq. 2) then
c
c  	  ::::: R-THETA
c
          do i = vol_l1, vol_h1
             ri = offset(1) + dx(1)*i
             ro = ri + dx(1)
             do j = vol_l2, vol_h2
                pi = offset(2) + dx(2)*j
                po = pi + dx(2)
                v = RZFACTOR*(ro**3 - ri**3)*(cos(pi)-cos(po))/3.0D0
                vol(i,j) = abs(v)
             enddo
          enddo

       end if

       end

c========================================================
c========================================================
        subroutine setvolpt(vol,ro, ri, po, pi, dx, coord)
        implicit none
        integer coord
        DOUBLE PRECISION dx(2)
        DOUBLE PRECISION     vol
        DOUBLE PRECISION     ro, po, pi
        DOUBLE PRECISION     ri

        DOUBLE PRECISION     RZFACTOR
        parameter (RZFACTOR = 2*3.1415926535897932D0)

        if(coord .eq. 0) then
           vol = (ro-ri)*dx(2)
        elseif(coord .eq. 1) then
           vol = 0.5D0*RZFACTOR*dx(2)*(ro**2 - ri**2)
           vol = abs(vol)
        elseif(coord .eq. 2) then
           vol = RZFACTOR*(ro**3-ri**3)*(cos(pi)-cos(po))/3.0D0
        else
            call bl_abort('bogus value of coord ... bndrylib::SETVOLPT')
        endif

        return
        end

c :: ----------------------------------------------------------
c :: SETDLOGA
c ::             Compute  d(log(A))/dr in each cell
c ::
c :: INPUTS / OUTPUTS:
c ::  dloga        <=  dloga array
c ::  dlo,dhi      => index limits of dloga array
c ::  offset       => shift to origin of computational domain
c ::  dx           => cell size
c ::  coord        => coordinate flag (0 = cartesian, 1 = RZ)
c :: ----------------------------------------------------------
c ::
       subroutine setdloga(dloga,dloga_l1, dloga_l2, dloga_h1, dloga_h2,
     &offset,dx,dir,coord)
       implicit none
       integer    dloga_l1, dloga_l2, dloga_h1, dloga_h2
       integer    coord
       DOUBLE PRECISION     dx(2), offset(2)
       DOUBLE PRECISION     dloga(dloga_l1:dloga_h1, dloga_l2:dloga_h2)
       integer dir

       integer    i, j
       DOUBLE PRECISION     rc, dlga, po, pi

       if (coord .eq. 0) then

          do j = dloga_l2, dloga_h2
             do i = dloga_l1, dloga_h1
                dloga(i,j) = 0.0D0
             end do
          end do

       else if( coord .eq. 1 ) then

          if (dir .eq. 0) then
             do i = dloga_l1, dloga_h1
                rc = offset(1) + dx(1)*(dble(i)+0.5d0)
                dlga = 1.d0/rc
                do j = dloga_l2, dloga_h2
                   dloga(i,j) = dlga
                end do
             end do
          else if (dir .eq. 1) then
             do i = dloga_l1, dloga_h1
                do j = dloga_l2, dloga_h2
                   dloga(i,j) = 0.0D0
                end do
             end do
          else
             call bl_abort('setdloga: illegal direction')
          end if

       else if( coord .eq. 2) then
          if (dir .eq. 0) then
             do i = dloga_l1, dloga_h1
                rc = offset(1) + dx(1)*(dble(i)+0.5d0)
                dlga = 2.d0/rc
                do j = dloga_l2, dloga_h2
                   dloga(i,j) = dlga
                enddo
             enddo
          else if (dir .eq. 1) then
             do i = dloga_l1, dloga_h1
                rc = offset(1) + dx(1)*(dble(i)+0.5d0)
                dlga = 1.d0/rc
                do j = dloga_l2, dloga_h2
                   pi = offset(2) + dx(2)*j
                   po = pi + dx(2)
                   dloga(i,j) = dlga/tan(0.5D0*(pi+po))
                enddo
             enddo
          else
          call bl_abort('setdloga: illegal coordinate system')
          endif
       end if

       end

c :: ----------------------------------------------------------
c :: SETAREA
c ::             Compute the area of given cell face
c ::
c :: INPUTS / OUTPUTS:
c ::  area        <=  area array
c ::  alo,ahi      => index limits of area array
c ::  offset       => shift to origin of computational domain
c ::  dx           => cell size
c ::  coord        => coordinate flag (0 =cartesian, 1 = RZ)
c :: ----------------------------------------------------------
c ::
       subroutine setarea(area,area_l1, area_l2, area_h1, area_h2,offset
     &,dx,dir,coord)
       implicit none
       integer    area_l1, area_l2, area_h1, area_h2
       integer    coord, dir
       DOUBLE PRECISION     dx(2), offset(2)
       DOUBLE PRECISION     area(area_l1:area_h1, area_l2:area_h2)

       integer    i, j
       DOUBLE PRECISION     rc, ri, ro, a, pi, po
       DOUBLE PRECISION     RZFACTOR
       parameter (RZFACTOR = 2.d0*3.14159265358979323846d0)

       if (coord .eq. 0) then
c
c         ::::: cartesian
c
          if (dir .eq. 0) then
             do j = area_l2, area_h2
                do i = area_l1, area_h1
                   area(i,j) = dx(2)
                end do
             end do
          else
             do j = area_l2, area_h2
                do i = area_l1, area_h1
                   area(i,j) = dx(1)
                end do
             end do
          end if

       else if (coord .eq. 1) then
c
c         ::::: R-Z
c
          if (dir .eq. 0) then
             do i = area_l1, area_h1
                ri = offset(1) + dx(1)*i
                a = abs(RZFACTOR*ri*dx(2))
                do j = area_l2, area_h2
                   area(i,j) = a
                end do
             end do
          else
             do i = area_l1, area_h1
                rc = offset(1) + dx(1)*(dble(i)+0.5d0)
                a = abs(dx(1)*RZFACTOR*rc)
                do j = area_l2, area_h2
                   area(i,j) = a
                end do
             end do
          end if

       elseif(coord .eq. 2) then
              if (dir .eq. 0) then
                 do i = area_l1, area_h1
                    ri = offset(1) + dx(1)*i
                    do j = area_l2, area_h2
                       pi = offset(2) + dx(2)*j
                       po = pi + dx(2)
                       a = RZFACTOR*ri*ri*(cos(pi)-cos(po))
                       area(i,j) = abs(a)
                    enddo
                 enddo
              elseif(dir .eq. 1) then
                 do i = area_l1, area_h1
                    ri = offset(1) + dx(1)*i
                    ro = ri + dx(1)
                    do j = area_l2, area_h2
                       pi = offset(2) + dx(2)*j
                       a = RZFACTOR*sin(pi)*(ro**2 - ri**2)/2.0D0
                       area(i,j) = abs(a)
                    enddo
                 enddo
              else
                 write(6,*)' bogus dir ', dir
                 call bl_abort(" ")
              endif

       end if

       end

