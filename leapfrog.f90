program leapfrog
    use, intrinsic :: iso_fortran_env
    use :: OMP_LIB
    implicit none
    real :: xs(106), ys(106), zs(106), xvs(106), yvs(106), zvs(106), ms(106)
    real :: pos(3,106), poshalf(3,106), vs(3,106), vsold(3,106)
    real :: as(3,106), add(3)
    real :: r
    real(real64) :: deltat, g
    integer :: i, j, n, a, v, p, nbodies
    logical :: isfile
    real(real64) :: start, end
    nbodies = 106
    deltat =  1d0/12
    g =  1.0

    open(1,file = 'initial_pos.txt',status = 'old')
    do i = 1, nbodies
        read(1,*) xs(i),ys(i),zs(i),xvs(i),yvs(i),zvs(i),ms(i)
        pos(1:3,i) = (/xs(i), ys(i), zs(i)/)
        vs(1:3, i) = (/xvs(i),yvs(i),zvs(i)/)
    enddo

    do j = 1, 5
        write(*,*) SQRT(xs(j)**2+ys(j)**2+zs(j)**2)
    enddo

as(1:3,1:nbodies) = 0
!look for textfiles w/ body positions and delete if they exist already
inquire(file="mercpos.txt",exist=isfile)
if (isfile) then
    open(UNIT=2, file="mercpos.txt",status="OLD")
    close(UNIT=2, status="DELETE")
end if

inquire(file="venuspos.txt",exist=isfile)
if (isfile) then
    open(UNIT=2, file="venuspos.txt",status="OLD")
    close(UNIT=2, status="DELETE")
end if


inquire(file="earthpos.txt",exist=isfile)
if (isfile) then
    open(UNIT=2, file="earthpos.txt",status="OLD")
    close(UNIT=2, status="DELETE")
end if


inquire(file="marspos.txt",exist=isfile)
if (isfile) then
    open(UNIT=2, file="marspos.txt",status="OLD")
    close(UNIT=2, status="DELETE")
end if

inquire(file="juppos.txt",exist=isfile)
if (isfile) then
    open(UNIT=2, file="juppos.txt",status="OLD")
    close(UNIT=2, status="DELETE")
end if

!open textfiles to write positions of bodies
!open(3, file='mercpos.txt',status='new')
!open(4, file='venuspos.txt',status='new')
!open(5, file='earthpos.txt',status='new')
!open(6, file='marspos.txt',status='new')
!open(7, file='juppos.txt',status='new')


!do leapfrog
call cpu_time(start)
do n=1, 100000
vsold =  vs
!take half step in position
    do i = 1, nbodies
        poshalf(1:3,i) = pos(1:3,i)+0.5*deltat*vs(1:3,i)
    enddo

!calculate accelerations on each body
    do a = 1, nbodies
        do i = 1, nbodies
        if (a .NE. i) then
            r = SQRT((poshalf(1,a)-poshalf(1,i))**2+(poshalf(2,a)-poshalf(2,i))**2+(poshalf(3,a)-poshalf(3,i))**2)
            add = (-g*ms(i)*(/poshalf(1,a)-poshalf(1,i),poshalf(2,a)-poshalf(2,i),poshalf(3,a)-poshalf(3,i)/) / r**3)
            
            as(1:3,a) = as(1:3,a)+add
        endif
        enddo
    enddo

!advance velocity
!$OMP PARALLEL
!$OMP DO    
do v= 1, nbodies
    vs(1:3,v) = vsold(1:3,v) +deltat*as(1:3,v)
enddo
!$OMP END DO
as(1:3,:) = 0
!advance positions


!$OMP DO
do p = 1, nbodies
    !write(*,*) 'hello from process', OMP_GET_THREAD_NUM()
    pos(1:3,p) = poshalf(1:3,p) + 0.5*deltat*vs(1:3,p)
enddo
!$OMP END DO
!$OMP END PARALLEL
    !write(3,*) pos(:,1), deltat*n, ms(1)
    !write(4,*) pos(:,2), deltat*n, ms(2)
    !write(5,*) pos(:,3), deltat*n, ms(3)
    !write(6,*) pos(:,4), deltat*n, ms(4)
    !write(7,*) pos(:,5), deltat*n, ms(5)
enddo
call cpu_time(end)


write(*,*) 'time elapsed', (end-start)
end program leapfrog