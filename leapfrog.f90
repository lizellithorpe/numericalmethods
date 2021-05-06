program leapfrog
    use, intrinsic :: iso_fortran_env
    use :: OMP_LIB
    implicit none
    real :: xs(106), ys(106), zs(106), xvs(106), yvs(106), zvs(106), ms(106)
    real :: pos(3,106), poshalf(3,106), vs(3,106), vsold(3,106)
    real :: as(3,106), add(3)
    real :: r
    real(real64) :: deltat, g
    integer :: i, j, n, a, v, p, nbodies, k, tsteps
    logical :: isfile
    real(real64) :: start, end, vel_end, vel_start, pos_end, pos_start
    character(100) :: fname
    nbodies = 6
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
do k=1,nbodies
    write(fname, fmt='(a,i1,a)') 'data_',k, '.txt'
    
    open(unit=1,file=fname, status='new')
    write(1,*) 0d0, 0d0, 0d0
    close(unit=1)
end do
!open textfiles to write positions of bodies
tsteps=10000

!do leapfrog
call cpu_time(start)
do n=1, tsteps
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
  
!advance velocity
call cpu_time(vel_start)
!$OMP PARALLEL
!$OMP DO    
    do v= 1, nbodies
        vs(1:3,v) = vsold(1:3,v) +deltat*as(1:3,v)
    enddo
!$OMP END DO
 call cpu_time(vel_end)
    as(1:3,:) = 0
    !advance positions
    
 call cpu_time(pos_start)   
!$OMP DO
    do p = 1, nbodies
        !write(*,*) 'hello from process', OMP_GET_THREAD_NUM()
        pos(1:3,p) = poshalf(1:3,p) + 0.5*deltat*vs(1:3,p)
    enddo
!$OMP END DO
call cpu_time(pos_end)
!$OMP END PARALLEL
if (MODULO(n,100) .EQ. 0) then 
    do k=1,nbodies
    write(fname, fmt='(a,i1,a)') 'data_',k, '.txt'
    
    open(unit=1,file=fname, status='old',position='append')
    write(1,*) pos(:,k), deltat*n, ms(k)
    close(unit=1)
    end do
endif 

enddo
call cpu_time(end)


write(*,*) 'time elapsed total:', (end-start)
write(*,*) 'time elapsed velocity update:', (vel_end-vel_start)
write(*,*) 'time elapses position update:', (pos_end-pos_start)
end program leapfrog