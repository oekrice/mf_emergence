!*******************************************************************************
MODULE output
!*******************************************************************************
! Contains tools for saving and printing data to the correct directories etc. Put separately jsut to stop things being messy
!*******************************************************************************
    USE shared_data
    USE netcdf

    IMPLICIT NONE

!*******************************************************************************

CONTAINS

SUBROUTINE diagnostics(diag_num)
    !Calculates some diagnostics and saves to netcdf file as for the triangle code, which was fairly neat (if i say so myself...). Should make for easy pythonning
    IMPLICIT NONE
    INTEGER:: diag_num, proc_test
    INTEGER:: i,j,k
    LOGICAL:: flag1
    character(len=100):: filename

    integer:: id_1, id_2, id_3, id_4, id_5, id_6, id_7, id_8, ncid, nd_id, nz_id

    real(num), dimension(:,:):: jx0(1:nx,1:ny,1:nz),jy0(1:nx,1:ny,1:nz),jz0(1:nx,1:ny,1:nz) !
    real(num), dimension(:,:):: j0(1:nx,1:ny,1:nz)
    real(num), dimension(:,:):: bx0(1:nx,1:ny,1:nz),by0(1:nx,1:ny,1:nz),bz0(1:nx,1:ny,1:nz) !
    real(num), dimension(:,:):: b0(1:nx,1:ny,1:nz)
    real(num), dimension(:,:):: ex0(1:nx,1:ny,1:nz),ey0(1:nx,1:ny,1:nz),ez0(1:nx,1:ny,1:nz) !
    real(num), dimension(:,:):: e0(1:nx,1:ny,1:nz)
    real(num), dimension(:,:):: lx0(1:nx,1:ny,1:nz),ly0(1:nx,1:ny,1:nz),lz0(1:nx,1:ny,1:nz) !
    real(num), dimension(:,:):: l0(1:nx,1:ny,1:nz)

    real(num), dimension(:,:):: time(0:nprocs-1), oflux(0:nprocs-1)
    real(num), dimension(:,:):: sumj(0:nprocs-1), sume(0:nprocs-1)
    real(num), dimension(:,:):: energy(0:nprocs-1), sumlf(0:nprocs-1)

    !Allocate space for the slice, using the global coordinates
    real(num), dimension(:,:,:):: bz0_global(1:nx_global, 1:ny_global, 1:nz_global)
    real(num), dimension(:,:,:):: bx0_global(1:nx_global, 1:ny_global, 1:nz_global)
    real(num), dimension(:,:,:):: l0_global(1:nx_global, 1:ny_global, 1:nz_global)

    real(num), dimension(:):: bx_slice(1:nz_global), lf_heights(1:nz_global)

    !Allocate diagnostic arrays
    if (diag_num == 0) then
        allocate(diag_time(0:ndiags-1))
        allocate(diag_oflux(0:ndiags-1)); allocate(diag_sumj(0:ndiags-1))
        allocate(diag_sume(0:ndiags-1))
        allocate(diag_avgj(0:ndiags-1)); allocate(diag_energy(0:ndiags-1))
        allocate(diag_maxlorentz(0:ndiags-1)); allocate(diag_avglorentz(0:ndiags-1))
        allocate(diag_nulls(0:ndiags-1)); allocate(diag_lfheights(0:ndiags-1, 1:nz_global))
        diag_oflux = 1e6; diag_sumj = 1e6; diag_avgj = 1e6; diag_energy = 1e6
        diag_maxlorentz = 1e6; diag_avglorentz = 1e6; diag_time = 1e6; diag_nulls = 1e6
    end if

    !CURRENT THINGS
    jx0 = 0.25_num*(jx(1:nx,0:ny-1,0:nz-1) + jx(1:nx,1:ny,0:nz-1) + jx(1:nx,0:ny-1,1:nz)  + jx(1:nx,0:ny-1,0:nz-1) )
    jy0 = 0.25_num*(jy(0:nx-1,1:ny,0:nz-1) + jy(1:nx,1:ny,0:nz-1) + jy(0:nx-1,1:ny,1:nz)  + jy(0:nx-1,1:ny,0:nz-1) )
    jz0 = 0.25_num*(jz(0:nx-1,0:ny-1,1:nz) + jz(1:nx,0:ny-1,1:nz) + jz(0:nx-1,1:ny,1:nz)  + jz(1:nx,1:ny,1:nz) )

    j0 = jx0**2 + jy0**2 + jz0**2

    !MAGNETIC FIELD THINGS
    bx0 = 0.5_num*(bx(0:nx-1,1:ny,1:nz) + bx(1:nx,1:ny,1:nz))
    by0 = 0.5_num*(by(1:nx,0:ny-1,1:nz) + by(1:nx,1:ny,1:nz))
    bz0 = 0.5_num*(bz(1:nx,1:ny,0:nz-1) + bz(1:nx,1:ny,1:nz))

    b0 = bx0**2 + by0**2 + bz0**2

    !ELECTRIC FIELD THINGS
    ex0 = 0.25_num*(ex(1:nx,0:ny-1,0:nz-1) + ex(1:nx,1:ny,0:nz-1) + ex(1:nx,0:ny-1,1:nz)  + ex(1:nx,0:ny-1,0:nz-1) )
    ey0 = 0.25_num*(ey(0:nx-1,1:ny,0:nz-1) + ey(1:nx,1:ny,0:nz-1) + ey(0:nx-1,1:ny,1:nz)  + ey(0:nx-1,1:ny,0:nz-1) )
    ez0 = 0.25_num*(ez(0:nx-1,0:ny-1,1:nz) + ez(1:nx,0:ny-1,1:nz) + ez(0:nx-1,1:ny,1:nz)  + ez(1:nx,1:ny,1:nz) )

    e0 = ex0**2 + ey0**2 + ez0**2

    !LORENTZ FORCE THINGS

    lx0 = jy0*bz0 - jz0*by0
    ly0 = jz0*bx0 - jx0*bz0
    lz0 = jx0*by0 - jy0*bx0

    l0 = lx0**2 + ly0**2 + lz0**2


    !TIME
    time(proc_num) = t
    CALL MPI_REDUCE(time(proc_num)/nprocs, diag_time(diag_num), 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, comm, ierr)

    !OPEN FLUX
    if (z_up < 0) then
    oflux(proc_num) = sum(abs(bz(1:nx,1:ny,nz)))*dx*dy
    else
    oflux(proc_num) = 0.0_num
    end if
    !diag_oflux(diag_num) = diag
    CALL MPI_REDUCE(oflux(proc_num), diag_oflux(diag_num), 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, comm, ierr)

    sumj(proc_num) = sum(sqrt(j0))*dx*dy*dz
    CALL MPI_REDUCE(sumj(proc_num), diag_sumj(diag_num), 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, comm, ierr)

    sume(proc_num) = sum(sqrt(e0))*dx*dy*dz
    CALL MPI_REDUCE(sume(proc_num), diag_sume(diag_num), 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, comm, ierr)

    sumlf(proc_num) = sum(sqrt(l0))*dx*dy*dz
    CALL MPI_REDUCE(sumlf(proc_num), diag_avglorentz(diag_num), 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, comm, ierr)

    diag_avglorentz(diag_num) = diag_avglorentz(diag_num)/volume_global

    diag_avgj(diag_num) = diag_sumj(diag_num)/volume_global

    energy(proc_num) = 0.5_num*sum(b0)*dx*dy*dz
    CALL MPI_REDUCE(energy(proc_num), diag_energy(diag_num), 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, comm, ierr)

    !Want the null height, but this depends on tricksy MPI stuff. Might have to be smart :(
    !MPI all the individual processes into bz0_global, one at a time
    bx0_global = 0.0_num; bz0_global = 0.0_num

    !Fill in the bit from the root process, this is easy.
    !What does MPI gather do?!?!

    !This gives the rank information so don't need to send it or anything. Good good. Why bz?!
    bx0_global(1+nx*x_rank:nx*(x_rank+1), 1+ny*y_rank:ny*(y_rank+1), 1+nz*z_rank:nz*(z_rank+1)) = bx0(1:nx,1:ny,1:nz)

    bz0_global(1+nx*x_rank:nx*(x_rank+1), 1+ny*y_rank:ny*(y_rank+1), 1+nz*z_rank:nz*(z_rank+1)) = bz0(1:nx,1:ny,1:nz)

    l0_global(1+nx*x_rank:nx*(x_rank+1), 1+ny*y_rank:ny*(y_rank+1), 1+nz*z_rank:nz*(z_rank+1)) = l0(1:nx,1:ny,1:nz)

    !Doing this one element at a time is EXTEMELY silly but it is the only thing which seems to work...

    do i = 1, nx_global
    do j = 1, ny_global
    do k = 1, nz_global
        call MPI_ALLREDUCE(bx0_global(i,j,k),bx0_global(i,j,k),1,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
        call MPI_ALLREDUCE(l0_global(i,j,k),l0_global(i,j,k),1,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
    end do
    end do
    end do

    if (proc_num == 0) then
        !DO ROPE HEIGHT
        !Take the slice of relevance
        bx_slice(1:nz_global) = bx0_global(nx_global/2, ny_global/2, 1:nz_global)
        !Centre of the rope will be the point at which (working down) positive will go to negative
        flag1 = .false.; diag_nulls(diag_num) = 0.0_num
        do k = nx_global, 2, -1
            if (bx_slice(k) > maxval(bx_slice)/4 .and. bx_slice(k) > 1e-6) flag1 = .true.

            if (flag1 .and. bx_slice(k-1) < 0.0 .and. bx_slice(k) > 0.0) then
                !Centre is between these two values
                diag_nulls(diag_num) = (zc_global(k)*abs(bx_slice(k-1)) + & 
		zc_global(k-1)*abs(bx_slice(k)))/(abs(bx_slice(k)) + abs(bx_slice(k-1)))
                print*, 'Flux rope found at height', diag_nulls(diag_num)
                exit

            end if
        end do
        !DO LORENTZ FORCE SLICES
        do k = 1, nx_global
            diag_lfheights(diag_num, k) = sum(l0_global(1:nx,1:ny,k))
        end do

    end if

    call MPI_BARRIER(comm, ierr)

    if (proc_num == 0) then
      print*, '______________________________________'
      !print*, 'Total current squared', t, diag_sumj(diag_num)
      print*, 'Time', diag_time(diag_num)
      print*, 'Open Flux', diag_oflux(diag_num)
      print*, 'Total Current', diag_sumj(diag_num)
      print*, 'Total Efield', diag_sume(diag_num)
      !print*, 'bz0 check', sum(bz0_global)

      !print*, 'Average Current', diag_avgj(diag_num)
      !print*, 'Magnetic Energy', diag_energy(diag_num)

    !ADMIN
    if (run_number < 10) then
        write (filename, "(A18,I1,A3)") "./diagnostics/run0", int(run_number), ".nc"
    else
        write (filename, "(A17,I2,A3)") "./diagnostics/run", int(run_number), ".nc"
    end if

    !Write to diagnostics file, using netcdf
    call try(nf90_create(trim(filename), nf90_clobber, ncid))
    call try(nf90_def_dim(ncid, 'ndiags', ndiags, nd_id))  !Make up fake dimensions here
    call try(nf90_def_dim(ncid, 'nz', nz_global, nz_id))  !Make up fake dimensions here


    call try(nf90_def_var(ncid, 'time', nf90_double, (/nd_id/), id_1))
    call try(nf90_def_var(ncid, 'openflux', nf90_double, (/nd_id/), id_2))
    call try(nf90_def_var(ncid, 'sumcurrent', nf90_double, (/nd_id/), id_3))
    call try(nf90_def_var(ncid, 'avgcurrent', nf90_double, (/nd_id/), id_4))
    call try(nf90_def_var(ncid, 'energy', nf90_double, (/nd_id/), id_5))
    call try(nf90_def_var(ncid, 'ropeheight', nf90_double, (/nd_id/), id_6))
    call try(nf90_def_var(ncid, 'avglorentz', nf90_double, (/nd_id/), id_7))
    call try(nf90_def_var(ncid, 'lfheights', nf90_double, (/nd_id, nz_id/), id_8))

    call try(nf90_enddef(ncid))

    call try(nf90_put_var(ncid, id_1, diag_time))
    call try(nf90_put_var(ncid, id_2, diag_oflux))
    call try(nf90_put_var(ncid, id_3, diag_sumj))
    call try(nf90_put_var(ncid, id_4, diag_avgj))
    call try(nf90_put_var(ncid, id_5, diag_energy))
    call try(nf90_put_var(ncid, id_6, diag_nulls))
    call try(nf90_put_var(ncid, id_7, diag_avglorentz))
    call try(nf90_put_var(ncid, id_8, diag_lfheights))

    call try(nf90_close(ncid))
    end if

    call MPI_BARRIER(comm, ierr)


END SUBROUTINE diagnostics

subroutine try(status)
! Catch error in reading netcdf fild.
INTEGER, INTENT(IN):: status

if (status /= NF90_noerr) THEN
    PRINT*,TRIM(ADJUSTL(NF90_STRERROR(status)))
end if

end subroutine try


SUBROUTINE save_snap(snap_num)
    !Exports the magnetic field at this plot_num to an appropriate netcdf file
    IMPLICIT NONE

    CHARACTER(LEN =64):: output_filename
    INTEGER:: snap_num, proc_write
    INTEGER:: ncid, vid
    INTEGER:: xs_id, ys_id, zs_id
    INTEGER:: xc_id, yc_id, zc_id
    INTEGER:: bx_id, by_id, bz_id
    INTEGER:: jx_id, jy_id, jz_id
    INTEGER:: ex_id, ey_id, ez_id

    if (snap_num < 10) then
        write (output_filename, "(A27,A3,I1,A3)") trim(data_directory), "000", snap_num, ".nc"
    else if (snap_num < 100) then
        write (output_filename, "(A27,A2,I2,A3)") trim(data_directory), "00", snap_num, ".nc"
    else if (snap_num < 1000) then
        write (output_filename, "(A27,A1,I3,A3)") trim(data_directory), "0", snap_num, ".nc"
    else if (snap_num < 10000) then
        write (output_filename, "(A27,I4,A3)") trim(data_directory), snap_num, ".nc"
    end if

    if (proc_num == 0) then
    call try(nf90_create(trim(output_filename), nf90_clobber, ncid))

    call try(nf90_def_dim(ncid, 'xs', nx_global+1, xs_id))
    call try(nf90_def_dim(ncid, 'ys', ny_global+1, ys_id))
    call try(nf90_def_dim(ncid, 'zs', nz_global+1, zs_id))

    call try(nf90_def_dim(ncid, 'xc', nx_global, xc_id))
    call try(nf90_def_dim(ncid, 'yc', ny_global, yc_id))
    call try(nf90_def_dim(ncid, 'zc', nz_global, zc_id))

    call try(nf90_def_var(ncid, 'bx', nf90_double, (/xs_id ,yc_id, zc_id/), bx_id))
    call try(nf90_def_var(ncid, 'by', nf90_double, (/xc_id ,ys_id, zc_id/), by_id))
    call try(nf90_def_var(ncid, 'bz', nf90_double, (/xc_id ,yc_id, zs_id/), bz_id))

    call try(nf90_def_var(ncid, 'jx', nf90_double, (/xc_id ,ys_id, zs_id/), jx_id))
    call try(nf90_def_var(ncid, 'jy', nf90_double, (/xs_id ,yc_id, zs_id/), jy_id))
    call try(nf90_def_var(ncid, 'jz', nf90_double, (/xs_id ,ys_id, zc_id/), jz_id))

    call try(nf90_def_var(ncid, 'ex', nf90_double, (/xc_id ,ys_id, zs_id/), ex_id))
    call try(nf90_def_var(ncid, 'ey', nf90_double, (/xs_id ,yc_id, zs_id/), ey_id))
    call try(nf90_def_var(ncid, 'ez', nf90_double, (/xs_id ,ys_id, zc_id/), ez_id))

    call try(nf90_enddef(ncid))
    call try(nf90_close(ncid))

    end if
    call MPI_BARRIER(comm,ierr)

    !Each process writes data in turn

    do proc_write = 0 ,nprocs-1
        call MPI_BARRIER(comm,ierr)

        if (proc_num == proc_write) then
            call try(nf90_open(trim(output_filename), nf90_write, ncid))

            call try(nf90_inq_varid(ncid, 'bx', vid))
            call try(nf90_put_var(ncid, vid, bx(0:nx,1:ny,1:nz), &
            start = (/x_rank*nx+1,y_rank*ny+1,z_rank*nz+1/),count = (/nx+1,ny,nz/)))

            call try(nf90_inq_varid(ncid, 'by', vid))
            call try(nf90_put_var(ncid, vid, by(1:nx,0:ny,1:nz), &
            start = (/x_rank*nx+1,y_rank*ny+1,z_rank*nz+1/),count = (/nx,ny+1,nz/)))

            call try(nf90_inq_varid(ncid, 'bz', vid))
            call try(nf90_put_var(ncid, vid, bz(1:nx,1:ny,0:nz), &
            start = (/x_rank*nx+1,y_rank*ny+1,z_rank*nz+1/),count = (/nx,ny,nz+1/)))

            call try(nf90_inq_varid(ncid, 'jx', vid))
            call try(nf90_put_var(ncid, vid, jx(1:nx,0:ny,0:nz), &
            start = (/x_rank*nx+1,y_rank*ny+1,z_rank*nz+1/),count = (/nx,ny+1,nz+1/)))

            call try(nf90_inq_varid(ncid, 'jy', vid))
            call try(nf90_put_var(ncid, vid, jy(0:nx,1:ny,0:nz), &
            start = (/x_rank*nx+1,y_rank*ny+1,z_rank*nz+1/),count = (/nx+1,ny,nz+1/)))

            call try(nf90_inq_varid(ncid, 'jz', vid))
            call try(nf90_put_var(ncid, vid, jz(0:nx,0:ny,1:nz), &
            start = (/x_rank*nx+1,y_rank*ny+1,z_rank*nz+1/),count = (/nx+1,ny+1,nz/)))

            call try(nf90_inq_varid(ncid, 'ex', vid))
            call try(nf90_put_var(ncid, vid, ex(1:nx,0:ny,0:nz), &
            start = (/x_rank*nx+1,y_rank*ny+1,z_rank*nz+1/),count = (/nx,ny+1,nz+1/)))

            call try(nf90_inq_varid(ncid, 'ey', vid))
            call try(nf90_put_var(ncid, vid, ey(0:nx,1:ny,0:nz), &
            start = (/x_rank*nx+1,y_rank*ny+1,z_rank*nz+1/),count = (/nx+1,ny,nz+1/)))

            call try(nf90_inq_varid(ncid, 'ez', vid))
            call try(nf90_put_var(ncid, vid, ez(0:nx,0:ny,1:nz), &
            start = (/x_rank*nx+1,y_rank*ny+1,z_rank*nz+1/),count = (/nx+1,ny+1,nz/)))

            call try(nf90_close(ncid))

        end if
        call MPI_BARRIER(comm,ierr)

    end do


    call mpi_barrier(comm, ierr)
    if (proc_num == 0) print*, 'Saved snapshot number', snap_num, ' at time', t, 'to file ', output_filename

    return


END SUBROUTINE save_snap


!*******************************************************************************
END MODULE output
!********************************************************************************
