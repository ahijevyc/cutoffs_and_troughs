module netcdf_routines_mod

! ifort -L${NETCDF}/lib -lnetcdf -lnetcdff -I${NETCDF}/include average_netcdf_files_parallel.f90 -o average_netcdf_files_parallel.x
! mpif90 -L${NETCDF}/lib -lnetcdf -lnetcdff -I${NETCDF}/include average_netcdf_files_parallel.f90 -o average_netcdf_files_parallel.x

! Need to overwrite variables in file filenameout, which must be present in run directory

use netcdf

implicit none

integer, parameter  :: r_single = selected_real_kind(6)  ! single precision
integer, parameter  :: r_double = selected_real_kind(15) ! double precision
integer, parameter  :: i_byte   = selected_int_kind(1)   ! byte integer
integer, parameter  :: i_short  = selected_int_kind(4)   ! short integer
integer, parameter  :: i_long   = selected_int_kind(8)   ! long integer
integer, parameter  :: i_kind   = i_long                 ! default integer
integer, parameter  :: r_kind   = r_double               ! default real

!!!!!

contains

subroutine open_netcdf(fname,ncfileid)
   character(len=*), intent(in) :: fname
   integer, intent(out) :: ncfileid

   integer :: ncstatus

   ncstatus = nf90_open(path=trim(adjustl(fname)),mode=nf90_nowrite,ncid=ncfileid)  ! open file
   if ( ncstatus .eq. 0 ) then
!     write(*,fmt='(a)') 'opened '//trim(adjustl(fname))//' for reading'
!     write(*,fmt='(a,i8)') 'fileid = ',ncfileid
   else
      write(*,fmt='(a)') 'error reading '//trim(adjustl(fname))
!     call stop2(31) ! stop
   endif

   return
end subroutine open_netcdf

subroutine close_netcdf(fname,ncfileid)
   character(len=*), intent(in) :: fname
   integer, intent(in) :: ncfileid

   integer :: ncstatus

   ncstatus = nf90_close(ncfileid) ! close file
   if ( ncstatus .ne. 0 ) then
      write(*,fmt='(a)') 'error closing '//trim(adjustl(fname))
!     call stop2(31) ! stop
   endif
end subroutine close_netcdf

subroutine get_netcdf_dims( cdfid, var, dims, ndims, do_output )
  implicit none

  integer, intent(in) :: cdfid
  character (len=*), intent(in) :: var
  integer, intent(inout), dimension(4) :: dims
  integer, intent(inout) :: ndims
  logical, intent(in) :: do_output

  integer :: rcode
  character (len=80) :: varname
  integer :: natts, dimids(10)
  integer :: i, ncvarid, xtype

  rcode = nf90_inq_varid( cdfid, trim(adjustl(var)), ncvarid)
  rcode = nf90_Inquire_Variable(cdfid, ncvarid, varname, xtype, ndims, dimids, natts)

  dims(:) = -1

  do i=1,ndims
    rcode = nf90_inquire_dimension( cdfid, dimids(i), len=dims(i) )
    if ( do_output ) write(6,*) ' dimension ',i,dims(i)
  enddo

end subroutine get_netcdf_dims

subroutine get_netcdf_var_1d_real(fileid,variable,dims,output)

   integer, intent(in) :: fileid, dims(4)
   character(len=*), intent(in) :: variable
   real(r_single), intent(inout), dimension(dims(1)) :: output

   integer :: ncstatus, ncvarid, istatus

   istatus = 0
   ncstatus = nf90_inq_varid(fileid,trim(adjustl(variable)),ncvarid)
   if ( ncstatus /= 0 ) istatus = istatus + ncstatus
   ncstatus = nf90_get_var(fileid,ncvarid,output)
   if ( ncstatus /= 0 ) istatus = istatus + ncstatus

   if ( istatus /= 0 ) then
      write(*,*) 'Error reading '//trim(adjustl(variable))
  !   call mpi_finalize(iret)
      stop
   endif

end subroutine get_netcdf_var_1d_real

!subroutine check(status)
!   integer, intent ( in) :: status
!   if(status /= nf90_noerr) then 
!      print *, trim(nf90_strerror(status))
!      stop "Stopped"
!   end if
!end subroutine check 

!------------

end module netcdf_routines_mod
