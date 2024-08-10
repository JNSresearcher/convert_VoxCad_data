! This program generates domain images from a geoPHYS array on an unstructured grid in .vtk format for ParaView

! Fortran code created by J.Sochor  ( https://github.com/JNSresearcher) 

subroutine writeVtk ( nPHYS, Ncells, sdx, sdy, sdz, delta, geoPHYS, files)
implicit none

integer,      intent(IN):: nPHYS, &               ! number of physical domains
                           sdx,sdy,sdz            ! number of cells along X,  Y and  Z

integer,      intent(IN):: Ncells(nPHYS)          ! array of number of cells in each physical domain
real(8),      intent(IN):: delta(3)               ! grid spacing along X, Y and Z
character(*), intent(IN):: files                  ! output files name
integer(1),   intent(IN):: geoPHYS(sdx,sdy,sdz)   ! array with physical domain numbers

! --------------working variable ------------ !
character (len = 1)  ci
character (len = 4)  ch_sd
character (len = 30) buf1, buf2
integer:: i,j,k,m, n,numcells,  nim,njm,nkm,nip,njp,nkp, ios, n_new
real(8):: s,sxm,sym,szm
    
! --- start of calculation --- !
  numcells=sum(Ncells)

  ios=0
  open(newunit=n_new,file=trim(files)//'.vtk',  iostat=ios)  
    
  if (ios/=0) then
    print *,'error! Could not open the outfile'//trim(files)//'.vtk';   stop
  end if
    
! HEADER OUTPUT
  write (n_new, '(a )') "# vtk DataFile Version 3.0"//new_line(ci)//"out data result"//new_line(ci)//"ASCII" 
  write (n_new, '(a )' )"DATASET UNSTRUCTURED_GRID"
  write (buf1,'(i8)')   numcells*8 
  write (n_new, '(a )' ) "POINTS "//trim(adjustl(buf1))//" double" 

  ! output of point coordinates
  do m = 1, nPHYS                                     
      do k=1,sdz;    do j=1, sdy;    do i=1, sdx
        if (geoPHYS(i,j,k) == m)  then 
  
            ! point 0 x = i; y=j ; z=k
              sxm = real(i,8)*delta(1) - delta(1)
              sym = real(j,8)*delta(2) - delta(2) 
              szm = real(k,8)*delta(3) - delta(3) 
              write(n_new, * ) sxm, sym, szm
              
            ! point 1 x = i+1; y=j ; z=k
              sxm = real(i+1,8)*delta(1) - delta(1)
              sym = real(j,8)*delta(2) - delta(2) 
              szm = real(k,8)*delta(3) - delta(3) 
              write(n_new, * ) sxm, sym, szm
              
            ! point 2 x = i; y=j+1 ; z=k
              sxm = real(i,8)*delta(1) - delta(1)
              sym = real(j+1,8)*delta(2) - delta(2) 
              szm = real(k,8)*delta(3) - delta(3) 
              write(n_new, * ) sxm, sym, szm
              
            ! point 3 x = i+1; y=j+1 ; z=k
              sxm = real(i+1,8)*delta(1) - delta(1)
              sym = real(j+1,8)*delta(2) - delta(2) 
              szm = real(k,8)*delta(3) - delta(3) 
              write(n_new, * ) sxm, sym, szm
              
            ! point 4 x = i; y=j ; z=k+1
              sxm = real(i,8)*delta(1) - delta(1)
              sym = real(j,8)*delta(2) - delta(2) 
              szm = real(k+1,8)*delta(3) - delta(3) 
              write(n_new, * ) sxm, sym, szm
              
            ! point 5 x = i+1; y=j ; z=k+1
              sxm = real(i+1,8)*delta(1) - delta(1)
              sym = real(j,8)*delta(2) - delta(2) 
              szm = real(k+1,8)*delta(3) - delta(3) 
              write(n_new, * ) sxm, sym, szm
              
            ! point 6 x = i; y=j+1 ; z=k+1
              sxm = real(i,8)*delta(1) - delta(1)
              sym = real(j+1,8)*delta(2) - delta(2) 
              szm = real(k+1,8)*delta(3) - delta(3) 
              write(n_new, * ) sxm, sym, szm
              
            ! point 7 x = i+1; y=j+1 ; z=k+1
              sxm = real(i+1,8)*delta(1) - delta(1)
              sym = real(j+1,8)*delta(2) - delta(2) 
              szm = real(k+1,8)*delta(3) - delta(3) 
              write(n_new, * ) sxm, sym, szm
          endif
    enddo; enddo; enddo
  enddo
  write(n_new,*)  

  write(buf1,'(i8)')   numcells
  write(buf2,'(i8)')   numcells*9
  write(n_new, '(a )' ) "CELLS "//trim(adjustl(buf1))//" "//trim(adjustl(buf2))
  do i = 0, numcells - 1
      write(n_new, *) 8, 8*i+[0, 1, 2, 3, 4, 5, 6, 7]
  enddo
  write(n_new, *)  

  write(n_new, '(a )' ) "CELL_TYPES "//trim(adjustl(buf1))
  do i=1,numcells
      write(n_new, *) 11
  enddo
  write(n_new,*)

  !   SCALAR FIELD WRITE
  write(n_new, '(a )' ) "CELL_DATA "//trim(adjustl(buf1))
  write(n_new, '(a )' ) "SCALARS "//'Scalar_field '//"long" //" 1"//new_line(ci)//"LOOKUP_TABLE default"
    ! write(n_new)"SCALARS "//'Scalar_field '//"float" //" 1"//new_line(ci)//"LOOKUP_TABLE default"
    
  do i=1,nPHYS
    write(n_new, *) (i, m=1, Ncells(i))
  enddo
   
  close(n_new)

end subroutine writeVtk


subroutine writeVtk_bin ( nPHYS, Ncells, sdx, sdy, sdz, delta, geoPHYS, files)
implicit none

integer,      intent(IN):: nPHYS, &              ! number of physical domains
                           sdx,sdy,sdz           ! number of cells along X,  Y and  Z

integer,      intent(IN):: Ncells(nPHYS)         ! array of number of cells in each physical domain
real(8),      intent(IN):: delta(3)              ! grid spacing along X, Y and Z
character(*), intent(IN):: files                 ! output file name
integer(1),   intent(IN):: geoPHYS(sdx,sdy,sdz)  ! array with physical domain numbers

! working variable 
character (len = 1)  ci
character (len = 4)  ch_sd
character (len = 30) buf1, buf2
integer:: i,j,k,m, n,numcells,  nim,njm,nkm,nip,njp,nkp, ios, n_new
real(8):: s,sxm,sym,szm
    
! --- start of calculation --- !
  numcells=sum(Ncells)

  ios=0
  open(newunit=n_new,file=trim(files)//'_bin.vtk',form="UNFORMATTED",&
  access="STREAM",  convert="big_endian", iostat=ios)  

    
  if (ios/=0) then
    print *,'fatal error! Could not open the outfile result.';   stop
  end if
    
  ! HEADER OUTPUT
  write (n_new) "# vtk DataFile Version 3.0"//new_line(ci)//"out data result"//new_line(ci)//"BINARY" //new_line(ci)
  write(n_new)"DATASET UNSTRUCTURED_GRID"//new_line(ci)
  write(buf1,'(i8)')   numcells*8  !
  write(n_new) "POINTS "//trim(adjustl(buf1))//" double" //new_line(ci)

  ! output of point coordinates
  do m = 1, nPHYS
      do k=1,sdz;    do j=1, sdy;    do i=1, sdx
        if (geoPHYS(i,j,k) == m)  then 
  
          ! point 0 x = i; y=j ; z=k
            sxm = real(i,8)*delta(1) - delta(1)
            sym = real(j,8)*delta(2) - delta(2) 
            szm = real(k,8)*delta(3) - delta(3) 
            write(n_new ) sxm, sym, szm
            
          ! point 1 x = i+1; y=j ; z=k
            sxm = real(i+1,8)*delta(1) - delta(1)
            sym = real(j,8)*delta(2) - delta(2) 
            szm = real(k,8)*delta(3) - delta(3) 
            write(n_new) sxm, sym, szm
            
          ! point 2 x = i; y=j+1 ; z=k
            sxm = real(i,8)*delta(1) - delta(1)
            sym = real(j+1,8)*delta(2) - delta(2) 
            szm = real(k,8)*delta(3) - delta(3) 
            write(n_new) sxm, sym, szm
            
          ! point 3 x = i+1; y=j+1 ; z=k
            sxm = real(i+1,8)*delta(1) - delta(1)
            sym = real(j+1,8)*delta(2) - delta(2) 
            szm = real(k,8)*delta(3) - delta(3) 
            write(n_new) sxm, sym, szm
            
          ! point 4 x = i; y=j ; z=k+1
            sxm = real(i,8)*delta(1) - delta(1)
            sym = real(j,8)*delta(2) - delta(2) 
            szm = real(k+1,8)*delta(3) - delta(3) 
            write(n_new ) sxm, sym, szm
            
          ! point 5 x = i+1; y=j ; z=k+1
            sxm = real(i+1,8)*delta(1) - delta(1)
            sym = real(j,8)*delta(2) - delta(2) 
            szm = real(k+1,8)*delta(3) - delta(3) 
            write(n_new ) sxm, sym, szm
            
          ! point 6 x = i; y=j+1 ; z=k+1
            sxm = real(i,8)*delta(1) - delta(1)
            sym = real(j+1,8)*delta(2) - delta(2) 
            szm = real(k+1,8)*delta(3) - delta(3) 
            write(n_new) sxm, sym, szm
            
          ! point 7 x = i+1; y=j+1 ; z=k+1
            sxm = real(i+1,8)*delta(1) - delta(1)
            sym = real(j+1,8)*delta(2) - delta(2) 
            szm = real(k+1,8)*delta(3) - delta(3) 
            write(n_new) sxm, sym, szm
          endif
      enddo; enddo; enddo
  enddo
  write(n_new)  new_line(ci)

  write(buf1,'(i8)')   numcells  !
  write(buf2,'(i8)')   numcells*9  !
  write(n_new) "CELLS "//trim(adjustl(buf1))//" "//trim(adjustl(buf2))//new_line(ci)
  do i=0,numcells-1
    write(n_new) 8, 8*i+[0, 1, 2, 3, 4, 5, 6, 7]
  enddo
  write(n_new)   new_line(ci)

  write(n_new ) "CELL_TYPES "//trim(adjustl(buf1))//new_line(ci)
  do i=1,numcells
    write(n_new) 11
  enddo
  write(n_new)  new_line(ci)

   !       SCALAR FIELD WRITE
  write(n_new) "CELL_DATA "//trim(adjustl(buf1))//new_line(ci)
  write(n_new)"SCALARS "//'Scalar_field '//"long" //" 1"//new_line(ci)//"LOOKUP_TABLE default"//new_line(ci)
    ! write(n_new)"SCALARS "//'Scalar_field '//"float" //" 1"//new_line(ci)//"LOOKUP_TABLE default"//new_line(ci)
    
  do i=1,nPHYS
    write(n_new) (i, m=1, Ncells(i))
  enddo
  write(n_new)  new_line(ci)
   
  close(n_new)

end subroutine writeVtk_bin
