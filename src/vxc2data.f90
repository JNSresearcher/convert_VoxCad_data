! This is a program for converting VoxCad format into text data.

! Fortran code created by J.Sochor  ( https://github.com/JNSresearcher ) 

program vxc2data
implicit none

! ============= INPUT DATA =========
! input data is contained in a file in.vxc

! ============= OUTPUT DATA =========
! Grid spacing along X, Y and Z
REAL(8) :: delta(3) 

! Number of cells along X,  Y and  Z
integer ::  sdx, sdy,  sdz

! Total number of space cells:  Cells = sdx*sdy*sdz
integer :: Cells

integer :: nPHYS,     &   ! number of physical domains
           nPHYS_env, &   ! number of environments domains
           nPHYS_glob     ! total number of domains

! array with physical domain numbers
integer(1), allocatable  :: geoPHYS(:,:,:)

! array of number of cells in each physical domain
integer,   allocatable  :: Ncells(:)

! physical domain parameters
REAL(8),  allocatable  :: valPHYS(:)

! physical domain type and name
character (LEN=6), allocatable  :: typPHYS(:), namePHYS(:)

!output file name
character(:), allocatable ::  files


!============ WORKING VARIABLE ==========
real(8) delta0  ! Lattice Dimension
integer(1), allocatable::    v(:) 
character (len = 10) ch_e
character (len = 50) words(20)
character (len = 74) letter
character(len=:),allocatable :: st, ch

integer m1, m2,  i,j,k,  L, m, n, k0, k1, k2, ios
integer lst, neww, kp, ier, idev

!  when reading data the variable "uncompress" is set to .true if the data is zlib compressed
logical uncompress

! function to convert string to number taking into account prefix
REAL(8)  numeric
external numeric

files = '' ! name for .vtk and .txt

! the position number of the character string "letter" corresponds to the domain number
letter='123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz'

open(1,file='in.vxc')

!  ---------- start data entry -----------!

k0=0;k1=0;k2=0;
m=0 !  m - it is counter CDATA - not USE in do while convert_strings0

convert_strings0: &
do while ( 1==1)

    idev = 1
    call readline(st,ch,ier,idev)             !============================= 1 idev= 1

    if(ier == -1) exit convert_strings0
    
    if (st(1:6)=='</VXC>' ) exit
    lst = len_trim(st)
    
    m1=index(st,'<Lattice_Dim>')
    if (m1>=1) then
        m2=index(st,'</Lattice_Dim>')
        ch_e=trim(st(m1+13:m2-1))
        delta0=numeric(ch_e)
        ! write (*, '( a,g10.3 )') 'Lattice Dimension=',delta0
    endif
    
    m1=index(st,'<X_Dim_Adj>')
    if (m1>=1) then
        m2=index(st,'</X_Dim_Adj>')
        ch_e=trim(st(m1+11:m2-1))
        delta(1)=numeric(ch_e) * delta0
        ! write (*,'( a,g10.3 )') 'Grid spacing along X=',REAL(delta(1))
    endif
    
    m1=index(st,'<Y_Dim_Adj>')
    if (m1>=1) then
        m2=index(st,'</Y_Dim_Adj>')
        ch_e=trim(st(m1+11:m2-1))
        delta(2)=numeric(ch_e) * delta0
        ! write (*,'( a,g10.3 )') 'Grid spacing along Y=',REAL(delta(2))
    endif
    
    m1=index(st,'<Z_Dim_Adj>')
    if (m1>=1) then
        m2=index(st,'</Z_Dim_Adj>')
        ch_e=trim(st(m1+11:m2-1))
        delta(3)=numeric(ch_e) * delta0
        ! write (*,'( a,g10.3 )') 'Grid spacing along Z=',REAL(delta(3))
    endif

    if (index(st,'Structure Compression="ZLIB"') >= 1) then
        uncompress = .true.
    elseif (index(st,'Structure Compression="ASCII_READABLE"') >= 1) then
        uncompress = .false.
    endif
    
    m1=index(st,'<X_Voxels>')
    if (m1>=1) then
        m2=index(st,'</X_Voxels>')
        ch_e=st(m1+10:m2-1)
        read(ch_e,'(i3)')sdx
        ! write (*,'( a,i3 )') 'Number of cells along X=',sdx
    endif

    m1=index(st,'<Y_Voxels>')
    if (m1>=1) then
        m2=index(st,'</Y_Voxels>')
        ch_e=st(m1+10:m2-1)
        read(ch_e,'(i3)')sdy
        ! write (*,'( a,i3 )') 'Number of cells along Y=',sdy
    endif

    m1=index(st,'<Z_Voxels>')
    if (m1>=1 .and. m==0) then
        m2=index(st,'</Z_Voxels>')
        ch_e=st(m1+10:m2-1)
        read(ch_e,'(i3)')sdz
        ! write (*,'( a,i3/ )') 'Number of cells along Z=',sdz

        allocate(v(sdx*sdy*sdz), source=0_1) 

        read(1,'(a)',end=91,iostat=ios) ch_e ! line skip  <Data>
        if ( uncompress ) then
        ! START UNCOMPRESS
            print '(a )','start uncompress'
            ios=0
            open(2,file='compress.txt', iostat=ios) 
            
            idev=1
            do j=1,sdz
                call readline(st,ch, ier, idev)             !================== 2 dev=1
                k=index(st,'CDATA')
                k2=index(st,']]></')
                if (k2 == 0) then
                    stop ' Error: long string CDATA'
                endif
                k1=k+6
                ch=st(k1:k2-1)                        
                lst = len_trim(ch)
                write(2, '( a )') ch(1:lst)
                m = m + lst
            enddo
            close(2)

            call system('python uncompress_zlib.py', i)
            write(*, '(a,i3)') "Uncompress.txt created, status=", i 

            open(2,file='uncompress.txt', iostat=ios)
            st=''; m=0
            idev = 2
            do j=1,sdz
                call readline(st,ch, ier, idev)               !=======     3 idev= 2
                lst = len_trim(st)
                lst=lst-1
                do i=1,lst
                    v(m+i) = index(letter, st(i:i))
                enddo
                m = m + lst
            enddo
            close(2)
            
            ! for unix
            call get_environment_variable("HOME", st)
            if (len_trim(st) > 0)   &   ! for unix
            call system('rm compress.txt; rm uncompress.txt ',i)
            
            ! for windows
            call get_environment_variable("USERPROFILE", st)
            if (len_trim(st) > 0)   &   ! for windows
            call system('del compress.txt; del uncompress.txt ',i)
            
            print *, "delete files uncompress.txt and compress.txt, status: ", i 
        else
            idev=1
            do j=1,sdz
                call readline(st,ch,ier, idev)              !--------------4  idev=1
                k=index(st,'CDATA')
                k2=index(st,']]></')
                if (k2 == 0)   stop ' Error: long string CDATA'
                k1=k+6
                ch=st(k1:k2-1)
                lst = len_trim(ch)
                do i=1,lst
                    v(m+i) = index(letter, ch(i:i))
                enddo
                m = m + lst
            enddo
        endif
    endif

91 if (ios/=0) exit

end do convert_strings0
!=============================

nPHYS=maxval(v)                     ! number of physical domains only
Cells = sdx*sdy*sdz

j=0; k=1
do i=1,Cells
    if (v(i) == 0) then
        j = j + 1                    ! new cell
        if (j == 50000000) then      !3500 2200  300000  
            j=0                      ! reset cells counter
            k = k + 1                ! new domain
        endif
        v(i) = nPHYS + k             ! another one added for environment
    endif
enddo
if (j==0) k = k - 1                  ! there was no enlargement of the cells
nPHYS_env = k

nPHYS_glob = nPHYS+nPHYS_env

ALLOCATE ( typPHYS(nPHYS_glob), namePHYS(nPHYS_glob ) )

typPHYS=''; namePHYS='';

allocate ( valPHYS(nPHYS_glob), source=0.d0 )

! default param for ENVIRON
if (nPHYS_env /=0) then
    do j=nPHYS+1, nPHYS_glob
        typPHYS(j)   = 'R     '
        namePHYS(j)  = 'ENV   '
        valPHYS(j) = 1.d0 
        ! write (*, '(a)') 'created a default environment with parameters: '
        ! write (*, '(a,i3,a/)') 'nPHYS= ',j,' valPHYS= 1.0 typPHYS= "R", namePHYS= "ENV"'
    enddo
endif


rewind(1)

! ------------filling arrays-----------!
kp=0; ios=0;

idev = 1
convert_strings1: &
do while ( 1==1)
    call readline(st,ch, ier, idev)              !-------------- 5
    if (st(1:6)=='</VXC>' ) exit
    k0=index(st,'<Name>')
    if ( k0 > 0 ) then 
        kp = kp + 1      ! count all lines with <Name>
        m2=index(st,'</')
        st=st(k0+6:m2-1)
        i=1
        do
            k = index( st(i:),'=')
            if (k/=0) then
                i=i+k-1; st(i:i) = ' '
            else
                exit
            endif
        enddo

        call Upp (st, st)
        call string2words(trim(st),words,neww)

        do i=2,neww
            if ( words(i)(1:1)=='Z' .and. kp <= nPHYS ) then  ! Key word Z 
                valPHYS(kp) = numeric( words(i+1) )
                typPHYS(kp) = 'R'
                namePHYS(kp) = trim(words(1))
                                  
                IF (i+1 == neww) then
                    ! write (*, '( "nPHYS= ",i3,  " valPHYS= ",e10.3,  " typPHYS= ", a,   " namPHYS= ", a  )') &
                                  ! kp ,  valPHYS(kp),       trim(typPHYS(kp)),   trim(namePHYS(kp))
                    CYCLE ! it is last word
                endif

            elseif ( index( words(i),'ENVIRON') /=0  )  then  !   ! Key word ENVIRON
                do j=i+1, neww-1
                    if     (  words(j)(1:1)=='Z' ) then 
                        valPHYS(nPHYS_glob) = numeric(words(j + 1))
                    endif
                    
                    if     (  index (words(j),'FILE') /=0 ) then 
                        files = trim ( words(j + 1))
                    endif
                    
                enddo
                namePHYS(nPHYS_glob) = trim(words(1))
                
                ! write (*, '( "new environment with parameters:",/ "nPHYS= ",i3,  " valPHYS= ",e10.3,  " typPHYS= ", a,   " namPHYS= ", a  )') &
                        ! nPHYS_glob ,  valPHYS(nPHYS_glob),       trim(typPHYS(nPHYS_glob)),   trim(namePHYS(nPHYS_glob))
            endif  

        enddo  ! words

    endif ! <Name>

92  if (ios/=0) exit

end do convert_strings1
    
close(1)

allocate (geoPHYS(sdx,sdy,sdz), source=0_1  )
geoPHYS = reshape(v, shape =(/sdx,sdy,sdz/) ,  order = (/ 1, 2, 3/) )

allocate (Ncells(nPHYS), source=0)

do m = 1, nPHYS
    do k=1,sdz;    do j=1, sdy;    do i=1, sdx
        if (geoPHYS(i,j,k) == m) Ncells(m) = Ncells(m) + 1
    enddo; enddo; enddo
enddo

if (files == '') files='out'

! --------- Data output -------------- !

open(10,file=files//'.txt')

write (10,'( a,g10.3 )') 'Lattice Dimension=',delta0
write (10,'( a,g10.3 )') 'Grid spacing along X=',REAL(delta(1))
write (10,'( a,g10.3 )') 'Grid spacing along Y=',REAL(delta(2))
write (10,'( a,g10.3 )') 'Grid spacing along Z=',REAL(delta(3))
write (10,'( a,i3 )')  'Number of cells along X=',sdx
write (10,'( a,i3 )')  'Number of cells along Y=',sdy
write (10,'( a,i3/ )') 'Number of cells along Z=',sdz

write (10,'( "Cells= ",i9, / "nPHYS= ", i3, " nPHYS_env= ", i3 ," nPHYS_glob=nPHYS+nPHYS_env=", i3/ )') &
            Cells ,          nPHYS,          nPHYS_env,          nPHYS_glob
            
do i=1,nPHYS
write (10, '( "numPHYS= ",i3,  " valPHYS= ",e10.3,  " typPHYS= ", a,   " namPHYS= ", a  )') &
                                  i ,  valPHYS(i),       trim(typPHYS(i)),   trim(namePHYS(i))
enddo

write (10, '( "created environment with parameters:",/ "numPHYS= ",i3,  " valPHYS= ",e10.3,  " typPHYS= ", a,   " namPHYS= ", a  )') &
               nPHYS_glob ,  valPHYS(nPHYS_glob),       trim(typPHYS(nPHYS_glob)),   trim(namePHYS(nPHYS_glob))

if (Cells < 1000) then
    write (10,'( a)')''
    write (10,'( a)') '3d matrix of geometric location of physical domains:'
    do k=1,sdz
        write (10,'( a)')''
        write (10,'( a,i3)')'          k=' ,k
        write (10,'( a )') ' j ^' 
        do j=sdy,1,-1
            write (10,'(i2, a,$ )')  j, ' |'
            write (10,'(  *(i5) )')(geoPHYS(i,j,k),i=1,sdx)
        enddo
        write (10,'( a, $)')'      ' 
        write (10,'(  *( a ) )' )  ( '-----',i=1,sdx),'>'
        write (10,'( a, $)')'    i= ' 
        write (10,'(  *( i2,a ) )' )  ( i,'   ',i=1,sdx)
    enddo
    call writeVtk     ( nPHYS, Ncells, sdx, sdy, sdz, delta, geoPHYS, files)
    write(*,'( a)') 'converted data saved in files: '//files//'.txt and '//files//'.vtk'
    
else
    call writeVtk_bin ( nPHYS, Ncells, sdx, sdy, sdz, delta, geoPHYS, files)
    write(*,'( a)') 'converted data saved in files: '//files//'.txt and '//files//'_bin.vtk'
endif


deallocate (v)
close(10)


!========================================

contains


subroutine readline(line,ch,ier, idev)

    character(len=:),allocatable,intent(out) :: line ,ch
    integer,intent(out)                      :: ier
    integer                                  :: idev
    character :: buffer

    line=''
    ier=0
    do
        read( idev, '(a)', advance='no', iostat=ier ) buffer
        if ( ier /= 0 ) then
            exit
        else 
            line = line // buffer
        endif
    end do
    
    ch = line
end subroutine readline

end program vxc2data





