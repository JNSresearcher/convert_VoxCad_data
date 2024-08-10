subroutine string2words(Ligne,words,num_words)
! converting string to array of words

character (len = *), intent(IN) :: Ligne
character (len = *), intent(OUT):: words(*) 

integer Start
integer Finish
integer num_words
INTEGER length 

    length= len_trim(Ligne)
    Start = 0;Finish=0;num_words=0
    do j = 1, length
        if (Ligne(j:j) == ' ' .or. Ligne(j:j)== '	') then
            if (Start > 0) then
                num_words = num_words + 1
                words(num_words) = Ligne(Start : Finish)
                Start = 0
            endif
        else if (Start == 0) then
            Start = j;  Finish = j
        else
            Finish = Finish + 1
        endif
    end do
    if (Start > 0) then
        num_words = num_words + 1
        words(num_words) = Ligne(Start : Finish)
        Start = 0
    endif

end subroutine string2words

SUBROUTINE UPP (st1, st2)
! convert string to uppercase

    CHARACTER (LEN=*),  INTENT(in) :: st1
    CHARACTER (LEN=*),  INTENT(out) :: st2
    CHARACTER (LEN=*),  PARAMETER :: lsym = 'abcdefghijklmnopqrstuvwxyz'
    CHARACTER (LEN=*),  PARAMETER :: usym = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    st2 = st1
    DO j=1,LEN_TRIM(st1)
        k = INDEX(lsym,st1(j:j))
        IF (k > 0) st2(j:j) = usym(k:k)
    END DO
END SUBROUTINE UPP



REAL(8) function numeric(SA)
! Converting a string to a number using decimal prefixes.
! (Spice-like format used in electrical circuit simulation programs).

! Accepted decimal prefix formats in this program  :
! 'f',   femto    1e-15
! 'p',   pico     1e-12
! 'n',   nano     1e-9
! 'u',   micro    1e-6
! 'm',   milli    1e-3
! 'c'    centi    1e-2
! 'h'    hecto    1e+2
! 'k',   kilo     1e+3
! 'meg'  mega     1e+6
! 'g',   giga     1e+9
! 't',   tera     1e+12
! 'pet'  peta     1e+15

! some features:
! - Instead of a decimal point, on can use prefix: 1.3k = 1k3
! - can put exponential format at the end 
! - spaces inside are ignored
! - prefixes are automatically converted to uppercase during conversion

! Fortran code created by J.Sochor  ( https://github.com/JNSresearcher/convert_prefix ) 

implicit none

    CHARACTER (LEN=*),   PARAMETER :: lsym = 'abcdefghijklmnopqrstuvwxyz'
    CHARACTER (LEN=*),   PARAMETER :: usym = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

!                                            1   2   3    4    5   6   7    8    9   10
    character*1 :: ch, prefixes(10)=(/'M','K','U','N','P','G','T','F','C','H'/) 
    character (len = *) sa
    character*20 inter
    REAL(8) multiplier   ! factor
    integer i,j,k,l,m,n,lm,lp
    
    DO j=1,LEN_TRIM(SA)
        k = INDEX(lsym,SA(j:j))
        IF (k > 0) SA(j:j) = usym(k:k)
    END DO
    
    multiplier = 1.0_8;lp=0;lm=0
    m=index(SA,',')
    if (m /= 0) SA(m:m)='.'

    do i=1,9
        l=index(SA,prefixes(i))
        if (l /= 0) then
            m=1
            if ( i == 1 ) then
                lm=index(SA, 'MEG')
                if (lm /= 0) then
                    multiplier=1.E+06_8
                else
                    multiplier=1.E-3_8
                endif
                exit
            elseif ( i == 2 ) then
                multiplier=1.E+03_8
                exit
            elseif ( i == 3 ) then
                multiplier=1.E-06_8
                exit
            elseif ( i == 4 ) then
                multiplier=1.E-09_8
                exit
            elseif ( i == 5 ) then
                lm=index(SA, 'PET')
                if (lm /= 0) then
                    multiplier=1.E+15_8
                else
                    multiplier=1.E-12_8
                endif
                exit
            elseif ( i == 6 ) then
                multiplier=1.E+09_8
                exit
            elseif ( i == 7 ) then
                multiplier=1.E+12_8
                exit
            elseif ( i == 8 ) then
                multiplier=1.E-15_8
                exit
            elseif ( i == 9 ) then
                multiplier=1.E-2_8
                exit
            elseif ( i == 10 ) then
                multiplier=1.E+2_8
                exit
            endif
        endif
    enddo

    k=index(SA,'.');j=len_trim(SA)
    if ( (m /= 0).AND.(k==0) )  then
        SA(l:l)='.'
        if (lm /= 0) then
            SA(lm+1:lm+2)='  '
            sa=sa(1:lm)//sa(lm+3:)
        endif
    else
        if (lm /= 0) then
            SA(lm:lm+2)='   '
            sa=sa(1:lm)//sa(lm+3:)
        endif
    endif
    
    j=len_trim(SA)
    i=index(SA,'E')
    if (i == 0) then
        do i=1,j
            ch=SA(i:i)
            n=ichar(ch)
            if ((n.GT.57).or.(n.LT.48)) then
                if     (ch == '.') then
                    cycle
                elseif (ch == '-') then
                    cycle
                else
                    SA(i:i) = ' '
                endif
            else
            endif
        enddo
    endif

    write (inter,'(A20)') SA(1:j)
    !print*,SA(1:j)
    read (inter,'(BN,G20.0)') numeric
    numeric = multiplier * numeric

end function numeric

