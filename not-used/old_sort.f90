! Currently this code is taken from NR. It should be replaced
! by GPL code as soon as possible?
module sort
  use types
  use warnings_and_errors
  use assertions
  implicit none
  private

  INTEGER, PARAMETER :: I4B = SELECTED_INT_KIND(9)

  interface icomp_xchg
     module procedure icomp_xchg_dp, icomp_xchg_i4b
  end interface
  
  
  interface indexx
     module procedure indexx_dp, indexx_i4b
  end interface
  
  public :: indexx

contains

  ! arranges for index to corresond to order of ascending arr?
  SUBROUTINE indexx_dp(arr,index)
    IMPLICIT NONE
    REAL(DP), DIMENSION(:), INTENT(IN) :: arr
    INTEGER(I4B), DIMENSION(:), INTENT(OUT) :: index
    INTEGER(I4B), PARAMETER :: NN=15, NSTACK=50
    REAL(DP) :: a
    INTEGER(I4B) :: n,k,i,j,indext,jstack,l,r
    INTEGER(I4B), DIMENSION(NSTACK) :: istack
    n=assert_eq(size(index),size(arr),'indexx_sp')
    index=arth(1,1,n)
    jstack=0
    l=1
    r=n
    do
       if (r-l < NN) then
          do j=l+1,r
             indext=index(j)
             a=arr(indext)
             do i=j-1,l,-1
                if (arr(index(i)) <= a) exit
                index(i+1)=index(i)
             end do
             index(i+1)=indext
          end do
          if (jstack == 0) RETURN
          r=istack(jstack)
          l=istack(jstack-1)
          jstack=jstack-2
       else
          k=(l+r)/2
          call swap(index(k),index(l+1))
          call icomp_xchg(arr,index(l),index(r))
          call icomp_xchg(arr,index(l+1),index(r))
          call icomp_xchg(arr,index(l),index(l+1))
          i=l+1
          j=r
          indext=index(l+1)
          a=arr(indext)
          do
             do
                i=i+1
                if (arr(index(i)) >= a) exit
             end do
             do
                j=j-1
                if (arr(index(j)) <= a) exit
             end do
             if (j < i) exit
             call swap(index(i),index(j))
          end do
          index(l+1)=index(j)
          index(j)=indext
          jstack=jstack+2
          if (jstack > NSTACK) call wae_error('indexx: NSTACK too small')
          if (r-i+1 >= j-l) then
             istack(jstack)=r
             istack(jstack-1)=i
             r=j-1
          else
             istack(jstack)=j-1
             istack(jstack-1)=l
             l=i
          end if
       end if
    end do
  END SUBROUTINE indexx_dp

  SUBROUTINE indexx_i4b(iarr,index)
    IMPLICIT NONE
    INTEGER(I4B), DIMENSION(:), INTENT(IN) :: iarr
    INTEGER(I4B), DIMENSION(:), INTENT(OUT) :: index
    INTEGER(I4B), PARAMETER :: NN=15, NSTACK=50
    INTEGER(I4B) :: a
    INTEGER(I4B) :: n,k,i,j,indext,jstack,l,r
    INTEGER(I4B), DIMENSION(NSTACK) :: istack
    n=assert_eq(size(index),size(iarr),'indexx_sp')
    index=arth(1,1,n)
    jstack=0
    l=1
    r=n
    do
       if (r-l < NN) then
          do j=l+1,r
             indext=index(j)
             a=iarr(indext)
             do i=j-1,1,-1
                if (iarr(index(i)) <= a) exit
                index(i+1)=index(i)
             end do
             index(i+1)=indext
          end do
          if (jstack == 0) RETURN
          r=istack(jstack)
          l=istack(jstack-1)
          jstack=jstack-2
       else
          k=(l+r)/2
          call swap(index(k),index(l+1))
          call icomp_xchg(iarr,index(l),index(r))
          call icomp_xchg(iarr,index(l+1),index(r))
          call icomp_xchg(iarr,index(l),index(l+1))
          i=l+1
          j=r
          indext=index(l+1)
          a=iarr(indext)
          do
             do
                i=i+1
                if (iarr(index(i)) >= a) exit
             end do
             do
                j=j-1
                if (iarr(index(j)) <= a) exit
             end do
             if (j < i) exit
             call swap(index(i),index(j))
          end do
          index(l+1)=index(j)
          index(j)=indext
          jstack=jstack+2
          if (jstack > NSTACK) call wae_error('indexx: NSTACK too small')
          if (r-i+1 >= j-l) then
             istack(jstack)=r
             istack(jstack-1)=i
             r=j-1
          else
             istack(jstack)=j-1
             istack(jstack-1)=l
             l=i
          end if
       end if
    end do
  END SUBROUTINE indexx_i4b

  !BL
  SUBROUTINE icomp_xchg_dp(arr,i,j)
    real(dp),     intent(in)    :: arr(:)
    INTEGER(I4B), INTENT(INOUT) :: i,j
    INTEGER(I4B) :: swp
    if (arr(j) < arr(i)) then
       swp=i
       i=j
       j=swp
    end if
  END SUBROUTINE icomp_xchg_dp
  !BL
  SUBROUTINE icomp_xchg_i4b(iarr,i,j)
    integer(i4B), intent(in)    :: iarr(:)
    INTEGER(I4B), INTENT(INOUT) :: i,j
    INTEGER(I4B) :: swp
    if (iarr(j) < iarr(i)) then
       swp=i
       i=j
       j=swp
    end if
  END SUBROUTINE icomp_xchg_i4b

  FUNCTION arth(first,increment,n)
    INTEGER(I4B), INTENT(IN) :: first,increment,n
    INTEGER(I4B), DIMENSION(n) :: arth
    INTEGER(I4B) :: k,k2,temp
    INTEGER(I4B), PARAMETER :: NPAR_ARTH=16,NPAR2_ARTH=8
    if (n > 0) arth(1)=first
    if (n <= NPAR_ARTH) then
       do k=2,n
          arth(k)=arth(k-1)+increment
       end do
    else
       do k=2,NPAR2_ARTH
          arth(k)=arth(k-1)+increment
       end do
       temp=increment*NPAR2_ARTH
       k=NPAR2_ARTH
       do
          if (k >= n) exit
          k2=k+k
          arth(k+1:min(k2,n))=temp+arth(1:min(k,n-k))
          temp=temp+temp
          k=k2
       end do
    end if
  END FUNCTION arth

  SUBROUTINE swap(a,b)
    INTEGER(I4B), INTENT(INOUT) :: a,b
    INTEGER(I4B) :: dum
    dum=a
    a=b
    b=dum
  END SUBROUTINE swap

end module sort
