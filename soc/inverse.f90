  subroutine inv(ndim,Amat)
     ! A pack for Lapack subroutine zgesv to get the inverse of 
     ! a generic complex matrix


     implicit none

     integer,parameter :: dp=8

     integer           :: i
     integer           :: info

     integer,intent(in):: ndim

!    IPIV   : INTEGER. Array, DIMENSION at least max(1, n). The pivot indices that define
!    the permutation matrix P; row i of the matrix was interchanged with
!    row ipiv(i). Corresponds to the single precision factorization (if info=
!    0 and iter ≥ 0) or the double precision factorization (if info= 0 and
!    iter < 0).
     integer,allocatable   :: ipiv(:)


     complex(dp),parameter :: zone=(1.0d0,0.0d0)



!    Amat  :
!    Overwritten by the factors L and U from the factorization of A = P*L*U;
!    the unit diagonal elements of L are not stored.
!    If iterative refinement has been successfully used (info= 0 and iter≥
!    0), then A is unchanged.
!    If double precision factorization has been used (info= 0 and iter <
!    0), then the array A contains the factors L and U from the factorization
!    A = P*L*U; the unit diagonal elements of L are not stored.
     complex(dp),intent(inout):: Amat(ndim,ndim)

!    Bmat  :
!    Overwritten by the solution matrix X for dgesv, sgesv,zgesv,zgesv.
     complex(dp),allocatable :: Bmat(:,:)


     allocate(ipiv(ndim))
     allocate(Bmat(ndim,ndim))

     ipiv=0

     ! unit matrix
     Bmat= (0d0, 0d0)
     do i=1,ndim
        Bmat(i,i)= zone
     enddo

     call zgesv(ndim,ndim,Amat,ndim,ipiv,Bmat,ndim,info)

     if(info.ne.0)print *,'something wrong with zgesv'

     Amat=Bmat
     
     return
  end subroutine inv 

  subroutine inv_r(ndim,Amat)
     ! A pack for Lapack subroutine zgesv to get the inverse of 
     ! a generic real matrix


     implicit none

     integer,parameter :: dp=8

     integer           :: i
     integer           :: info

     integer,intent(in):: ndim

!    IPIV   : INTEGER. Array, DIMENSION at least max(1, n). The pivot indices that define
!    the permutation matrix P; row i of the matrix was interchanged with
!    row ipiv(i). Corresponds to the single precision factorization (if info=
!    0 and iter ≥ 0) or the double precision factorization (if info= 0 and
!    iter < 0).
     integer,allocatable   :: ipiv(:)


!    Amat  :
!    Overwritten by the factors L and U from the factorization of A = P*L*U;
!    the unit diagonal elements of L are not stored.
!    If iterative refinement has been successfully used (info= 0 and iter≥
!    0), then A is unchanged.
!    If double precision factorization has been used (info= 0 and iter <
!    0), then the array A contains the factors L and U from the factorization
!    A = P*L*U; the unit diagonal elements of L are not stored.
     real(dp),intent(inout):: Amat(ndim,ndim)

!    Bmat  :
!    Overwritten by the solution matrix X for dgesv, sgesv,zgesv,zgesv.
     real(dp),allocatable :: Bmat(:,:)


     allocate(ipiv(ndim))
     allocate(Bmat(ndim,ndim))

     ipiv=0

     ! unit matrix
     Bmat= 0d0 
     do i=1,ndim
        Bmat(i,i)= 1d0
     enddo

     call dgesv(ndim,ndim,Amat,ndim,ipiv,Bmat,ndim,info)

     if(info.ne.0)print *,'something wrong with dgesv'

     Amat=Bmat
     
     return
  end subroutine inv_r
