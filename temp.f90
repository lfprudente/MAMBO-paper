module myproblem

use modprobdata, only: problem,seed

implicit none

  ! SUBROUTINES
  public :: inip, evalf, evalg, evalh

contains

    subroutine inip(n,m,x,l,u,strconvex,scaleF,checkder)

    implicit none

    ! SCALAR ARGUMENTS
    integer, intent(out) :: n,m
    logical, intent(out) :: scaleF,checkder
    
    ! ARRAY ARGUMENTS
    real(kind=8), allocatable, intent(out) :: x(:),l(:),u(:)
    logical, allocatable, intent(out) :: strconvex(:)
    
    ! LOCAL SCALARS
    integer :: i,allocerr
    real, parameter :: pi = 3.1415927
    real(kind=8), parameter :: zero = 0.0d0, one = 1.0d0
    real(kind=8) :: a,b
    
    ! FUNCTIONS
    real(kind=8) :: drand
        
    ! ----------------------------------------------------------------------

    ! AP1
    ! Example 1 of "A modified Quasi-Newton method for vector optimization problem"

    if ( problem == 'AP1' ) then 
    
        ! Number of variables
        
        n = 2
        
        ! Number of objectives
        
        m = 3
        
        allocate(x(n),l(n),u(n),strconvex(m),stat=allocerr)
        if ( allocerr .ne. 0 ) then
             write(*,*) 'Allocation error in main program'
             stop
        end if
        
        ! Strictly strconvex?
        
        strconvex(1) = .false.
        strconvex(2) = .true.
        strconvex(3) = .true.
        
        ! Box constraints
    
        l(:) = - 1.0d1
        u(:) =   1.0d1
        
        ! Initial point
        
        do i = 1,n
            x(i) = l(i) + ( u(i) - l(i) ) * drand(seed)
        end do 
        
        ! Scale the problem?
            
        scaleF   = .true.

        ! Check derivatives?

        checkder = .false.
        
        return
    end if
    
! ----------------------------------------------------------------------

    ! AP2
    ! Example 2 of "A modified Quasi-Newton method for vector optimization problem"

    if ( problem == 'AP2' ) then 
    
        ! Number of variables
        
        n = 1
        
        ! Number of objectives
        
        m = 2
        
        allocate(x(n),l(n),u(n),strconvex(m),stat=allocerr)
        if ( allocerr .ne. 0 ) then
             write(*,*) 'Allocation error in main program'
             stop
        end if
        
        ! Strictly strconvex?
        
        strconvex(1) = .true.
        strconvex(2) = .true.
        
        ! Box constraints
    
        l(:) = - 1.0d2
        u(:) =   1.0d2
        
        ! Initial point
        
        a = - 1.0d2
        b =   1.0d2
        
        do i = 1,n
            x(i) = a + ( b - a ) * drand(seed)
        end do 
        
        ! Scale the problem?
            
        scaleF   = .true.

        ! Check derivatives?

        checkder = .false.
    end if
    
! ----------------------------------------------------------------------

    ! AP3
    ! Example 3 of "A modified Quasi-Newton method for vector optimization problem"

    if ( problem == 'AP3' ) then 
    
        ! Number of variables
        
        n = 2
        
        ! Number of objectives
        
        m = 2
        
        allocate(x(n),l(n),u(n),strconvex(m),stat=allocerr)
        if ( allocerr .ne. 0 ) then
             write(*,*) 'Allocation error in main program'
             stop
        end if
        
        ! Strictly strconvex?
        
        strconvex(1) = .false.
        strconvex(2) = .false.
        
        ! Box constraints
        
        l(:) = - 1.0d2
        u(:) =   1.0d2
        
        ! Initial point
        
        a = - 1.0d2
        b =   1.0d2
        
        do i = 1,n
            x(i) = a + ( b - a ) * drand(seed)
        end do 
        
        ! Scale the problem?
            
        scaleF   = .true.

        ! Check derivatives?

        checkder = .false.
        
        return
    end if
    
! ----------------------------------------------------------------------

    ! AP4
    ! Exemple 4 of "A modified Quasi-Newton method for vector optimization problem"
    
    if ( problem == 'AP4' ) then 
    
        ! Number of variables
        
        n = 3
        
        ! Number of objectives
        
        m = 3
        
        allocate(x(n),l(n),u(n),strconvex(m),stat=allocerr)
        if ( allocerr .ne. 0 ) then
             write(*,*) 'Allocation error in main program'
             stop
        end if
        
        ! Strictly strconvex?
        
        strconvex(1) = .false.
        strconvex(2) = .true.
        strconvex(3) = .true.
        
        ! Box constraints
    
        l(:) = - 1.0d1
        u(:) =   1.0d1
        
        ! Initial point
        
        a = - 1.0d1
        b =   1.0d1
        
        do i = 1,n
            x(i) = a + ( b - a ) * drand(seed)
        end do 
        
        ! Scale the problem?
            
        scaleF   = .true.

        ! Check derivatives?

        checkder = .false.
        
        return
    end if
    
! ----------------------------------------------------------------------

    !  BK1
    !  A Review of Multiobjective Test Problems and a Scalable Test Problem Toolkit
    
    if ( problem == 'BK1' ) then 
    
        ! Number of variables
        
        n = 2
        
        ! Number of objectives
        
        m = 2
        
        allocate(x(n),l(n),u(n),strconvex(m),stat=allocerr)
        if ( allocerr .ne. 0 ) then
             write(*,*) 'Allocation error in main program'
             stop
        end if
        
        ! Strictly strconvex?
        
        strconvex(1) = .true.
        strconvex(2) = .true.
        
        ! Box constraints
    
        l(:) = - 5.0d0
        u(:) =   1.0d1
        
        ! Initial point
        
        a = - 5.0d0
        b =   1.0d1
        
        do i = 1,n
            x(i) = a + ( b - a ) * drand(seed)
        end do 
        
        ! Scale the problem?
            
        scaleF   = .true.

        ! Check derivatives?

        checkder = .false.
        
        return
    end if	
        
    end subroutine inip

    !***********************************************************************
    !***********************************************************************

    subroutine evalf(n,x,f,ind)
    
    implicit none

    ! SCALAR ARGUMENTS
    integer, intent(in) :: n,ind
    real(kind=8), intent(out) :: f
    
    ! ARRAY ARGUMENTS
    real(kind=8), intent(in) :: x(n)
    
    ! LOCAL SCALARS
    integer :: i,k
    real, parameter :: pi = 3.141592653589793
    real(kind=8) :: faux,A1,A2,A3,B1,B2,a,b,t,y
    
    ! LOCAL ARRAYS
    real(kind=8) :: M(3,3),p(3)
    
! ----------------------------------------------------------------------

  ! AP1: Exemple 1 of "A modified Quasi-Newton method for vector optimization problem"

    if ( problem == 'AP1' ) then 
            
        if ( ind == 1 ) then
            f = 0.25d0 * ( ( x(1) - 1.0d0 ) ** 4 + 2.0d0 * ( x(2) - 2.0d0 ) ** 4 )
            return
        end if
        
        if ( ind == 2 ) then
            f = exp( ( x(1) + x(2) ) / 2.0d0 ) + x(1) ** 2 + x(2) ** 2
            return
        end if	
        
        if ( ind == 3 ) then
            f = 1.0d0/6.0d0 * ( exp( - x(1) ) + 2.0d0 * exp( - x(2) ) )
            return
        end if	
            
    end if		
    
! ----------------------------------------------------------------------

  ! AP2: Exemple 2 of "A modified Quasi-Newton method for vector optimization problem"

    if ( problem == 'AP2' ) then 
            
        if ( ind == 1 ) then
            f = x(1) ** 2 - 4.0d0
            return
        end if
        
        if ( ind == 2 ) then
            f = ( x(1) - 1.0d0 ) ** 2
            return
        end if	
            
    end if		
    
! ----------------------------------------------------------------------

  ! AP3: Exemple 3 of "A modified Quasi-Newton method for vector optimization problem"

    if ( problem == 'AP3' ) then 
            
        if ( ind == 1 ) then
            f = 0.25d0 * ( ( x(1) - 1.0d0 ) ** 4 + 2.0d0 * ( x(2) - 2.0d0 ) ** 4 )
            return
        end if
        
        if ( ind == 2 ) then
            f = ( x(2) - x(1) ** 2 ) ** 2 + ( 1.0d0 - x(1) ) ** 2
            return
        end if	
            
    end if	
    
! ----------------------------------------------------------------------

 ! AP4: Exemple 4 of "A modified Quasi-Newton method for vector optimization problem"

    if ( problem == 'AP4' ) then 
            
        if ( ind == 1 ) then
            f = 1.0d0/9.0d0 * ( ( x(1) - 1.0d0 ) ** 4 + 2.0d0 * ( x(2) - 2.0d0 ) ** 4 &
            + 3.0d0 * ( x(3) - 3.0d0 ) ** 4 )
            return
        end if
        
        if ( ind == 2 ) then
            f = exp( ( x(1) + x(2) + x(3) ) / 3.0d0 ) + x(1) ** 2 + x(2) ** 2 + x(3) ** 2
            return
        end if	
        
        if ( ind == 3 ) then
            f = 1.0d0/1.2d1 * ( 3.0d0 * exp( -x(1) ) + 4.0d0 * exp( - x(2) ) + 3.0d0 * exp( -x(3) ) )
            return
        end if	
            
    end if		
    
! ----------------------------------------------------------------------	

    !  BK1
    !  A Review of Multiobjective Test Problems and a Scalable Test Problem Toolkit
    
    if ( problem == 'BK1' ) then 
        
        if ( ind == 1 ) then
            f = x(1) ** 2 + x(2) ** 2
            return
        end if
        
        if ( ind == 2 ) then
            f = ( x(1) - 5.0d0 ) ** 2 + ( x(2) - 5.0d0 ) ** 2
            return	
        end if
                    
    end if					
            
    end subroutine evalf	

    !***********************************************************************
    !***********************************************************************

    subroutine evalg(n,x,g,ind)
    
    implicit none
    
    ! SCALAR ARGUMENTS
    integer, intent(in) :: n,ind
    
    ! ARRAY ARGUMENTS
    real(kind=8), intent(in)  :: x(n)
    real(kind=8), intent(out) :: g(n)
    
    ! LOCAL SCALARS
    
    integer :: i,j,k
    real, parameter :: pi = 3.141592653589793
    real(kind=8) :: faux,gaux1,gaux2,A1,A2,A3,B1,B2,a,b,t
    
    ! LOCAL ARRAYS
    real(kind=8) :: M (3,3),p(3)
    
! ----------------------------------------------------------------------

! AP1: Exemple 1 of "A modified Quasi-Newton method for vector optimization problem"

    if ( problem == 'AP1' ) then 
            
        if ( ind == 1 ) then
            g(1) = ( x(1) - 1.0d0 ) ** 3 
            g(2) = 2.0d0 * ( x(2) - 2.0d0 ) ** 3
            return
        end if
        
        if ( ind == 2 ) then
            g(1) = 0.5d0 * exp( ( x(1) + x(2) ) / 2.0d0 ) + 2.0d0 * x(1)
            g(2) = 0.5d0 * exp( ( x(1) + x(2) ) / 2.0d0 ) + 2.0d0 * x(2)
            return
        end if	
        
        if ( ind == 3 ) then
            g(1) = - 1.0d0/6.0d0 * exp( - x(1) ) 
            g(2) = - 1.0d0/3.0d0 * exp( - x(2) ) 
            return
        end if	
        
    end if	
    
! ----------------------------------------------------------------------

! AP2: Exemple 2 of "A modified Quasi-Newton method for vector optimization problem"

    if ( problem == 'AP2' ) then 
            
        if ( ind == 1 ) then
            g(1) = 2.0d0 * x(1)
            return
        end if
        
        if ( ind == 2 ) then
            g(1) = 2.0d0 * ( x(1) - 1.0d0 )
            return
        end if	
            
    end if		
    
    
! ----------------------------------------------------------------------

! AP3: Exemple 3 of "A modified Quasi-Newton method for vector optimization problem"

    if ( problem == 'AP3' ) then 
            
        if ( ind == 1 ) then
            g(1) = ( x(1) - 1.0d0 ) ** 3 
            g(2) = 2.0d0 * ( x(2) - 2.0d0 ) ** 3
            return
        end if
        
        if ( ind == 2 ) then
            g(1) = - 4.0d0 * x(1) * ( x(2) - x(1) ** 2 ) - 2.0d0 * ( 1.0d0 - x(1) )
            g(2) = 2.0d0 * ( x(2) - x(1) ** 2 )
            return
        end if	
            
    end if	
    
! ----------------------------------------------------------------------

! AP4: Exemple 4 of "A modified Quasi-Newton method for vector optimization problem"

    if ( problem == 'AP4' ) then 
        
        if ( ind == 1 ) then
            g(1) = 4.0d0/9.0d0 * ( x(1) - 1.0d0 ) ** 3 
            g(2) = 8.0d0/9.0d0 * ( x(2) - 2.0d0 ) ** 3 
            g(3) = 1.2d1/9.0d0 * ( x(3) - 3.0d0 ) ** 3 
            return
        end if
        
        if ( ind == 2 ) then
            g(1) = 1.0d0/3.0d0 * exp( ( x(1) + x(2) + x(3) ) / 3.0d0 ) + 2.0d0 * x(1) 
            g(2) = 1.0d0/3.0d0 * exp( ( x(1) + x(2) + x(3) ) / 3.0d0 ) + 2.0d0 * x(2) 
            g(3) = 1.0d0/3.0d0 * exp( ( x(1) + x(2) + x(3) ) / 3.0d0 ) + 2.0d0 * x(3) 
            return
        end if	
        
        if ( ind == 3 ) then
            g(1) = - 1.0d0/4.0d0 * exp( -x(1) ) 
            g(2) = - 1.0d0/3.0d0 * exp( -x(2) ) 
            g(3) = - 1.0d0/4.0d0 * exp( -x(3) ) 
            return
        end if	
            
    end if	
    
! ----------------------------------------------------------------------	

    !  BK1

    if ( problem == 'BK1' ) then 
            
        if ( ind == 1 ) then
            do i = 1,n
                g(i) = 2.0d0 * x(i)
            end do
            return
        end if
        
        if ( ind == 2 ) then
            do i = 1,n
                g(i) = 2.0d0 * ( x(i) - 5.0d0 ) 
            end do
            return	
        end if
            
    end if		

    end subroutine evalg
    
! ******************************************************************
! ******************************************************************

subroutine evalh(n,x,H,ind)

    implicit none

    ! SCALAR ARGUMENTS
    integer, intent(in) :: n,ind
    
    ! ARRAY ARGUMENTS
    real(kind=8), intent(in)  :: x(n)
    real(kind=8), intent(out) :: H(n,n)
    
    ! LOCAL SCALAR
    real, parameter :: pi = 3.141592653589793
    integer :: i,j
    real(kind=8) :: t,a,b,c,A1,A2,A3,B1,B2,gaux1,gaux2
    
    ! LOCAL ARRAYS
    real(kind=8) :: M(3,3),p(3)
    
    
! ----------------------------------------------------------------------

! AP1: Exemple 1 of "A modified Quasi-Newton method for vector optimization problem"

    if ( problem == 'AP1' ) then 
            
        if ( ind == 1 ) then
            H(:,:) = 0.0d0
            
            H(1,1) = 3.0d0 * ( x(1) - 1.0d0 ) ** 2
            H(2,2) = 6.0d0 * ( x(2) - 2.0d0 ) ** 2
            return
        end if
        
        if ( ind == 2 ) then			
            H(:,:) = 0.0d0
            
            H(1,1) = 0.5d0 * exp( ( x(1) + x(2) ) / 2.0d0 ) / 2.0d0 + 2.0d0
            H(1,2) = 0.5d0 * exp( ( x(1) + x(2) ) / 2.0d0 ) / 2.0d0 
            H(2,1) = H(1,2)
            H(2,2) = 0.5d0 * exp( ( x(1) + x(2) ) / 2.0d0 ) / 2.0d0 + 2.0d0
            return
        end if	
        
        if ( ind == 3 ) then
            H(:,:) = 0.0d0
            
            H(1,1) = 1.0d0/6.0d0 * exp( - x(1) )
            H(2,2) = 1.0d0/3.0d0 * exp( - x(2) ) 
            return
        end if	
        
    end if	
    
! ----------------------------------------------------------------------

! AP2: Exemple 2 of "A modified Quasi-Newton method for vector optimization problem"

    if ( problem == 'AP2' ) then 
            
        if ( ind == 1 ) then
            H(1,1) = 2.0d0
            return
        end if
        
        if ( ind == 2 ) then
            H(1,1) = 2.0d0
            return
        end if	
            
    end if		
    
    
! ----------------------------------------------------------------------

! AP3: Exemple 3 of "A modified Quasi-Newton method for vector optimization problem"

    if ( problem == 'AP3' ) then 
            
        if ( ind == 1 ) then
            !g(1) = ( x(1) - 1.0d0 ) ** 3 
            !g(2) = 2.0d0 * ( x(2) - 2.0d0 ) ** 3
            
            H(:,:) = 0.0d0
            H(1,1) = 3.0d0 * ( x(1) - 1.0d0 ) ** 2
            H(2,2) = 6.0d0 * ( x(2) - 2.0d0 ) ** 2
            return
        end if
        
        if ( ind == 2 ) then
!			g(1) = - 4.0d0 * x(1) * ( x(2) - x(1) ** 2 ) - 2.0d0 * ( 1.0d0 - x(1) )
!			g(2) = 2.0d0 * ( x(2) - x(1) ** 2 )

            H(1,1) = - 4.0d0  * ( x(2) - x(1) ** 2 ) + 8.0d0 * x(1) ** 2 + 2.0d0
            H(1,2) = - 4.0d0 * x(1)
            H(2,1) = H(1,2)
            H(2,2) = 2.0d0
            return
        end if	
            
    end if	
    
! ----------------------------------------------------------------------

! AP4: Exemple 4 of "A modified Quasi-Newton method for vector optimization problem"

    if ( problem == 'AP4' ) then 
        
        if ( ind == 1 ) then
            H(:,:) = 0.0d0
            H(1,1) = 1.2d1/9.0d0 * ( x(1) - 1.0d0 ) ** 2
            H(2,2) = 2.4d1/9.0d0 * ( x(2) - 2.0d0 ) ** 2 
            H(3,3) = 3.6d1/9.0d0 * ( x(3) - 3.0d0 ) ** 2
            return
        end if
        
        if ( ind == 2 ) then
            t = 1.0d0/9.0d0 * exp( ( x(1) + x(2) + x(3) ) / 3.0d0 )
            H(:,:) = t
            forall (i = 1:n) H(i,i) = H(i,i) + 2.0d0
            return
        end if	
        
        if ( ind == 3 ) then
            H(:,:) = 0.0d0
            H(1,1) = 1.0d0/4.0d0 * exp( -x(1) ) 
            H(2,2) = 1.0d0/3.0d0 * exp( -x(2) ) 
            H(3,3) = 1.0d0/4.0d0 * exp( -x(3) ) 
            return
        end if	
            
    end if	
    
! ----------------------------------------------------------------------	

    !  BK1

    if ( problem == 'BK1' ) then 
            
        if ( ind == 1 ) then
            H(:,:) = 0.0d0
            H(1,1) = 2.0d0 
            H(2,2) = 2.0d0
            return
        end if
        
        if ( ind == 2 ) then
            H(:,:) = 0.0d0
            H(1,1) = 2.0d0 
            H(2,2) = 2.0d0
            return	
        end if
            
    end if		
    
  end subroutine evalh

end module myproblem