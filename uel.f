 ************************************************************************
!
! User element for large elastic deformation in 2D using 4 node
! isoparametric element with full (4pt) gauss integration. This is for
! plane strain only.
!
! Displacements and chemical potential are the nodal variables
! 
!     
!              A eta (=xi_2)
!  4-node      |
!   quad       |Face 3
!        4-----------3
!        |     |     |
!        |     |     |
!  Face 4|     ------|---> xi (=xi_1)
!        |           | Face2
!        |           |
!        1-----------2
!          Face 1
!
***********************************************************************
!
!     State Variables
!     --------------------------------------------------------------
!     Global SDV's (used for visualization)
!       1) polymer volume fraction (phi)
!       2) concentration (c)
!
!     Material Properties Vector
!     --------------------------------------------------------------
!     Gshear = props(1) ! Shear modulus
!     Kbulk  = props(2) ! Bulk modulus
!     chi    = props(3) ! Chi parameter
!     D      = props(4) ! Coefficient of permeability
!     mu0    = props(5) ! Chemical potential of pure fluid
!     Vmol   = props(6) ! Volume of a mole of fluid particles
!     Rgas   = props(7) ! Universal gas constant
!     theta  = props(8) ! Absolute temperature
!     phi0   = props(9) ! Initial polymer volume fraction
!     nlSdv  = jprops(1) ! Number of local sdv's per integ pt
!     ngSdv  = jprops(2) ! Number of global sdv's per integ pt
!
!***********************************************************************

      module global

      integer numElem,ElemOffset,err

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Set the number of UEL elements used here, must be greater than
      ! or equal to the value in input file
      parameter(numElem=10000)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Set the offset here for UVARM plotting, must match input file!
      parameter(ElemOffset=10000)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      real*8, allocatable :: globalSdv(:,:,:)

      end module global

***********************************************************************

      SUBROUTINE UVARM(UVAR,DIRECT,T,TIME,DTIME,CMNAME,ORNAME,
     1 NUVARM,NOEL,NPT,LAYER,KSPT,KSTEP,KINC,NDI,NSHR,COORD,
     2 JMAC,JMATYP,MATLAYO,LACCFLA)
     
      use global
     
      include 'ABA_PARAM.INC'

      CHARACTER*80 CMNAME,ORNAME
      CHARACTER*3 FLGRAY(15)
      DIMENSION UVAR(NUVARM),DIRECT(3,3),T(3,3),TIME(2)
      DIMENSION ARRAY(15),JARRAY(15),JMAC(*),JMATYP(*),COORD(*)

      uvar(1) = globalSdv(noel-ElemOffset,npt,1)
      uvar(2) = ((1/uvar(1))-1)/(1.8d-5)

      return
      end subroutine uvarm

****************************************************************************

      SUBROUTINE UEL(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     +     PROPS,NPROPS,coords,MCRD,NNODE,Uall,DUall,Vel,Accn,JTYPE,
     +     TIME,DTIME,KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,
     +     PREDEF,NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,
     +     NJPROP,PERIOD)

      use global
*
      IMPLICIT NONE
*
*     VARIABLES DEFINED IN UEL, PASSED BACK TO ABAQUS
*
      REAL(8) :: RHS,AMATRX,SVARS,ENERGY
*
*     VARIABLES PASSED INTO UEL 
*
      REAL(8) :: PROPS,coords,Uall,DUall,Vel,Accn,TIME,
     1  DTIME,PARAMS,ADLMAG,PREDEF,DDLMAG,PNEWDT,PERIOD
      INTEGER :: NDOFEL,NRHS,NSVARS,NPROPS,MCRD,NNODE,JTYPE,KSTEP,KINC,
     1  JELEM,NDLOAD,JDLTYP,NPREDF,LFLAGS,MLVARX,MDLOAD,JPROPS,NJPROP
*
      DIMENSION RHS(MLVARX,*),AMATRX(NDOFEL,NDOFEL),PROPS(*),
     1 SVARS(*),ENERGY(8),coords(MCRD,NNODE),Uall(NDOFEL),
     2 DUall(MLVARX,*),Vel(NDOFEL),Accn(NDOFEL),TIME(2),PARAMS(*),
     3 JDLTYP(MDLOAD,*),ADLMAG(MDLOAD,*),DDLMAG(MDLOAD,*),
     4 PREDEF(2,NPREDF,NNODE),LFLAGS(*),JPROPS(*)

      integer lenJobName,lenOutDir,nDim,nInt,nIntS
      character*256 jobName,outDir,fileName



      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      parameter(nInt=4)  ! number of volume integration pionts
      parameter(nIntS=1) ! number of surface integration points
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




      !----------------------------------------------------------------
      ! 
      ! Perform initial checks
      !
      !
      ! Open the debug/error message file
      !
      call getJobName(jobName,lenJobName)
      call getOutDir(outDir,lenOutDir)
      fileName = outDir(1:lenOutDir)//'\aaMSGS_'//
     +     jobName(1:lenJobName)//'.dat'
      open(unit=80,file=fileName,status='unknown')


      ! Checking the procedure type
      !
      if((lflags(1).eq.64).or.(lflags(1).eq.65).or.
     +     (lflags(1).eq.72).or.(lflags(1).eq.73)) then
         !
         ! all is good
         !
      else
         write(*,*) 'Abaqus does not have the right procedure'
         write(*,*) 'go back and chekc the procedure type'
         write(*,*) 'lflags(1)=',lflags(1)
         write(80,*) 'Abaqus does not have the right procedure'
         write(80,*) 'go back and chekc the procedure type'
         write(80,*) 'lflags(1)=',lflags(1)
         call xit
      endif


      ! Checking large deformation
      !
      if(lflags(2).eq.0) then
         !
         ! lflags(2)=0 -> small disp.
         ! lflags(2)=1 -> large disp.
         !
         write(*,*) 'Abaqus thinks you are doing'
         write(*,*) 'a small displacement analysis'
         write(*,*) 'go in and set nlgeom=yes'
         write(80,*) 'Abaqus thinks you are doing'
         write(80,*) 'a small displacement analysis'
         write(80,*) 'go in and set nlgeom=yes'
         call xit
      endif


      ! Check to see if you are doing a general
      !  step or a linear purturbation step
      !
      if(lflags(4).eq.1) then
         !
         ! lflags(4)=0 -> general step
         ! lflags(4)=1 -> linear purturbation step
         !
         write(*,*) 'Abaqus thinks you are doing'
         write(*,*) 'a linear purturbation step'
         write(80,*) 'Abaqus thinks you are doing'
         write(80,*) 'a linear purturbation step'
         call xit         
      endif


      ! Do nothing if a ``dummy'' step
      !
      if(dtime.eq.0.0) return
      !
      ! Done with initial checks
      !
      !----------------------------------------------------------------



      !----------------------------------------------------------------
      ! 
      ! Call the paricular element to perform the analysis
      !
      if(jtype.eq.1) then
         nDim = 2
         call UPE41(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     +        PROPS,NPROPS,coords,MCRD,NNODE,Uall,DUall,Vel,Accn,JTYPE,
     +        TIME,DTIME,KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,
     +        PREDEF,NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,
     +        NJPROP,PERIOD,
     +        nDim,nInt,nIntS)         !
         !
      else
         !
         ! We have a problem...
         !
         write(*,*) 'Element type not supported, jtype=',jtype
         write(80,*) 'Element type not supported, jtype=',jtype
         call xit
         !
      endif
      !----------------------------------------------------------------


      return
      end subroutine uel

************************************************************************

      subroutine UPE41(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     +     PROPS,NPROPS,coords,MCRD,NNODE,Uall,DUall,Vel,Accn,JTYPE,
     +     TIME,DTIME,KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,
     +     PREDEF,NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,
     +     NJPROP,PERIOD,
     +     nDim,nInt,nIntS)

      use global
*
      IMPLICIT NONE
*
*     VARIABLES DEFINED IN UEL, PASSED BACK TO ABAQUS
*
      REAL(8) :: RHS,AMATRX,SVARS,ENERGY
*
*     VARIABLES PASSED INTO UEL 
*
      REAL(8) :: PROPS,coords,Uall,DUall,Vel,Accn,TIME,
     1  DTIME,PARAMS,ADLMAG,PREDEF,DDLMAG,PNEWDT,PERIOD
      INTEGER :: NDOFEL,NRHS,NSVARS,NPROPS,MCRD,NNODE,JTYPE,KSTEP,KINC,
     1  JELEM,NDLOAD,JDLTYP,NPREDF,LFLAGS,MLVARX,MDLOAD,JPROPS,NJPROP
*
      DIMENSION RHS(MLVARX,*),AMATRX(NDOFEL,NDOFEL),PROPS(*),
     1 SVARS(*),ENERGY(8),coords(MCRD,NNODE),Uall(NDOFEL),
     2 DUall(MLVARX,*),Vel(NDOFEL),Accn(NDOFEL),TIME(2),PARAMS(*),
     3 JDLTYP(MDLOAD,*),ADLMAG(MDLOAD,*),DDLMAG(MDLOAD,*),
     4 PREDEF(2,NPREDF,NNODE),LFLAGS(*),JPROPS(*)

      real*8 u(nNode,2),du(nNode,ndofel),thetaNew(nNode)
      real*8 thetaOld(nNode),dtheta(nNode),muNew(nNode)
      real*8 muOld(nNode),dMU(nNode),uNew(nNode,ndofel)
      real*8 uOld(nNode,ndofel),u_t(nNode,ndofel),v(nNode,2)
      real*8 coordsC(mcrd,nNode)

      integer i,j,k,l,m,n,nIntPt,nDim,intpt,pOrder,a1,b1,a11,b11,face
      integer nInt,ii,jj,a12,b12,pe,nSdv,stat,q,nIntV,nIntPtV,p,ngSdv
      integer nlSdv,kk,lenJobName,lenOutDir,faceFlag,nIntS

      real*8 Iden(3,3),Le,theta0,phi0,Ru(2*nNode,1),Rc(nNode,1),body(3)
      real*8 Kuu(2*nNode,2*nNode),Kcc(nNode,nNode),sh0(nNode),detMapJ0
      real*8 dshxi(nNode,2),dsh0(nNode,2),dshC0(nNode,2),detMapJ0C,Vmol
      real*8 Fc_tau(3,3),Fc_t(3,3),detFc_tau,detFc_t,w(nInt),ds,flux
      real*8 sh(nNode),detMapJ,phi_t,dsh(nNode,2),detMapJC,phiLmt,umeror
      real*8 dshC(nNode,2),mu_tau,mu_t,dMUdX(2,1),dMUdt,F_tau(3,3),DmDJ
      real*8 F_t(3,3),detF_tau,xi(nInt,2),detF,TR_tau(3,3),T_tau(3,3)
      real*8 SpTanMod(3,3,3,3),phi_tau,dPdt,DphiDmu,DphidotDmu,Mfluid
      real*8 Smat(3,1),Bmat(3,2*nNode),BodyForceRes(2*nNode,1),Qmat(4,4)
      real*8 DmDmu,Gmat(4,2*nNode),G0mat(4,2*nNode),Amat(4,4),wS(nIntS)
      real*8 xLocal(nIntS),yLocal(nIntS),Kcu(nNode,2*nNode),detF_t
      real*8 Kuc(2*nNode,nNode),Nvec(1,nNode),ResFac,AmatUC(3,1),TanFac
      real*8 SpUCMod(3,3),SpCUMod(3,3,3),SpCUModFac(3,3),AmatCU(2,4)

      real*8 zero,one,two,half,Pi,three,third,pow
      parameter(zero=0.d0,one=1.d0,two=2.d0,half=0.5d0,Pi=3.141592653d0,
     +     three=3.d0,third=1.d0/3.d0)


      ! Get element parameters
      !
      nlSdv = jprops(1) ! number of local sdv's per integ point
      ngSdv = jprops(2) ! number of global sdv's per integ point


      ! Allocate memory for the globalSdv's used for viewing
      !  results on the dummy mesh
      !
      if(.not.allocated(globalSdv)) then
         !
         ! allocate memory for the globalSdv's
         !
         ! numElem needs to be set in the MODULE
         ! nInt needs to be set in the UEL
         !
         stat=0
c         allocate(globalSdv(numElem,nInt,ngSdv))
c         deallocate(globalSdv)
         allocate(globalSdv(numElem,nInt,ngSdv),stat=err)
         if(stat.ne.0) then
            write(*,*) '//////////////////////////////////////////////'
            write(*,*) 'error when allocating globalSdv'
            write(*,*) '//////////////////////////////////////////////'
            write(*,*) '   stat=',stat
            write(*,*) '  ngSdv=',ngSdv
            write(*,*) '   nInt=',nInt
            write(*,*) 'numElem=',numElem
            write(*,*) '  nNode=',nNode
            write(*,*) 'lbound(globalSdv)',lbound(globalSdv)
            write(*,*) 'ubound(globalSdv)',ubound(globalSdv)
            write(*,*) '//////////////////////////////////////////////'
            write(80,*) '//////////////////////////////////////////////'
            write(80,*) 'error when allocating globalSdv'
            write(80,*) '//////////////////////////////////////////////'
            write(80,*) '   stat=',stat
            write(80,*) '  ngSdv=',ngSdv
            write(80,*) '   nInt=',nInt
            write(80,*) 'numElem=',numElem
            write(80,*) '  nNode=',nNode
            write(80,*) 'lbound(globalSdv)=',lbound(globalSdv)
            write(80,*) 'ubound(globalSdv)=',ubound(globalSdv)
            write(80,*) '//////////////////////////////////////////////'
            call xit
         endif
         write(*,*) '-------------------------------------------------'
         write(*,*) '----------- globalSDV ALLOCATED -----------------'
         write(*,*) '-------------------------------------------------'
         write(*,*) '---------- YOU PUT NUMBER OF ELEMENTS -----------'
         write(*,*) '---------- numElem=',numElem
         write(*,*) '---------- UPE4 ELEMENTS ------------------------'
         write(*,*) '-------------------------------------------------'
         write(*,*) '---------- YOU PUT NUMBER OF POINTS -------------'
         write(*,*) '---------- nInt =',nInt
         write(*,*) '---------- nIntS=',nIntS
         write(*,*) '-------------------------------------------------'
         write(*,*) '---------- YOU PUT NUMBER OF SDVs ---------------'
         write(*,*) '---------- ngSdv=',ngSdv
         write(*,*) '-------------------------------------------------'
      endif


      ! Identity tensor
      !
      call onem(Iden)


      ! Obtain initial conditions
      !
      theta0 = props(8)
      phi0   = props(9)


      ! Initialize the residual and tangent matrices to zero.
      !
      Ru  = zero
      Rc = zero
      Kuu = zero
      Kcc = zero
      Kuc = zero
      Kcu = zero
      Energy = zero


      ! Body forces
      !
      body(1:3) = zero


      ! Obtain nodal displacements and chemical potentials
      !
      k = 0
      do i=1,nNode
         do j=1,nDim
            k = k + 1
            u(i,j) = Uall(k)
            du(i,j) = DUall(k,1)
            uOld(i,j) = u(i,j) - du(i,j)
         enddo
         k = k + 1
         muNew(i) = Uall(k)
         dMU(i) = DUall(k,1)
         muOld(i) = muNew(i) - dMU(i)
      enddo


      ! Obtain current nodal coordinates
      !
      do i=1,nNode
         do j=1,nDim
            coordsC(j,i) = coords(j,i) + u(i,j)
         enddo
      enddo


      ! Impose any time-stepping changes on the increments of
      !  chemical potential or displacement if you wish
      !
      ! chemical potential increment
      !
      do i=1,nNode
         if(dabs(dMU(i)).gt.1.d6) then
            pnewdt = 0.5
            return
         endif
      enddo
      !
      ! displacement increment, based on element diagonal
      !
      Le = dsqrt(((coordsC(1,1)-coordsC(1,3))**two) + 
     +     ((coordsC(2,1)-coordsC(2,3))**two))
      !
      do i=1,nNode
         do j=1,nDim
            if(dabs(du(i,j)).gt.10.0*Le) then
               pnewdt = 0.5
               return
            endif
         enddo
      enddo


      !----------------------------------------------------------------
      !
      ! Obtain shape functions and their local gradients at the element
      !  centriod, that means xi=eta=zeta=0.0, and nIntPt=1
      !
      if(nNode.eq.4) then
         call calcShape2DLinear(1,zero,1,sh0,dshxi)
      else
         write(*,*) 'Incorrect number of nodes: nNode.ne.4'
         write(80,*) 'Incorrect number of nodes: nNode.ne.4'
         call xit
      endif


      ! Map shape functions from local to global reference coordinate system
      !
      if(mcrd.eq.2) then
         call mapShape2Da(nNode,dshxi,coords,dsh0,detMapJ0,stat)
         if(stat.eq.0) then
            pnewdt = 0.5
            return
         endif
      elseif(mcrd.eq.3) then
         call mapShape2D(nNode,dshxi,coords,dsh0,detMapJ0,stat)
         if(stat.eq.0) then
            pnewdt = 0.5
            return
         endif
      else
         ! big problem
         write(*,*) 'Unexpected error, mcrd=',mcrd
         write(80,*) 'Unexpected error, mcrd=',mcrd
         call xit
      endif


      ! Map shape functions from local to global current coordinate system
      !
      if(mcrd.eq.2) then
         call mapShape2Da(nNode,dshxi,coordsC,dshC0,detMapJ0C,stat)
         if(stat.eq.0) then
            pnewdt = 0.5
            return
         endif
      elseif(mcrd.eq.3) then
         call mapShape2D(nNode,dshxi,coordsC,dshC0,detMapJ0C,stat)
         if(stat.eq.0) then
            pnewdt = 0.5
            return
         endif
      else
         ! big problem
         write(*,*) 'Unexpected error, mcrd=',mcrd
         write(80,*) 'Unexpected error, mcrd=',mcrd
         call xit
      endif

      ! Calculate the deformation gradient at the element centriod
      !
      Fc_tau = Iden
      Fc_t = Iden
      do i=1,nDim
         do j=1,nDim
            do k=1,nNode
               Fc_tau(i,j) = Fc_tau(i,j) + dsh0(k,j)*u(k,i)
               Fc_t(i,j) = Fc_t(i,j) + dsh0(k,j)*uOld(k,i)
            enddo
         enddo
      enddo
      !
      ! modify for plane-strain
      !
      Fc_tau(3,3) = one
      Fc_t(3,3) = one
      !
      ! 2D plane-strain implementation detF
      !
      detFc_t = Fc_t(1,1)*Fc_t(2,2) - Fc_t(1,2)*Fc_t(2,1)
      detFc_tau = Fc_tau(1,1)*Fc_tau(2,2) - Fc_tau(1,2)*Fc_tau(2,1)
      !
      !----------------------------------------------------------------




      !----------------------------------------------------------------
      ! Begin the loop over body integration points
      !
      ! Obtain integration point local coordinates and weights
      !
      if(nNode.eq.4) then
         !
         ! gauss integration for a rectangular element
         !
         if(nInt.eq.4) then
            call xint2D4pt(xi,w,nIntPt) ! 4-pt integration, nInt=4 above
         else
            write(*,*) 'Invalid number of int points, nInt=',nInt
            write(80,*) 'Invalid number of int points, nInt=',nInt
            call xit
         endif
      else
         write(*,*) 'Incorrect number of nodes: nNode.ne.4'
         write(80,*) 'Incorrect number of nodes: nNode.ne.4'
         call xit
      endif


      ! Loop over integration points
      !
      jj = 0 ! jj is used for tracking the state variables
      do intpt=1,nIntPt


         ! Obtain state variables from previous increment
         !
         if((kinc.le.1).and.(kstep.eq.1)) then
            !
            ! this is the first increment, of the first step
            !  give initial conditions
            !
            phi_t  = phi0
            !
         else
            !
            ! this is not the first increment, read old values
            !
            phi_t  = svars(1+jj)
            !
         endif


         ! Obtain shape functions and their local gradients
         !
         if(nNode.eq.4) then
            call calcShape2DLinear(nIntPt,xi,intpt,sh,dshxi)
         else
            write(*,*) 'Incorrect number of nodes: nNode.ne.4'
            write(80,*) 'Incorrect number of nodes: nNode.ne.4'
            call xit
         endif
         

         ! Map shape functions from local to global reference coordinate system
         !
         if(mcrd.eq.2) then
            call mapShape2Da(nNode,dshxi,coords,dsh,detMapJ,stat)
            if(stat.eq.0) then
               pnewdt = 0.5
               return
            endif
         elseif(mcrd.eq.3) then
            call mapShape2D(nNode,dshxi,coords,dsh,detMapJ,stat)
            if(stat.eq.0) then
               pnewdt = 0.5
               return
            endif
         else
            ! big problem
            write(*,*) 'Unexpected error, mcrd=',mcrd
            write(80,*) 'Unexpected error, mcrd=',mcrd
            call xit
         endif


         ! Map shape functions from local to global current coordinate system
         !
         if(mcrd.eq.2) then
            call mapShape2Da(nNode,dshxi,coordsC,dshC,detMapJC,stat)
            if(stat.eq.0) then
               pnewdt = 0.5
               return
            endif
         elseif(mcrd.eq.3) then
            call mapShape2D(nNode,dshxi,coordsC,dshC,detMapJC,stat)
            if(stat.eq.0) then
               pnewdt = 0.5
               return
            endif
         else
            ! big problem
            write(*,*) 'Unexpected error, mcrd=',mcrd
            write(80,*) 'Unexpected error, mcrd=',mcrd
            call xit
         endif


         ! Obtain the chemical potential and its derivative's at 
         !  this intPt at the begining and end of the incrment
         !
         mu_tau = zero
         mu_t = zero
         dMUdt = zero
         dMUdX = zero
         do k=1,nNode
            mu_tau = mu_tau + muNew(k)*sh(k)
            mu_t   = mu_t + muOld(k)*sh(k)
            if(ndload.gt.0) then
              if (u(k,2).gt.1.d-3) then
                pow = -u(k,2)*1.d3
              else
                pow = 0
              endif
            else
              pow = 0
            endif
            do i=1,nDim
              !
              ! Implement cut-off
              !
              dMUdX(i,1) = dMUdX(i,1) + muNew(k)*dshC(k,i)*exp(pow)
            enddo
         enddo
         dMUdt = (mu_tau - mu_t)/dtime


         ! Obtain, and modify the deformation gradient at this integration
         !  point.
         !
         F_tau = Iden
         F_t = Iden
         do i=1,nDim
            do j=1,nDim
               do k=1,nNode
                  F_tau(i,j) = F_tau(i,j) + dsh(k,j)*u(k,i)
                  F_t(i,j) = F_t(i,j) + dsh(k,j)*uOld(k,i)
               enddo
            enddo
         enddo
         !
         ! modify F(3,3) for plane-strain 
         !
         F_tau(3,3) = one
         F_t(3,3) = one
         !
         call mdet(F_tau,detF)


         !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
         !
         ! Perform the constitutive time integration at this integ. point 
         !
         call integ1(props,nprops,dtime,
     +        F_tau,mu_tau,phi_t,theta0,
     +        T_tau,SpTanMod,
     +        phi_tau,dPdt,DphiDmu,DphidotDmu,
     +        Mfluid,DmDmu,DmDJ,Vmol,
     +        SpUCMod,SpCUModFac,jelem,intPt)
         !
         !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

         ! Save the state variables at this integ point
         !  at the end of the increment
         !
         svars(1+jj) = phi_tau
         jj = jj + nlSdv ! setup for the next intPt


         ! Save the state variables at this integ point in the
         !  global array used for plotting field output
         !
         globalSdv(jelem,intPt,1) = phi_tau   ! polymer volume fraction


         ! Time stepping algorithim based on the constitutive response.
         !  Here based on the change in the polymer volume fraction change.
         !
         phiLmt = 0.005d0
         umeror = dabs((phi_tau - phi_t)/phiLmt)
         if(umeror.le.half) then
            pnewdt = 1.5d0
         elseif((umeror.gt.half).and.(umeror.le.0.8d0)) then
            pnewdt = 1.25d0
         elseif((umeror.gt.0.8d0).and.(umeror.le.1.25d0)) then
            pnewdt = 0.75d0 
         else
            pnewdt = half
         endif


         ! Compute/update the displacement residual vector
         !
         Smat(1,1) = T_tau(1,1)
         Smat(2,1) = T_tau(2,2)
         Smat(3,1) = T_tau(1,2)
         !
         Bmat = zero
         do kk=1,nNode
            Bmat(1,1+nDim*(kk-1)) = dshC(kk,1)
            Bmat(2,2+nDim*(kk-1)) = dshC(kk,2)
            Bmat(3,1+nDim*(kk-1)) = dshC(kk,2)
            Bmat(3,2+nDim*(kk-1)) = dshC(kk,1)
         enddo
         !
         BodyForceRes = zero
         do kk=1,nNode
            BodyForceRes(1+nDim*(kk-1),1) = sh(kk)*body(1)
            BodyForceRes(2+nDim*(kk-1),1) = sh(kk)*body(2)
         enddo
         !
         Ru = Ru + detmapJC*w(intpt)*
     +        (
     +        -matmul(transpose(Bmat),Smat)
     +        + BodyForceRes
     +        )



         ! Compute/update the chemical potential residual vector
         !
         do kk=1,nNode
            Nvec(1,kk) = sh(kk)
         enddo
         !
         ResFac = (dPdt)/(detF*Vmol*phi_tau*phi_tau)
         !
         Rc = Rc + detmapJC*w(intpt)*
     +        (
     +        transpose(Nvec)*ResFac - Mfluid*matmul(dshC,dMUdX)
     +        )



         ! Compute/update the displacement tangent matrix
         !
         Gmat = zero
         do kk=1,nNode
            Gmat(1,1+nDim*(kk-1)) = dshC(kk,1)
            Gmat(2,2+nDim*(kk-1)) = dshC(kk,1)
            Gmat(3,1+nDim*(kk-1)) = dshC(kk,2)
            Gmat(4,2+nDim*(kk-1)) = dshC(kk,2)
         enddo

         G0mat = zero
         do kk=1,nNode
            G0mat(1,1+nDim*(kk-1)) = dshC0(kk,1)
            G0mat(2,2+nDim*(kk-1)) = dshC0(kk,1)
            G0mat(3,1+nDim*(kk-1)) = dshC0(kk,2)
            G0mat(4,2+nDim*(kk-1)) = dshC0(kk,2)
         enddo

         Amat = zero
         Amat(1,1) = SpTanMod(1,1,1,1)
         Amat(1,2) = SpTanMod(1,1,2,1)
         Amat(1,3) = SpTanMod(1,1,1,2)
         Amat(1,4) = SpTanMod(1,1,2,2)
         Amat(2,1) = SpTanMod(2,1,1,1)
         Amat(2,2) = SpTanMod(2,1,2,1)
         Amat(2,3) = SpTanMod(2,1,1,2)
         Amat(2,4) = SpTanMod(2,1,2,2)
         Amat(3,1) = SpTanMod(1,2,1,1)
         Amat(3,2) = SpTanMod(1,2,2,1)
         Amat(3,3) = SpTanMod(1,2,1,2)
         Amat(3,4) = SpTanMod(1,2,2,2)
         Amat(4,1) = SpTanMod(2,2,1,1)
         Amat(4,2) = SpTanMod(2,2,2,1)
         Amat(4,3) = SpTanMod(2,2,1,2)
         Amat(4,4) = SpTanMod(2,2,2,2)



         Qmat = zero
         Qmat(1,1) = half*(Amat(1,1)+Amat(1,4)) - half*T_tau(1,1)
         Qmat(2,1) = half*(Amat(2,1)+Amat(2,4)) - half*T_tau(1,2)
         Qmat(3,1) = half*(Amat(3,1)+Amat(3,4)) - half*T_tau(1,2)
         Qmat(4,1) = half*(Amat(4,1)+Amat(4,4)) - half*T_tau(2,2)
         Qmat(1,4) = half*(Amat(1,1)+Amat(1,4)) - half*T_tau(1,1)
         Qmat(2,4) = half*(Amat(2,1)+Amat(2,4)) - half*T_tau(1,2)
         Qmat(3,4) = half*(Amat(3,1)+Amat(3,4)) - half*T_tau(1,2)
         Qmat(4,4) = half*(Amat(4,1)+Amat(4,4)) - half*T_tau(2,2)
            

         ! Tangent modulus
         Kuu = Kuu + detMapJC*w(intpt)*
     +           (
     +           matmul(matmul(transpose(Gmat),Amat),Gmat)
     +           )



         ! Compute/update the chemical potential tangent matrix
         !
         TanFac = (one/(detF*Vmol*phi_tau**two))*
     +        (two*(dPdt/phi_tau)*DphiDmu - DphidotDmu)
         !
         Kcc = Kcc + detmapJC*w(intPt)*
     +        (
     +        TanFac*matmul(transpose(Nvec),Nvec)
     +        + Mfluid*matmul(dshC,transpose(dshC))
     +        + DmDmu*matmul(matmul(dshC,dMUdX),Nvec)
     +        )

         

         ! Compute/update the chemical potential - displacement tangent matrix
         !
         SpCUMod = zero
         do i=1,nDim
            do k=1,nDim
               do l=1,nDim
                  SpCUMod(i,k,l) = SpCUMod(i,k,l)
     +                 + dMUdX(k,1)*SpCUModFac(i,l)
               enddo
            enddo
         enddo
         !
         AmatCU = zero
         AmatCU(1,1) = SpCUMod(1,1,1)
         AmatCU(1,2) = SpCUMod(1,2,1)
         AmatCU(1,3) = SpCUMod(1,1,2)
         AmatCU(1,4) = SpCUMod(1,2,2)
         AmatCU(2,1) = SpCUMod(2,1,1)
         AmatCU(2,2) = SpCUMod(2,2,1)
         AmatCU(2,3) = SpCUMod(2,1,2)
         AmatCU(2,4) = SpCUMod(2,2,2)
         !
         Kcu = Kcu - detMapJC*w(intpt)*
     +        (
     +        matmul(matmul(dshC,AmatCU),Gmat)
     +        )


         ! Compute/update the displacement - chemical potential tangent matrix
         !
         AmatUC = zero
         AmatUC(1,1) = SpUCMod(1,1)
         AmatUC(2,1) = SpUCMod(2,2)
         AmatUC(3,1) = SpUCMod(1,2)
         !
         Kuc = Kuc + detMapJC*w(intpt)*
     +        (
     +        matmul(matmul(transpose(Bmat),AmatUC),Nvec)
     +        )

      enddo
      !
      ! End the loop over body integration points
      !----------------------------------------------------------------




      !----------------------------------------------------------------
      ! Return Abaqus the RHS vector and the Stiffness matrix.
      !
      call AssembleElement(nDim,nNode,nDofEl,
     +     Ru,Rc,Kuu,Kuc,Kcu,Kcc,
     +     rhs,amatrx)
      !
      ! End return of RHS and AMATRX
      !----------------------------------------------------------------


      return
      end subroutine UPE41

************************************************************************

      subroutine integ1(props,nprops,dtime,
     +     F_tau,mu_tau,phi_t,theta,
     +     T_tau,SpTanMod,
     +     phi_tau,dPdt,DphiDmu,DphidotDmu,
     +     Mfluid,DmDmu,DmDJ,Vmol,
     +     SpUCMod,SpCUModFac,jelem,intPt)

      ! This subroutine computes everything required for the time integration
      ! of the problem.
      !
      ! Inputs:
      !  1) material parameters, props(nprops)
      !  2) time increment, dtime
      !  3) deformation gradient, F_tau(3,3)
      !  4) chemical potential, mu_tau
      !  5) old polymer volume fraction, phi_t
      !  6) temperature, theta
      !
      ! Outputs:
      !  1) Cauchy stress, T_tau(3,3)
      !  2) spatial tangent modulus, SpTanMod(3,3,3,3)
      !  3) polymer volume fraction, phi_tau
      !  4) time rate of polymer volume fraction, dPdt
      !  5) derivative of the phi with mu, DphiDmu
      !  6) derivative of the time rate of phi with mu, DphidotDmu
      !  7) scalar fluid permeability, Mfluid
      !  8) derivative of permeability with chemical potential, DmDmu
      !  9) volume of a mole of fluid, Vmol
      ! 10) displacement - chemical potential modulus terms
      ! 11) chemical potential - displacement modulus terms

      use global

      implicit none

      integer i,j,k,l,m,n,nprops,nargs,stat,jelem,intPt,ee
      parameter(nargs=8)

      real*8 Iden(3,3),props(nprops),F_tau(3,3),phi_tau,mu_tau,phi_t
      real*8 theta,TR_tau(3,3),T_tau(3,3),dTRdF(3,3,3,3),Gshear,Kbulk
      real*8 spTanMod(3,3,3,3),chi,D,mu0,Vmol,Rgas,detF,FinvT(3,3),dPdt
      real*8 B_tau(3,3),trB_tau,C_tau(3,3),trC_tau,args(nargs),detFe
      real*8 deltaMU,DphiDmu,dPdt_per,dPdt_m,DphidotDmu,Mfluid,Finv(3,3)
      real*8 phi_per,phi_m,dtime,DmDmu,DphiDJ,SpUCMod(3,3),DmDphi,DmDJ
      real*8 SpCUModFac(3,3),detFs,ec,RT,kp

      real*8 zero,one,two,three,third,half
      parameter(zero=0.d0,one=1.d0,two=2.d0,three=3.d0,third=1.d0/3.d0,
     +     half=1.d0/2.d0)


      ! Identity tensor
      !
      call onem(Iden)
      !
      ! Incorporate effect of concentration on Gshear, Kbulk
      !
      if(globalSdv(jelem,intPt,1).gt.0.65) then
        kp = globalSdv(jelem,intPt,1)
      else
        kp = 0.65
      endif
      ec = (-0.867*dlog((1-kp)/
     +       (kp*props(6)))+9.182)


      ! Obtain material properties
      !
      Gshear = props(1)*ec*1.d5
      Kbulk  = props(2)*ec*1.d5
      chi    = props(3)
      D      = 1.7d-11
      mu0    = props(5)
      Vmol   = props(6)
      Rgas   = props(7)
      !print*, Gshear, Kbulk


      ! Compute the inverse of F, its determinant, and its transpose
      !
      call matInv3D(F_tau,Finv,detF,stat)
      if(stat.eq.0) then
         write(*,*) 'Problem: detF.lt.zero'
         call xit
      endif
      FinvT = transpose(Finv)


      ! Compute the left Cauchy-Green tensor and its trace
      !
      B_tau = matmul(F_tau,transpose(F_tau))
      trB_tau = B_tau(1,1) + B_tau(2,2) + B_tau(3,3)


      ! Compute the right Cauchy-Green tensor and its trace
      !
      C_tau = matmul(transpose(F_tau),F_tau)
      trC_tau = C_tau(1,1) + C_tau(2,2) + C_tau(3,3)


      ! Compute the polymer volume fraction
      !
      args(1)  = mu_tau
      args(2)  = mu0
      args(3)  = Rgas
      args(4)  = theta
      args(5)  = chi
      args(6)  = Vmol
      args(7)  = Kbulk
      args(8)  = detF
      args(9)  = Gshear
      args(10) = props(1)
      args(11) = props(2)
      call solvePhi1(phi_tau,args,nargs,phi_t,F_tau)


      ! Compute the elastic volume ratio, detFe
      !
      detFe = detF*phi_tau


      ! Compute the swelling volume ratio, detFs
      !
      detFs = one/phi_tau


      ! Compute the time rate of the polymer volume fraction using
      !  a finite difference in time
      !
      dPdt = (phi_tau - phi_t)/dtime

      ! Useful
      RT = Rgas*theta


      ! Compute the derivative of the polymer volume fraction with
      !  respect to the chemical potential.  Computed via implicit
      !  differentiation on the chemical potential equation.
      !
      DphiDmu =  (one/RT)/
     +     ((one/(phi_tau - one) + one
     +     + two*chi*phi_tau) + half*(trC_tau - three 
     +     - two*dlog(detF))*props(1)*(-0.867)*Vmol/(((one
     +     -phi_tau)**2)*RT)+ props(2)*(-0.867)*(Vmol**2)*
     +     dlog(detF*phi_tau)*((dlog(detF*phi_tau)-2)*(
     +     -(phi_tau**2)+phi_tau+one)*half
     +     +one/phi_tau)/(((1-phi_tau)**2)*RT) 
     +     +Kbulk*Vmol*(dlog(detF*phi_tau)-one)/(RT*phi_tau))



      ! Compute the derivative of the polymer volume fraction with
      !  respect to the chemical potential.  Computed via implicit
      !  differentiation on the chemical potential equation.
      !
      DphiDJ = (one/RT)*(two)/
     +     ((one/(phi_tau - one) + one
     +     + two*chi*phi_tau) + half*(trC_tau - three 
     +     - two*dlog(detF))*props(1)*(-0.867)*Vmol/(((one
     +     -phi_tau)**2)*RT)+ props(2)*(-0.867)*(Vmol**2)*
     +     dlog(detF*phi_tau)*((dlog(detF*phi_tau)-2)*(
     +     -(phi_tau**2)+phi_tau+one)*half
     +     +one/phi_tau)/(((1-phi_tau)**2)*RT) 
     +     +Kbulk*Vmol*(dlog(detF*phi_tau)-one)/(RT*phi_tau))

      ! Compute the perturbation on the chemical potential
      !
      if(dabs(mu_tau).gt.one) then
         deltaMU = dabs(mu_tau)*1.d-8
      else
         deltaMU = 1.d-8
      endif


      ! Compute a perturbed polymer volume fraction
      !
      args(1)  = mu_tau + deltaMU
      args(2)  = mu0
      args(3)  = Rgas
      args(4)  = theta
      args(5)  = chi
      args(6)  = Vmol
      args(7)  = Kbulk
      args(8)  = detF
      call solvePhi1(phi_per,args,nargs,phi_t,F_tau)


      ! Compute a perturbed polymer volume fraction
      !
      args(1)  = mu_tau - deltaMU
      args(2)  = mu0
      args(3)  = Rgas
      args(4)  = theta
      args(5)  = chi
      args(6)  = Vmol
      args(7)  = Kbulk
      args(8)  = detF
      call solvePhi1(phi_m,args,nargs,phi_t,F_tau)


      ! Compute the derivative of the time rate of change of the 
      !  polymer volume fraction with respect to the chemical potential
      !
      dPdt_per = (phi_per - phi_t)/dtime
      dPdt_m   = (phi_m - phi_t)/dtime
      DphidotDmu = (dPdt_per - dPdt_m)/(two*deltaMU)


      ! Compute the fluid permeability at this integ. point
      !
      Mfluid = (D*(one/phi_tau - one))/(detF*Vmol*Rgas*theta)


      ! Compute the tangents of the fluid mobility
      !
      DmDmu = DmDphi*DphiDmu
      DmDJ  = DmDphi*DphiDJ



      ! Compute the Cauchy stress
      !
      T_tau = (Gshear*(B_tau-Iden) + detFs*Kbulk*dlog(detFe)*Iden)/detF

      
      ! Compute the 1st Piola stress
      !
      TR_tau = Gshear*(F_tau - FinvT) + detFs*Kbulk*dlog(detFe)*FinvT


      ! Compute dTRdF, the so-called material tangent modulus
      !
      dTRdF = zero
      do i=1,3
         do j = 1,3
            do k = 1,3
               do l = 1,3          
                  dTRdF(i,j,k,l) = dTRdF(i,j,k,l)
     +                 + Gshear*Iden(i,k)*Iden(j,l)
     +                 + Gshear*Finv(l,i)*Finv(j,k)
     +                 + detFs*Kbulk*Finv(j,i)*Finv(l,k)
     +                 - detFs*Kbulk*dlog(detFe)*Finv(l,i)*Finv(j,k)
               enddo
            enddo
         enddo
      enddo
      !
      ! Calculate the so-called spatial tangent modulus, based
      !  on the push forward of the material tangent modulus
      !
      SpTanMod = zero
      do i=1,3
         do j=1,3
            do k=1,3
               do l=1,3
                  do m=1,3
                     do n=1,3
                        SpTanMod(i,j,k,l) = SpTanMod(i,j,k,l) + 
     +                       (dTRdF(i,m,k,n)*F_tau(j,m)*F_tau(l,n))/detF
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo


      ! Compute the displacement - chemical potential modulus
      !
      SpUCMod = (Kbulk/(detFe*phi_tau))*Iden*DphiDmu


      ! Compute the chemical potential - displacement modulus
      !
      SpCUModFac = Mfluid*Iden

      return
      end subroutine integ1

************************************************************************

      subroutine solvePhi1(root,args,nargs,rootOld,Fmat)

      ! This subroutine will numerically solve for the polymer
      !  volume fraction based on the current osmotic pressure
      !  and the previous state.  See numerical recipies RTSAFE.

      implicit none

      integer maxit,j,nargs

      real*8 xacc,f,df,fl,fh,xl,xh,x1,x2,swap,root,dxold,one
      real*8 dx,args(nargs),zero,rootOld,temp,rootMax,rootMin
      real*8 Fmat(3,3)

      parameter(maxit=50)
      parameter(xacc=1.d-6,zero=0.d0,one=1.d0)

      rootMax = 0.9999d0 ! corresponds to nearly 100% dry polymer
      rootMin = 0.05d0   ! corresponds to nearly 100% fluid

      x1 = rootMin
      x2 = rootMax
      call phiFunc1(x1,FL,DF,args,nargs,Fmat)
      call phiFunc1(x2,FH,DF,args,nargs,Fmat)

      if(fl*fh.ge.zero) then
         root = rootOld
         write(*,*) 'FYI, root not bracketed on phi'
         write(*,*) 'fl=',fl
         write(*,*) 'fh=',fh
         write(*,*) 'rootOld=',rootOld
         write(80,*) 'FYI, the root is not bracketed on phi'
         write(80,*) 'fl=',fl
         write(80,*) 'fh=',fh
         write(80,*) 'rootOld=',rootOld

         write(*,*) 'mu =',args(1)
         write(*,*) 'mu0=',args(2)
         write(*,*) 'Rgas=',args(3)
         write(*,*) 'theta=',args(4)
         write(*,*) 'chi=',args(5)
         write(*,*) 'Vmol=',args(6)
         write(*,*) 'Kbulk=',args(7)
         write(*,*) 'detF=',args(8)

         call xit
         return
      endif

C
C		ORIENT THE SEARCH SO THAT F(XL) < 0.
C
      IF( FL .LT. 0.D0 ) THEN
         XL = X1
         XH = X2
      ELSE
         XH = X1
         XL = X2
         SWAP = FL
         FL = FH
         FH = SWAP
      END IF
C
C		INITIALIZE THE GUESS FOR THE ROOT, THE ''STEP SIZE
C		BEFORE LAST'', AND THE LAST STEP
C
      if(rootOld.lt.rootMin) rootOld = rootMin
      if(rootOld.gt.rootMax) rootOld = rootMax
      ROOT = rootOld !0.5D0 *( X1 + X2)
      DXOLD = DABS(X2 - X1)
      DX = DXOLD
      
      call phiFunc1(root,F,DF,args,nargs,Fmat)

C
C			LOOP OVER ALLOWED ITERATIONS
C
      DO 10 J = 1,MAXIT
C
C			BISECT IF NEWTON OUT OF RANGE, OR NOT DECREASING
C			FAST ENOUGH.
C
         IF( ((ROOT-XH)*DF - F)*((ROOT - XL)*DF -F) .GE. 0.D0
     +        .OR. DABS(2.D0*F) .GT. DABS(DXOLD*DF) ) THEN

            DXOLD = DX
            DX = 0.5D0*(XH-XL)
            ROOT = XL + DX
            IF( XL .EQ. ROOT ) THEN
C
C			CHANGE IN ROOT IS NEGLIGIBLE
C
               RETURN
            END IF

         ELSE
C
C			NEWTON STEP IS ACCEPTABLE. TAKE IT.
C
            DXOLD = DX
            DX = F/DF
            TEMP = ROOT
            ROOT = ROOT - DX
            IF( TEMP .EQ. ROOT) THEN
C
C			 CHANGE IN ROOT IS NEGLIGIBLE
C
               RETURN
            END IF

         END IF
C
C		CONVERVEGENCE CRITERION
C
         IF( DABS(DX) .LT. XACC) RETURN

C
C			THE ONE NEW FUNCTION EVALUATION PER ITERATION
C
         call phiFunc1(root,F,DF,args,nargs,Fmat)

C
C		MAINTAIN THE BRACKET ON THE ROOT
C
         IF( F .LT. 0.D0) THEN
            XL = ROOT
            FL = F
         ELSE
            XH = ROOT
            FH = F
         END IF

 10   CONTINUE

      WRITE(*,'(/1X,A)') 'solvePhi EXCEEDING MAXIMUM ITERATIONS'
      WRITE(80,'(/1X,A)') 'solvePhi EXCEEDING MAXIMUM ITERATIONS'

      return
      end subroutine solvePhi1

****************************************************************************

      subroutine phiFunc1(phi,fr,dfr,args,nargs,Fmat)

      ! This subroutine serves as the function we would like to solve for
      !  the polymer volume fraction by finding phi such that ``f=0''

      use global

      implicit none

      integer nargs,NeoHookean,Langevin,material
      parameter(NeoHookean=1,Langevin=2)

      real*8 args(nargs),fr,dfr,mu,mu0,Rgas,theta,chi,Vmol,Gshear,Kbulk
      real*8 detF,phi,RT,C(3,3),trC,Fmat(3,3),G0,K0,pnewdt,time(2)

      real*8 zero,one,two,three,third,half
      parameter(zero=0.d0,one=1.d0,two=2.d0,three=3.d0,third=1.d0/3.d0)
      parameter(half=1.d0/2.d0)

      
      ! Obtain relevant quantities
      !
      mu     = args(1)
      mu0    = args(2)
      Rgas   = args(3)
      theta  = args(4)
      chi    = args(5)
      Vmol   = args(6)
      Kbulk  = args(7)
      detF   = args(8)
      Gshear = args(9)
      G0     = args(10)
      K0     = args(11)


      ! Compute the useful quantity
      !
      RT = Rgas*theta

      C = matmul(transpose(Fmat),Fmat)
      trC = C(1,1) + C(2,2) + C(3,3)


      ! Compute the residual
      !
      fr = (mu0 - mu)/RT
     +     + (dlog(one - phi) + phi + chi*phi*phi)
     +     + half*(trC-three-two*dlog(detF))*G0*
     +     (-0.867)*Vmol*phi/(RT*(1-phi))+half*Kbulk
     +     *Vmol*((dlog(detF*phi))**two)/RT+half*((dlog(
     +     detF*phi))**two)*K0*(-0.867)*(Vmol**two)/(
     +     (one-phi)*RT)-Vmol*Kbulk*dlog(detF*phi)/RT



      ! Compute the tangent
      !
      if(phi.gt.0.999d0) then
         dfr = zero
      else
         dfr = (one/(phi - one) + one
     +        +two*chi*phi) + half*(trC - three 
     +        -two*dlog(detF))*G0*(-0.867)*Vmol/(((one
     +        -phi)**2)*RT) + K0*(-0.867)*(Vmol**2)*
     +        dlog(detF*phi)*((dlog(detF*phi)-2)*(
     +        -(phi**2)+phi+one)*half
     +        +one/phi)/(((1-phi)**2)*RT) 
     +        +Kbulk*Vmol*(dlog(detF*phi)-one)/(RT*phi)
      endif

      return
      end subroutine phiFunc1


************************************************************************
************************************************************************
************************************************************************
************************************************************************

      subroutine AssembleElement(nDim,nNode,ndofel,
     +     Ru,Rc,Kuu,Kuc,Kcu,Kcc,
     +     rhs,amatrx)
      
      !
      ! Subroutine to assemble the local elements residual and tangent
      !

      implicit none

      integer i,j,k,l,m,n,A11,A12,B11,B12,nDim,nNode,nDofEl,nDofN

      real*8 Ru(nDim*nNode,1),Rc(nNode,1),Kuu(nDim*nNode,nDim*nNode)
      real*8 Kcc(nNode,nNode),Kuc(nDim*nNode,nNode),rhs(ndofel,1)
      real*8 Kcu(nNode,nDim*nNode),amatrx(ndofel,ndofel)


      ! Total number of degrees of freedom per node
      !
      nDofN = nDofEl/nNode


      ! init
      !
      rhs(:,1) = 0.d0
      amatrx = 0.d0

      if(nDim.eq.2) then
         !
         ! Assemble the element level residual
         !
         do i=1,nNode
            A11 = nDofN*(i-1)+1
            A12 = nDim*(i-1)+1
            !
            ! displacement
            !
            rhs(A11,1) = Ru(A12,1)
            rhs(A11+1,1) = Ru(A12+1,1)
            !
            ! chemical potential
            !
            rhs(A11+2,1) = Rc(i,1)
         enddo
         !
         ! Assemble the element level tangent matrix
         !
         do i=1,nNode
            do j=1,nNode
               A11 = nDofN*(i-1)+1
               A12 = nDim*(i-1)+1
               B11 = nDofN*(j-1)+1
               B12 = nDim*(j-1)+1
               !
               ! displacement
               !
               amatrx(A11,B11) = Kuu(A12,B12)
               amatrx(A11,B11+1) = Kuu(A12,B12+1)
               amatrx(A11+1,B11) = Kuu(A12+1,B12)
               amatrx(A11+1,B11+1) = Kuu(A12+1,B12+1)
               !
               ! chemical potential
               !
               amatrx(A11+2,B11+2) = Kcc(i,j)
               !
               ! displacement - chemical potential
               !
               amatrx(A11,B11+2) = Kuc(A12,j)
               amatrx(A11+1,B11+2) = Kuc(A12+1,j)
               !
               ! chemical potential - displacement
               !
               amatrx(A11+2,B11) = Kcu(i,B12)
               amatrx(A11+2,B11+1) = Kcu(i,B12+1)
               !
            enddo
         enddo
         !
      else
         write(*,*) 'How did you get nDim=',nDim
         call xit
      endif

      return
      end subroutine AssembleElement

!****************************************************************************
!     Element subroutines
!****************************************************************************

      subroutine xint2D4pt(xi,w,nIntPt)
      !
      ! This subroutine will get the integration point locations
      !  and corresponding gauss quadrature weights for 2D elements
      !  using 4 gauss points for integration
      !
      !  xi(nIntPt,2): xi,eta coordinates for the integration pts
      !  w(nIntPt):    corresponding integration weights
      !
      implicit none
      !
      integer nIntPt,nDim
      !
      real*8 xi(4,2), w(4)


      ! Initialize
      !
      w = 0.d0
      xi = 0.d0


      ! Number of Gauss points
      !
      nIntPt = 4


      ! Gauss weights
      !
      w(1) = 1.d0
      w(2) = 1.d0
      w(3) = 1.d0
      w(4) = 1.d0
      

      ! Gauss pt locations in master element
      !
      xi(1,1) = -dsqrt(1.d0/3.d0)
      xi(1,2) = -dsqrt(1.d0/3.d0)
      xi(2,1) = dsqrt(1.d0/3.d0)
      xi(2,2) = -dsqrt(1.d0/3.d0)
      xi(3,1) = -dsqrt(1.d0/3.d0)
      xi(3,2) = dsqrt(1.d0/3.d0)
      xi(4,1) = dsqrt(1.d0/3.d0)
      xi(4,2) = dsqrt(1.d0/3.d0)


      return
      end subroutine xint2D4pt

************************************************************************

      subroutine xintSurf2D1pt(face,xLocal,yLocal,w)

      ! This subroutine will get the integration point locations
      !  and corresponding gauss quadrature weights for 2D elements
      !  using 2 gauss points for surface integration
      !
      !  xLocal(nIntPt): x coordinates for the integration pts
      !  yLocal(nIntPt): y coordinates for the integration pts
      !  w(nIntPt):    corresponding integration weights

      implicit none

      integer face

      real*8 xLocal(1),yLocal(1),w(1),zero,one,two
      parameter(zero=0.d0,one=1.d0,two=2.d0)


      ! Gauss weights
      !
      w(1) = two
      

      ! Gauss pt locations in master element
      !
      if(face.eq.1) then
         xLocal(1) = zero
         yLocal(1) = -one
      elseif(face.eq.2) then
         xLocal(1) = one
         yLocal(1) = zero
      elseif(face.eq.3) then
         xLocal(1) = zero
         yLocal(1) = one
      elseif(face.eq.4) then
         xLocal(1) = -one
         yLocal(1) = zero
      else
         write(*,*) 'face.ne.1,2,3,4'
         write(80,*) 'face.ne.1,2,3,4'
         call xit
      endif

      end subroutine xintSurf2D1pt

************************************************************************

      subroutine calcShape2DLinear(nIntPt,xi_int,intpt,sh,dshxi)
      !
      ! Calculate the shape functions and their derivatives at the
      ! given integration point in the master element


      ! Calculate the shape functions and their derivatives at the
      ! given integration point in the master element
      !
      !                          eta
      !   4-----------3          |
      !   |           |          |
      !   |           |          |
      !   |           |          |
      !   |           |          |
      !   |           |          O--------- xi
      !   1-----------2        origin at center
      !
      !
      ! sh(i) = shape function of node i at the intpt.
      ! dshxi(i,j) = derivative wrt j direction of shape fn of node i
      !
      implicit none
      !
      integer intpt,nDim,nIntPt
      !
      real*8 xi_int(nIntPt,2),sh(4),dshxi(4,2),xi,eta
      !
      real*8 zero,one,fourth
      parameter(zero=0.d0,one=1.d0,fourth=1.d0/4.d0)
      

      ! Location in the master element
      !
      xi = xi_int(intpt,1)
      eta = xi_int(intpt,2)
      
      
      ! The shape functions
      !
      sh(1) = fourth*(one - xi)*(one - eta)
      sh(2) = fourth*(one + xi)*(one - eta)
      sh(3) = fourth*(one + xi)*(one + eta)
      sh(4) = fourth*(one - xi)*(one + eta)
      
      
      ! The first derivatives
      !
      dshxi(1,1) = -fourth*(one - eta)
      dshxi(1,2) = -fourth*(one - xi)
      dshxi(2,1) = fourth*(one - eta)
      dshxi(2,2) = -fourth*(one + xi)
      dshxi(3,1) = fourth*(one + eta)
      dshxi(3,2) = fourth*(one + xi)
      dshxi(4,1) = -fourth*(one + eta)
      dshxi(4,2) = fourth*(one - xi)
      

      return
      end subroutine calcShape2DLinear

!************************************************************************


      subroutine computeSurf(xLocal,yLocal,face,coords,sh,ds)

      ! This subroutine computes the shape functions, derivatives
      !  of shape functions, and the length ds, so that one can
      !  do the numerical integration on the boundary for fluxes 
      !  on the 4-node quadrilateral elements

      implicit none

      integer face

      real*8 xLocal,yLocal,ds,dshxi(4,2),sh(4),dXdXi,dXdEta,dYdXi
      real*8 dYdEta,one,coords(2,4),fourth,shape,normal(2,1)
      parameter(one=1.d0,fourth=1.d0/4.d0)

      sh(1) = fourth*(one - xLocal)*(one - yLocal)
      sh(2) = fourth*(one + xLocal)*(one - yLocal)
      sh(3) = fourth*(one + xLocal)*(one + yLocal)
      sh(4) = fourth*(one - xLocal)*(one + yLocal)
      
      dshxi(1,1) = -fourth*(one - yLocal)
      dshxi(1,2) = -fourth*(one - xLocal)
      dshxi(2,1) = fourth*(one - yLocal)
      dshxi(2,2) = -fourth*(one + xLocal)
      dshxi(3,1) = fourth*(one + yLocal)
      dshxi(3,2) = fourth*(one + xLocal)
      dshxi(4,1) = -fourth*(one + yLocal)
      dshxi(4,2) = fourth*(one - xLocal)

      dXdXi = dshxi(1,1)*coords(1,1)+dshxi(2,1)*coords(1,2)
     +     + dshxi(3,1)*coords(1,3)+dshxi(4,1)*coords(1,4)
      dXdEta = dshxi(1,2)*coords(1,1)+dshxi(2,2)*coords(1,2)
     +     + dshxi(3,2)*coords(1,3)+dshxi(4,2)*coords(1,4)
      dYdXi = dshxi(1,1)*coords(2,1)+dshxi(2,1)*coords(2,2)
     +     + dshxi(3,1)*coords(2,3)+dshxi(4,1)*coords(2,4)
      dYdEta = dshxi(1,2)*coords(2,1)+dshxi(2,2)*coords(2,2)
     +     + dshxi(3,2)*coords(2,3)+dshxi(4,2)*coords(2,4)


      ! Jacobian of the mapping
      !
      if((face.eq.2).or.(face.eq.4)) then
         ds = dsqrt(dXdEta*dXdEta + dYdEta*dYdEta)
      elseif((face.eq.1).or.(face.eq.3)) then
         ds = dsqrt(dXdXi*dXdXi + dYdXi*dYdXi)
      else
         write(*,*) 'never should get here'
         call xit
      endif


      ! Surface normal, outward pointing in this case. Useful for
      !  ``follower'' type loads. The normal is referential or spatial
      !  depending on which coords were supplied to this subroutine
      !  (NOT fully tested)
      !
      if((face.eq.2).or.(face.eq.4)) then
         normal(1,1) = dYdEta/dsqrt(dXdEta*dXdEta + dYdEta*dYdEta)
         normal(2,1) = -dXdEta/dsqrt(dXdEta*dXdEta + dYdEta*dYdEta)
         if(face.eq.4) normal = -normal
      elseif((face.eq.1).or.(face.eq.3)) then
         normal(1,1) = dYdXi/dsqrt(dXdXi*dXdXi + dYdXi*dYdXi)
         normal(2,1) = -dXdXi/dsqrt(dXdXi*dXdXi + dYdXi*dYdXi)
         if(face.eq.3) normal = -normal
      else
         write(*,*) 'never should get here'
         call xit
      endif

      return
      end subroutine computeSurf

************************************************************************

      subroutine mapShape2D(nNode,dshxi,coords,dsh,detMapJ,stat)
      !
      ! Map derivatives of shape fns from xi-eta-zeta domain
      !  to x-y-z domain.
      !
      implicit none
      !
      integer i,j,k,nNode,ieror,stat
      !
      real*8 dshxi(nNode,2),dsh(nNode,2),coords(3,nNode),mapJ(2,2),
     +  mapJ_inv(2,2),detmapJ
      !
      real*8 zero,one,two,half,fourth,eighth
      parameter(zero=0.d0,one=1.d0,two=2.d0,half=0.5d0,fourth=0.25d0,
     +     eighth=1.d0/8.d0)


      ! Calculate the mapping Jacobian matrix:
      !
      mapJ = zero
      do i=1,2
        do j=1,2
          do k=1,nNode
              mapJ(i,j) = mapJ(i,j) + dshxi(k,i)*coords(j,k)
          end do
        end do
      end do


      ! Calculate the inverse and the determinant of Jacobian
      !
      call matInv2D(mapJ,mapJ_inv,detMapJ,stat)
      if(stat.eq.0) then
         write(*,*) 'Problem: detF.lt.zero in mapShape2D'
         call xit
      endif


      ! Calculate first derivatives wrt x, y, z
      !
      dsh = transpose(matmul(mapJ_inv,transpose(dshxi)))
      

      return
      end subroutine mapShape2D

!*************************************************************************

      subroutine mapShape2Da(nNode,dshxi,coords,dsh,detMapJ,stat)
      !
      ! Map derivatives of shape fns from xi-eta-zeta domain
      !  to x-y-z domain.
      !
      ! This subroutine is exactly the same as the regular mapShape2D
      !  with the exception that coords(2,nNode) here and coords(3,nNode)
      !  in the regular.  I have noticed that a "heat transfer" and 
      !  "static" step uses MCRD=2, but for "coupled-temperature-displacement"
      !  you will get MCRD=3, even for a plane analysis.
      !
      implicit none
      !
      integer i,j,k,nNode,ieror,stat
      !
      real*8 dshxi(nNode,2),dsh(nNode,2),coords(2,nNode),mapJ(2,2),
     +  mapJ_inv(2,2),detmapJ
      !
      real*8 zero,one,two,half,fourth,eighth
      parameter(zero=0.d0,one=1.d0,two=2.d0,half=0.5d0,fourth=0.25d0,
     +     eighth=1.d0/8.d0)


      ! Calculate the mapping Jacobian matrix:
      !
      mapJ = zero
      do i=1,2
        do j=1,2
          do k=1,nNode
              mapJ(i,j) = mapJ(i,j) + dshxi(k,i)*coords(j,k)
          end do
        end do
      end do


      ! Calculate the inverse and the determinant of Jacobian
      !
      call matInv2D(mapJ,mapJ_inv,detMapJ,stat)
      if(stat.eq.0) then
         write(*,*) 'Problem: detF.lt.zero in mapShape2Da'
         call xit
      endif


      ! Calculate first derivatives wrt x, y, z
      !
      dsh = transpose(matmul(mapJ_inv,transpose(dshxi)))
      

      return
      end subroutine mapShape2Da

!****************************************************************************
!     Utility subroutines
!****************************************************************************

      subroutine matInv3D(A,A_inv,det_A,istat)
      !
      ! Returns A_inv, the inverse and det_A, the determinant
      ! Note that the det is of the original matrix, not the
      ! inverse
      !
      implicit none
      !
      integer istat
      !
      real*8 A(3,3),A_inv(3,3),det_A,det_A_inv


      istat = 1
      
      det_A = A(1,1)*(A(2,2)*A(3,3) - A(3,2)*A(2,3)) -
     +        A(2,1)*(A(1,2)*A(3,3) - A(3,2)*A(1,3)) +
     +        A(3,1)*(A(1,2)*A(2,3) - A(2,2)*A(1,3))
      
      if (det_A .le. 0.d0) then
        write(*,*) 'WARNING: subroutine matInv3D:'
        write(*,*) 'WARNING: det of mat=',det_A
        istat = 0
        return
      end if
          
      det_A_inv = 1.d0/det_A
        
      A_inv(1,1) = det_A_inv*(A(2,2)*A(3,3)-A(3,2)*A(2,3))
      A_inv(1,2) = det_A_inv*(A(3,2)*A(1,3)-A(1,2)*A(3,3))
      A_inv(1,3) = det_A_inv*(A(1,2)*A(2,3)-A(2,2)*A(1,3))
      A_inv(2,1) = det_A_inv*(A(3,1)*A(2,3)-A(2,1)*A(3,3))
      A_inv(2,2) = det_A_inv*(A(1,1)*A(3,3)-A(3,1)*A(1,3))
      A_inv(2,3) = det_A_inv*(A(2,1)*A(1,3)-A(1,1)*A(2,3))
      A_inv(3,1) = det_A_inv*(A(2,1)*A(3,2)-A(3,1)*A(2,2))
      A_inv(3,2) = det_A_inv*(A(3,1)*A(1,2)-A(1,1)*A(3,2))
      A_inv(3,3) = det_A_inv*(A(1,1)*A(2,2)-A(2,1)*A(1,2))
      

      return
      end subroutine matInv3D

!****************************************************************************

      subroutine matInv2D(A,A_inv,det_A,istat)
      !
      ! Returns A_inv, the inverse, and det_A, the determinant
      ! Note that the det is of the original matrix, not the
      ! inverse
      !
      implicit none
      !
      integer istat
      !
      real*8 A(2,2),A_inv(2,2),det_A,det_A_inv

      
      istat = 1
      
      det_A = A(1,1)*A(2,2) - A(1,2)*A(2,1)
        
      if (det_A .le. 0.d0) then
        write(*,*) 'WARNING: subroutine matInv2D:'
        write(*,*) 'WARNING: det of mat=',det_A
        istat = 0
        return
      end if
            
      det_A_inv = 1.d0/det_A
          
      A_inv(1,1) =  det_A_inv*A(2,2)
      A_inv(1,2) = -det_A_inv*A(1,2)
      A_inv(2,1) = -det_A_inv*A(2,1)
      A_inv(2,2) =  det_A_inv*A(1,1)


      return
      end subroutine matInv2D

!****************************************************************************

      subroutine mdet(A,det)
      !
      ! This subroutine calculates the determinant
      ! of a 3 by 3 matrix [A]
      !
      implicit none
      !
      real*8  A(3,3),det


      det = A(1,1)*A(2,2)*A(3,3) 
     +	  + A(1,2)*A(2,3)*A(3,1)
     +	  + A(1,3)*A(2,1)*A(3,2)
     +	  - A(3,1)*A(2,2)*A(1,3)
     +	  - A(3,2)*A(2,3)*A(1,1)
     +	  - A(3,3)*A(2,1)*A(1,2)


      return
      end subroutine mdet
	
!****************************************************************************

      subroutine onem(A)
      !
      ! This subroutine stores the identity matrix in the
      ! 3 by 3 matrix [A]
      !
      implicit none
      !
      integer i,j
      !
      real*8 A(3,3)


      do i=1,3
         do J=1,3
	    if (i .eq. j) then
              A(i,j) = 1.0
            else
              A(i,j) = 0.0
            end if
         end do
      end do


      return
      end subroutine onem

****************************************************************************
