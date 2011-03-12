!  F77Test.f90 
!
!  FUNCTIONS:
!  PoseEst      - Entry point of console application.
!
!

!****************************************************************************
!
!  PROGRAM: PoseEst
!
!  PURPOSE:  Test the SoftPOSIT pose estimation algorithm
!
!****************************************************************************


    program PoseEst

    implicit none
    
!    interface           ! for external procedures
!      subroutine SoftPOSIT(rot, trans, imagePts, worldPts, beta0, noiseStd, initRot, initTrans, focalLength, center)
!        use mkl95_blas, only: syr, axpy, symv
!        use mkl95_lapack, only: sytrf, sytri, gesvd
!        integer, intent(in) :: imagePts(:,:), center(2)
!        real, intent(in) :: worldPts(:,:), beta0, noiseStd, initRot(3,3), initTrans(3), focalLength
!        real, intent(out) :: rot(3,3), trans(3)
!      end subroutine SoftPOSIT
!    end interface
    
    integer :: imagePtsData(4,2), centerData(2)
    real :: worldPts(4,3), beta0, noiseStd, initRot(3,3), initTrans(3), focalLength, rot(3,3), trans(3)
    
    !integer*8 :: clock1, clock2
    real :: clock1, clock2
    !real :: freq
    
    integer :: count = 0
    
    imagePtsData = reshape( (/ 299, 641, 330, 522, 249, 254, 55, 436 /), (/ 4, 2 /) )
    worldPts = reshape( (/ -3.75, 7.50, -3., 3., 0., 0., -5.0, 5.0, 0.5, 2.75, -2.0, -2.0 /), (/ 4, 3 /) )
    beta0 = 2.0E-4
    noiseStd = 0.0
    initRot = reshape( (/1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0/), (/ 3, 3 /) )
    initTrans = (/0.0, 0.0, 30.0/)
    focalLength = 982.1
    centerData = (/376, 240/)
    
    !freq = getcpufrequency()
    !freq = 3.0D9
    !call getcpuclocks(clock1)
    call cpu_time(clock1)
    do count = 1, 1000
    call SoftPOSIT(rot, trans, imagePtsData, worldPts, beta0, noiseStd, initRot, initTrans, focalLength, centerData)
    !call getcpuclocks(clock2)
    end do
    call cpu_time(clock2)
    
    print *, 'ImagePts = ', 4
    print *, 'WorldPts = ', 4
    print *, 'Rot = '
    print *, rot
    print *, 'Trans = '
    print *, trans
    print *, 'Eplased time = ', (clock2-clock1)







    contains




!!! Subroutine of SoftPOSIT algorithm

    
    subroutine SoftPOSIT(rot, trans, imagePts, worldPts, beta0, noiseStd, initRot, initTrans, focalLength, center)
    use mkl95_blas, only: syr, axpy, symv, ger
    use mkl95_lapack, only: sytrf, sytri, gesvd
    implicit none
    integer, intent(in) :: imagePts(:,:), center(2)
    real, intent(in) :: worldPts(:,:), beta0, noiseStd, initRot(3,3), initTrans(3), focalLength
    real, intent(out) :: rot(3,3), trans(3)
    
!   imagePts: M * 2 matrix
!   worldPts: N * 3 matrix
    integer, parameter :: nbImagePts = 4
    integer, parameter :: nbWorldPts = 4
    
    integer :: wait
    

    
    real :: alpha
    real :: maxDelta
    real, parameter :: betaFinal = 0.5                 ! Terminate iteration when beta == betaFinal.
    real, parameter :: betaUpdate = 1.05               ! Update rate on beta.
    real, parameter :: epsilon0 = 0.01                 ! Used to initialize assignement matrix.

    integer, parameter :: maxCount = 2
    integer, parameter :: minBetaCount = 20
    
    !integer :: nbImagePts
    !integer :: nbWorldPts

    integer :: minNbPts
    integer :: maxNbPts
    
    real :: scale
    
    real, save :: centeredImage(nbImagePts,2)
    
    real, save :: imageOnes(nbImagePts)
!    real :: worldOnes(nbWorldPts, 1) = 1.0  !not used

!   homogeneousWorldPts = [worldPts, worldOnes];    
    real, save :: homogeneousWorldPts(nbWorldPts, 4)
    

    
!    wk = homogeneousWorldPts * [rot(3,:)/trans(3), 1]';
    real, save :: wk(nbWorldPts)
    
    real :: r1T(4),  r2T(4)
    real :: r1Tprev(4),  r2Tprev(4)
    
    integer :: betaCount = 0;
    logical :: poseConverged = .FALSE.
    logical :: assignConverged = .FALSE.
    logical :: foundPose = .FALSE.
    real :: beta
    
    real, save :: assignMat(nbImagePts+1,nbWorldPts+1)
    
    real, save :: projectedU(nbWorldPts)
    real, save :: projectedV(nbWorldPts)
    real, save :: replicatedProjectedU(nbImagePts, nbWorldPts), replicatedProjectedV(nbImagePts, nbWorldPts)
    
    real, save  :: wkxj(nbImagePts,nbWorldPts), wkyj(nbImagePts,nbWorldPts)
    real, save :: distMat(nbImagePts, nbWorldPts)
    
!    integer :: numMatchPts
    real :: sumNonslack
    real, save :: summedByColAssign(nbWorldPts)
    
    real :: sumSkSkT(4,4) = 0.0
    integer :: j, k
    integer :: inv_fail
    integer :: ipiv(4)
    
    integer :: count
    
    real :: weightedUi(4) = 0.0
    real :: weightedVi(4) = 0.0
    real :: weightedTemp(4) = 0.0
    
    real :: svdMat(3,2), svdS(2), svdU(3,3), svdVt(2,2), A(3,2)
    real :: delta
    


!!!!!!!!!!  End of initialization  !!!!!!!!!!!!!!!
    
    !print *, 'start SoftPOSIT...'
    
    alpha = 9.21*noiseStd**2 + 1
    maxDelta = sqrt(alpha)/2.0        ! Max allowed error per world point.
    
    !nbImagePts = size(imagePts, 1)    ! Number of image points (M below).
    !nbWorldPts = size(worldPts, 1)    ! Number of world points (N below).
    
    minNbPts = min(nbImagePts, nbWorldPts)
    maxNbPts = max(nbImagePts, nbWorldPts)
    
    scale = 1.0/(maxNbPts + 1)

    
    centeredImage(:,1) = (imagePts(:,1) - center(1))/focalLength
    centeredImage(:,2) = (imagePts(:,2) - center(2))/focalLength
    
   
    imageOnes = 1.0
    
    
    homogeneousWorldpts(:,1:3) = worldPts
    homogeneousWorldpts(:,4) = 1.0
    
    rot = initRot;
    trans = initTrans;
    
    
    wk = matmul(homogeneousWorldPts, (/rot(3,:)/trans(3), 1.0/))
!    wk = matmul(homogeneousWorldPts, (/1.0, 1.0, 1.0, 1.0/))
!   Need to use (/a, b, c, d/) for const arrays
    
    r1T = (/rot(1,:)/trans(3), trans(1)/trans(3)/);
    r2T = (/rot(2,:)/trans(3), trans(2)/trans(3)/);
    
    betaCount = 0;
    poseConverged = .FALSE.
    assignConverged = .FALSE.
    foundPose = .FALSE.
    beta = beta0

!   The assignment matrix includes a slack row and slack column.
    
    
    assignMat = 1.0 + epsilon0;

    
!   Do some initialization outside the loop
    

   
    
    do while (beta < betaFinal .AND. (.NOT.assignConverged))

      !print *,'Start while loop...', assignMat
      projectedU = matmul(homogeneousWorldPts, r1T)          
      projectedV = matmul(homogeneousWorldPts, r2T)
      
      
      
      replicatedProjectedU = spread(projectedU, 1, nbImagePts)
      replicatedProjectedV = spread(projectedV, 1, nbImagePts)
      
      
      
      !wkxj = dot_product(centeredImage(:,1), wk)
      !wkyj = dot_product(centeredImage(:,2), wk)
      wkxj = 0
      wkyj = 0
      
      call ger(wkxj, centeredImage(:,1), wk)
      call ger(wkyj, centeredImage(:,2), wk)
      
      
      
      
      distMat = focalLength**2 * ((replicatedProjectedU - wkxj)**2 + (replicatedProjectedV - wkyj)**2)
      
      
!    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    % Use the softAssign algorithm to compute the best assignment of image to
!    % world points based on the computed distances.  The use of alpha 
!    % determines when to favor assignments instead of use of slack.
!    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      assignMat(1:nbImagePts, 1:nbWorldPts) = scale*exp(-beta*(distMat - alpha))
      assignMat(1:nbImagePts+1, nbWorldPts+1) = scale
      assignMat(nbImagePts+1, 1:nbWorldPts+1) = scale
      
      !print *, distMat

          
      call sinkhornImp(assignMat)    ! My "improved" Sinkhorn.
     
      
      !print *, 'After sinkhorm Imp', assignMat
      !read *, wait
      
!    % About how many matching model points do we currently have?

!      numMatchPts = numMatches(assignMat)

      sumNonslack = sum(assignMat(1:nbImagePts,1:nbWorldPts))
      
!    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    % Use POSIT to calculate the pose that minimizes the objective function
!    % (equation (10), the weighted sum of the differences between corrected
!    % image points and SOP's of the corresponding world points), given the
!    % current assignment.  The pose parameters and the w[i] are iteratively
!    % updated until some convergence criterion is satisfied.
!    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          
      summedByColAssign = sum(assignMat(1:nbImagePts, 1:nbWorldPts), 1)
               
      sumSkSkT = 0.0
      
      do k = 1, nbWorldPts
        call syr(sumSkSkT, homogeneousWorldPts(k,:), 'U', summedByColAssign(k))
        !sumSkSkT = sumSkSkT + summedByColAssign(k) * homogeneousWorldPts(k,:)' * homogeneousWorldPts(k,:);
      end do

     
      call sytrf(sumSkSkT, 'U', ipiv)
      call sytri(sumSkSkT, ipiv, 'U', inv_fail)  ! inv_fail = 0: succeed
      
     
      poseConverged = .FALSE.                              ! Initialize for POSIT loop.
      count = 0
      
      r1Tprev = r1T
      r2Tprev = r2T
      
      
      weightedUi = 0.0
      weightedVi = 0.0
      
      do j = 1, nbImagePts
            do k = 1, nbWorldPts
                weightedTemp = assignMat(j,k) * wk(k) * homogeneousWorldPts(k,:)
                call axpy(weightedTemp, weightedUi, centeredImage(j,1))
                call axpy(weightedTemp, weightedVi, centeredImage(j,2))
                
                ! weightedUi = weightedUi + assignMat(j,k) * wk(k) * centeredImage(j,1) * homogeneousWorldPts(k,:)';
                ! weightedVi = weightedVi + assignMat(j,k) * wk(k) * centeredImage(j,2) * homogeneousWorldPts(k,:)';
                
            end do
      end do
      
        
      call symv(sumSkSkT, weightedUi, r1T, 'U', 1.0, 0.0)
      call symv(sumSkSkT, weightedVi, r2T, 'U', 1.0, 0.0)
       
 
      !r1T = objectMat * weightedUi;               % M
      !r2T = objectMat * weightedVi;              % N
        
      ! Always use this        if 1  % Chang & Tsai calculation of R and T.
      
      svdMat(:,1) = r1T(1:3)
      svdMat(:,2) = r2T(1:3)
        
      call gesvd(svdMat, svdS, svdU, svdVt)
      
      ! [U, S, V] = svd([r1T(1:3)'; r2T(1:3)']');
      A = matmul(matmul(svdU, reshape((/1, 0, 0, 0, 1, 0/), (/3, 2/))), svdVt)
      ! A = U * [1 0; 0 1; 0 0] * V';
      
      rot(:,1) = A(:,1)
      rot(:,2) = A(:,2)
      !r1 = A(:,1);
      !r2 = A(:,2);
      rot(:,3) = cross_prod(A(:,1), A(:,2))
      !r3 = cross(r1,r2);
      trans(3) = 2.0 / (svdS(1) + svdS(2))
      trans(1) = r1T(4) * trans(3)
      trans(2) = r2T(4) * trans(3)
      !Tz = 2 / (S(1,1) + S(2,2));
      !Tx = r1T(4) * Tz;
      !Ty = r2T(4) * Tz;
      
      !r3T= [r3; Tz]; 
      r1T = (/rot(:,1), trans(1)/) / trans(3)
      r2T = (/rot(:,2), trans(2)/) / trans(3)
      !r1T = [r1; Tx]/Tz;
      !r2T = [r2; Ty]/Tz;
    
      wk = matmul(homogeneousWorldPts, (/rot(:,3), trans(3)/)) / trans(3)
      !wk = homogeneousWorldPts * r3T /Tz;
      
      delta = sqrt(sum(assignMat(1:nbImagePts, 1:nbWorldPts) * distMat)/nbWorldPts)
      poseConverged = delta < maxDelta
      count = count + 1                         !% Number of POSIT iterations.
!            % Update the "annealing temperature" and determine if the assignments
!    % have converged.
      beta = betaUpdate * beta
      betaCount = betaCount + 1   !% Number of deterministic annealing iterations.
      assignConverged = poseConverged .AND. betaCount > minBetaCount

!    % Form the best estimates for the translation vector and rotation matrix.
!    trans = [Tx; Ty; Tz];
!    rot = [r1'; r2'; r3'];

!%     if ismember(dispLevel,[5,6])
!%         disp(' '); disp('Current transformation:'); rot, trans
!%     end

!    % Has the pose converged?
      foundPose = delta < maxDelta .AND. betaCount > minBetaCount 
      
    end do
    
    end subroutine SoftPOSIT
    
    !contains
      subroutine sinkhornImp(M)
        implicit none
        real, dimension(:,:), intent(inout) :: M
        
        integer :: wait
        
        integer, parameter :: iMaxIterSinkhorn=60          ! In PAMI paper
        real, parameter :: fEpsilon2 = 0.001            ! Used in termination of Sinkhorn Loop.
        
        integer, parameter :: nbCols = 5, nbRows = 5
        
        integer :: iNumSinkIter = 0
        !integer :: nbRows, nbCols
        real :: fMdiffsum
        
        real, save :: Mprev(nbRows, nbCols)
        real, save :: McolSums(nbCols), MrowSums(nbRows)
        real, save :: McolSumsRep(nbRows, nbCols), MrowSumsRep(nbRows, nbCols)
        
        integer :: i
        integer :: rowPosmax = 0
        integer, save :: posmax(nbCols-1, 2)
        real, save :: ratios(nbCols-1, 2)
        
        !nbRows = size(M,1)
        !nbCols = size(M,2)
        fMdiffSum = fEpsilon2 + 1
                

        
        call maxPosRatio(rowPosmax, posmax, ratios, M)
        

        
        
        
        do while (abs(fMdiffSum) > fEpsilon2 .AND. iNumSinkIter < iMaxIterSinkhorn)
    
          Mprev = M;  ! Save M from previous iteration to test for loop termination

    !% Col normalization (except outlier row - do not normalize col slacks 
    !% against each other)
          McolSums = sum(M, 1)  
          McolSums(nbCols) = 1;  !% Don't normalize slack col terms against each other.
          McolSumsRep = spread(McolSums, 1, nbRows)
          M = M / McolSumsRep


    !% Fix values in the slack column.
          do i = 1, rowPosmax
            M(posmax(i,1),nbCols) = ratios(i,1)*M(posmax(i,1),posmax(i,2))
          end do

    !% Row normalization (except outlier row - do not normalize col slacks 
    !% against each other)
          MrowSums = sum(M, 2)  !% Column vector.
          MrowSums(nbRows) = 1  !% Don't normalize slack row terms against each other.
          MrowSumsRep = spread(MrowSums, 2, nbCols)  
          M = M / MrowSumsRep

    !% Fix values in the slack row.
          do i = 1, rowPosmax
            M(nbRows,posmax(i,2)) = ratios(i,2)*M(posmax(i,1),posmax(i,2))
          end do

          iNumSinkIter = iNumSinkIter + 1
          fMdiffSum = sum(abs(M - Mprev))
          
          !print *, McolSumsRep
          !print *, M
      
        end do
           
      end subroutine sinkhornImp
      
      subroutine maxPosRatio(rowPosmax, posmax, ratios, assignMat)
        implicit none
        real, dimension(:,:), intent(in) :: assignMat
        integer, intent(out) :: rowPosmax
        integer, dimension(:,:), intent(out) :: posmax
        real, dimension(:,:), intent(out) :: ratios
      
        integer :: nrows, ncols, nimgpnts, nmodpnts
        real :: vmax
        integer :: imax
        
        integer :: k
        
        rowPosmax = 0 

        nrows = size(assignMat,1)
        ncols = size(assignMat,2)
        nimgpnts  = nrows - 1
        nmodpnts = ncols - 1
        
        do k = 1, nmodpnts
        
          vmax = maxval(assignMat(:,k))                  !% max value in column k.
          imax = maxloc(assignMat(:,k),1)
        
          if (imax == nrows)  cycle                  !% slack value is maximum in this column.

          if ( all(vmax > assignMat(imax,1:k-1) .AND. vmax > assignMat(imax,k+1:ncols)) ) then
              rowPosmax = rowPosmax + 1
              posmax(rowPosmax,:) = (/imax, k/) != [pos; [imax, k]];     % This value is maximal in its row & column.

    !% Compute the ratios to row and column slack values.
    !rr = assignMat(imax,ncols)/assignMat(imax,k);
    !cr = assignMat(nrows,k)/assignMat(imax,k);
              ratios(rowPosmax,:) = (/assignMat(imax,ncols), assignMat(nrows,k)/) / assignMat(imax,k)  ![rr cr]];
          end if
        end do
        
      end subroutine maxPosRatio
   

      function cross_prod(a, b)
        implicit none
        real, dimension(3) :: cross_prod
        real, dimension(3), intent(in) :: a, b
        
        cross_prod(1) = a(2) * b(3) - a(3) * b(2)
        cross_prod(2) = a(3) * b(1) - a(1) * b(3)
        cross_prod(3) = a(1) * b(2) - a(2) * b(1)
        
      end function cross_prod
    
   
    
    end program PoseEst
!   End of program
