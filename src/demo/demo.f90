!****************************************************************************
!
!  PROGRAM: PoseEst
!
!  PURPOSE:  Test the SoftPOSIT pose estimation algorithm
!
!****************************************************************************

    program PoseEst
    implicit none

    integer, parameter :: nbImagePts = 5, nbWorldPts = 6
    integer :: imagePtsData(nbImagePts,2), centerData(2), foundPose
    real :: worldPts(nbWorldPts,3), beta0, noiseStd, initRot(3,3), initTrans(3), focalLength, rot(3,3), trans(3)
    real :: clock1, clock2
    integer :: k
    
    imagePtsData = reshape( (/ 612, 486, 567, 476, 441, 117, 145, 206,&
234, 329 /), (/ nbImagePts, 2 /) )
    worldPts = reshape( (/ -3.75, 7.50, -3., 3., 0., 0., 0., 0., -5.0, 5.0, &
2.25, -2.25, 0.5, 2.75, -2.0, -2.0, -0.75, -0.75 /), (/ nbWorldPts, 3 /) )
    beta0 = 2.0E-4
    noiseStd = 10.0
    initRot = reshape( (/0.0, 1.0, 0.0, -1.0, 0.0, 0.0, 0.0, &
0.0, 1.0/), (/ 3, 3 /) )
    initTrans = (/0.0, 0.0, 30.0/)
    focalLength = 982.1
    centerData = (/376, 240/)
    
    !freq = getcpufrequency()
    !freq = 3.0D9
    !call getcpuclocks(clock1)
    call cpu_time(clock1)
    call softposit(rot, trans, foundPose, imagePtsData, worldPts, 5, 6, beta0, noiseStd,&
initRot, initTrans, focalLength, centerData)
    !call getcpuclocks(clock2)
    call cpu_time(clock2)

    ! Print the results
    print *, '#ImagePts = ', nbImagePts
    print *, '#WorldPts = ', nbWorldPts
    
    print *, 'ImagePts = '
    do k = 1, nbImagePts
       print *, imagePtsData(k,:)
    end do

    print *, 'WorldPts = '
    do k = 1, nbWorldPts
       print *, worldPts(k,:)
    end do
    
    print *, 'Rot = '
    do k = 1,3
       print *, rot(k,:)
    end do
    print *, 'Trans = '
    print *, trans
    print *, 'Eplased time = ', (clock2-clock1)

    end program PoseEst
!   End of program
