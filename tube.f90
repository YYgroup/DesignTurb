! Shiying Xiong, Peking University, 2019, ver-1.0 Generation of knotted fields
! Weiyu Shen, Peking University, 2021, ver-1.1 fix construction method of twisted cases
! Weiyu Shen, Peking University, 2021 Dec, ver-2.0 add sigmaFunc/etaFunc for bursting cases and fix bugs
! Weiyu Shen, Peking University, 2022 Jan, ver-2.1 add twist of different vortex surfaces 
! Weiyu Shen, Peking University, 2022 Apr, ver-2.2 unify the output format
! Weiyu Shen, Peking University, 2022 Jul, ver-3.0 update construction algorithm to closed flux tubes with arbitrary centerlines, vortex surface field, tube thickness distribution and local twist rate distribution
! Weiyu Shen, Peking University, 2022 Nov, ver-4.0 update construction algorithm to accommodate centerline across boundaries
! Weiyu Shen, Peking University, 2022 Nov, ver-4.1 use cubic splines to fit the scatter centerline
! Weiyu Shen, Peking University, 2023 Feb, ver-4.2 use cubic splines to fit multiple scatter centerlines
! Weiyu Shen, Peking University, 2023 Mar, ver-4.3 fix bugs
! Weiyu Shen, Peking University, 2023 Mar, ver-4.4 update construction algorithm for stop-points and straight lines
! Weiyu Shen, Peking University, 2023 Mar, ver-4.5 fix bugs
! Weiyu Shen, Peking University, 2023 Mar, ver-4.6 add the output for structure functions
! Weiyu Shen, Peking University, 2023 Mar, ver-4.7 fix bugs; add the calculation of the total length, total writhe, average vortex length density and average vortex volume ratio; update sigmaFunc/etaFunc for arc-length parameter
! Weiyu Shen, Peking University, 2023 Apr, ver-4.8 fix bugs
! Weiyu Shen, Peking University, 2023 Apr, ver-4.9 fix bugs; adaptive nt, ntc
! Weiyu Shen, Peking University, 2023 Apr, ver-4.10 change output field from dat file to plt file
! Weiyu Shen, Peking University, 2023 Jun, ver-4.11 add the 7box mode
! Weiyu Shen, Peking University, 2024 Apr, ver-4.12 add tpFunc
! Weiyu Shen, Peking University, 2024 Apr, ver-4.13 add Lx Ly Lz
! Weiyu Shen, Peking University, 2024 May, ver-4.14 adjust boundary processing
 
program turbvsf
    implicit none
    include 'mpif.h'
    include "fftw3.f"   
    integer status(MPI_STATUS_SIZE)
    integer, parameter :: nx = 128, ny = 128, nz = 128
    real * 8, parameter :: Gamma = 1.d0, nu = 1.d0, yt = 1.d0
    real * 8, parameter :: GammaR = 0.d0
    real * 8 eta
    real * 8 sigma, sigmax, dsigmads, pi, Rtube, pi2, pdx
    integer, parameter :: igst = 3
    real * 8 xstart, ystart, zstart, lx, ly, lz, lt
    real * 8 deviation, helicity, Etot, Diss, Enstrophy
    character * 200 name
    integer i, j, k, nzp, tii, tij, tnumber, ii, norder, nt, ntc, ntf
    real * 8 dx, dy, dz, dt, dv, dtc, dtf
    real * 8 tp, Rkappa, t0, t1, t2, t3, t4, t5, t6, t7, t8, s0, s1, s2, sb
    real * 8, dimension(3) :: xczeta, Tzeta, Nzeta, Bzeta, flag, meshxyz
    real * 8, dimension(3) :: czeta, dczeta, ddczeta, dddczeta, dcddczeta
    real * 8, dimension(3) :: ci0, ci1, ci2, dci0, dci1, tci0, tci1, cib
    real * 8 ndci0, ndci1
    real * 8 zeta, ndczeta, ndcddczeta, kappazeta, tauzeta, nxczeta
    integer ierr
    integer id, id_l, id_r, nproc
    double precision starttime, endtime
    integer * 8 planxf, planxb, planyf, planyb, planzf, planzb !!! fft plans
    real * 8, allocatable :: meshx(:,:,:), meshy(:,:,:), meshz(:,:,:)
    integer, allocatable :: kx(:), ky(:), kz(:), k2(:,:,:)
    real * 8, allocatable :: vorx(:,:,:), vory(:,:,:), vorz(:,:,:), phiv(:,:,:)
    real * 8, allocatable :: velx(:,:,:), vely(:,:,:), velz(:,:,:)
    double complex, allocatable :: specx(:, :, :), specy(:, :, :), specz(:, :, :)
    real, allocatable :: data_box(:,:,:,:)
    character*40, allocatable :: varname(:)
    real * 8, allocatable :: c(:, :), dc(:, :), ddc(:, :), dddc(:, :), dcddc(:, :)
    real * 8, allocatable :: kappa(:), tau(:), ndc(:), ndcddc(:), vc(:), Sn(:), r(:)
	integer nline, npoint, npointall, npointstart
	integer, allocatable :: npointlist(:)
	real * 8, allocatable :: cx(:), cy(:), cz(:), cze(:), ckx(:), cky(:), ckz(:), cxall(:), cyall(:), czall(:)
    real * 8 rho, theta, costh, sinth
    integer ntm, m
    real * 8 writhe, twist, length, slk, dist, cdx, cdy, cdz, dist2, writhetot, lengthtot, beta1, beta
    call MPI_Init(ierr)
    call MPI_Comm_rank(MPI_COMM_WORLD,id,ierr)
    call MPI_Comm_size(MPI_COMM_WORLD,nproc,ierr)
    starttime = MPI_WTIME()
    id_l = mod(nproc + id - 1,nproc)
    id_r= mod(id + 1,nproc)
    nzp = nz / nproc
    if(id==0) print * , 'nproc=', nproc
    if(id==0) print * , 'nzp=', nzp
	
    pi = 4.0d0 * datan(1.0d0)
    pi2 = 8.0d0 * datan(1.0d0)
    xstart = - pi
    ystart = - pi
    zstart = - pi
    lx = pi2
    ly = pi2
    lz = pi2
    lt = pi2
    dx = lx / nx
    dy = ly / ny
    dz = lz / nz
    dv = dx * dy * dz
    
    allocate(k2(nx,ny,nzp))
    allocate(kx(nx))
    allocate(ky(ny))
    allocate(kz(nzp))
	allocate(Sn(nx/2+1))
	allocate(r(nx/2+1))
    allocate(meshx(nx,ny,nzp))
    allocate(meshy(nx,ny,nzp))
    allocate(meshz(nx,ny,nzp))
    allocate(vorx(nx,ny,nzp))
    allocate(vory(nx,ny,nzp))
    allocate(vorz(nx,ny,nzp))
    allocate(velx(nx,ny,nzp))
    allocate(vely(nx,ny,nzp))
    allocate(velz(nx,ny,nzp))
    allocate(specx(nx,ny,nzp))
    allocate(specy(nx,ny,nzp))
    allocate(specz(nx,ny,nzp))
    allocate(phiv(1-igst:nx+igst,1-igst:ny+igst,1-igst:nzp+igst))
    call set_fft_plans(nx,ny,nz,planxf,planxb,planyf,planyb,planzf,planzb)
    call initialize_mesh(nx,ny,nzp,xstart,ystart,zstart,dx,dy,dz,meshx,meshy,meshz,id)
     
	
	if(id==0)then
         open(90,file='./QTCenterline.dat',status='unknown')
         read(90,*) nline
         print *,'==================================='
         print *,'number of lines:',nline
	endif
	call mpi_bcast(nline,1,MPI_INTEGER,0,mpi_comm_world,ierr)
	allocate(npointlist(nline))
	if(id==0)then
		 do i = 1, nline
			read(90,*) npointlist(i)
			print *,'==================================='
			print *,'line:',i
			print *,'npoint:',npointlist(i)
		 end do
	endif
	call mpi_bcast(npointlist,nline,MPI_INTEGER,0,mpi_comm_world,ierr)
	npointall = 0
	do i = 1, nline
		npointall = npointall + npointlist(i)
	end do
	allocate(cxall(npointall))
	allocate(cyall(npointall))
	allocate(czall(npointall))
	if(id==0)then
		 do i = 1, npointall
		 	read(90,*) cxall(i)
		 end do
         do i = 1, npointall
		 	read(90,*) cyall(i)
		 end do
		 do i = 1, npointall
		 	read(90,*) czall(i)
		 end do
         print *,'==================================='
         print *,' '
		 close(90)
    endif 
	call mpi_bcast(cxall,npointall,MPI_REAL8,0,mpi_comm_world,ierr)
	call mpi_bcast(cyall,npointall,MPI_REAL8,0,mpi_comm_world,ierr)
	call mpi_bcast(czall,npointall,MPI_REAL8,0,mpi_comm_world,ierr)
	
	npointstart = 0	
	vorx = 0.d0
    vory = 0.d0
    vorz = 0.d0
    phiv = 0.d0
	lengthtot = 0.d0
	writhetot = 0.d0
	!!!!!!!!!!
	if (id==0) then
		open(31, file='./output/curvedata.dat', status='unknown')
	endif

	do ii = 1, nline
	
	nt = 2 * npointlist(ii)
	ntc = 2 * npointlist(ii)
	ntf = 100
	dtc = lt / ntc
    dtf = lt / ntc / ntf
	
	npoint = npointlist(ii)
	allocate(cx(npoint))
	allocate(cy(npoint))
	allocate(cz(npoint))
	allocate(cze(npoint))
	allocate(ckx(npoint))
	allocate(cky(npoint))
	allocate(ckz(npoint))

	if (ii > 1) then
		npointstart = npointstart + npointlist(ii-1)
	endif
	do i = 1, npointlist(ii)
		cx(i) = cxall(npointstart+i)
		cy(i) = cyall(npointstart+i)
		cz(i) = czall(npointstart+i)
	end do

	
	!!!!!!get cze: initial zeta
	cze(1)=0.0d0
	do i = 2, npoint
		cdx=abs(cx(i)-cx(i-1))
		cdy=abs(cy(i)-cy(i-1))
		cdz=abs(cz(i)-cz(i-1))
		if (cdx > lx/2) cdx=cdx-lx
		if (cdy > ly/2) cdy=cdy-ly
		if (cdz > lz/2) cdz=cdz-lz
		dist = dsqrt(cdx**2.0d0+cdy**2.0d0+cdz**2.0d0)
		cze(i) = cze(i-1) + dist
	end do
	cdx=abs(cx(1)-cx(npoint))
	cdy=abs(cy(1)-cy(npoint))
	cdz=abs(cz(1)-cz(npoint))
	if (cdx > lx/2) cdx=cdx-lx
	if (cdy > ly/2) cdy=cdy-ly
	if (cdz > lz/2) cdz=cdz-lz
	dist = dsqrt(cdx**2.0d0+cdy**2.0d0+cdz**2.0d0)
	cze = cze * pi2 / (cze(npoint)+dist)
		
	!!!!!!!!get ckx cky ckz
	do i = 1, npoint
		if (i==1) then
			cdx=cx(2)-cx(npoint)
			cdy=cy(2)-cy(npoint)
			cdz=cz(2)-cz(npoint)
			dist = cze(2)-cze(npoint)+pi2
		else if (i==npoint) then
			cdx=cx(1)-cx(npoint-1)
			cdy=cy(1)-cy(npoint-1)
			cdz=cz(1)-cz(npoint-1)
			dist = cze(1)-cze(npoint-1)+pi2
		else
			cdx=cx(i+1)-cx(i-1)
			cdy=cy(i+1)-cy(i-1)
			cdz=cz(i+1)-cz(i-1)
			dist = cze(i+1)-cze(i-1)
		endif
		if (cdx > lx/2) cdx=cdx-lx
		if (cdy > ly/2) cdy=cdy-ly
		if (cdz > lz/2) cdz=cdz-lz
		if (cdx < -lx/2) cdx=cdx+lx
		if (cdy < -ly/2) cdy=cdy+ly
		if (cdz < -lz/2) cdz=cdz+lz
		ckx(i)=cdx/dist
		cky(i)=cdy/dist
		ckz(i)=cdz/dist
	end do
	
	allocate(c(nt, 3))
    allocate(dc(nt, 3))
    allocate(ddc(nt, 3))
    allocate(dddc(nt, 3))
    allocate(dcddc(nt, 3))
    allocate(vc(nt))
    allocate(ndc(nt))
    allocate(ndcddc(nt))
	allocate(kappa(nt))
    allocate(tau(nt))
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	m = 0
    
	dt = lt / nt
    do i = 1, nt
        tp = (i - 1.0d0) * dt
        call curve(m,tp,c(i, :),npoint,cx,cy,cz,cze,ckx,cky,ckz,lx,ly,lz)
		if (id==0) then
			write(31,*) c(i, 1), c(i, 2), c(i, 3)
		endif
    enddo
	
    call ptboundary(nt, dt, c, dc,lx,ly,lz)
    call pt(nt, dt, dc, ddc)
    call pt(nt, dt, ddc, dddc)
    call cross_t(nt, dc, ddc, dcddc)
    call norm_t(nt, dc, ndc)
    call norm_t(nt, dcddc, ndcddc)
    call dot_t(nt, dcddc, dddc, vc)
    tau = vc / (ndcddc ** 2.d0)
    writhe = 0.0d0
	kappa = ndcddc / (ndc ** 3.d0)
	Rkappa = 1.d0 / maxval(real(kappa))
    do i = 1, nt
        do j = 1, nt
            t1 = (c(i, 1) - c(j, 1)) * (dc(i, 2) * &
                &dc(j, 3) - dc(i, 3) * dc(j, 2))
            t2 = (c(i, 2) - c(j, 2)) * (dc(i, 3) * &
                &dc(j, 1) - dc(i, 1) * dc(j, 3))
            t3 = (c(i, 3) - c(j, 3)) * (dc(i, 1) * &
                &dc(j, 2) - dc(i, 2) * dc(j, 1))
            t4 = (c(i, 1) - c(j, 1)) ** 2.d0
            t5 = (c(i, 2) - c(j, 2)) ** 2.d0
            t6 = (c(i, 3) - c(j, 3)) ** 2.d0
            tp = t4 + t5 + t6
            if(tp /=0) then
                writhe = writhe + (t1 + t2 + t3) * dt / (tp ** (3.0d0 / 2.0d0))
            endif
        enddo
    enddo
	length = sum(ndc) * dt
	lengthtot = lengthtot + length
    writhe = writhe * dt  / 4.0d0 / pi
    twist = sum(tau * ndc) * dt / 2.0d0 / pi
	writhetot = writhetot + writhe
	slk = floor(writhe + twist + 0.5)
	if(id==0) then
		print *, '================'
		print *, 'curve number = ', ii
        print *, 'length = ', length
        print *, 'writhe = ', writhe
		print *, 'twist = ', twist
		print *, 'self-linking = ', slk
    endif
    deallocate(c)
    deallocate(dc)
    deallocate(ddc)
    deallocate(dddc)
    deallocate(dcddc)
    deallocate(vc)
    deallocate(ndc)
    deallocate(ndcddc)
	deallocate(kappa)
    deallocate(tau)
    	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
    endtime=MPI_WTIME()
    if(id==0) print *,'initialize time is',endtime-starttime
    starttime=MPI_WTIME()
    
	zeta=0.0d0
    call sigmaFunc(zeta,length,sigma,sigmax)
	Rtube = 5.0d0 * sigmax
    pdx = 0.0d0
 
	!!!!!!!!!! Face boundary processing : do m = 0, 6 ; Full boundary processing : do m = 0, 26
	!!!!!!!!!! If no boundary processing is performed, please comment out the m loop
	do m = 0, 6!26
    do k = 1, nzp
        do j = 1, ny
            do i = 1, nx
                tnumber = 0
                do tii = 1, ntc
                    call curve(m,dtc * (tii - 1), ci0,npoint,cx,cy,cz,cze,ckx,cky,ckz,lx,ly,lz)
                    call curve(m,dtc * tii, ci1,npoint,cx,cy,cz,cze,ckx,cky,ckz,lx,ly,lz)
                    call curve(m,dtc * (tii + 1), ci2,npoint,cx,cy,cz,cze,ckx,cky,ckz,lx,ly,lz)
                    call cminus(ci1,ci0,dci0,lx,ly,lz)
					call cminus(ci2,ci1,dci1,lx,ly,lz)
                    call norm_zeta(dci0,ndci0)
                    call norm_zeta(dci1,ndci1)
                    tci0 = dci0 / ndci0
                    tci1 = dci1 / ndci1
					call distflag(ci1,ci0,flag,lx,ly,lz)
					if (flag(1) /= 0.0d0 .or. flag(2) /= 0.0d0 .or. flag(3) /= 0.0d0) then
						dist = dsqrt((meshx(i, j, k) - ci0(1)) ** 2 + &
                            &(meshy(i, j, k) - ci0(2)) ** 2 + &
                            &(meshz(i, j, k) - ci0(3)) ** 2)
						dist2 = dsqrt((meshx(i, j, k) - ci1(1)) ** 2 + &
                            &(meshy(i, j, k) - ci1(2)) ** 2 + &
                            &(meshz(i, j, k) - ci1(3)) ** 2)
						if (dist<=dist2) then
							ci1 = ci1 - flag
						else
							ci0 = ci0 + flag
						endif
					endif
                    t0 = (meshx(i, j, k) - ci0(1)) * tci0(1) + &
                        &(meshy(i, j, k) - ci0(2)) * tci0(2) +  &
                        &(meshz(i, j, k) - ci0(3)) * tci0(3)
                    t1 = (meshx(i, j, k) - ci1(1)) * tci1(1) + &
                        &(meshy(i, j, k) - ci1(2)) * tci1(2) +  &
                        &(meshz(i, j, k) - ci1(3)) * tci1(3)
                    if (t0 >= 0 .and. t1 < 0) then
                        rho = dsqrt((meshx(i, j, k) - ci0(1)) ** 2 + &
                            &(meshy(i, j, k) - ci0(2)) ** 2 + &
                            &(meshz(i, j, k) - ci0(3)) ** 2 - t0 ** 2)
                        if (rho < Rtube) then
                            do tij = - 500, ntf + 500
                                s0 = dtc * (tii - 1) + dtf * (tij - 1.d0)
                                s1 = dtc * (tii - 1) + dtf * (tij + 0.d0)
                                s2 = dtc * (tii - 1) + dtf * (tij + 1.d0)
                                call curve(m,s0, ci0,npoint,cx,cy,cz,cze,ckx,cky,ckz,lx,ly,lz)
                                call curve(m,s1, ci1,npoint,cx,cy,cz,cze,ckx,cky,ckz,lx,ly,lz)
                                call curve(m,s2, ci2,npoint,cx,cy,cz,cze,ckx,cky,ckz,lx,ly,lz)
                                call cminus(ci1,ci0,dci0,lx,ly,lz)
								call cminus(ci2,ci1,dci1,lx,ly,lz)
                                call norm_zeta(dci0,ndci0)
                                call norm_zeta(dci1,ndci1)
                                tci0 = dci0 / ndci0
                                tci1 = dci1 / ndci1
								call distflag(ci1,ci0,flag,lx,ly,lz)
								if (m > 0) then
									if (flag(1) /= 0.0d0 .or. flag(2) /= 0.0d0 .or. flag(3) /= 0.0d0) then
										t0 = -1.0d0
										t1 = 1.0d0
									else
										t0 = (meshx(i, j, k) - ci0(1)) * tci0(1) + &
											&(meshy(i, j, k) - ci0(2)) * tci0(2) +  &
											&(meshz(i, j, k) - ci0(3)) * tci0(3)
										t1 = (meshx(i, j, k) - ci1(1)) * tci1(1) + &
											&(meshy(i, j, k) - ci1(2)) * tci1(2) +  &
											&(meshz(i, j, k) - ci1(3)) * tci1(3)
									endif
								else
									if (flag(1) /= 0.0d0 .or. flag(2) /= 0.0d0 .or. flag(3) /= 0.0d0) then
										dist = dsqrt((meshx(i, j, k) - ci0(1)) ** 2 + &
											&(meshy(i, j, k) - ci0(2)) ** 2 + &
											&(meshz(i, j, k) - ci0(3)) ** 2)
										dist2 = dsqrt((meshx(i, j, k) - ci1(1)) ** 2 + &
											&(meshy(i, j, k) - ci1(2)) ** 2 + &
											&(meshz(i, j, k) - ci1(3)) ** 2)
										if (dist<=dist2) then
											ci1 = ci1 - flag
										else
											ci0 = ci0 + flag
										endif
									endif
									t0 = (meshx(i, j, k) - ci0(1)) * tci0(1) + &
										&(meshy(i, j, k) - ci0(2)) * tci0(2) +  &
										&(meshz(i, j, k) - ci0(3)) * tci0(3)
									t1 = (meshx(i, j, k) - ci1(1)) * tci1(1) + &
										&(meshy(i, j, k) - ci1(2)) * tci1(2) +  &
										&(meshz(i, j, k) - ci1(3)) * tci1(3)
								endif			
                                if (t0 >= 0.d0 .and. t1 < 0.d0) then
                                    tnumber = tnumber + 1
                                    ! if(tnumber > 1) print *, tnumber
                                    t1 = (meshx(i, j, k) - ci1(1)) * tci0(1) + &
                                        &(meshy(i, j, k) - ci1(2)) * tci0(2) +  &
                                        &(meshz(i, j, k) - ci1(3)) * tci0(3)
                                    zeta = (s0 * t0 - s1 * t1) / ndci0
                                    call curve(m,zeta, czeta,npoint,cx,cy,cz,cze,ckx,cky,ckz,lx,ly,lz)
                                    call d1curve(m,zeta, dczeta,npoint,cx,cy,cz,cze,ckx,cky,ckz,lx,ly,lz)
                                    call d2curve(m,zeta, ddczeta,npoint,cx,cy,cz,cze,ckx,cky,ckz,lx,ly,lz)
                                    call d3curve(m,zeta, dddczeta,npoint,cx,cy,cz,cze,ckx,cky,ckz,lx,ly,lz)
                                    call cross_zeta(dczeta,ddczeta,dcddczeta)
                                    call norm_zeta(dczeta,ndczeta)
                                    call norm_zeta(dcddczeta,ndcddczeta)
                                    Tzeta = dczeta / ndczeta
                                    Bzeta = dcddczeta / ndcddczeta
                                    call cross_zeta(Bzeta,Tzeta,Nzeta)
                                    kappazeta = ndcddczeta / (ndczeta ** 3.d0)
                                    call dot_zeta(dcddczeta,dddczeta,t4) 
                                    tauzeta = t4 / (ndcddczeta ** 2.d0)
                                    meshxyz(1) = meshx(i, j, k)
                                    meshxyz(2) = meshy(i, j, k)
                                    meshxyz(3) = meshz(i, j, k)
									call cminus(meshxyz,czeta,xczeta,lx,ly,lz)
                                    call norm_zeta(xczeta,nxczeta)
                                    call dot_zeta(xczeta,Nzeta,t4) 
                                    costh = t4 / (nxczeta + 1.0d-15)
                                    call dot_zeta(xczeta,Bzeta,t4)
                                    sinth = t4 / (nxczeta + 1.0d-15)
									theta = datan2(sinth,costh)
                                    rho = nxczeta
									call sigmaFunc(zeta,length,sigma,sigmax)
                                    !tp = dexp(- ((rho - pdx) ** 2.d0) / 2.d0 / (sigma ** 2.0d0))/ (pi2 * (sigma ** 2.0d0))
									call tpFunc(zeta,tp,rho,pdx,sigma,pi2,theta)
                                    t5 = - sinth * rho / (1.d0 - kappazeta * rho * costh)
                                    t6 = costh * rho / (1.d0 - kappazeta * rho * costh)
									call dsigmadsFunc(m,zeta,length,dsigmads,npoint,cx,cy,cz,cze,ckx,cky,ckz,lx,ly,lz)
									t7 = costh * rho * dsigmads / sigma / (1.d0 - kappazeta * rho * costh)
									t8 = sinth * rho * dsigmads / sigma / (1.d0 - kappazeta * rho * costh)
									phiv(i, j, k) = phiv(i, j, k) + (dexp(- ((rho - pdx) ** 2.d0) / 2.d0 / (sigma ** 2.0d0)))
									call etaFunc(zeta,length,phiv(i, j, k),eta)
                                    vorx(i,j,k) = vorx(i,j,k) + tp * Gamma * (Tzeta(1) + &
                                        & eta * (t5 * Nzeta(1) + t6 * Bzeta(1)) + (t7 * Nzeta(1) + t8 * Bzeta(1)) )
                                    vory(i,j,k) = vory(i,j,k) + tp * Gamma * (Tzeta(2) + &
                                        & eta * (t5 * Nzeta(2) + t6 * Bzeta(2)) + (t7 * Nzeta(2) + t8 * Bzeta(2)) )
                                    vorz(i,j,k) = vorz(i,j,k) + tp * Gamma * (Tzeta(3) + &
                                        & eta * (t5 * Nzeta(3) + t6 * Bzeta(3)) + (t7 * Nzeta(3) + t8 * Bzeta(3)) )
                                endif
                            enddo
                        endif
                    endif
                enddo
            enddo
        enddo
    enddo
	enddo
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	deallocate(cx)
	deallocate(cy)
	deallocate(cz)
	deallocate(cze)
	deallocate(ckx)
	deallocate(cky)
	deallocate(ckz)
	enddo
	
	if (id == 0) then
		close(31)
	endif
	
    endtime=MPI_WTIME()
    if(id==0) print *,'construction time is',endtime-starttime
    starttime=MPI_WTIME()
	
	beta1 = 0.0d0
	do k = 1, nzp
        do j = 1, ny
            do i = 1, nx
			 if (phiv(i,j,k) > dexp(-12.5d0)) beta1=beta1+1.0d0
			enddo
		enddo	
	enddo
	call MPI_REDUCE(beta1,beta,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
	beta = beta / nx / ny / nz
	
    call wavenumber(nx,ny,nz,nzp,kx,ky,kz,k2,id)
    call cross_vector(vorx,vory,vorz,velx,vely,velz,-1,nx,ny,nzp,&
       &kx,ky,kz,k2,planxf,planyf,planzf,planxb,planyb,planzb,id,nproc)
    if(id==0) then
       open(49, file='./output/spectrum.dat',status='unknown')
    endif
    call spectrum(nx,ny,nzp,velx,vely,velz,49,k2,id,nproc,planxf,planyf,planzf)
    if(id==0) then
       close(49)
    endif
    call wall_phiv(nx,ny,nzp,phiv,igst,id_l,id_r)
    call wall_phiv(nx,ny,nzp,phiv,igst,id_l,id_r)
    call wall_phiv(nx,ny,nzp,phiv,igst,id_l,id_r)
    call cal_deviation(vorx,vory,vorz,phiv,nx,ny,nz,nzp,deviation,dx,dy,dz,id,igst)
    allocate(data_box(nx,ny,nzp,8))
    data_box(:,:,:,1)=meshx
    data_box(:,:,:,2)=meshy
    data_box(:,:,:,3)=meshz
    data_box(:,:,:,4)=vorx
    data_box(:,:,:,5)=vory
    data_box(:,:,:,6)=vorz
    data_box(:,:,:,7)=phiv(1:nx,1:ny,1:nzp)
    data_box(:,:,:,8)=vorx * velx + vory * vely + vorz * velz
	!data_box(:,:,:,8)=dsqrt(vorx * vorx + vory * vory + vorz * vorz)
    name='./output/vorticity.plt'   
    allocate(varname(8))
    varname(1)='x'
    varname(2)='y'
    varname(3)='z'
    varname(4)='w1'
    varname(5)='w2'
    varname(6)='w3'
    varname(7)='phiv'
    varname(8)='h'
    call output_v(name,varname,nx,ny,nzp,data_box,8,id,nproc)
    deallocate(data_box)
    deallocate(varname)

    allocate(data_box(nx,ny,nzp,8))
    data_box(:,:,:,1)=meshx
    data_box(:,:,:,2)=meshy
    data_box(:,:,:,3)=meshz
    data_box(:,:,:,4)=velx
    data_box(:,:,:,5)=vely
    data_box(:,:,:,6)=velz
    data_box(:,:,:,7)=phiv(1:nx,1:ny,1:nzp)
    data_box(:,:,:,8)=vorx * velx + vory * vely + vorz * velz
    name='./output/velocity.plt'   
    allocate(varname(8))
    varname(1)='x'
    varname(2)='y'
    varname(3)='z'
    varname(4)='v1'
    varname(5)='v2'
    varname(6)='v3'
    varname(7)='phiv'
    varname(8)='h'
    call output_v(name,varname,nx,ny,nzp,data_box,8,id,nproc)
    deallocate(data_box)
    deallocate(varname)

    call fourier_forward(velx,specx,nx,ny,nz,planxf,planyf,planzf,id,nproc)
    call fourier_forward(vely,specy,nx,ny,nz,planxf,planyf,planzf,id,nproc)
    call fourier_forward(velz,specz,nx,ny,nz,planxf,planyf,planzf,id,nproc)
    allocate(data_box(nx,ny,nzp,6))
    data_box(:,:,:,1)=real(specx)
    data_box(:,:,:,2)=real(specy)
    data_box(:,:,:,3)=real(specz)
    data_box(:,:,:,4)=aimag(specx)
    data_box(:,:,:,5)=aimag(specy)
    data_box(:,:,:,6)=aimag(specz)
    name='./output/specVelocity.dat'   
    allocate(varname(6))
    varname(1)='rvx'
    varname(2)='rvy'
    varname(3)='rvz'
    varname(4)='ivx'
    varname(5)='ivy'
    varname(6)='ivz'
    call output_v(name,varname,nx,ny,nzp,data_box,6,id,nproc)
    deallocate(data_box)
    deallocate(varname)

    call cal_helicity(vorx,vory,vorz,velx,vely,velz,nx,ny,nz,nzp,helicity,dv,id)
    call cal_helicity(velx,vely,velz,velx,vely,velz,nx,ny,nz,nzp,Etot,dv,id)
    Etot = Etot / 2.d0
    call cal_helicity(vorx,vory,vorz,vorx,vory,vorz,nx,ny,nz,nzp,Diss,dv,id)
    Enstrophy = Diss / 2.0d0
	Diss = Diss / lx / ly / lz
    call destroy_fft_plans(planxf,planxb,planyf,planyb,planzf,planzb)
    if(id==0) then
        open(19, file='./output/parameter.dat ',status='unknown')
            write(19, *) 'nx, ny, nz = ', nx, ny, nz
            write(19, *) 'Total length of curve = ', lengthtot
			write(19, *) 'Average vortex length density = ', lengthtot / lx / ly / lz
			write(19, *) 'Average vortex volume ratio = ', beta
            write(19, *) 'Minimum radius of curvature = ', Rkappa
            write(19, *) 'Radius of vortex tube = ', Rtube
            write(19, *) 'Maximal standard deviation of flux distribution = ', sigmax
            write(19, *) 'VSF deviation of vorticity magnitude = ', deviation
            write(19, *) 'Total helicity = ', helicity
            write(19, *) 'Total writhe = ', writhetot
            write(19, *) 'Gamma = ', Gamma, GammaR
            write(19, *) 'eta = ', eta
            write(19, *) 'Etot = ', Etot
            write(19, *) 'Diss/nu = ', Diss
			write(19, *) 'Enstrophy = ', Enstrophy			
            write(19, *) '(helicity - (Gamma ** 2) * (writhetot + eta * lengthtot / pi2)) / helicity  = ', &
                &(helicity - (Gamma ** 2) * (writhetot + eta * lengthtot / pi2)) / helicity
        close(19)
    endif
	norder = 2
	call getStructureFunc(Sn,r,norder,velx,vely,velz,lx,ly,lz,nx,ny,nz,nzp,id)
	if(id==0) then
        open(20, file='./output/structureFunc_2.dat ',status='unknown')
            do i=1,nx/2+1
				write(20, *) r(i), Sn(i)
			enddo
        close(20)
    endif
	norder = 3
	call getStructureFunc(Sn,r,norder,velx,vely,velz,lx,ly,lz,nx,ny,nz,nzp,id)
	if(id==0) then
        open(21, file='./output/structureFunc_3.dat ',status='unknown')
            do i=1,nx/2+1
				write(21, *) r(i), Sn(i)
			enddo
        close(21)
    endif
	norder = 4
	call getStructureFunc(Sn,r,norder,velx,vely,velz,lx,ly,lz,nx,ny,nz,nzp,id)
	if(id==0) then
        open(22, file='./output/structureFunc_4.dat ',status='unknown')
            do i=1,nx/2+1
				write(22, *) r(i), Sn(i)
			enddo
        close(22)
    endif
	norder = 5
	call getStructureFunc(Sn,r,norder,velx,vely,velz,lx,ly,lz,nx,ny,nz,nzp,id)
	if(id==0) then
        open(23, file='./output/structureFunc_5.dat ',status='unknown')
            do i=1,nx/2+1
				write(23, *) r(i), Sn(i)
			enddo
        close(23)
    endif
    if(id==0) then
        print *, 'nx, ny, nz = ', nx, ny, nz
        print *, 'Total length of curve = ', lengthtot
		print *, 'Average vortex length density = ', lengthtot / lx / ly / lz
		print *, 'Average vortex volume ratio = ', beta
        print *, 'Minimum radius of curvature = ', Rkappa
        print *, 'Radius of vortex tube = ', Rtube
        print *, 'Maximal standard deviation of flux distribution = ', sigmax
        print *, 'Total writhe = ', writhetot
        print *, 'Gamma = ', Gamma, GammaR
        print *, 'eta = ', eta
        print *, 'VSF deviation of vorticity magnitude = ', deviation
        print *, 'Total helicity = ', helicity
        print *, 'expected helicity =', (Gamma ** 2) * (writhetot + eta * lengthtot / pi2)
        print *, 'Etot = ', Etot
        print *, 'Diss/nu = ', Diss
		print *, 'Enstrophy = ', Enstrophy
        print *, '(helicity - (Gamma ** 2) * (writhetot + eta * lengthtot / pi2)) / helicity  = ', &
        &(helicity - (Gamma ** 2) * (writhetot + eta * lengthtot / pi2)) / helicity
    endif
    deallocate(k2)
    deallocate(kx)
    deallocate(ky)
    deallocate(kz)
	deallocate(Sn)
	deallocate(r)
    deallocate(meshx)
    deallocate(meshy)
    deallocate(meshz)
    deallocate(vorx)
    deallocate(vory)
    deallocate(vorz)
    deallocate(velx)
    deallocate(vely)
    deallocate(velz)
    deallocate(phiv)
	deallocate(cxall)
	deallocate(cyall)
	deallocate(czall)
	deallocate(npointlist)
    endtime=MPI_WTIME()
    if(id==0) print *,'output time is',endtime-starttime
    call MPI_Finalize(ierr)
end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
