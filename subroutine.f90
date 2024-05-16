subroutine tpFunc(zeta,tp,rho,pdx,sigma,pi2,theta)
    real * 8 zeta,tp,rho,pdx,sigma,pi,pi2,theta,h,f,omega0,Omega,OmegaPrime,omega2,t
	pi = 4.0d0 * datan(1.0d0)
    pi2 = 8.0d0 * datan(1.0d0)
	h = 0.5d0
	nus = 4.0d-5
	! Guassian
	!tp = dexp(- ((rho - pdx) ** 2.d0) / 2.d0 / (sigma ** 2.0d0))/ (pi2 * (sigma ** 2.0d0))
	t = 0.4d0
	!
	!f = h * rho**2.0d0/(2.0d0*sigma**2.0d0) * dexp(-rho**2.0d0/(2.0d0*sigma**2.0d0) ) / (2.0d0*(1.0d0+h)*pi*(2.0d0*sigma**2.0d0))
	!omega0 = dexp(-rho**2.0d0/(2.0d0*sigma**2.0d0)) / ((1.0d0+h)*pi*(2.0d0*sigma**2.0d0)) + 2.0d0 * f
	!Omega = (1.0d0+h-(1.0d0+h+h*(rho**2.0d0/(2.0d0*sigma**2.0d0))) * dexp(-rho**2.0d0/(2.0d0*sigma**2.0d0))) / (2.0d0*(1.0d0+h)*pi*(2.0d0*sigma**2.0d0))
	!OmegaPrime = (omega0-2.0d0*Omega)/rho
	!tp = omega0 + 2.0d0 * f * dcos(2.0d0*(theta - Omega * t)) * dexp(-4 * OmegaPrime**2 * nus * t**3 / 3)
	tp = dexp(- ((rho - pdx) ** 2.d0) / 2.d0 / (sigma ** 2.0d0))/ (pi2 * (sigma ** 2.0d0))
end subroutine

subroutine getStructureFunc(Sn,r,norder,velx,vely,velz,lx,ly,lz,nx,ny,nz,nzp,id)
	implicit none
    include 'mpif.h'
    real * 8 dr, lx, ly, lz
	integer i, j, k, m, p, norder
	integer nx, ny, nzp, nz, ierr, id
	real*8, dimension (nx,ny,nzp) :: velx, vely, velz
	real*8, dimension (nx/2+1) :: Sn1, Sn, r
	dr=lx/nx
	do m=1,nx/2+1
		r(m)=(m-1)*dr
	enddo
	Sn1 = 0.0d0
	do k=1,nzp
		do j=1,ny
			do i=1,nx
				do m=1, nx/2+1
					p = i+(m-1)
					if (p > nx) p = p - nx
					Sn1(m) = Sn1(m) + (velx(p,j,k) - velx(i,j,k))**norder
				enddo
			enddo
		enddo
	enddo
	
	call MPI_REDUCE(Sn1,Sn,nx/2+1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
	
	Sn = Sn / (nx*ny*nz)

end subroutine

subroutine sigmaFunc(zeta,length,sigma,sigmax)
    real * 8 zeta, sigma, sigma0, sigmax, pi, pi2, beta, length, omega
	pi = 4.0d0 * datan(1.0d0)
    pi2 = 8.0d0 * datan(1.0d0)
    sigma0 = 3.0d-2
	beta = 1.5d0
	!!!!!!!asy
	! if(zeta<-pi/2.0d0) then
	   ! sigma=sigma0 * (2.0d0 + dsin(2.0d0 * zeta - pi / 2.0d0))
    ! else
       ! if(zeta<pi/2.0d0) then
	      ! sigma=sigma0 * (2.0d0 - dsin(zeta))
	   ! else
          ! sigma=sigma0
       ! endif
    ! endif
	!!!!!!!sy
	! if(zeta<0.0d0) then
	   ! sigma=sigma0 * ((beta+1.0d0)/2.0d0 + (beta-1.0d0)/2.0d0 * dsin(2.0d0 * zeta - pi / 2.0d0))
    ! else
       ! sigma=sigma0
    ! endif
	!!!!!!!vortex ring
	!sigma=sigma0
	omega = floor(length / pi2 + 0.5d0)
	sigma=sigma0 * (1.0d0+ 0.5d0*dsin(omega*(zeta)))
	sigmax = sigma0 * beta
end subroutine

subroutine etaFunc(zeta,length,phivijk,eta)
    real * 8 eta, phivijk, zeta, sigma, Tw, pi, pi2, length
	pi = 4.0d0 * datan(1.0d0)
    pi2 = 8.0d0 * datan(1.0d0)
    Tw = 0.0d0
    !sigma = 1.0d0 / 1.6d0 / dsqrt(pi2)
	!if(zeta<0.0d0) then
    !    eta=(dsqrt(pi2)*Tw)/sigma * dexp( - (zeta ** 2.d0) / 2.d0 / (sigma ** 2.0d0))
    !else
    !    eta=(dsqrt(pi2)*Tw)/sigma * dexp( - ((pi2-zeta) ** 2.d0) / 2.d0 / (sigma ** 2.0d0))
    !endif
	eta = 0.0d0
	!eta=10.0d0*(3.0d0*phivijk-1.0d0)*dsin((zeta+pi)/2.0d0)
	!eta = Tw*phivijk
	!eta = Tw*phivijk*dcos(zeta/2.0d0)
	!eta = Tw*dsin(pi*phivijk)*dcos(zeta/2.0d0)
end subroutine

subroutine getsplinecoe(a,b,c,d,zeta1,zeta2,x1,x2,dx1,dx2,lx)
    real * 8 a,b,c,d,zeta1,zeta2,x1,x2,x1p,x2p,dx1,dx2,lx
	x1p = x1
	x2p = x2
	if (x1-x2 > lx/2) x2p=x2+lx
	if (x1-x2 < -lx/2) x1p=x1+lx
	a = -2.0d0/(zeta1-zeta2)*x1p + 2.0d0/(zeta1-zeta2)*x2p + dx1 + dx2
	b = 3.0d0*(zeta1+zeta2)/(zeta1-zeta2)*x1p - 3.0d0*(zeta1+zeta2)/(zeta1-zeta2)*x2p + (-zeta1-2.0d0*zeta2)*dx1 + (-2.0d0*zeta1-zeta2)*dx2
	c = -6.0d0*zeta1*zeta2/(zeta1-zeta2)*x1p + 6.0d0*zeta1*zeta2/(zeta1-zeta2)*x2p + (2.0d0*zeta1*zeta2+zeta2**2.0d0)*dx1 + (zeta1**2.0d0+2.0d0*zeta1*zeta2)*dx2
	d = (3.0d0*zeta1*zeta2**2.0d0-zeta2**3.0d0)/(zeta1-zeta2)*x1p + (zeta1**3.0d0-3.0d0*zeta1**2.0d0*zeta2)/(zeta1-zeta2)*x2p - zeta1*zeta2**2.0d0*dx1 - zeta1**2.0d0*zeta2*dx2
	a = a /(zeta1-zeta2)**2.0d0
	b = b /(zeta1-zeta2)**2.0d0
	c = c /(zeta1-zeta2)**2.0d0
	d = d /(zeta1-zeta2)**2.0d0
end subroutine

subroutine linetocurve(a,b,c,zeta1,zeta2,x1,x2)
    real * 8 a,b,c,d,zeta1,zeta2,x1,x2,pi,pi2
	pi = 4.0d0 * datan(1.0d0)
    pi2 = 8.0d0 * datan(1.0d0)
	a = (zeta2-zeta1)/1.0d10
	b = pi2/(zeta2-zeta1)
	c = -pi*(3.0d0*zeta1+zeta2)/2.0d0/(zeta2-zeta1)
end subroutine

subroutine curve(m,zeta0,c3,npoint,cx,cy,cz,cze,ckx,cky,ckz,lx,ly,lz)
	real * 8 zeta, zeta0, pi, pi2, a, b, c, d, flag, lx, ly, lz
    real * 8, dimension(3) :: c3
	integer npoint, i, m
	real * 8, dimension(npoint) :: cx, cy, cz, cze, ckx, cky, ckz
	pi = 4.0d0 * datan(1.0d0)
    pi2 = 8.0d0 * datan(1.0d0)
	zeta = dmod(zeta0,pi2)
	if (zeta<0.0d0) zeta=zeta+pi2
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!
	if (cze(npoint) <= zeta) then
	    flag=0.0d0
		call getsplinecoe(a,b,c,d,cze(npoint),pi2,cx(npoint),cx(1),ckx(npoint),ckx(1),lx)
		if ((abs(a)<1.0d-10) .and. (abs(b)<1.0d-10)) flag=1.0d0
		c3(1)=a*zeta**3+b*zeta**2+c*zeta+d
		call getsplinecoe(a,b,c,d,cze(npoint),pi2,cy(npoint),cy(1),cky(npoint),cky(1),ly)
		if ((abs(a)<1.0d-10) .and. (abs(b)<1.0d-10) .and. (flag==1.0d0)) flag=2.0d0
		c3(2)=a*zeta**3+b*zeta**2+c*zeta+d
		call getsplinecoe(a,b,c,d,cze(npoint),pi2,cz(npoint),cz(1),ckz(npoint),ckz(1),lz)
		if ((abs(a)<1.0d-10) .and. (abs(b)<1.0d-10) .and. (flag==2.0d0)) flag=3.0d0
		c3(3)=a*zeta**3+b*zeta**2+c*zeta+d
		if (flag==3.0d0) then
			! print *,'straight line'
			call linetocurve(a,b,c,cze(npoint),pi2,cx(npoint),cx(1))
			c3(1)=c3(1)+a*dsin(b*zeta+c)+a
			call linetocurve(a,b,c,cze(npoint),pi2,cy(npoint),cy(1))
			c3(2)=c3(2)+a*dsin(b*zeta+c)+a
			call linetocurve(a,b,c,cze(npoint),pi2,cz(npoint),cz(1))
			c3(3)=c3(3)+a*dsin(b*zeta+c)+a
		endif
	else
		do i = 1, npoint-1
		if ((cze(i) <= zeta) .and. (cze(i+1) > zeta)) then
			flag=0.0d0
			call getsplinecoe(a,b,c,d,cze(i),cze(i+1),cx(i),cx(i+1),ckx(i),ckx(i+1),lx)
			if ((abs(a)<1.0d-10) .and. (abs(b)<1.0d-10)) flag=1.0d0
			c3(1)=a*zeta**3+b*zeta**2+c*zeta+d
			call getsplinecoe(a,b,c,d,cze(i),cze(i+1),cy(i),cy(i+1),cky(i),cky(i+1),ly)
			if ((abs(a)<1.0d-10) .and. (abs(b)<1.0d-10) .and. (flag==1.0d0)) flag=2.0d0
			c3(2)=a*zeta**3+b*zeta**2+c*zeta+d
			call getsplinecoe(a,b,c,d,cze(i),cze(i+1),cz(i),cz(i+1),ckz(i),ckz(i+1),lz)
			if ((abs(a)<1.0d-10) .and. (abs(b)<1.0d-10) .and. (flag==2.0d0)) flag=3.0d0
			c3(3)=a*zeta**3+b*zeta**2+c*zeta+d
			if (flag==3.0d0) then
				! print *,'straight line'
				call linetocurve(a,b,c,cze(i),cze(i+1),cx(i),cx(i+1))
				c3(1)=c3(1)+a*dsin(b*zeta+c)+a
				call linetocurve(a,b,c,cze(i),cze(i+1),cy(i),cy(i+1))
				c3(2)=c3(2)+a*dsin(b*zeta+c)+a
				call linetocurve(a,b,c,cze(i),cze(i+1),cz(i),cz(i+1))
				c3(3)=c3(3)+a*dsin(b*zeta+c)+a
			endif
		endif
		end do
	endif
	! if (flag<3.0d0) print *,flag
	if (c3(1) >= lx/2) c3(1)=c3(1)-lx
	if (c3(2) >= ly/2) c3(2)=c3(2)-ly
	if (c3(3) >= lz/2) c3(3)=c3(3)-lz
	if (c3(1) < -lx/2) c3(1)=c3(1)+lx
	if (c3(2) < -ly/2) c3(2)=c3(2)+ly
	if (c3(3) < -lz/2) c3(3)=c3(3)+lz
	!!!!!!!!!!!!!!!!!!!
	if (m == 1) then
		c3(1) = c3(1)
		c3(2) = c3(2)
		c3(3) = c3(3) + lz
	endif
	if (m == 2) then
		c3(1) = c3(1)
		c3(2) = c3(2)
		c3(3) = c3(3) - lz
	endif
	if (m == 3) then
		c3(1) = c3(1)
		c3(2) = c3(2) + ly
		c3(3) = c3(3)
	endif
	if (m == 4) then
		c3(1) = c3(1)
		c3(2) = c3(2) - ly
		c3(3) = c3(3)
	endif
	if (m == 5) then
		c3(1) = c3(1) + lx
		c3(2) = c3(2)
		c3(3) = c3(3)
	endif
	if (m == 6) then
		c3(1) = c3(1) - lx
		c3(2) = c3(2)
		c3(3) = c3(3)
	endif
	!!!!!!
	if (m == 7) then
		c3(1) = c3(1)
		c3(2) = c3(2) + ly
		c3(3) = c3(3) + lz
	endif
	if (m == 8) then
		c3(1) = c3(1)
		c3(2) = c3(2) + ly
		c3(3) = c3(3) - lz
	endif
	if (m == 9) then
		c3(1) = c3(1)
		c3(2) = c3(2) - ly
		c3(3) = c3(3) + lz
	endif
	if (m == 10) then
		c3(1) = c3(1)
		c3(2) = c3(2) - ly
		c3(3) = c3(3) - lz
	endif
	if (m == 11) then
		c3(1) = c3(1) + lx
		c3(2) = c3(2)
		c3(3) = c3(3) + lz
	endif
	if (m == 12) then
		c3(1) = c3(1) + lx
		c3(2) = c3(2)
		c3(3) = c3(3) - lz
	endif
	if (m == 13) then
		c3(1) = c3(1) + lx
		c3(2) = c3(2) + ly
		c3(3) = c3(3)
	endif
	if (m == 14) then
		c3(1) = c3(1) + lx
		c3(2) = c3(2) - ly
		c3(3) = c3(3)
	endif
	if (m == 15) then
		c3(1) = c3(1) - lx
		c3(2) = c3(2)
		c3(3) = c3(3) + lz
	endif
	if (m == 16) then
		c3(1) = c3(1) - lx
		c3(2) = c3(2)
		c3(3) = c3(3) - lz
	endif
	if (m == 17) then
		c3(1) = c3(1) - lx
		c3(2) = c3(2) + ly
		c3(3) = c3(3)
	endif
	if (m == 18) then
		c3(1) = c3(1) - lx
		c3(2) = c3(2) - ly
		c3(3) = c3(3)
	endif
	!!!!!!
	if (m == 19) then
		c3(1) = c3(1) + lx
		c3(2) = c3(2) + ly
		c3(3) = c3(3) + lz
	endif
	if (m == 20) then
		c3(1) = c3(1) + lx
		c3(2) = c3(2) + ly
		c3(3) = c3(3) - lz
	endif
	if (m == 21) then
		c3(1) = c3(1) + lx
		c3(2) = c3(2) - ly
		c3(3) = c3(3) + lz
	endif
	if (m == 22) then
		c3(1) = c3(1) + lx
		c3(2) = c3(2) - ly
		c3(3) = c3(3) - lz
	endif
	if (m == 23) then
		c3(1) = c3(1) - lx
		c3(2) = c3(2) + ly
		c3(3) = c3(3) + lz
	endif
	if (m == 24) then
		c3(1) = c3(1) - lx
		c3(2) = c3(2) + ly
		c3(3) = c3(3) - lz
	endif
	if (m == 25) then
		c3(1) = c3(1) - lx
		c3(2) = c3(2) - ly
		c3(3) = c3(3) + lz
	endif
	if (m == 26) then
		c3(1) = c3(1) - lx
		c3(2) = c3(2) - ly
		c3(3) = c3(3) - lz
	endif
end subroutine

subroutine cminus(c3r,c3l,dc3,lx,ly,lz)
    real * 8 cdx, cdy, cdz, lx, ly, lz
    real * 8, dimension(3) :: dc3, c3r, c3l
	cdx = c3r(1) - c3l(1)
	cdy = c3r(2) - c3l(2)
	cdz = c3r(3) - c3l(3)
	if (cdx > lx/2) cdx=cdx-lx
	if (cdy > ly/2) cdy=cdy-ly
	if (cdz > lz/2) cdz=cdz-lz
	if (cdx < -lx/2) cdx=cdx+lx
	if (cdy < -ly/2) cdy=cdy+ly
	if (cdz < -lz/2) cdz=cdz+lz
	dc3(1) = cdx
	dc3(2) = cdy
	dc3(3) = cdz
end subroutine

subroutine distflag(c3r,c3l,flag,lx,ly,lz)
    real * 8 cdx, cdy, cdz, lx, ly, lz
    real * 8, dimension(3) :: c3r, c3l, flag
	flag(1) = 0.0d0
	flag(2) = 0.0d0
	flag(3) = 0.0d0
	cdx = c3r(1) - c3l(1)
	cdy = c3r(2) - c3l(2)
	cdz = c3r(3) - c3l(3)
	if (cdx > lx/2) flag(1) = lx
	if (cdy > ly/2) flag(2) = ly
	if (cdz > lz/2) flag(3) = lz
	if (cdx < -lx/2) flag(1) = -lx
	if (cdy < -ly/2) flag(2) = -ly
	if (cdz < -lz/2) flag(3) = -lz
end subroutine

subroutine dsigmadsFunc(m,zeta,length,dsigmads,npoint,cx,cy,cz,cze,ckx,cky,ckz,lx,ly,lz)
    real * 8 zeta, dsigmads, sigmal, sigmar, sigmax, dzeta, ds2, length, lx, ly, lz
	real * 8, dimension(3) :: czeta, czetal, czetar, dczetal, dczetar
	integer m, npoint
	real * 8, dimension(npoint) :: cx, cy, cz, cze, ckx, cky, ckz
	dzeta = 1.0d-4
	call sigmaFunc(zeta-dzeta,length,sigmal,sigmax)
	call sigmaFunc(zeta+dzeta,length,sigmar,sigmax)
	call curve(m,zeta, czeta,npoint,cx,cy,cz,cze,ckx,cky,ckz,lx,ly,lz)
	call curve(m,zeta-dzeta, czetal,npoint,cx,cy,cz,cze,ckx,cky,ckz,lx,ly,lz)
	call curve(m,zeta+dzeta, czetar,npoint,cx,cy,cz,cze,ckx,cky,ckz,lx,ly,lz)
	call cminus(czeta,czetal,dczetal,lx,ly,lz)
	call cminus(czetar,czeta,dczetar,lx,ly,lz)
	ds2=dsqrt(dczetal(1)**2.d0+dczetal(2)**2.d0+dczetal(3)**2.d0)+dsqrt(dczetar(1)**2.d0+dczetar(2)**2.d0+dczetar(3)**2.d0)
	dsigmads=(sigmar-sigmal)/ds2
end subroutine

subroutine d1curve(m,zeta,dc3,npoint,cx,cy,cz,cze,ckx,cky,ckz,lx,ly,lz)
    real * 8 zeta, dzeta, cdx, cdy, cdz, lx, ly, lz
    real * 8, dimension(3) :: dc3, c3r, c3l
	integer m, npoint
	real * 8, dimension(npoint) :: cx, cy, cz, cze, ckx, cky, ckz
    dzeta = 1.0d-4
    call curve(m,zeta + dzeta, c3r,npoint,cx,cy,cz,cze,ckx,cky,ckz,lx,ly,lz)
    call curve(m,zeta - dzeta, c3l,npoint,cx,cy,cz,cze,ckx,cky,ckz,lx,ly,lz)
	call cminus(c3r,c3l,dc3,lx,ly,lz)
    dc3 = dc3 / 2.d0 / dzeta
end subroutine

subroutine d2curve(m,zeta,dc3,npoint,cx,cy,cz,cze,ckx,cky,ckz,lx,ly,lz)
    real * 8 zeta, dzeta, lx, ly, lz
    real * 8, dimension(3) :: dc3, c3r, c3l
	integer m, npoint
	real * 8, dimension(npoint) :: cx, cy, cz, cze, ckx, cky, ckz
    dzeta = 1.0d-4
    call d1curve(m,zeta + dzeta, c3r,npoint,cx,cy,cz,cze,ckx,cky,ckz,lx,ly,lz)
    call d1curve(m,zeta - dzeta, c3l,npoint,cx,cy,cz,cze,ckx,cky,ckz,lx,ly,lz)
    dc3 = (c3r - c3l) / 2.d0 / dzeta    
end subroutine

subroutine d3curve(m,zeta,dc3,npoint,cx,cy,cz,cze,ckx,cky,ckz,lx,ly,lz)
    real * 8 zeta, dzeta, lx, ly, lz
    real * 8, dimension(3) :: dc3, c3r, c3l
	integer m, npoint
	real * 8, dimension(npoint) :: cx, cy, cz, cze, ckx, cky, ckz
    dzeta = 1.0d-4
    call d2curve(m,zeta + dzeta, c3r,npoint,cx,cy,cz,cze,ckx,cky,ckz,lx,ly,lz)
    call d2curve(m,zeta - dzeta, c3l,npoint,cx,cy,cz,cze,ckx,cky,ckz,lx,ly,lz)
    dc3 = (c3r - c3l) / 2.d0 / dzeta    
end subroutine

subroutine norm_zeta(c3,nc3)
    real * 8 nc3
    real * 8, dimension(3) :: c3
    nc3 = dsqrt(c3(1) ** 2 + c3(2) ** 2 + c3(3) ** 2)
end subroutine

subroutine cross_zeta(vec1,vec2,vec3)
    real * 8, dimension(3) :: vec1, vec2, vec3
    vec3(1) = vec1(2) * vec2(3) - vec1(3) * vec2(2)
    vec3(2) = vec1(3) * vec2(1) - vec1(1) * vec2(3)
    vec3(3) = vec1(1) * vec2(2) - vec1(2) * vec2(1)
end subroutine

subroutine dot_zeta(vec1,vec2,sc)
    real * 8, dimension(3) :: vec1, vec2
    real * 8 sc
    sc = vec1(1) * vec2(1) + vec1(2) * vec2(2) + vec1(3) * vec2(3)
end subroutine

subroutine norm_t(n,vec,normv)
    integer n
    real * 8, dimension(n,3) :: vec
    real * 8, dimension(n) :: normv
    normv = dsqrt(vec(:, 1) * vec(:, 1) + vec(:, 2) * vec(:, 2) + vec(:, 3) * vec(:, 3))
end subroutine

subroutine cross_t(n,vec1,vec2,vec3)
    integer n
    real * 8, dimension(n,3) :: vec1, vec2, vec3
    vec3(:, 1) = vec1(:, 2) * vec2(:, 3) - vec1(:, 3) * vec2(:, 2)
    vec3(:, 2) = vec1(:, 3) * vec2(:, 1) - vec1(:, 1) * vec2(:, 3)
    vec3(:, 3) = vec1(:, 1) * vec2(:, 2) - vec1(:, 2) * vec2(:, 1)
end subroutine

subroutine dot_t(n,vec1,vec2,sc)
    integer n
    real * 8, dimension(n,3) :: vec1, vec2
    real * 8, dimension(n) :: sc
    sc = vec1(:, 1) * vec2(:, 1) + vec1(:, 2) * vec2(:, 2) + vec1(:, 3) * vec2(:, 3)
end subroutine

subroutine pt(n,dt,vec,dvec)
    integer n
    real * 8 dt
    real * 8, dimension(n,3) :: vec, dvec, vecr, vecl
    vecr(1 : n - 1, :) = vec(2 : n, :)
    vecr(n, :) = vec(1, :)
    vecl(2 : n, :) = vec(1 : n - 1, :)
    vecl(1, :) = vec(n, :)
    dvec = (vecr - vecl) / 2.d0 / dt
end subroutine

subroutine ptboundary(n,dt,vec,dvec,lx,ly,lz)
    integer n, i
    real * 8 dt, lx, ly, lz
    real * 8, dimension(n,3) :: vec, dvec, vecr, vecl
    vecr(1 : n - 1, :) = vec(2 : n, :)
    vecr(n, :) = vec(1, :)
    vecl(2 : n, :) = vec(1 : n - 1, :)
    vecl(1, :) = vec(n, :)
    !dvec = (vecr - vecl) / 2.d0 / dt
	do i = 1, n
		call cminus(vecr(i, :),vecl(i, :),dvec(i, :),lx,ly,lz)
	enddo
	dvec = dvec / 2.d0 / dt
end subroutine

subroutine set_fft_plans(nx, ny, nz, planxf, planxb, planyf, planyb, planzf, planzb)
    implicit none
    include "fftw3.f"
    integer*8 :: planxf, planyf, planzf, planxb, planyb, planzb
    integer nx, ny, nz
    double complex,dimension(nx) :: tempx
    double complex,dimension(ny) :: tempy
    double complex,dimension(nz) :: tempz
    call dfftw_plan_dft_1d(planxf, nx, tempx, tempx, FFTW_FORWARD, FFTW_ESTIMATE)
    call dfftw_plan_dft_1d(planxb, nx, tempx, tempx, FFTW_BACKWARD, FFTW_ESTIMATE)
    call dfftw_plan_dft_1d(planyf, ny, tempy, tempy, FFTW_FORWARD, FFTW_ESTIMATE)
    call dfftw_plan_dft_1d(planyb, ny, tempy, tempy, FFTW_BACKWARD, FFTW_ESTIMATE)
    call dfftw_plan_dft_1d(planzf, nz, tempz, tempz, FFTW_FORWARD, FFTW_ESTIMATE)
    call dfftw_plan_dft_1d(planzb, nz, tempz, tempz, FFTW_BACKWARD, FFTW_ESTIMATE)
end subroutine set_fft_plans

subroutine destroy_fft_plans(planxf, planxb, planyf, planyb, planzf, planzb)
    implicit none
    include "fftw3.f"
    integer*8 :: planxf, planyf, planzf, planxb, planyb, planzb   
    call dfftw_destroy_plan(planxf)
    call dfftw_destroy_plan(planxb)
    call dfftw_destroy_plan(planyf)
    call dfftw_destroy_plan(planyb)
    call dfftw_destroy_plan(planzf)
    call dfftw_destroy_plan(planzb)
end subroutine destroy_fft_plans

SUBROUTINE wavenumber(nx,ny,nz,nzp,kx,ky,kz,k2,id)
    implicit none
    integer nx, ny, nz, nzp, i, j, k, id
    integer, DIMENSION (nx) :: kx
    integer, DIMENSION (ny) :: ky
    integer, DIMENSION (nzp) :: kz
    integer, dimension (nx,ny,nzp) :: k2
    do i=1,nx
        kx(i)=mod(i-1+nx/2,nx)-nx/2
    end do
    do j=1,ny
        ky(j)=mod(j-1+ny/2,ny)-ny/2
    end do
    do k=1,nzp
        kz(k)=mod(k-1+nz/2+id*nzp,nz)-nz/2
    end do
    do k=1,nzp
        do j=1,ny
            do i=1,nx
                k2(i,j,k)=kx(i)**2+ky(j)**2+kz(k)**2
            enddo
        enddo
    enddo
end subroutine wavenumber

subroutine initialize_mesh(nx, ny, nzp, xstart, ystart, zstart, dx, dy, dz, meshx, meshy, meshz, id)
    implicit none
    integer nx, ny, nzp, i, j, k, id
    real * 8 xstart, ystart ,zstart, dx, dy, dz
    real * 8, dimension(nx,ny,nzp) :: meshx, meshy, meshz 
    do i = 1, nx
        meshx(i, :, :) = (i - 1.d0) * dx + xstart
    enddo
    do j = 1, ny
        meshy(:, j, :) = (j - 1.d0) * dy + ystart
    enddo
    do k = 1, nzp
        meshz(:, :, k) = (k - 1.d0 + nzp * id) * dz + zstart
    enddo
end subroutine initialize_mesh 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine transpose_yz(mx, my, mz, nproc, spec, spectemp, id)
    implicit none
    include 'mpif.h'
    real ierr
    integer mx, my, mz, nproc, i, j, k, imp, mzp, myp, id
    integer, dimension (nproc) :: counts, displs
    double complex, dimension (mx, my, mz / nproc) :: spec
    double complex, dimension (mx, mz, my / nproc) :: spectemp
    double complex, dimension (mx, my / nproc, mz) :: spectemp1
    double complex, dimension (mx, my / nproc, mz / nproc) :: specb
    mzp = mz / nproc
    myp = my / nproc
    do i = 1, nproc
        counts(i) = mx * myp * mzp
        displs(i) = (i - 1) * counts(i)
    enddo
    do imp = 1, nproc
        do k = 1, mzp
            do j = 1, myp
                specb(:, j, k) = spec(:, j + (imp - 1) * myp, k)
            enddo
        enddo
        call mpi_gatherv(specb, counts(id + 1), MPI_DOUBLE_COMPLEX, spectemp1,& 
            &counts,displs, MPI_DOUBLE_COMPLEX,imp-1,mpi_comm_world,ierr)
    enddo
    do j = 1, myp
        do k = 1, mz
            spectemp(:, k, j) = spectemp1(:, j, k)
        enddo
    enddo
end subroutine transpose_yz
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! 3D forward fast fourier transformation
!! mpi for z direction, 
!!    number of blocks: nproc=nz/nzp 
!!    id = 0, 1, ... , nproc-1
!! input:
!!    phy: real, dimension (nx,ny,nzp)
!!    mesh size: nx, ny, nz
!!	  fft plans: planxf, planyf, planzf
!! output:
!!	  spec: complex,dimension (nx,ny,nzp)	
subroutine fourier_forward(phy,spec,nx,ny,nz,planxf,planyf,planzf,id,nproc)
    implicit none
    include 'mpif.h'
    include "fftw3.f"
    integer*8 planxf, planyf, planzf
    !!!MPI
    integer :: id, nproc, ierr
    integer nx, ny, nz, nzp, nyp, i, j, k
    real * 8, dimension (nx,ny,nz/nproc) :: phy
    double complex, dimension (nx,ny,nz/nproc) :: spec
    double complex, dimension (nx) :: tempx
    double complex, dimension (ny) :: tempy
    double complex, dimension (nz) :: tempz
    double complex, dimension (nx,nz,ny/nproc) :: spectemp
    nyp=ny/nproc
    nzp=nz/nproc
    do k=1,nzp
        do j=1,ny
        	do i=1,nx
                tempx(i)=dcmplx(phy(i,j,k)/nx+0.0d0,0.0d0)
            enddo
            call dfftw_execute_dft(planxf,tempx,tempx)
            spec(:,j,k)=tempx
        enddo
    enddo
    do k=1,nzp
        do i=1,nx
            tempy=spec(i,:,k)/ny
            call dfftw_execute_dft(planyf,tempy,tempy)
            spec(i,:,k)=tempy
        enddo
    enddo
    call transpose_yz(nx,ny,nz,nproc,spec,spectemp,id)
    do j=1,nyp
        do i=1,nx
            tempz=spectemp(i,:,j)/nz
            call dfftw_execute_dft(planzf,tempz,tempz)
            spectemp(i,:,j)=tempz
        enddo
    enddo
    call transpose_yz(nx,nz,ny,nproc,spectemp,spec,id)
end subroutine fourier_forward


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! 3D backward fast fourier transformation
!! mpi for z direction, 
!!    number of blocks: nproc=nz/nzp 
!!    id = 0, 1, ... , nproc-1
!! input:
!!    phy: real, dimension (nx,ny,nzp)
!!    mesh size: nx, ny, nz
!!	  fft plans: planxf, planyf, planzf
!! output:
!!	  spec: complex,dimension (nx,ny,nzp)	
subroutine fourier_backward(phy,spec,nx,ny,nz,planxb,planyb,planzb,id,nproc)
    implicit none
    include 'mpif.h'
    include "fftw3.f"
    integer*8 planxb, planyb, planzb
    integer :: id, nproc, ierr
    integer nx, ny, nz, nzp, nyp, i, j, k
    real * 8, dimension (nx,ny,nz/nproc) :: phy
    double complex, dimension (nx,ny,nz/nproc) :: spec
    double complex, dimension (nx) :: tempx
    double complex, dimension (ny) :: tempy
    double complex, dimension (nz) :: tempz
    double complex, dimension (nx,nz,ny/nproc) :: spectemp
    double complex, allocatable :: specall(:,:,:)
    double complex, allocatable :: specalltemp(:,:,:)
    if(id==0) then
        allocate(specall(nx,ny,nz))
        allocate(specalltemp(nx,nz,ny))
    endif
    nyp=ny/nproc
    nzp=nz/nproc
    do k=1,nzp
        do j=1,ny
            tempx=spec(:,j,k)
            call dfftw_execute_dft(planxb,tempx,tempx)
            spec(:,j,k)=tempx
        enddo
    enddo
    do k=1,nzp
        do i=1,nx
            tempy=spec(i,:,k)
            call dfftw_execute_dft(planyb,tempy,tempy)
            spec(i,:,k)=tempy
        enddo
    enddo
    call transpose_yz(nx,ny,nz,nproc,spec,spectemp,id)
    do j=1,nyp
        do i=1,nx
            tempz=spectemp(i,:,j)
            call dfftw_execute_dft(planzb,tempz,tempz)
            spectemp(i,:,j)=tempz
        enddo
    enddo
    call transpose_yz(nx,nz,ny,nproc,spectemp,spec,id)
    phy=dreal(spec)
end subroutine fourier_backward

subroutine dx_dy_dz_dp_dm(phy,dphy,switch_d,nx,ny,nzp,kx,ky,kz,k2,planxf,planyf,planzf,&
    &planxb,planyb,planzb,id,nproc)
    implicit none
    integer*8 planxf,planyf,planzf,planxb,planyb,planzb
    integer nx, ny, nz, nzp, i, j, k
    integer, DIMENSION (nx) :: kx
    integer, DIMENSION (ny) :: ky
    integer, DIMENSION (nzp) :: kz
    integer, dimension (nx,ny,nzp) :: k2
    integer :: id, nproc, switch_d
    real * 8, dimension (nx,ny,nzp) :: phy, dphy
    double complex, dimension (nx,ny,nzp) :: spec
    nz=nzp*nproc
    call fourier_forward(phy,spec,nx,ny,nz,planxf,planyf,planzf,id,nproc)
    if(switch_d==1) then
        do i=1,nx
            spec(i,:,:)=spec(i,:,:)*dcmplx(0.0d0,kx(i)+0.0d0)
        enddo
    elseif(switch_d==2) then
        do j=1,ny
            spec(:,j,:)=spec(:,j,:)*dcmplx(0.0d0,ky(j)+0.0d0)
        enddo
    elseif(switch_d==3) then
        do k=1,nzp
            spec(:,:,k)=spec(:,:,k)*dcmplx(0.0d0,kz(k)+0.0d0)
        enddo
    elseif(switch_d==6) then
        spec=-k2*spec
    elseif(switch_d==-6) then
        spec=-spec/k2
        if(id==0) then
            spec(1,1,1)=dcmplx(0.0d0,0.0d0)
        endif
    endif
    call fourier_backward(dphy,spec,nx,ny,nz,planxb,planyb,planzb,id,nproc)
end subroutine dx_dy_dz_dp_dm

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!if (switch_d==-1) (dphy1,dphy2,dphy3) is divergence free
subroutine cross_vector(phy1,phy2,phy3,dphy1,dphy2,dphy3,switch_d,nx,ny,nzp,kx,ky,kz,k2,planxf,planyf,planzf,&
    &planxb,planyb,planzb,id,nproc)
    implicit none
    integer*8 planxf,planyf,planzf,planxb,planyb,planzb
    integer nx, ny, nz, nzp, i, j, k
    integer, DIMENSION (nx) :: kx
    integer, DIMENSION (ny) :: ky
    integer, DIMENSION (nzp) :: kz
    integer, dimension (nx,ny,nzp) :: k2
    integer :: id, nproc, switch_d
    real * 8, dimension (nx,ny,nzp) :: phy1, phy2, phy3, dphy1, dphy2, dphy3
    double complex, dimension (nx,ny,nzp) :: spec1, spec2, spec3, spec4, spec5, spec6
    nz=nzp*nproc
    call fourier_forward(phy1,spec1,nx,ny,nz,planxf,planyf,planzf,id,nproc)
    call fourier_forward(phy2,spec2,nx,ny,nz,planxf,planyf,planzf,id,nproc)
    call fourier_forward(phy3,spec3,nx,ny,nz,planxf,planyf,planzf,id,nproc)
    if(switch_d==1) then
    	do k=1,nzp
    		do j=1,ny
    			do i=1,nx
                    spec4(i,j,k)=dcmplx(0.0d0,1.0d0)*(ky(j)*spec3(i,j,k)-kz(k)*spec2(i,j,k))
                    spec5(i,j,k)=dcmplx(0.0d0,1.0d0)*(kz(k)*spec1(i,j,k)-kx(i)*spec3(i,j,k))
                    spec6(i,j,k)=dcmplx(0.0d0,1.0d0)*(kx(i)*spec2(i,j,k)-ky(j)*spec1(i,j,k))
    			enddo
    		enddo
    	enddo
    elseif(switch_d==-1) then
    	do k=1,nzp
    		do j=1,ny
    			do i=1,nx
                    spec4(i,j,k)=dcmplx(0.0d0,1.0d0)*(ky(j)*spec3(i,j,k)-kz(k)*spec2(i,j,k))/k2(i,j,k)
                    spec5(i,j,k)=dcmplx(0.0d0,1.0d0)*(kz(k)*spec1(i,j,k)-kx(i)*spec3(i,j,k))/k2(i,j,k)
                    spec6(i,j,k)=dcmplx(0.0d0,1.0d0)*(kx(i)*spec2(i,j,k)-ky(j)*spec1(i,j,k))/k2(i,j,k)
    			enddo
    		enddo
    	enddo
    	if(id==0) then
            spec4(1,1,1)=dcmplx(0.0d0,0.0d0)
            spec5(1,1,1)=dcmplx(0.0d0,0.0d0)
            spec6(1,1,1)=dcmplx(0.0d0,0.0d0)
        endif
    endif
    do k=1,nzp
        do j=1,ny
            do i=1,nx
                if(k2(i,j,k) > (nx**2 + ny**2 + nz**2) / 27) then
                    spec4(i,j,k) = dcmplx(0.0d0,0.0d0)
                    spec5(i,j,k) = dcmplx(0.0d0,0.0d0)
                    spec6(i,j,k) = dcmplx(0.0d0,0.0d0)
                endif
            enddo
        enddo
    enddo
    call fourier_backward(dphy1,spec4,nx,ny,nz,planxb,planyb,planzb,id,nproc)
    call fourier_backward(dphy2,spec5,nx,ny,nz,planxb,planyb,planzb,id,nproc)
    call fourier_backward(dphy3,spec6,nx,ny,nz,planxb,planyb,planzb,id,nproc)
end subroutine cross_vector

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine divergence(phy1,phy2,phy3,dphy,nx,ny,nzp,kx,ky,kz,planxf,planyf,planzf,&
    &planxb,planyb,planzb,id,nproc)
    implicit none
    integer*8 planxf,planyf,planzf,planxb,planyb,planzb
    integer nx, ny, nz, nzp, i, j, k
    integer, DIMENSION (nx) :: kx
    integer, DIMENSION (ny) :: ky
    integer, DIMENSION (nzp) :: kz
    integer :: id, nproc
    real * 8, dimension (nx,ny,nzp) :: phy1, phy2, phy3, dphy
    double complex, dimension (nx,ny,nzp) :: spec1, spec2, spec3
    nz=nzp*nproc
    call fourier_forward(phy1,spec1,nx,ny,nz,planxf,planyf,planzf,id,nproc)
    call fourier_forward(phy2,spec2,nx,ny,nz,planxf,planyf,planzf,id,nproc)
    call fourier_forward(phy3,spec3,nx,ny,nz,planxf,planyf,planzf,id,nproc)
    do k=1,nzp
    	do j=1,ny
    		do i=1,nx
                spec1(i,j,k)=dcmplx(0.0d0,1.0d0)*(kx(i)*spec1(i,j,k)+ky(j)*spec2(i,j,k)+kz(k)*spec3(i,j,k))
    		enddo
    	enddo
    enddo
    call fourier_backward(dphy,spec1,nx,ny,nz,planxb,planyb,planzb,id,nproc)
end subroutine divergence

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine gradient(phy,dphy1,dphy2,dphy3,nx,ny,nzp,kx,ky,kz,planxf,planyf,planzf,&
    &planxb,planyb,planzb,id,nproc)
    implicit none
    integer nx, ny, nz, nzp, i, j, k
    integer*8 planxf,planyf,planzf,planxb,planyb,planzb
    integer, DIMENSION (nx) :: kx
    integer, DIMENSION (ny) :: ky
    integer, DIMENSION (nzp) :: kz
    integer :: id, nproc
    real * 8, dimension (nx,ny,nzp) :: dphy1, dphy2, dphy3, phy
    double complex, dimension (nx,ny,nzp) :: spec, spec1, spec2, spec3
    nz=nzp*nproc
    call fourier_forward(phy,spec,nx,ny,nz,planxf,planyf,planzf,id,nproc)
    do k=1,nzp
    	do j=1,ny
    		do i=1,nx
                spec1(i,j,k)=dcmplx(0.0d0,1.0d0)*kx(i)*spec(i,j,k)
                spec2(i,j,k)=dcmplx(0.0d0,1.0d0)*ky(j)*spec(i,j,k)
                spec3(i,j,k)=dcmplx(0.0d0,1.0d0)*kz(k)*spec(i,j,k)
    		enddo
    	enddo
    enddo
    call fourier_backward(dphy1,spec1,nx,ny,nz,planxb,planyb,planzb,id,nproc)
    call fourier_backward(dphy2,spec2,nx,ny,nz,planxb,planyb,planzb,id,nproc)
    call fourier_backward(dphy3,spec3,nx,ny,nz,planxb,planyb,planzb,id,nproc)
end subroutine gradient

subroutine cross_product(nx,ny,nzp,phy1,phy2,phy3,phy11,phy12,phy13,phy21,phy22,phy23)
    implicit none
    integer nx, ny, nzp
    real * 8, dimension (nx,ny,nzp) :: phy1, phy2, phy3, phy11, phy12, phy13, phy21, phy22, phy23
    phy21=phy2*phy13-phy3*phy12
    phy22=phy3*phy11-phy1*phy13
    phy23=phy1*phy12-phy2*phy11
end subroutine cross_product

subroutine output_v(name,varname,nx,ny,nzp,data_box,nbox,id,nproc)
    implicit none
    include 'mpif.h'
    character*200 name
    integer nx, ny, nzp, nz, id, nproc, ierr, i, nbox
    integer, dimension (nproc) :: counts,displs
    real, dimension(nx,ny,nzp,nbox) :: data_box
    real, dimension(nx,ny,nzp) :: temp
    real, allocatable :: v_box(:,:,:,:)
    character*40, dimension (nbox) :: varname
    if(id==0) then
        nz=nzp*nproc
        allocate(v_box(nx,ny,nz,nbox))
    endif
    do i=1,nproc
        counts(i)=nx*ny*nzp
        displs(i)=(i-1)*counts(i)
    enddo
    do i=1,nbox
        temp=data_box(:,:,:,i)
        call mpi_gatherv(temp,counts(id+1),MPI_real,v_box(:,:,:,i),counts,displs,&
            &MPI_real,0,mpi_comm_world,ierr)
    enddo
    if(id==0) call output_3d_tecplot_bin(v_box,nx,ny,nz,nbox,name,varname)
    if(id==0) deallocate(v_box)
end subroutine output_v

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! output binary data for tecplot
!! input:
!!    v_box: real, dimension (nx,ny,nz,nbox)
!!    nbox: number of output varible
!!    mesh size: nx, ny, nz
!!	  name: character*200 name of the output data
!!    varname: character*40, dimension (nbox) names of varible
subroutine output_3d_tecplot_bin(v_box,nx,ny,nz,nbox,name,varname)
    implicit none
    integer nx, ny, nz, nbox
    real, dimension (nx,ny,nz,nbox) :: v_box
    real*8, dimension (nbox) :: min_value, max_value
    real*4 ZONEMARKER,EOHMARKER
    integer len,i
    character*40 Title,var
    character*200 name
    character*40, dimension (nbox) :: varname
    character*40 Zonename
    character(40) instring
    ZONEMARKER= 299.0
    EOHMARKER = 357.0
    do i=1,nbox
        min_value(i)=minval(v_box(:,:,:,i))
        max_value(i)=maxval(v_box(:,:,:,i))
    enddo
    open(unit=99,file=name,form="BINARY")
    !I. The header section.
    !1.1 Magic number, Version number
    write(99) "#!TDV112"
    !1.2. Integer value of 1.
    write(99) 1
    !1.3. Title and variable names.
    !Filetype
    write(99) 0
    !1.3.1. The TITLE.
    Title=""
    call dumpstring(Title)
    !1.3.2 Number of variables (NumVar) in the datafile.
    write(99) nbox
    !1.3.3 Variable names. N = L[1] + L[2] + .... L[NumVar]
    do i=1,nbox
        call dumpstring(varname(i))
    enddo
    !1.4. Zones
    !Zone marker. Value = 299.0
    write(99) ZONEMARKER
    !Zone name.
    Zonename="ZONE 001"
    call dumpstring(Zonename)
    !ParentZone
    write(99) -1
    !StrandID
    write(99) -1
    !solution time
    write(99) 0
    write(99) 0
    !not used
    write(99) -1
    !ZoneType  
    write(99) 0
    !DataPacking 0=Block, 1=Point
    write(99) 0
    !Specify Var Location. 0 = Don't specify, all data is located at the nodes. 1 = Specify
    write(99) 0
    !Number of user defined face neighbor connections (value >= 0)
    write(99) 0
    !IMax,JMax,KMax
    write(99) nx
    write(99) ny
    write(99) nz
    !1=Auxiliary name/value pair to follow   0=No more Auxiliar name/value pairs.
    write(99) 0
    !I HEADER OVER
    !EOHMARKER, value=357.0
    write(99) EOHMARKER
    !II. Data section
    !2.1 zone
    write(99) Zonemarker
    !variable data format, 1=Float, 2=Double, 3=LongInt, 4=ShortInt, 5=Byte, 6=Bit                
    do i=1,nbox
        write(99) 1                                                               
    enddo
    !Has variable sharing 0 = no, 1 = yes.       
    write(99) 0
    !Has passive variables 0 = no, 1 = yes.       
    write(99) 0                                  
    !Zone number to share connectivity list with (-1 = no sharing).
    write(99) -1
    !min value
    !max value
    do i=1,nbox
        write(99) min_value(i)
        write(99) max_value(i)                                                               
    enddo
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Zone Data. Each variable is in data format asspecified above.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    write(99) v_box
    close(99)
end subroutine output_3d_tecplot_bin
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine dumpstring(instring)
    !!!for binary output
    character(40) instring
    integer len
    len=LEN_TRIM(instring)
    do i=1,len
        ii=ICHAR(instring(i:i))
        write(99) ii
    enddo
    write(99) 0
end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine cal_deviation(vorx,vory,vorz,phiv,nx,ny,nz,nzp,deviation,dx,dy,dz,id,igst)
    implicit none
    include 'mpif.h'
    integer i, j, k
    real * 8 dx, dy, dz
    integer nx, ny, nzp, nz, igst, ierr, id
    real*8 deviation, deviation1
    real * 8 eps
    real * 8 dxphi,dyphi,dzphi,ome_dot_graph,graph,vorti
    real * 8, dimension (nx,ny,nzp) :: vorx,vory,vorz
    real * 8 phiv(1-igst:nx+igst,1-igst:ny+igst,1-igst:nzp+igst)
    eps=0.1d-10
    deviation1=0.d0
    do k=1,nzp
        do j=1,ny
            do i=1,nx
                dxphi=(9.0d0*phiv(i-2,j,k)-phiv(i-3,j,k)-45.0d0*phiv(i-1,j,k)+45.0d0*phiv(i+1,j,k)-9.0d0*phiv(i+2,j,k)+phiv(i+3,j,k))/60.0d0/dx
                dyphi=(9.0d0*phiv(i,j-2,k)-phiv(i,j-3,k)-45.0d0*phiv(i,j-1,k)+45.0d0*phiv(i,j+1,k)-9.0d0*phiv(i,j+2,k)+phiv(i,j+3,k))/60.0d0/dy
                dzphi=(9.0d0*phiv(i,j,k-2)-phiv(i,j,k-3)-45.0d0*phiv(i,j,k-1)+45.0d0*phiv(i,j,k+1)-9.0d0*phiv(i,j,k+2)+phiv(i,j,k+3))/60.0d0/dz
                graph=dsqrt(dxphi**2.d0+dyphi**2.d0+dzphi**2.d0)
                vorti=dsqrt((vorx(i,j,k))**2.d0+(vory(i,j,k))**2.d0+(vorz(i,j,k))**2.d0)
                ome_dot_graph=vorx(i,j,k)*dxphi+vory(i,j,k)*dyphi+vorz(i,j,k)*dzphi
                deviation1=deviation1+dabs(ome_dot_graph)/(vorti*graph+eps)/nx/ny/nz
            enddo
        enddo
    enddo   
    call MPI_REDUCE(deviation1,deviation,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
end subroutine cal_deviation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine cal_helicity(vorx,vory,vorz,velx,vely,velz,nx,ny,nz,nzp,helicity,dv,id)
    implicit none
    include 'mpif.h'
    integer i, j, k
    real*8 dv
    integer nx, ny, nzp, nz, ierr, id
    real*8 helicity, helicity1
    real*8, dimension (nx,ny,nzp) :: vorx, vory, vorz, velx, vely, velz
    helicity1=0.d0
    do k=1,nzp
        do j=1,ny
            do i=1,nx
                helicity1=helicity1+vorx(i, j, k) * velx(i, j, k) + &
                &vory(i, j, k) * vely(i, j, k) + vorz(i, j, k) * velz(i, j, k)
            enddo
        enddo
    enddo
    call MPI_REDUCE(helicity1,helicity,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
    helicity = helicity * dv
end subroutine cal_helicity
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Set the value on the ghost points for phiv
subroutine wall_phiv(nx,ny,nzp,phiv,igst,id_l,id_r)
    implicit none
    include 'mpif.h'
    integer ierr, i, j
    integer nx, ny, nzp, igst, nxy, id_l, id_r
    integer status(MPI_STATUS_SIZE)
    real*8 phiv(1-igst:nx+igst,1-igst:ny+igst,1-igst:nzp+igst)
    if(igst/=0) then
        nxy=(nx+2*igst)*(ny+2*igst)*igst
        do i=1,igst
            phiv(nx+i,:,:)=phiv(i,:,:)
            phiv(-igst+i,:,:)=phiv(nx-igst+i,:,:)
        enddo
        do j=1,igst
            phiv(:,ny+j,:)=phiv(:,j,:)
            phiv(:,-igst+j,:)=phiv(:,ny-igst+j,:)
        enddo
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!z-direction, need communication between nodes
        !!!Communication with each other as on a ring. 
        call mpi_sendrecv(phiv(:,:,1:igst), nxy, MPI_DOUBLE_PRECISION,id_l,99,&
                    phiv(:,:,nzp+1:nzp+igst), nxy, MPI_DOUBLE_PRECISION,id_r,99,mpi_comm_world, status, ierr)
        call mpi_sendrecv(phiv(:,:,nzp-igst+1:nzp), nxy, MPI_DOUBLE_PRECISION,id_r,99,&
                    phiv(:,:,1-igst:0), nxy, MPI_DOUBLE_PRECISION,id_l,99,mpi_comm_world, status, ierr)
        call mpi_barrier(mpi_comm_world,ierr)
    endif
end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine spectrum(nx,ny,nzp,velx,vely,velz,num_data,k2,id,nproc,planxf,planyf,planzf)
    implicit none
    integer*8 planxf, planyf, planzf
    integer nx, ny, nzp, nz, num_data
    real * 8, dimension (nx,ny,nzp) :: velx, vely, velz
    integer id, nproc
    integer, dimension (nx,ny,nzp) :: k2
    double complex, dimension (nx,ny,nzp) :: specx, specy, specz
    real, dimension (nx,ny,nzp) :: spec
    nz = nzp * nproc
    call fourier_forward(velx,specx,nx,ny,nz,planxf,planyf,planzf,id,nproc)
    call fourier_forward(vely,specy,nx,ny,nz,planxf,planyf,planzf,id,nproc)
    call fourier_forward(velz,specz,nx,ny,nz,planxf,planyf,planzf,id,nproc)
    spec = real(specx*conjg(specx) + specy*conjg(specy) + specz*conjg(specz)) / 2.0
    call get_spectrum(nx,ny,nzp,nz,spec,num_data,k2,id)
end subroutine spectrum


subroutine get_spectrum(nx,ny,nzp,nz,spec,num_data,k2,id)
    implicit none
    include 'mpif.h'
    real ierr
    integer nx, ny, nz, nzp, nproc, num_data, nek, i, id
    real, dimension (nx,ny,nzp) :: spec
    integer, dimension (nx,ny,nzp) :: k2, ik2
    real spectrum_k1, spectrum_k
    nek = int(sqrt((nx**2+ny**2+nz**2)/27.0))
    ik2 = int(sqrt(k2+0.00000001)+0.5)
    do i=1,nek
        spectrum_k1=sum(spec,mask=(ik2.eq.i))
        call MPI_REDUCE(spectrum_k1,spectrum_k,1,MPI_real,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        if(id==0) write(num_data,*) i, spectrum_k
    end do
end subroutine get_spectrum