PROGRAM Model
REAL, dimension(50,50,50) :: eps,Power
REAL, dimension(50,50,50,2) :: flux,oldflux,guessflux
REAL, dimension(2,2) :: sigs
REAL, dimension(2) :: rnubar,b
REAL, dimension(2) :: siga,rnusigf,sigt,chi,D,rlsquared
REAL :: Abar,rmubar
REAL :: epsk,A,tolerin,tolerout,totalnew,totalold,T,dex2,deltax,dey2,deltay,dez2,deltaz,aL,aB,aR,aT,aU,aD,x,y,z
REAL :: rkeff,rkold,Ef,Powermult,totpower,counter,PApower,analyt
CHARACTER :: dummy
INTEGER :: itflag,outflag,itergs,itouter,M,o,i,j,k
M = 50
itouter = 0
itergs = 0
! Read in values from the input file
OPEN(1,FILE='input.dta')
OPEN(2,FILE='output.dta')
READ(1,'(A)') dummy
READ(1,'(A)') dummy
READ(1,'(A)') dummy
READ(1,*) siga(1)
READ(1,'(A)') dummy
READ(1,*) rnusigf(1)
READ(1,'(A)') dummy
READ(1,*) sigt(1)
READ(1,'(A)') dummy
READ(1,*) sigs(1,1)
READ(1,'(A)') dummy
READ(1,*) chi(1)
READ(1,'(A)') dummy
READ(1,'(A)') dummy
READ(1,*) siga(2)
READ(1,'(A)') dummy
READ(1,*) rnusigf(2)
READ(1,'(A)') dummy
READ(1,*) sigt(2)
READ(1,'(A)') dummy
READ(1,*) sigs(2,2)
READ(1,'(A)') dummy
READ(1,*) sigs(1,2)
READ(1,'(A)') dummy
READ(1,*) chi(2)
READ(1,'(A)') dummy
READ(1,*) Abar
READ(1,'(A)') dummy
READ(1,*) aL
READ(1,'(A)') dummy
READ(1,*) aR
READ(1,'(A)') dummy
READ(1,*) aT
READ(1,'(A)') dummy
READ(1,*) aB
READ(1,'(A)') dummy
READ(1,*) aU
READ(1,'(A)') dummy
READ(1,*) aD
READ(1,'(A)') dummy
READ(1,*) x
READ(1,'(A)') dummy
READ(1,*) y
READ(1,'(A)') dummy
READ(1,*) z
READ(1,'(A)') dummy
READ(1,*) rnubar(1)
READ(1,'(A)') dummy
READ(1,*) rnubar(2)
rmubar = 2.0/(3.0*Abar)
! Energy per fission in joules
Ef = 195.6*1.602E-13
!write(*,*) ''
!write(*,*) 'mubar of fuel =', rmubar(1)
!write(*,*) ''
!write(*,*) 'mubar of water =', rmubar(2)
!write(*,*) ''
!write(*,*) 'mubar of CR =', rmubar(3)
deltax = x/M
deltay = y/M
deltaz = z/M
dex2 = deltax*deltax
dey2 = deltay*deltay
dez2 = deltaz*deltaz
D(1) = 1.0/(3.0*(sigt(1)-(rmubar*(sigs(1,1)+sigs(1,2)))))
D(2) = 1.0/(3.0*(sigt(2)-(rmubar*sigs(2,2))))
rlsquared(1) = D(1)/(sigt(1)-sigs(1,1))
rlsquared(2) = D(2)/(sigt(2)-sigs(2,2))
!DO i = 1,2
!	DO j = 1,3
!		write(*,*) ''
!		write(*,*) 'D group', i, 'material', j, '=', D(i,j)
!		write(*,*) ''
!		write(*,*) 'L squared group', i, 'material', j, '=', rlsquared(i,j)
!	END DO
!END DO
!DO i = 1,2
!	DO j = 1,3
!		write(*,*) ''
!		write(*,*) 'nusigmaf of group', i, 'material', j, '=', rnusigf(i,j)
!		write(*,*) ''
!		write(*,*) 'chi of group 1', i, 'material', j, '=', chi(i,j)
!	END DO
!END DO
write(*,*) '3D Bare Takeda Reactor Criticality Search'
write(2,*) '3D Bare Takeda Reactor Criticality Search'
write(*,*) ''
! Initializing epsilon and coefficient values
A = 0.0
b = 0.0
eps = 0.0
tolerin = 0.0001
tolerout = 0.000001
! Guess a keffective and a flux profile
rkeff = 1.0
flux = 1000.0
outflag = 1
DO WHILE (outflag == 1)
	itouter = itouter + 1
	! Set new and old fluxes and ks equal to one another
	rkold = rkeff
	guessflux = flux
	! Start inner iterations
	! Start the energy group loop to converge each flux value.
	DO o = 1,2
		! Reset itflag to one for the next group
		itflag = 1
		DO WHILE (itflag == 1)
			itergs = itergs + 1
			oldflux = flux
			DO k = 1,M
				DO j = 1,M
					DO i = 1,M
						T = 0.0
						! compute the source term at the node of interest, depending on energy group.
						IF (o == 1) THEN
							b(o) = -(chi(o)/(D(o)*rkeff))*(rnusigf(1)*guessflux(i,j,k,1) + rnusigf(2)*guessflux(i,j,k,2))
						ELSE
							b(o) = -(chi(o)/(D(o)*rkeff))*(rnusigf(1)*guessflux(i,j,k,1)+rnusigf(2)*guessflux(i,j,k,2))-(flux(i,j,k,1)*sigs(1,2)/D(o))
						ENDIF
						! Start constructing equation for flux coefficient at the node of interest.
							A = -(2.0/dex2)-(2.0/dey2)-(2.0/dez2)-(1/rlsquared(o))
						! Albedo factors are Left/Right for x, Bottom/Top for y, Down/Up for z.
						! Construct x dimension portion of the coefficients.
						IF (i == 1) THEN
							A = A-((1.0-aL)/(D(o)*deltax*(1.0+aL)))
							T = T+((2.0/dex2)*flux(i+1,j,k,o))
						ELSEIF (i == M) THEN
							A = A-((1.0-aR)/(D(o)*deltax*(1.0+aR)))
							T = T+((2.0/dex2)*flux(i-1,j,k,o))
						ELSE
							T = T+((flux(i+1,j,k,o)+flux(i-1,j,k,o))/dex2)
						ENDIF
						! Construct y dimension coefficients.
						IF (j == 1) THEN
							A = A-((1.0-aB)/(D(o)*deltay*(1.0+aB)))
							T = T+((2.0/dey2)*flux(i,j+1,k,o))
						ELSEIF (j == M) THEN
							A = A-((1.0-aT)/(D(o)*deltay*(1.0+aT)))
							T = T+((2.0/dey2)*flux(i,j-1,k,o))
						ELSE
							T = T+((flux(i,j+1,k,o)+flux(i,j-1,k,o))/dey2)
						ENDIF
						! Construct z dimension coefficients.
						IF (k == 1) THEN
							A = A-((1.0-aD)/(D(o)*deltaz*(1.0+aD)))
							T = T+((2.0/dez2)*flux(i,j,k+1,o))
						ELSEIF (k == M) THEN
							A = A-((1.0-aU)/(D(o)*deltaz*(1.0+aU)))
							T = T+((2.0/dez2)*flux(i,j,k-1,o))
						ELSE
							T = T+((flux(i,j,k+1,o)+flux(i,j,k-1,o))/dez2)
						ENDIF
						!Apply Gauss Seidel forward substitution
						flux(i,j,k,o) = (b(o)-T)/A
						! create corresponding epsilon value
						IF (oldflux(i,j,k,o) > 0) THEN
							eps(i,j,k) = ABS(flux(i,j,k,o)-oldflux(i,j,k,o))/oldflux(i,j,k,o)
						ELSE
							eps(i,j,k) = ABS(flux(i,j,k,o)-oldflux(i,j,k,o))
						ENDIF
					END DO
				END DO
			END DO
			itflag = 0
			IF (maxval(eps) > tolerin) THEN
				itflag = 1
			END IF
			!write(*,*) 'b is', b(60,60,1)
			!write(*,*) 'flux is', flux(30,30,2)
		END DO
	END DO
	totalnew = 0.0
	totalold = 0.0
	! Compute source term
	DO i = 1,M
		DO j = 1,M
			DO k = 1,M
				totalnew = totalnew + ((rnusigf(1)*flux(i,j,k,1))+(rnusigf(2)*flux(i,j,k,2)))
				totalold = totalold + ((rnusigf(1)*guessflux(i,j,k,1))+(rnusigf(2)*guessflux(i,j,k,2)))
			END DO
		END DO
	END DO
	! Use source term computation to produce a new k.
	rkeff = rkold*(totalnew/totalold)
	outflag = 0
	! Check for K convergence
	IF (rkold > 0) THEN 
		epsk = ABS((rkeff-rkold))/rkold
	ELSE
		epsk = ABS(rkeff-rkold)
	ENDIF
	IF (epsk > tolerout) THEN
		outflag = 1
	ENDIF
END DO
Power(i,j,k) = 0.0
totpower = 0.0
DO i = 1,M
	DO j = 1,M
		DO k = 1,M
			Power(i,j,k) = ((rnusigf(1)*flux(i,j,k,1)/rnubar(1)) + (rnusigf(2)*flux(i,j,k,2)/rnubar(2)))*Ef*deltax*deltay*deltaz
			totpower = totpower + Power(i,j,k)
		END DO
	END DO
END DO
Powermult = 125000.0/totpower
Power = Powermult*Power
totpower = 125000.0
flux = Powermult*flux
! Count the amount of fuel nodes
counter = 0.0
DO i = 1,M
	DO j = 1,M
		DO k = 1,M
			counter = counter + 1.0
		END DO
	END DO
END DO
write (*,*) 'number of fuel nodes =', counter
PApower = maxval(Power)/(totpower/counter)
!DO i = 1,M
!	write(2,*) flux(i,:,1)
!END DO
!write(2,*) ''
!DO i = 1,M
!	write(2,*) flux(i,:,2)
!END DO
analyt = rnusigf(2)*sigs(1,2)
analyt = analyt/siga(2)
analyt = rnusigf(1)-(siga(1)+sigs(1,2))+analyt
analyt = (3.0*D(1)*(3.14159**2.0))/analyt
analyt = SQRT(analyt) - (2.0*D(1))
WRITE(*,*) 'analytical side length is', analyt
WRITE(*,*) ''
WRITE(*,*) 'keff is', rkeff
WRITE(*,*) ''
WRITE(*,*) 'Total number of inner iterations =', itergs
WRITE(*,*) ''
WRITE(*,*) 'Total number of outer iterations =', itouter
WRITE(*,*) ''
WRITE(*,*) 'Peak flux with total 1/8th reactor power of 125 KW =', maxval(flux), 'n/cm^2'
WRITE(*,*) ''
WRITE(*,*) 'Peak to Average Power =', PApower
WRITE(2,*) ''
WRITE(2,*) 'keff is', rkeff
WRITE(2,*) ''
WRITE(2,*) 'Total number of inner iterations =', itergs
WRITE(2,*) ''
WRITE(2,*) 'Total number of outer iterations =', itouter
WRITE(2,*) ''
WRITE(2,*) 'Peak flux with total 1/8th reactor power of 125 KW =', maxval(flux), 'n/cm^2'
WRITE(2,*) ''
WRITE(2,*) 'Peak to Average Power =', PApower
close(1)
close(2)
END PROGRAM Model
