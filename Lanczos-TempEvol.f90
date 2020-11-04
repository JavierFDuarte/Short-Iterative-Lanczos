Program lanczos
 implicit none
 integer, parameter                         :: dp = kind(1.d0)
 integer                                    :: N,i,k
 complex(dp),dimension(:,:),allocatable     :: H,D,Z
 real(dp),dimension(:),allocatable          :: eigenval
 real(dp)                                   :: s,pi=2.d0*asin(1.d0),re,im
 complex(dp),dimension(:),allocatable       :: nu, Psi_s
 real(dp),dimension(:),allocatable          :: realalpha,imagalpha,beta
 !Variables para dsyev
 integer                                    :: information, lworkk
 complex(dp), dimension(:), allocatable     :: workk
 real(dp),dimension(:),allocatable          :: rworkk


 !Dimension Espacio Krylov, N
 N=4

!Delta de Tiempo de la Evolución
 s = 0.8

 allocate(H(1:N,1:N),D(1:N,1:N),Z(1:N,1:N))
 allocate(eigenval(1:N),nu(1:N),Psi_s(1:N))
 allocate(realalpha(0:N-1),imagalpha(0:N-1),beta(1:N-1))
 !Para zheev
 lworkk = 2*N-1
 allocate(workk(1:lworkk))
 allocate(rworkk(1:3*n-2))

 !Construcción Hamiltoniano en base Krylov
  open(unit=10,file="coeficientesAlpha.txt",status="old")
  open(unit=11,file="coeficientesBeta.txt",status="old")
  read(10,*)realalpha(0)
  do i = 1, N-1
    read(10,*)realalpha(i),imagalpha(i)
    read(11,*)beta(i)
  end do
  close(10)
  close(11)
realalpha = 0.25_dp*realalpha
imagalpha = 0.25_dp*imagalpha
beta = 0.25_dp*beta
H = (0._dp,0._dp)

H(1,1) = complex(realalpha(0),imagalpha(0))
H(2,1) = complex(beta(1),0._dp)
 do i = 2 ,N-1
   H(i,i) = complex(realalpha(i-1),imagalpha(i-1))
   H(i-1,i) = complex(beta(i-1),0._dp)
   H(i+1,i) = complex(beta(i),0._dp)
 end do
 H(N,N) = complex(realalpha(N-1),imagalpha(N-1))
 H(N-1,N) = complex(beta(N-1),0._dp)
   !CHECK Lectura Coeficinetes
  do i=0,N-1
      write(*,*) "alpha:",realalpha(i),imagalpha(i)
  end do
  do i=1,N-1
      write(*,*) "beta:",beta(i)
  end do


 !CHECK H
 do i=1,N
    write(*,*) H(i,:)
 end do

!Diagonalización
 Z = H
 !call cheevd('V','U',N,Z,N,eigenval,workk,lworkk,rworkk,lrworkk,iworkk,liworkk,information)


 call zheev("V","U",N,Z,N,eigenval,workk,lworkk,rworkk,information)
 
 if (int(workk(1)) /= lworkk .and. information /= 0)then
  lworkk = int(workk(1))
  deallocate(workk)
  allocate(workk(1:lworkk))
  call zheev("V","U",N,Z,N,eigenval,workk,lworkk,rworkk,information)
 end if

 if (information==0) then
 	write(*,*)"------------------------------------------"
 	write(*,*) "Info = 0 - Funcionó bien"
  write(*,*) "lwork óptimo = ", workk(1)
  
 else
  write(*,*)"Info != 0 - La concha de la lora"
 end if

 	write(*,*)"------------------------------------------"
 	write(*,*) "Los autovalores son:"
   
  do i = 1, n
 	  write(*,*) "lambda",i,"=",eigenval(i)
 	end do
 	write(*,*)"------------------------------------------"
 	write(*,*) "La matriz de autovectores es:"
 	do i = 1,n
 		write(*,*) z(i,:)
 	end do
 	write(*,*)"------------------------------------------"





 do i=1,N
    re= cos(eigenval(i)*s)
    im= -sin(eigenval(i)*s)
    nu(i) = complex(re,im)
!    write(*,*) re, im , lambda(i) !Check
 end do

 Psi_s = 0._dp
 do i=1,N
    do k=1,N
        Psi_s(i) = Psi_s(i) + Z(i,k)*nu(k)*conjg(Z(1,k))
    end do
 end do


 !Check Evolución
 open(unit=12, file="CoefEvTemp.txt", status="unknown", action="write")
 write(12,*) " recoef                         imcoef"
 do i=1,N
    write(12,*) REAL(Psi_s(i),8),"    ", DIMAG(Psi_s(i))
 end do

 close(12)
End Program
