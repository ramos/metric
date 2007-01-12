!
! PROGRAM to compute the metric in the space of solutions of the 
!         Bogomolny equations of the Abelian Higgs Model.
!
! Copyright (C) 2006  Alberto Ramos <alberto@martin.ft.uam.es>
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program; if not, write to the Free Software
! Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
! USA
! 

! $ v. 1.0; Released: 10/01/2007; $

Program Metrica 

  USE NumTypes
  USE Constants
  USE Geometry
  USE Error
  USE Fourier
  USE Bases
  USE Time

  Integer, Parameter :: MaxOrder = 51
  Integer :: Nqs = 2, Norders = 20

  Integer :: Nterm = 20, Terminos(MaxOrder)
  Integer, Allocatable :: Multi(:)
  Complex (kind=DPC), Allocatable :: C(:), g(:,:,:), Wzero(:), WzeroRel(:)
  Real (kind=DP) :: Dnorm, tbg, tnd, BCR, BCI
  Real (kind=DP), Allocatable :: Factor(:)
  Complex (kind=DPC) :: Wfac
  Complex (kind=DPC), Allocatable :: B(:,:) !->For ST2 
  Character (len=100) :: dirbase, FileSave, FileMulti

  Type (Fourier_Serie_2D) :: Chif, Prod, &!For Stage 2
       & Lcont, Aux
  Type (Fourier_Serie_2D), Allocatable :: hf(:), &
       & Pf(:), deltaf(:), Acum2(:), Aparcial(:), &
       & Powhf(:,:), & ! For Stage 2
       & Hdot(:,:), Lf(:,:), Acum(:), & ! For Stage 3
       & Lmc(:)
  

  Interface 
     Recursive Function Factorial(N1) Result (Fac)
       
       USE NumTypes

       Integer, Intent (in) :: N1
       Real (kind=DP) :: Fac
  
     End Function Factorial
  End Interface

  Write(stderr,*)
  Write(stderr,*)'Stage 1: Reading input and needed data:'
  Write(stderr,*)'======================================='

  Read(*,*)Nstart1, Nstart2
  Read(*,*)Nqs, Norders, Nterm, dirbase

  ! Allocate space for all the components...
  Allocate(C(Nqs), Wzero(Nqs), WzeroRel(Nqs))
  Allocate(hf(Norders), Pf(Norders), deltaf(Norders), &
       & Acum2(Norders), Aparcial(Nqs), Powhf(Norders,Norders))
  Allocate(Factor(Norders), Multi(Norders))


  CALL Init_Geometry(qq=Nqs)

  Do I = 1, q-1
     Read(*,'(2ES33.25)')Wzero(I)
  End Do
  
  If (q > 1) Then
     Wzero(q) = -Sum(Wzero(1:q-1))
  Else
     Wzero(q) = Cmplx(0.0_DP, 0.0_DP, kind=DPC)
  End If
  WzeroRel = Wzero
  Do I = 1, q
     Wzero(I) = l1/2.0_DP * &
          & Cmplx(0.5_DP+Real(Wzero(I),kind=DP), &
          &       taui*(0.5_DP+Aimag(Wzero(I))), kind=DPC)
  End Do

  C = Cbase(Wzero)
  Dnorm = Sqrt(Sum(Abs(C(:))**2))
  C = C/Dnorm

  Write(stderr,'(1A,1I4)')' Flux (a.k.a q, first chern class): ', q
  Write(stderr,'(1A,1I7)')' Orders to compute: ', Norders
  Write(stderr,'(1A,1I4)')' Number of terms used for the Fourier series: ', Nterm
  Write(stderr,*)
  Write(stderr,'(1A)')' Positions of the zeros (relative to the center of the Torus): '
  Do I = 1, q
     Write(stderr,'(2ES33.25)')WzeroRel(I)
  End Do
  Write(stderr,*)

  ! Set values to Terminos(51)
  Data Terminos /0, 1, 2, 3, 5, 7, 11, 15, 22, 30, 42, 56, 77, 101,&
       & 135, 176, 231, 297, 385, 490, 627, 792, 1002, 1255, 1575, &
       & 1958, 2436, 3010, 3718, 4565, 5604, 6842, 8349, 10143, &
       & 12310, 14883, 17977, 21637, 26015, 31185, 37338, 44583, &
       & 53174, 63261, 75175, 89134, 105558, 124754, 147273, 173525, &
       & 204226/

  ! Write a small summary in the dir, with date 
  ! and other interesting things
  Write(FileSave,'(1A)')'summary_metric.dat'
  FileSave = Trim(Trim(dirbase) // '/' // FileSave)

  Open (Unit = 69, File = Trim(FileSave))

  Write(69,'(1A)')'# '
  Write(69,'(2A)')'# SUMMARY of computed data. ', &
       & asctime(gettime())
  Write(69,'(1A)')'# '
  Write(69,'(1A,1I4)')'# Flux (a.k.a q, first chern class): ', q
  Write(69,'(1A,1I7)')'# Orders computed: ', Norders
  Write(69,'(1A,1I4)')'# Number of terms used for the Fourier series: ', Nterm
  Write(69,*)
  Write(69,'(1A)')'# Positions of the zeros (relative to the center of the Torus): '
  Do I = 1, q
     Write(69,'(2ES33.25)')WzeroRel(I)
  End Do

  Write(69,*)
  Write(69,'(1A)')'# Values of C (position in moduli space in the coordinates of CP): '
  Do I = 1, q
     Write(69,'(1A2,2ES33.25)')'# ', C(I)
  End Do
  

  ! Init the fourier Series
  CALL Init_Serie(Chif, Nterm)
  CALL Init_Serie(Prod, Nterm)
  
  Do I = 1, q
     CALL Init_Serie(Aparcial(I), Nterm)
  End Do
    
  Do I = 1, Norders
     CALL Init_Serie(hf(I), Nterm)
     CALL Init_Serie(Deltaf(I), Nterm)
     CALL Init_Serie(Acum2(I), Nterm)
     CALL Init_Serie(Pf(I), Nterm)
     Do K = 1, Norders
        CALL Init_Serie(Powhf(I,K), Nterm)
     End Do
  End Do

  ! Set the displacements and
  ! calculate the initial function.
  Do N1 = -Nterm, Nterm
     Do N2 = -Nterm, Nterm
        Chif%Coef(N1, N2) = (0.0_DP, 0.0_DP)
        Do I = 1, q
           Do J = 1, q
              Chif%Coef(N1,N2) = Chif%Coef(N1,N2) + &
                   & L(i, j, N1, N2)*Conjg(C(i))*C(j)
           End Do
        End Do
     End Do
  End Do
  
  I = 1
  ! Set "maually" the first order.
  Factor(1) = Real(1.0_DP/Chif%Coef(0,0), kind=DP)
  
  Deltaf(1) = -0.5_DP * Factor(1) * Chif
  Deltaf(1)%Coef(0,0) = (0.0_DP, 0.0_DP)

  Pf(1)%Coef = (0.0_DP, 0.0_DP)
  Pf(1)%Coef(0,0) = Cmplx(Factor(1), kind=DPC)

  
  Do N1 = -Nterm, Nterm
     Do N2 = -Nterm, Nterm
        If ( (N1 == 0) .and. (N2 == 0) ) Then
           hf(1)%Coef(0,0) = (0.0_DP, 0.0_DP)
        Else
           hf(1)%Coef(N1,N2) = - Deltaf(1)%Coef(N1,N2) / xi(N1,N2)
        End If
     End Do
  End Do
  
  Write(FileSave, '(1A9,1I2.2,1A7,1I2.2,1A6,1I1.1,1A4)')'factor:O=',I&
       &,':Nterm=',Nterm ,':flux=',q,'.dat'
  FileSave = Trim(Trim(dirbase) // '/' // FileSave)
  Open (Unit=99, File = FileSave)
  Write(99,'(1ES33.25)')Factor(1)
  Close(99)
  
  Write(FileSave,'(1A9,1I2.2,1A7,1I2.2,1A6,1I1.1,1A4)')'deltaf:O=',I&
       &,':Nterm=',Nterm ,':flux=',q,'.dat'
  FileSave = Trim(Trim(dirbase) // '/' // FileSave)
  CALL Save_Serie(Deltaf(I), Trim(FileSave))
  
  Write(FileSave,'(1A5,1I2.2,1A7,1I2.2,1A6,1I1.1,1A4)')'hf:O=',I&
       &,':Nterm=',Nterm ,':flux=',q,'.dat'
  FileSave = Trim(Trim(dirbase) // '/' // FileSave)
  CALL Save_Serie(hf(I), Trim(FileSave))
  
  Powhf(1,1) = hf(1)
  Do I2 = 2, Int((Norders-1))
     Powhf(1,I2) = hf(1) * Powhf(1,I2-1)
  End Do
  
  Write(69,'(1A)')'# '
  Write(69,'(1A)')'# Times of the computation (PART I):'
  Write(69,'(1A)')'# '
  ! If we have to read data from HD, its time to do it
  If (Nstart1 > 1) Then
     Do I = 2, Nstart1
        CALL CPU_Time(tbg)
        Write(FileSave,'(1A5,1I2.2,1A7,1I2.2,1A6,1I1.1,1A4)')'Pf:O=',I&
             &,':Nterm=',Nterm ,':flux=',q,'.dat'
        FileSave = Trim(Trim(dirbase) // '/' // FileSave)
        CALL Read_Serie(Pf(I), Trim(FileSave))
        
        Write(FileSave, '(1A9,1I2.2,1A7,1I2.2,1A6,1I1.1,1A4)')'factor:O=',I&
             &,':Nterm=',Nterm ,':flux=',q,'.dat'
        FileSave = Trim(Trim(dirbase) // '/' // FileSave)
        Open (Unit=99, File = FileSave)
        Read(99,'(1ES33.25)')Factor(I)
        Close(99)
        
        Write(FileSave,'(1A9,1I2.2,1A7,1I2.2,1A6,1I1.1,1A4)')'deltaf:O=',I&
             &,':Nterm=',Nterm ,':flux=',q,'.dat'
        FileSave = Trim(Trim(dirbase) // '/' // FileSave)
        CALL Read_Serie(Deltaf(I), Trim(FileSave))
        
        Write(FileSave,'(1A5,1I2.2,1A7,1I2.2,1A6,1I1.1,1A4)')'hf:O=',I&
             &,':Nterm=',Nterm ,':flux=',q,'.dat'
        FileSave = Trim(Trim(dirbase) // '/' // FileSave)
        CALL Read_Serie(hf(I), Trim(FileSave))

        Powhf(I,1) = hf(I)
        Do Npow = 2, Int((Norders-1)/I)
           Powhf(I, Npow) = hf(I) * Powhf(I,Npow-1)
        End Do

        Acum2(I)%Coef = Cmplx(0.0_DP, kind=DPC)
        Do J = 1, I-2
           Acum2(I) = Acum2(I) + Factor(I-J)*Acum2(J+1)
        End Do
        Acum2(I) = (1.0_DP/Factor(1))*(Pf(I) - Acum2(I))
        Acum2(I)%Coef(0,0) = Acum2(I)%Coef(0,0) - Factor(I)/Factor(1)

        CALL CPU_Time(tnd)
        Write(stderr,'(1A,1I5)')'Reading order: ', I
        Write(69,'(1A,1I5,1A,1ES8.2,1A)')'Reading order: ', I, &
             & '(',tnd-tbg,'sec)'
     End Do
  End If

  Write(stderr,*)
  Write(stderr,*)'Stage 2: Computing h, Factor:'
  Write(stderr,*)'==============================='
  
  ! Now start the iteration until Norders is reached.
  Do I = Max(2,Nstart1+1), Norders
     ! This loop must be done for each
     ! order. I is the order.
     CALL CPU_Time(tbg)
     Write(FileMulti, '(1A15,1I2.2,1A4)')'combinatoria/O=', I-1, '.dat'     
     Open (Unit = 34, File = Trim(FileMulti), ACTION="READ")     
     
     Acum2(I)%Coef = (0.0_DP, 0.0_DP)
     Do J = 1, Terminos(I)
        Read(34,*)(Multi(K), K=1, Norders-1)
        Ntot = Sum(Multi)
        
        CALL Unit(Prod,Nterm)
        Do K = 1, I-1
           If (Multi(K) /= 0) Then
              Prod = (1.0_DP/Factorial(Multi(K))) * Prod * &
                   & Powhf(K,Multi(K))
           End If
        End Do
        Prod = (-2.0_DP) ** Ntot * Prod
        
        Acum2(I) = Acum2(I) + Prod
     End Do
     close(34)
     
     Prod%Coef = (0.0_DP, 0.0_DP)
     Do J = 1, I-1
        Prod = Prod + Factor(I-J) * Acum2(J+1)
     End Do
     
     Pf(I) = Prod
     Prod = Prod * Chif
     
     Factor(I) = -Real(Prod%Coef(0,0),kind=DP) / Real(Chif%Coef(0,0),kind=DP)
     Deltaf(I) = -0.5_DP * ( (Factor(I)*Chif) + Prod )

     Pf(I)%Coef(0,0) = Pf(I)%Coef(0,0) + Factor(I)
     Write(FileSave,'(1A5,1I2.2,1A7,1I2.2,1A6,1I1.1,1A4)')'Pf:O=',I&
          &,':Nterm=',Nterm ,':flux=',q,'.dat'
     FileSave = Trim(Trim(dirbase) // '/' // FileSave)
     CALL Save_Serie(Pf(I), Trim(FileSave))

     
     Do N1 = -Nterm, Nterm
        Do N2 = -Nterm, Nterm
           If ( (N1 == 0) .and. (N2 == 0) ) Then
              hf(I)%Coef(0,0) = (0.0_DP, 0.0_DP)
           Else
              hf(I)%Coef(N1,N2) = hf(I-1)%Coef(N1,N2) - &
                   & Deltaf(I)%Coef(N1,N2) / xi(N1,N2)
           End If
        End Do
     End Do
     
     Write(FileSave, '(1A9,1I2.2,1A7,1I2.2,1A6,1I1.1,1A4)')'factor:O=',I&
          &,':Nterm=',Nterm ,':flux=',q,'.dat'
     FileSave = Trim(Trim(dirbase) // '/' // FileSave)
     Open (Unit=99, File = FileSave)
     Write(99,'(1ES33.25)')Factor(I)
     Close(99)
     
     Write(FileSave,'(1A9,1I2.2,1A7,1I2.2,1A6,1I1.1,1A4)')'deltaf:O=',I&
          &,':Nterm=',Nterm ,':flux=',q,'.dat'
     FileSave = Trim(Trim(dirbase) // '/' // FileSave)
     CALL Save_Serie(Deltaf(I), Trim(FileSave))
     
     Write(FileSave,'(1A5,1I2.2,1A7,1I2.2,1A6,1I1.1,1A4)')'hf:O=',I&
          &,':Nterm=',Nterm ,':flux=',q,'.dat'
     FileSave = Trim(Trim(dirbase) // '/' // FileSave)
     CALL Save_Serie(hf(I), Trim(FileSave))
     
     Powhf(I,1) = hf(I)
     Do Npow = 2, Int((Norders-1)/I)
        Powhf(I, Npow) = hf(I) * Powhf(I,Npow-1)
     End Do
     
     CALL CPU_Time(tnd)
     Write(stderr, '(1A,1I5,1A)')'Order: ', I, ' done'
     Write(69,'(1A,1I5,1A,1ES8.2,1A)')'Order: ', I, &
          & ' done (',tnd-tbg,'sec)'
  End Do


  Write(69,'(1A)')'# '
  Write(69,'(1A)')'# Times of the computation (PART II):'
  Write(69,'(1A)')'# '

  Write(0,*)
  Write(0,*)'Stage 3: Computing \dot H_i, B_i:'
  Write(0,*)'==================================='
  

  ! Allocate space for stage 2 and init the new Series
  ALLOCATE(   B(q,0:Norders))
  ALLOCATE(Hdot(q,Norders))
  ALLOCATE(   Lf(q,q))

  CALL Init_Serie(Lcont, Nterm)
  CALL Init_Serie(Aux  , Nterm)
  Do J1 = 1, q
     Do J2 = 1, q
        CALL Init_Serie(Lf(J1,J2), Nterm)
     End Do
  End Do
  Do I = 1, Norders
     Do J2 = 1, q
        CALL Init_Serie(Hdot(J2,I), Nterm)
     End Do
  End Do

  ! First of all set zero order things.
  Lcont%Coef = (0.0_DP, 0.0_DP)
  Do J1 = 1, q
     Do J2 = 1, q
        Do N1 = -Nterm, Nterm
           Do N2 = -Nterm, Nterm
              Lf(J1,J2)%Coef(N1,N2) = L(J1,J2,N1,N2)
           End Do
        End Do
     End Do
  End Do
  
  Do J1 = 1, q
     Do J2 = 1, q
        Lcont = Lcont + Conjg(C(J1))*C(J2) * Lf(J1,J2)
     End Do
  End Do
  
  
  ! Set and save in HD B(:,0) manually.
  Dnorm = Abs(C(1))**2 + Abs(C(2))**2
  I = 0 ! The order
  Do J1 = 1, q
     B(J1,0) = Conjg(C(J1))/Dnorm
     
     Write(FileSave,'(1A3,1I1.1,1A3,1I2.2,1A7,1I2.2,1A6,1I1.1,1A4)') &
          & 'Bf:',J1,':O=',I,':Nterm=',Nterm ,':flux=',q,'.dat'
     FileSave = Trim(Trim(dirbase) // '/' // FileSave)
     Open (Unit=99, File = FileSave, Action = "WRITE")
     Write(99,'(1ES33.25)')B(J1,0)
     Close(99)
  End Do
  
  
  ! Set and save in HD H(I,1).
  I = 1 ! The order
  Do J1 = 1, q
     Aux%Coef = (0.0_DP, 0.0_DP)
     Do J2 = 1, q
        Aux = Aux + Lf(J2,J1)*Conjg(C(J2))
        !          Write(*,*)J1,J2,Aux%Coef(1,2)
     End Do
     Aux = Aux - Lcont*B(J1,0)
     !       Aux = Aux - (Conjg(C(J1))/Dnorm)*Lcont
     !       Write(*,*)J1,J2,Aux%Coef(1,2)       
     
     Do Nx = -Nterm, Nterm
        Do Ny = -Nterm, Nterm
           If ( (Nx /= 0) .or. (Ny /= 0) ) Then
              Hdot(J1,1)%Coef(Nx,Ny) = Aux%Coef(Nx,Ny)/(Dnorm*xi(Nx,Ny))
           Else
              Hdot(J1,1)%Coef(Nx,Ny) = (0.0_DP, 0.0_DP)
           End If
        End Do
     End Do
     
     Write(FileSave,'(1A6,1I1.1,1A3,1I2.2,1A7,1I2.2,1A6,1I1.1,1A4)') &
          & 'hdotf:',J1,':O=',I,':Nterm=',Nterm ,':flux=',q,'.dat'
     FileSave = Trim(Trim(dirbase) // '/' // FileSave)
     CALL Save_Serie(Hdot(J1,1), Trim(FileSave))
  End Do
  
  ! Read the computed things
  Do I = 2, Nstart2
     CALL CPU_Time(tbg)
     Do J1 = 1, q
        Write(FileSave,'(1A6,1I1.1,1A3,1I2.2,1A7,1I2.2,1A6,1I1.1,1A4)') &
             & 'hdotf:',J1,':O=',I,':Nterm=',Nterm ,':flux=',q,'.dat'
        FileSave = Trim(Trim(dirbase) // '/' // FileSave)
        CALL Read_Serie(Hdot(J1,I), Trim(FileSave))
        
        Write(FileSave,'(1A3,1I1.1,1A3,1I2.2,1A7,1I2.2,1A6,1I1.1,1A4)') &
             & 'Bf:',J1,':O=',I-1,':Nterm=',Nterm ,':flux=',q,'.dat'
        FileSave = Trim(Trim(dirbase) // '/' // FileSave)
        Open (Unit=99, File = FileSave, Action = "READ")
        Read(99,'(1ES33.25)')BCR
        Read(99,'(1ES33.25)')BCI
        Close(99)

        B(J1,I-1) = Cmplx(BCR, BCI, kind=DPC)
     End Do
     CALL CPU_Time(tnd)
     Write(stderr,'(1A,1I5)')'Reading order: ', I
     Write(69,'(1A,1I5,1A,1ES8.2,1A)')'Reading order: ', I, &
          & '(',tnd-tbg,'sec)'
  End Do


  ! Now we have to iterate the perturbative process.
  ! In each loop we calculate:
  !  - Hdot(:,I)
  !  - B(:,I-1)
  ! and save them to HD.
  Do I = Max(2, Nstart2), Norders
     CALL CPU_Time(tbg)

     Do J1 = 1, q
        Aux%Coef = (0.0_DP,0.0_DP)
        
        Aux = -B(J1,0)*Pf(I)
        Do M = 1, I-2
           Aux = Aux - Pf(I-M)*Hdot(J1,M) - &
                & Pf(I-M)*B(J1,M)
        End Do
        Aux = Aux - Factor(1)*Hdot(J1,I-1)
        Aux = Aux * Lcont
        
        Do J2 = 1, q
           Aux = Aux + Pf(I)*Lf(J2,J1)*Conjg(C(J2))
        End Do
        
        B(J1,I-1) = Aux%Coef(0,0)/( Factor(1)*Lcont%Coef(0,0) )
        
        Aux = Aux - Factor(1)*B(J1,I-1)*Lcont
!        Write(stderr,*)I, Aux%Coef(0,0)!, Aux%Coef(1,1)
        
        Do Nx = -Nterm, Nterm
           Do Ny = -Nterm, Nterm
              If( (Nx == 0) .and. (Ny == 0) ) Then
                 Hdot(J1,I)%Coef(Nx,Ny) = (0.0_DP,0.0_DP)
              Else
                 Hdot(J1,I)%Coef(Nx,Ny) = Aux%Coef(Nx,Ny)/xi(Nx,Ny)&
                      & + Hdot(J1,I-1)%Coef(Nx,Ny)
              End If
           End Do
        End Do
        
        !          Write(stderr,'(2I3,6ES10.3)')I,J1,Aux%Coef(1,1)/xi(1,1)&
        !               &,Hdot(J1,I)%Coef(1,1),&
        !               & Hdot(J1,I-1)%Coef(1,1)
        
        
        Write(FileSave,'(1A6,1I1.1,1A3,1I2.2,1A7,1I2.2,1A6,1I1.1,1A4)') &
             & 'hdotf:',J1,':O=',I,':Nterm=',Nterm ,':flux=',q,'.dat'
        FileSave = Trim(Trim(dirbase) // '/' // FileSave)
        CALL Save_Serie(Hdot(J1,I), Trim(FileSave))
        
        Write(FileSave,'(1A3,1I1.1,1A3,1I2.2,1A7,1I2.2,1A6,1I1.1,1A4)') &
             & 'Bf:',J1,':O=',I-1,':Nterm=',Nterm ,':flux=',q,'.dat'
        FileSave = Trim(Trim(dirbase) // '/' // FileSave)
        Open (Unit=99, File = FileSave, Action = "WRITE")
        Write(99,'(1ES33.25)')B(J1,I-1)
        Close(99)
        
     End Do

     CALL CPU_Time(tnd)
     Write(stderr, '(1A,1I5,1A)')'Order: ', I, ' done'
     Write(69,'(1A,1I5,1A,1ES8.2,1A)')'Order: ', I, &
          & ' done (',tnd-tbg,'sec)'
  End Do
  
  Write(0,*)
  Write(0,*)'Stage 4: Computing the metric:'
  Write(0,*)'================================'

  Write(69,'(1A)')'# '
  Write(69,'(1A)')'# Times of the computation (PART III):'
  Write(69,'(1A)')'# '
  

  ! Allocate and init things for stg3
  Allocate(g(Norders,q,q))
  g = (0.0_DP, 0.0_DP)
  ALLOCATE(  Lmc(q))

  Do I = 1, q
     CALL Init_Serie(Lmc(I), Nterm)
  End Do

  Do J1 = 1, q
     Lmc(J1)%Coef = (0.0_DP,0.0_DP)
     Do J2 = 1, q
        Lmc(J1) = Lmc(J1) + Lf(J1,J2) * C(J2)
     End Do
  End Do


  ! Calculates the metric
  Do I = Max(1,Nstart2), Norders
     CALL CPU_Time(tbg)
     Do J1 = 1, q
        Do J2 = 1, q
           Aux = Pf(I)*Lf(J1,J2)
           g(I,J1,J2) = Aux%Coef(0,0)
!           if ((J1*J2 == 1).and.(I==2)) Write(0,'(1A3,2ES33.25)')'1: ', Aux%Coef(0,0)
           Aux = - Pf(I)*B(J2,0)*Lmc(J1)
!           if ((J1*J2 == 1).and.(I==2)) Write(0,'(1A3,2ES33.25)')'2: ', Aux%Coef(0,0)
           g(I,J1,J2) = g(I,J1,J2) + Aux%Coef(0,0)
           Do M = 1, I-1
              Aux = - Pf(I-M) * B(J2,M) * Lmc(J1) 
              g(I,J1,J2) = g(I,J1,J2) + Aux%Coef(0,0)
!              if ((J1*J2 == 1).and.(I==2)) Write(0,'(1A3,2ES33.25)')'3: ', Aux%Coef(0,0)
              Aux = - Pf(I-M) * Hdot(J2,M) * Lmc(J1)
              g(I,J1,J2) = g(I,J1,J2) + Aux%Coef(0,0)
!              if ((J1*J2 == 1).and.(I==2)) Write(0,'(1A3,2ES33.25)')'4&
!                   &: ', Aux%Coef(0,0)
           End Do
           
!           g(I,J1,J2) = g(I,J1,J2) + Aux%Coef(0,0)
        End Do
     End Do

     ! Save the metric in HD with a "nice" format
     Write(FileSave,'(1A9,1I2.2,1A7,1I2.2,1A6,1I1.1,1A4)') &
          & 'metric:O=',I,':Nterm=',Nterm ,':flux=',q,'.dat'
     FileSave = Trim(Trim(dirbase) // '/' // Trim(FileSave))
     Open (Unit = 99, File = Trim(FileSave), Action = 'WRITE')
     Do I1 = 1, q
        Write(99,'(100ES33.25)')(g(I,I1,J), J = 1, q)
     End Do
     Close(99)
     CALL CPU_Time(tnd)
     Write(stderr, '(1A,1I5,1A)')'Order: ', I, ' done'
     Write(69,'(1A,1I5,1A,1ES8.2,1A)')'Order: ', I, &
          & ' done (',tnd-tbg,'sec)'
  End Do

!  Write(*,'(8ES33.25,500ES33.25)')Wzero(1), C(1), &
!       & C(2), C(2)/C(1), g(:,1,1) + g(:,2,2)

  Close(69)
  
  Stop
End Program Metrica

! ****************************
! *
Recursive Function Factorial(N1) Result (Fac)
! *
! ****************************
! * Returns the factorial of a 
! * integer number.
! ****************************
  
  USE NumTypes

  Integer, Intent (in) :: N1
  Real (kind=DP) :: Fac
  
  
  If (N1 == 0) Then 
     Fac = 1.0_DP
  Else
     Fac = Real(N1,kind=DP) * Factorial(N1 - 1)
  End If
  
  Return
End Function Factorial
