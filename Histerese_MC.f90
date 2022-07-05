!! Artigo referência: Diversity Enabling Equilibration: Disorder and the Ground State in Artificial Spin Ice
!!                         Zoe Budrikis, Paolo Politi and R. L. Stamps
!!
!!
!!    ListStreamPlot -> data = Partition[#, 2] & /@ data;


module ranmod1
   integer,parameter::i4b=selected_int_kind(9)
   integer,parameter::k15=selected_int_kind(15)
   integer,parameter::dp=selected_real_kind(8)
end module ranmod1
  
module ranmod2
   use ranmod1
   real(dp)::U(97),cc,cd,cm
   integer(k15)::i97,j97
end module ranmod2
  
module ranutil
   use ranmod1
   integer(i4b)::rand_inst=0
   integer(i4b)::rand_feedback=1
  contains
   subroutine initrandom(i,i2)
    use ranmod1
    implicit none
    integer(k15),optional,intent(in)::i,i2
    integer(k15)::seed_in,kl,ij
    character(len=10)::fred
    real(dp)::klr
    
    if (present(i)) then
     seed_in = i
    else
     seed_in = -1
    end if
    if (seed_in /= -1) then
     if (present(i2)) then
      kl = i2
      if (i2 > 30081) stop 'i2 > 30081'
     else
      kl = 9373
     end if
     ij = i
    else
     call system_clock(count=ij)
     ij = mod(ij+rand_inst*100,31328)
     call date_and_time(time=fred)
     read(fred,'(e10.3)')klr
     kl = mod(int(klr*1000),30081)
    end if
    call rmarin(ij,kl)
   end subroutine initrandom
   
   subroutine rmarin(ij,kl)
    use ranmod1
    use ranmod2
    implicit none
    integer(k15),intent(in)::ij,kl
    integer(k15)::i,j,k,l,ii,jj,m
    real(dp)::s,t
    
    if (ij<0 .or. ij>31328 .or. kl<0 .or. kl>30081) then
     stop 'Erro ao iniciar o gerador de números aleatórios'
    end if
    i = mod(ij/177,177)+2
    j = mod(ij,177)+2
    k = mod(kl/169,178)+1
    l = mod(kl,169)
    do ii = 1,97
     s = 0.0d0
     t = 0.5d0
     do jj = 1,24
      m = mod(mod(i*j,179)*k,179)
      i = j
      j = k
      k = m
      l = mod(53*l+1,169)
      if (mod(l*m,64)>=32) then
       s = s + t
      end if
      t = 0.5d0*t
     end do
     u(ii) = s
    end do
    cc = 362436.0d0/16777216.0d0
    cd = 7654321.0d0/16777216.0d0
    cm = 16777213.0d0/16777216.0d0
    i97 = 97
    j97 = 33
    
    return
   end subroutine rmarin
   
   function ranmar() result(fn_val)
    use ranmod1
    use ranmod2
    implicit none
    real(dp) :: fn_val
    real(dp) :: uni
    
    uni =u(i97) -u(j97)
    if ( uni < 0.0d0 ) uni = uni + 1.0d0
    u(i97) = uni
    i97 = i97 - 1
    if (i97 == 0) i97 = 97
    j97 = j97 - 1
    if (j97 == 0) j97 = 97
    cc = cc - cd
    if ( cc < 0.0d0 ) cc = cc + cm
    uni = uni - cc
    if ( uni < 0.0d0 ) uni = uni + 1.0d0
    fn_val = uni
    
    return
   end function ranmar
   
   real(dp) function rand_real(a,b)
    use ranmod1
    implicit none
    real(dp),intent(in)::a,b
    
    rand_real = (b-a)*ranmar() + a
   end function rand_real
   
   integer(i4b) function rand_int(i,j)
    use ranmod1
    implicit none
    integer(i4b),intent(in)::i,j
    
    rand_int = int(j*ranmar()) + i
   end function rand_int

   function gaussrnd()
      implicit none
      real(dp) :: gaussrnd
      real(dp) :: fac,v1,v2,r
      real(dp), save :: gset
      integer, save :: iset=0
  
      if (iset==0) then ! Create a new RN
         r=100.0
         do while (r>1.0)
            v1 = 2.0*ranmar()-1.0
            v2 = 2.0*ranmar()-1.0
            r = v1*v1+v2*v2
         end do
         fac = sqrt(-2.0*log(r)/r)
         gset = v1*fac
         gaussrnd = v2*fac
         iset = 1
      else ! Use the 2nd NR from the previous call
         gaussrnd = gset
         iset = 0
      endif
      return
    end function gaussrnd
   
end module ranutil

  !-----------------------------------------------------------------!
  !                          MARSAGLIA                              !
  ! Usage example:                                                  !
  !    use ranutil                                                  !
  !    ...                                                          !
  !    integer :: i                                                 !
  !    real(rk) :: x,g                                              !
  !    ...                                                          !                                           !
  !    call initrandom()                                            !
  !    x=ranmar()                                                   !
  !    x=rand_real(25.,59.)                                         !
  !    i=rand_int(0,100)                                            !
  !    ...                                                          !
  !                                                                 !
  !-----------------------------------------------------------------!

module var_global
   real(8), parameter :: pi = 4.0d0*atan(1.0d0)
   integer :: Ns,N_viz,N_campo,ncor,Nmssf
   integer, dimension(:), allocatable :: S,Nviz,jviz,cor
   real(8), dimension(:), allocatable :: Aij,Bi,Hi,hcut
   real(8), dimension(:), allocatable :: x,y,mx,my
   real(8), dimension(:), allocatable :: Bx,By,theta,ModB
   integer,dimension(:,:),allocatable :: cortab
   real(8) :: BMax,delB,ExtMax,theta_c,D
   character(600) :: dir1,dir2
   integer :: kk
end module var_global

subroutine inicial
   use var_global
   use ranutil, only : initrandom
   implicit none
   integer :: qual
   logical :: direx

   !!Ler input

   call initrandom()

   print*, "Histerese = 0; Campo Rotativo = 1"
   read(*,*) qual
   print*, "Número para MSSF:"; read(*,*) Nmssf

   call ler_input(1)
   call ler_config(2)
   call ler_Aij(3)
   call ler_Nviz(4)
   call ler_cortab(5)
   !!Calculo dos campos iniciais
   allocate(Bi(Ns),Hi(Ns))
   call campo_dip
   call campo_flip

   if (qual == 0) then
      call campo_zee_H
   else if (qual == 1) then
      call campo_zee
   else
      print*, "Erro"
      stop
   end if
   
   !!Diretórios

   write(dir1,"('Resultados/')")
   inquire(file=dir1,exist=direx)
   if (direx .eqv. .false.) then
      call system("mkdir " // trim(dir1))
   end if

   write(dir1,"('Resultados/His/')")
   inquire(file=dir1,exist=direx)
   if (direx .eqv. .false.) then
      call system("mkdir " // trim(dir1))
   end if

   write(dir1,"('Resultados/Rot/')")
   inquire(file=dir1,exist=direx)
   if (direx .eqv. .false.) then
      call system("mkdir " // trim(dir1))
   end if

   if (qual == 0) then
      kk = 0
      direx = .true.
      do while (direx .eqv. .true.)
         kk = kk + 1
         write(dir2,"('Resultados/His/Simulacao_',I3.3,'/')") kk
         inquire(file=dir2,exist=direx)
      end do
   else
      kk = 0
      direx = .true.
      do while (direx .eqv. .true.)
         kk = kk + 1
         write(dir2,"('Resultados/Rot/Simulacao_',I3.3,'/')") kk
         inquire(file=dir2,exist=direx)
      end do
   end if
   call system ("mkdir " // trim(dir2))

   return

end subroutine inicial

subroutine ler_input(unit_in)
   use var_global, only : Ns,N_viz,BMax,delB,ExtMax,theta_c,D
   implicit none
   integer, intent(in) :: unit_in
   real(8) :: a
   
   open(unit=unit_in,file='input.dat',status='old',action='read')
   ! read(unit_in,*) D
   ! read(unit_in,*) Ns
   ! read(unit_in,*) N_viz
   
   read(unit_in,*) D
   read(unit_in,*) Ns
   read(unit_in,*) N_viz
   read(unit_in,*) a
   read(unit_in,*) a
   read(unit_in,*) a
   read(unit_in,*) a
   read(unit_in,*) a
   read(unit_in,*) a
   read(unit_in,*) a
   read(unit_in,*) a
   read(unit_in,*) a
   read(unit_in,*) a
   read(unit_in,*) a
   read(unit_in,*) a
   read(unit_in,*) a
   read(unit_in,*) a
   read(unit_in,*) BMax
   read(unit_in,*) delB
   read(unit_in,*) ExtMax
   read(unit_in,*) theta_c
   close(unit_in)

   return
end subroutine ler_input

subroutine ler_config(unit_in)
   use ranutil, only : ranmar
   use var_global, only : Ns,S,x,y,mx,my,cor
   implicit none
   integer, intent(in) :: unit_in
   integer :: i
   real(8) :: z

   allocate(S(Ns),x(Ns),y(Ns),mx(Ns),my(Ns),cor(Ns))
   open(unit=unit_in,file='config0.xyz',status='old',action='read')
   read(unit_in,*) i
   read(unit_in,*)
   do i = 1,Ns
      read(unit_in,*) x(i),y(i),mx(i),my(i),cor(i)
   end do
   close(unit_in)

   print*, 'Config dentro'

   S(:) = 1
   do i = 1,Ns
      z = ranmar()
      if (z .lt. 0.5d0) then
         S(i) = -1
      end if
      
   end do

   return
end subroutine ler_config

subroutine ler_Aij(unit_in)
   use var_global, only : N_viz,Aij,jviz
   implicit none
   integer, intent(in) :: unit_in
   integer :: i

   allocate(Aij(N_viz),jviz(N_viz))
   open(unit=unit_in,file='Aij.dat',status='old',action='read')
   do i = 1,N_viz
      read(unit_in,*) jviz(i),Aij(i)
   end do
   close(unit_in)

   return
end subroutine ler_Aij

subroutine ler_Nviz(unit_in)
   use var_global, only : Ns,Nviz
   implicit none
   integer, intent(in) :: unit_in
   integer :: i

   allocate(Nviz(0:Ns))
   open(unit=unit_in,file='Nviz.dat',status='old',action='read')

   Nviz(0) = 0
   do i = 1,Ns
      read(unit_in,*) Nviz(i)
   end do

   close(unit_in)
   return
end subroutine ler_Nviz

subroutine ler_cortab(unit_in)
   use var_global, only : ncor,cortab
   implicit none
   integer,intent(in) :: unit_in
   integer :: i

   open(unit=unit_in,file='cor.dat',status='old',action='read')
   read(unit_in,*) ncor
   allocate(cortab(-1:1,ncor))
   cortab = 0
   do i = 1,ncor
      read(unit_in,*) cortab(1,i),cortab(-1,i)
   end do
   close(unit_in) 

end subroutine ler_cortab

subroutine campo_dip
   use var_global, only : Ns,S,Nviz,jviz,Aij,Bi
   implicit none
   integer :: i,j,k

   Bi(:) = 0.0d0
   do i = 1,Ns
      do k = Nviz(i-1)+1,Nviz(i)
         j = jviz(k)
         Bi(i) = Bi(i) + S(j)*Aij(k)
      end do
   end do

   return
end subroutine campo_dip

subroutine campo_zee_H
   use var_global, only : N_campo,ExtMax,Bx,By,pi,theta_c,ModB
   implicit none
   integer :: i
   real(8) :: dB,Bmod

   N_campo = 10**5!256
   !print*, 'Entre com N_campo'
   !read(6,*) N_campo


   allocate(Bx(N_campo),By(N_campo),ModB(N_campo))
   dB = 2.0d0*ExtMax/real(N_campo-1,8)

   do i = 1,N_campo/2
      Bmod = -ExtMax + 2.0d0*dB*(i-1)
      ModB(i) = Bmod
      Bx(i) = Bmod*cos(theta_c*pi/180.0d0)
      By(i) = Bmod*sin(theta_c*pi/180.0d0)
   end do

   do i = N_campo/2+1,N_campo
      Bmod = ExtMax - 2.0d0*dB*(i-1-N_campo/2)
      ModB(i) = Bmod
      Bx(i) = Bmod*cos(theta_c*pi/180.0d0)
      By(i) = Bmod*sin(theta_c*pi/180.0d0)
   end do

   return
end subroutine campo_zee_H

subroutine campo_zee
   use var_global, only : N_campo,ExtMax,Bx,By,pi,theta,ModB
   implicit none
   integer :: i

   N_campo = 10**7!256
   print*, 'Entre com N_campo'
   read(*,*) N_campo

   allocate(Bx(1:N_campo),By(1:N_campo),theta(1:N_campo),ModB(1:N_campo))

   do i = 1,N_campo
      theta(i) = i*pi/180.0d0
      Bx(i) = ExtMax*cos(i*pi/180.0d0)
      By(i) = ExtMax*sin(i*pi/180.0d0)
      ModB(i) = sqrt(Bx(i)**2 + By(i)**2)
   end do

   return
end subroutine campo_zee

subroutine campo_flip
   use ranutil, only : gaussrnd
   use var_global, only : Ns,hcut
   implicit none
   integer :: i

   allocate(hcut(Ns))
   !do i = 1,Ns
   !   hcut(i) = delB*gaussrnd() + Bmax
   !end do
   open(50,file='hc.dat',status='old',action='read')
   do i = 1,Ns
      read(50,*) hcut(i)
   end do
   close(50)

   return
end subroutine campo_flip

subroutine abrir_fechar(i,a,b,c)
   use var_global, only : dir2
   implicit none
   integer,intent(in) :: i,a,b,c

   if (i .eq. 1) then
      open(unit = a,file=trim(dir2) // trim("conf.xyz"))
      open(unit = b,file=trim(dir2) // trim("mag.dat"))
      open(unit = c,file=trim(dir2) // trim("campo.xyz"))

      open(100,file=trim(dir2) // trim("Ising.dat"))

   else
      close(a)
      close(b)
      close(c)
   end if

   return
end subroutine abrir_fechar

subroutine budrikis(try)
   use ranutil, only : rand_int
   use var_global, only : Ns,S,Bi,Hi,hcut
   implicit none
   integer, intent(out) :: try
   integer :: SetPos(Ns)
   integer :: i
   real(8) :: E_try

   !!Percorre toda a rede determinando todo o conjunto de spins que são flipáveis!!
   try = 0
   do i = 1,Ns
      E_try = S(i)*(Bi(i) - Hi(i))
      if (E_try .gt. hcut(i)) then !.and. E_try.gt.maxE
         !maxE = E_try
         !j = i
         try = try + 1    !! Como documentado no artigo de Zoe Budrikis em 2012
         SetPos(try) = i
      end if
   end do
   !call update(j)
   !!Flipa um spin aleatório pertencente ao conjunto e atualiza os campos!!
   if (try .ne. 0) then
      i = rand_int(1,try)
      call update(SetPos(i))
   end if

   return
end subroutine budrikis

subroutine update(i)
   use var_global, only : S,Nviz,Aij,Bi,jviz
   implicit none
   integer,intent(in) :: i
   integer :: j,k
   real(8) :: dB

   !! Atualiza o Campo de i nos outros Spin!!
   do k = Nviz(i-1)+1,Nviz(i)
      j = jviz(k)
      dB = -2.0d0*S(i)*Aij(k)
      Bi(j) = Bi(j) + dB
   end do
   !! Atualiza o sentido do Spin i
   S(i) = -S(i)
   
   return
end subroutine update

subroutine config(u_conf)
   use var_global, only : Ns,x,y,mx,my,S,Bi,cor,cortab
   implicit none
   integer, intent(in) :: u_conf
   integer :: i

   write(u_conf,*) Ns
   write(u_conf,*) " "
   do i = 1,Ns
      write(u_conf,*) x(i),y(i),S(i)*mx(i),S(i)*my(i),S(i)*Bi(i),cortab(S(i),cor(i))
   end do

   return
end subroutine config

subroutine config_ising(u_conf)
   use var_global, only: Ns,S
   implicit none
   integer,intent(in) :: u_conf
   integer :: i

   do i = 1,Ns
      write(u_conf,*) (S(i)+1)/2
   end do

   return
end subroutine config_ising

subroutine mag(ibx,iby,MB)
   use var_global
   implicit none
   real(8), intent(in) :: ibx,iby,MB
   real(8) :: Magx,Magy

   Magx = sum(S(:)*mx(:))
   Magy = sum(S(:)*my(:))
   write(30,*) ibx,Magx/real(Ns,8),iby,Magy/real(Ns,8),MB,(Magx*cos(theta_c*pi/180.0d0) + Magy*sin(theta_c*pi/180.0d0))/real(Ns,8)

   return
end subroutine mag

subroutine mssf(unit_in,flag)
   use var_global, only: Ns,x,y,mx,my,S,Nmssf,dir2
   implicit none
   integer, intent(in) :: unit_in,flag
   integer :: i

   if (flag==1) then
      open(unit = unit_in,file=trim(dir2) // trim("mssf.dat"))
      write(unit_in,*) Ns,Nmssf
      do i = 1,Ns
         write(unit_in,*) x(i),y(i),mx(i),my(i)
      end do
   else if (flag==2) then
      do i = 1,Ns
         write(unit_in,*) (S(i)+1)/2
      end do
   else
      close(unit_in)
   end if

   return
end subroutine mssf

subroutine one_rev
   use var_global, only : Hi,Bx,By,mx,my,N_campo
   implicit none
   integer :: i,try

   do i = 1,N_campo
      try = 1
      Hi(:) = Bx(i)*mx(:) + By(i)*my(:)
      do while (try>0)
         call budrikis(try)
         !call config(40)
      end do
   end do

end subroutine one_rev

subroutine simu_rev
   use var_global, only : N_campo,Hi,Bx,By,mx,my,ModB,y,Nmssf
   implicit none
   integer :: i,try,N
   real(8) :: dM

   dM = 360.0d0/real(N_campo)
   N = N_campo/Nmssf

   do i = 1,N_campo
      try = 1
      Hi(:) = Bx(i)*mx(:) + By(i)*my(:)
      
      do while (try>0)
         call budrikis(try)
         call config_ising(100)
         !call config(40)
      end do

      call config(20)
      call mag(Bx(i),By(i),ModB(i))

      !if (mod(i,N)==0) then
         call mssf(50,2)
      !   j = j + 1
      !end if

      write(40,*) 1
      write(40,*) ModB(i)
      write(40,*) 0.,maxval(y)+1.5d0,Bx(i)/abs(ModB(i)),By(i)/abs(ModB(i))

      call flush()
   end do

end subroutine simu_rev

program main
   implicit none
   integer :: i

   call inicial
   call abrir_fechar(1,20,30,40)
   call mssf(50,1)

   ! Primeira revolução do Campo - Não salva resultados !
   call one_rev

   do i = 1,2
      call simu_rev
      print*, i
   end do

   call abrir_fechar(0,20,30,40)
   call mssf(50,3)

end program main
