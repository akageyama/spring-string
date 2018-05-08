
module basic
  implicit none
  integer, parameter :: dp = kind(1.d0)
  real(dp), parameter :: pi = 3.14159265358979323846_dp
  real(dp), parameter :: twopi = 2*pi
  integer, parameter :: imax = 501     !  チューブの長さ方向。
  integer, parameter :: jmax = 3        !  断面の正多角形の頂点数
  !
  !      /\
  !   a /  \ a
  !    /    \
  !    ------
  !      a: 正三角形の一辺の長さ    
  real(dp), parameter :: a = 1.0_dp          ! 長さの規格化単位
  real(dp), parameter :: sc = 1.0_dp         ! バネ定数
  real(dp), parameter :: mass = 1.0_dp       ! 質点の質量
  real(dp), parameter :: resis = 0.01_dp     ! 速度に比例する摩擦係数
  real(dp), parameter :: edge_tension = 1.e-4_dp   ! 両端の張力
                                      ! 1.e-3は強すぎた (imax=31)
  real(dp) :: v_twist                        ! 両端の質点の(回転)速度
  real(dp) :: time                           ! 時刻
  integer :: nloop                           ! シミュレーションループ数
  real(dp) :: tau                            ! 特性時間。
  real(dp) :: dt                             ! 時間刻み
  real(dp), parameter :: dt_factor = 0.001_dp ! 積分ファクター。充分小さい数
  ! チューブの両端を次の角速度で(互いに反対方向に)捻る。
  real(dp) :: dphidt
  !
  type tube_ 
    ! チューブの各質点の位置 (x,y,zの３成分)
    real(dp), dimension(imax,jmax,3) :: pos
    ! チューブの各質点の速度 (x,y,zの３成分)
    real(dp), dimension(imax,jmax,3) :: vel
    ! チューブの両端のjmax-個の質点の初期位置(角度)
    real(dp), dimension(jmax) :: phi_btm_init, phi_top_init
    ! チューブの両端の初期高さ
    real(dp) :: z_btm_init, z_top_init
    ! チューブの半径。この半径の円上に正多角形(jmax-角形)の頂点が乗る。
    ! 各頂点間の距離が丁度 a になるようにチューブの半径をとる。
    real(dp) :: rad 
    ! バネのベクトル。 (i,j)番目の質点に向かうベクトル
    real(dp), dimension(imax-1,jmax,3,3) :: bond_v  ! bond vector
    !                      左の3は東、東北、北西方向のバネに対応
    !                      右の3はベクトルの(x,y,z)成分
    !
    !         (i+1,j)     (i+1,j+1)
    !             o         o
    !              \       /
    !               \3   2/
    !                \   /  1
    !  (i,j-1)o------- o -------o (i,j+1)              
    !                (i,j)
    !                /   \
    !               /     \                          z方向
    !              /       \                           |
    !             o         o                          |
    !          (i-1,j)    (i-1,j+1)                    +---> x,y方向
  end type tube_
!
  type(tube_) :: tube
!
end module basic

module init
  use basic
  implicit none
  private
  public :: init_param, init_tube, init_check, init_perturb
contains
!
  subroutine init_check
!    if(jmax < 6) then
!      print *, '*** ERROR.  jmax must >= 6. '
!      stop
!    end if
  end subroutine init_check
!  
  subroutine init_param
    tube%rad = a / sqrt(2*(1.0_dp-cos(twopi/jmax)))    ! チューブの半径
    tau = sqrt(mass/sc) * imax      ! 波がチューブを伝わる時間(概算)
    dphidt = (twopi/tau)/ 60        ! tau 時間に 1/60 回転
    v_twist = tube%rad * dphidt     ! 捻りの回転速度
    dt = min(0.001_dp, dt_factor/v_twist)          ! (初期の)時間刻み
  end subroutine init_param
!  
  subroutine init_tube
    call position
    call velocity
    tube%z_btm_init = tube%pos(   1,1,3)
    tube%z_top_init = tube%pos(imax,1,3)
  end subroutine init_tube
!
  subroutine position
    integer :: i, j
    real(dp), parameter :: dphi = twopi / jmax
    real(dp), parameter :: dphi2 = dphi / 2
    real(dp), parameter :: dphi4 = dphi / 4
    real(dp), dimension(jmax) :: phi
    real(dp) :: dz    ! 正三角形の高さ(少し傾いていることに注意)

    dz = a * sqrt(1.0_dp - (sin(dphi4)/sin(dphi2))**2 )

    do i = 1 , imax
      do j = 1 , jmax
        if ( mod(i-1,2) == 0 ) then
          phi(j) = dphi * (j-1)
        else
          phi(j) = dphi * (j-1) + dphi2
        end if
      end do
      ! 両端の各質点の初期位置を記憶しておく
      if (i==1) then
        tube%phi_btm_init(:) = phi(:)
      else if (i==imax) then
        tube%phi_top_init(:) = phi(:)
      end if
      ! 初期位置の設定
      tube%pos(i,:,1) = tube%rad * cos(phi(:))
      tube%pos(i,:,2) = tube%rad * sin(phi(:))
      tube%pos(i,:,3) = dz * (i-1)
    end do
  end subroutine position
!
  subroutine velocity
    ! 初速度はゼロ
    tube%vel(:,:,:) = 0.0_dp
  end subroutine velocity
!
  subroutine init_perturb
    real(dp) :: amp        ! amplitude
    integer :: j, m        ! m=mode
    real(dp) :: phi, ran, z_l
    real(dp), dimension(imax) :: dx, dy, z, phase

    amp = a * 0.001_dp

    z = tube%pos(:,1,3) - tube%pos(1,1,3)
    z_l = z(imax)

    do m = 1 , 10
      call random_number(ran)
      phi = twopi * ran
      print *, ' - perturbation: amp=', amp,' m=',m,' phi=',phi
      phase = m*z*pi/z_l
      dx = amp * cos(phi) * sin(phase)
      dy = amp * sin(phi) * sin(phase)
      do j = 1 , jmax
        tube%pos(:,j,1) = tube%pos(:,j,1) + dx
        tube%pos(:,j,2) = tube%pos(:,j,2) + dy
      end do
    end do
  end subroutine init_perturb
end module init


module spring
  use basic
  implicit none
  private
  public :: spring_force, spring_bonds
contains
!
  subroutine spring_force(pos,force)
    real(dp), dimension(imax,jmax,3), intent(in) :: pos
    real(dp), dimension(imax,jmax,3), intent(out) :: force
    call spring_bonds(pos)
    call bulk(force)
    call boundary_with_tension(force)
  end subroutine spring_force
!
  function f(b)
    ! 二つの質点間に働くバネ力
    real(dp), dimension(3) :: f
    real(dp), dimension(3), intent(in) :: b
    real(dp), dimension(3) :: unit_b
    real(dp) :: babs, dx

    babs = sqrt(dot_product(b,b))
    unit_b = b / babs
    dx = a - babs        ! バネの縮み
    f = unit_b * dx
  end function f
!
  !         (i+1,j)     (i+1,j+1)
  !             o         o
  !              \       /
  !               \3   2/
  !                \   / bindex=1
  !  (i,j-1)o------- o -------o (i,j+1)
  !                (i,j)
  !                /   \
  !               /     \
  !              /       \
  !             o         o
  !          (i-1,j)    (i-1,j+1)
!
  subroutine spring_bonds(pos)
    real(dp), dimension(imax,jmax,3), intent(in) :: pos
    ! 全てのバネの伸びと方向
    !
    !               
    !        ------j=3-----4------5------6-------    i=4
    !                               \   /
    !        --j=3------4------5------6-----7----    i=3
    !                    \   /
    !        ------j=3-----4------5------6-------    i=2
    !                     /  \
    !        --j=3------4------5------6-----7----    i=1
    !
    integer :: i, j, jp1, jm1
    do i = 1 , imax-1 , 2        ! 奇数
      do j = 1 , jmax
        jp1 = j+1; if(j==jmax) jp1 = 1
        jm1 = j-1; if(j==1)    jm1 = jmax
        tube%bond_v(i,j,1,:) = pos(i,j,:) - pos(i,  jp1,:)
        tube%bond_v(i,j,2,:) = pos(i,j,:) - pos(i+1,j,  :)
        tube%bond_v(i,j,3,:) = pos(i,j,:) - pos(i+1,jm1,:)
      end do
    end do
    do i = 2 , imax-1 , 2        ! 偶数
      do j = 1 , jmax
        jp1 = j+1; if(j==jmax) jp1 = 1
        tube%bond_v(i,j,1,:) = pos(i,j,:) - pos(i,  jp1,:)
        tube%bond_v(i,j,2,:) = pos(i,j,:) - pos(i+1,jp1,:)
        tube%bond_v(i,j,3,:) = pos(i,j,:) - pos(i+1,j,  :)
      end do
    end do
  end subroutine spring_bonds
!
  subroutine bulk(force)
    real(dp), dimension(imax,jmax,3), intent(out) :: force
    ! チューブの各質点に掛かる力の計算
    integer :: i, j, jp1, jm1
    real(dp), dimension(3) :: f1, f2, f3
    
    do i = 2 , imax-1 , 2        ! 偶数
      do j = 1 , jmax
        jp1 = j+1; if(j==jmax) jp1 = 1
        jm1 = j-1; if(j==1)    jm1 = jmax
        f1 = sc*(f(tube%bond_v(i,j,1,:))-f(tube%bond_v(  i,jm1,1,:)))
        f2 = sc*(f(tube%bond_v(i,j,2,:))-f(tube%bond_v(i-1,  j,2,:)))
        f3 = sc*(f(tube%bond_v(i,j,3,:))-f(tube%bond_v(i-1,jp1,3,:)))
        force(i,j,:) = f1 + f2 + f3
      end do                               
    end do

    do i = 3 , imax-1 , 2        ! 奇数
      do j = 1 , jmax
        jm1 = j-1; if(j==1) jm1 = jmax
        f1 = sc*(f(tube%bond_v(i,j,1,:))-f(tube%bond_v(  i,jm1,1,:)))
        f2 = sc*(f(tube%bond_v(i,j,2,:))-f(tube%bond_v(i-1,jm1,2,:)))
        f3 = sc*(f(tube%bond_v(i,j,3,:))-f(tube%bond_v(i-1,  j,3,:)))
        force(i,j,:) = f1 + f2 + f3
      end do                               
    end do
  end subroutine bulk
!
  subroutine boundary_with_tension(force)
    real(dp), dimension(imax,jmax,3), intent(out) :: force

    integer :: i, j, jp1, jm1
    ! 「底」境界の質点にかかる力
    i = 1
    do j = 1 , jmax
      force(i,j,:) = sc * ( f(tube%bond_v( i, j, 2, :))   & 
                          + f(tube%bond_v( i, j, 3, :)) ) 
    end do
    force(i,:,1) = 0.0_dp
    force(i,:,2) = 0.0_dp
    force(i,:,3) = sum(force(i,:,3)) / jmax - edge_tension

    ! 「天井」境界の質点にかかる力
    i = imax
    if (mod(i,2)==0) then
      do j = 1 , jmax
        jp1 = j+1; if(j==jmax) jp1 = 1
        force(i,j,:) = sc * ( - f(tube%bond_v(i-1,   j, 2, :))      &
                              - f(tube%bond_v(i-1, jp1, 3, :))  )
      end do                               
    else
      do j = 1 , jmax
        jm1 = j-1; if(j==1) jm1 = jmax
        force(i,j,:) = sc * ( - f(tube%bond_v(i-1, jm1, 2, :))      &
                              - f(tube%bond_v(i-1,   j, 3, :))  )
      end do                               
    end if
    force(i,:,1) = 0.0_dp
    force(i,:,2) = 0.0_dp
    force(i,:,3) = sum(force(i,:,3)) / jmax + edge_tension
  end subroutine boundary_with_tension
end module spring

module self
  use basic
  implicit none
  private
  public :: self_interact, self_set_centers, self_near_table
  integer, parameter :: imax2 = imax/4                ! 概算
  integer, parameter :: imax3 = imax + imax2
  real(dp), dimension(-imax2:imax3,3) :: centers
  real(dp), dimension(2:imax-1,-imax2:imax3,3) :: repuls_f
  real(dp), parameter :: r0 = 2*a           ! ぶつかり始める距離
  real(dp), parameter :: rfar = 20*a        ! 充分遠い距離
  logical, dimension(imax,-imax2:imax3) :: isnear ! 近いか？
contains
  subroutine self_near_table
    integer :: i1, i2
    
    isnear = .false.                  ! 初期化

    ! チューブ自身
    do i1 = 1 , imax
      do i2 = i1 , imax
        isnear(i1,i2) = check(centers(i1,:), centers(i2,:))
      end do
    end do
 
    do i1 = 1 , imax
      do i2 = 1 , i1-1
        isnear(i1,i2) = isnear(i2,i1)
      end do
    end do
    
    do i1 = 1 , imax
      ! 下端よりも下にある仮想的な棒との距離
      do i2 = -imax2 , 0
        isnear(i1,i2) = check(centers(i1,:), centers(i2,:))
      end do
      ! 上端よりも上にある仮想的な棒との距離
      do i2 = imax+1 , imax3
        isnear(i1,i2) = check(centers(i1,:), centers(i2,:))
      end do
    end do
  end subroutine self_near_table
!
  function check(pos1, pos2)
    logical :: check
    real(dp), intent(in), dimension(3) :: pos1, pos2
    real(dp) :: dx, dy, dz
    dx = abs(pos1(1)-pos2(1))
    dy = abs(pos1(2)-pos2(2))
    dz = abs(pos1(3)-pos2(3))
    if (dx <= rfar .and. dy <= rfar .and. dz <= rfar) then
      check = .true.
    else
      check = .false.
    end if
  end function check
!
  subroutine self_interact(force)
    real(dp), dimension(imax,jmax,3), intent(out) :: force

    real(dp) :: total_repuls_fx, total_repuls_fy, total_repuls_fz
    integer :: i
    real(dp), dimension(2:imax-1) :: r

    do i = 2 , imax-1
      r(i) = sqrt(dot_product(centers(i,1:2),centers(i,1:2)))
    end do
 
    if (maxval(r) < 2*a) return      ! 未だまっすぐ

    call repulsions

    do i = 2 , imax-1
      total_repuls_fx = sum(repuls_f(i,:,1))
      total_repuls_fy = sum(repuls_f(i,:,2))
      total_repuls_fz = sum(repuls_f(i,:,3))
      force(i,:,1) = force(i,:,1) + total_repuls_fx
      force(i,:,2) = force(i,:,2) + total_repuls_fy
      force(i,:,3) = force(i,:,3) + total_repuls_fz
    end do
  end subroutine self_interact
!
  subroutine self_set_centers(pos)
    real(dp), dimension(imax,jmax,3), intent(in) :: pos

    integer :: i

    do i = 1 , imax  
      centers(i,1) = sum(pos(i,:,1)) / jmax     ! x
      centers(i,2) = sum(pos(i,:,2)) / jmax     ! y
      centers(i,3) = sum(pos(i,:,3)) / jmax     ! z
    end do

    do i = 0 , -imax2 , -1           ! 底境界点よりも下にある仮想的な棒
      centers(i,1) = 0.0_dp
      centers(i,2) = 0.0_dp
      centers(i,3) = centers(1,3) + a*(i-1)
    end do

    do i = imax+1 , imax3            ! 天井境界点より上にある仮想的な棒
      centers(i,1) = 0.0_dp
      centers(i,2) = 0.0_dp
      centers(i,3) = centers(imax,3) + a*(i-imax)
    end do
  end subroutine self_set_centers
!
  subroutine repulsions
    real(dp) :: z_top, z_btm

    repuls_f(:,:,:) = 0.0_dp             ! 初期化

    call itself
    call symmetry                        ! 作用-反作用の法則
    
    z_top = centers(imax,3)
    z_btm = centers(   1,3)

    ! 輪くぐり現象の禁止
    if ( maxval(centers(1:imax-1,3)) > z_top ) call edge('top')
    if ( minval(centers(2:imax,  3)) < z_btm ) call edge('btm')
  end subroutine repulsions
!
  subroutine itself
    ! チューブ自身の相互作用
    integer, parameter :: iskip = 5                ! 近所の定義
    integer :: i1, i2

    do i1 = 2 , imax-1
      do i2 = i1+1 , imax-1
        if (i2-i1 <= iskip) cycle                  ! 御近所どうしは反発しない
        if (isnear(i1,i2)) then
          repuls_f(i1,i2,:) = force(centers(i1,:), centers(i2,:))
        end if
      end do
    end do
  end subroutine itself
!
  subroutine symmetry
    integer :: i1, i2
    do i1 = 2 , imax-1
      do i2 = 2 , i1-1
        if (isnear(i1,i2)) then
          repuls_f(i1,i2,:) = -repuls_f(i2,i1,:)      ! 反作用
        end if
      end do
    end do
  end subroutine symmetry
!
  subroutine edge(which)
    character(len=3), intent(in) :: which

    integer :: i1, i2

    if (which=='btm') then
      ! 底境界が輪くぐりをしないように
      do i1 = 2 , imax-1
        do i2 = -imax2 , 0
          if (isnear(i1,i2)) then
            repuls_f(i1,i2,:) = force(centers(i1,:), centers(i2,:))
          end if
        end do
      end do
    end if
 
    if (which=='top') then
      ! 天井境界が輪くぐりをしないように
      do i1 = 2 , imax-1
        do i2 = imax+1 , imax3
          if (isnear(i1,i2)) then
            repuls_f(i1,i2,:) = force(centers(i1,:), centers(i2,:))
          end if
        end do
      end do
    end if
  end subroutine edge
!
  function force(pos1, pos2)
    real(dp), dimension(3) :: force
    real(dp), dimension(3), intent(in) :: pos1, pos2

    real(dp), dimension(3) :: vect, unit_v
    real(dp) :: dist

    vect = pos1 - pos2                         ! 相手から自分への向き 
    dist = sqrt(dot_product(vect,vect))        ! 距離
    unit_v = vect / dist                       ! 単位ベクトル
    force = unit_v * amp(dist)
  end function force
!
  function amp(dist)
    real(dp) :: amp
    real(dp), intent(in) :: dist

    amp = (r0/dist)**13
  end function amp
end module self

module motion
  use basic; use spring; use self
  implicit none
  private
  public :: motion_advance
  real(dp), dimension(imax,jmax,3), save :: dvdt, dpdt
  real(dp), dimension(imax,jmax,3) :: vel2, pos2
  real(dp), dimension(imax,jmax,3) :: force
contains
!
  subroutine motion_advance
    ! 方程式：   dv/dt = ( force  - resis * mass * v ) / mass
    !            dx/dt = v
    ! 手法：     2次ルンゲ＝クッタ法  (今、積分精度はそれほど高い必要はない)
    logical, save :: first = .true.

    if(first) then
      first = .false.
    else
      call reset_dt
    end if

    call step_1
    call step_2
  end subroutine motion_advance
!
  subroutine step_1
    integer, parameter :: nloop_ref = 2.0_dp / dt_factor
    real(dp) :: dt2, resismass
 
    dt2 = dt / 2
    resismass = resis * mass

    ! -- RK, 1st step --
    call spring_force(tube%pos, force)
    call self_set_centers(tube%pos)
    if (mod(nloop,nloop_ref)==0) call self_near_table
    call self_interact(force)

    dvdt(:,:,:) = ( force(:,:,:) - resismass * tube%vel(:,:,:) ) / mass
    dpdt(:,:,:) = tube%vel(:,:,:)

    vel2(:,:,:) = tube%vel(:,:,:) + dt2 * dvdt(:,:,:)
    pos2(:,:,:) = tube%pos(:,:,:) + dt2 * dpdt(:,:,:)

    call twist(time+dt/2, pos2)
  end subroutine step_1
!
  subroutine step_2
    ! -- RK, 2nd step --
    real(dp) :: resismass
    resismass = resis * mass
    call spring_force(pos2, force)
    call self_set_centers(pos2)
    call self_interact(force)

    dvdt(:,:,:) = ( force(:,:,:) - resismass * vel2(:,:,:) ) / mass
    dpdt(:,:,:) = vel2(:,:,:)

    tube%vel(:,:,:) = tube%vel(:,:,:) + dt * dvdt(:,:,:)
    tube%pos(:,:,:) = tube%pos(:,:,:) + dt * dpdt(:,:,:)

    call twist(time+dt, tube%pos)
  end subroutine step_2
!
  subroutine twist(t, pos)
    real(dp) :: t              ! 時刻
    real(dp), dimension(imax,jmax,3), intent(out) :: pos
    ! 両端で反対方向に半分ずつ捻る
    real(dp) :: theta 
    real(dp), dimension(jmax) :: phi_btm, phi_top
    theta = dphidt * t
    phi_btm(:) = tube%phi_btm_init(:) + theta
    phi_top(:) = tube%phi_top_init(:) - theta
    pos(   1,:,1) = tube%rad * cos(phi_btm(:))
    pos(   1,:,2) = tube%rad * sin(phi_btm(:))
    pos(imax,:,1) = tube%rad * cos(phi_top(:))
    pos(imax,:,2) = tube%rad * sin(phi_top(:))
  end subroutine twist
!
  subroutine reset_dt
    real(dp) :: ff_max, vv_max
    ff_max = maxval(abs(dpdt))
    vv_max = maxval(abs(dvdt))
    vv_max = max(vv_max,v_twist)
    dt = min(dt_factor/ff_max, dt_factor/vv_max)
  end subroutine reset_dt
end module motion


module diagno
  use basic; use spring
  implicit none
  private
  public :: diagno_first, diagno_memo, diagno_position
  integer, parameter :: ids_memo = 10
  integer, parameter :: ids_xyz1 = 11
  integer, parameter :: ids_xyz2 = 12
  integer, parameter :: ids_xyz3 = 13
contains
!
  subroutine diagno_first
    print *, '--------------------------------------------'
    print *, '   imax = ', imax
    print *, '   jmax = ', jmax
    print *, '   mass = ', mass
    print *, '      a = ', a
    print *, '  resis = ', resis
    print *, '    tau = ', tau
    print *, '     dt = ', dt
    print *, ' dphidt = ', dphidt
    print *, '    rad = ', tube%rad
    print *, '  z_btm = ', tube%z_btm_init
    print *, '  z_top = ', tube%z_top_init
    print *, '--------------------------------------------'
    write(ids_memo,*) '--------------------------------------------'
    write(ids_memo,*) '   imax = ', imax
    write(ids_memo,*) '   jmax = ', jmax
    write(ids_memo,*) '   mass = ', mass
    write(ids_memo,*) '      a = ', a
    write(ids_memo,*) '  resis = ', resis
    write(ids_memo,*) '    tau = ', tau
    write(ids_memo,*) '     dt = ', dt
    write(ids_memo,*) ' dphidt = ', dphidt
    write(ids_memo,*) '    rad = ', tube%rad
    write(ids_memo,*) '  z_btm = ', tube%z_btm_init
    write(ids_memo,*) '  z_top = ', tube%z_top_init
    write(ids_memo,*) '--------------------------------------------'
  end subroutine diagno_first
!  
  subroutine diagno_memo
    real(dp) :: ke  ! kinetic energy
    real(dp) :: pe  ! potential energy

    ke = 0.5_dp * mass * sum(tube%vel(:,:,:)**2)
    if (nloop==0) call spring_bonds(tube%pos)
    pe = potential()
    
    print *, nloop, time, dt, dphidt*time/twopi*2, ke, pe
    write(ids_memo,*) nloop, time, dt, dphidt*time/twopi*2, ke, pe
  end subroutine diagno_memo
!
  function potential()
    real(dp) :: potential
    integer :: i, j, n
    real(dp) :: length

    potential = 0.0_dp
    do i = 1 , imax-1
      do j = 1 , jmax
        do n = 1 , 3        ! バネの3つの方向
          length = sqrt(dot_product(tube%bond_v(i,j,n,:),    &
                                tube%bond_v(i,j,n,:)))
          potential = potential + (length-a)**2
        end do
      end do
    end do
    potential = 0.5_dp * sc * potential
  end function potential
!
  subroutine diagno_position
    integer, save :: counter = 0
    character(len=4) :: ctr

    counter = counter + 1

    if (counter>=5000) then
      print *,'*** ERROR: too many output files.'
      stop
    end if

    ctr = i2c(counter)

    call printline(1,ids_xyz1)
    call printline(2,ids_xyz2)
    call printline(3,ids_xyz3)
 
    print *,' position data counter = ', ctr
  end subroutine diagno_position
!
  subroutine printline(j,ids)
    integer, intent(in) :: j, ids
    integer :: i

    write(ids,*)'#  '
    write(ids,*)'# nloop = ',nloop, '  time = ', time
    do i = 1 , imax 
      if(mod(i,2)==1) then
        write(ids,*) tube%pos(i,j,1), tube%pos(i,j,2), tube%pos(i,j,3)
      else
        write(ids,*) '## ',tube%pos(i,j,1), tube%pos(i,j,2), tube%pos(i,j,3)
      end if
    end do
  end subroutine printline
!
  function i2c(i)
    integer, intent(in) :: i
    character(len=4) :: i2c
    write(i2c,fmt='(i4)') i
    if(i2c(1:1)==' ') i2c(1:1) = '0'
    if(i2c(2:2)==' ') i2c(2:2) = '0'
    if(i2c(3:3)==' ') i2c(3:3) = '0'
    if(i2c(4:4)==' ') i2c(4:4) = '0'
  end function i2c
end module diagno


program twisting_tube
  use basic; use init; use motion; use spring; use diagno; use self
  implicit none

  real(dp) :: time_prev_diag
  real(dp) :: diag_factor = 1.0_dp/20
  integer :: ctr = 0
  
  time = 0.0_dp
  nloop = 0
  time_prev_diag = 0.0_dp

  call init_check
  call init_param
  call init_tube
  call diagno_first
  call diagno_memo
  call diagno_position
  call init_perturb

! do while ( nloop < 10000000 )
  do while ( nloop < 200000 )
    call motion_advance
    time = time + dt
    nloop = nloop + 1
    if( (time-time_prev_diag) > tau*diag_factor ) then
      ctr = ctr + 1
      if(ctr==10) then
        call diagno_position
        ctr = 0
      end if
      call diagno_memo
      time_prev_diag = time
    end if
  end do
end program twisting_tube
