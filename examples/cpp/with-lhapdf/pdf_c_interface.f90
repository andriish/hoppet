module pdf_c_interface
  ! Canonical LHAPDF C interop layer for the example programs.
  use iso_c_binding, only: c_char, c_int, c_double, c_null_char, c_intptr_t, c_ptr, c_null_ptr, c_loc, c_f_pointer
  use types, only: dp
  implicit none
  private

  !-------------------------------------------------------------------
  ! Derived type returned by the struct-based locator interface.
  !-------------------------------------------------------------------
  type, bind(C), public :: pdf_location
    character(kind=c_char) :: pdfsets_path(1024)
    character(kind=c_char) :: pdfdata_path(1024)
    character(kind=c_char) :: pdfinfo_path(1024)
    integer(c_int) :: valid
  end type pdf_location

  !> C-compatible container for PDF data points. The array data lives in
  !> module-level target storage; this struct only carries opaque pointer
  !> slots (plus sizes) so it can be BIND(C).
  type, bind(C) :: pdfdat_t
    integer(c_int) :: n = 0
    integer(c_int) :: nQ = 0
    integer(c_int) :: nFlav = 0
    type(c_ptr) :: x_ptr = c_null_ptr
    type(c_ptr) :: q_ptr = c_null_ptr
    type(c_ptr) :: f_ptr = c_null_ptr
    type(c_ptr) :: flav_ptr = c_null_ptr
  end type pdfdat_t

  !-------------------------------------------------------------------
  ! C-compatible container for selected PDF metadata fields. This type is
  ! intentionally kept as a plain BIND(C) struct with opaque pointer slots for
  ! the dynamically allocated arrays. The C/C++ side only preserves the layout
  ! and does not dereference these slots directly.
  !-------------------------------------------------------------------
  type, bind(C) :: pdf_meta_t
    character(kind=c_char) :: SetDesc(1024) = ''
    integer(c_int) :: SetIndex = 0
    character(kind=c_char) :: Authors(1024) = ''
    character(kind=c_char) :: Reference(1024) = ''
    character(kind=c_char) :: Format(1024) = ''
    character(kind=c_char) :: DataVersion(1024) = ''
    integer(c_int) :: NumMembers = 0
    integer(c_int) :: Particle = 0
    integer(c_intptr_t) :: Flavors = 0_c_intptr_t
    integer(c_int) :: OrderQCD = 0
    character(kind=c_char) :: FlavorScheme(1024) = ''
    integer(c_int) :: NumFlavors = 0
    character(kind=c_char) :: ErrorType(1024) = ''
    real(c_double) :: XMin = 0.0_c_double
    real(c_double) :: XMax = 0.0_c_double
    real(c_double) :: QMin = 0.0_c_double
    real(c_double) :: QMax = 0.0_c_double
    real(c_double) :: MZ = 0.0_c_double
    real(c_double) :: MUp = 0.0_c_double
    real(c_double) :: MDown = 0.0_c_double
    real(c_double) :: MStrange = 0.0_c_double
    real(c_double) :: MCharm = 0.0_c_double
    real(c_double) :: MBottom = 0.0_c_double
    real(c_double) :: MTop = 0.0_c_double
    real(c_double) :: AlphaS_MZ = 0.0_c_double
    integer(c_int) :: AlphaS_OrderQCD = 0
    character(kind=c_char) :: AlphaS_Type(1024) = ''
    integer(c_intptr_t) :: AlphaS_Qs = 0_c_intptr_t
    integer(c_intptr_t) :: AlphaS_Vals = 0_c_intptr_t
  end type pdf_meta_t

  public :: pdf_locate_c, pdf_get_values_c, pdf_read_grid_c, read_pdfdat, get_pdfdat_point
  public :: pdfdat_t, pdf_meta_t, pdf_read_metadata_c

  real(c_double), allocatable, target, save :: store_x(:)
  real(c_double), allocatable, target, save :: store_q(:)
  real(c_double), allocatable, target, save :: store_f(:,:,:)
  integer(c_int), allocatable, target, save :: store_flav(:)
  real(c_double), allocatable, target, save :: store_alpha_qs(:)
  real(c_double), allocatable, target, save :: store_alpha_vals(:)
  integer(c_int), save :: n_alpha_points = 0

contains

  subroutine copy_c_string(dest, src)
    character(kind=c_char), intent(out) :: dest(:)
    character(len=*), intent(in) :: src
    integer :: i, n, sz

    sz = size(dest)
    do i = 1, sz
      dest(i) = c_null_char
    end do
    n = min(len_trim(src), sz-1)
    do i = 1, n
      dest(i) = src(i:i)
    end do
    if (sz > 0) dest(min(n+1, sz)) = c_null_char
  end subroutine copy_c_string

  function pdf_locate_c(pdfname, pdfname_len) result(loc) bind(C)
    character(kind=c_char), intent(in) :: pdfname(*)
    integer(c_int), value :: pdfname_len
    type(pdf_location) :: loc
    character(len=:), allocatable :: f_pdfname
    character(len=:), allocatable :: info_path, dat_path, index_file
    integer :: ios, unit, i, n_dirs, pos, start, slash_pos
    integer :: dummy_id, dir_idx, err_code
    character(len=512) :: line, set, path, env, token
    logical :: found
    character(len=256) :: dirs(100)
    character(len=512) :: cand

    f_pdfname = ''
    do i = 1, pdfname_len
      if (pdfname(i) == c_null_char) exit
      f_pdfname = f_pdfname // pdfname(i)
    end do

    err_code = 0
    call get_environment_variable('LHAPDF_DATA_PATH', env, status=ios)
    n_dirs = 0
    if (ios == 0 .and. len_trim(env) > 0) then
      start = 1
      do
        pos = index(env(start:), ':')
        if (pos == 0) then
          token = env(start:)
          if (len_trim(token) > 0) then
            n_dirs = n_dirs + 1
            dirs(n_dirs) = trim(token)
          end if
          exit
        else
          token = env(start:start+pos-2)
          if (len_trim(token) > 0) then
            n_dirs = n_dirs + 1
            dirs(n_dirs) = trim(token)
          end if
          start = start + pos
        end if
      end do
    end if
    if (n_dirs == 0) then
      n_dirs = 2
      dirs(1) = '/usr/share/LHAPDF/'
      dirs(2) = '/usr/local/share/LHAPDF/'
    end if
    print *, '[DEBUG] PDF locator will search', n_dirs, 'directories:'
    do i = 1, n_dirs
      print *, '  ', trim(dirs(i))
    end do

    index_file = ''
    do i = 1, n_dirs
      cand = trim(dirs(i)) // '/pdfsets.index'
      inquire(file=cand, exist=found)
      if (found) then
        index_file = cand
        exit
      end if
    end do
    if (len_trim(index_file) == 0) then
      err_code = 1
    end if
    if (err_code == 0) then
      print *, '[DEBUG] Found pdfsets.index at', trim(index_file)
      open(newunit=unit, file=index_file, status='old', action='read', iostat=ios)
      if (ios /= 0) then
        err_code = 2
      else
        found = .false.
        do while (.true.)
          read(unit, '(A)', iostat=ios) line
          if (ios /= 0) exit
          if (len_trim(line) == 0) cycle
          if (line(1:1) == '#') cycle
          read(line, *, iostat=ios) dummy_id, set, dir_idx
          if (ios /= 0) then
            print *, '[DEBUG] failed to parse line, skipping'
            cycle
          end if
          if (dir_idx >= 1 .and. dir_idx <= n_dirs) then
            path = trim(dirs(dir_idx))
          else
            path = ''
          end if
          if (trim(set) == trim(f_pdfname)) then
            found = .true.
            print *, '[DEBUG] MATCH FOUND for set', trim(f_pdfname)
            exit
          end if
        end do
        close(unit)
        if (.not. found) err_code = 3
      end if
    end if

    if (err_code == 0) then
      info_path = trim(path) // '/' // trim(f_pdfname) // '/' // trim(f_pdfname) // '.info'
      dat_path  = trim(path) // '/' // trim(f_pdfname) // '/' // trim(f_pdfname) // '_0000.dat'
      print *, '[DEBUG] Resolved info_path =', trim(info_path)
      print *, '[DEBUG] Resolved dat_path  =', trim(dat_path)

      call copy_c_string(loc%pdfinfo_path, info_path)
      call copy_c_string(loc%pdfdata_path, dat_path)
      call copy_c_string(loc%pdfsets_path, '')
      if (len_trim(info_path) > 0) then
        slash_pos = index(info_path, '/', .true.)
        if (slash_pos > 0) then
          call copy_c_string(loc%pdfsets_path, info_path(1:slash_pos-1))
        end if
      end if
      loc%valid = 0
    else
      call copy_c_string(loc%pdfinfo_path, '')
      call copy_c_string(loc%pdfdata_path, '')
      call copy_c_string(loc%pdfsets_path, '')
      loc%valid = err_code
    end if
  end function pdf_locate_c

  function pdf_read_metadata_c(fname, flen) result(meta) bind(C)
    character(kind=c_char), intent(in) :: fname(*)
    integer(c_int), value :: flen
    type(pdf_meta_t) :: meta
    character(len=:), allocatable :: f_name, key, val
    integer :: i, ierr, ios, unit, pos, nentry
    character(len=4096) :: line, keybuf, valbuf

    f_name = ''
    do i = 1, flen
      if (fname(i) == c_null_char) exit
      f_name = f_name // fname(i)
    end do

    meta = pdf_meta_t()
    if (allocated(store_alpha_qs)) deallocate(store_alpha_qs)
    if (allocated(store_alpha_vals)) deallocate(store_alpha_vals)
    n_alpha_points = 0
    open(newunit=unit, file=trim(f_name), status='old', action='read', iostat=ios)
    if (ios /= 0) return

    nentry = 0
    do
      read(unit, '(A)', iostat=ios) line
      if (ios /= 0) exit
      if (len_trim(line) == 0) cycle

      keybuf = ''
      valbuf = ''
      pos = index(line, ':')
      if (pos == 0) cycle

      keybuf = line(:pos-1)
      valbuf = line(pos+1:)
      keybuf = trim(adjustl(keybuf))
      valbuf = trim(adjustl(valbuf))
      if (len_trim(keybuf) == 0) cycle

      nentry = nentry + 1
      key = trim(keybuf)
      val = trim(valbuf)
      if ((trim(key) == 'AlphaS_Qs' .or. trim(key) == 'AlphaS_Vals') .and. index(trim(val), '[') /= 0 .and. index(trim(val), ']') == 0) then
        do
          read(unit, '(A)', iostat=ios) line
          if (ios /= 0) exit
          val = trim(val) // ' ' // trim(line)
          if (index(trim(line), ']') /= 0) exit
        end do
      end if
      print *, '[RETURN CHECK] i=', nentry, ' key=', trim(key), ' val=', trim(val)
      select case (key)
      case('SetDesc')
        call copy_c_string(meta%SetDesc, val)
      case('SetIndex')
        read(val, *, iostat=ierr) meta%SetIndex
      case('Authors')
        call copy_c_string(meta%Authors, val)
      case('Reference')
        call copy_c_string(meta%Reference, val)
      case('Format')
        call copy_c_string(meta%Format, val)
      case('DataVersion')
        call copy_c_string(meta%DataVersion, val)
      case('NumMembers')
        read(val, *, iostat=ierr) meta%NumMembers
      case('Particle')
        read(val, *, iostat=ierr) meta%Particle
      case('OrderQCD')
        read(val, *, iostat=ierr) meta%OrderQCD
      case('FlavorScheme')
        call copy_c_string(meta%FlavorScheme, val)
      case('NumFlavors')
        read(val, *, iostat=ierr) meta%NumFlavors
      case('ErrorType')
        call copy_c_string(meta%ErrorType, val)
      case('XMin')
        read(val, *, iostat=ierr) meta%XMin
      case('XMax')
        read(val, *, iostat=ierr) meta%XMax
      case('QMin')
        read(val, *, iostat=ierr) meta%QMin
      case('QMax')
        read(val, *, iostat=ierr) meta%QMax
      case('MZ')
        read(val, *, iostat=ierr) meta%MZ
      case('MUp')
        read(val, *, iostat=ierr) meta%MUp
      case('MDown')
        read(val, *, iostat=ierr) meta%MDown
      case('MStrange')
        read(val, *, iostat=ierr) meta%MStrange
      case('MCharm')
        read(val, *, iostat=ierr) meta%MCharm
      case('MBottom')
        read(val, *, iostat=ierr) meta%MBottom
      case('MTop')
        read(val, *, iostat=ierr) meta%MTop
      case('AlphaS_MZ')
        read(val, *, iostat=ierr) meta%AlphaS_MZ
      case('AlphaS_OrderQCD')
        read(val, *, iostat=ierr) meta%AlphaS_OrderQCD
      case('AlphaS_Type')
        call copy_c_string(meta%AlphaS_Type, val)
      case('AlphaS_Qs')
        call parse_real_tokens(val, store_alpha_qs, n_alpha_points)
        if (n_alpha_points > 0) then
          meta%AlphaS_Qs = transfer(c_loc(store_alpha_qs(1)), 0_c_intptr_t)
        end if
      case('AlphaS_Vals')
        call parse_real_tokens(val, store_alpha_vals, n_alpha_points)
        if (n_alpha_points > 0) then
          meta%AlphaS_Vals = transfer(c_loc(store_alpha_vals(1)), 0_c_intptr_t)
        end if
      end select
    end do
    close(unit)
  end function pdf_read_metadata_c

  subroutine parse_real_tokens(line, vals, nvals)
    character(len=*), intent(in) :: line
    real(dp), allocatable, intent(out) :: vals(:)
    integer, intent(out) :: nvals
    integer :: i, j, ios, ntotal
    character(len=64) :: token
    real(dp) :: tmp

    nvals = 0
    i = 1
    do while (i <= len_trim(line))
      if (line(i:i) == ' ' .or. line(i:i) == char(9) .or. line(i:i) == ',') then
        i = i + 1
        cycle
      end if
      j = i
      do while (j <= len_trim(line))
        if (line(j:j) == ' ' .or. line(j:j) == char(9) .or. line(j:j) == ',') exit
        j = j + 1
      end do
      token = line(i:j-1)
      read(token, *, iostat=ios) tmp
      if (ios == 0) nvals = nvals + 1
      i = j
    end do

    if (nvals == 0) then
      allocate(vals(0))
      return
    end if

    allocate(vals(nvals))
    ntotal = 0
    i = 1
    do while (i <= len_trim(line))
      if (line(i:i) == ' ' .or. line(i:i) == char(9) .or. line(i:i) == ',') then
        i = i + 1
        cycle
      end if
      j = i
      do while (j <= len_trim(line))
        if (line(j:j) == ' ' .or. line(j:j) == char(9) .or. line(j:j) == ',') exit
        j = j + 1
      end do
      token = line(i:j-1)
      read(token, *, iostat=ios) tmp
      if (ios == 0) then
        ntotal = ntotal + 1
        vals(ntotal) = tmp
      end if
      i = j
    end do
  end subroutine parse_real_tokens

  subroutine parse_int_tokens(line, vals, nvals)
    character(len=*), intent(in) :: line
    integer, allocatable, intent(out) :: vals(:)
    integer, intent(out) :: nvals
    integer :: i, j, ios, ntotal
    character(len=64) :: token
    integer :: tmp

    nvals = 0
    i = 1
    do while (i <= len_trim(line))
      if (line(i:i) == ' ' .or. line(i:i) == char(9) .or. line(i:i) == ',') then
        i = i + 1
        cycle
      end if
      j = i
      do while (j <= len_trim(line))
        if (line(j:j) == ' ' .or. line(j:j) == char(9) .or. line(j:j) == ',') exit
        j = j + 1
      end do
      token = line(i:j-1)
      read(token, *, iostat=ios) tmp
      if (ios == 0) nvals = nvals + 1
      i = j
    end do

    if (nvals == 0) then
      allocate(vals(0))
      return
    end if

    allocate(vals(nvals))
    ntotal = 0
    i = 1
    do while (i <= len_trim(line))
      if (line(i:i) == ' ' .or. line(i:i) == char(9) .or. line(i:i) == ',') then
        i = i + 1
        cycle
      end if
      j = i
      do while (j <= len_trim(line))
        if (line(j:j) == ' ' .or. line(j:j) == char(9) .or. line(j:j) == ',') exit
        j = j + 1
      end do
      token = line(i:j-1)
      read(token, *, iostat=ios) tmp
      if (ios == 0) then
        ntotal = ntotal + 1
        vals(ntotal) = tmp
      end if
      i = j
    end do
  end subroutine parse_int_tokens

  function read_pdfdat(filename, ierr) result(data)
    character(len=*), intent(in) :: filename
    integer, intent(out), optional :: ierr
    type(pdfdat_t) :: data

    integer :: unit, ios, i, ncount, nrow, j, nf
    character(len=2000000) :: line
    real(dp), allocatable :: tmp_x(:), tmp_q(:), tmp_f(:,:,:), tmp_vals(:)
    integer, allocatable :: tmp_flav(:), tmp_ints(:)

    data = pdfdat_t()

    write(*,*) 'DEBUG: read_pdfdat - opening file', trim(filename)
    open(newunit=unit, file=filename, status='old', action='read', iostat=ios)
    if (ios /= 0) then
       write(*,*) 'ERROR: read_pdfdat - failed to open file', trim(filename), 'iostat =', ios
       if (present(ierr)) ierr = ios
       return
    end if

    read(unit, '(A)', iostat=ios) line
    if (ios /= 0) then
       write(*,*) 'ERROR: read_pdfdat - missing PdfType header'
       if (present(ierr)) ierr = -1
       close(unit)
       return
    end if
    read(unit, '(A)', iostat=ios) line
    if (ios /= 0) then
       write(*,*) 'ERROR: read_pdfdat - missing Format header'
       if (present(ierr)) ierr = -2
       close(unit)
       return
    end if
    read(unit, '(A)', iostat=ios) line
    if (ios /= 0) then
       write(*,*) 'ERROR: read_pdfdat - missing separator line'
       if (present(ierr)) ierr = -3
       close(unit)
       return
    end if

    do
       read(unit, '(A)', iostat=ios) line
       if (ios /= 0) then
          write(*,*) 'ERROR: read_pdfdat - EOF before x-grid (ios=', ios, ')'
          if (present(ierr)) ierr = -4
          close(unit)
          return
       end if
       if (len_trim(line) == 0) cycle
       exit
    end do

    write(*,*) 'DEBUG: read_pdfdat - raw x-grid line:', trim(line(1:min(len_trim(line), 200)))
    call parse_real_tokens(line, tmp_vals, ncount)
    if (ncount == 0) then
       write(*,*) 'ERROR: read_pdfdat - failed to parse x-grid line'
       if (present(ierr)) ierr = -5
       close(unit)
       return
    end if
    allocate(tmp_x(ncount))
    tmp_x = tmp_vals
    data%n = ncount

    do
       read(unit, '(A)', iostat=ios) line
       if (ios /= 0) then
          write(*,*) 'ERROR: read_pdfdat - EOF before Q-grid (ios=', ios, ')'
          if (present(ierr)) ierr = -6
          close(unit)
          return
       end if
       if (len_trim(line) == 0) cycle
       exit
    end do

    write(*,*) 'DEBUG: read_pdfdat - raw Q-grid line:', trim(line(1:min(len_trim(line), 200)))
    call parse_real_tokens(line, tmp_vals, ncount)
    if (ncount == 0) then
       write(*,*) 'ERROR: read_pdfdat - failed to parse Q-grid line'
       if (present(ierr)) ierr = -7
       close(unit)
       return
    end if
    allocate(tmp_q(ncount))
    tmp_q = tmp_vals
    data%nQ = ncount

    do
       read(unit, '(A)', iostat=ios) line
       if (ios /= 0) then
          allocate(tmp_flav(0))
          exit
       end if
       if (len_trim(line) == 0) cycle
       exit
    end do

    if (ios == 0) then
       call parse_int_tokens(line, tmp_flav, ncount)
       if (ncount == 0) then
          allocate(tmp_flav(0))
       end if
    end if

    data%nFlav = max(size(tmp_flav), 1)
    allocate(tmp_f(data%nQ, data%n, data%nFlav))
    tmp_f = 0.0_dp

    i = 0
    do while (i < data%n * data%nQ)
       read(unit, '(A)', iostat=ios) line
       if (ios /= 0) then
          write(*,*) 'ERROR: read_pdfdat - reached EOF while reading data rows; expected', data%n * data%nQ, 'entries but got', i
          if (present(ierr)) ierr = -8
          close(unit)
          return
       end if
       if (len_trim(line) == 0) cycle

       call parse_real_tokens(line, tmp_vals, nrow)
       if (nrow == 0) cycle

       i = i + 1
       nf = min(nrow, data%nFlav)
       tmp_f(mod(i-1, data%nQ) + 1, (i-1)/data%nQ + 1, 1:nf) = tmp_vals(1:nf)
    end do

    close(unit)

    if (i /= data%n * data%nQ) then
       write(*,*) 'WARNING: read_pdfdat - expected', data%n * data%nQ, 'grid entries but read', i, 'entries'
    end if

    if (allocated(store_x)) deallocate(store_x)
    if (allocated(store_q)) deallocate(store_q)
    if (allocated(store_f)) deallocate(store_f)
    if (allocated(store_flav)) deallocate(store_flav)
    allocate(store_x(data%n), store_q(data%nQ), store_f(data%nQ, data%n, data%nFlav), store_flav(data%nFlav))

    store_x = real(tmp_x, c_double)
    store_q = real(tmp_q, c_double)
    store_f = real(tmp_f, c_double)
    if (allocated(tmp_flav) .and. size(tmp_flav) > 0) then
       store_flav = int(tmp_flav, c_int)
    else
       store_flav = 0_c_int
    end if

    data%x_ptr = c_loc(store_x(1))
    data%q_ptr = c_loc(store_q(1))
    data%f_ptr = c_loc(store_f(1,1,1))
    data%flav_ptr = c_loc(store_flav(1))

    if (present(ierr)) ierr = 0
  end function read_pdfdat

  function pdf_read_grid_c(fname, flen) result(grid) bind(C)
    character(kind=c_char), intent(in) :: fname(*)
    integer(c_int), value :: flen
    type(pdfdat_t) :: grid
    character(len=:), allocatable :: f_name
    integer :: i, ierr

    f_name = ''
    do i = 1, flen
      if (fname(i) == c_null_char) exit
      f_name = f_name // fname(i)
    end do

    grid = read_pdfdat(f_name, ierr)
  end function pdf_read_grid_c

  subroutine get_pdfdat_point(data, i, xout, fout, ierr, iflav)
    type(pdfdat_t), intent(in) :: data
    integer, intent(in) :: i
    real(dp), intent(out) :: xout, fout
    integer, intent(out), optional :: ierr
    integer, intent(in), optional :: iflav
    integer :: jflav
    real(c_double), pointer :: xarr(:), farr(:,:,:)
    if (i < 1 .or. i > data%n) then
      if (present(ierr)) ierr = 1
      return
    end if
    jflav = 1
    if (present(iflav)) jflav = iflav
    call c_f_pointer(data%x_ptr, xarr, [data%n])
    xout = real(xarr(i), dp)
    if (data%nQ >= 1 .and. jflav >= 1 .and. jflav <= data%nFlav) then
      call c_f_pointer(data%f_ptr, farr, [data%nQ, data%n, data%nFlav])
      fout = real(farr(1, i, jflav), dp)
    else
      fout = 0.0_dp
    end if
    if (present(ierr)) ierr = 0
  end subroutine get_pdfdat_point

  real(c_double) function pdf_get_alphas_c(meta, Q) result(val) bind(C)
    type(pdf_meta_t), value, intent(in) :: meta
    real(c_double), value :: Q
    type(c_ptr) :: qptr, vptr
    real(c_double), pointer :: qarr(:), varr(:)
    integer(c_int) :: n
    integer :: i
    real(dp) :: qval, q0, q1, q2, f0, f1, f2, l0, l1, l2

    val = 0.0_c_double
    if (meta%AlphaS_Qs == 0_c_intptr_t .or. meta%AlphaS_Vals == 0_c_intptr_t .or. n_alpha_points < 2) then
      val = meta%AlphaS_MZ
      return
    end if

    qptr = transfer(meta%AlphaS_Qs, c_null_ptr)
    vptr = transfer(meta%AlphaS_Vals, c_null_ptr)
    n = n_alpha_points
    call c_f_pointer(qptr, qarr, [n])
    call c_f_pointer(vptr, varr, [n])

    qval = real(Q, dp)
    if (qval <= qarr(1)) then
      val = real(varr(1), c_double)
      return
    else if (qval >= qarr(n)) then
      val = real(varr(n), c_double)
      return
    end if

    if (n == 2) then
      q1 = real(qarr(1), dp)
      q2 = real(qarr(2), dp)
      f1 = real(varr(1), dp)
      f2 = real(varr(2), dp)
      if (q2 == q1) then
        val = real(f1, c_double)
      else
        l1 = (qval - q1) / (q2 - q1)
        val = real((1.0_dp - l1) * f1 + l1 * f2, c_double)
      end if
      return
    end if

    do i = 2, n-1
      if (qval <= qarr(i)) then
        q0 = real(qarr(max(1, i-1)), dp)
        q1 = real(qarr(i), dp)
        q2 = real(qarr(min(n, i+1)), dp)
        f0 = real(varr(max(1, i-1)), dp)
        f1 = real(varr(i), dp)
        f2 = real(varr(min(n, i+1)), dp)
        if (qval <= q0) then
          q0 = real(qarr(1), dp)
          q1 = real(qarr(2), dp)
          q2 = real(qarr(3), dp)
          f0 = real(varr(1), dp)
          f1 = real(varr(2), dp)
          f2 = real(varr(3), dp)
        else if (qval >= q2) then
          q0 = real(qarr(n-2), dp)
          q1 = real(qarr(n-1), dp)
          q2 = real(qarr(n), dp)
          f0 = real(varr(n-2), dp)
          f1 = real(varr(n-1), dp)
          f2 = real(varr(n), dp)
        end if
        l0 = (qval - q1) * (qval - q2) / ((q0 - q1) * (q0 - q2))
        l1 = (qval - q0) * (qval - q2) / ((q1 - q0) * (q1 - q2))
        l2 = (qval - q0) * (qval - q1) / ((q2 - q0) * (q2 - q1))
        val = real(l0 * f0 + l1 * f1 + l2 * f2, c_double)
        return
      end if
    end do

    val = real(varr(n), c_double)
  end function pdf_get_alphas_c

  subroutine pdf_get_values_c(grid, x, Q, flv, val) bind(C)
    type(pdfdat_t), value, intent(in) :: grid
    real(c_double), value :: x, Q
    integer(c_int), value :: flv
    real(c_double), intent(out) :: val
    integer :: n, nQcols, nflav, ix_lo, ix_hi, iq_lo, iq_hi, j, jflav, pdgcode
    real(c_double), pointer :: xarr(:), qarr(:), farr(:,:,:)
    integer(c_int), pointer :: flavarr(:)
    real(dp) :: xval, qval, x1, x2, q1, q2
    real(dp) :: f11, f12, f21, f22, tx, tq

    n = grid%n
    nQcols = grid%nQ
    nflav = grid%nFlav
    if (n == 0 .or. nQcols == 0 .or. nflav == 0) then
      val = 0.0_c_double
      return
    end if

    call c_f_pointer(grid%x_ptr, xarr, [n])
    call c_f_pointer(grid%q_ptr, qarr, [nQcols])
    call c_f_pointer(grid%f_ptr, farr, [nQcols, n, nflav])
    call c_f_pointer(grid%flav_ptr, flavarr, [nflav])

    pdgcode = flv
    if (flv == 0) pdgcode = 21

    jflav = 0
    do j = 1, nflav
      if (flavarr(j) == pdgcode) then
        jflav = j
        exit
      end if
    end do
    if (jflav == 0) then
      val = 0.0_c_double
      return
    end if

    xval = real(x, dp)
    qval = real(Q, dp)

    ix_lo = 1
    ix_hi = 1
    if (xval <= xarr(1)) then
      ix_hi = min(2, n)
    else if (xval >= xarr(n)) then
      ix_lo = max(1, n-1)
      ix_hi = n
    else
      do j = 1, n-1
        if (xval >= xarr(j) .and. xval <= xarr(j+1)) then
          ix_lo = j
          ix_hi = j + 1
          exit
        end if
      end do
    end if

    iq_lo = 1
    iq_hi = 1
    if (qval <= qarr(1)) then
      iq_hi = min(2, nQcols)
    else if (qval >= qarr(nQcols)) then
      iq_lo = max(1, nQcols-1)
      iq_hi = nQcols
    else
      do j = 1, nQcols-1
        if (qval >= qarr(j) .and. qval <= qarr(j+1)) then
          iq_lo = j
          iq_hi = j + 1
          exit
        end if
      end do
    end if

    if (ix_lo == ix_hi .and. iq_lo == iq_hi) then
      val = real(farr(iq_lo, ix_lo, jflav), c_double)
      return
    end if

    x1 = xarr(ix_lo)
    x2 = xarr(ix_hi)
    q1 = qarr(iq_lo)
    q2 = qarr(iq_hi)

    f11 = farr(iq_lo, ix_lo, jflav)
    f12 = farr(iq_hi, ix_lo, jflav)
    f21 = farr(iq_lo, ix_hi, jflav)
    f22 = farr(iq_hi, ix_hi, jflav)

    if (x2 == x1) then
      val = real(f11, c_double)
      return
    end if
    if (q2 == q1) then
      val = real(f11, c_double)
      return
    end if

    tx = (xval - x1) / (x2 - x1)
    tq = (qval - q1) / (q2 - q1)

    val = real((1.0_dp - tx) * (1.0_dp - tq) * f11 + &
               (1.0_dp - tx) * tq * f12 + &
               tx * (1.0_dp - tq) * f21 + &
               tx * tq * f22, c_double)
  end subroutine pdf_get_values_c

end module pdf_c_interface
