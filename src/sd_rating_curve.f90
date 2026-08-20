      subroutine sd_rating_curve (i)
      
      use sd_channel_module
      use channel_velocity_module
      use maximum_data_module
      !use hydrograph_module
      
      implicit none      

      integer, intent (in) :: i     !none           |counter  
      integer :: i_dep              !none           |counter
      integer :: ifp_dep            !none           |counter
      
      real :: a                     !m^2            |cross-sectional area of channel
      real :: b                     !m              |bottom width of channel
      real :: p                     !m              |wetting perimeter
      real :: rh                    !m              |hydraulic radius
      real :: qman                  !m^3/s or m/s   |flow rate or flow velocity
      real :: dep                   !               |
      real :: a_bf                  !               |
      real :: a_fp                  !m^2            |cross-sectional area on the flood plain
      real :: p_bf                  !               |
      real :: vol_bf                !m^3            |channel volume at bankfull depth
      real :: vel                   !               |
      real :: frac_abov

      b = sd_ch(i)%chw - 2. * sd_ch(i)%chd * sd_ch(i)%chss
      !! check if bottom width (b) is < 0
      if (b <= 0.) then
        b = .5 * sd_ch(i)%chw
        b = Max(0., b)
        if (sd_ch(i)%chd > 1.e-9) then
          sd_ch(i)%chss = (sd_ch(i)%chw - b) / (2. * sd_ch(i)%chd)
        else
          sd_ch(i)%chss = 0.
        end if
      end if
      ch_rcurv(i)%wid_btm = b
      
      !! compute rating curve at 0.1 and 1.0 times bankfull depth
      do i_dep = 1, 2
        if (i_dep == 1) dep = 0.1 * sd_ch(i)%chd
        if (i_dep == 2) dep = sd_ch(i)%chd
        !! chss is horizontal change per unit vertical change
        p = b + 2. * dep * Sqrt(1. + sd_ch(i)%chss ** 2)
        a = b * dep + sd_ch(i)%chss * dep ** 2
        if (p > 1.e-9 .and. a > 1.e-9) then
          rh = a / p
        else
          rh = 0.
        end if
        ch_rcurv(i)%elev(i_dep)%dep = dep
        ch_rcurv(i)%elev(i_dep)%wet_perim = p
        ch_rcurv(i)%elev(i_dep)%xsec_area = a
        
        ch_rcurv(i)%elev(i_dep)%top_wid = b + 2. * dep * sd_ch(i)%chss
        !! m2 = m * km * 1000 m/km
        ch_rcurv(i)%elev(i_dep)%surf_area = ch_rcurv(i)%elev(i_dep)%top_wid * sd_ch(i)%chl * 1000.
        ch_rcurv(i)%elev(i_dep)%vol = a * sd_ch(i)%chl * 1000.
        ch_rcurv(i)%elev(i_dep)%vol_ch = ch_rcurv(i)%elev(i_dep)%vol
        ch_rcurv(i)%elev(i_dep)%vol_fp = 0.
        
        ch_rcurv(i)%elev(i_dep)%flo_rate = Qman(a, rh, sd_ch(i)%chn, sd_ch(i)%chs)
        vel = Qman(1., rh, sd_ch(i)%chn, sd_ch(i)%chs)
        if (vel > 1.e-9) then
          ch_rcurv(i)%elev(i_dep)%ttime = sd_ch(i)%chl / (3.6 * vel)
        else
          ch_rcurv(i)%elev(i_dep)%ttime = 0.
        end if
        
        !! save bankfull depth and area for flood plain calculations
        if (i_dep == 2) then
          p_bf = p
          a_bf = a
          vol_bf = ch_rcurv(i)%elev(i_dep)%vol_ch
        end if
      end do
        
      !! compute rating curve at 1.2 and 2.0 times bankfull depth (flood plain)
      do i_dep = 1, 2
        !! dep = depth above bankfull
        if (i_dep == 1) frac_abov = 0.2
        if (i_dep == 2) frac_abov = 1.
        dep = frac_abov * sd_ch(i)%chd
        !! flood plain perimeter - p^2 = dep^2 + width^2
        p = p_bf + 2. * Sqrt(dep ** 2 * (1. + 1. / (sd_ch(i)%fps ** 2)))
        !! fps is vertical change per unit horizontal change
        !! total area includes water above the bankfull channel and flood-plain wedges
        a_fp = dep ** 2 / sd_ch(i)%fps
        a = a_bf + sd_ch(i)%chw * dep + a_fp
        if (p > 1.e-9 .and. a > 1.e-9) then
          rh = a / p
        else
          rh = 0.
        end if
        ifp_dep = i_dep + 2
        ch_rcurv(i)%elev(ifp_dep)%dep = (1. + frac_abov) * sd_ch(i)%chd
        ch_rcurv(i)%elev(ifp_dep)%wet_perim = p
        ch_rcurv(i)%elev(ifp_dep)%xsec_area = a
        ch_rcurv(i)%elev(ifp_dep)%top_wid = sd_ch(i)%chw + 2. * dep / sd_ch(i)%fps
        !! m2 = m * km * 1000 m/km
        ch_rcurv(i)%elev(ifp_dep)%surf_area = ch_rcurv(i)%elev(ifp_dep)%top_wid * sd_ch(i)%chl * 1000.
        !! routing stores the channel up to bankfull; all excess water is flood-plain storage
        ch_rcurv(i)%elev(ifp_dep)%vol_ch = vol_bf
        ch_rcurv(i)%elev(ifp_dep)%vol_fp = (a - a_bf) * sd_ch(i)%chl * 1000.
        ch_rcurv(i)%elev(ifp_dep)%vol = ch_rcurv(i)%elev(ifp_dep)%vol_ch +  &
                                                    ch_rcurv(i)%elev(ifp_dep)%vol_fp
        ch_rcurv(i)%elev(ifp_dep)%flo_rate = Qman(a, rh, sd_ch(i)%fpn, sd_ch(i)%chs)
        vel = Qman(1., rh, sd_ch(i)%fpn, sd_ch(i)%chs)
        if (vel > 1.e-9) then
          ch_rcurv(i)%elev(ifp_dep)%ttime = sd_ch(i)%chl / (3.6 * vel)
        else
          ch_rcurv(i)%elev(ifp_dep)%ttime = 0.
        end if
      end do

      return
      end subroutine sd_rating_curve
