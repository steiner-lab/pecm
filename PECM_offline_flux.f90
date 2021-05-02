program fluxcode

! MODEL VERSION: PECM1.1
! Modified by Yingxiao, read 15 diff species 02/05/20
! Modified by Yingxiao on 24/03/20, add precipitation
! Final Code Update: Matthew Wozniak, 08/02/16
! Last ALS update: 1 June 2015 - substantial updates to data input; creation of netcdf output, phenological calculations
! Based on earlier version (beldread.f90, partially developed by MAZ)
! UPDATE: now includes normalization of emission factor (normalized to Gaussian curve) (06/2015)
! UPDATE: phenology now parameterized with sDOY and eDOY versus PYAAT (5/2016)
! UPDATE: final emissions output is aggregated by PFT (DBL, ENL, Grasses, Ragweed)
! UPDATE: Ragweed landcover is given by CLM crop (corn and soybean) PFTs and CLM urban level 3. Parameterization is drawn from Clay et al. 2006 (crop ragweed) and Katz et al. 2014 (urban ragweed) ragweed densities. 
! Crops are split into tilled vs. no-till (fractions given by USDA)
! UPDATE: Ragweed crop parameterization reduced by a factor of 100 based on comparison of test-run RegCM-pollen simulated surface counts versus observed ragweed counts. Grass
! UPDATE: Grasses' emission factor has been reduced to minimum from literature for C3 and x10^-1 the minimum for C4

use netcdf

implicit none

logical :: plac_flag
integer :: nx,ny,nyr,ndy,ndy_nextyr,nspec,nbeld,nclm,cnt,cnt_yrs,cnt_mnths,nlev,yrst,yrend,ndytot,ngra,nrag,ntree,cnt_begin,cnt_end     ! mcwoz added nlev, yrst, yrend, ndytot, nclm, ngra, nrag
integer :: idbeld,idLat,idLong,idyear,idday,idclm,idpft,idurb,idc3,idc4,idurblev ! mcwoz added idc3, idc4, idurblev
integer :: idcru,idTmp,status,dimval,nlats,idcpc,idPr
integer :: len, wid
integer :: phen_dimids(4),phen_count(4)
integer :: i,ii,iii,iv, iv2,jj,id_2d(2),start(4),start_all(5),count(4),count_all(5),start1(4),count1(4),rcode,pcoutid,ncidout ! mcwoz added jj
integer :: dimid1,dimid2,dimid3,dimid4,dimid5,dimids(4),dimids_all(5),dimid6,dimid7
integer :: lonoutid,latoutid,dayoutid,yroutid,specoutid,timeoutid,levoutid,lcoutid,sdoyoutid,edoyoutid,duroutid,dbloutid,enloutid,graoutid,ragoutid,fluxoutid ! mcwoz added dbl-,enl-,gra-,ragoutid
real,allocatable :: sDOY_m(:),sDOY_b(:),eDOY_m(:),eDOY_b(:),ef(:)
double precision,allocatable   :: tval(:,:,:),pft(:,:,:),urban(:,:,:)
real,allocatable :: xlat(:,:),xlon(:,:),grass(:,:)
integer,allocatable   :: ispec(:),idspec(:),ttest(:,:,:),lattest(:,:)
integer,allocatable   :: iyr(:),idy(:),idcrop(:),iddbl(:),idenl(:) ! mcwoz added idcrop
real, allocatable  :: pyaat_sum(:,:,:),pyaat(:,:,:),pyaat_sum_adj(:,:,:),pyaat_adj(:,:,:)    ! pyaat = prior year annual average temperature
real, allocatable :: frac_area(:,:,:),beld(:,:,:),sDOY(:,:,:,:),eDOY(:,:,:,:),dur(:,:,:,:),pc(:,:,:,:),dbl_poll(:,:,:,:),enl_poll(:,:,:,:),gra_poll(:,:,:,:),rag_poll(:,:,:,:)
real, allocatable :: flux_all(:,:,:,:,:)
real, allocatable :: pr(:,:,:), cdur(:)
real, dimension(103,2) :: stat_llind
real :: sigma, mu, gphen,sigma1,sigma2,mu1,mu2,tot_flx,ef_d ! mcwoz added tot_flx, ef_d
real :: vf = 0.10 ! for application of flowering probability (vitality) to emission factor
character(len=256) :: longname(8), units(8), lc_option
character(4),allocatable :: name(:), year(:)
yrst = 1 ! Index for start year of emissions
yrend = 10! Index for end year of emissions
nyr = yrend-yrst + 1 ! Number of years for this run
units(4) = 'days since 1995-01-01 00:00:00' ! First time index refers to 2003-01-01 00:00:00 with these units

! Indicate below which land cover specification to use ("PFT" = PFT-level from CLM, "TAX" taxon-level from BELD)
lc_option = "TAX"

if (lc_option(1:3) .eq. "TAX") then
! Read input file scalars
      open(unit=20,file='/nfs/alsteine/shared_data/PECM/input_parameters/emission_params_beld.in',status='unknown')
      read(20,*) nx              ! Number of grid points in x direction (longitude)
      read(20,*) ny              ! Number of grid points in y direction (latitude)
      read(20,*) ndy             ! number of days per year (calendar definition)
      read(20,*) nbeld           ! number of tree taxa in the BELD species database (currently 11)
      read(20,*) ngra            ! number of grass taxa (currently 2: C3 and C4)
      read(20,*) nrag             ! number of ragweed taxa (currently 1)
      read(20,*) idc3            ! index of c3 non-arctic grasses pft in CLM species data
      read(20,*) idc4            ! index of c4 grasses pft in CLM species data
      read(20,*) idurblev        ! index of the urban level used for ragweed parameterization
      read(20,*) nspec           ! number of pollen genus/type to include in model
      read(20,*) nlev            ! number of vertical levels in emissions model

allocate(year(nyr))
year = (/'1995','1996','1997','1998','1999','2000','2001','2002','2003','2004'/)

! count the total number of days in the emissions calculation accounting for leap days
ndytot = 0 
do iv = yrst,yrend
  if (year(iv) .eq. '1996' .or. year(iv) .eq. '2000' .or. year(iv) .eq. '2004') then 
      ndytot = ndytot + 366
   else
      ndytot = ndytot + 365
   end if
end do

!--- Allocate memory for arrays
allocate(xlon(nx,ny))     ! note CRU dims are different
allocate(xlat(nx,ny))
allocate(iyr(nyr))
allocate(idy(ndytot))
allocate(idspec(nspec))
allocate(frac_area(nx,ny,nspec))
allocate(beld(nx,ny,nbeld))
allocate(tval(nx,ny,96))  ! note revised format for CRU data and that CRU dims are different
allocate(pyaat(nx,ny,nyr))
allocate(pyaat_sum(nx,ny,nyr))
allocate(pyaat_adj(nx,ny,nyr))
allocate(pyaat_sum_adj(nx,ny,nyr))
allocate(sDOY_m(nspec))
allocate(sDOY_b(nspec))
allocate(eDOY_m(nspec))
allocate(eDOY_b(nspec))
allocate(ef(nspec))
allocate(sDOY(nx,ny,nyr,nspec))
allocate(dur(nx,ny,nyr,nspec))
allocate(eDOY(nx,ny,nyr,nspec))
!allocate(pc(nx,ny,nlev,ndytot))
allocate(dbl_poll(nx,ny,nlev,ndytot))
allocate(pr(nx,ny,ndytot))
allocate(cdur(nspec))
allocate(enl_poll(nx,ny,nlev,ndytot))
allocate(gra_poll(nx,ny,nlev,ndytot))
allocate(rag_poll(nx,ny,nlev,ndytot))
allocate(flux_all(nx,ny,nlev,ndytot,nspec))
allocate(name(nspec))
allocate(ttest(nx,ny,144))
allocate(lattest(ny,nx))
allocate(urban(nx,ny,3))
allocate(grass(ny,nx))
allocate(pft(nx,ny,25))
allocate(idcrop(4)) ! currently just corn and soybean crops

!cdur = (/37,37,33,43,29,33,47,24,30,44,29,27,53,28,51/)
! set idy to integer counting array counting from 1 to the number of days of emissions
idy = (/ (i, i = 1,ndytot) /)
iyr = (/ (i, i = 1,nyr) /)

! read in rest of input file
      read(20,*) (idcrop(i),i=1,4)               ! read in crop pft indices (2 corn, 2 soybean)
      read(20,*) (name(i),i=1,nspec)             ! read the names of the emitting taxa
      read(20,*) (ef(i),i=1,nspec)               ! read in emission factors (grains/m2; RW is grains/plant) 
      read(20,*) (sDOY_m(i),i=1,nspec)           ! read in sDOY regression slope 
      read(20,*) (sDOY_b(i),i=1,nspec)           ! read in sDOY regression intercept 
      read(20,*) (eDOY_m(i),i=1,nspec)            ! read in eDOY regression slope 
      read(20,*) (eDOY_b(i),i=1,nspec)            ! read in eDOY intercept

! print the phenological parameters to screen
do i=1,nspec
write(*,*) name(i),ef(i),sDOY_m(i),sDOY_b(i),eDOY_m(i),eDOY_b(i)
end do

!---Read in BELD3 25km res lat/lon file---------------
status = nf90_open(path = "/nfs/alsteine/shared_data/PECM/input_data/merged_beld_25km.nc", mode = nf90_nowrite, ncid = idbeld)

do i=1,nbeld
  status = nf90_inq_varid(idbeld, name(i), idspec(i))
  status = nf90_get_var(idbeld, idspec(i), beld(:,:,i))
end do

frac_area(:,:,1:nbeld) = beld(:,:,:)

status = nf90_close(idbeld)

! read in CLM land cover data (PFTS for grasses and crops, urban landcover)
status = nf90_open(path = "/nfs/alsteine/shared_data/PECM/input_data/CLM45_surface_emisvars.nc", mode = nf90_nowrite, ncid = idclm)
status = nf90_inq_varid(idclm, 'pft_2d', idpft)
status = nf90_get_var(idclm, idpft, pft)
status = nf90_inq_varid(idclm, 'urb_2d', idurb)
status = nf90_get_var(idclm, idurb, urban)

! Duplicate the elm landcover                                               
!do i=1,ny
!do ii=1,nx
   !frac_area(ii,i,12) = frac_area(ii,i,11)         
!end do
!end do
frac_area(:,:,12) = frac_area(:,:,11)

frac_area(:,:,:) = 0.01*frac_area(:,:,:) ! convert from percent to fraction (m2 veg/m2 total)
pft(:,:,:) = 0.01*pft(:,:,:)
urban(:,:,:) = 0.01*urban(:,:,:)

! grass landcover
nclm = ngra
frac_area(:,:,13) = real(pft(:,:,idc3),4)
frac_area(:,:,14) = real(pft(:,:,idc4),4)
!frac_area(:,:,13) = 1.                                                                  
!frac_area(:,:,14) = 1.
!write(*,*) size(pft,1),size(pft,2),size(pft,3)
!write(*,*) idc3,idc4

! Ragweed landcover; NOTE: in plants/m2, NOT fractional area
frac_area(:,:,15) = real(0.01*(2.*0.75 + 10.*0.25)*(pft(:,:,idcrop(1)) + pft(:,:,idcrop(2))) + 0.01*(2.*0.55 + 10.*0.45)*(pft(:,:,idcrop(3)) + pft(:,:,idcrop(4))) + 0.1*0.5*urban(:,:,idurblev),4)
! (t_density*till + nt_density*no-till)*(corn) + (t_density*till + nt_density*no-till)*(soybean) + urb_density*urban
! now with crop density reduced by factor of 100, urban reduced by factor of 10
   
do i = 1,nx
   do ii = 1,ny
      do iii = 13,15
         if (frac_area(i,ii,iii).lt.0) then
            frac_area(i,ii,iii) = 0.
         end if
      end do
   end do
end do

end if ! ends "if ('TAX')" block


if (lc_option(1:3) .eq. "PFT") then

! Read input file scalars
      open(unit=20,file='/nfs/alsteine/shared_data/PECM/input_parameters/emission_params_pft.in',status='unknown')
      read(20,*) nx              ! Number of grid points in x direction (longitude)
      read(20,*) ny              ! Number of grid points in y direction (latitude)
      read(20,*) ndy             ! number of days per year (calendar definition)
      read(20,*) ntree           ! number of tree plant functional types (currently 2: DBL, ENL)
      read(20,*) ngra            ! number of grass taxa
      read(20,*) nrag             ! number of ragweed taxa (currently 1)
      read(20,*) idc3            ! index of c3 non-arctic grasses pft in CLM species data
      read(20,*) idc4            ! index of c4 grasses pft in CLM species data
      read(20,*) idurblev        ! index of the urban level used for ragweed parameterization
      read(20,*) nspec           ! number of pollen genus/type to include in model
      read(20,*) nlev            ! number of vertical levels in emissions model

allocate(year(nyr))
year = (/'1995','1996','1997','1998','1999','2000','2001','2002','2003','2004'/)

! count the total number of days in the emissions calculation accounting for leap days
ndytot = 0 
do iv = yrst,yrend
  if (year(iv) .eq. '1996' .or. year(iv) .eq. '2000' .or. year(iv) .eq. '2004') then 
      ndytot = ndytot + 366
   else
      ndytot = ndytot + 365
   end if
end do


!--- Allocate memory for arrays
allocate(xlon(nx,ny))     ! note CRU dims are different
allocate(xlat(nx,ny))
allocate(iyr(nyr))
allocate(idy(ndytot))
allocate(idspec(nspec))
allocate(frac_area(nx,ny,nspec))
allocate(tval(nx,ny,96))  ! note revised format for CRU data and that CRU dims are different
allocate(pyaat(nx,ny,nyr))
allocate(pyaat_sum(nx,ny,nyr))
allocate(pyaat_adj(nx,ny,nyr))
allocate(pyaat_sum_adj(nx,ny,nyr))
allocate(sDOY_m(nspec))
allocate(sDOY_b(nspec))
allocate(eDOY_m(nspec))
allocate(eDOY_b(nspec))
allocate(ef(nspec))
allocate(sDOY(nx,ny,nyr,nspec))
allocate(dur(nx,ny,nyr,nspec))
allocate(eDOY(nx,ny,nyr,nspec))
!allocate(pc(nx,ny,nlev,ndytot))
allocate(dbl_poll(nx,ny,nlev,ndytot))
allocate(pr(nx,ny,ndytot))
allocate(enl_poll(nx,ny,nlev,ndytot))
allocate(gra_poll(nx,ny,nlev,ndytot))
allocate(rag_poll(nx,ny,nlev,ndytot))
allocate(flux_all(nx,ny,nlev,ndytot,nspec))
allocate(name(nspec))
allocate(ttest(nx,ny,144))
allocate(lattest(ny,nx))
allocate(urban(nx,ny,3))
allocate(grass(ny,nx))
allocate(pft(nx,ny,25))
allocate(idcrop(4)) ! currently just corn and soybean crops
allocate(iddbl(3))
allocate(idenl(2))

! set idy to integer counting array counting from 1 to the number of days of emissions
idy = (/ (i, i = 1,ndytot) /)
iyr = (/ (i, i = 1,nyr) /)
  write(*,*) "ndytot, nyr 2"
  write(*,*) ndytot, nyr
! read in rest of input file
      read(20,*) (idcrop(i),i=1,4)               ! read in crop pft indices
      read(20,*) (iddbl(i),i=1,3)
      read(20,*) (idenl(i),i=1,2)                ! recently updated to include boreal ENF
      read(20,*) (name(i),i=1,nspec)             ! read the names of the emitting taxa
      read(20,*) (ef(i),i=1,nspec)               ! read in emission factors (grains/m2; RW is grains/plant) 
      read(20,*) (sDOY_m(i),i=1,nspec)           ! read in sDOY regression slope 
      read(20,*) (sDOY_b(i),i=1,nspec)           ! read in sDOY regression intercept 
      read(20,*) (eDOY_m(i),i=1,nspec)            ! read in eDOY regression slope 
      read(20,*) (eDOY_b(i),i=1,nspec)            ! read in eDOY intercept

! print the phenological parameters to screen
do i=1,nspec
write(*,*) name(i),ef(i),sDOY_m(i),sDOY_b(i),eDOY_m(i),eDOY_b(i)
end do

! read in CLM land cover data (PFTS for grasses and crops, urban landcover)
status = nf90_open(path = "/nfs/alsteine/shared_data/PECM/input_data/CLM45_surface_emisvars.nc", mode = nf90_nowrite, ncid = idclm)
status = nf90_inq_varid(idclm, 'pft_2d', idpft)
status = nf90_get_var(idclm, idpft, pft)
status = nf90_inq_varid(idclm, 'urb_2d', idurb)
status = nf90_get_var(idclm, idurb, urban)

pft(:,:,:) = 0.01*pft(:,:,:)
urban(:,:,:) = 0.01*urban(:,:,:)

! tree landcover
frac_area(:,:,1) = pft(:,:,iddbl(1)) + pft(:,:,iddbl(2)) + pft(:,:,iddbl(3))
frac_area(:,:,2) = pft(:,:,idenl(1)) + pft(:,:,idenl(2))

! grass landcover
nclm = ngra
frac_area(:,:,3) = real(pft(:,:,idc3),4)
frac_area(:,:,4) = real(pft(:,:,idc4),4)
!frac_area(:,:,13) = 1.                                                                  
!frac_area(:,:,14) = 1.
!write(*,*) size(pft,1),size(pft,2),size(pft,3)
!write(*,*) idc3,idc4

! Ragweed landcover; NOTE: in plants/m2, NOT fractional area
frac_area(:,:,5) = real(0.01*(2.*0.75 + 10.*0.25)*(pft(:,:,idcrop(1)) + pft(:,:,idcrop(2))) + 0.01*(2.*0.55 + 10.*0.45)*(pft(:,:,idcrop(3)) + pft(:,:,idcrop(4))) + 0.1*0.5*urban(:,:,idurblev),4)
! (t_density*till + nt_density*no-till)*(corn) + (t_density*till + nt_density*no-till)*(soybean) + urb_density*urban
! now with crop density reduced by factor of 100
   
do i = 1,nx
   do ii = 1,ny
      do iii = 1,5
         if (frac_area(i,ii,iii).lt.0) then
            frac_area(i,ii,iii) = 0.
         end if
      end do
   end do
end do

end if ! ends "if ('PFT')" block

! officially close landcover files
status = nf90_close(idclm)

! officially close landcover files
status = nf90_close(idclm)

!---Read in daily precipitation ----------------

status = nf90_open(path = "/alsteine2/yingxz/PECM/input_data/CMIP6/GFDL_ESM4/CMIP6_noaa_hist_pr_day_1995-2014.nc", mode = nf90_nowrite, ncid = idcpc)
status = nf90_inq_varid(idcpc, "pr", idPr)
status = nf90_get_var(idcpc, idPr, pr, start=(/1,1,1/),count=(/nx,ny,ndytot/))



!---Read in latitude values ----------

status = nf90_inq_varid(idcpc, "xlat", idLat)
status = nf90_get_var(idcpc, idLat, xlat, start=(/1,1/),count=(/nx,ny/))
status = nf90_get_var(idcru, idLat, xlat)

!---Read in longitude values ----------

status = nf90_inq_varid(idcpc, "xlon", idLong)
status = nf90_get_var(idcpc, idLong, xlon, start=(/1,1/),count=(/nx,ny/))


if(status /= nf90_noerr) write(*,*)"netcdf error = 0 ",trim(nf90_strerror(status))
!write(*,*) tval(101,101,1:24)
status = nf90_close(idcpc)

!---Read in daily temperature  ----------------

status = nf90_open(path = "/nfs/alsteine/yingxz/PECM/input_data/CMIP6/GFDL_ESM4/CMIP6_noaa_hist_tas_aat_1994-2013.nc", mode = nf90_nowrite, ncid = idcru)
status = nf90_inq_varid(idcru, "aat", idTmp)
status = nf90_get_var(idcru, idTmp, pyaat, start=(/1,1,1/),count=(/nx,ny,nyr/))  
status = nf90_inq_varid(idcru, "aajt", idTmp)
status = nf90_get_var(idcru, idTmp, pyaat_adj, start=(/1,1,1/),count=(/nx,ny,nyr/))


!---Read in latitude values ----------

status = nf90_inq_varid(idcru, "xlat", idLat)
status = nf90_get_var(idcru, idLat, xlat, start=(/1,1/),count=(/nx,ny/))
!status = nf90_get_var(idcru, idLat, xlat)

!---Read in longitude values ----------

status = nf90_inq_varid(idcru, "xlon", idLong)
status = nf90_get_var(idcru, idLong, xlon, start=(/1,1/),count=(/nx,ny/))


if(status /= nf90_noerr) write(*,*)"netcdf error = 1 ",trim(nf90_strerror(status))
!write(*,*) tval(101,101,1:24)
status = nf90_close(idcru)


! on/off - calculate pyaat from monthly temperature data - need to read in monthly temperature
! as 'tval' for this to work
if (.false.) then
write(*,*) "false"
!---Calculate prior year annual average temperature (spring taxa; start in January)
pyaat_sum = -9999.
pyaat = -9999.
do i = 1,nx        !Loop over all longitudes
do ii = 1,ny     !Loop over all latitudes
   cnt_yrs = 1 ! reset count years
   cnt_mnths = 1 ! count of good months (months without missing data)
   do iv = 13,240  ! year loop, starting on July 1 1999  and ending early (2018) 
      if (pyaat_sum(i,ii,cnt_yrs) .eq. -9999. .and. tval(i,ii,iv).lt.100) then
         pyaat_sum(i,ii,cnt_yrs) = 0
         pyaat_sum(i,ii,cnt_yrs) = pyaat_sum(i,ii,cnt_yrs) + tval(i,ii,iv)
         cnt_mnths = cnt_mnths + 1 ! update count of good months
      else if (pyaat_sum(i,ii,cnt_yrs) .ne. -9999. .and. tval(i,ii,iv).lt.100) then
         pyaat_sum(i,ii,cnt_yrs) = pyaat_sum(i,ii,cnt_yrs) + tval(i,ii,iv)
         cnt_mnths = cnt_mnths + 1 ! update count of good months
      else
      	! do nothing
      end if
      if(mod(iv,12).eq.0) then
         if (pyaat_sum(i,ii,cnt_yrs).ne.-9999) then
         	pyaat(i,ii,cnt_yrs) = pyaat_sum(i,ii,cnt_yrs)/cnt_mnths ! to get annual averages
         end if
         cnt_yrs = cnt_yrs + 1 ! update count years
         cnt_mnths = 1 ! reset count months
      end if
   end do ! end time loop
end do
end do

!---Calculate prior year annual average temperature (winter taxa; start in July)
pyaat_sum_adj = -9999.
pyaat_adj = -9999.
do i = 1,nx        !Loop over all longitudes
do ii = 1,ny     !Loop over all latitudes
   cnt_yrs = 1 ! reset count years
   cnt_mnths = 1 ! count of good months (months without missing data)
   do iv = 7,234  ! year loop, starting on July 1 1999  and ending early (2018) 
      if (pyaat_sum_adj(i,ii,cnt_yrs) .eq. -9999. .and. tval(i,ii,iv).lt.100) then
         pyaat_sum_adj(i,ii,cnt_yrs) = 0
         pyaat_sum_adj(i,ii,cnt_yrs) = pyaat_sum_adj(i,ii,cnt_yrs) + tval(i,ii,iv)
         cnt_mnths = cnt_mnths + 1 ! update count of good months
	  end if
      if (pyaat_sum_adj(i,ii,cnt_yrs) .ne. -9999. .and. tval(i,ii,iv).lt.100) then
         pyaat_sum_adj(i,ii,cnt_yrs) = pyaat_sum_adj(i,ii,cnt_yrs) + tval(i,ii,iv)
         cnt_mnths = cnt_mnths + 1
      end if
      if(mod(iv,6).eq.0.and.mod(iv,12).ne.0) then  ! add in extra conditional
         if (pyaat_sum_adj(i,ii,cnt_yrs).ne.-9999.) then
         	pyaat_adj(i,ii,cnt_yrs) = pyaat_sum_adj(i,ii,cnt_yrs)/cnt_mnths ! to get annual averages
         end if
         cnt_yrs = cnt_yrs + 1 ! update count years
         cnt_mnths = 1 ! reset count months
      end if
   end do ! end time loop
end do
end do

end if ! on/off - calculate pyaat from monthy temperature data

!---Calculate Pollen Emissions Start Day of Year, end day of year and Duration-------

if (.true.) then
do i=1,nspec
   if(i.eq.2 .or. i.eq.4 .or. i.eq.7)  then   !Alnus, Cupressaceae, Pinus
     sDOY(:,:,:,i) = NINT((sDOY_m(i)*pyaat_adj(:,:,:))+sDOY_b(i))
     eDOY(:,:,:,i)  = NINT((eDOY_m(i)*pyaat_adj(:,:,:))+eDOY_b(i))
     !     write(*,*) sDOY(69,74,:,i),dur(69,74,:,i)
   !else if(i.eq.2) then !for Alnus, set the same constant for start date and end date
   !  sDOY(:,:,:,i) = NINT((sDOY_m(i)*pyaat_adj(:,:,:))+sDOY_b(i))
   !  eDOY(:,:,:,i)  = NINT((eDOY_m(i)*pyaat_adj(:,:,:))+sDOY_b(i))
   else
     sDOY(:,:,:,i) = NINT((sDOY_m(i)*pyaat(:,:,:))+sDOY_b(i))
     eDOY(:,:,:,i)  = NINT((eDOY_m(i)*pyaat(:,:,:))+eDOY_b(i))
   end if
   dur(:,:,:,i) = eDOY(:,:,:,i) - sDOY(:,:,:,i) !calculate duration from phenological dates
   do ii = 1,nx
      do iii = 1,ny
         do iv = 1,nyr
            if(i.eq.2 .or. i.eq.4 .or. i.eq.7)  then
               if (pyaat_adj(ii,iii,iv) .lt. -9998) then
                  sDOY(ii,iii,iv,i) = -9999.
                  dur(ii,iii,iv,i) = -9999.
                  eDOY(ii,iii,iv,i) = -9999.
               end if
            else
               if (pyaat(ii,iii,iv) .lt. -9998) then
                  sDOY(ii,iii,iv,i) = -9999.
                  dur(ii,iii,iv,i) = -9999.
                  eDOY(ii,iii,iv,i) = -9999.
               end if
            end if
            if (dur(ii,iii,iv,i).le.0) then
               !write(*,*) 'duration is zero or negative'
                dur(ii,iii,iv,i) = -9999. 
            end if
         end do
      end do
    end do
end do
end if

! Due to memory array issues, calculate pollen in species loop, then write netCDF output for each species type
!---Output new netcdf file containing pollen count data for each grid cell for each day (2003-2010)---

! netcdf-output write parameters
  start = (/1,1,1,1/)
  start_all = (/1,1,1,1,1/)
!  count = (/nx,ny,ndy,nyr,1/)   ! note that will set start in the loop below
  count = (/nx,ny,nlev,ndytot/)
  count_all = (/nx,ny,nlev,ndytot,nspec/)
! dimension and variable attributes for netcdf-output
  longname(1) = 'longitude'
  longname(2) = 'latitude'
  longname(3) = 'level index'
  longname(4) = 'time'
  !longname(5) = 'pollen type' 
  longname(5) = 'pollen emission potential (flux)'

  units(1) = 'degrees_east'
  units(2) = 'degrees_north'
  units(3) = 'level'
  units(4) = 'days since 1995-01-01 00:00:00' ! set at beginning of code for convenience!
  units(5) = 'grains m-2 day-1'

! Output NetCDF file create
  rcode = nf90_create(path = "/nfs/alsteine/yingxz/emis_out/Obs/pollen_emissions_2010-2019_25km_speciated-calc_taxon.nc",cmode = nf90_64bit_offset,initialsize = 0,ncid = ncidout)

! NetCDF file variable definition
  rcode = nf90_redef(ncidout)

  rcode = nf90_def_dim(ncidout,'jx', count(1),dimid1)
  rcode = nf90_def_dim(ncidout,'iy', count(2),dimid2)
  rcode = nf90_def_dim(ncidout,'lev',  count(3),dimid3)
  rcode = nf90_def_dim(ncidout,'time', count(4),dimid4)
  rcode = nf90_def_dim(ncidout,'spec', nspec,dimid5)
  rcode = nf90_def_dim(ncidout,'year',nyr,dimid6)
  if(rcode /= nf90_noerr) write(*,*)"netcdf error = 2",trim(nf90_strerror(rcode))
  write(*,*) "error after define dim"
  
  id_2d(1) = dimid1
  id_2d(2) = dimid2
!  rcode = nf90_def_dim(ncidout,'lon', id_2d,dimid6)
!  rcode = nf90_def_dim(ncidout,'lat', id_2d,dimid7)
!  rcode = nf90_def_dim(ncidout,'spec',nf90_unlimited,dimid5)
  write(*,*) "ndytot, nyr 3"
  write(*,*) ndytot, nyr
       ! The dimids array is used to pass the dimids of the dimensions of
       ! the netCDF variables. Both of the netCDF variables we are creating
       ! share the same four dimensions. In Fortran, the unlimited
       ! dimension must come last on the list of dimids.

  dimids = (/ dimid1,dimid2,dimid3,dimid4/) 
  dimids_all = (/ dimid1,dimid2,dimid3,dimid4,dimid5/)  
  phen_dimids = (/ dimid1,dimid2,dimid6,dimid5/) 
  phen_count = (/ nx,ny,nyr,nspec /)

! M. Wozniak edits: changed output variables to those of RegCM chemical emissions format
  rcode = nf90_def_var(ncidout,'lon'   ,nf90_real,id_2d,lonoutid)
  rcode = nf90_def_var(ncidout,'lat'   ,nf90_real,id_2d,latoutid)
  rcode = nf90_def_var(ncidout,'lev'   ,nf90_real,dimid3,levoutid)
  !rcode = nf90_def_var(ncidout,'day'   ,nf90_real,dimid3,dayoutid)
  rcode = nf90_def_var(ncidout,'time'  ,nf90_real,dimid4,timeoutid)
  rcode = nf90_def_var(ncidout,'spec'  ,nf90_real,dimid5,specoutid)
  rcode = nf90_def_var(ncidout,'year'  ,nf90_real,dimid6,yroutid)
  !rcode = nf90_def_var(ncidout,'yr'    ,nf90_real,dimid4,yroutid)
  !rcode = nf90_def_var(ncidout,'nspec' ,nf90_real,dimid5,speoutid)
  !rcode = nf90_def_var(ncidout,'pc'    ,nf90_real,dimids,pcoutid)
  !rcode = nf90_def_var(ncidout,'POLLEN',nf90_real,dimids,polloutid)
  rcode = nf90_def_var(ncidout,'clm_lc',nf90_double,(/dimid1,dimid2/),lcoutid)
  rcode = nf90_def_var(ncidout,'sDOY',nf90_real,phen_dimids,sdoyoutid)
  rcode = nf90_def_var(ncidout,'eDOY',nf90_real,phen_dimids,edoyoutid)
  rcode = nf90_def_var(ncidout,'dur',nf90_real,phen_dimids,duroutid)
  rcode = nf90_def_var(ncidout,'DBL_POLL',nf90_real,dimids,dbloutid)
  rcode = nf90_def_var(ncidout,'ENL_POLL',nf90_real,dimids,enloutid)
  rcode = nf90_def_var(ncidout,'GRA_POLL',nf90_real,dimids,graoutid)
  rcode = nf90_def_var(ncidout,'RAG_POLL',nf90_real,dimids,ragoutid)
  rcode = nf90_def_var(ncidout,'FLUX_ALL',nf90_real,dimids_all,fluxoutid)
  !rcode = nf90_def_var(ncidout,'clm_lc',nf90_double,(/dimid1,dimid2/),polloutid)
  !rcode = nf90_def_var(ncidout,'sDOY',nf90_real,(/dimid1,dimid2,dimid3,dimid4/),sdoyoutid)
  if(rcode /= nf90_noerr) write(*,*)"netcdf error = 3",trim(nf90_strerror(rcode))
  write(*,*) "error after define var"

  rcode = nf90_put_att(ncidout,lonoutid, 'long_name',trim(longname(1)))
  rcode = nf90_put_att(ncidout,latoutid, 'long_name',trim(longname(2)))
  !rcode = nf90_put_att(ncidout,dayoutid, 'long_name',trim(longname(3)))
  !rcode = nf90_put_att(ncidout,yroutid,  'long_name',trim(longname(4)))
  !rcode = nf90_put_att(ncidout,speoutid, 'long_name',trim(longname(5)))
  !rcode = nf90_put_att(ncidout,pcoutid,  'long_name',trim(longname(6)))
  rcode = nf90_put_att(ncidout,timeoutid, 'long_name',trim(longname(4)))
  rcode = nf90_put_att(ncidout,levoutid, 'long_name',trim(longname(3)))
  !rcode = nf90_put_att(ncidout,polloutid, 'long_name',trim(longname(5)))
  rcode = nf90_put_att(ncidout,dbloutid, 'long_name',trim(longname(5)))
  rcode = nf90_put_att(ncidout,enloutid, 'long_name',trim(longname(5)))
  rcode = nf90_put_att(ncidout,graoutid, 'long_name',trim(longname(5)))
  rcode = nf90_put_att(ncidout,ragoutid, 'long_name',trim(longname(5)))
  rcode = nf90_put_att(ncidout,fluxoutid, 'long_name',trim(longname(5)))
  if(rcode /= nf90_noerr) write(*,*)"netcdf error =4 ",trim(nf90_strerror(rcode))
  write(*,*) "error after put att long_name"

  rcode = nf90_put_att(ncidout,lonoutid,'units',trim(units(1)))
  rcode = nf90_put_att(ncidout,latoutid,'units',trim(units(2)))
  !rcode = nf90_put_att(ncidout,dayoutid,'units',trim(units(3)))
  !rcode = nf90_put_att(ncidout,yroutid, 'units',trim(units(4)))
!  rcode = nf90_put_att(ncidout,speoutid,'units',trim(units(5)))
  !rcode = nf90_put_att(ncidout,pcoutid, 'units',trim(units(6)))
  rcode = nf90_put_att(ncidout,timeoutid,'units',trim(units(4)))
  rcode = nf90_put_att(ncidout,levoutid,'units',trim(units(3)))
  !rcode = nf90_put_att(ncidout,polloutid,'units',trim(units(5)))
  rcode = nf90_put_att(ncidout,dbloutid,'units',trim(units(5)))
  rcode = nf90_put_att(ncidout,enloutid,'units',trim(units(5)))
  rcode = nf90_put_att(ncidout,graoutid,'units',trim(units(5)))
  rcode = nf90_put_att(ncidout,ragoutid,'units',trim(units(5)))
  rcode = nf90_put_att(ncidout,fluxoutid,'units',trim(units(5)))
  if(rcode /= nf90_noerr) write(*,*)"netcdf error = 5",trim(nf90_strerror(rcode))
  write(*,*) "error after put att units"

  rcode = nf90_enddef(ncidout)  !end definition mode

! NetCDF file data addition, basic dimension information
  !rcode = nf90_put_var(ncidout,dayoutid,idy,  start(3:3),count(3:3))
  !rcode = nf90_put_var(ncidout,yroutid, iyr,  start(4:4),count(4:4))
  rcode = nf90_put_var(ncidout,lonoutid,xlon, start(1:2),count(1:2))
  rcode = nf90_put_var(ncidout,latoutid,xlat, start(1:2),count(1:2))
  rcode = nf90_put_var(ncidout,timeoutid,idy,  start(4:4),count(4:4))
  rcode = nf90_put_var(ncidout,levoutid,(/1/),start(3:3),count(3:3))
  rcode = nf90_put_var(ncidout,specoutid,(/ (i,i=1,nspec) /),start(4:4),phen_count(4:4))
  rcode = nf90_put_var(ncidout,yroutid,iyr,start(3:3),phen_count(3:3))
  rcode = nf90_put_var(ncidout,lcoutid,frac_area(:,:,13) + frac_area(:,:,14),start(1:2),count(1:2))
  rcode = nf90_put_var(ncidout,sdoyoutid,sDOY,(/1,1,1,1/),phen_count)
  rcode = nf90_put_var(ncidout,sdoyoutid,sDOY,(/1,1,1,1/),phen_count)
  rcode = nf90_put_var(ncidout,edoyoutid,eDOY,(/1,1,1,1/),phen_count)
  rcode = nf90_put_var(ncidout,duroutid,dur,(/1,1,1,1/),phen_count)
  if(rcode /= nf90_noerr) write(*,*)"netcdf error = 6",trim(nf90_strerror(rcode))
  write(*,*) "error after put var"
!---Calculate pollen count-----
  write(*,*) "flux_all"
  write(*,*) size(flux_all)
  write(*,*) shape(flux_all)
wid = 3 ! M. Wozniak edit: gaussian width adjustment (allows one to tune duration of emissions based on count phenology)
!pc = 0.
dbl_poll = 0.
enl_poll = 0.
gra_poll = 0.
rag_poll = 0.
flux_all = 0.
write(*,*) 'starting emissions calculation loop...'
do iv2 = 1,nspec ! -2 cuts off the grasses
plac_flag = .true. ! write out which PFT this taxon's emissions will be added to
cnt = 1 ! now counting days for total period of emissions

! initialize all values to zero
  !pc = 0.   ! reinitialize for next count, comment out if doing sum
  write(*,*) 'Calculating taxon...',iv2,name(iv2)
    do i = 1,nx
       !write(*,*) 'X coord ', i
    do ii = 1,ny
       !write(*,*) 'Y coord ', ii
       cnt = 1
        do iv = yrst,yrend
           !write(*,*) 'Year ', iv 
           !if (year(iv) .eq. '2084' .or. year(iv) .eq. '2088' .or. year(iv) .eq. '2092'.or. year(iv) .eq. '2096' .or. year(iv) .eq. '2100') then
           !if (year(iv) .eq. '2040' .or. year(iv) .eq. '2044' .or. year(iv) .eq. '2048'.or. year(iv) .eq. '2052' .or. year(iv) .eq. '2056') then
           !if (year(iv) .eq. '1996' .or. year(iv) .eq. '2000' .or. year(iv) .eq. '2004' .or. year(iv) .eq. '2008' .or. year(iv) .eq. '2012') then
           !if (year(iv) .eq. '2000' .or. year(iv) .eq. '2004' .or. year(iv) .eq. '2008'.or. year(iv) .eq. '2012' .or. year(iv) .eq. '2016') then
           !if (year(iv) .eq. '1996' .or. year(iv) .eq. '2000' .or. year(iv) .eq. '2004' ) then
           !if (year(iv) .eq. '2008' .or. year(iv) .eq. '2012') then
           !   ndy = 366
           !else
              ndy = 365
           !end if
           tot_flx = 0
           cnt_begin = cnt
           if (dur(i,ii,iv,iv2).gt.0) then
              !write(*,*) 'duration is greater than zero, writing emissions'
           if(.false. .and. iv.lt.nyr) then
                mu1 = sDOY(i,ii,iv,iv2) + dur(i,ii,iv,iv2)/2.
                mu2 = (ndy+sDOY(i,ii,iv+1,iv2)) + dur(i,ii,iv+1,iv2)/2.
                sigma1    = dur(i,ii,iv,iv2)/wid
                sigma2    = dur(i,ii,iv+1,iv2)/wid
                do jj = 1,ndy
                   tot_flx = tot_flx + exp(-(real(jj)-mu1)**2./(2.*sigma1**2))+exp(-(real(jj)-mu2)**2./(2.*sigma2**2))
                end do
                ef_d = ef(iv2)/tot_flx ! daily emissions flux is division of annual emission factor over total yearly flux                                                     
           else
                mu = sDOY(i,ii,iv,iv2) + dur(i,ii,iv,iv2)/2.
                sigma = dur(i,ii,iv,iv2)/wid
                do jj = 1,ndy
                     tot_flx = tot_flx + exp(-(real(jj)-mu)**2./(2.*sigma**2))
                end do
                ef_d = ef(iv2)/tot_flx*2   !set the daily production factor twice as large as before becasue of CO2
           end if
        do iii = 1,ndy
!********************************
!    add precipitation here
!   When there is rainfall, the emission would be 0
!******************************** 
!           write(*,*) iv2,i,ii,iv,iii
          if(.false. .and. iv.lt.nyr) then
             !write(*,*) 'negative sDOY. calculating double gaussian.'
            gphen = exp(-(real(iii)-mu1)**2./(2.*sigma1**2))+exp(-(real(iii)-mu2)**2./(2.*sigma2**2))
           ! pc(i,ii,iv,iii) = ef(iv2)*frac_area(i,ii,iv2)*gphen
            if (iv2.eq.4 .or. iv2.eq.7) then ! ENL taxa
              if(plac_flag) write(*,*) 'adding to ENL PFT, negative sDOY'
              if(plac_flag) write(*,*) 'Year ', iv
               enl_poll(i,ii,1,cnt) = enl_poll(i,ii,1,cnt) + ef_d*frac_area(i,ii,iv2)*gphen
            end if
            if (iv2.eq.13 .or. iv2.eq.14) then ! grass taxa
               if(plac_flag) write(*,*) 'adding to grass PFT, negative sDOY'
               if(plac_flag) write(*,*) 'Year ', iv
               gra_poll(i,ii,1,cnt) = gra_poll(i,ii,1,cnt) + ef_d*frac_area(i,ii,iv2)*gphen
            end if
            if (iv2.eq.15) then ! ragweed
               if(plac_flag) write(*,*) 'adding to ragweed PFT, negative sDOY'
               if(plac_flag) write(*,*) 'Year ', iv
               rag_poll(i,ii,1,cnt) = rag_poll(i,ii,1,cnt) + ef_d*frac_area(i,ii,iv2)*gphen
            end if
            if (iv2.eq.1 .or. iv2.eq.2 .or. iv2.eq.3 .or. iv2.eq.5 .or. iv2.eq.6 .or. iv2.eq.8 .or. iv2.eq.9 .or. iv2.eq.10 .or. iv2.eq.11 .or. iv2.eq.12) then! DBL taxa
               if(plac_flag) write(*,*) 'adding to DBL PFT, negative sDOY'
               if(plac_flag) write(*,*) 'Year ', iv
               dbl_poll(i,ii,1,cnt) = dbl_poll(i,ii,1,cnt) + ef_d*frac_area(i,ii,iv2)*gphen
            end if
            flux_all(i,ii,1,cnt,iv2) = ef_d*frac_area(i,ii,iv2)*gphen
            plac_flag = .false.
            
            !pc(i,ii,1,cnt) = pc(i,ii,1,cnt) + ef_d*frac_area(i,ii,iv2)*gphen*vf
!           if(i.eq.67.and.ii.eq.74) then      ! test for high CUPR grid cell
!               write(*,*) iv,iii,ef(iv2),frac_area(i,ii,iv2),gphen, pc(i,ii,iii,iv)
!           endif
          else
           gphen = exp(-(real(iii)-mu)**2./(2.*sigma**2))
           ! PFT-level pollen counts
           !write(*,*) 'positive sDOY. calculating single gaussian.'
           if (iv2.eq.4 .or. iv2.eq.7)then ! ENL taxa
             if(plac_flag) write(*,*) 'adding to ENL PFT'
             if(plac_flag) write(*,*) 'Year ', iv  
             enl_poll(i,ii,1,cnt) = enl_poll(i,ii,1,cnt) + ef_d*frac_area(i,ii,iv2)*gphen
           end if
           if (iv2.eq.13 .or. iv2.eq.14) then ! grass taxa
               if(plac_flag) write(*,*) 'adding to grass PFT'
               if(plac_flag) write(*,*) 'Year ', iv
               gra_poll(i,ii,1,cnt) = gra_poll(i,ii,1,cnt) + ef_d*frac_area(i,ii,iv2)*gphen
           end if 
           if (iv2.eq.15) then ! ragweed
               if(plac_flag) write(*,*) 'adding to ragweed PFT'
               if(plac_flag) write(*,*) 'Year ', iv
               rag_poll(i,ii,1,cnt) = rag_poll(i,ii,1,cnt) + ef_d*frac_area(i,ii,iv2)*gphen
           end if
           if (iv2.eq.1 .or. iv2.eq.2 .or. iv2.eq.3 .or. iv2.eq.5 .or. iv2.eq.6 .or. iv2.eq.8 .or. iv2.eq.9 .or. iv2.eq.10 .or. iv2.eq.11 .or. iv2.eq.12) then! DBL taxa
              if(plac_flag) write(*,*) 'adding to DBL PFT'
              if(plac_flag) write(*,*) 'Year ', iv
              dbl_poll(i,ii,1,cnt) = dbl_poll(i,ii,1,cnt) + ef_d*frac_area(i,ii,iv2)*gphen
           end if
           flux_all(i,ii,1,cnt,iv2) = ef_d*frac_area(i,ii,iv2)*gphen
           plac_flag = .false.
           ! pc(i,ii,iv,iii) = ef(iv2)*frac_area(i,ii,iv2)*gphen
           !pc(i,ii,1,cnt) = pc(i,ii,1,cnt) + ef_d*frac_area(i,ii,iv2)*gphen*vf
           if (pr(i,ii,cnt) .gt. 5.00) then !When daily precipitation is larger than 5mm/day, set the pollen emission as 0
          	enl_poll(i,ii,1,cnt) = 0.
          	gra_poll(i,ii,1,cnt) = 0.
          	rag_poll(i,ii,1,cnt) = 0.
          	dbl_poll(i,ii,1,cnt) = 0.
          	flux_all(i,ii,1,cnt,iv2) = 0.
           end if
          end if ! if iv.lt.nyr for emissions calculation and placement

          cnt = cnt + 1 ! update the day for pc
         ! if(pc(i,ii,iv,iii) .eq. sqrt(-1) ) then
          !   pc(i,ii,iv,iii) = 0
           !  write(*,*) "I have changed nan to zero."
          !end if
          
        enddo !day loop
           
        else
            cnt_end = cnt_begin+ndy-1
            flux_all(i,ii,1,cnt_begin:cnt_end,iv2) = -9999.
            
            cnt = cnt_end+1
            ! When duration is negative, set the emission as missing value
        end if ! if duration not equal to zero
        enddo ! year loop
    enddo ! latitude loop
    enddo ! longitude loop
! NetCDF file data addition
!  start(5) = iv2 
!  rcode = nf90_put_var(ncidout,pcoutid,pc,start=start,count=count)

!  if(rcode /= nf90_noerr) write(*,*)"netcdf error = ",trim(nf90_strerror(rcode))


  end do          ! end iv2 species loop

  !rcode = nf90_put_var(ncidout,polloutid,pc,start=start,count=count)
  rcode = nf90_put_var(ncidout,dbloutid,dbl_poll,start=start,count=count)
  rcode = nf90_put_var(ncidout,enloutid,enl_poll,start=start,count=count)
  rcode = nf90_put_var(ncidout,graoutid,gra_poll,start=start,count=count)
  rcode = nf90_put_var(ncidout,ragoutid,rag_poll,start=start,count=count)
  rcode = nf90_put_var(ncidout,fluxoutid,flux_all,start=start_all,count=count_all)
  !start1 = (/1,1,1,1/)
  !count1 = (/nx,ny,1,nspec/)
  !rcode = nf90_put_var(ncidout,sdoyoutid,sDOY,start=start1,count=count1)

  if(rcode /= nf90_noerr) write(*,*)"netcdf error = ",trim(nf90_strerror(rcode))

  rcode = nf90_close(ncidout)
  

end program fluxcode