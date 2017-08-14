### Helper functions for calc_balanced.jl
### Author: Annette Foerster (foerster@hawaii.edu)
### 3/2/2017
### julia 0.5


module Thermo

import NetCDF
using NetCDF
using DataStructures
using DataArrays
using Grid


export  uv_to_vrvt,vrvt_to_uv, calc_mean_single,calc_mean_single_RTZ, calc_profile_vars, calc_balanced_profile,
        calc_horizontal_derivatives, calc_vertical_derivative, dict_to_array, read_nc_var,
        write_additional_variables, create_mask_nan, create_mask_value, regrid_r_to_xy, regrid_rtz_to_xyz, BCnan, InterpLinear, CoordInterpGrid

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function uv_to_vrvt(u,v,x,y,centerX,centerY)
###Convert u and v into vr and vt
  println("Convert to vr/vt ...")

  d1,d2,d3,d4 = size(u)

  # Initialize
  vr      = similar(u); fill!(vr, -999)
  vt      = similar(u); fill!(vt, -999)
  radius  = similar(u); fill!(radius, -999)
  azimuth = similar(u); fill!(azimuth, -999)

  for k in 1:1:d3
    for j in 1:1:d2
      for i in 1:1:d1
        xdist = x[i] - centerX
        ydist = y[j] - centerY
        radius[i,j,k]  = sqrt(xdist^2+ydist^2)
        azimuth[i,j,k] = atan2(ydist,xdist)
      if radius[i,j,k] > 0
        if u[i,j,k] != -999
          sn =  ydist/radius[i,j,k]
          cs =  xdist/radius[i,j,k]
        vr[i,j,k] =  u[i,j,k]*cs + v[i,j,k]*sn
        vt[i,j,k] = -u[i,j,k]*sn + v[i,j,k]*cs

      end
      end
      end
    end
  end
  return vr,vt,radius,azimuth
end

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function vrvt_to_uv(vr,vt,radius,azimuth)
###Convert vr and vt into u and v
  println("Convert to u/v ...")

  d1,d2,d3,d4 = size(vr)

  # Initialize
  u    = similar(vr); fill!(u, -999)
  v    = similar(vr); fill!(v, -999)
  x    = similar(vr); fill!(x, -999)
  y    = similar(vr); fill!(y, -999)

  for k in 1:1:d3
    for j in 1:1:d2
      for i in 1:1:d1
        x[i,j,k] = radius[i] * cos(azimuth[j]./180.0.*pi)
        y[i,j,k] = radius[i] * sin(azimuth[j]./180.0.*pi)

        if vr[i,j,k] != -999
          u[i,j,k] = vr[i,j,k] * cos(azimuth[j]./180.0.*pi)  - vt[i,j,k] * sin(azimuth[j]./180.0.*pi)
          v[i,j,k] = vr[i,j,k] * sin(azimuth[j]./180.0.*pi)  + vt[i,j,k] * cos(azimuth[j]./180.0.*pi)
        end
      end
    end
  end
  return u,v,x,y
end

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function calc_mean_single(vt,x,y,z,centerX,centerY)
  println("Calculate mean state in rz-coordinates ...")
  d1,d2,d3 = size(vt)
  dx = x[2]-x[1]
  dy = y[2]-y[1]
  #ds = int((dx+dy)/2.0)  ### julia v 0.3.12
  ds = convert(Int64,round((dx+dy)/2.0))
  if ds==0
    ds = 1
  end

  radius  = similar(vt); fill!(radius, -999)
  azimuth = similar(vt); fill!(azimuth, -999)
  maxradius = 0

  for k in 1:d3
    for j in 1:d2
      for i in 1:d1
        xdist = x[i] - centerX
        ydist = y[j] - centerY
        radius[i,j,k]  = sqrt(xdist^2+ydist^2)
        azimuth[i,j,k] = atan2(ydist,xdist)
        if radius[i,j,k] > maxradius
          maxradius = radius[i,j,k]
        end
      end
    end
  end

  ###Azimuthal mean
  #radii = [1:int(maxradius/ds)]*ds   #### julia v 0.3.12
  radii = ds:convert(Int64,round(maxradius/ds))*ds

  vt_mean = OrderedDict()
  for k in z
    k_idx = convert(Float64,k)
    vt_mean[k_idx] = OrderedDict()
    for r in radii
      r_idx = convert(Float64,r)
      vt_mean[k_idx][r_idx] = -999
    end
  end

  for k in 1:d3
    for r in radii
      sum_vt = 0
      hwt = 0
      for j in 1:d2
        for i in 1:d1
          raddist = abs(radius[i,j,k]-r)
          if (raddist < ds)
            weight = ds-raddist
            if (vt[i,j,k] != -999 && !isnan(vt[i,j,k]))     # need to handle missing data better!
              hwt += weight
              sum_vt += vt[i,j,k]*weight
            end
          end
        end
      end

      if hwt > 0
        vt_mean[convert(Float64,z[k])][convert(Float64,r)] = sum_vt/hwt
      end
    end
  end

  return vt_mean,radii
end
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function calc_mean_single_RTZ(vt,radius,azimuth,z)
  println("Calculate mean state in rz-coordinates ...")

  d1,d2,d3 = size(vt)

###Azimuthal mean
  radii = radius

  vt_mean = OrderedDict()
  for k in z
    k_idx = convert(Float64,k)
    vt_mean[k_idx] = OrderedDict()
    for r in radii
      r_idx = convert(Float64,r)
      vt_mean[k_idx][r_idx] = -999
    end
  end

  for k in 1:d3
    for i in 1:d1
      r = radii[i]
      sum_vt = 0
      hwt = 0
      for j in 1:d2
        if (vt[i,j,k] != -999 && !isnan(vt[i,j,k]))     # need to handle missing data better!
          hwt += 1
          sum_vt += vt[i,j,k]
        end
      end

      if hwt > 0
        vt_mean[convert(Float64,z[k])][convert(Float64,r)] = sum_vt/hwt
      end

    end
  end

  return vt_mean,radii
end

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function calc_profile_vars(z,DBZ,RHOA,T,THETA,QV)
    ### Based on thermo/calc_aux_vars_dbz
    ### Modified to work for a 1D array instead of a 3D array
    println("Calculate auxiliary profile variables ...")
    d3 = length(DBZ)

  pi       = similar(DBZ); fill!(pi, -999)
  thetarho = similar(DBZ); fill!(thetarho, -999)
  qr       = similar(DBZ); fill!(qr, -999)

  for k in 1:1:d3
    height = z[k]
    dbz = DBZ[k]
    dbz = dbz>50.0 ? 50.0 : DBZ[k]
    zz=10.0^(dbz*0.1)
    rainmass = (zz / 14630.0)^(0.6905)
    icemass =  (zz / 670.0)^(0.5587)
    mixed_dbz = 20.0
    rain_dbz = 30.0
    hlow = 5.0
    melting_zone = 1.0
    hhi = hlow + melting_zone
    if dbz > mixed_dbz && dbz <= rain_dbz
      weightr=(dbz-mixed_dbz)/(rain_dbz - mixed_dbz)
      weights=1.-weightr
      icemass=(rainmass*weightr+icemass*weights)/(weightr+weights)
    elseif dbz > 30
      icemass=rainmass
    end
    precipmass = rainmass*(hhi-height)/melting_zone + icemass*(height-hlow)/melting_zone
    height < hlow ? precipmass = rainmass : nothing
    height > hhi  ? precipmass = icemass  : nothing

    qr[k] = precipmass/RHOA[k]
    pi[k] = T[k]/THETA[k]
    thetarho[k] = THETA[k] * (1 + QV[k]/1000.0)/(1 + QV[k]/1000.0 + qr[k]/1000.0)
  end
  return pi,thetarho,qr
end

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function calc_balanced_profile(radii,alt,vt_mean::OrderedDict,dvtdz_mean::OrderedDict,pi_mean::OrderedDict,tr_mean::OrderedDict,f::Number)
  ### In contrast to calc_balanced, this uses a 1-d vector for the environmental pressure and thetarhobar profiles
  ### Can handle input containing -999s
  println("Calculated balanced state ...")

  shape0 = size(radii)[1]
  shape2 = size(alt)[1]
  dz = 1000.0 * (alt[2] - alt[1])
  dr = 1000.0 * (radii[2] - radii[1])


  pibar = Dict()
  thetarhobar = Dict()
  for k in alt
    k_idx = Float64(k)
    pibar[k_idx] = OrderedDict()
    thetarhobar[k_idx] = OrderedDict()
    for r in radii
      r_idx = Float64(r)
      pibar[k_idx][r_idx] = -999
      thetarhobar[k_idx][r_idx] = -999
    end
  end

  for i in 1:(shape0-1)
    for k in 1:(shape2-1)
      pi_z = 0.0
      lnthetarho = 0.0
      nan_flag = false

      for ii in i:shape0
        z = (1000.0 * (alt[k]-alt[1]) + pi_z)/dz +1    # plus one because julia indexing starts with 1 instead of 0
        kk = convert(Int64,floor(z))
        dk = z-kk
        ### Check if any -999 in the integration path
        if (vt_mean[Float64(alt[kk])][Float64(radii[ii])]==-999  || vt_mean[Float64(alt[kk+1])][Float64(radii[ii])]==-999 ||
            dvtdz_mean[Float64(alt[kk])][Float64(radii[ii])]==-999  || dvtdz_mean[Float64(alt[kk+1])][Float64(radii[ii])]==-999)
               nan_flag = true
               break
        else
          vbar = (1.0 - dk)* vt_mean[Float64(alt[kk])][Float64(radii[ii])] + (dk)* vt_mean[Float64(alt[kk+1])][Float64(radii[ii])]
          dvbardz = (1.0 - dk)*dvtdz_mean[Float64(alt[kk])][Float64(radii[ii])] + (dk)*dvtdz_mean[Float64(alt[kk+1])][Float64(radii[ii])]
          radius = 1000.0 * radii[ii]

          if (ii==i || ii==shape0)
            pi_z += 0.5*dr*(vbar*vbar/radius + f*vbar)/9.81
            lnthetarho += 0.5*dr*1.0e-5*dvbardz*(2*vbar/radius + f)/9.81
          else
            pi_z += dr*(vbar*vbar/radius + f*vbar)/9.81
            lnthetarho += dr*1.0e-5*dvbardz*(2*vbar/radius + f)/9.81
          end
        end
      end

      if nan_flag
        pibar[Float64(alt[k])][Float64(radii[i])]  = NaN
        thetarhobar[Float64(alt[k])][Float64(radii[i])]  = NaN

      else
        z = (1000.0 * (alt[k]-alt[1]) + pi_z)/dz +1
        kk = convert(Int64,floor(z))
        dk = z-kk
        ii = shape0
        #abs(pi_z) > dz ? println("Error! DeltaZ larger than dz, BoundsError very likely!") : nothing
        pibar[Float64(alt[k])][Float64(radii[i])]       = (1.0 - dk)*pi_mean[Float64(alt[kk])] + (dk)*pi_mean[Float64(alt[kk+1])]
        thetarhobar[Float64(alt[k])][Float64(radii[i])] =  exp(log((1.0 - dk)*tr_mean[Float64(alt[kk])] + (dk)*tr_mean[Float64(alt[kk+1])]) - lnthetarho)
      end
    end
  end
  return pibar, thetarhobar
end


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function calc_horizontal_derivatives(pibar,x,y)
  println("Calculate horizontal derivatives ...")

  dpibardx = similar(pibar); fill!(dpibardx, -999.0)
  dpibardy = similar(pibar); fill!(dpibardy, -999.0)

  #d1,d2,d3,d4 = size(pibar)
  d1,d2,d3 = size(pibar)
  dx = (x[2]-x[1])
  dy = (y[2]-y[1])

  for k in 1:d3
    for j in 1:d2
      for i in 1:d1
        if i == 1
          if pibar[i,j,k] != -999.0 && pibar[i+1,j,k] != -999.0
            dpibardx[i,j,k] = (pibar[i+1,j,k]-pibar[i,j,k])/dx
          end
        elseif i == d1
          if pibar[i,j,k] != -999.0 && pibar[i-1,j,k] != -999.0
            dpibardx[i,j,k] = (pibar[i,j,k]-pibar[i-1,j,k])/dx
          end
        else
          if pibar[i-1,j,k] != -999.0 && pibar[i+1,j,k] != -999.0
            dpibardx[i,j,k] = (pibar[i+1,j,k]-pibar[i-1,j,k])/(2*dx)
          end
        end

        if j == 1
          if pibar[i,j,k] != -999.0 && pibar[i,j+1,k] != -999.0
            dpibardy[i,j,k] = (pibar[i,j+1,k]-pibar[i,j,k])/dy
          end
        elseif j == d2
          if pibar[i,j,k] != -999.0 && pibar[i,j-1,k] != -999.0
            dpibardy[i,j,k] = (pibar[i,j,k]-pibar[i,j-1,k])/dy
          end
        else
          if pibar[i,j-1,k] != -999.0 && pibar[i,j+1,k] != -999.0
            dpibardy[i,j,k] = (pibar[i,j+1,k]-pibar[i,j-1,k])/(2*dy)
          end
        end
      end
    end
  end


  return dpibardx,dpibardy
end

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function calc_vertical_derivative(pibar,z)
  println("Calculate vertical derivative ...")

  dpibardz = similar(pibar); fill!(dpibardz, -999.0)
  d1,d2,d3 = size(pibar)
  dz = (z[2]-z[1])

  for k in 1:d3
    for j in 1:d2
      for i in 1:d1
        if k == 1
          if pibar[i,j,k] != -999.0 && pibar[i,j,k+1] != -999.0
            dpibardz[i,j,k] = (pibar[i,j,k+1]-pibar[i,j,k])/dz
          end
        elseif k == d3
          if pibar[i,j,k] != -999.0 && pibar[i,j,k-1] != -999.0
            dpibardz[i,j,k] = (pibar[i,j,k]-pibar[i,j,k-1])/dz
          end
        else
          if pibar[i,j,k-1] != -999.0 && pibar[i,j,k+1] != -999.0
            dpibardz[i,j,k] = (pibar[i,j,k+1]-pibar[i,j,k-1])/(2*dz)
          end
        end
      end
    end
  end


  return dpibardz
end

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dict_to_array(vt_mean::OrderedDict)
    z = []; r = []
    for (key,value) in vt_mean
        append!(z,key)
    end

    for (key,value) in vt_mean[z[1]]
        append!(r,key)
    end

    vt_array = zeros(Float64,(length(z),length(r)))
    for i in 1:1:length(z)
        for j in 1:1:length(r)
            vt_array[i,j] = vt_mean[z[i]][r[j]]
        end
    end
    return vt_array,z,r
end

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function read_nc_var(filename::String,varnames)
###Read in Cartesian nc-file
  #println("Read in nc file ...")
  if typeof(varnames) <: String
    res = ncread(filename, string(varnames))
    return res
  elseif typeof(varnames) <: Array
    data = OrderedDict()
    for varname in varnames
      data[string(varname)] = ncread(filename, string(varname))
    end
    return collect(values(data))
  else
    println("Can't read file ... exit()")
    exit()
  end
end

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function write_additional_variables(filename,x,y,z,lat,lon,varnames,datadict)
###Write to nc-file
    ### data has to be in datadict["varname"]
  println("Write to nc file ...")

  ncvars = NcVar[]
  xatts = Dict("long_name" => "x (longitude)", "units" => "km", "missing_value" => -999.0, "_FillValue" => -999.0)
  yatts = Dict("long_name" => "y (latitude)",  "units" => "km", "missing_value" => -999.0, "_FillValue" => -999.0)
  zatts = Dict("long_name" => "z (altitude)",  "units" => "km", "missing_value" => -999.0, "_FillValue" => -999.0)
  tatts = Dict("long_name" => "time",  "units" => "s", "missing_value" => -999.0, "_FillValue" => -999.0)
  lonatts = Dict("long_name" => "longitude", "units" => "degrees_east", "missing_value" => -999.0, "_FillValue" => -999.0)
  latatts = Dict("long_name" => "latitude", "units" => "degrees_north", "missing_value" => -999.0, "_FillValue" => -999.0)

  t = [1.0]
  x_dim = NcDim("x",x,xatts)
  y_dim = NcDim("y",y,yatts)
  z_dim = NcDim("z",z,zatts)
  t_dim = NcDim("t",t,tatts)

  lonvar = NetCDF.NcVar("lon",[x_dim],atts=lonatts,t=Float64); push!(ncvars,lonvar);
  latvar = NetCDF.NcVar("lat",[y_dim],atts=latatts,t=Float64); push!(ncvars,latvar);

  for varname in varnames
    atts  = Dict("long_name" => varname, "units" => "???", "missing_value" => -999.0, "_FillValue" => -999.0)
    push!(ncvars,NcVar(varname,[x_dim,y_dim,z_dim,t_dim],atts=atts,t=Float64))
  end

  nc = NetCDF.create(filename,ncvars)

  NetCDF.putvar(nc,"lon",lon);
  NetCDF.putvar(nc,"lat",lat);

  for varname in varnames
        NetCDF.putvar(nc,varname,datadict[varname])
  end

  NetCDF.close(nc)
  return 0
end

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function create_mask_nan(a::Array)
    b =falses(size(a))
    dims = size(a)
    if (dims[end] == 1 && length(dims) >1)
      rank = length(dims)-1
    else
      rank = length(dims)
    end

    if rank == 1
        for i in 1:1:dims[1]
            b[i] = (isnan(a[i]) ? true : false)
        end
    elseif rank == 2
        for i in 1:1:dims[1]
            for j in 1:1:dims[2]
                b[i,j] = (isnan(a[i,j]) ? true : false)
            end
        end
    elseif rank == 3
        for i in 1:1:dims[1]
            for j in 1:1:dims[2]
                for k in 1:1:dims[3]
                    b[i,j,k] = (isnan(a[i,j,k]) ? true : false)
                end
            end
        end
    elseif rank == 4
        for i in 1:1:dims[1]
            for j in 1:1:dims[2]
                for k in 1:1:dims[3]
                  for l in 1:1:dims[4]
                    b[i,j,k,l] = (isnan(a[i,j,k,l]) ? true : false)
                  end
                end
            end
        end
    else
        println("Can't create mask, array has too many dimensions ... error!")
    end

    return b
end

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function create_mask_value(a::Array,value::Number)
    b =falses(size(a))
    dims = size(a)
    rank = length(dims)
    if rank == 1
        for i in 1:1:dims[1]
            b[i] = (a[i]==value ? true : false)
        end
    elseif rank == 2
        for i in 1:1:dims[1]
            for j in 1:1:dims[2]
                b[i,j] = (a[i,j]==value ? true : false)
            end
        end
    elseif rank == 3
        for i in 1:1:dims[1]
            for j in 1:1:dims[2]
                for k in 1:1:dims[3]
                    b[i,j,k] = (a[i,j,k]==value ? true : false)
                end
            end
        end
    elseif rank ==4
        for i in 1:1:dims[1]
            for j in 1:1:dims[2]
                for k in 1:1:dims[3]
                    for l in 1:1:dims[4]
                      b[i,j,k,l] = (a[i,j,k,l]==value ? true : false)
                    end
                end
            end
        end
    else
        println("Can't create mask, array has too many dimensions ... error!")
    end

    return b
end

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function regrid_r_to_xy(r,field1D,x,y)
  ### calculate r values corresponding to the xy values

  ### Interpolate in R-space to the r value of the cartesian grid
  r_range = r[1]:(r[2]-r[1]):r[end]
  ri = CoordInterpGrid(r_range, float(field1D), BCnan, InterpLinear)


  field_interp = ones(Float64,(length(x),length(y))) .* -999
  for (idx1,i) in enumerate(x)
    for (idx2,j) in enumerate(y)
      ### calculate r value corresponding to the xy values
      r_cartesian = sqrt(i.^2 .+ j.^2)
      field_interp[idx1,idx2] = ri[r_cartesian]
    end
  end

  return field_interp
end

end

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function create_lat_lon_arrays(center_lat,center_lon,x,y)
    ### This uses equations from samurai_lineartrack.pl
    lat  = similar(x); fill!(lat, NaN)
    lon  = similar(x); fill!(lon, NaN)

    latrad = center_lat *pi/180.0
    fac_lat = 111.13209-0.56605*cos(2.0*latrad)+0.00012*cos(4.0*latrad)-0.000002*cos(6.0*latrad);
    fac_lon = 111.41513*cos(latrad)-0.09455*cos(3.0*latrad)+0.00012*cos(5.0*latrad);

    for i in 1:1:length(x)
        lon[i] = center_lon + x[i]/fac_lon
        lat[i] = center_lat + y[i]/fac_lat;
    end

    return lat,lon
end

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function regrid_rtz_to_xyz(r,t,z,field,x,y,bc,interp)
    ### Copied from regrid.jl
    field_new = zeros(Float64,(length(x),length(y),length(z)))
    for (index,alt) in enumerate(z)
        field_new[:,:,index] = regrid_rt_to_xy(r,t,field[:,:,index],x,y,bc,interp)
    end
    return field_new
end

function regrid_rt_to_xy(r,t,field,x,y,bc,interp)
  ### Copied from regrid.jl
  ### calculate rt values corresponding to the xy values
  x_2d,y_2d = grid_2d(x,y)
  r_cartesian = sqrt(x_2d.^2 .+ y_2d.^2)
  az_cartesian = atan2(y_2d,x_2d)./pi .*180.0

  for i in 1:1:length(x)
      for j in 1:1:length(y)
          if az_cartesian[i,j] < 0.0
              az_cartesian[i,j] = az_cartesian[i,j] + 360.0
          end
      end
  end

  ### Interpolate in R-Az-space to the r and az values of the cartesian grid
  r_range = r[1]:(r[2]-r[1]):r[end]
  t_range = t[1]:(t[2]-t[1]):t[end]
  field_interp = interp2D_non_regular_grid(r_range,t_range,field,r_cartesian,az_cartesian,bc,interp)
  return field_interp
end

function interp2D_non_regular_grid(x,y,field,x_new,y_new,bc,interp)
  ### x_new and y_new have to be 2D arrays
  ###interp options: InterpNearest, InterpLinear, InterpQuadratic, InterpCubic
  ###bc options: BCnil,BCnan,BCreflect,BCperiodic,BCnearest,BCfill

    field_interp = CoordInterpGrid((x,y), field, bc, interp);
    field_new = zeros(Float64,size(x_new))

    ### Interpolate to new grid
    for i in 1:1:size(x_new)[1]
      for j in 1:1:size(x_new)[2]
        field_new[i,j] = field_interp[x_new[i,j],y_new[i,j]]
      end
    end
    return field_new
end

function grid_2d(x,y)
  ### same as meshgrid
    dim1 = length(x)
    dim2 = length(y)
    x_2d = ones(Float64,(dim1,dim2))
    y_2d = ones(Float64,(dim1,dim2))
    for i in 1:1:dim1
       x_2d[i,:] = x[i]
    end
    for j in 1:1:dim2
        y_2d[:,j] = y[j]
    end
    return x_2d,y_2d
end

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
