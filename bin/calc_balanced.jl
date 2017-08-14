### calc_balanced.jl
### Author: Annette Foerster (foerster@hawaii.edu)
### 3/2/2017
### julia 0.5
### calls module Thermo.jl

push!(LOAD_PATH, "./")
using Thermo

using NetCDF
using DataStructures
using DataArrays


### -----------  TO MODIFY -------------------------

filein = "./1023_pass1/samurai_RTZ_analysis.nc"
fileout = "./20151023_pass1_balanced_1730.nc"

f = 0.0000585   ### Coriolis parameter
paramBL = false   ### Option to "parameterize" the reference state definition in the boundary layer by setting the mean
                  ### tangential wind speed as constant with height below 2km altitude (taking on the value at 2km)


grid = "RTZ"  ### "XYZ" or "RTZ", gridding of the input file

# Modify the following two lines if grid == XYZ
centerX = 0  ### usually zero, unless storm center is different from domain center
centerY = 0  ### usually zero, unless storm center is different from domain center


# Modify the following four lines if grid == RTZ
lat_center = 18.157      ### latitude of domain center
lon_center = -105.262    ### longitude of domain center
x_rtz = -75.0:1.0:75.0  ### horizontal dimensions of XYZ output grid
y_rtz = -75.0:1.0:75.0  ### horizontal dimensions of XYZ output grid


### ----------- END MODIFY -------------------------


### Read in data

dbz = ncread(filein,"DBZ");
u = ncread(filein,"U");
v = ncread(filein,"V");
ww = ncread(filein,"W");
p = ncread(filein,"P");
T = ncread(filein,"T");
theta = ncread(filein,"THETA");
qv = ncread(filein,"QV");
rhoa = ncread(filein,"RHOA");

if grid == "XYZ"
    x = ncread(filein,"x");
    y = ncread(filein,"y");
    z = ncread(filein,"altitude");
    lat = ncread(filein,"latitude");
    lon = ncread(filein,"longitude");
elseif grid == "RTZ"
    radius = ncread(filein,"radius");
    azimuth = ncread(filein,"theta");
    z = ncread(filein,"altitude");
    dvdz = ncread(filein,"DVDZ");
else
  println("Can't read input file. Exit")
  exit()
end

###Calc aux variables
if grid == "XYZ"

  ### reduce domain size, so that azimuthal average of pressure and temperature
  ### are calculated using data from almost an entire ring everywhere in the domain
  offset_from_boundary = convert(Int64,floor(0.1*length(x)))   ### in gridpoints, 0.1 somewhat arbitrary, corresponds to 60km for x=75km

  x = x[offset_from_boundary+1:end-offset_from_boundary]
  y = y[offset_from_boundary+1:end-offset_from_boundary]
  lat = lat[offset_from_boundary+1:end-offset_from_boundary]
  lon = lon[offset_from_boundary+1:end-offset_from_boundary]

  dbz = dbz[offset_from_boundary+1:end-offset_from_boundary,offset_from_boundary+1:end-offset_from_boundary,:,:]
  u = u[offset_from_boundary+1:end-offset_from_boundary,offset_from_boundary+1:end-offset_from_boundary,:,:]
  v = v[offset_from_boundary+1:end-offset_from_boundary,offset_from_boundary+1:end-offset_from_boundary,:,:]
  ww = ww[offset_from_boundary+1:end-offset_from_boundary,offset_from_boundary+1:end-offset_from_boundary,:,:]
  p = p[offset_from_boundary+1:end-offset_from_boundary,offset_from_boundary+1:end-offset_from_boundary,:,:]
  T = T[offset_from_boundary+1:end-offset_from_boundary,offset_from_boundary+1:end-offset_from_boundary,:,:]
  theta = theta[offset_from_boundary+1:end-offset_from_boundary,offset_from_boundary+1:end-offset_from_boundary,:,:]
  qv = qv[offset_from_boundary+1:end-offset_from_boundary,offset_from_boundary+1:end-offset_from_boundary,:,:]
  rhoa = rhoa[offset_from_boundary+1:end-offset_from_boundary,offset_from_boundary+1:end-offset_from_boundary,:,:]


  vr,vt,radius,azimuth = uv_to_vrvt(u,v,x,y,centerX,centerY)

  dvtdz = calc_vertical_derivative(vt,z) ./1000.0 .*100000.0  # newly calculated -> this one to use in the future
  dvtdz = DataArray(dvtdz,create_mask_value(dvtdz,-99900.0))
  dvtdz = convert(Array,dvtdz,-999);

  ###Calculate mean state
  vt_mean,radii = calc_mean_single(vt,x,y,z,centerX,centerY);
  dvtdz_mean,radii = calc_mean_single(dvtdz,x,y,z,centerX,centerY);
  p_mean,radii = calc_mean_single(p,x,y,z,centerX,centerY);
  T_mean,radii = calc_mean_single(T,x,y,z,centerX,centerY);
  theta_mean,radii = calc_mean_single(theta,x,y,z,centerX,centerY);
  rhoa_mean,radii = calc_mean_single(rhoa,x,y,z,centerX,centerY);
  dbz_mean,radii = calc_mean_single(dbz,x,y,z,centerX,centerY);
  qv_mean,radii = calc_mean_single(qv,x,y,z,centerX,centerY);

elseif grid == "RTZ"
  vt = v
  vr = u

  dvtdz = dvdz ./100000.0
  dvtdz = DataArray(dvtdz,create_mask_value(dvtdz,-0.00999))
  dvtdz = convert(Array,dvtdz,-999);

  ###Calculate mean state
  vt_mean,radii = calc_mean_single_RTZ(vt,radius,azimuth,z);
  dvtdz_mean,radii = calc_mean_single_RTZ(dvtdz,radius,azimuth,z);
  p_mean,radii = calc_mean_single_RTZ(p,radius,azimuth,z);
  T_mean,radii = calc_mean_single_RTZ(T,radius,azimuth,z);
  theta_mean,radii = calc_mean_single_RTZ(theta,radius,azimuth,z);
  rhoa_mean,radii = calc_mean_single_RTZ(rhoa,radius,azimuth,z);
  dbz_mean,radii = calc_mean_single_RTZ(dbz,radius,azimuth,z);
  qv_mean,radii = calc_mean_single_RTZ(qv,radius,azimuth,z);

else
  println("Encountered problem. Exit")
  exit()
end



### "Parameterize" the boundary layer, i.e. set vt constant and dvtdz to 0 in the BL (lowest 2km)
if paramBL == true
  z_keys = sort([key for (key,value) in vt_mean])
  r_keys = sort([key for (key,value) in vt_mean[z_keys[1]]])
  index_bl_top = findfirst(x->x>=2.0,z_keys)

  for k in 1:1:index_bl_top
     vt_mean[z_keys[k]] = vt_mean[z_keys[index_bl_top]]
     for r in r_keys
        dvtdz_mean[z_keys[k]][r] = 0.0
     end
  end
end

### Convert dictionaries into arrays (initially for easier plotting, not sure if it does anything else)
vt_array,z_array,radii_array = dict_to_array(vt_mean);
dvtdz_array,z_array,radii_array = dict_to_array(dvtdz_mean);
p_array,z_array,radii_array = dict_to_array(p_mean);
T_array,z_array,radii_array = dict_to_array(T_mean);
theta_array,z_array,radii_array = dict_to_array(theta_mean);
rhoa_array,z_array,radii_array = dict_to_array(rhoa_mean);
dbz_array,z_array,radii_array = dict_to_array(dbz_mean);
qv_array,z_array,radii_array = dict_to_array(qv_mean);

### Replace -999s with NaN for easier plotting (otherwise colorbar goes down to -999)
vt_array = DataArray(vt_array,create_mask_value(vt_array,-999.0))
vt_array = convert(Array,vt_array,NaN);

dvtdz_array = DataArray(dvtdz_array,create_mask_value(dvtdz_array,-999.0))
dvtdz_array = convert(Array,dvtdz_array,NaN);

p_array = DataArray(p_array,create_mask_value(p_array,-999.0))
p_array = convert(Array,p_array,NaN);

T_array = DataArray(T_array,create_mask_value(T_array,-999.0))
T_array = convert(Array,T_array,NaN);

theta_array = DataArray(theta_array,create_mask_value(theta_array,-999.0))
theta_array = convert(Array,theta_array,NaN);

rhoa_array = DataArray(rhoa_array,create_mask_value(rhoa_array,-999.0))
rhoa_array = convert(Array,rhoa_array,NaN);

dbz_array = DataArray(dbz_array,create_mask_value(dbz_array,-999.0))
dbz_array = convert(Array,dbz_array,NaN);

qv_array = DataArray(qv_array,create_mask_value(qv_array,-999.0))
qv_array = convert(Array,qv_array,NaN);

### Get profiles to anchor inward integration
if grid == "XYZ"
  index_boundary = findfirst(radii,abs(x[1]))
elseif grid == "RTZ"
   index_boundary = length(radii)
else
  println("Encountered problem. Exit")
  exit()
end

p_profile_boundary = p_array[:,index_boundary]
T_profile_boundary = T_array[:,index_boundary]
theta_profile_boundary = theta_array[:,index_boundary]
qv_profile_boundary = qv_array[:,index_boundary]
dbz_profile_boundary = dbz_array[:,index_boundary]
rhoa_profile_boundary = rhoa_array[:,index_boundary]


### Calc pi and tr for vertical profile
pi_profile,thetarho_profile,qr_profile = calc_profile_vars(z,dbz_profile_boundary,rhoa_profile_boundary,
            T_profile_boundary,theta_profile_boundary,qv_profile_boundary);

### Save profile in ordered dict
pi_profile_dict = OrderedDict()
tr_profile_dict = OrderedDict()
qr_profile_dict = OrderedDict()

for (index,k) in enumerate(z)
    pi_profile_dict[Float64(k)] = pi_profile[index]
    tr_profile_dict[Float64(k)] = thetarho_profile[index]
    qr_profile_dict[Float64(k)] = qr_profile[index]
end

###Calculate balanced state
pibar_balanced, thetarhobar_balanced = calc_balanced_profile(radii[1:index_boundary],
z,vt_mean,dvtdz_mean,pi_profile_dict,tr_profile_dict,f);


### Convert Dict to Array for plotting purposes
pibar_bal = ones(Float64,(length(z),length(radii[1:index_boundary]))) .* -999
trbar_bal = ones(Float64,(length(z),length(radii[1:index_boundary]))) .* -999

for (index1,k) in enumerate(z)
    for (index2,ri) in enumerate(radii[1:index_boundary])
        pibar_bal[index1,index2] = pibar_balanced[Float64(k)][Float64(ri)]
        trbar_bal[index1,index2] = thetarhobar_balanced[Float64(k)][Float64(ri)]
    end
end

pibar_bal = DataArray(pibar_bal,create_mask_value(pibar_bal,-999.0))
pibar_bal = convert(Array,pibar_bal,NaN);
trbar_bal = DataArray(trbar_bal,create_mask_value(trbar_bal,-999.0))
trbar_bal = convert(Array,trbar_bal,NaN);


if grid == "XYZ"
  ### Map back to x,y,z grid
  pibar_xy = similar(u); fill!(pibar_xy,-999)
  trbar_xy = similar(u); fill!(trbar_xy,-999)
  for k in 1:1:length(z)
      pibar_xy[:,:,k] = regrid_r_to_xy(radii[1:index_boundary],pibar_bal[k,:],x,y)
      trbar_xy[:,:,k] = regrid_r_to_xy(radii[1:index_boundary],trbar_bal[k,:],x,y)
  end

  ### Calc horizontal gradients of pibar
  dpibdx, dpibdy = calc_horizontal_derivatives(pibar_xy,x,y);

elseif grid == "RTZ"
  ### Map to x,y,z grid
  pibar_xy = ones(Float64, (length(x_rtz),length(y_rtz),length(z))) .* 999.0
  trbar_xy = ones(Float64, (length(x_rtz),length(y_rtz),length(z))) .* 999.0
  for k in 1:1:length(z)
      pibar_xy[:,:,k] = regrid_r_to_xy(radii[1:index_boundary],pibar_bal[k,:],x_rtz,y_rtz)
      trbar_xy[:,:,k] = regrid_r_to_xy(radii[1:index_boundary],trbar_bal[k,:],x_rtz,y_rtz)
  end

  ### Calc horizontal gradients of pibar
  dpibdx, dpibdy = calc_horizontal_derivatives(pibar_xy,x_rtz,y_rtz);

else
  println("Encountered problem. Exit")
  exit()
end

#### Write to nc file to use for thermo retrieval input
data_out = Dict()

dpibdx = DataArray(dpibdx,create_mask_nan(dpibdx));
data_out["dpibdx"] = convert(Array,dpibdx,-999);

dpibdy = DataArray(dpibdy,create_mask_nan(dpibdy));
data_out["dpibdy"] = convert(Array,dpibdy,-999);

trbar_xy = DataArray(trbar_xy,create_mask_nan(trbar_xy));
data_out["trb"] = convert(Array,trbar_xy,-999);

pibar_xy = DataArray(pibar_xy,create_mask_nan(pibar_xy));
data_out["pib"] = convert(Array,pibar_xy,-999);


if grid == "XYZ"
  data_out["lon"] = ncread(filein,"longitude")[offset_from_boundary+1:end-offset_from_boundary]
  data_out["lat"] = ncread(filein,"latitude")[offset_from_boundary+1:end-offset_from_boundary]
  data_out["x"] = ncread(filein,"x")[offset_from_boundary+1:end-offset_from_boundary]
  data_out["y"] = ncread(filein,"y")[offset_from_boundary+1:end-offset_from_boundary]
  data_out["z"] = ncread(filein,"altitude")
  data_out["u"] = ncread(filein,"U")[offset_from_boundary+1:end-offset_from_boundary,offset_from_boundary+1:end-offset_from_boundary,:,:]
  data_out["v"] = ncread(filein,"V")[offset_from_boundary+1:end-offset_from_boundary,offset_from_boundary+1:end-offset_from_boundary,:,:]
  data_out["w"] = ncread(filein,"W")[offset_from_boundary+1:end-offset_from_boundary,offset_from_boundary+1:end-offset_from_boundary,:,:]
  data_out["dudx"] = ncread(filein,"DUDX")[offset_from_boundary+1:end-offset_from_boundary,offset_from_boundary+1:end-offset_from_boundary,:,:]
  data_out["dvdx"] = ncread(filein,"DVDX")[offset_from_boundary+1:end-offset_from_boundary,offset_from_boundary+1:end-offset_from_boundary,:,:]
  data_out["dwdx"] = ncread(filein,"DWDX")[offset_from_boundary+1:end-offset_from_boundary,offset_from_boundary+1:end-offset_from_boundary,:,:]
  data_out["dudy"] = ncread(filein,"DUDY")[offset_from_boundary+1:end-offset_from_boundary,offset_from_boundary+1:end-offset_from_boundary,:,:]
  data_out["dvdy"] = ncread(filein,"DVDY")[offset_from_boundary+1:end-offset_from_boundary,offset_from_boundary+1:end-offset_from_boundary,:,:]
  data_out["dwdy"] = ncread(filein,"DWDY")[offset_from_boundary+1:end-offset_from_boundary,offset_from_boundary+1:end-offset_from_boundary,:,:]
  data_out["dudz"] = ncread(filein,"DUDZ")[offset_from_boundary+1:end-offset_from_boundary,offset_from_boundary+1:end-offset_from_boundary,:,:]
  data_out["dvdz"] = ncread(filein,"DVDZ")[offset_from_boundary+1:end-offset_from_boundary,offset_from_boundary+1:end-offset_from_boundary,:,:]
  data_out["dwdz"] = ncread(filein,"DWDZ")[offset_from_boundary+1:end-offset_from_boundary,offset_from_boundary+1:end-offset_from_boundary,:,:]

elseif grid == "RTZ"

  data_out["x"] = x_rtz
  data_out["y"] = y_rtz
  data_out["z"] = ncread(filein,"altitude")

  ### Add lat/lon
  data_out["lat"],data_out["lon"] = create_lat_lon_arrays(lat_center,lon_center,x_rtz,y_rtz)

  ### Add u,v
  vt = ncread(filein,"V")
  vr = ncread(filein,"U")
  radius = ncread(filein,"radius")
  azimuth = ncread(filein,"theta")

  u,v,tmp_x,tmp_y = vrvt_to_uv(vr,vt,radius,azimuth)
  w = ncread(filein,"W")

  ### Regrid to XYZ
  data_out["u"] =  regrid_rtz_to_xyz(radius,azimuth,data_out["z"],u,data_out["x"],data_out["y"],BCnan,InterpLinear)
  data_out["v"] =  regrid_rtz_to_xyz(radius,azimuth,data_out["z"],v,data_out["x"],data_out["y"],BCnan,InterpLinear)
  data_out["w"] =  regrid_rtz_to_xyz(radius,azimuth,data_out["z"],w,data_out["x"],data_out["y"],BCnan,InterpLinear)

   #data_out["u"][findin(data_out["u"],NaN)] = -999.0
   #data_out["v"][findin(data_out["v"],NaN)] = -999.0
   #data_out["w"][findin(data_out["w"],NaN)] = -999.0

   data_out["u"] = DataArray(data_out["u"],create_mask_nan(data_out["u"]));
   data_out["u"] = convert(Array,data_out["u"],-999);
   data_out["v"] = DataArray(data_out["v"],create_mask_nan(data_out["v"]));
   data_out["v"] = convert(Array,data_out["v"],-999);
   data_out["w"] = DataArray(data_out["w"],create_mask_nan(data_out["w"]));
   data_out["w"] = convert(Array,data_out["w"],-999);


### Calculate spatial derivatives of all three wind components
  data_out["dudx"], data_out["dudy"] = calc_horizontal_derivatives(data_out["u"],data_out["x"],data_out["y"])
  data_out["dudz"] = calc_vertical_derivative(data_out["u"],data_out["z"])
  data_out["dvdx"], data_out["dvdy"] = calc_horizontal_derivatives(data_out["v"],data_out["x"],data_out["y"])
  data_out["dvdz"] = calc_vertical_derivative(data_out["v"],data_out["z"])
  data_out["dwdx"], data_out["dwdy"] = calc_horizontal_derivatives(data_out["w"],data_out["x"],data_out["y"])
  data_out["dwdz"] = calc_vertical_derivative(data_out["w"],data_out["z"])

  data_out["x"] = convert(Array,x_rtz)
  data_out["y"] = convert(Array,y_rtz)

else
  println("Encountered problem. Exit")
  exit()
end

varnames = ["u","v","w","dudx","dvdx","dwdx","dudy","dvdy","dwdy",
    "dudz","dvdz","dwdz","dpibdx","dpibdy","trb","pib"]

for varname in varnames
  data_out[varname] = reshape(data_out[varname],(length(data_out["x"]),length(data_out["y"]),length(data_out["z"]),1))
end

write_additional_variables(fileout,data_out["x"],data_out["y"],data_out["z"],data_out["lat"],data_out["lon"],varnames,data_out)
