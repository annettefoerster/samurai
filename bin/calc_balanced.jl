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
import PyPlot
plt = PyPlot


### -----------  TO MODIFY -------------------------

filein = "./samurai_XYZ_analysis.nc"
fileout = "./samurai_balanced.nc"

paramBL = false

centerX = 0
centerY = 0

f = 0.0000585   ### Coriolis parameter


### ----------- END MODIFY -------------------------


### Read in data

dims = ["x","y","altitude","latitude","longitude"]
x,y,z,lat,lon = read_nc_var(filein,dims);

### added dbz instead of qr
vars = ["DBZ","U","V","W","P","T","THETA","QV","RHOA"]
dbz,u,v,ww,p,T,theta,qv,rhoa = read_nc_var(filein,vars);

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


###Calc aux variables
vr,vt,radius,azimuth = uv_to_vrvt(u,v,x,y,centerX,centerY) ## Wspd2 can be included to check if coordinate transformation worked

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
index_boundary = findfirst(radii,abs(x[1]))
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

### Map back to x,y,z grid 
pibar_xy = similar(u); fill!(pibar_xy,-999)
trbar_xy = similar(u); fill!(trbar_xy,-999)
for k in 1:1:length(z)
    pibar_xy[:,:,k] = regrid_r_to_xy(radii[1:index_boundary],pibar_bal[k,:],x,y)
    trbar_xy[:,:,k] = regrid_r_to_xy(radii[1:index_boundary],trbar_bal[k,:],x,y)
end


### Calc horizontal gradients of pibar
dpibdx, dpibdy = calc_horizontal_derivatives(pibar_xy,x,y);

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


varnames = ["u","v","w","dudx","dvdx","dwdx","dudy","dvdy","dwdy",
    "dudz","dvdz","dwdz","dpibdx","dpibdy","trb","pib"]

write_additional_variables(fileout,data_out["x"],data_out["y"],data_out["z"],data_out["lat"],data_out["lon"],varnames,data_out)
