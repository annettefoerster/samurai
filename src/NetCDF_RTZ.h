 /*
 *  NetCDF.h
 *  samurai
 *
 *  Created by Annette Foerster.
 *  Based on example code from netCDF library
 *  Copyright 2013 Michael Bell. All rights reserved.
 *
 */

#ifndef NETCDF_RTZ_H
#define NETCDF_RTZ_H

#include "NetCDF.h"
#include <netcdfcpp.h>
#include <QString>
	 
class NetCDF_RTZ : public NetCDF  
{

public:
	NetCDF_RTZ();
  ~NetCDF_RTZ();
  
  int readNetCDF(const char* filename);
	double getValue(const int &i,const int &j,const int &k,const QString &varName);	
	double getDerivative(const int &i,const int &j,const int &k, const QString &var, const int &der);
	double calc_A(const int &i,const int &j,const int &k);
	double calc_B(const int &i,const int &j,const int &k);
	double calc_C(const int &i,const int &j,const int &k);
	double calc_D(const int &i,const int &j,const int &k);
        double calc_E(const int &i,const int &j,const int &k);  
  
protected:
  int NRADIUS, NTHETA;
  
  float* radius;
	float* theta;
	float* altitude;	
	float* u;
	float* v;
	float* w;
	float* dudr;
	float* dvdr;
	float* dwdr;
	float* dudt;
	float* dvdt;
	float* dwdt;
	float* dudz;
	float* dvdz;
	float* dwdz;
	float* rhoa;
	float* pibar;
	float* thetarhobar;
	float* pip;
	float* thetarhop;
	float* vbar;
	float* vp;
	
	QString varName;  
	
};


#endif
