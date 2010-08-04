/*
 *  VarDriver3D.h
 *  samurai
 *
 *  Created by Michael Bell on 4/12/08.
 *  Copyright 2008 Michael Bell. All rights reserved.
 *
 */

#ifndef VARDRIVER3D_H
#define VARDRIVER3D_H

#include "VarDriver.h"
#include "BSpline.h"
#include "Observation.h"
#include "CostFunction3D.h"
#include "CostFunctionAnalytic.h"
#include "MetObs.h"
#include "TCcenter.h"
#include <iostream>
#include <vector>
#include <QHash>
#include <QDir>
#include <QList>
#include <QString>

using namespace std;

class VarDriver3D : public VarDriver
{

public:
	VarDriver3D();
	~VarDriver3D();
	// ESMF type calls
	bool initialize();
	bool run();
	bool finalize();
	
private:
	typedef BSplineBase<real> SplineBase;
	typedef BSpline<real> SplineD;
	
	// Common methods
	void preProcessMetObs();
	bool loadMetObs();
	bool loadBGfromFile();
	bool bilinearMish();
	real bilinearField(real xPos, real yPos, int var);
	
	vector<real> x;
	vector<real> y;
	unsigned int* RnumGridpts;
	vector<Observation> obVector;
	int bc;
	double yincr;
	double xincr;
	double zincr;
	double CQTOL;
	int maxIter;
	real zLevel;
	
	vector<real>** BG;
	vector<real>** BGsave;
	
	// Passable variables
	real* bgB;
	real* bgU;
	real* obs;
	const real* iaScalar;
	const real* iaVector;
	const real* jaScalar;
	const real* jaVector;
	real* ia;
	real* ja;
	real imin, imax, jmin, jmax, kmin, kmax;
	int idim;
	int jdim;
	int kdim;
	// Cost Functions
	CostFunction3D* cost3D;

};

#endif