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
#include "CostFunctionXYZ.h"
#include "CostFunctionXYP.h"
#include "CostFunctionRTZ.h"
#include "MetObs.h"
#include "FrameCenter.h"
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
	virtual bool initialize(const QDomElement& configuration);
	virtual bool run();
	virtual bool finalize();
        virtual bool validateXMLconfig();
	
protected:
	typedef BSplineBase<real> SplineBase;
	typedef BSpline<real> SplineD;

	// Common methods
	bool preProcessMetObs();
	bool loadMetObs();
	bool loadBGfromFile();
	int loadBackgroundObs();
	int loadBackgroundCoeffs();
	bool adjustBackground(const int& bStateSize);
	void updateAnalysisParams(const int& iteration);

	QList<real> bgIn;
	QList<Observation> obVector;
	int maxIter;

	// Cost Functions
	CostFunction3D* obCost3D;
	CostFunction3D* bgCost3D;

	// Variables passed to Cost function
	real* bgB;
	real* bgU;
	real* bgWeights;
	real* obs;
	real* bgObs;
	real imin, imax, jmin, jmax, kmin, kmax;
	real iincr;
	real jincr;
	real kincr;
	int idim;
	int jdim;
	int kdim;
    int runMode;

    enum RunModes {
        XYZ,
        RTZ
    };
};

#endif
