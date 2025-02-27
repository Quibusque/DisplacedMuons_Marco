#ifndef PROPAGATE_UTILITIES_H
#define PROPAGATE_UTILITIES_H

#include <limits>
#include <vector>

#include "DataFormats/GeometrySurface/interface/Cylinder.h"
#include "DataFormats/GeometrySurface/interface/Plane.h"
#include "DataFormats/GeometrySurface/interface/Surface.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/GeometryVector/interface/GlobalVector.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "TMatrixF.h"
#include "TrackingTools/GeomPropagators/interface/Propagator.h"
#include "TrackingTools/TrajectoryParametrization/interface/CartesianTrajectoryError.h"
#include "TrackingTools/TrajectoryParametrization/interface/CurvilinearTrajectoryError.h"
#include "TrackingTools/TrajectoryParametrization/interface/GlobalTrajectoryParameters.h"
#include "TrackingTools/TrajectoryState/interface/FreeTrajectoryState.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"

typedef std::pair<TrajectoryStateOnSurface, double> TsosPath;

enum kPropagationSurface {
    RECO_OUTSIDE_CMS = -4,
    GEN_OUTSIDE_CMS = -3,
    OUTSIDE_DELTAR = -2,
    NONE = -1,
    CYLINDER = 1,
    MAX_Z = 2,
    MIN_Z = 3
};

void markMinimalValues(const TMatrixF& matrix, TMatrixF& boolMatrix);

bool propagateToCylinder(Float_t radius, Float_t minZ, Float_t maxZ, const FreeTrajectoryState& fts,
                         const Propagator* propagatorAlong, const Propagator* propagatorOpposite,
                         TsosPath& tsosPath);

bool propagateToZPlane(Float_t maxRadius, Float_t planeZ, const FreeTrajectoryState& fts,
                       const Propagator* propagatorAlong, const Propagator* propagatorOpposite,
                       TsosPath& tsosPath, bool log = false);

kPropagationSurface propagateToCommonSurface(FreeTrajectoryState genFTS,
                                             FreeTrajectoryState recoFTS, TsosPath& genTsosPath,
                                             TsosPath& recoTsosPath, const MagneticField* magField,
                                             const Propagator* propagatorAlong,
                                             const Propagator* propagatorOpposite);

#endif