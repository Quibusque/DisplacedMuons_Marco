#include "propagate_utilities.h"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <limits>
#include <memory>
#include <string>
#include <vector>

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/GeometrySurface/interface/Cylinder.h"
#include "DataFormats/GeometrySurface/interface/Plane.h"
#include "DataFormats/GeometrySurface/interface/Surface.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/GeometryVector/interface/GlobalVector.h"
#include "DataFormats/Math/interface/AlgebraicROOTObjects.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "Math/Error.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "TFile.h"
#include "TH1F.h"
#include "TLorentzVector.h"
#include "TMatrixF.h"
#include "TTree.h"
#include "TrackPropagation/SteppingHelixPropagator/interface/SteppingHelixPropagator.h"
#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"
#include "TrackingTools/TrajectoryParametrization/interface/CartesianTrajectoryError.h"
#include "TrackingTools/TrajectoryParametrization/interface/CurvilinearTrajectoryError.h"
#include "TrackingTools/TrajectoryParametrization/interface/GlobalTrajectoryParameters.h"
#include "TrackingTools/TrajectoryState/interface/FreeTrajectoryState.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"

bool propagateToCylinder(Float_t radius, Float_t minZ, Float_t maxZ, const FreeTrajectoryState& fts,
                         const Propagator* propagatorAlong, const Propagator* propagatorOpposite,
                         TsosPath& tsosPath, bool checkFinalZ) {
    const Surface::RotationType dummyRot;
    Cylinder::CylinderPointer theTargetCylinder =
        Cylinder::build(Surface::PositionType(0., 0., 0.), dummyRot, radius);

    // The track should be propagated alongMomentum if the starting point is inside the cylinder
    // and oppositeToMomentum if the starting point is outside the cylinder
    bool isInsideInitial =
        fts.position().perp() < radius && fts.position().z() >= minZ && fts.position().z() <= maxZ;

    const Propagator* selectedPropagator = isInsideInitial ? propagatorAlong : propagatorOpposite;
    tsosPath = selectedPropagator->propagateWithPath(fts, *theTargetCylinder);

    if (!tsosPath.first.isValid()) {
        return false;
    }

    if (checkFinalZ) {
        bool withinZRange = tsosPath.first.globalPosition().z() >= minZ &&
                            tsosPath.first.globalPosition().z() <= maxZ;
        return withinZRange;
    }

    return true;
}

bool propagateToZPlane(Float_t maxRadius, Float_t planeZ, const FreeTrajectoryState& fts,
                       const Propagator* propagatorAlong, const Propagator* propagatorOpposite,
                       TsosPath& tsosPath, bool checkFinalRadius) {
    const Surface::RotationType dummyRot;
    Plane::PlanePointer theTargetPlane =
        Plane::build(Surface::PositionType(0., 0., planeZ), dummyRot);

    // The track should be propagated alongMomentum if going towards the plane
    // and oppositeToMomentum if going away from the plane
    Float_t posZ = fts.position().z();
    Float_t pZ = fts.momentum().z();

    bool isGoingTowardsPlane = (posZ < planeZ && pZ > 0) || (posZ > planeZ && pZ < 0);
    const Propagator* selectedPropagator =
        isGoingTowardsPlane ? propagatorAlong : propagatorOpposite;

    tsosPath = selectedPropagator->propagateWithPath(fts, *theTargetPlane);

    if (!tsosPath.first.isValid()) {
        return false;
    }
    if (checkFinalRadius) {
        bool withinMaxRadius = tsosPath.first.globalPosition().perp() <= maxRadius;
        return withinMaxRadius;
    }
    return true;
}

PropagationSurface findAndPropagateToOptimalSurface(FreeTrajectoryState fts,
                                                    GlobalTrajectoryParameters& finalParams,
                                                    const MagneticField* magField,
                                                    const Propagator* propagatorAlong,
                                                    const Propagator* propagatorOpposite) {
    // Check initial position to be within CMS bounds
    Float_t posZ = fts.position().z();
    Float_t posR = fts.position().perp();
    if (posZ < PropagationConstants::MIN_CMS_Z || posZ > PropagationConstants::MAX_CMS_Z ||
        posR > PropagationConstants::MAX_CMS_CYLINDER_RADIUS) {
        return PropagationConstants::GEN_OUTSIDE_CMS;
    }
    TsosPath tsosPath = TsosPath();
    PropagationSurface optimalSurface = PropagationConstants::GEN_PROPAGATION_FAIL;

    // First try with propagating to cylinder
    Float_t radius = PropagationConstants::CYLINDER.radius;
    Float_t minZ = PropagationConstants::CYLINDER.minZ;
    Float_t maxZ = PropagationConstants::CYLINDER.maxZ;

    bool genPropagationGood =
        propagateToCylinder(radius, minZ, maxZ, fts, propagatorAlong, propagatorOpposite, tsosPath);
    if (genPropagationGood) {
        optimalSurface = PropagationConstants::CYLINDER;
    } else {  // If the cylinder propagation fails, try to propagate to the z planes
        tsosPath = TsosPath();
        Float_t pZ = fts.momentum().z();
        bool prioritizeMaxZ = pZ > 0;
        if (prioritizeMaxZ) {
            Float_t maxRadius = PropagationConstants::POS_ENDCAP.radius;
            Float_t maxZ = PropagationConstants::POS_ENDCAP.maxZ;
            genPropagationGood = propagateToZPlane(maxRadius, maxZ, fts, propagatorAlong,
                                                   propagatorOpposite, tsosPath);
            if (genPropagationGood) {
                optimalSurface = PropagationConstants::POS_ENDCAP;
            } else {
                tsosPath = TsosPath();
                maxRadius = PropagationConstants::NEG_ENDCAP.radius;
                Float_t minZ = PropagationConstants::NEG_ENDCAP.minZ;
                genPropagationGood = propagateToZPlane(maxRadius, minZ, fts, propagatorAlong,
                                                       propagatorOpposite, tsosPath);
                if (genPropagationGood) {
                    optimalSurface = PropagationConstants::NEG_ENDCAP;
                }
            }
        } else {
            Float_t maxRadius = PropagationConstants::NEG_ENDCAP.radius;
            Float_t minZ = PropagationConstants::NEG_ENDCAP.minZ;
            genPropagationGood = propagateToZPlane(maxRadius, minZ, fts, propagatorAlong,
                                                   propagatorOpposite, tsosPath);
            if (genPropagationGood) {
                optimalSurface = PropagationConstants::NEG_ENDCAP;
            } else {
                tsosPath = TsosPath();
                Float_t maxRadius = PropagationConstants::POS_ENDCAP.radius;
                Float_t maxZ = PropagationConstants::POS_ENDCAP.maxZ;
                genPropagationGood = propagateToZPlane(maxRadius, maxZ, fts, propagatorAlong,
                                                       propagatorOpposite, tsosPath);
                if (genPropagationGood) {
                    optimalSurface = PropagationConstants::POS_ENDCAP;
                }
            }
        }
    }

    if (genPropagationGood) {
        finalParams =
            GlobalTrajectoryParameters(tsosPath.first.globalPosition(),
                                       tsosPath.first.globalMomentum(), fts.charge(), magField);
    }

    return optimalSurface;
}

bool propagateToSurface(FreeTrajectoryState fts, TsosPath& tsosPath,
                        const PropagationSurface propagationSurface, const MagneticField* magField,
                        const Propagator* propagatorAlong, const Propagator* propagatorOpposite) {
    if (propagationSurface == PropagationConstants::CYLINDER) {
        Float_t radius = PropagationConstants::CYLINDER.radius;
        Float_t minZ = PropagationConstants::CYLINDER.minZ;
        Float_t maxZ = PropagationConstants::CYLINDER.maxZ;
        bool checkZRange = false;
        return propagateToCylinder(radius, minZ, maxZ, fts, propagatorAlong, propagatorOpposite,
                                   tsosPath, checkZRange);
    } else if (propagationSurface == PropagationConstants::POS_ENDCAP) {
        Float_t maxZ = PropagationConstants::POS_ENDCAP.maxZ;
        Float_t maxRadius = PropagationConstants::POS_ENDCAP.radius;
        bool checkFinalRadius = false;
        return propagateToZPlane(maxRadius, maxZ, fts, propagatorAlong, propagatorOpposite,
                                 tsosPath, checkFinalRadius);
    } else if (propagationSurface == PropagationConstants::NEG_ENDCAP) {
        Float_t minZ = PropagationConstants::NEG_ENDCAP.minZ;
        Float_t maxRadius = PropagationConstants::NEG_ENDCAP.radius;
        bool checkFinalRadius = false;
        return propagateToZPlane(maxRadius, minZ, fts, propagatorAlong, propagatorOpposite,
                                 tsosPath, checkFinalRadius);
    }
    return false;
}
