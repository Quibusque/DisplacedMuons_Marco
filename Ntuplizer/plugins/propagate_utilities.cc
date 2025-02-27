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

void markMinimalValues(const TMatrixF& matrix, TMatrixF& boolMatrix) {
    unsigned int nrows = matrix.GetNrows();
    unsigned int ncols = matrix.GetNcols();

    // Initialize the boolean matrix with false values
    boolMatrix.ResizeTo(nrows, ncols);
    boolMatrix.Zero();

    // Find the minimal value in each column and set it to true in the boolean matrix
    for (unsigned int col = 0; col < ncols; ++col) {
        float minVal = std::numeric_limits<float>::infinity();
        int minRow = -1;
        for (unsigned int row = 0; row < nrows; ++row) {
            if (matrix(row, col) < minVal) {
                minVal = matrix(row, col);
                minRow = row;
            }
        }
        if (minRow != -1) {
            boolMatrix(minRow, col) = 1.0;
        }
    }

    // If a row has more than one true value keep the lowest one and set the others to false
    for (unsigned int row = 0; row < nrows; ++row) {
        std::vector<unsigned int> trueCols;
        for (unsigned int col = 0; col < ncols; ++col) {
            if (boolMatrix(row, col) == 1.0) {
                trueCols.push_back(col);
            }
        }
        if (trueCols.size() > 1) {
            float minVal = std::numeric_limits<float>::infinity();
            unsigned int minCol = -1;
            for (unsigned int col : trueCols) {
                if (matrix(row, col) < minVal) {
                    minVal = matrix(row, col);
                    minCol = col;
                }
            }
            for (unsigned int col : trueCols) {
                if (col != minCol) {
                    boolMatrix(row, col) = 0.0;
                }
            }
        }
    }

    // Need to store the columns and rows that still need to be filled
    //  i.e. empty columns and rows
    std::vector<unsigned int> unfilledCols;
    for (unsigned int col = 0; col < ncols; ++col) {
        bool hasTrue = false;
        for (unsigned int row = 0; row < nrows; ++row) {
            if (boolMatrix(row, col) == 1.0) {
                hasTrue = true;
                break;
            }
        }
        if (!hasTrue) {
            unfilledCols.push_back(col);
        }
    }
    std::vector<unsigned int> filledRows;
    for (unsigned int row = 0; row < nrows; ++row) {
        for (unsigned int col = 0; col < ncols; ++col) {
            if (boolMatrix(row, col) == 1.0) {
                filledRows.push_back(row);
                break;
            }
        }
    }

    // Submatrix for the remaining rows and columns, apply the function recursively
    if (!unfilledCols.empty() && filledRows.size() < nrows) {
        TMatrixF subMatrix(nrows - filledRows.size(), unfilledCols.size());
        TMatrixF subBoolMatrix(nrows - filledRows.size(), unfilledCols.size());

        unsigned int subRow = 0;
        for (unsigned int row = 0; row < nrows; ++row) {
            if (std::find(filledRows.begin(), filledRows.end(), row) == filledRows.end()) {
                unsigned int subCol = 0;
                for (unsigned int col : unfilledCols) {
                    subMatrix(subRow, subCol) = matrix(row, col);
                    subCol++;
                }
                subRow++;
            }
        }

        markMinimalValues(subMatrix, subBoolMatrix);

        // Fill the values in the main matrix with the values from the submatrix
        subRow = 0;
        for (unsigned int row = 0; row < nrows; ++row) {
            if (std::find(filledRows.begin(), filledRows.end(), row) == filledRows.end()) {
                unsigned int subCol = 0;
                for (unsigned int col : unfilledCols) {
                    boolMatrix(row, col) = subBoolMatrix(subRow, subCol);
                    subCol++;
                }
                subRow++;
            }
        }
    }
}

bool propagateToCylinder(Float_t radius, Float_t minZ, Float_t maxZ, const FreeTrajectoryState& fts,
                         const Propagator* propagatorAlong, const Propagator* propagatorOpposite,
                         TsosPath& tsosPath) {
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

    // Final propagated position counts only if it is within the z range
    bool withinZRange =
        tsosPath.first.globalPosition().z() >= minZ && tsosPath.first.globalPosition().z() <= maxZ;
    return withinZRange;
}

bool propagateToZPlane(Float_t maxRadius, Float_t planeZ, const FreeTrajectoryState& fts,
                       const Propagator* propagatorAlong, const Propagator* propagatorOpposite,
                       TsosPath& tsosPath, bool log = false) {
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

    // Final propagated position counts only if within the max radius
    bool withinMaxRadius = tsosPath.first.globalPosition().perp() <= maxRadius;
    return withinMaxRadius;
}

kPropagationSurface propagateToCommonSurface(FreeTrajectoryState genFTS,
                                             FreeTrajectoryState recoFTS, TsosPath& genTsosPath,
                                             TsosPath& recoTsosPath, const MagneticField* magField,
                                             const Propagator* propagatorAlong,
                                             const Propagator* propagatorOpposite) {
    Float_t cylinderRadius = 420.0;
    Float_t maxRadius = 700.0;
    Float_t minZ = -700.0;
    Float_t maxZ = 700.0;

    // Check if the gen and the reco are inside the propagation cylinder (this is not to stop the
    // code but just a check)
    if (genFTS.position().perp() > cylinderRadius || genFTS.position().z() > maxZ ||
        genFTS.position().z() < minZ) {
        std::cout << "Gen states is outside the propagation cylinder, with coordinates\n"
                  << "R = " << genFTS.position().perp() << " and z = " << genFTS.position().z()
                  << std::endl;
    } else {
        std::cout << "Gen states is inside the propagation cylinder" << std::endl;
    }
    if (recoFTS.position().perp() > cylinderRadius || recoFTS.position().z() > maxZ ||
        recoFTS.position().z() < minZ) {
        std::cout << "Reco states is outside the propagation cylinder, with coordinates\n"
                  << "R = " << recoFTS.position().perp() << " and z = " << recoFTS.position().z()
                  << std::endl;
    } else {
        std::cout << "Reco states is inside the propagation cylinder" << std::endl;
    }

    //Sanity check - initial states within the max dimensions of CMS
    Float_t maxCMSCylinderRadius = 800.0;
    Float_t maxCMSZ = 1200.0;
    Float_t minCMSZ = -1200.0;
    if (genFTS.position().perp() > maxCMSCylinderRadius || genFTS.position().z() > maxCMSZ ||
        genFTS.position().z() < minCMSZ) {
        return GEN_OUTSIDE_CMS;
    }
    if (recoFTS.position().perp() > maxCMSCylinderRadius || recoFTS.position().z() > maxCMSZ ||
        recoFTS.position().z() < minCMSZ) {
        return RECO_OUTSIDE_CMS;
    }

    // First try with propagating both to the cylinder
    bool genPropagationGood = propagateToCylinder(cylinderRadius, minZ, maxZ, genFTS,
                                                  propagatorAlong, propagatorOpposite, genTsosPath);
    bool recoPropagationGood = propagateToCylinder(
        cylinderRadius, minZ, maxZ, recoFTS, propagatorAlong, propagatorOpposite, recoTsosPath);
    if (genPropagationGood && recoPropagationGood) {
        return CYLINDER;
    }
    genTsosPath = TsosPath();
    recoTsosPath = TsosPath();

    // If the cylinder propagation fails, try to propagate to the z planes
    Float_t pZ = genFTS.momentum().z();
    bool prioritizeMaxZ = pZ > 0;
    if (prioritizeMaxZ) {
        genPropagationGood = propagateToZPlane(maxRadius, maxZ, genFTS, propagatorAlong,
                                               propagatorOpposite, genTsosPath);
        recoPropagationGood = propagateToZPlane(maxRadius, maxZ, recoFTS, propagatorAlong,
                                                propagatorOpposite, recoTsosPath);
        if (genPropagationGood && recoPropagationGood) {
            return MAX_Z;
        }
        genTsosPath = TsosPath();
        recoTsosPath = TsosPath();
        genPropagationGood = propagateToZPlane(maxRadius, minZ, genFTS, propagatorAlong,
                                               propagatorOpposite, genTsosPath);
        recoPropagationGood = propagateToZPlane(maxRadius, minZ, recoFTS, propagatorAlong,
                                                propagatorOpposite, recoTsosPath);
        if (genPropagationGood && recoPropagationGood) {
            return MIN_Z;
        }
        genTsosPath = TsosPath();
        recoTsosPath = TsosPath();
    } else {
        genPropagationGood = propagateToZPlane(maxRadius, minZ, genFTS, propagatorAlong,
                                               propagatorOpposite, genTsosPath);
        recoPropagationGood = propagateToZPlane(maxRadius, minZ, recoFTS, propagatorAlong,
                                                propagatorOpposite, recoTsosPath);
        if (genPropagationGood && recoPropagationGood) {
            return MIN_Z;
        }
        genTsosPath = TsosPath();
        recoTsosPath = TsosPath();
        genPropagationGood = propagateToZPlane(maxRadius, maxZ, genFTS, propagatorAlong,
                                               propagatorOpposite, genTsosPath);
        recoPropagationGood = propagateToZPlane(maxRadius, maxZ, recoFTS, propagatorAlong,
                                                propagatorOpposite, recoTsosPath);
        if (genPropagationGood && recoPropagationGood) {
            return MAX_Z;
        }
    }

    return NONE;
}