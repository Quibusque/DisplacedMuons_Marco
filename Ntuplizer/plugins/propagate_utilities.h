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
#include "propagation_definitions.h"
typedef std::pair<TrajectoryStateOnSurface, double> TsosPath;

void markMinimalValues(const TMatrixF& matrix, TMatrixF& boolMatrix);

/**
 * @brief Propagate the fts to the cylinder defined by the radius, minZ, and maxZ. Final
 * state is in tsosPath.
 *
 * @param radius the radius of the cylinder
 * @param minZ the minimum z value of the cylinder (final state z must be greater than this)
 * @param maxZ the maximum z value of the cylinder (final state z must be less than this)
 * @param fts the initial state
 * @param propagatorAlong
 * @param propagatorOpposite
 * @param tsosPath the final state (modified by the function)
 * @param checkFinalZ if true, the final state z must be less than maxZ
 * @return true if the propagation was successful, false otherwise
 */
bool propagateToCylinder(Float_t radius, Float_t minZ, Float_t maxZ, const FreeTrajectoryState& fts,
                         const Propagator* propagatorAlong, const Propagator* propagatorOpposite,
                         TsosPath& tsosPath, bool checkFinalZ = true);
/**
 * @brief Propagate the fts to the z plane defined by planeZ. Final state is in tsosPath.
 *
 * @param maxRadius the maximum radial distance (final state r must be less than this)
 * @param planeZ the z value of the plane
 * @param fts the initial state
 * @param propagatorAlong
 * @param propagatorOpposite
 * @param tsosPath the final state (modified by the function)
 * @param checkFinalRadius if true, the final state r must be less than maxRadius
 * @return true if the propagation was successful, false otherwise
 */
bool propagateToZPlane(Float_t maxRadius, Float_t planeZ, const FreeTrajectoryState& fts,
                       const Propagator* propagatorAlong, const Propagator* propagatorOpposite,
                       TsosPath& tsosPath, bool checkFinalRadius = true);
/**
 * @brief Propagate the fts to the optimal surface. Final state is in finalParams. Surface
 * is returned.
 *
 * First attempt is always to barrel, if that fails endcap is attempted.
 *
 * @param fts the initial state
 * @param finalParams the final state (modified by the function)
 * @param magField
 * @param propagatorAlong
 * @param propagatorOpposite
 * @return PropagationSurface the surface that the fts was propagated to
 */
PropagationSurface findAndPropagateToOptimalSurface(FreeTrajectoryState fts,
                                                    GlobalTrajectoryParameters& finalParams,
                                                    const MagneticField* magField,
                                                    const Propagator* propagatorAlong,
                                                    const Propagator* propagatorOpposite);
/**
 * @brief Propagate the fts to the surface at targetVertex. Final state is in tsosPath.
 *
 * The propagation is "forced" i.e. no constraints are applied to the final z while
 * propagating to barrel and no constraints are applied to the final r while propagating
 * to endcap. This is to ensure that the propagation is successful to go to that surface even
 * in edge cases.
 *
 * @param fts the initial state
 * @param tsosPath the final state (modified by the function)
 * @param targetVertex the target vertex where the surface is located
 * @param propagationSurface the kind of surface to propagate to (cylinder or plane)
 * @param magField
 * @param propagatorAlong
 * @param propagatorOpposite
 * @return true if the propagation was successful, false otherwise
 */
bool propagateToSurface(FreeTrajectoryState fts, TsosPath& tsosPath, GlobalPoint targetVertex,
                        const PropagationSurface propagationSurface, const MagneticField* magField,
                        const Propagator* propagatorAlong, const Propagator* propagatorOpposite);

#endif