#ifndef DataFormats_ClusterGeometrySoA_interface_ClusterGeometrysSoA_h
#define DataFormats_ClusterGeometrySoA_interface_ClusterGeometrysSoA_h

#include "DataFormats/SoATemplate/interface/SoALayout.h"

GENERATE_SOA_LAYOUT(ClusterGeometrysLayout,
                    SOA_COLUMN(uint32_t, clusterIds),
                    SOA_COLUMN(float, pitchX),
                    SOA_COLUMN(float, pitchY),
                    SOA_COLUMN(float, thickness),
                    SOA_COLUMN(float, tanLorentzAngles))

using ClusterGeometrysSoA = ClusterGeometrysLayout<>;
using ClusterGeometrysSoAView = ClusterGeometrysSoA::View;
using ClusterGeometrysSoAConstView = ClusterGeometrysSoA::ConstView;

#endif  // DataFormats_ClusterGeometrySoA_interface_ClusterGeometrysSoA_h
