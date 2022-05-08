/******************************************************************************
*                 SOFA, Simulation Open-Framework Architecture                *
*                    (c) 2006 INRIA, USTL, UJF, CNRS, MGH                     *
*                                                                             *
* This program is free software; you can redistribute it and/or modify it     *
* under the terms of the GNU Lesser General Public License as published by    *
* the Free Software Foundation; either version 2.1 of the License, or (at     *
* your option) any later version.                                             *
*                                                                             *
* This program is distributed in the hope that it will be useful, but WITHOUT *
* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or       *
* FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License *
* for more details.                                                           *
*                                                                             *
* You should have received a copy of the GNU Lesser General Public License    *
* along with this program. If not, see <http://www.gnu.org/licenses/>.        *
*******************************************************************************
* Authors: Alireza Montazeri                                                  *
*                                                                             *
* Contact information: alireza.montazeri9675@gmail.com                        *
******************************************************************************/
#ifndef SOFA_TEARPLUGIN_H
#define SOFA_TEARPLUGIN_H

#include <config.h>

#include <sofa/core/behavior/ForceField.h>
#include <sofa/helper/RandomGenerator.h>
#include <sofa/defaulttype/VecTypes.h>
#include <sofa/core/topology/BaseMeshTopology.h>
#include <SofaMiscFem/TriangularFEMForceField.h>
#include <sofa/simulation/AnimateBeginEvent.h>
#include <sofa/simulation/AnimateEndEvent.h>
#include <SofaBaseTopology/TriangleSetTopologyModifier.h>
#include <SofaBaseTopology/TriangleSetGeometryAlgorithms.h>
// #include <SofaBaseTopology/TriangleSetTopologyAlgorithms.h>
#include <sofa/helper/ColorMap.h>
#include <sofa/type/RGBAColor.h>
#include <sofa/core/visual/VisualParams.h>

#include <sofa/helper/map.h>


namespace sofa::component::forcefield
{

template<class DataTypes>
class TearPlugin : public core::behavior::ForceField<DataTypes>
{
public:
    SOFA_CLASS(SOFA_TEMPLATE(TearPlugin, DataTypes), SOFA_TEMPLATE(core::behavior::ForceField, DataTypes));

    typedef core::behavior::ForceField<DataTypes> Inherit;
    typedef typename DataTypes::Deriv Deriv;
    typedef typename DataTypes::VecCoord VecCoord;
    typedef typename DataTypes::VecDeriv VecDeriv;
    typedef core::objectmodel::Data<VecCoord> DataVecCoord;
    typedef core::objectmodel::Data<VecDeriv> DataVecDeriv;

    typedef sofa::component::forcefield::TriangularFEMForceField<defaulttype::Vec3Types>::TriangleInformation triangleInfo;
    typedef core::topology::BaseMeshTopology::TriangleID TriangleID;
    typedef core::topology::BaseMeshTopology::EdgeID EdgeID;
    typedef core::topology::BaseMeshTopology::PointID PointID;
    typedef core::topology::BaseMeshTopology::Triangle Triangle;
    typedef core::topology::BaseMeshTopology::TrianglesAroundVertex TrianglesAroundVertex;
    typedef core::topology::BaseMeshTopology::EdgesInTriangle EdgesInTriangle;

    typedef sofa::core::topology::BaseMeshTopology::Edge Edge;
    typedef sofa::core::topology::BaseMeshTopology::EdgesAroundVertex EdgesAroundVertex;

    typedef typename DataTypes::Real Real;
    
public:
    Data<Real> d_tearThreshold;
    Data<bool> d_draw;
    double tearStep;

protected:    
    /// Component constructor
    TearPlugin();
    /// Pointer to the current topology
    sofa::core::topology::BaseMeshTopology* m_topology;
    sofa::component::forcefield::TriangularFEMForceField<defaulttype::Vec3Types>* triangleFF;
    sofa::component::topology::TriangleSetTopologyModifier* triangleMod;
    sofa::component::topology::TriangleSetGeometryAlgorithms<defaulttype::Vec3Types>* triangleGeo;

    core::topology::BaseMeshTopology::TrianglesAroundVertex trianglesAroundPointA;
    core::topology::BaseMeshTopology::TrianglesAroundVertex trianglesAroundPointB;
    Coord a;
    Coord b;
    TriangleID ind_ta;
    TriangleID ind_tb;

    std::vector<type::Vector3> pointVertices;
    std::vector<type::Vector3> lineVertices;
    std::vector<type::Vector3> triangleVertices;
    bool shouldDraw;

    TrianglesAroundVertex triAroundPa;
    TrianglesAroundVertex triAroundPb;

    EdgesAroundVertex edgesAroundPa,trueEdgeAroundPa;
    EdgesAroundVertex edgesAroundPb,trueEdgeAroundPb;

    unsigned int lastPa,lastPb;
public:
    /// Init function
    void init() override;

    /// Forces addition for explicit and implicit integration schemes
    virtual void addForce (const core::MechanicalParams* params, DataVecDeriv& currentForce, const DataVecCoord& currentPosition, const DataVecDeriv& currentVelocities) override;

    /// Forces addition for implicit integration schemes
    virtual void addDForce(const core::MechanicalParams* /*mparams*/, DataVecDeriv& /*d_df*/ , const DataVecDeriv& /*d_dx*/) override {}

    virtual SReal getPotentialEnergy(const core::MechanicalParams* /*params*/, const DataVecCoord& /*x*/) const override { return 0; } // Keep it simpl

    void draw(const core::visual::VisualParams* vparams) override;
protected:
    bool firstCut;
    unsigned int stepCnt;
};

} // namespace sofa::component::forcefield

#ifndef SOFA_TEARPLUGIN_CPP
extern template class SOFATEARING_API TearPlugin<defaulttype::Vec3Types>;
#endif

#endif // SOFA_TEARPLUGIN_H
