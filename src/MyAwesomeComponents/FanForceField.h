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
* Authors: The SOFA Team and external contributors (see Authors.txt)          *
*                                                                             *
* Contact information: contact@sofa-framework.org                             *
******************************************************************************/
#ifndef SOFA_COMPONENT_FORCEFIELD_FANFORCEFIELD_H
#define SOFA_COMPONENT_FORCEFIELD_FANFORCEFIELD_H

#include <MyAwesomeComponents/config.h>

#include <sofa/core/behavior/ForceField.h>
#include <sofa/helper/RandomGenerator.h>
#include <sofa/defaulttype/VecTypes.h>
#include <sofa/core/topology/BaseMeshTopology.h>
#include <SofaMiscFem/TriangularFEMForceField.h>
#include <sofa/simulation/AnimateBeginEvent.h>
#include <sofa/simulation/AnimateEndEvent.h>
#include <SofaBaseTopology/TriangleSetTopologyModifier.h>
#include <SofaBaseTopology/TriangleSetGeometryAlgorithms.h>
#include <SofaBaseTopology/TriangleSetTopologyAlgorithms.h>
#include <sofa/helper/ColorMap.h>
#include <sofa/defaulttype/RGBAColor.h>
#include <sofa/core/visual/VisualParams.h>

#include <sofa/helper/map.h>


namespace sofa::component::forcefield
{

template<class DataTypes>
class FanForceField : public core::behavior::ForceField<DataTypes>
{
public:
    SOFA_CLASS(SOFA_TEMPLATE(FanForceField, DataTypes), SOFA_TEMPLATE(core::behavior::ForceField, DataTypes));

    typedef core::behavior::ForceField<DataTypes> Inherit;
    typedef typename DataTypes::Deriv Deriv;
    typedef typename DataTypes::VecCoord VecCoord;
    typedef typename DataTypes::VecDeriv VecDeriv;
    typedef core::objectmodel::Data<VecCoord> DataVecCoord;
    typedef core::objectmodel::Data<VecDeriv> DataVecDeriv;

    typedef sofa::component::forcefield::TriangularFEMForceField<defaulttype::Vec3Types>::TriangleInformation triangleInfo;
    typedef typename Coord::value_type   Real;
    typedef core::topology::BaseMeshTopology::TriangleID TriangleID;
    typedef core::topology::BaseMeshTopology::EdgeID EdgeID;
    typedef core::topology::BaseMeshTopology::PointID PointID;
    typedef core::topology::BaseMeshTopology::Triangle Triangle;
    typedef core::topology::BaseMeshTopology::TrianglesAroundVertex TrianglesAroundVertex;
    typedef sofa::helper::fixed_array<EdgeID,3> EdgesInTriangle;
    typedef sofa::core::topology::BaseMeshTopology::Edge Edge;
    typedef sofa::helper::vector<EdgeID>			EdgesAroundVertex;
    
public:
    Data<Real> d_tearThreshold;
    double tearStep;

protected:    
    /// Component constructor
    FanForceField();
    /// Pointer to the current topology
    sofa::core::topology::BaseMeshTopology* m_topology;
    sofa::component::forcefield::TriangularFEMForceField<defaulttype::Vec3Types>* triangleFF;
    sofa::component::topology::TriangleSetTopologyModifier* triangleMod;
    sofa::component::topology::TriangleSetTopologyAlgorithms<defaulttype::Vec3Types>* triangleAlg;
    sofa::component::topology::TriangleSetGeometryAlgorithms<defaulttype::Vec3Types>* triangleGeo;

    sofa::helper::vector<TriangleID> trianglesAroundLastVertex;
    Coord a;
    Coord b;
    TriangleID ind_ta;
    TriangleID ind_tb;

    std::vector<sofa::defaulttype::Vector3> lineVertices;
    std::vector<sofa::defaulttype::Vector3> triangleVertices;
    bool isDraw;

    TrianglesAroundVertex triAroundPa;
    TrianglesAroundVertex triAroundPb;
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

#ifndef SOFA_COMPONENT_FORCEFIELD_FANFORCEFIELD_CPP
extern template class MYAWESOMECOMPONENTS_API FanForceField<defaulttype::Vec3Types>;
#endif

#endif // SOFA_COMPONENT_FORCEFIELD_FANFORCEFIELD_H
