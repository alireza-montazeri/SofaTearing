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
#ifndef SOFA_COMPONENT_FORCEFIELD_FANFORCEFIELD_INL
#define SOFA_COMPONENT_FORCEFIELD_FANFORCEFIELD_INL

#include <MyAwesomeComponents/FanForceField.h>

namespace sofa::component::forcefield
{

template<class DataTypes>
FanForceField<DataTypes>::FanForceField()
    : d_tearThreshold(initData(&d_tearThreshold, "tearThreshold", ""))
{
    // Nothing more is done here
}


template<class DataTypes>
void FanForceField<DataTypes>::init()
{
    m_topology = this->getContext()->getMeshTopology(); // get the mesh topology to access to the points
    m_topology->getContext()->get(triangleFF);
    m_topology->getContext()->get(triangleMod);
    m_topology->getContext()->get(triangleAlg);
    m_topology->getContext()->get(triangleGeo);

    isTear = false;
    stepCnt=0;
    
    Inherit::init(); // call parent class init
}


template<class DataTypes>
void FanForceField<DataTypes>::addForce(const core::MechanicalParams* /*params*/, DataVecDeriv& currentForce, const DataVecCoord& /*currentPosition*/, const DataVecDeriv& /*currentVelocities*/)
{
    stepCnt++;
    if(stepCnt%5==0)
    {
        Real maxStressOfAll = 0;
        TriangleID maxStressTriangleID=0;
        sofa::helper::vector<TriangleID> triangleToBeRemoved;

        sofa::helper::vector<triangleInfo> triangleInfo = *(triangleFF->triangleInfo.beginEdit());

        for(unsigned int i=0;i<trianglesAroundLastVertex.size();i++)
        {
            if(triangleInfo[trianglesAroundLastVertex.at(i)].maxStress > maxStressOfAll)
            {
                maxStressOfAll = fabs(triangleInfo[trianglesAroundLastVertex.at(i)].maxStress);
                maxStressTriangleID = trianglesAroundLastVertex.at(i);
            }
        }

        if(maxStressOfAll < d_tearThreshold.getValue())
        {
            size_t nbTriangle = m_topology->getNbTriangles();

            for (unsigned int i = 0 ; i < nbTriangle ; ++i)
            {
                if(fabs(triangleInfo[i].maxStress) > maxStressOfAll)
                {
                    maxStressOfAll = fabs(triangleInfo[i].maxStress);
                    maxStressTriangleID = i;
                }
            }
        }
    
        DataTypes::Coord principalStressDirection;
        sofa::defaulttype::Vec<3,double> triangleNormal;

        if(maxStressOfAll > d_tearThreshold.getValue())
        {
            if(trianglesAroundLastVertex.size() == 0)
            {
                a = triangleGeo->computeTriangleCenter(maxStressTriangleID);
            }
            else
            {
                a = b;
            }
            triangleNormal = triangleGeo->computeTriangleNormal(maxStressTriangleID);
            principalStressDirection  = triangleInfo[maxStressTriangleID].principalStressDirection;

            
            
            triangleToBeRemoved.push_back(maxStressTriangleID);
            triangleMod->removeItems(triangleToBeRemoved);
        }
        else
        {
            maxStressTriangleID = 0;
        }
        trianglesAroundLastVertex = m_topology->getTrianglesAroundVertex(10);
        dmsg_warning() << "maxStressOfAll: " << maxStressOfAll << " maxStressID: " <<  trianglesAroundVertex;
    }
}

} // namespace sofa::component::forcefield

#endif // SOFA_COMPONENT_FORCEFIELD_FANFORCEFIELD_INL



