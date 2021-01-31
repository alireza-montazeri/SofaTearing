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
#include <ctime>


namespace sofa::component::forcefield
{

template<class DataTypes>
FanForceField<DataTypes>::FanForceField()
    : d_force(initData(&d_force, "force", "applied force to all points"))
    , d_randForceMinCoeff(initData(&d_randForceMinCoeff, "randForceMinCoeff", ""))
    , d_randForceMaxCoeff(initData(&d_randForceMaxCoeff, "randForceMaxCoeff", ""))
    , d_randForceCoeffChangeProba(initData(&d_randForceCoeffChangeProba, "randForceCoeffChangeProba", ""))
    // و d_tearThreshold(initData(&d_tearThreshold, "tearThreshold", ""))
{
    // Nothing more is done here
}


template<class DataTypes>
void FanForceField<DataTypes>::init()
{
    m_topology = this->getContext()->getMeshTopology(); // get the mesh topology to access to the points
    m_topology->getContext()->get(triangleFF);
    m_topology->getContext()->get(triangleMod);
    isTear = false;
    stepCnt=0;
    // m_topology->getContext()->get(triangleAlg);
    // m_topology->getContext()->get(triangleGeo);

    // m_randomGenerator.initSeed( (unsigned int)time(NULL) ); // init random generator
    // m_randForceCoeff = 1.0; // init random force coefficient

    Inherit::init(); // call parent class init
}


template<class DataTypes>
void FanForceField<DataTypes>::addForce(const core::MechanicalParams* /*params*/, DataVecDeriv& currentForce, const DataVecCoord& /*currentPosition*/, const DataVecDeriv& /*currentVelocities*/)
{
    // float randProba = m_randomGenerator.random<float>(0, 1);
    // if( randProba < d_randForceCoeffChangeProba.getValue() )
    // {
    //     m_randForceCoeff = m_randomGenerator.random<float>(d_randForceMinCoeff.getValue(), d_randForceMaxCoeff.getValue()); // generating new random force coefficient
    // }

    // sofa::helper::WriteAccessor<core::objectmodel::Data< VecDeriv> > force = currentForce; // create writer on the current force
    // for(int i = 0 ; i < m_topology->getNbPoints() ; i++)
    // {
    //     force[i] += d_force.getValue() * m_randForceCoeff; // Add asked force randomized with coeff
    // }
    stepCnt++;
    if(stepCnt%5==0)
    {
        size_t nbTriangle = m_topology->getNbTriangles();
        Real maxStressOfAll = 0;
        TriangleID maxStressID=0;
        sofa::helper::vector<TriangleID> triangleToBeRemoved;
        sofa::helper::vector<triangleInfo> triangleInfo = *(triangleFF->triangleInfo.beginEdit());
        for ( unsigned int i = 0 ; i < nbTriangle ; ++i)
        {
            if(fabs(triangleInfo[i].maxStress) > maxStressOfAll)
            {
                maxStressOfAll = fabs(triangleInfo[i].maxStress);
                maxStressID = i;
            }
        }

        if(maxStressOfAll > 100)
        {
            triangleToBeRemoved.push_back(maxStressID);
            triangleMod->removeItems(triangleToBeRemoved);
        }
        else
        {
            maxStressID = 0;
        }
        dmsg_warning() << "maxStressOfAll: " << maxStressOfAll << " maxStressID: " <<  maxStressID;
    }
}

template< class DataTypes>
void FanForceField< DataTypes >::handleEvent(core::objectmodel::Event *event)
{
    msg_error() << "test";
    if (dynamic_cast< sofa::simulation::AnimateBeginEvent *>(event))
    {
        msg_error() << "test";
    }
    if (dynamic_cast< sofa::simulation::AnimateEndEvent *>(event))
    {
        // size_t nbTriangle = m_topology->getNbTriangles();
        // Real maxStressOfAll = 0;
        // TriangleID maxStressID=0;
        // sofa::helper::vector<TriangleID> triangleToBeRemoved;
        // sofa::helper::vector<triangleInfo> triangleInfo = *(triangleFF->triangleInfo.beginEdit());
        // for ( unsigned int i = 0 ; i < nbTriangle ; ++i)
        // {
        //     if(fabs(triangleInfo[i].maxStress) > maxStressOfAll)
        //     {
        //         maxStressOfAll = fabs(triangleInfo[i].maxStress);
        //         maxStressID = i;
        //     }
        // }
        // if(maxStressOfAll > 1)
        // {
        //     triangleToBeRemoved.push_back(maxStressID);
        //     triangleMod->removeTrianglesWarning(triangleToBeRemoved);
        // }
        // triangleToBeRemoved.push_back(100);
        // triangleMod->removeTrianglesWarning(triangleToBeRemoved);
        msg_error() << "test";
    }
}

} // namespace sofa::component::forcefield

#endif // SOFA_COMPONENT_FORCEFIELD_FANFORCEFIELD_INL



