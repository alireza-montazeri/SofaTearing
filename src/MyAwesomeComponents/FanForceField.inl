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

    firstCut = true;
    stepCnt=0;
    tearStep=0.00001;
    
    isDraw = false;
    Inherit::init(); // call parent class init
}


template<class DataTypes>
void FanForceField<DataTypes>::addForce(const core::MechanicalParams* /*params*/, DataVecDeriv& currentForce, const DataVecCoord& /*currentPosition*/, const DataVecDeriv& /*currentVelocities*/)
{
    stepCnt++;
    // dmsg_warning() << "debug: " << m_topology->getTrianglesAroundVertex(100);
    if(stepCnt%5==0)
    {
        bool isNeighber = true;
        Real maxStressOfAll = 0;
        TriangleID maxStressTriangleID=0;

        sofa::helper::vector<triangleInfo> triangleInfo = *(triangleFF->triangleInfo.beginEdit());

        // for(unsigned int i=0;i<trianglesAroundLastVertex.size();i++)
        // {
        //     if(fabs(triangleInfo[trianglesAroundLastVertex.at(i)].maxStress) > maxStressOfAll)
        //     {
        //         maxStressOfAll = fabs(triangleInfo[trianglesAroundLastVertex.at(i)].maxStress);
        //         maxStressTriangleID = trianglesAroundLastVertex.at(i);
        //     }
        // }
        
        size_t nbTriangle = m_topology->getNbTriangles();

        if(maxStressOfAll < d_tearThreshold.getValue())
        {
            for (unsigned int i = 0 ; i < nbTriangle ; ++i)
            {
                if(fabs(triangleInfo[i].maxStress) > maxStressOfAll)
                {
                    maxStressOfAll = fabs(triangleInfo[i].maxStress);
                    maxStressTriangleID = i;                        
                }
            }
        }
        if(maxStressOfAll > d_tearThreshold.getValue())
        {
            DataTypes::Coord stressDirection;
            sofa::defaulttype::Vec<3,Real> stressDirectionVec;
            sofa::defaulttype::Vec<3,double> triangleNormal;
            DataTypes::Coord ref;
            

            if(firstCut)
            {
                ref = triangleGeo->computeTriangleCenter(maxStressTriangleID);
            }
            else
            {
                ref = b;
            }

            triangleNormal = triangleGeo->computeTriangleNormal(maxStressTriangleID);
            stressDirection  = triangleInfo[maxStressTriangleID].principalStressDirection;
            stressDirectionVec[0] = (Real)(stressDirection[0]); 
            stressDirectionVec[1] = (Real)(stressDirection[1]); 
            stressDirectionVec[2] = (Real)(stressDirection[2]);            
            sofa::defaulttype::Vec<3,Real> intersectedLine = stressDirectionVec.cross(triangleNormal);
            // Coord d = intersectedLine*2.5;
            // lineVertices.push_back(sofa::defaulttype::Vector3(ref-d));
            // lineVertices.push_back(sofa::defaulttype::Vector3(ref+d));
            Coord p[3];
            triangleGeo->getTriangleVertexCoordinates(maxStressTriangleID, p);
            triangleVertices.push_back(sofa::defaulttype::Vector3(p[0]));
            triangleVertices.push_back(sofa::defaulttype::Vector3(p[1]));
            triangleVertices.push_back(sofa::defaulttype::Vector3(p[2]));

            unsigned int nbPoints=0;
            for(unsigned int i=0;i<3;i++)
            {
                for(unsigned int j=i+1;j<3;j++)
                {
                    Coord edgeDirection = p[j]-p[i];
                    double denum;
                    if((denum=intersectedLine[1]*edgeDirection[0]-intersectedLine[0]*edgeDirection[1])> 0.0001)
                    {
                        double s = (intersectedLine[0]*p[i][1]-intersectedLine[1]*p[i][0]+intersectedLine[1]*ref[0]-intersectedLine[0]*ref[1])/denum;
                        if(s>=0 && s<=1)
                        {
                            double t = (edgeDirection[0]*s+p[i][0]-ref[0])/intersectedLine[0];
                            if((t*intersectedLine[2]+ref[2])-(s*edgeDirection[2]+p[i][2]) < 0.0001)
                            {
                                if(nbPoints == 0)
                                {
                                    a[0] = intersectedLine[0]*t + ref[0];
                                    a[1] = intersectedLine[1]*t + ref[1];
                                    a[2] = intersectedLine[2]*t + ref[2];
                                }
                                else if(nbPoints == 1)
                                {
                                    b[0] = intersectedLine[0]*t + ref[0];
                                    b[1] = intersectedLine[1]*t + ref[1];
                                    b[2] = intersectedLine[2]*t + ref[2];
                                }
                                nbPoints++;
                            }
                        }
                    }                        
                }
            }
            lineVertices.push_back(sofa::defaulttype::Vector3(a));
            lineVertices.push_back(sofa::defaulttype::Vector3(b));
            isDraw = true;
            dmsg_warning() << "nbPoints: " << nbPoints;


            const Triangle t = m_topology->getTriangle(maxStressTriangleID);
            for(unsigned int i=0;i<3;i++)
            {
                TrianglesAroundVertex triAroundVer = m_topology->getTrianglesAroundVertex(t[i]);
                for(unsigned int j=0;j<triAroundVer.size();j++)
                {
                    if(triAroundVer[j] != maxStressTriangleID)
                    {
                        const bool is_tested = false;
                        unsigned int indTest = 0;
                        if(triangleGeo->isPointInsideTriangle( triAroundVer[j], is_tested, a, indTest))
                        {
                            dmsg_warning() << "point a is also on triangle " << triAroundVer[j];
                        }
                        if(triangleGeo->isPointInsideTriangle( triAroundVer[j], is_tested, b, indTest))
                        {
                            dmsg_warning() << "point b is also on triangle " << triAroundVer[j];
                        }
                    }
                }
            }

            // Vector3 minAABB, maxAABB;
            // const bool is_tested = false;
            // unsigned int indTest = 0;
            // bool isPointInTriangle = true;
            // unsigned int cnt = 0;
            // while(isPointInTriangle)
            // {
            //     cnt++;
            //     b[0] = intersectedLine[0]*cnt*tearStep + a[0];
            //     b[1] = intersectedLine[1]*cnt*tearStep + a[1];
            //     b[2] = intersectedLine[2]*cnt*tearStep + a[2];

            //     isPointInTriangle = triangleGeo->isPointInsideTriangle( ind_ta, is_tested, b, indTest);
            //     // triangleGeo->computeTriangleAABB(ind_ta, minAABB, maxAABB);
            //     // for (int j = 0 ; j < 3 ; j++)
            //     // {
            //     //     if (b[j] < minAABB[j] || b[j] > maxAABB[j])
            //     //     {
            //     //         isPointInTriangle = false;
            //     //         break;
            //     //     }
            //     // }
            // }
            // dmsg_warning() << "cnt: " << cnt << msgendl;
            // std::vector<unsigned int> triIndices;
            // triIndices.clear();
            // for ( unsigned int i = 0 ; i < nbTriangle ; i++)
            // {
            //     triangleGeo->computeTriangleAABB(i, minAABB, maxAABB);
            //     bool isBPointInAABB = true;
            //     bool isAPointInAABB = true;
            //     for (int j = 0 ; j < 3 ; j++)
            //     {
            //         if (b[j] < minAABB[j] || b[j] > maxAABB[j])
            //         {
            //             isBPointInAABB = false;
            //         }
            //         if (a[j] < minAABB[j] || a[j] > maxAABB[j])
            //         {
            //             isAPointInAABB = false;
            //         }
            //     }

            //     if (isBPointInAABB)
            //     {
            //         ind_tb = i;
            //         // triIndices.push_back(i);
            //     }
            //     if(isAPointInAABB)
            //     {
            //         ind_ta = i;
            //     }
            // }
            
            // dmsg_warning() << "indexB: " << ind_tb  << " indexA: " << ind_ta << msgendl;
            
            unsigned int a_last = core::topology::BaseMeshTopology::InvalidID;
            unsigned int b_last = core::topology::BaseMeshTopology::InvalidID;
            sofa::helper::vector<sofa::core::topology::TopologyObjectType> topoPath_list;
            sofa::helper::vector<unsigned int> indices_list;
            sofa::helper::vector<Vec<3, double> > coords2_list;
            ind_ta = maxStressTriangleID;
            ind_tb = maxStressTriangleID;

            if(nbPoints == 2)
            {
                bool isPathOk = triangleGeo->computeIntersectedObjectsList(
                            a_last, a, b, ind_ta, ind_tb,
                            topoPath_list, indices_list,
                            coords2_list);

                dmsg_warning() << "isPathOk:  " << isPathOk;
                if(isPathOk)
                {
                    sofa::helper::vector<unsigned int> new_edges;

                    //Split triangles to create edges along a path given as a the list of existing edges and triangles crossed by it.
                    triangleAlg->SplitAlongPath(a_last, a, b_last, b,
                            topoPath_list, indices_list, coords2_list,
                            new_edges, 0.1, 0.25);

                    sofa::helper::vector<unsigned int> new_points;
                    sofa::helper::vector<unsigned int> end_points;
                    bool reachBorder = false;

                    //Duplicates the given edges
                    triangleAlg->InciseAlongEdgeList(new_edges,
                            new_points, end_points, reachBorder);

                    if (!end_points.empty())
                    {
                        a_last = end_points.back();
                        // trianglesAroundLastVertex.clear();     
                        // trianglesAroundLastVertex = m_topology->getTrianglesAroundVertex(a_last);
                    }

                    // firstCut = false;

                    triangleMod->propagateTopologicalChanges();
                    // notify the end for the current sequence of topological change events
                    triangleMod->notifyEndingEvent();

                    triangleMod->propagateTopologicalChanges();
                }
            }
        }
        else
        {
            maxStressTriangleID = 0;
        }
    }
}

template<class DataTypes>
void FanForceField<DataTypes>::draw(const core::visual::VisualParams* vparams)
{
    if(isDraw)
    {
        sofa::defaulttype::RGBAColor color = sofa::defaulttype::RGBAColor(1,0,0,1);
        vparams->drawTool()->drawLines(lineVertices,3,color);
        color = sofa::defaulttype::RGBAColor(0,0,1,1);
        vparams->drawTool()->drawTriangles(triangleVertices,color);
        lineVertices.clear();
        triangleVertices.clear();
        isDraw = false;
    }
}


} // namespace sofa::component::forcefield

#endif // SOFA_COMPONENT_FORCEFIELD_FANFORCEFIELD_INL



