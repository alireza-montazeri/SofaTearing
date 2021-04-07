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
    if(stepCnt%10==0)
    {
        bool isNeighber = true;
        Real maxStressOfAll = 0;
        TriangleID maxStressTriangleID=0;
        int pointSide=0;

        sofa::helper::vector<triangleInfo> triangleInfo = *(triangleFF->triangleInfo.beginEdit());

        trueEdgeAroundPa.clear();
        trueEdgeAroundPb.clear();

        for(unsigned int i=0;i<edgesAroundPa.size();i++)
        {
            sofa::helper::vector<TriangleID> tri = m_topology->getTrianglesAroundEdge(edgesAroundPa[i]);
            if(tri.size()>1)
            {
                for(unsigned int j=0;j<tri.size();j++)
                {
                    if(fabs(triangleInfo[tri[j]].maxStress) > maxStressOfAll)
                    {
                        maxStressOfAll = fabs(triangleInfo[tri[j]].maxStress);
                        maxStressTriangleID = tri[j];
                        pointSide = 1;
                    }
                }
            }
            trueEdgeAroundPa.push_back(edgesAroundPa[i]);
        }

        for(unsigned int i=0;i<edgesAroundPb.size();i++)
        {
            sofa::helper::vector<TriangleID> tri = m_topology->getTrianglesAroundEdge(edgesAroundPb[i]);
            if(tri.size()>1)
            {
                for(unsigned int j=0;j<tri.size();j++)
                {
                    if(fabs(triangleInfo[tri[j]].maxStress) > maxStressOfAll)
                    {
                        maxStressOfAll = fabs(triangleInfo[tri[j]].maxStress);
                        maxStressTriangleID = tri[j];
                        pointSide = 2;
                    }
                }
            }
            trueEdgeAroundPb.push_back(edgesAroundPb[i]);
        }

        if(firstCut)
        {
            size_t nbTriangle = m_topology->getNbTriangles();
            for (unsigned int i = 0 ; i < nbTriangle ; ++i)
            {
                if(fabs(triangleInfo[i].maxStress) > maxStressOfAll)
                {
                    maxStressOfAll = fabs(triangleInfo[i].maxStress);
                    maxStressTriangleID = i;   
                    pointSide = 0;                   
                }
            }
        }

        if(maxStressOfAll > d_tearThreshold.getValue())
        {
            DataTypes::Coord stressDirection;
            sofa::defaulttype::Vec<3,Real> stressDirectionVec;
            sofa::defaulttype::Vec<3,double> triangleNormal;

            triangleNormal = triangleGeo->computeTriangleNormal(maxStressTriangleID);
            stressDirection  = triangleInfo[maxStressTriangleID].principalStressDirection;
            stressDirectionVec[0] = (Real)(stressDirection[0]); 
            stressDirectionVec[1] = (Real)(stressDirection[1]); 
            stressDirectionVec[2] = (Real)(stressDirection[2]);            
            sofa::defaulttype::Vec<3,Real> intersectedLine = stressDirectionVec.cross(triangleNormal);
            intersectedLine /= sqrt(intersectedLine[0]*intersectedLine[0]+intersectedLine[1]*intersectedLine[1]+intersectedLine[2]*intersectedLine[2]);

            Coord p[3];
            triangleGeo->getTriangleVertexCoordinates(maxStressTriangleID, p);
            triangleVertices.push_back(sofa::defaulttype::Vector3(p[0]));
            triangleVertices.push_back(sofa::defaulttype::Vector3(p[1]));
            triangleVertices.push_back(sofa::defaulttype::Vector3(p[2]));

            double maxDotProduct=0,dotProduct;
            EdgeID firstEdgeId=0, secondEdgeId=0;

            if(pointSide == 1)
            {
                for(unsigned int i=0;i<trueEdgeAroundPa.size();i++)
                {
                    Coord edgePointsPos[2];
                    triangleGeo->getEdgeVertexCoordinates(trueEdgeAroundPa[i], edgePointsPos);
                    Coord edgeDirection = edgePointsPos[1]-edgePointsPos[0];
                    edgeDirection /= sqrt(edgeDirection[0]*edgeDirection[0]+edgeDirection[1]*edgeDirection[1]+edgeDirection[2]*edgeDirection[2]);
                    dotProduct = edgeDirection[0]*intersectedLine[0] + edgeDirection[1]*intersectedLine[1] + edgeDirection[2]*intersectedLine[2];
                    if(fabs(dotProduct) > maxDotProduct)
                    {
                        firstEdgeId = trueEdgeAroundPa[i];
                        maxDotProduct = fabs(dotProduct);
                    }
                }
            }
            else if(pointSide == 2)
            {
                for(unsigned int i=0;i<trueEdgeAroundPb.size();i++)
                {
                    Coord edgePointsPos[2];
                    triangleGeo->getEdgeVertexCoordinates(trueEdgeAroundPb[i], edgePointsPos);
                    Coord edgeDirection = edgePointsPos[1]-edgePointsPos[0];
                    edgeDirection /= sqrt(edgeDirection[0]*edgeDirection[0]+edgeDirection[1]*edgeDirection[1]+edgeDirection[2]*edgeDirection[2]);
                    dotProduct = edgeDirection[0]*intersectedLine[0] + edgeDirection[1]*intersectedLine[1] + edgeDirection[2]*intersectedLine[2];
                    if(fabs(dotProduct) > maxDotProduct)
                    {
                        firstEdgeId = trueEdgeAroundPb[i];
                        maxDotProduct = fabs(dotProduct);
                    }
                }
            }
            else if(pointSide == 0)
            {
                unsigned int ver1,ver2;

                for(unsigned int i=0;i<3;i++)
                {
                    for(unsigned int j=i+1;j<3;j++)
                    {
                        Coord edgeDirection = p[j]-p[i];
                        edgeDirection /= sqrt(edgeDirection[0]*edgeDirection[0]+edgeDirection[1]*edgeDirection[1]+edgeDirection[2]*edgeDirection[2]);
                        dotProduct = edgeDirection[0]*intersectedLine[0] + edgeDirection[1]*intersectedLine[1] + edgeDirection[2]*intersectedLine[2];
                        if(fabs(dotProduct) > maxDotProduct)
                        {
                            ver1 = i; ver2 = j;
                            maxDotProduct = fabs(dotProduct);
                        }                     
                    }
                }

                Triangle t = m_topology->getTriangle(maxStressTriangleID);
                EdgesInTriangle edges = m_topology->getEdgesInTriangle(maxStressTriangleID);
                for(unsigned int i=0;i<3;i++)
                {
                    Edge edge = m_topology->getEdge(edges[i]);
                    if(edge[0] == t[ver1] || edge[0] == t[ver2])
                    {
                        if(edge[1] == t[ver1] || edge[1] == t[ver2])
                        {
                            firstEdgeId = edges[i];
                        }
                    }
                }
            }

            maxDotProduct = 0;
            Edge edge = m_topology->getEdge(firstEdgeId);
            if(pointSide == 1)
            {
                unsigned int p;
                if(edge[0] == lastPa)       p = edge[1];
                else if(edge[1] == lastPa)  p = edge[0];

                EdgesAroundVertex edgesAroundVer = m_topology->getEdgesAroundVertex(p);
                for(unsigned int j=0;j<edgesAroundVer.size();j++)
                {
                    if(edgesAroundVer[j] != firstEdgeId)
                    {
                        Coord edgePointsPos[2];
                        triangleGeo->getEdgeVertexCoordinates(edgesAroundVer[j], edgePointsPos);
                        Coord edgeDirection = edgePointsPos[1]-edgePointsPos[0];
                        edgeDirection /= sqrt(edgeDirection[0]*edgeDirection[0]+edgeDirection[1]*edgeDirection[1]+edgeDirection[2]*edgeDirection[2]);
                        dotProduct = edgeDirection[0]*intersectedLine[0] + edgeDirection[1]*intersectedLine[1] + edgeDirection[2]*intersectedLine[2];
                        if(fabs(dotProduct) > maxDotProduct && m_topology->getTrianglesAroundEdge(edgesAroundVer[j]).size()>1)
                        {
                            secondEdgeId = edgesAroundVer[j];
                            maxDotProduct = fabs(dotProduct);
                        }
                    }
                }
            }
            else if(pointSide == 2)
            {
                unsigned int p;
                if(edge[0] == lastPb)       p = edge[1];
                else if(edge[1] == lastPb)  p = edge[0];

                EdgesAroundVertex edgesAroundVer = m_topology->getEdgesAroundVertex(p);
                for(unsigned int j=0;j<edgesAroundVer.size();j++)
                {
                    if(edgesAroundVer[j] != firstEdgeId)
                    {
                        Coord edgePointsPos[2];
                        triangleGeo->getEdgeVertexCoordinates(edgesAroundVer[j], edgePointsPos);
                        Coord edgeDirection = edgePointsPos[1]-edgePointsPos[0];
                        edgeDirection /= sqrt(edgeDirection[0]*edgeDirection[0]+edgeDirection[1]*edgeDirection[1]+edgeDirection[2]*edgeDirection[2]);
                        dotProduct = edgeDirection[0]*intersectedLine[0] + edgeDirection[1]*intersectedLine[1] + edgeDirection[2]*intersectedLine[2];
                        if(fabs(dotProduct) > maxDotProduct && m_topology->getTrianglesAroundEdge(edgesAroundVer[j]).size()>1)
                        {
                            secondEdgeId = edgesAroundVer[j];
                            maxDotProduct = fabs(dotProduct);
                        }
                    }
                }
            }
            else if(pointSide == 0)
            {
                EdgesAroundVertex edgesAroundVer = m_topology->getEdgesAroundVertex(edge[0]);
                for(unsigned int j=0;j<edgesAroundVer.size();j++)
                {
                    if(edgesAroundVer[j] != firstEdgeId)
                    {
                        Coord edgePointsPos[2];
                        triangleGeo->getEdgeVertexCoordinates(edgesAroundVer[j], edgePointsPos);
                        Coord edgeDirection = edgePointsPos[1]-edgePointsPos[0];
                        edgeDirection /= sqrt(edgeDirection[0]*edgeDirection[0]+edgeDirection[1]*edgeDirection[1]+edgeDirection[2]*edgeDirection[2]);
                        dotProduct = edgeDirection[0]*intersectedLine[0] + edgeDirection[1]*intersectedLine[1] + edgeDirection[2]*intersectedLine[2];
                        if(fabs(dotProduct) > maxDotProduct && m_topology->getTrianglesAroundEdge(edgesAroundVer[j]).size()>1)
                        {
                            secondEdgeId = edgesAroundVer[j];
                            maxDotProduct = fabs(dotProduct);
                        }
                    }
                }
            }


            if(secondEdgeId!=0 && firstEdgeId!=0)
            {
                sofa::helper::vector<EdgeID> inciesEdge;
                inciesEdge.clear();

                inciesEdge.push_back(firstEdgeId);
                inciesEdge.push_back(secondEdgeId);

                sofa::helper::vector<PointID> new_points;
                sofa::helper::vector<PointID> end_points;
                bool reachBorder = false;

                bool incieseOk = triangleAlg->InciseAlongEdgeList(inciesEdge, new_points, end_points, reachBorder);
                dmsg_warning() << "end_points " << end_points;   

                if (!end_points.empty())
                {
                    firstCut = false;
                    if(pointSide == 1)
                    {
                        edgesAroundPa.clear();
                        lastPa = end_points.back();
                        edgesAroundPa = m_topology->getEdgesAroundVertex(lastPa);
                    }
                    else if(pointSide == 2)
                    {
                        edgesAroundPb.clear();
                        lastPb = end_points.back();
                        edgesAroundPb = m_topology->getEdgesAroundVertex(lastPb);
                    }
                    else if(pointSide == 0)
                    {
                        edgesAroundPa.clear();
                        lastPa = end_points.back();
                        edgesAroundPa = m_topology->getEdgesAroundVertex(lastPa);
                        if(!end_points.empty())
                        {
                            edgesAroundPb.clear();
                            lastPb = end_points.back();
                            edgesAroundPb = m_topology->getEdgesAroundVertex(lastPb);
                        }
                    }
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



