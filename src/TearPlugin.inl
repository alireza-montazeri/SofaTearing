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
#ifndef SOFA_TEARPLUGIN_INL
#define SOFA_TEARPLUGIN_INL

#include <TearPlugin.h>

namespace sofa::component::forcefield
{

template<class DataTypes>
TearPlugin<DataTypes>::TearPlugin()
    : d_tearThreshold(initData(&d_tearThreshold, "tearThreshold", ""))
    , d_draw(initData(&d_draw, false, "draw", "draw"))
{
    // Nothing more is done here
}


template<class DataTypes>
void TearPlugin<DataTypes>::init()
{
    m_topology = this->getContext()->getMeshTopology(); // get the mesh topology to access to the points
    m_topology->getContext()->get(triangleFF);
    m_topology->getContext()->get(triangleMod);
    // m_topology->getContext()->get(triangleAlg);
    m_topology->getContext()->get(triangleGeo);

    firstCut = true;
    stepCnt=0;
    tearStep=0.00001;
    
    shouldDraw = false;
    Inherit::init(); // call parent class init
}


template<class DataTypes>
void TearPlugin<DataTypes>::addForce(const core::MechanicalParams* /*params*/, DataVecDeriv& currentForce, const DataVecCoord& /*currentPosition*/, const DataVecDeriv& /*currentVelocities*/)
{
    stepCnt++;
    // dmsg_warning() << "debug: " << m_topology->getTrianglesAroundVertex(100);
    if(stepCnt%10==0)
    {
        bool isNeighber = true;
        Real maxStressOfAll = 0;
        TriangleID maxStressTriangleID=0;
        int pointSide=0;

        std::vector<triangleInfo> triangleInfo = *(triangleFF->triangleInfo.beginEdit());

        trueEdgeAroundPa.clear();
        trueEdgeAroundPb.clear();

        lineVertices.clear();
        triangleVertices.clear();

        for(unsigned int i=0;i<edgesAroundPa.size();i++)
        {
            std::vector<TriangleID> tri = m_topology->getTrianglesAroundEdge(edgesAroundPa[i]);
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
            std::vector<TriangleID> tri = m_topology->getTrianglesAroundEdge(edgesAroundPb[i]);
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
            sofa::type::Vec<3,Real> stressDirectionVec;
            sofa::type::Vec<3,double> triangleNormal;

            triangleNormal = triangleGeo->computeTriangleNormal(maxStressTriangleID);
            stressDirection  = triangleInfo[maxStressTriangleID].principalStressDirection;
            stressDirectionVec[0] = (Real)(stressDirection[0]); 
            stressDirectionVec[1] = (Real)(stressDirection[1]); 
            stressDirectionVec[2] = (Real)(stressDirection[2]);            
            sofa::type::Vec<3,Real> intersectedLine = stressDirectionVec.cross(triangleNormal);
            intersectedLine = intersectedLine.normalized();

            Coord p[3];
            triangleGeo->getTriangleVertexCoordinates(maxStressTriangleID, p);
            triangleVertices.push_back(sofa::type::Vector3(p[0]));
            triangleVertices.push_back(sofa::type::Vector3(p[1]));
            triangleVertices.push_back(sofa::type::Vector3(p[2]));

            double maxDotProduct=0,dotProduct;
            EdgeID firstEdgeId=0, secondEdgeId=0;
            sofa::type::Vec<3, Real> edgeDirectionVec;

            if(pointSide == 1)
            {
                for(unsigned int i=0;i<trueEdgeAroundPa.size();i++)
                {
                    if (m_topology->getTrianglesAroundEdge(trueEdgeAroundPa[i]).size() > 1)
                    {
                        Coord edgePointsPos[2];
                        triangleGeo->getEdgeVertexCoordinates(trueEdgeAroundPa[i], edgePointsPos);
                        Coord edgeDirection = edgePointsPos[1] - edgePointsPos[0];

                        edgeDirectionVec[0] = (Real)(edgeDirection[0]);
                        edgeDirectionVec[1] = (Real)(edgeDirection[1]);
                        edgeDirectionVec[2] = (Real)(edgeDirection[2]);
                        edgeDirectionVec = edgeDirectionVec.normalized();
                        dotProduct = edgeDirectionVec * intersectedLine;

                        if (fabs(dotProduct) > maxDotProduct)
                        {
                            firstEdgeId = trueEdgeAroundPa[i];
                            maxDotProduct = fabs(dotProduct);
                        }
                    }
                }
            }
            else if(pointSide == 2)
            {
                for(unsigned int i=0;i<trueEdgeAroundPb.size();i++)
                {
                    if (m_topology->getTrianglesAroundEdge(trueEdgeAroundPb[i]).size() > 1)
                    {
                        Coord edgePointsPos[2];
                        triangleGeo->getEdgeVertexCoordinates(trueEdgeAroundPb[i], edgePointsPos);
                        Coord edgeDirection = edgePointsPos[1] - edgePointsPos[0];

                        edgeDirectionVec[0] = (Real)(edgeDirection[0]);
                        edgeDirectionVec[1] = (Real)(edgeDirection[1]);
                        edgeDirectionVec[2] = (Real)(edgeDirection[2]);
                        edgeDirectionVec = edgeDirectionVec.normalized();
                        dotProduct = edgeDirectionVec * intersectedLine;

                        if (fabs(dotProduct) > maxDotProduct)
                        {
                            firstEdgeId = trueEdgeAroundPb[i];
                            maxDotProduct = fabs(dotProduct);
                        }
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
                        edgeDirectionVec[0] = (Real)(edgeDirection[0]);
                        edgeDirectionVec[1] = (Real)(edgeDirection[1]);
                        edgeDirectionVec[2] = (Real)(edgeDirection[2]);
                        edgeDirectionVec = edgeDirectionVec.normalized();
                        dotProduct = edgeDirectionVec * intersectedLine;

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
                    if(edgesAroundVer[j] != firstEdgeId && m_topology->getTrianglesAroundEdge(edgesAroundVer[j]).size() > 1)
                    {
                        Coord edgePointsPos[2];
                        triangleGeo->getEdgeVertexCoordinates(edgesAroundVer[j], edgePointsPos);
                        Coord edgeDirection = edgePointsPos[1]-edgePointsPos[0];
                        edgeDirectionVec[0] = (Real)(edgeDirection[0]);
                        edgeDirectionVec[1] = (Real)(edgeDirection[1]);
                        edgeDirectionVec[2] = (Real)(edgeDirection[2]);
                        edgeDirectionVec = edgeDirectionVec.normalized();
                        dotProduct = edgeDirectionVec * intersectedLine;

                        if(fabs(dotProduct) > maxDotProduct)
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
                    if(edgesAroundVer[j] != firstEdgeId && m_topology->getTrianglesAroundEdge(edgesAroundVer[j]).size() > 1)
                    {
                        Coord edgePointsPos[2];
                        triangleGeo->getEdgeVertexCoordinates(edgesAroundVer[j], edgePointsPos);
                        Coord edgeDirection = edgePointsPos[1]-edgePointsPos[0];
                        edgeDirectionVec[0] = (Real)(edgeDirection[0]);
                        edgeDirectionVec[1] = (Real)(edgeDirection[1]);
                        edgeDirectionVec[2] = (Real)(edgeDirection[2]);
                        edgeDirectionVec = edgeDirectionVec.normalized();
                        dotProduct = edgeDirectionVec * intersectedLine; 
                        
                        if(fabs(dotProduct) > maxDotProduct)
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
                    if(edgesAroundVer[j] != firstEdgeId && m_topology->getTrianglesAroundEdge(edgesAroundVer[j]).size() > 1)
                    {
                        Coord edgePointsPos[2];
                        triangleGeo->getEdgeVertexCoordinates(edgesAroundVer[j], edgePointsPos);
                        Coord edgeDirection = edgePointsPos[1]-edgePointsPos[0];
                        edgeDirectionVec[0] = (Real)(edgeDirection[0]);
                        edgeDirectionVec[1] = (Real)(edgeDirection[1]);
                        edgeDirectionVec[2] = (Real)(edgeDirection[2]);
                        edgeDirectionVec = edgeDirectionVec.normalized();
                        dotProduct = edgeDirectionVec * intersectedLine;
                        
                        if(fabs(dotProduct) > maxDotProduct)
                        {
                            secondEdgeId = edgesAroundVer[j];
                            maxDotProduct = fabs(dotProduct);
                        }
                    }
                }
            }


            if(secondEdgeId!=0 && firstEdgeId!=0)
            {
                std::vector<EdgeID> inciesEdge;
                inciesEdge.clear();

                inciesEdge.push_back(firstEdgeId);
                inciesEdge.push_back(secondEdgeId);

                sofa::type::vector<PointID> new_points;
                sofa::type::vector<PointID> end_points;
                bool reachBorder = false;

                bool incieseOk = triangleGeo->InciseAlongEdgeList(inciesEdge, new_points, end_points, reachBorder);
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

// draw for standard types (i.e Vec<1,2,3>)
template<class DataTypes>
void TearPlugin<DataTypes>::draw(const core::visual::VisualParams* vparams)
{
     if(d_draw.getValue())
     {
         sofa::type::RGBAColor color = sofa::type::RGBAColor(1,0,0,1);
         vparams->drawTool()->drawLines(lineVertices, 5, color);
         color = sofa::type::RGBAColor(0,0,1,1);
         vparams->drawTool()->drawTriangles(triangleVertices, color);
         //lineVertices.clear();
         //triangleVertices.clear();
         shouldDraw = false;
     }
}


} // namespace sofa::component::forcefield

#endif // SOFA_TEARPLUGIN_INL



