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
    if(stepCnt%100==0)
    {
        Real maxStressOfAll = 0;
        TriangleID maxStressTriangleID=0;
        int pointSide=0;

        std::vector<triangleInfo> triangleInfo = *(triangleFF->triangleInfo.beginEdit());

        trueEdgeAroundPa.clear();
        trueEdgeAroundPb.clear();

        //pointVertices.clear();
        //lineVertices.clear();
        //triangleVertices.clear();

        size_t nbTriangle = m_topology->getNbTriangles();
        for (unsigned int i = 0; i < nbTriangle; ++i)
        {
            if (fabs(triangleInfo[i].maxStress) > maxStressOfAll)
            {
                maxStressOfAll = fabs(triangleInfo[i].maxStress);
                maxStressTriangleID = i;
                pointSide = 0;
            }
        }

        if (maxStressOfAll > d_tearThreshold.getValue() && firstCut)
        {
            //firstCut = false;
            DataTypes::Coord stressDirection;
            sofa::type::Vec<3, Real> stressDirectionVec;
            sofa::type::Vec<3, double> triangleNormal;

            triangleNormal = triangleGeo->computeTriangleNormal(maxStressTriangleID);
            stressDirection = triangleInfo[maxStressTriangleID].principalStressDirection;
            stressDirectionVec[0] = (Real)(stressDirection[0]);
            stressDirectionVec[1] = (Real)(stressDirection[1]);
            stressDirectionVec[2] = (Real)(stressDirection[2]);
            sofa::type::Vec<3, Real> intersectedLine = stressDirectionVec.cross(triangleNormal);
            intersectedLine = intersectedLine.normalized();

            Coord p[3];
            triangleGeo->getTriangleVertexCoordinates(maxStressTriangleID, p);
            /*triangleVertices.push_back(sofa::type::Vector3(p[0]));
            triangleVertices.push_back(sofa::type::Vector3(p[1]));
            triangleVertices.push_back(sofa::type::Vector3(p[2]));*/

            unsigned int ver1, ver2, ver3, t2v3;
            sofa::type::Vec<3, Real> edgeDirectionVec;
            double maxDotProduct = 0, dotProduct;
            for (unsigned int i = 0; i < 3; i++)
            {
                for (unsigned int j = i + 1; j < 3; j++)
                {
                    Coord edgeDirection = p[j] - p[i];
                    edgeDirectionVec[0] = (Real)(edgeDirection[0]);
                    edgeDirectionVec[1] = (Real)(edgeDirection[1]);
                    edgeDirectionVec[2] = (Real)(edgeDirection[2]);
                    edgeDirectionVec = edgeDirectionVec.normalized();
                    dotProduct = edgeDirectionVec * intersectedLine;

                    if (fabs(dotProduct) > maxDotProduct)
                    {
                        ver1 = i; ver2 = j; ver3 = (unsigned int)fabs((int)i + (int)j - (int)3);
                        maxDotProduct = fabs(dotProduct);
                    }
                }
            }

            const Triangle& t1 = m_topology->getTriangle(maxStressTriangleID);
            EdgeID edgeID = m_topology->getEdgeIndex(t1[ver1], t1[ver2]);
            std::vector<TriangleID> tAroundEdge = m_topology->getTrianglesAroundEdge(edgeID);

            if (tAroundEdge.size() == 2)
            {
                type::vector< type::vector< Index > > p_ancestors(1);
                sofa::type::vector< type::vector< double > > p_baryCoefs(1);
                for (int i = 0; i < p_ancestors.size(); i++)
                {
                    auto& ancestor = p_ancestors[i];
                    ancestor.resize(3);
                    ancestor[0] = t1[0];
                    ancestor[1] = t1[1];
                    ancestor[2] = t1[2];
                }
                for (int i = 0; i < p_baryCoefs.size(); i++)
                {
                    type::vector<double>& baryCoef = p_baryCoefs[i];
                    baryCoef.resize(3);
                    baryCoef.fill(0);
                    baryCoef[ver1] = 0.5;
                    baryCoef[ver2] = 0.5;
                }

                triangleMod->addPoints(1, p_ancestors, p_baryCoefs, true);
                triangleMod->addPoints(1, p_ancestors, p_baryCoefs, true);
                PointID newPointIndex1 = m_topology->getNbPoints() - 1;
                PointID newPointIndex2 = m_topology->getNbPoints() - 2;
                dmsg_warning() << "newPointIndex1: " << newPointIndex1;
                dmsg_warning() << "newPointIndex2: " << newPointIndex2;
                const Coord& NewPos1 = triangleGeo->getPointPosition(newPointIndex1);
                const Coord& NewPos2 = triangleGeo->getPointPosition(newPointIndex2);
                dmsg_warning() << "NewPos1: " << NewPos1;
                dmsg_warning() << "NewPos2: " << NewPos2;

                //Draw new point
                const Coord& pPos1 = triangleGeo->getPointPosition(newPointIndex1);
                pointVertices.push_back(sofa::type::Vector3(pPos1)); 
                const Coord& pPos2 = triangleGeo->getPointPosition(newPointIndex2);
                pointVertices.push_back(sofa::type::Vector3(pPos2));
                //Draw max stress triangle
                /*const Coord& tPos1 = triangleGeo->getPointPosition(t1[0]);
                triangleVertices.push_back(tPos1);
                const Coord& tPos2 = triangleGeo->getPointPosition(t1[1]);
                triangleVertices.push_back(tPos2);
                const Coord& tPos3 = triangleGeo->getPointPosition(t1[2]);
                triangleVertices.push_back(tPos3);*/

                const Triangle& t2 = m_topology->getTriangle((tAroundEdge[0] == maxStressTriangleID) ? tAroundEdge[1] : tAroundEdge[0]);
                if ((t2[0] == t1[ver1] && t2[1] == t1[ver2]) || (t2[0] == t1[ver2] && t2[1] == t1[ver1]))       t2v3 = 2;
                else if ((t2[1] == t1[ver1] && t2[2] == t1[ver2]) || (t2[1] == t1[ver2] && t2[2] == t1[ver1]))  t2v3 = 0;
                else if ((t2[0] == t1[ver1] && t2[2] == t1[ver2]) || (t2[0] == t1[ver2] && t2[2] == t1[ver1]))  t2v3 = 1;

                
                type::vector<Triangle> newTriangle;
                newTriangle.resize(4);
                newTriangle[0][0] = t1[ver3]; newTriangle[0][1] = t1[ver2]; newTriangle[0][2] = newPointIndex1;
                newTriangle[1][1] = t1[ver3]; newTriangle[1][0] = t1[ver1]; newTriangle[1][2] = newPointIndex1;
                newTriangle[2][1] = t2[t2v3]; newTriangle[2][0] = t1[ver2]; newTriangle[2][2] = newPointIndex2;
                newTriangle[3][0] = t2[t2v3]; newTriangle[3][1] = t1[ver1]; newTriangle[3][2] = newPointIndex2;
                dmsg_warning() << "newTriangle: " << newTriangle;

                triangleMod->addTriangles(newTriangle);

                for (int i = 0; i < 4; i++)
                {
                    const Triangle& newT1 = m_topology->getTriangle(nbTriangle+i);
                    const Coord& newTPos1 = triangleGeo->getPointPosition(newT1[0]);
                    triangleVertices.push_back(newTPos1);
                    const Coord& newTPos2 = triangleGeo->getPointPosition(newT1[1]);
                    triangleVertices.push_back(newTPos2);
                    const Coord& newTPos3 = triangleGeo->getPointPosition(newT1[2]);
                    triangleVertices.push_back(newTPos3);
                }

               /* dmsg_warning() << "Start cut: ";
                type::vector<EdgeID> inciesEdge;
                inciesEdge.resize(2);
                inciesEdge[0] = m_topology->getEdgeIndex(t1[ver1], newPointIndex);
                inciesEdge[1] = m_topology->getEdgeIndex(t1[ver2], newPointIndex);
                dmsg_warning() << "insices edge " << inciesEdge;

                sofa::type::vector<PointID> new_points;
                sofa::type::vector<PointID> end_points;
                bool reachBorder = false;

                bool incieseOk = triangleGeo->InciseAlongEdgeList(inciesEdge, new_points, end_points, reachBorder);
                dmsg_warning() << "end_points " << end_points;*/

                triangleMod->removeItems(tAroundEdge);
            }
        }
        else
        {
            maxStressTriangleID = 0;
        }
        

        /*for(unsigned int i=0;i<edgesAroundPa.size();i++)
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
        }*/
    }
}

// draw for standard types (i.e Vec<1,2,3>)
template<class DataTypes>
void TearPlugin<DataTypes>::draw(const core::visual::VisualParams* vparams)
{
     if(d_draw.getValue())
     {
         sofa::type::RGBAColor color = sofa::type::RGBAColor(1, 0, 0, 1);
         vparams->drawTool()->drawPoints(pointVertices, 5, color);
         color = sofa::type::RGBAColor(1,0,0,1);
         vparams->drawTool()->drawLines(lineVertices, 2, color);
         color = sofa::type::RGBAColor(0,0,1,1);
         vparams->drawTool()->drawTriangles(triangleVertices, color);
         //lineVertices.clear();
         //triangleVertices.clear();
         shouldDraw = false;
     }
}


} // namespace sofa::component::forcefield

#endif // SOFA_TEARPLUGIN_INL



