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
    if(stepCnt%10==0 && !reachEnd)
    {
        /*if (reachEnd)
        {
            firstCut = true;
            trueEdgeAroundPa.clear();
            trueEdgeAroundPb.clear();
        }*/

        Real maxStressOfAll = 0;
        TriangleID maxStressTriangleID=0;
        int pointSide=0;

        std::vector<triangleInfo> triangleInfo = *(triangleFF->triangleInfo.beginEdit());

        //pointVertices.clear();
        //lineVertices.clear();
        //triangleVertices.clear();
        for (unsigned int i = 0; i < trueEdgeAroundPa.size(); i++)
        {
            std::vector<TriangleID> tri = m_topology->getTrianglesAroundEdge(trueEdgeAroundPa[i]);
            if (tri.size() > 1)
            {
                for (unsigned int j = 0; j < tri.size(); j++)
                {
                    if (fabs(triangleInfo[tri[j]].maxStress) > maxStressOfAll)
                    {
                        maxStressOfAll = fabs(triangleInfo[tri[j]].maxStress);
                        maxStressTriangleID = tri[j];
                        pointSide = 1;
                    }
                }
            }
        }

        for (unsigned int i = 0; i < trueEdgeAroundPb.size(); i++)
        {
            std::vector<TriangleID> tri = m_topology->getTrianglesAroundEdge(trueEdgeAroundPb[i]);
            if (tri.size() > 1)
            {
                for (unsigned int j = 0; j < tri.size(); j++)
                {
                    if (fabs(triangleInfo[tri[j]].maxStress) > maxStressOfAll)
                    {
                        maxStressOfAll = fabs(triangleInfo[tri[j]].maxStress);
                        maxStressTriangleID = tri[j];
                        pointSide = 2;
                    }
                }
            }
        }

        if (firstCut)
        {
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
        }

        if (maxStressOfAll > d_tearThreshold.getValue())
        {
            firstCut = false;
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

            sofa::type::Vec<3, Real> edgeDirectionVec;
            double maxDotProduct = 0, dotProduct;
            EdgeID inciesEdge;

            if (pointSide == 0)
            {
                EdgesInTriangle eInTri = m_topology->getEdgesInTriangle(maxStressTriangleID);
                for (unsigned int i = 0; i < 3; i++)
                {
                    Coord edgesPos[2];
                    triangleGeo->getEdgeVertexCoordinates(eInTri[i], edgesPos);

                    Coord edgeDirection = p[1] - p[0];
                    edgeDirectionVec[0] = (Real)(edgeDirection[0]);
                    edgeDirectionVec[1] = (Real)(edgeDirection[1]);
                    edgeDirectionVec[2] = (Real)(edgeDirection[2]);
                    edgeDirectionVec = edgeDirectionVec.normalized();
                    dotProduct = edgeDirectionVec * intersectedLine;

                    if (fabs(dotProduct) > maxDotProduct)
                    {
                        inciesEdge = eInTri[i];
                        maxDotProduct = fabs(dotProduct);
                    }
                }
            }
            else
            {
                EdgesAroundVertex edges;
                if (pointSide == 1)         edges = trueEdgeAroundPa;
                else if (pointSide == 2)    edges = trueEdgeAroundPb;

                for (int i = 0; i < edges.size(); i++)
                {
                    Coord edgesPos[2];
                    triangleGeo->getEdgeVertexCoordinates(edges[i], edgesPos);

                    Coord edgeDirection = edgesPos[1] - edgesPos[0];
                    edgeDirectionVec[0] = (Real)(edgeDirection[0]);
                    edgeDirectionVec[1] = (Real)(edgeDirection[1]);
                    edgeDirectionVec[2] = (Real)(edgeDirection[2]);
                    edgeDirectionVec = edgeDirectionVec.normalized();
                    dotProduct = edgeDirectionVec * intersectedLine;

                    if (fabs(dotProduct) > maxDotProduct)
                    {
                        inciesEdge = edges[i];
                        maxDotProduct = fabs(dotProduct);
                    }
                }
            }
            Edge pointInInciesEdge = m_topology->getEdge(inciesEdge);
            std::vector<TriangleID> tAroundInciesEdge = m_topology->getTrianglesAroundEdge(inciesEdge);
            const Triangle& t1 = m_topology->getTriangle(tAroundInciesEdge[0]);
            const Triangle& t2 = m_topology->getTriangle(tAroundInciesEdge[1]);
            PointID t1v1, t1v2, t1v3, t2v3;
            for (int i = 0; i < 3; i++)
            {
                if (t1[i] == pointInInciesEdge[0]) t1v1 = i;
                if (t1[i] == pointInInciesEdge[1]) t1v2 = i;
                if (t1[i] != pointInInciesEdge[0] && t1[i] != pointInInciesEdge[1])
                    t1v3 = i;
            }
            for (int i = 0; i < 3; i++)
            {
                if (t2[i] != pointInInciesEdge[0] && t2[i] != pointInInciesEdge[1])
                    t2v3 = i;
            }

            type::vector< type::vector< Index > > p_ancestors(1);
            sofa::type::vector< type::vector< double > > p_baryCoefs(1);
            auto& ancestor = p_ancestors[0];
            ancestor.resize(3);
            ancestor[0] = t1[0];
            ancestor[1] = t1[1];
            ancestor[2] = t1[2];

            type::vector<double>& baryCoef = p_baryCoefs[0];
            baryCoef.resize(3);
            baryCoef.fill(0);
            baryCoef[t1v1] = 0.5;
            baryCoef[t1v2] = 0.5;

            triangleMod->addPoints(1, p_ancestors, p_baryCoefs, true);
            PointID newPointIndex = m_topology->getNbPoints() - 1;
            dmsg_warning() << "newPointIndex1: " << newPointIndex;
            sofa::type::Vec<3, Real> NewPos = triangleGeo->getPointPosition(newPointIndex);
            dmsg_warning() << "NewPos1: " << NewPos;
            //Draw new point
            //pointVertices.push_back(sofa::type::Vector3(NewPos));

            //Draw max stress triangle
            /*const Coord& tPos1 = triangleGeo->getPointPosition(t1[0]);
            triangleVertices.push_back(tPos1);
            const Coord& tPos2 = triangleGeo->getPointPosition(t1[1]);
            triangleVertices.push_back(tPos2);
            const Coord& tPos3 = triangleGeo->getPointPosition(t1[2]);
            triangleVertices.push_back(tPos3);*/

              
            type::vector<Triangle> newTriangle;
            newTriangle.resize(4);
            newTriangle[0][0] = t1[t1v3]; newTriangle[0][1] = pointInInciesEdge[0]; newTriangle[0][2] = newPointIndex;
            newTriangle[1][1] = t1[t1v3]; newTriangle[1][0] = pointInInciesEdge[1]; newTriangle[1][2] = newPointIndex;
            newTriangle[2][1] = t2[t2v3]; newTriangle[2][0] = pointInInciesEdge[0]; newTriangle[2][2] = newPointIndex;
            newTriangle[3][0] = t2[t2v3]; newTriangle[3][1] = pointInInciesEdge[1]; newTriangle[3][2] = newPointIndex;
            dmsg_warning() << "newTriangle: " << newTriangle;

            triangleMod->addTriangles(newTriangle);
            triangleMod->removeItems(tAroundInciesEdge);

            dmsg_warning() << "Start cut: ";
            type::vector<EdgeID> inciesEdges;
            inciesEdges.resize(2);
            inciesEdges[0] = m_topology->getEdgeIndex(pointInInciesEdge[0], newPointIndex);
            inciesEdges[1] = m_topology->getEdgeIndex(pointInInciesEdge[1], newPointIndex);
            dmsg_warning() << "insices edge " << inciesEdges;
            /*const Coord& newTPos1 = triangleGeo->getPointPosition(pointInInciesEdge[0]);
            lineVertices.push_back(newTPos1);
            const Coord& newTPos2 = triangleGeo->getPointPosition(newPointIndex);
            lineVertices.push_back(newTPos2);
            const Coord& newTPos3 = triangleGeo->getPointPosition(pointInInciesEdge[1]);
            lineVertices.push_back(newTPos3);*/

            sofa::type::vector<PointID> new_points;
            sofa::type::vector<PointID> end_points;
            bool reachBorder = false;

            bool incieseOk = triangleGeo->InciseAlongEdgeList(inciesEdges, new_points, end_points, reachBorder);
            dmsg_warning() << "end_points " << end_points << " new_points " << new_points;
            if (end_points.size() < 1)
            {
                reachEnd = true;
                dmsg_warning() << ">>>>>>>>>>>>>>>>>>>>>>>>>> REACH END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<";
            }
            if (incieseOk && !reachEnd)
            {
                if (pointSide == 0)
                {
                    lastPa = end_points[0];
                    trueEdgeAroundPa.clear();

                    type::vector<EdgeID> forbiddenEdge;
                    forbiddenEdge.push_back(m_topology->getEdgeIndex(t1[t1v3], lastPa));
                    forbiddenEdge.push_back(m_topology->getEdgeIndex(t2[t2v3], lastPa));
                    forbiddenEdge.push_back(m_topology->getEdgeIndex(new_points[0], lastPa));
                    forbiddenEdge.push_back(m_topology->getEdgeIndex(newPointIndex, lastPa));

                    edgesAroundPa = m_topology->getEdgesAroundVertex(lastPa);
                    for (int i = 0; i < edgesAroundPa.size(); i++)
                    {
                        bool exist = false;
                        for (int j = 0; j < forbiddenEdge.size(); j++)
                        {
                            if (edgesAroundPa[i] == forbiddenEdge[j])
                                exist = true;
                        }
                        if (!exist)
                            trueEdgeAroundPa.push_back(edgesAroundPa[i]);
                    }
                    dmsg_warning() << "trueEdgeAroundPa " << trueEdgeAroundPa;


                    lastPb = end_points[1];
                    trueEdgeAroundPb.clear();

                    forbiddenEdge.clear();
                    forbiddenEdge.push_back(m_topology->getEdgeIndex(t1[t1v3], lastPb));
                    forbiddenEdge.push_back(m_topology->getEdgeIndex(t2[t2v3], lastPb));
                    forbiddenEdge.push_back(m_topology->getEdgeIndex(new_points[0], lastPb));
                    forbiddenEdge.push_back(m_topology->getEdgeIndex(newPointIndex, lastPb));

                    edgesAroundPb = m_topology->getEdgesAroundVertex(lastPb);
                    for (int i = 0; i < edgesAroundPb.size(); i++)
                    {
                        bool exist = false;
                        for (int j = 0; j < forbiddenEdge.size(); j++)
                        {
                            if (edgesAroundPb[i] == forbiddenEdge[j])
                                exist = true;
                        }
                        if (!exist)
                            trueEdgeAroundPb.push_back(edgesAroundPb[i]);
                    }
                    dmsg_warning() << "trueEdgeAroundPa " << trueEdgeAroundPb;
                }
                else if (pointSide == 1)
                {
                    lastPa = end_points[0];
                    trueEdgeAroundPa.clear();

                    type::vector<EdgeID> forbiddenEdge;
                    forbiddenEdge.push_back(m_topology->getEdgeIndex(t1[t1v3], lastPa));
                    forbiddenEdge.push_back(m_topology->getEdgeIndex(t2[t2v3], lastPa));
                    forbiddenEdge.push_back(m_topology->getEdgeIndex(newPointIndex, lastPa));

                    sofa::type::Vec<3, Real> checkPos = triangleGeo->getPointPosition(new_points[0]);
                    unsigned int newPointIdx;
                    if ((checkPos - NewPos).norm() == 0) { newPointIdx = 0; dmsg_warning() << "is Zero"; }
                    else                                 newPointIdx = 1;
                    forbiddenEdge.push_back(m_topology->getEdgeIndex(new_points[newPointIdx], lastPa));

                    edgesAroundPa = m_topology->getEdgesAroundVertex(lastPa);
                    for (int i = 0; i < edgesAroundPa.size(); i++)
                    {
                        bool exist = false;
                        for (int j = 0; j < forbiddenEdge.size(); j++)
                        {
                            if (edgesAroundPa[i] == forbiddenEdge[j])
                                exist = true;
                        }
                        if (!exist)
                            trueEdgeAroundPa.push_back(edgesAroundPa[i]);
                    }
                    dmsg_warning() << "trueEdgeAroundPa " << trueEdgeAroundPa;
                }
                else if (pointSide == 2)
                {
                    lastPb = end_points[0];
                    trueEdgeAroundPb.clear();

                    type::vector<EdgeID> forbiddenEdge;
                    forbiddenEdge.push_back(m_topology->getEdgeIndex(t1[t1v3], lastPb));
                    forbiddenEdge.push_back(m_topology->getEdgeIndex(t2[t2v3], lastPb));
                    forbiddenEdge.push_back(m_topology->getEdgeIndex(newPointIndex, lastPb));

                    sofa::type::Vec<3, Real> checkPos = triangleGeo->getPointPosition(new_points[0]);
                    unsigned int newPointIdx;
                    if ((checkPos - NewPos).norm() == 0) { newPointIdx = 0; dmsg_warning() << "is Zero"; }
                    else                                    newPointIdx = 1;
                    forbiddenEdge.push_back(m_topology->getEdgeIndex(new_points[newPointIdx], lastPb));

                    edgesAroundPb = m_topology->getEdgesAroundVertex(lastPb);
                    for (int i = 0; i < edgesAroundPb.size(); i++)
                    {
                        bool exist = false;
                        for (int j = 0; j < forbiddenEdge.size(); j++)
                        {
                            if (edgesAroundPb[i] == forbiddenEdge[j])
                                exist = true;
                        }
                        if (!exist)
                            trueEdgeAroundPb.push_back(edgesAroundPb[i]);
                    }
                    dmsg_warning() << "trueEdgeAroundPa " << trueEdgeAroundPb;
                }
                
                    
                /*TrianglesAroundVertex tAround = m_topology->getTrianglesAroundEdge(forbiddenEdge[0]);
                forbiddenTriangle.push_back(tAround[0]); forbiddenTriangle.push_back(tAround[1]);
                tAround = m_topology->getTrianglesAroundEdge(forbiddenEdge[1]);
                forbiddenTriangle.push_back(tAround[0]); forbiddenTriangle.push_back(tAround[1]);

                for (int i = 0; i < AllTriAroundP.size(); i++)
                {
                    bool exist = false;
                    for (int j = 0; j < forbiddenTriangle.size(); j++)
                    {
                        if (AllTriAroundP[i] == forbiddenTriangle[j])
                            exist = true;
                    }
                    if (!exist)
                    {
                        trianglesAroundPointA.push_back(AllTriAroundP[i]);
                    }
                }*/
            }

            /*std::vector<TriangleID> tAroundFirstP = m_topology->getTrianglesAroundVertex(newPointIndex);
            std::vector<TriangleID> tAroundSecondP = m_topology->getTrianglesAroundVertex(new_points[0]);
                
            for (int i = 0; i < tAroundFirstP.size(); i++)
            {
                const Triangle& newT1 = m_topology->getTriangle(tAroundFirstP[i]);
                const Coord& newTPos1 = triangleGeo->getPointPosition(newT1[0]);
                triangleVertices.push_back(newTPos1);
                const Coord& newTPos2 = triangleGeo->getPointPosition(newT1[1]);
                triangleVertices.push_back(newTPos2);
                const Coord& newTPos3 = triangleGeo->getPointPosition(newT1[2]);
                triangleVertices.push_back(newTPos3);
            }

            for (int i = 0; i < tAroundSecondP.size(); i++)
            {
                const Triangle& newT1 = m_topology->getTriangle(tAroundSecondP[i]);
                const Coord& newTPos1 = triangleGeo->getPointPosition(newT1[0]);
                triangleVertices.push_back(newTPos1);
                const Coord& newTPos2 = triangleGeo->getPointPosition(newT1[1]);
                triangleVertices.push_back(newTPos2);
                const Coord& newTPos3 = triangleGeo->getPointPosition(newT1[2]);
                triangleVertices.push_back(newTPos3);
            }*/

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
         color = sofa::type::RGBAColor(0,1,0,1);
         vparams->drawTool()->drawLines(lineVertices, 2, color);

         size_t nbTri = triangleVertices.size() / 3;
         for (int i = 0; i < nbTri; i++)
         {
             std::vector<type::Vector3> drawTri;
             drawTri.push_back(triangleVertices[3 * i + 0]);
             drawTri.push_back(triangleVertices[3 * i + 1]);
             drawTri.push_back(triangleVertices[3 * i + 2]);

             color = sofa::type::RGBAColor(i / (float)nbTri, 0, 1 - (i / (float)nbTri), 1);
             vparams->drawTool()->drawTriangles(drawTri, color);
         }
         //lineVertices.clear();
         //triangleVertices.clear();
         shouldDraw = false;
     }
}


} // namespace sofa::component::forcefield

#endif // SOFA_TEARPLUGIN_INL



