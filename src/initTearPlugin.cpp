/******************************************************************************
*       SOFA, Simulation Open-Framework Architecture, development version     *
*                (c) 2006-2017 INRIA, USTL, UJF, CNRS, MGH                    *
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
#include <config.h>

namespace sofa
{
namespace component
{

extern "C" {

SOFATEARING_API
void initExternalModule()
{
    static bool first = true;
    if (first)
    {
        first = false;
    }
}

SOFATEARING_API
const char* getModuleName()
{
    return "SofaTearing";
}

SOFATEARING_API
const char* getModuleVersion()
{
    return "0.1";
}

SOFATEARING_API
const char* getModuleLicense()
{
    return "LGPL";
}

SOFATEARING_API
const char* getModuleDescription()
{
    return "SOFA tearing plugin";
}

SOFATEARING_API
const char* getModuleComponentList()
{
    // string containing the names of the classes provided by the plugin
    return "Tear";
}

} // extern "C"

} // namespace component
} // namespace sofa
