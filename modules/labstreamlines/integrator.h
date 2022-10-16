/*********************************************************************
 *  Author  : Himangshu Saikia
 *  Init    : Wednesday, September 20, 2017 - 12:04:15
 *
 *  Project : KTH Inviwo Modules
 *
 *  License : Follows the Inviwo BSD license model
 *********************************************************************
 */

#pragma once

#include <inviwo/core/common/inviwo.h>
#include <inviwo/core/datastructures/geometry/basicmesh.h>
#include <inviwo/core/datastructures/volume/volumeram.h>
#include <labstreamlines/labstreamlinesmoduledefine.h>
#include <labutils/scalarvectorfield.h>

namespace inviwo {

class IVW_MODULE_LABSTREAMLINES_API Integrator {

public:
    // Construction / Deconstruction
public:
    Integrator() {}
    virtual ~Integrator() = default;

    // Methods
public:

    // Add a point to a mesh
    static void drawPoint(const dvec2& p, const vec4& color, IndexBufferRAM* indexBuffer,
                          std::vector<BasicMesh::Vertex>& vertices);
    // Add a point to a polyline, assumes that the indexBuffer uses Strip Connectivity
    static void drawNextPointInPolyline(const dvec2& p, const vec4& color,
                                        IndexBufferRAM* indexBuffer,
                                        std::vector<BasicMesh::Vertex>& vertices);
    // Add a line segment to a mesh, assuming no connectivity in the indexBuffer
    static void drawLineSegment(const dvec2& v1, const dvec2& v2, const vec4& color,
                                IndexBufferRAM* indexBuffer,
                                std::vector<BasicMesh::Vertex>& vertices);


    // TODO: Implement the methods below (one integration step with either Euler or
    // Runge-Kutte of 4th order integration method)
    // Pass any other properties that influence the integration process
    // Examples would be the stepsize, inegreation direction, ...
    // Integration step with  Euler
    static dvec2 Euler(const VectorField2& vectorField, const dvec2& position, float stepSize, int dir, bool norm);
    static dvec2 RK4(const VectorField2& vectorField, const dvec2& position, float stepSize, int dir, bool norm);
    static std::tuple<dvec2, bool> RK4_(const VectorField2& vectorField, const dvec2& position, float stepSize, int dir,
                                            bool norm, bool zeros, bool stopSlow, float velocity);
    static int StreamLines(bool Rk4type, const VectorField2& vectorField, const dvec2& startPoint,
        bool direction, float stepSize, bool normalized, bool stopSteps, int numIter,
        bool stopArc, float arc, bool boundary, bool zeros, bool stopSlow, float velocity,
        bool displayPoints, IndexBufferRAM* indexBufferPoints, IndexBufferRAM* indexBuffer, std::vector<BasicMesh::Vertex>& vertices);
    
    static std::list<vec2> StreamLine(const VectorField2& vectorField, const  dvec2& startPoint, float stepSize, bool stopArc, float arc);
    static int Separatrix(const VectorField2& vectorField, const  dvec2& startPoint, const  dvec2& criticalPoint, int dir, float stepSize,
                            IndexBufferRAM* indexBufferlines, std::vector<BasicMesh::Vertex>& vertices);
    static int SeparatrixBoundary(const VectorField2& vectorField, const  dvec2& startPoint, float stepSize,
                            IndexBufferRAM* indexBufferlines, std::vector<BasicMesh::Vertex>& vertices);
};

}  // namespace inviwo
