/*********************************************************************
 *  Author  : Himangshu Saikia
 *  Init    : Wednesday, September 20, 2017 - 12:04:15
 *
 *  Project : KTH Inviwo Modules
 *
 *  License : Follows the Inviwo BSD license model
 *********************************************************************
 */

#include <labstreamlines/integrator.h>


namespace inviwo {

// TODO: Implement a single integration step here

dvec2 Integrator::Euler(const VectorField2& vectorField, const dvec2& position, float stepSize, int dir, bool norm)
 {
     dvec2 interpolationValue = dir * vectorField.interpolate(position);
     if (norm) {
         interpolationValue = glm::normalize(interpolationValue);
     }
     interpolationValue[0] = stepSize * interpolationValue[0]; 
     interpolationValue[1] = stepSize * interpolationValue[1];
     dvec2 newPosition = position + interpolationValue; 

     return newPosition; 
}

dvec2 Integrator::RK4(const VectorField2& vectorField, const dvec2& position, float stepSize, int dir, bool norm)
{
    dvec2 v1 = dir * vectorField.interpolate(position);
    if (norm && (glm::length(v1) > 1e-4)) {
        v1 = glm::normalize(v1);
    }

    dvec2 s_v1 = v1;
    s_v1[0] = stepSize/2 * v1[0];
    s_v1[1] = stepSize/2 * v1[1];
    /*if (norm) {
        s_v1 = glm::normalize(s_v1);
    }*/
    dvec2 v2 = dir * vectorField.interpolate(position + s_v1);
    if (norm && (glm::length(v2) > 1e-4)) {
        v2 = glm::normalize(v2);
    }

    dvec2 s_v2 = v2;
    s_v2[0] = stepSize/2 * v2[0];
    s_v2[1] = stepSize/2 * v2[1];
    /*if (norm) {
        s_v2 = glm::normalize(s_v2);
    }*/
    dvec2 v3 = dir * vectorField.interpolate(position + s_v2);
    if (norm && (glm::length(v3) > 1e-4)) {
        v3 = glm::normalize(v3);
    }

    dvec2 s_v3 = v3;
    s_v3[0] = stepSize * v3[0];
    s_v3[1] = stepSize * v3[1];
    /*if (norm) {
        s_v3 = glm::normalize(s_v3);
    }*/
    dvec2 v4 = dir * vectorField.interpolate(position + s_v3);
    if (norm && (glm::length(v4) > 1e-4)) {
        v4 = glm::normalize(v4);
    }

    v1 = v1 / 6;
    v2 = v2 / 3;
    v3 = v3 / 3;
    v4 = v4 / 6;
    dvec2 v = (v1 + v2 + v3 + v4); 
    if (norm && (glm::length(v) > 1e-4)) {
        v = glm::normalize(v); 
    }
    
    dvec2 newPosition = position; 
    newPosition[0] = position[0] + stepSize * (v[0]); 
    newPosition[1] = position[1] + stepSize * (v[1]);
    
    return newPosition; 
}


std::tuple<dvec2, bool> Integrator::RK4_(const VectorField2& vectorField, const dvec2& position, float stepSize, int dir, 
                                           bool norm, bool zeros, bool stopSlow, float velocity)
{
    dvec2 newPosition = position;
    bool stop = false;
    double episilon = 1e-4;
    dvec2 v1 = dir * vectorField.interpolate(position);
    /*if (norm && (glm::length(v1) > episilon)) {
        v1 = glm::normalize(v1);
    }*/

    dvec2 s_v1 = v1;
    s_v1[0] = stepSize / 2 * v1[0];
    s_v1[1] = stepSize / 2 * v1[1];
    dvec2 v2 = dir * vectorField.interpolate(position + s_v1);
    //if (norm && (glm::length(v2) > episilon)) {
    //    v2 = glm::normalize(v2);
    //}

    dvec2 s_v2 = v2;
    s_v2[0] = stepSize / 2 * v2[0];
    s_v2[1] = stepSize / 2 * v2[1];
    dvec2 v3 = dir * vectorField.interpolate(position + s_v2);
    //if (norm && (glm::length(v3) > episilon)) {
    //    v3 = glm::normalize(v3);
    //}

    dvec2 s_v3 = v3;
    s_v3[0] = stepSize * v3[0];
    s_v3[1] = stepSize * v3[1];
    dvec2 v4 = dir * vectorField.interpolate(position + s_v3);
    //if (norm && (glm::length(v4) > episilon)) {
    //    v4 = glm::normalize(v4);
    //}

    v1 = v1 / 6;
    v2 = v2 / 3;
    v3 = v3 / 3;
    v4 = v4 / 6;
    dvec2 v = (v1 + v2 + v3 + v4);

    newPosition[0] = position[0] + stepSize * (v[0]);
    newPosition[1] = position[1] + stepSize * (v[1]);

    // Check if final vector magnitude is too small
    if (zeros && (glm::length(v) < episilon)) {
        stop = true;
    }
    
    // Check if velocity is too low
    if(stopSlow && (glm::length(newPosition-position) < velocity)){
        stop = true;
    }

    // Check if final vector is big enough to normalize
    if (norm && (glm::length(v) >= episilon)) {
        v = glm::normalize(v);
    }

    newPosition[0] = position[0] + stepSize * (v[0]);
    newPosition[1] = position[1] + stepSize * (v[1]);

    return std::tuple(newPosition, stop);
}
void Integrator::drawPoint(const dvec2& p, const vec4& color, IndexBufferRAM* indexBuffer,
                           std::vector<BasicMesh::Vertex>& vertices) {
    indexBuffer->add(static_cast<std::uint32_t>(vertices.size()));
    vertices.push_back({vec3(p[0], p[1], 0), vec3(0, 0, 1), vec3(p[0], p[1], 0), color});
}

// Alias for draw point
void Integrator::drawNextPointInPolyline(const dvec2& p, const vec4& color,
                                         IndexBufferRAM* indexBuffer,
                                         std::vector<BasicMesh::Vertex>& vertices) {
    Integrator::drawPoint(p, color, indexBuffer, vertices);
}

void Integrator::drawLineSegment(const dvec2& v1, const dvec2& v2, const vec4& color,
                                 IndexBufferRAM* indexBuffer,
                                 std::vector<BasicMesh::Vertex>& vertices) {
    indexBuffer->add(static_cast<std::uint32_t>(vertices.size()));
    vertices.push_back({vec3(v1[0], v1[1], 0), vec3(0, 0, 1), vec3(v1[0], v1[1], 0), color});
    indexBuffer->add(static_cast<std::uint32_t>(vertices.size()));
    vertices.push_back({vec3(v2[0], v2[1], 0), vec3(0, 0, 1), vec3(v2[0], v2[1], 0), color});
}

int Integrator::StreamLines(bool Rk4type, const VectorField2& vectorField, const  dvec2& startPoint,
                            bool direction, float stepSize, bool normalized, bool stopSteps, int numIter,
                            bool stopArc, float arc, bool boundary, bool zeros, bool stopSlow, float velocity,
                            bool displayPoints, IndexBufferRAM* indexBufferPoints, IndexBufferRAM* indexBufferlines, std::vector<BasicMesh::Vertex>& vertices) {
    vec4 color;
    bool stop = false;
    
    dvec2 newPoint = vec2(0.0, 0.0);
    dvec2 currentPoint = startPoint;
    float length = 0.0f;
    int steps = 0;
    for (int i = 0; i < 200; i++) {
        
        // Find next point 
        if(Rk4type) {
            color = vec4(0, 0, 0, 1); // Black
            std::tuple<dvec2, bool> rk4Results;
            if (direction) // true means forward
                rk4Results = Integrator::RK4_(vectorField, currentPoint, stepSize, +1, normalized, zeros, stopSlow, velocity);
                //newPoint = Integrator::RK4(vectorField, currentPoint, stepSize, +1, normalized);
            else
                rk4Results = Integrator::RK4_(vectorField, currentPoint, stepSize, -1, normalized, zeros, stopSlow, velocity);
                //newPoint = Integrator::RK4(vectorField, currentPoint, stepSize, -1, normalized);
            newPoint = std::get<0>(rk4Results);
            stop = std::get<1>(rk4Results);
        }
        else{

            color = vec4(0.5, 0.5, 0.5, 1); // Grey
            if (direction) // true means forward
                newPoint = Integrator::Euler(vectorField, currentPoint, stepSize, +1, normalized);
            else
                newPoint = Integrator::Euler(vectorField, currentPoint, stepSize, -1, normalized);
        }

        // break if number of steps
        if (stopSteps && i >= numIter) {
            break;
        }

        // Break if arc length exceeded
        if (stopArc && length > arc)
        {
            break;
        }

        // Break if outside boundary
        if (boundary && !vectorField.isInside(newPoint)) {
            break;
        }

        //// Break if magnitude of newpoint is nearly zero
        //if (zeros && glm::length(newPoint) < 1e-4)
        //{
        //    break;
        //}

        //// Break if velocity is too slow
        //if (stopSlow && glm::length(newPoint - currentPoint) < velocity ) {
        //    break;
        //}
        
        // Break if vector is nearly zero or if velocity too slow
        if (stop) {
            break;
        }

        // Draw points and  line to next point
        if (displayPoints) {
            Integrator::drawPoint(newPoint, color, indexBufferPoints, vertices);
        }
        Integrator::drawLineSegment(currentPoint, newPoint, color, indexBufferlines, vertices);

        // increment for next step
        length = length + glm::length(newPoint - currentPoint);
        currentPoint = newPoint;
        steps++;
        
    }
    return steps;
}



std::list<vec2> Integrator::StreamLine(const VectorField2& vectorField, const  dvec2& startPoint, float stepSize, bool stopArc, float arc) {
    
    dvec2 newPoint = vec2(0.0, 0.0);
    dvec2 currentPoint = startPoint;
    bool stop = false;
    std::list<vec2> streamlinePoints; 
    int numIter = 500; 
    //double arc = numIter*stepSize; 
    float length = 0.0f;

    for(int i = 0; i < numIter; i++){ // FORWARD/RIGHT Curve
        std::tuple<dvec2, bool> rk4Results;
        rk4Results = Integrator::RK4_(vectorField, currentPoint, stepSize, +1, true, true, true, 1e-4);
        newPoint = std::get<0>(rk4Results);
        stop = std::get<1>(rk4Results);
        //newPoint = Integrator::RK4(vectorField, currentPoint, stepSize, +1, true);

         // Break if outside boundary
        if (!vectorField.isInside(newPoint)) {
            break;
        }

         // Break if arc length exceeded
        if (stopArc && length > arc)
        {
            break;
        }

        // Break if magnitude of newpoint is nearly zero or slow velcity
        if (stop)
        {
            break;
        }

        streamlinePoints.push_back(newPoint); 
        length = length + glm::length(newPoint - currentPoint);
        currentPoint = newPoint;
    }

    streamlinePoints.push_front(startPoint); // Add start point in middle
    
    currentPoint = startPoint;
    length = 0.0f;
    stop = false;

    for(int i = 0; i < numIter; i++){ // REVERSE/LEFT Curve
        std::tuple<dvec2, bool> rk4Results;
        rk4Results = Integrator::RK4_(vectorField, currentPoint, stepSize, -1, true, true, true, 1e-4);
        newPoint = std::get<0>(rk4Results);
        stop = std::get<1>(rk4Results);
        //newPoint = Integrator::RK4(vectorField, currentPoint, stepSize, -1, true);

         // Break if outside boundary
        if (!vectorField.isInside(newPoint)) {
            break;
        }

         // Break if arc length exceeded
        if (stopArc && length > arc)
        {
            break;
        }

        // Break if magnitude of newpoint is nearly zero
        if (stop)
        {
            break;
        }

        streamlinePoints.push_front(newPoint); 
        length = length + glm::length(newPoint - currentPoint);
        currentPoint = newPoint;
    }
    return streamlinePoints;
}

int Integrator::Separatrix(const VectorField2& vectorField, const  dvec2& startPoint, const  dvec2& criticalPoint, int dir, float stepSize,
                            IndexBufferRAM* indexBufferlines, std::vector<BasicMesh::Vertex>& vertices) {
    
    vec4 color = vec4(1, 1, 1, 1); // White
    bool stop = false;
    
    dvec2 newPoint = vec2(0.0, 0.0);
    dvec2 currentPoint = startPoint;
    float length = 0.0f;
    int steps = 0;
    Integrator::drawLineSegment(criticalPoint, startPoint, color, indexBufferlines, vertices);
    for (int i = 0; i < 200; i++) {
        
        // Find next point 
        std::tuple<dvec2, bool> rk4Results;
        rk4Results = Integrator::RK4_(vectorField, currentPoint, stepSize, dir, false, true, false, 0.001);
        newPoint = std::get<0>(rk4Results);
        stop = std::get<1>(rk4Results);

        // Break if outside boundary
        if (!vectorField.isInside(newPoint)) {
            break;
        }
        
        // Break if vector is nearly zero or if velocity too slow
        if (stop) {
            break;
        }

        // Draw line to next point
        Integrator::drawLineSegment(currentPoint, newPoint, color, indexBufferlines, vertices);

        // increment for next step
        length = length + glm::length(newPoint - currentPoint);
        currentPoint = newPoint;
        steps++;
        
    }
    
    return 0;
}

int Integrator::SeparatrixBoundary(const VectorField2& vectorField, const  dvec2& startPoint, float stepSize,
                                    IndexBufferRAM* indexBufferlines, std::vector<BasicMesh::Vertex>& vertices) {

    vec4 color = vec4(1, 1, 1, 1); // White
    bool stop = false;

    dvec2 newPoint = vec2(0.0, 0.0);
    dvec2 currentPoint = startPoint;
    float length = 0.0f;
    int steps = 0;
    for (int i = 0; i < 200; i++) {

        // Find next point 
        std::tuple<dvec2, bool> rk4Results;
        rk4Results = Integrator::RK4_(vectorField, currentPoint, stepSize, +1, false, true, false, 0.001);
        newPoint = std::get<0>(rk4Results);
        stop = std::get<1>(rk4Results);

        // Break if outside boundary
        if (!vectorField.isInside(newPoint)) {
            break;
        }

        // Break if vector is nearly zero or if velocity too slow
        if (stop) {
            break;
        }

        // Draw line to next point
        Integrator::drawLineSegment(currentPoint, newPoint, color, indexBufferlines, vertices);

        // increment for next step
        length = length + glm::length(newPoint - currentPoint);
        currentPoint = newPoint;
        steps++;

    }

    currentPoint = startPoint;
    length = 0.0f;
    stop = false;

    for (int i = 0; i < 200; i++) {

        // Find next point 
        std::tuple<dvec2, bool> rk4Results;
        rk4Results = Integrator::RK4_(vectorField, currentPoint, stepSize, -1, false, true, false, 0.001);
        newPoint = std::get<0>(rk4Results);
        stop = std::get<1>(rk4Results);

        // Break if outside boundary
        if (!vectorField.isInside(newPoint)) {
            break;
        }

        // Break if vector is nearly zero or if velocity too slow
        if (stop) {
            break;
        }

        // Draw line to next point
        Integrator::drawLineSegment(currentPoint, newPoint, color, indexBufferlines, vertices);

        // increment for next step
        length = length + glm::length(newPoint - currentPoint);
        currentPoint = newPoint;
        steps++;

    }

    return 0;
}

}  // namespace inviwo
