/*********************************************************************
 *  Author  : Anke Friederici
 *
 *  Project : KTH Inviwo Modules
 *
 *  License : Follows the Inviwo BSD license model
 **********************************************************************/

#include <inviwo/core/datastructures/geometry/basicmesh.h>
#include <inviwo/core/datastructures/volume/volumeram.h>
#include <labstreamlines/integrator.h>
#include <labutils/scalarvectorfield.h>
#include <labtopo/topology.h>
#include <labtopo/utils/gradients.h>

namespace inviwo {

const vec4 Topology::ColorsCP[6] = {
    vec4(1, 1, 0, 1),    // Saddle - Yellow
    vec4(1, 0, 0, 1),    // AttractingNode - Red
    vec4(0, 0, 1, 1),    // RepellingNode - Blue
    vec4(0.5, 0, 1, 1),  // AttractingFocus - Purple
    vec4(1, 0.5, 0, 1),  // RepellingFocus - Orange
    vec4(0, 1, 0, 1)     // Center - Green
};

// The Class Identifier has to be globally unique. Use a reverse DNS naming scheme
const ProcessorInfo Topology::processorInfo_{
    "org.inviwo.Topology",    // Class identifier
    "Vector Field Topology",  // Display name
    "KTH Lab",                // Category
    CodeState::Experimental,  // Code state
    Tags::None,               // Tags
};

const ProcessorInfo Topology::getProcessorInfo() const { return processorInfo_; }

Topology::Topology()
    : Processor()
    , inData("inData")
    , outMesh("meshOut")
    , meshBBoxOut("meshBBoxOut")
// TODO: Initialize additional properties
// propertyName("propertyIdentifier", "Display Name of the Propery",
// default value (optional), minimum value (optional), maximum value (optional), increment
// (optional)); propertyIdentifier cannot have spaces
    ,propThresholdFactor("thresholdFactor", "Threshold Factor", 0.1, 0.001, 1, 0.001)
    , propThresholdBoundary("thresholdBoundary", "Threshold Boundary", 0.1, 0.001, 1, 0.001)
{
    // Register Ports
    addPort(outMesh);
    addPort(inData);
    addPort(meshBBoxOut);

    // TODO: Register additional properties
    // addProperty(propertyName);
    addProperty(propThresholdFactor); 
    addProperty(propThresholdBoundary);
}



void Topology::process() {
    // Get input
    if (!inData.hasData()) {
        return;
    }
    auto vol = inData.getData();

    // Retreive data in a form that we can access it
    const VectorField2 vectorField = VectorField2::createFieldFromVolume(vol);

    // Add a bounding box to the mesh
    const dvec2& BBoxMin = vectorField.getBBoxMin();
    const dvec2& BBoxMax = vectorField.getBBoxMax();
    const dvec2 cellSize = vectorField.getCellSize();
    auto bboxMesh = std::make_shared<BasicMesh>();
    std::vector<BasicMesh::Vertex> bboxVertices;
    auto indexBufferBBox = bboxMesh->addIndexBuffer(DrawType::Lines, ConnectivityType::Strip);
    // Bounding Box vertex 0
    vec4 black = vec4(0, 0, 0, 1);
    Integrator::drawNextPointInPolyline(BBoxMin, black, indexBufferBBox.get(), bboxVertices);
    Integrator::drawNextPointInPolyline(vec2(BBoxMin[0], BBoxMax[1]), black, indexBufferBBox.get(),
                                        bboxVertices);
    Integrator::drawNextPointInPolyline(BBoxMax, black, indexBufferBBox.get(), bboxVertices);
    Integrator::drawNextPointInPolyline(vec2(BBoxMax[0], BBoxMin[1]), black, indexBufferBBox.get(),
                                        bboxVertices);
    // Connect back to the first point, to make a full rectangle
    indexBufferBBox->add(static_cast<std::uint32_t>(0));
    bboxMesh->addVertices(bboxVertices);
    meshBBoxOut.setData(bboxMesh);

    // Initialize mesh, vertices and index buffers for seperatrices
    auto mesh = std::make_shared<BasicMesh>();
    std::vector<BasicMesh::Vertex> vertices;
    // Either add all line segments to this index buffer (one large buffer, two consecutive points
    // make up one line), or use several index buffers with connectivity type strip.
    auto indexBufferSeparatrices = mesh->addIndexBuffer(DrawType::Lines, ConnectivityType::None);
    // auto indexBufferSeparatrices = mesh->addIndexBuffer(DrawType::Lines,
    // ConnectivityType::Strip);

    auto indexBufferPoints = mesh->addIndexBuffer(DrawType::Points, ConnectivityType::None);

    // TODO: Compute the topological skeleton of the input vector field.
    // Find the critical points and color them according to their type.
    // Integrate all separatrices.

    size2_t dims = vectorField.getNumVerticesPerDim();
    dvec2 threshold = dvec2(propThresholdFactor, propThresholdFactor);        // CHECK LATER
    double stepSize = 0.05;

    std::cout << "BBoxMin_x: " << BBoxMin[0] << " BBoxMin_y: " << BBoxMin[1] << " dims_x: " << dims[0] << " dims_y: " << dims[1] << "\n";

    // Find Critical Points
    std::list<dvec2> criticalPoints; 
    for (size_t j = 0; j < dims[1]-1; ++j) {
        for (size_t i = 0; i < dims[0]-1; ++i) {
            dvec2 f00, f01, f10, f11, p00, p10, p01, p11;
            std::tuple<dvec2, bool> zeroPoint; 
            f00 = vectorField.getValueAtVertex({ i, j });
            f10 = vectorField.getValueAtVertex({ i + 1, j });
            f01 = vectorField.getValueAtVertex({ i, j + 1 });
            f11 = vectorField.getValueAtVertex({ i + 1, j + 1 }); 
            p00 = getVertexPos(i, j, BBoxMin, cellSize);
            p10 = getVertexPos(i + 1, j, BBoxMin, cellSize);
            p01 = getVertexPos(i, j + 1, BBoxMin, cellSize);
            p11 = getVertexPos(i + 1, j + 1, BBoxMin, cellSize);
            //Integrator::drawPoint(p00, vec4(0,0,0,1), indexBufferPoints.get(),vertices);                     
            if(checkSign(f00, f10, f01, f11)){
                //zeroPoint = domainDecomposition(vectorField, threshold, p00, p10, p01, p11);
                zeroPoint = domainDecomposition2(vectorField, threshold, p00, p10, p01, p11);
                if(std::get<1>(zeroPoint)){
                    criticalPoints.push_back(std::get<0>(zeroPoint));
                    std::cout << "\n" << criticalPoints.back();
                }
            }
        }
    }
    
    LogProcessorInfo("Critical: "<< criticalPoints.size());

   
    // Find all boundary switch points
    std::list<dvec2> boundarySwitchPoints;
    // Bottom & top
    for (size_t j = 0; j < dims[1]; j = j + dims[1] - 1) {
        for (size_t i = 0; i < dims[0]-1; i = i + 1) {
            dvec2 f0, f1, p0, p1;
            std::tuple<dvec2, bool> zeroPoint;
            f0 = vectorField.getValueAtVertex({ i, j });
            f1 = vectorField.getValueAtVertex({ i + 1, j });
            p0 = getVertexPos(i, j, BBoxMin, cellSize);
            p1 = getVertexPos(i + 1, j, BBoxMin, cellSize);
            if (checkSignBoundary(f0, f1)) {
                zeroPoint = domainDecompositionBoundary(vectorField, propThresholdBoundary, p0, p1, true);
                if (std::get<1>(zeroPoint)) {
                    boundarySwitchPoints.push_back(std::get<0>(zeroPoint));
                    LogProcessorInfo("\n" << boundarySwitchPoints.back());
                }
            }
        }
    }
    // Left & Right
    for (size_t j = 0; j < dims[1]-1; j = j + 1) {
        for (size_t i = 0; i < dims[0]; i = i + dims[0] - 1) {
            dvec2 f0, f1, p0, p1;
            std::tuple<dvec2, bool> zeroPoint;
            f0 = vectorField.getValueAtVertex({ i, j });
            f1 = vectorField.getValueAtVertex({ i, j + 1 });
            p0 = getVertexPos(i, j, BBoxMin, cellSize);
            p1 = getVertexPos(i, j + 1, BBoxMin, cellSize);
            if (checkSignBoundary(f0, f1)) {
                zeroPoint = domainDecompositionBoundary(vectorField, propThresholdBoundary, p0, p1, false);
                if (std::get<1>(zeroPoint)) {
                    boundarySwitchPoints.push_back(std::get<0>(zeroPoint));
                    LogProcessorInfo("\n" << boundarySwitchPoints.back());
                }
            }
        }
    }

    LogProcessorInfo("Boundary Switch: " << boundarySwitchPoints.size());

    // Mark all boundary switch points
    for (auto it : boundarySwitchPoints) {
        Integrator::drawPoint(it, vec4(0.8, 0.8, 0.8, 1), indexBufferPoints.get(), vertices);
        Integrator::SeparatrixBoundary(vectorField, it, stepSize, indexBufferSeparatrices.get(), vertices);
    }

    // Accessing the colors
    vec4 colorCenter = ColorsCP[static_cast<int>(TypeCP::Center)];
    std::vector<std::tuple<dvec2, int>> criticalPointsType(criticalPoints.size());

    // Analysis of critical point and draw Separatices for critical points

    int counter = 0;
    for (auto it:criticalPoints){
        // Computing the jacobian at a position
        dmat2 jacobian = vectorField.derive(it);
        // Doing the eigen analysis
        auto EigenResult = util::eigenAnalysis(jacobian);
        vec2 im = EigenResult.eigenvaluesIm;
        vec2 re = EigenResult.eigenvaluesRe;
        mat2 eigenvector = EigenResult.eigenvectors;
        dvec2 eivec1 = dvec2(eigenvector[0][0], eigenvector[0][1]);
        dvec2 eivec2 = dvec2(eigenvector[1][0], eigenvector[1][1]);
        //("eigenvector: " << eigenvector << "\n");
        if(im[0] == 0 && im[1] == 0 ){
            // repelling node
            if(re[0] > 0 && re[1] > 0){
                criticalPointsType[counter] = std::tuple(it, 2);
            }
            // attracting node
            else if(re[0] < 0 && re[1] < 0){
                criticalPointsType[counter] = std::tuple(it, 1);
            }
            // saddle point
            else{
                criticalPointsType[counter] = std::tuple(it, 0);
                //Draw separatrices
                std::vector<dvec2> seedPoints(4);
                seedPoints[0] = it + stepSize * eivec1;
                seedPoints[1] = it - stepSize * eivec1;
                seedPoints[2] = it + stepSize * eivec2;
                seedPoints[3] = it - stepSize * eivec2;
                if (re[0] < 0) { // eigen1 corresponds to incoming
                    Integrator::Separatrix(vectorField, seedPoints[0], it, -1, stepSize, indexBufferSeparatrices.get(), vertices);
                    Integrator::Separatrix(vectorField, seedPoints[1], it, -1, stepSize, indexBufferSeparatrices.get(), vertices);
                    Integrator::Separatrix(vectorField, seedPoints[2], it, +1, stepSize, indexBufferSeparatrices.get(), vertices);
                    Integrator::Separatrix(vectorField, seedPoints[3], it, +1, stepSize, indexBufferSeparatrices.get(), vertices);
                }
                else {// eigen2 corresponds to incoming
                    Integrator::Separatrix(vectorField, seedPoints[0], it, +1, stepSize, indexBufferSeparatrices.get(), vertices);
                    Integrator::Separatrix(vectorField, seedPoints[1], it, +1, stepSize, indexBufferSeparatrices.get(), vertices);
                    Integrator::Separatrix(vectorField, seedPoints[2], it, -1, stepSize, indexBufferSeparatrices.get(), vertices);
                    Integrator::Separatrix(vectorField, seedPoints[3], it, -1, stepSize, indexBufferSeparatrices.get(), vertices);
                }
                
            }
        }
        else {
            // repelling focus
            if(re[0] > 0 && re[1] > 0){
                criticalPointsType[counter] = std::tuple(it, 4);

            }
            // attracting focus
            else if(re[0] < 0 && re[1] < 0){
                criticalPointsType[counter] = std::tuple(it, 3);
            }
            // center
            else{
                criticalPointsType[counter] = std::tuple(it, 5);
            }
        }
        counter++;
    }
    
    // Mark all critical points
    for(int i = 0; i < criticalPoints.size(); i++){
         Integrator::drawPoint(std::get<0>(criticalPointsType[i]), ColorsCP[std::get<1>(criticalPointsType[i])], indexBufferPoints.get(),vertices);
    }

    





    mesh->addVertices(vertices);
    outMesh.setData(mesh);
}

void Topology::drawLineSegment(const dvec2& v1, const dvec2& v2, const vec4& color,
                               IndexBufferRAM* indexBuffer,
                               std::vector<BasicMesh::Vertex>& vertices) {
    indexBuffer->add(static_cast<std::uint32_t>(vertices.size()));
    vertices.push_back({vec3(v1[0], v1[1], 0), vec3(0, 0, 1), vec3(v1[0], v1[1], 0), color});
    indexBuffer->add(static_cast<std::uint32_t>(vertices.size()));
    vertices.push_back({vec3(v2[0], v2[1], 0), vec3(0, 0, 1), vec3(v2[0], v2[1], 0), color});
}

bool Topology::checkSign(dvec2 f00, dvec2 f10, dvec2 f01, dvec2 f11){

    // x-axis 
    if ((f00[0] > 0 && (f10[0] < 0 || f01[0] < 0 || f11[0] < 0)) || 
        (f00[0] < 0 && (f10[0] > 0 || f01[0] > 0 || f11[0] > 0))){
            // y-axis
            if ((f00[1] > 0 && (f10[1] < 0 || f01[1] < 0 || f11[1] < 0)) || 
                (f00[1] < 0 && (f10[1] > 0 || f01[1] > 0 || f11[1] > 0))){
                    return true;        // sign-change
                }
    }
    return false;       // no sign-change
}

bool Topology::checkSignBoundary(dvec2 f0, dvec2 f1) {

    // x-axis 
    if ((f0[0] > 0 && f1[0] < 0) || (f0[0] < 0 && f1[0] > 0)) {
        // y-axis
        if ((f0[1] > 0 && f1[1] < 0) || (f0[1] < 0 && f1[1] > 0)) {
            return true;        // sign-change
        }
    }
    return false;       // no sign-change
}

dvec2 Topology::getVertexPos(int i, int j, dvec2 bBoxMin, dvec2 cellSize) {

    return dvec2(bBoxMin[0] + (i * cellSize[0]), bBoxMin[1] + (j * cellSize[1]));
}

std::tuple<dvec2, bool> Topology::domainDecomposition(const VectorField2& vectorField, dvec2 threshold, dvec2 p00, dvec2 p10, dvec2 p01, dvec2 p11){
    
    dvec2 f00, f10, f01, f11;
    dvec2 new_p00, new_p10, new_p01, new_p11; 
    dvec2 bottom_mid = dvec2((p00[0] +p10[0])/2, p00[1]);
    dvec2 top_mid = dvec2((p01[0]+p11[0])/2, p01[1]);
    dvec2 left_mid = dvec2(p00[0], (p00[1]+p01[1])/2);
    dvec2 right_mid = dvec2(p10[0],(p10[1]+p11[1])/2);
    dvec2 mid_mid = dvec2((p00[0]+p10[0])/2, (p00[1]+p01[1])/2);

    if(abs(p00[0]-p10[0]) <= threshold[0] && abs(p00[1]-p01[1]) <= threshold[1]){
            return std::tuple(mid_mid, true);
    }

    else{
        // bottom left
        new_p00 = p00;
        new_p10 = bottom_mid;
        new_p01 = left_mid;
        new_p11 = mid_mid;
        f00 = vectorField.interpolate(new_p00);
        f10 = vectorField.interpolate(new_p10); 
        f01 = vectorField.interpolate(new_p01); 
        f11 = vectorField.interpolate(new_p11);
        if (checkSign(f00, f10, f01, f11)){
            return domainDecomposition(vectorField, threshold, new_p00, new_p10, new_p01, new_p11);
        }
        else{
            // bottom right
            new_p00 = bottom_mid;
            new_p10 = p10;
            new_p01 = mid_mid;
            new_p11 = right_mid;
            f00 = vectorField.interpolate(new_p00);
            f10 = vectorField.interpolate(new_p10); 
            f01 = vectorField.interpolate(new_p01); 
            f11 = vectorField.interpolate(new_p11);
            if (checkSign(f00, f10, f01, f11)){
                return domainDecomposition(vectorField, threshold, new_p00, new_p10, new_p01, new_p11);
            }  
            else {
                // top right
                new_p00 = mid_mid;
                new_p10 = right_mid;
                new_p01 = top_mid;
                new_p11 = p11;  
                f00 = vectorField.interpolate(new_p00);
                f10 = vectorField.interpolate(new_p10); 
                f01 = vectorField.interpolate(new_p01); 
                f11 = vectorField.interpolate(new_p11);
                if (checkSign(f00, f10, f01, f11)){
                    return domainDecomposition(vectorField, threshold, new_p00, new_p10, new_p01, new_p11);
                }  
                else{
                    // top left
                    new_p00 = left_mid;
                    new_p10 = mid_mid;
                    new_p01 = p01;
                    new_p11 = top_mid;
                    f00 = vectorField.interpolate(new_p00);
                    f10 = vectorField.interpolate(new_p10); 
                    f01 = vectorField.interpolate(new_p01); 
                    f11 = vectorField.interpolate(new_p11);
                    if (checkSign(f00, f10, f01, f11)){
                        return domainDecomposition(vectorField, threshold, new_p00, new_p10, new_p01, new_p11);
                    } 
                    else {
                        return std::tuple(mid_mid, false);
                    }
            
                }
            }
            
        }

    }

}


std::tuple<dvec2, bool> Topology::domainDecomposition2(const VectorField2& vectorField, dvec2 threshold, dvec2 p00, dvec2 p10, dvec2 p01, dvec2 p11) {

    dvec2 f00, f10, f01, f11;
    dvec2 new_p00, new_p10, new_p01, new_p11;
    dvec2 bottom_mid = dvec2((p00[0] + p10[0]) / 2, p00[1]);
    dvec2 top_mid = dvec2((p01[0] + p11[0]) / 2, p01[1]);
    dvec2 left_mid = dvec2(p00[0], (p00[1] + p01[1]) / 2);
    dvec2 right_mid = dvec2(p10[0], (p10[1] + p11[1]) / 2);
    dvec2 mid_mid = dvec2((p00[0] + p10[0]) / 2, (p00[1] + p01[1]) / 2);

    // Check if below threshold
    if (abs(p00[0] - p10[0]) <= threshold[0] && abs(p00[1] - p01[1]) <= threshold[1]) {
        return std::tuple(mid_mid, true);
    }

    else {
        std::tuple<dvec2, bool> criticalPoint;
        criticalPoint = std::tuple(mid_mid, false);


        // bottom left
        new_p00 = p00;
        new_p10 = bottom_mid;
        new_p01 = left_mid;
        new_p11 = mid_mid;
        f00 = vectorField.interpolate(new_p00);
        f10 = vectorField.interpolate(new_p10);
        f01 = vectorField.interpolate(new_p01);
        f11 = vectorField.interpolate(new_p11);
        if (checkSign(f00, f10, f01, f11)) {
            criticalPoint = domainDecomposition(vectorField, threshold, new_p00, new_p10, new_p01, new_p11);
        }

        // bottom right
        if (std::get<1>(criticalPoint) == false) {// critical point not found already
            new_p00 = bottom_mid;
            new_p10 = p10;
            new_p01 = mid_mid;
            new_p11 = right_mid;
            f00 = vectorField.interpolate(new_p00);
            f10 = vectorField.interpolate(new_p10);
            f01 = vectorField.interpolate(new_p01);
            f11 = vectorField.interpolate(new_p11);
            if (checkSign(f00, f10, f01, f11)) {
                criticalPoint = domainDecomposition(vectorField, threshold, new_p00, new_p10, new_p01, new_p11);
            }
        }

        // top right
        if (std::get<1>(criticalPoint) == false) {// critical point not found already
            new_p00 = mid_mid;
            new_p10 = right_mid;
            new_p01 = top_mid;
            new_p11 = p11;
            f00 = vectorField.interpolate(new_p00);
            f10 = vectorField.interpolate(new_p10);
            f01 = vectorField.interpolate(new_p01);
            f11 = vectorField.interpolate(new_p11);
            if (checkSign(f00, f10, f01, f11)) {
                criticalPoint = domainDecomposition(vectorField, threshold, new_p00, new_p10, new_p01, new_p11);
            }
        }

        // top left
        if (std::get<1>(criticalPoint) == false) {// critical point not found already

            new_p00 = left_mid;
            new_p10 = mid_mid;
            new_p01 = p01;
            new_p11 = top_mid;
            f00 = vectorField.interpolate(new_p00);
            f10 = vectorField.interpolate(new_p10);
            f01 = vectorField.interpolate(new_p01);
            f11 = vectorField.interpolate(new_p11);
            if (checkSign(f00, f10, f01, f11)) {
                criticalPoint =  domainDecomposition(vectorField, threshold, new_p00, new_p10, new_p01, new_p11);
            }
        }
         return criticalPoint;
    }

}

std::tuple<dvec2, bool> Topology::domainDecompositionBoundary(const VectorField2& vectorField, double threshold, dvec2 p0, dvec2 p1, bool x_axis) {

    dvec2 f0, f1;
    dvec2 new_p0, new_p1, mid;
    double diff;
    if (x_axis) {
        mid = dvec2((p0[0] + p1[0]) / 2, p0[1]);
        diff = abs(p0[0] - p1[0]) - threshold;

    }
    else
    {
       mid = dvec2(p0[0], (p0[1] + p1[1]) / 2);
       diff = abs(p0[1] - p1[1]) - threshold;
    }

    // Check if below threshold
    if (diff < 0) {
        return std::tuple(mid, true);
    }

    else {
        std::tuple<dvec2, bool> criticalPoint;
        criticalPoint = std::tuple(mid, false);

        // left/bootom
        new_p0 = p0;
        new_p1 = mid;
        f0 = vectorField.interpolate(new_p0);
        f1 = vectorField.interpolate(new_p1);
        if (checkSignBoundary(f0, f1)) {
            criticalPoint = domainDecompositionBoundary(vectorField, threshold, new_p0, new_p1, x_axis);
        }

        // right/top
        if (std::get<1>(criticalPoint) == false) {// critical point not found already
            new_p0 = mid;
            new_p1 = p1;
            f0 = vectorField.interpolate(new_p0);
            f1 = vectorField.interpolate(new_p1);
            if (checkSignBoundary(f0, f1)) {
                criticalPoint = domainDecompositionBoundary(vectorField, threshold, new_p0, new_p1, x_axis);
            }
        }
        return criticalPoint;
    }

}


}  // namespace inviwo
