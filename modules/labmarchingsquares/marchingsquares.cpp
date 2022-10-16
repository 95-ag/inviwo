/*********************************************************************
 *  Author  : Himangshu Saikia, Wiebke Koepp, Anke Friederici
 *  Init    : Monday, September 11, 2017 - 12:58:42
 *
 *  Project : KTH Inviwo Modules
 *
 *  License : Follows the Inviwo BSD license model
 *********************************************************************
 */

#include <labmarchingsquares/marchingsquares.h>
#include <inviwo/core/util/utilities.h>

namespace inviwo {

// The Class Identifier has to be globally unique. Use a reverse DNS naming scheme
const ProcessorInfo MarchingSquares::processorInfo_{
    "org.inviwo.MarchingSquares",  // Class identifier
    "Marching Squares",            // Display name
    "KTH Lab",                     // Category
    CodeState::Experimental,       // Code state
    Tags::None,                    // Tags
};

// Global variable
const int KERNEL_SIZE = 5;

const ProcessorInfo MarchingSquares::getProcessorInfo() const { return processorInfo_; }

MarchingSquares::MarchingSquares()
    : Processor()
    , inData("volumeIn")
    , meshIsoOut("meshIsoOut")
    , meshGridOut("meshGridOut")
    , propShowGrid("showGrid", "Show Grid")
    , propGridColor("gridColor", "Grid Lines Color", vec4(0.0f, 0.0f, 0.0f, 1.0f), vec4(0.0f),
                    vec4(1.0f), vec4(0.1f), InvalidationLevel::InvalidOutput,
                    PropertySemantics::Color)
    , propDeciderType("deciderType", "Decider Type")
    , propRandomSeed("seed", "Random Seed", 0, 0, std::mt19937::max())
    , propMultiple("multiple", "Iso Levels")
    , propIsoValue("isovalue", "Iso Value")
    , propIsoColor("isoColor", "Color", vec4(0.0f, 0.0f, 1.0f, 1.0f), vec4(0.0f), vec4(1.0f),
                   vec4(0.1f), InvalidationLevel::InvalidOutput, PropertySemantics::Color)
    , propNumContours("numContours", "Number of Contours", 1, 1, 50, 1)
    , propIsoTransferFunc("isoTransferFunc", "Colors", &inData)
    , propSmooth("smoothenInput", "Smoothen Input" )
    , propSigma("sigma", "Sigma") {
    // Register ports
    addPort(inData);
    addPort(meshIsoOut);
    addPort(meshGridOut);

    // Register properties
    addProperty(propShowGrid);
    addProperty(propGridColor);

    addProperty(propDeciderType);
    propDeciderType.addOption("asymptotic", "Asymptotic", 0);
    propDeciderType.addOption("random", "Random", 1);

    addProperty(propRandomSeed);
    propRandomSeed.setSemantics(PropertySemantics::Text);

    addProperty(propMultiple);

    propMultiple.addOption("single", "Single", 0);
    addProperty(propIsoValue);
    addProperty(propIsoColor);

    propMultiple.addOption("multiple", "Multiple", 1);
    addProperty(propNumContours);
    addProperty(propIsoTransferFunc);
    addProperty(propSmooth);
    addProperty(propSigma);



    // The default transfer function has just two blue points
    propIsoTransferFunc.get().clear();
    //propIsoTransferFunc.get().add(0.0f, vec4(0.5f, 0.0f, 1.0f, 1.0f)); //purple
    //propIsoTransferFunc.get().add(0.25f, vec4(0.0f, 0.0f, 1.0f, 1.0f)); //blue
    //propIsoTransferFunc.get().add(0.75f, vec4(0.0f, 1.0f, 1.0f, 1.0f)); //cyan
    //propIsoTransferFunc.get().add(1.0f, vec4(0.0f, 1.0f, 0.5f, 1.0f)); //light green
    propIsoTransferFunc.get().add(0.0f, vec4(0.5f, 0.0f, 1.0f, 1.0f)); //purple
    propIsoTransferFunc.get().add(0.5f, vec4(0.0f, 0.0f, 1.0f, 1.0f)); //blue
    propIsoTransferFunc.get().add(1.0f, vec4(0.0f, 1.0f, 1.0f, 1.0f)); //cyan
    propIsoTransferFunc.setCurrentStateAsDefault();

    util::hide(propGridColor, propRandomSeed, propNumContours, propIsoTransferFunc, propSigma);

    propDeciderType.onChange([this]() {
        if (propDeciderType.get() == 1) {
            util::show(propRandomSeed);
        } else {
            util::hide(propRandomSeed);
        }
    });

    // Show the grid color property only if grid is actually displayed
    propShowGrid.onChange([this]() {
        if (propShowGrid.get()) {
            util::show(propGridColor);
        } else {
            util::hide(propGridColor);
        }
    });

    // Show the sigma property only if smoothing is selected
    propSmooth.onChange([this]() {
        if (propSmooth.get()) {
            util::show(propSigma);
        }
        else {
            util::hide(propSigma);
        }
        });

    // Show options based on display of one or multiple iso contours
    propMultiple.onChange([this]() {
        if (propMultiple.get() == 0) {
            util::show(propIsoValue, propIsoColor);
            util::hide(propNumContours, propIsoTransferFunc);
        } else {
            //util::hide(propIsoValue);
            //util::show(propIsoColor, propNumContours);

            // TODO (Bonus): Comment out above if you are using the transfer function
            // and comment in below instead
            util::hide(propIsoValue, propIsoColor);
            util::show(propNumContours, propIsoTransferFunc);
        }
    });
}

void MarchingSquares::process() {
    if (!inData.hasData()) {
        return;
    }

    // Create a structured grid from the input volume
    auto vol = inData.getData();
    auto grid = ScalarField2::createFieldFromVolume(vol);

    // Extract the minimum and maximum value from the input data
    double minValue = grid.getMinValue();
    double maxValue = grid.getMaxValue();

    // Set the range for the isovalue to that minimum and maximum
    propIsoValue.setMinValue(minValue);
    propIsoValue.setMaxValue(maxValue);

    // You can print to the Inviwo console with Log-commands:
    LogProcessorInfo("This scalar field contains values between " << minValue << " and " << maxValue
                                                                  << ".");
    // You can also inform about errors and warnings:
    // LogProcessorWarn("I am warning about something"); // Will print warning message in yellow
    // LogProcessorError("I am letting you know about an error"); // Will print error message in red
    // (There is also LogNetwork...() and just Log...(), these display a different source,
    // LogProcessor...() for example displays the name of the processor in the workspace while
    // Log...() displays the identifier of the processor (thus with multiple processors of the
    // same kind you would not know which one the information is coming from

    // Get the definition of our structured grid with
    // - number of vertices in each dimension {nx, ny}
    const ivec2 nVertPerDim = grid.getNumVerticesPerDim();
    // - bounding box {xmin, ymin} - {xmax, ymax}
    const dvec2 bBoxMin = grid.getBBoxMin();
    const dvec2 bBoxMax = grid.getBBoxMax();
    // - cell size {dx, dy}
    const dvec2 cellSize = grid.getCellSize();

    // Values at the vertex positions can be accessed by the indices of the vertex
    // with index i ranging between [0, nx-1] and j in [0, ny-1]
    ivec2 ij = {0, 0};
    double valueAt00 = grid.getValueAtVertex(ij);
    LogProcessorInfo("The value at (0,0) is: " << valueAt00 << ".");

    // Initialize the output: mesh and vertices for the grid and bounding box
    auto gridmesh = std::make_shared<BasicMesh>();
    std::vector<BasicMesh::Vertex> gridvertices;

    auto indexBufferBBox = gridmesh->addIndexBuffer(DrawType::Lines, ConnectivityType::None);
    // bottomLeft to topLeft
    drawLineSegment(bBoxMin, vec2(bBoxMin[0], bBoxMax[1]), propGridColor.get(),
                    indexBufferBBox.get(), gridvertices);
    // topLeft to topRight
    drawLineSegment(vec2(bBoxMin[0], bBoxMax[1]), bBoxMax, propGridColor.get(),
                    indexBufferBBox.get(), gridvertices);
    // topRight to bottomRight
    drawLineSegment(bBoxMax, vec2(bBoxMax[0], bBoxMin[1]), propGridColor.get(),
                    indexBufferBBox.get(), gridvertices);
    // bottomRight to bottomLeft
    drawLineSegment(vec2(bBoxMax[0], bBoxMin[1]), bBoxMin, propGridColor.get(),
                    indexBufferBBox.get(), gridvertices);

    // Set the random seed to the one selected in the interface
    randGenerator.seed(static_cast<std::mt19937::result_type>(propRandomSeed.get()));
    // You can create a random sample between min and max with
    float minRand = 0.0;
    float maxRand = 1.0;
    float rand = randomValue(minRand, maxRand);
    LogProcessorInfo("The first random sample for seed " << propRandomSeed.get() << " between "
                                                         << minRand << " and " << maxRand << " is "
                                                         << rand << ".");

    // Properties are accessed with propertyName.get()
    if (propShowGrid.get()) {
        // TODO: Add grid lines of the given color

        // The function drawLineSegments creates two vertices at the specified positions,
        // that are placed into the Vertex vector defining our mesh.
        // An index buffer specifies which of those vertices should be grouped into to make up
        // lines/trianges/quads. Here two vertices make up a line segment.
        auto indexBufferGrid = gridmesh->addIndexBuffer(DrawType::Lines, ConnectivityType::None);

        // Draw vertical line segments on the x-axis
        for (int i = 1; i <= (bBoxMax[0] - bBoxMin[0]) / cellSize[0]; i++) {
            vec2 v1 = vec2(bBoxMin[0] + (i * cellSize[0]), bBoxMin[1]);
            vec2 v2 = vec2(bBoxMin[0] + (i * cellSize[0]), bBoxMax[1]);
            drawLineSegment(v1, v2, propGridColor.get(), indexBufferGrid.get(), gridvertices);
        }
        // Draw horizanrt line segments on the y-axis
        for (int i = 1; i <= (bBoxMax[1] - bBoxMin[1]) / cellSize[1]; i++) {
            vec2 v1 = vec2(bBoxMin[0], bBoxMin[1] + (i * cellSize[1]));
            vec2 v2 = vec2(bBoxMax[0], bBoxMin[1] + (i * cellSize[1]));
            drawLineSegment(v1, v2, propGridColor.get(), indexBufferGrid.get(), gridvertices);
        }
    }

    // Set the created grid mesh as output
    gridmesh->addVertices(gridvertices);
    meshGridOut.setData(gridmesh);

    // TODO (Bonus) Gaussian filter
    // Our input is const (i.e. cannot be altered), but you need to compute smoothed data and write
    // it somewhere
    // Create an editable structured grid with ScalarField2 smoothedField =
    // ScalarField2(nVertPerDim, bBoxMin, bBoxMax - bBoxMin); Values can be set with
    // smoothedField.setValueAtVertex({0, 0}, 4.2);
    // and read again in the same way as before
    // smoothedField.getValueAtVertex(ij);
    
    
    // Create an editable structured grid
    ScalarField2 smoothedField = ScalarField2(nVertPerDim, bBoxMin, bBoxMax - bBoxMin);

    // Do smooth if selected
    if (propSmooth.get()) {
        // Copy values to new scalar
        for (int i = 0; i < nVertPerDim[0]; i++) {
            for (int j = 0; j < nVertPerDim[1]; j++) {
                smoothedField.setValueAtVertex({ i, j }, grid.getValueAtVertex({ i, j }));
            }
        }

        // build gaussian kernal filter
        double sigma = propSigma.get();
        double filter[KERNEL_SIZE];
        gaussKernel(filter, sigma);

        //// 1D convolution on x-axis
        //for (int i = 0; i < nVertPerDim[0]; i++) {
        //    for (int j = KERNEL_SIZE/2; j < nVertPerDim[1]- KERNEL_SIZE/2; j++) {
        //        double inputs[KERNEL_SIZE];
        //        for (int k = j - KERNEL_SIZE/2; k < j + KERNEL_SIZE/2; k++) {
        //            int x = KERNEL_SIZE/2 + (k - j);
        //            inputs[x] = smoothedField.getValueAtVertex({i, k});
        //        }
        //        // smoothed values
        //        auto newVal = applyConv(inputs, filter);
        //        smoothedField.setValueAtVertex({ i, j }, newVal); 
        //    }
        //}

        //// 1D convolution on y-axis
        //for (int j = 0; j < nVertPerDim[1]; j++) {
        //    for (int i = KERNEL_SIZE/2; i < nVertPerDim[0] - KERNEL_SIZE/2; i++) {
        //        double inputs[KERNEL_SIZE];
        //        for (int k = i - KERNEL_SIZE/2; k < i + KERNEL_SIZE/2; k++) {
        //            int x = KERNEL_SIZE/2 + (k - i);
        //            inputs[x] = smoothedField.getValueAtVertex({ k, j });
        //        }
        //        // smoothed values
        //        auto newVal = applyConv(inputs, filter);
        //        smoothedField.setValueAtVertex({ i, j }, newVal);
        //    }
        //}

        // CNTR WITH EDGES
        // 1D convolution on x-axis
        for (int i = 0; i < nVertPerDim[0]; i++) {
            for (int j = 0; j < nVertPerDim[1]; j++) {
                double inputs[KERNEL_SIZE];
                int cntr = 0;
                for (int k = j - KERNEL_SIZE / 2; k < j + KERNEL_SIZE / 2; k++) {
                    if (k >= 0 && k < nVertPerDim[1]) {
                        int x = KERNEL_SIZE / 2 + (k - j);
                        inputs[x] = smoothedField.getValueAtVertex({ i, k });
                        cntr++;
                    }
                }
                // smoothed values
                double newVal = 0.0;
                //apply Gauss Filter/Convolution
                for (int x = 0; x < cntr; x++) {
                    std::cout << filter[x] << " " << inputs[x];
                    newVal += filter[x] * inputs[x];
                }
                newVal = newVal / cntr;
                smoothedField.setValueAtVertex({ i, j }, newVal);
            }
        }

        // 1D convolution on y-axis
        for (int j = 0; j < nVertPerDim[1]; j++) {
            for (int i = 0; i < nVertPerDim[0]; i++) {
                double inputs[KERNEL_SIZE];
                int cntr = 0;
                for (int k = i - KERNEL_SIZE / 2; k < i + KERNEL_SIZE / 2; k++) {
                    if (k >= 0 && k < nVertPerDim[0]) {
                        int x = KERNEL_SIZE / 2 + (k - i);
                        inputs[x] = smoothedField.getValueAtVertex({ k, j });
                        cntr++;
                    }
                }
                // smoothed values
                double newVal = 0.0;
                //apply Gauss Filter/Convolution
                for (int x = 0; x < cntr; x++) {
                    newVal += filter[x] * inputs[x];
                }
                newVal = newVal / cntr;
                smoothedField.setValueAtVertex({ i, j }, newVal);
            }
        }

        minValue = smoothedField.getMinValue();
        maxValue = smoothedField.getMaxValue();
    }

    // Initialize the output: mesh and vertices
    auto mesh = std::make_shared<BasicMesh>();
    std::vector<BasicMesh::Vertex> vertices;

    auto indexBufferIsolines = mesh->addIndexBuffer(DrawType::Lines, ConnectivityType::None);

    double isovalueArray[50];
    bool transferfunc = false;
    int n;
    vec4 color;

    if (propMultiple.get() == 0) {
        // TODO: Draw a single isoline at the specified isovalue (propIsoValue)
        // and color it with the specified color (propIsoColor)
        n = 1;
        isovalueArray[0] = propIsoValue;
    }

    else {
        // TODO: Draw the given number (propNumContours) of isolines between
        // the minimum and maximum value

        // Get array of isovalues for given number of contours
        n = propNumContours;
        for (int i = 0; i < n; i++) {
            isovalueArray[i] = minValue + ((i+1) * ((maxValue - minValue) /(n+1)));
        }
        
        // TODO (Bonus): Use the transfer function property to assign a color
        // The transfer function normalizes the input data and sampling colors
        // from the transfer function assumes normalized input, that means
        // vec4 color = propIsoTransferFunc.get().sample(0.0f);
        // is the color for the minimum value in the data
        // vec4 color = propIsoTransferFunc.get().sample(1.0f);
        // is the color for the maximum value in the data
        
        // UNCOMMENT if using transfer function
        transferfunc = true;

    }
    // Find and draw all isolines
    for( int l = 0; l < n; l++){
        // get isovalue
        double isovalue = isovalueArray[l];
        LogProcessorInfo("isovalues" << isovalue << "at" << n << ".");
        // get isovalue color
        if (transferfunc) { // normalized color from cyan to blue
            double normalized_isovalue = (isovalue - minValue) / (maxValue - minValue);
            color = propIsoTransferFunc.get().sample(normalized_isovalue);
        }
        else { // given color
            color = propIsoColor.get();
        }

        //Check each cell for isoline
        for (int i = 0; i < nVertPerDim[0] - 1; i++) {
            for (int j = 0; j < nVertPerDim[1] - 1; j++) {
                // get vertex values
                double f00, f01, f10, f11;
                if (propSmooth.get()) { // if smooothened input
                    f00 = smoothedField.getValueAtVertex({ i, j });
                    f10 = smoothedField.getValueAtVertex({ i + 1, j });
                    f01 = smoothedField.getValueAtVertex({ i, j + 1 });
                    f11 = smoothedField.getValueAtVertex({ i + 1, j + 1 });
                }
                else { // normal input
                    f00 = grid.getValueAtVertex({ i, j });
                    f10 = grid.getValueAtVertex({ i + 1, j });
                    f01 = grid.getValueAtVertex({ i, j + 1 });
                    f11 = grid.getValueAtVertex({ i + 1, j + 1 });
                }
                vec2 p00 = getVertexPos(i, j, bBoxMin, cellSize);
                vec2 p10 = getVertexPos(i + 1, j, bBoxMin, cellSize);
                vec2 p01 = getVertexPos(i, j + 1, bBoxMin, cellSize);
                vec2 p11 = getVertexPos(i + 1, j + 1, bBoxMin, cellSize);
                //check if isoline passes through cell
                if (isolineInCell(f00, f01, f10, f11, isovalue)) {
                    // find interpolation points
                    int pointCount = 0;
                    vec2 intersection[4];
                    if (isolineOnSide(f00, f10, isovalue)) { //bottom line
                        auto x = linearInterpolation(p00[0], p10[0], f00, f10, isovalue);
                        intersection[pointCount] = vec2(x, p00[1]);
                        pointCount++;
                    }
                    if (isolineOnSide(f10, f11, isovalue)) { //right line
                        auto y = linearInterpolation(p10[1], p11[1], f10, f11, isovalue);
                        intersection[pointCount] = vec2(p10[0], y);
                        pointCount++;
                    }
                    if (isolineOnSide(f01, f11, isovalue)) { //top line
                        auto x = linearInterpolation(p01[0], p11[0], f01, f11, isovalue);
                        intersection[pointCount] = vec2(x, p01[1]);
                        pointCount++;
                    }
                    if (isolineOnSide(f00, f01, isovalue)) { //left line
                        auto y = linearInterpolation(p00[1], p01[1], f00, f01, isovalue);
                        intersection[pointCount] = vec2(p00[0], y);
                        pointCount++;
                    }

                    //Connect interpolation points
                    if (pointCount == 2) { // Only one sign change
                        drawLineSegment(intersection[0], intersection[1],
                            color, indexBufferIsolines.get(), vertices);
                    }
                    if (pointCount == 4) { // Two sign changes 
                        if (propDeciderType == 0) { // Assympotic Decider
                            if (intersection[0][0] < intersection[2][0]) { // check x values of top and bottom cell lines
                                drawLineSegment(intersection[2], intersection[1],
                                    color, indexBufferIsolines.get(), vertices);
                                drawLineSegment(intersection[0], intersection[3],
                                    color, indexBufferIsolines.get(), vertices);
                            }
                            else {
                                drawLineSegment(intersection[2], intersection[3],
                                    color, indexBufferIsolines.get(), vertices);
                                drawLineSegment(intersection[1], intersection[0],
                                    color, indexBufferIsolines.get(), vertices);
                            }
                        }
                        else { // Random Decider
                            if (randomValue(0, 1) < 0.5) {
                                drawLineSegment(intersection[2], intersection[1],
                                    color, indexBufferIsolines.get(), vertices);
                                drawLineSegment(intersection[0], intersection[3],
                                    color, indexBufferIsolines.get(), vertices);
                            }
                            else {
                                drawLineSegment(intersection[2], intersection[3],
                                    color, indexBufferIsolines.get(), vertices);
                                drawLineSegment(intersection[1], intersection[0],
                                    color, indexBufferIsolines.get(), vertices);
                            }

                        }
                    }
                }
            }
        }
    }

    // Note: It is possible to add multiple index buffers to the same mesh,
    // thus you could for example add one for the grid lines and one for
    // each isoline
    // Also, consider to write helper functions to avoid code duplication
    // e.g. for the computation of a single iso contour



    mesh->addVertices(vertices);
    meshIsoOut.setData(mesh);
}

void MarchingSquares::gaussKernel(double kernel[KERNEL_SIZE], double sigma) {
    // intialize sum for normalizing
    double sum = 0.0;

    // compute kernal
    for(int x = -KERNEL_SIZE /2; x <= KERNEL_SIZE / 2; x++) {
        kernel[x + 2] = (1 / (sqrt(2 * M_PI) * sigma)) * exp(-(x * x) / (2 * sigma * sigma));
        sum += kernel[x + 2];
    }

    // normalize kernal
    for (int i = 0; i < KERNEL_SIZE; i++) {
        kernel[i] = kernel[i] / sum;
    }

    return;
}

double MarchingSquares::applyConv(double inputs[KERNEL_SIZE], double kernel[KERNEL_SIZE]){
    // smoothed values
    double smoothenValue = 0.0;
    
    //apply Gauss Filter/Convolution
    for (int i = 0; i < KERNEL_SIZE; i++) {
        smoothenValue += kernel[i] * inputs[i];
    }
    smoothenValue = smoothenValue / KERNEL_SIZE;

    return smoothenValue;
}

double MarchingSquares::linearInterpolation(double x0, double x1, double f0, double f1, double c) {

    return (x0 * (c - f1) / (f0 - f1)) + (x1 * (c - f0) / (f1 - f0));
}

vec2 MarchingSquares::getVertexPos(int i, int j, dvec2 bBoxMin, dvec2 cellSize) {

    return vec2(bBoxMin[0] + (i * cellSize[0]), bBoxMin[1] + (j * cellSize[1]));
}

bool  MarchingSquares::isolineInCell(double f00, double f01, double f10, double f11, double c) {
    double fmin = std::min(std::min(f00, f01), std::min(f10, f11));
    double fmax = std::max(std::max(f00, f01), std::max(f10, f11));
    return ((fmin <= c) && (c <= fmax));
}

bool  MarchingSquares::isolineOnSide(double f0, double f1, double c) {
    double fmin = std::min(f0, f1);
    double fmax = std::max(f0, f1);
    return ((fmin <= c) && (c <= fmax));
}

float MarchingSquares::randomValue(const float min, const float max) const {
    return min + uniformReal(randGenerator) * (max - min);
}

void MarchingSquares::drawLineSegment(const vec2& v1, const vec2& v2, const vec4& color,
                                      IndexBufferRAM* indexBuffer,
                                      std::vector<BasicMesh::Vertex>& vertices) {
    // Add first vertex
    indexBuffer->add(static_cast<std::uint32_t>(vertices.size()));
    // A vertex has a position, a normal, a texture coordinate and a color
    // we do not use normal or texture coordinate, but still have to specify them
    vertices.push_back({vec3(v1[0], v1[1], 0), vec3(0, 0, 1), vec3(v1[0], v1[1], 0), color});
    // Add second vertex
    indexBuffer->add(static_cast<std::uint32_t>(vertices.size()));
    vertices.push_back({vec3(v2[0], v2[1], 0), vec3(0, 0, 1), vec3(v2[0], v2[1], 0), color});
}

}  // namespace inviwo
