/*********************************************************************
 *  Author  : Himangshu Saikia
 *  Init    : Tuesday, September 19, 2017 - 15:08:33
 *
 *  Project : KTH Inviwo Modules
 *
 *  License : Follows the Inviwo BSD license model
 *********************************************************************
 */

#include <inviwo/core/interaction/events/mouseevent.h>
#include <inviwo/core/util/utilities.h>
#include <labstreamlines/integrator.h>
#include <labstreamlines/streamlineintegrator.h>
#include <labutils/scalarvectorfield.h>

namespace inviwo {

// The Class Identifier has to be globally unique. Use a reverse DNS naming
// scheme
const ProcessorInfo StreamlineIntegrator::processorInfo_{
    "org.inviwo.StreamlineIntegrator",  // Class identifier
    "Streamline Integrator",            // Display name
    "KTH Lab",                          // Category
    CodeState::Experimental,            // Code state
    Tags::None,                         // Tags
};

const ProcessorInfo StreamlineIntegrator::getProcessorInfo() const { return processorInfo_; }

StreamlineIntegrator::StreamlineIntegrator()
    : Processor()
    , inData("volIn")
    , meshOut("meshOut")
    , meshBBoxOut("meshBBoxOut")
    , propDisplayPoints("displayPoints", "Display Points", true)
    , propStartPoint("startPoint", "Start Point", vec2(0.5, 0.5), vec2(-1), vec2(1), vec2(0.1))
    , propSeedMode("seedMode", "Seeds")
    , propNumStepsTakenRK4("numstepstakenRK4", "RK4:Actual Steps", 0, 0, 100000)
    , propNumStepsTakenEuler("numstepstakenEuler", "Euler:Actual steps", 0, 0, 100000)
    , mouseMoveStart(
          "mouseMoveStart", "Move Start", [this](Event* e) { eventMoveStart(e); },
          MouseButton::Left, MouseState::Press | MouseState::Move)
// TODO: Initialize additional properties
// propertyName("propertyIdentifier", "Display Name of the Propery",
// default value (optional), minimum value (optional), maximum value (optional),
// increment (optional)); propertyIdentifier cannot have spaces
    , propFwdDir("fwdDir", "Forward Direction", true)
    , propStepSize("stepSize", "Step Size", 0.05, 0.005, 1, 0.005)
    , propNormalized("normalized", "Normalized Direction Field", false)
    , propStopSteps("stopSteps", "Stop after certain steps", false)
    , propNumSteps("numSteps", "No. of steps", 50, 1, 200, 1)
    , propStopArc("stopArc", "Stop after certain arc length", false)
    , propArc("arc", "Arc Length", 1.0, 0.1, 10.0, 0.01 )
    , propBounds("bounds", "Stop at boundary", false)
    , propZeros("zeros", "Stop at zero vector", false)
    , propStopSlow("stopSlow", "Stop at slow velocity", false)
    , propVelocity("velocity", "Velocity", 0.01, 0.001, 1.0, 0.0001)
    , propNumStreamLines("numStreamLines", "Stream Lines", 5, 1, 1000, 1)
    , propSeedUniform("seedUniform", "Seed Uniformly", false)
    , propNumXgrid("numXgrid", "X-grid", 2, 1, 50, 1)
    , propNumYgrid("numYgrid", "Y-grid", 2, 1, 50, 1)
    , propDistributed("distributed", "Distribute Seeds based on magintude", false)
    , propRk4("rk4", "Display RK4", true)
    , propEuler("euler", "Display Euler", true)
    , propRandomSeed("seed", "Random Seed", 0, 0, std::mt19937::max())
{
    // Register Ports
    addPort(inData);
    addPort(meshOut);
    addPort(meshBBoxOut);

    // Register Properties
    propSeedMode.addOption("one", "Single Start Point", 0);
    propSeedMode.addOption("multiple", "Multiple Seeds", 1);
    addProperty(propSeedMode);
    addProperty(propStartPoint);
    addProperty(propDisplayPoints);    
    addProperty(propRk4);
    addProperty(propEuler);
    addProperty(propNumStepsTakenRK4);
    propNumStepsTakenRK4.setReadOnly(true);
    propNumStepsTakenRK4.setSemantics(PropertySemantics::Text);
    addProperty(propNumStepsTakenEuler);
    propNumStepsTakenEuler.setReadOnly(true);
    propNumStepsTakenEuler.setSemantics(PropertySemantics::Text);
    addProperty(mouseMoveStart);

    // TODO: Register additional properties
    // addProperty(propertyName);
    addProperty(propStepSize);
    addProperty(propFwdDir);
    addProperty(propNormalized);
    addProperty(propStopSteps);
    addProperty(propNumSteps);
    addProperty(propStopArc);
    addProperty(propArc);
    addProperty(propBounds);
    addProperty(propZeros);
    addProperty(propStopSlow);
    addProperty(propVelocity);
    addProperty(propSeedUniform);
    addProperty(propRandomSeed);
    propRandomSeed.setSemantics(PropertySemantics::Text);
    addProperty(propNumStreamLines);
    addProperty(propNumXgrid);
    addProperty(propNumYgrid);
    addProperty(propDistributed);



    // Show properties for a single seed and hide properties for multiple seeds
    // (TODO)
    util::hide(propNumSteps, propArc, propVelocity);
    util::hide(propNumStreamLines, propRandomSeed, propSeedUniform, propDistributed);
    util::hide(propNumXgrid, propNumYgrid);
    propSeedMode.onChange([this]() {
        if (propSeedMode.get() == 0) {
            util::show(propStartPoint, mouseMoveStart, propNumStepsTakenRK4, propNumStepsTakenEuler);
            util::show(propRk4, propEuler, propFwdDir, propStepSize, propNormalized);
            util::show(propStopSteps, propStopArc, propBounds, propZeros, propStopSlow);
            util::hide(propNumSteps, propArc, propVelocity);
            util::hide(propNumStreamLines, propRandomSeed, propSeedUniform, propDistributed);
            util::hide(propNumXgrid, propNumYgrid);
        }
        else {
            util::hide(propStartPoint, mouseMoveStart, propNumStepsTakenRK4, propNumStepsTakenEuler);
            util::show(propRk4, propEuler, propFwdDir, propStepSize, propNormalized);
            util::hide(propStopSteps, propStopArc, propBounds, propZeros, propStopSlow);
            util::hide(propNumSteps, propArc, propVelocity);
            util::show(propNumStreamLines, propRandomSeed, propSeedUniform, propDistributed);
            util::hide(propNumXgrid, propNumYgrid);
        }
       
    });

    propStopSteps.onChange([this]() {
        if (propStopSteps.get() == true) {
            util::show(propNumSteps);
        }
        else {
            util::hide(propNumSteps);
        }
    });

    propStopArc.onChange([this]() {
        if (propStopArc.get() == true) {
            util::show(propArc);
        }
        else {
            util::hide(propArc);
        }
    });

    propStopSlow.onChange([this]() {
        if (propStopSlow.get() == true) {
            util::show(propVelocity);
        }
        else {
            util::hide(propVelocity);
        }
    });
    propSeedUniform.onChange([this]() {
        if (propSeedUniform.get() == true) {
            util::show(propNumXgrid, propNumYgrid);
            util::hide(propNumStreamLines, propRandomSeed, propDistributed);
        }
        else {
            util::hide(propNumXgrid, propNumYgrid);
            util::show(propNumStreamLines, propRandomSeed, propDistributed);
        }
    });
}

void StreamlineIntegrator::eventMoveStart(Event* event) {
    if (!inData.hasData()) return;
    auto mouseEvent = static_cast<MouseEvent*>(event);
    vec2 mousePos = mouseEvent->posNormalized();

    // Map to bounding box range
    mousePos[0] *= static_cast<float>(BBoxMax_[0] - BBoxMin_[0]);
    mousePos[1] *= static_cast<float>(BBoxMax_[1] - BBoxMin_[1]);
    mousePos += static_cast<vec2>(BBoxMin_);

    // Update starting point
    propStartPoint.set(mousePos);
    event->markAsUsed();
}

void StreamlineIntegrator::process() {
    // Get input
    if (!inData.hasData()) {
        return;
    }
    auto vol = inData.getData();

    // Retreive data in a form that we can access it
    auto vectorField = VectorField2::createFieldFromVolume(vol);
    BBoxMin_ = vectorField.getBBoxMin();
    BBoxMax_ = vectorField.getBBoxMax();
    const ivec2 nVertPerDim = vectorField.getNumVerticesPerDim();
    const dvec2 cellSize = vectorField.getCellSize();
    double minVal = glm::length(vectorField.getMinValue());
    double maxVal = glm::length(vectorField.getMaxValue());
    
    // The start point should be inside the volume (set maximum to the upper right corner)
    propStartPoint.setMinValue(BBoxMin_ - dvec2(1, 1));
    propStartPoint.setMaxValue(BBoxMax_ + dvec2(1, 1));

    auto bboxMesh = std::make_shared<BasicMesh>();
    std::vector<BasicMesh::Vertex> bboxVertices;

    // Make bounding box without vertex duplication, instead of line segments which duplicate
    // vertices, create line segments between each added points with connectivity type of the index
    // buffer
    auto indexBufferBBox = bboxMesh->addIndexBuffer(DrawType::Lines, ConnectivityType::Strip);
    // Bounding Box vertex 0
    vec4 black = vec4(0, 0, 0, 1);
    Integrator::drawNextPointInPolyline(BBoxMin_, black, indexBufferBBox.get(), bboxVertices);
    Integrator::drawNextPointInPolyline(vec2(BBoxMin_[0], BBoxMax_[1]), black,
                                        indexBufferBBox.get(), bboxVertices);
    Integrator::drawNextPointInPolyline(BBoxMax_, black, indexBufferBBox.get(), bboxVertices);
    Integrator::drawNextPointInPolyline(vec2(BBoxMax_[0], BBoxMin_[1]), black,
                                        indexBufferBBox.get(), bboxVertices);
    // Connect back to the first point, to make a full rectangle
    indexBufferBBox->add(static_cast<std::uint32_t>(0));
    bboxMesh->addVertices(bboxVertices);
    meshBBoxOut.setData(bboxMesh);

    auto mesh = std::make_shared<BasicMesh>();
    std::vector<BasicMesh::Vertex> vertices;

    auto indexBufferPoints = mesh->addIndexBuffer(DrawType::Points, ConnectivityType::None);


    // Set the random seed to the one selected in the interface
    randGenerator.seed(static_cast<std::mt19937::result_type>(propRandomSeed.get()));

    if (propSeedMode.get() == 0) {
        auto indexBufferEuler = mesh->addIndexBuffer(DrawType::Lines, ConnectivityType::Strip);
        auto indexBufferRK = mesh->addIndexBuffer(DrawType::Lines, ConnectivityType::Strip);
        vec2 startPoint = propStartPoint.get();
        // Clamp if outside boundary
        startPoint = vectorField.clampPositionToBBox(startPoint);

        
        
        // Draw start point
        if (propDisplayPoints.get() != 0)
            Integrator::drawPoint(startPoint, vec4(0, 0, 0, 1), indexBufferPoints.get(), vertices);

        // TODO: Create one stream line from the given start point

        // TODO: Use the propNumStepsTaken property to show how many steps have actually been
        // integrated This could be different from the desired number of steps due to stopping
        // conditions (too slow, boundary, ...)
        int stepsRk4 = 0, stepsEuler = 0;
        propNumStepsTakenRK4.set(stepsRk4);
        propNumStepsTakenEuler.set(stepsEuler);

        // TODO: Implement the Euler and Runge-Kutta of 4th order integration schemes
        // and then integrate forward for a specified number of integration steps and a given stepsize
        // (these should be additional properties of the processor)
        //Integrator::StreamLines(bool Rk4type, const VectorField2 & vectorField, const dvec2 & startPoint,
        //        bool direction, float stepSize, bool normalized, bool stopSteps, int numIter,
        //        bool stopArc, float arc, bool boundary, bool zeros, bool stopSlow, float velocity,
        //        bool displayPoints, IndexBufferRAM* indexBufferPoints, IndexBufferRAM * indexBuffer, std::vector<BasicMesh::Vertex>&vertices);
        if (propRk4.get()) {
            stepsRk4 = Integrator::StreamLines(true, vectorField, startPoint,
                propFwdDir, propStepSize, propNormalized, propStopSteps, propNumSteps,
                propStopArc, propArc, propBounds, propZeros, propStopSlow, propVelocity,
                propDisplayPoints.get() != 0, indexBufferPoints.get(), indexBufferRK.get(), vertices);
        }

        if (propEuler.get()) {
            stepsEuler = Integrator::StreamLines(false, vectorField, startPoint,
                propFwdDir, propStepSize, propNormalized, propStopSteps, propNumSteps,
                propStopArc, propArc, propBounds, propZeros, propStopSlow, propVelocity,
                propDisplayPoints.get() != 0, indexBufferPoints.get(), indexBufferEuler.get(), vertices);
        }
        
        propNumStepsTakenRK4.set(stepsRk4);
        propNumStepsTakenEuler.set(stepsEuler);

    } else {
        // TODO: Seed multiple stream lines either randomly or using a uniform grid
        // (TODO: Bonus, sample randomly according to magnitude of the vector field)
        
        int stepsRk4 = 0, stepsEuler = 0;
        if (propSeedUniform.get() != true ) {
            std::list<vec2> startPoints;
           if (propDistributed.get()) { // Distributed points
                std::list<double> magnitudes;
                std::list<int> prob;
                // find magnitudes at each pos 
                for (int i = 0; i < nVertPerDim[0] ; i++) {
                    for (int j = 0; j < nVertPerDim[1]; j++) {
                        magnitudes.push_back(glm::length(vectorField.getValueAtVertex({ i, j })));
                    }
                }
                // list them as integer probabilities
                while (!magnitudes.empty()) {
                    prob.push_back(int(1000 * ((magnitudes.front() - minVal) / (maxVal - minVal))));
                    magnitudes.pop_front();
                }
                std::discrete_distribution<int> distribution(prob.begin(), prob.end());

                for (int i = 0; i < propNumStreamLines; i++) {
                    // Sample from distribution and find point position
                    int list_pos = distribution(randGenerator);
                    int x = list_pos / nVertPerDim[1];
                    int y = list_pos % nVertPerDim[1];
                    vec2 startPoint = vectorField.getPositionAtVertex({ x, y });
                    // Randomly sample a point close to the found distribution point
                    startPoint[0] = StreamlineIntegrator::randomValue(startPoint[0] - cellSize[0]/2, startPoint[0] + cellSize[0]/2);
                    startPoint[1] = StreamlineIntegrator::randomValue(startPoint[1] - cellSize[1]/2, startPoint[1] + cellSize[1]/2);
                    startPoint = vectorField.clampPositionToBBox(startPoint);
                    startPoints.push_back(startPoint);
                }

            }
           else { // RANDOM POINTS
                for (int i = 0; i < propNumStreamLines; i++) {
                    vec2 startPoint;
                    startPoint[0] = StreamlineIntegrator::randomValue(BBoxMin_[0], BBoxMax_[0]);
                    startPoint[1] = StreamlineIntegrator::randomValue(BBoxMin_[1], BBoxMax_[1]);
                    startPoints.push_back(startPoint);
                }
           }

           // Draw N stream Lines
           vec2 startPoint;
           while (!startPoints.empty()) {
               startPoint = startPoints.front();
               startPoints.pop_front();
               if (propDisplayPoints.get() != 0)
                   Integrator::drawPoint(startPoint, vec4(0, 0, 0, 1), indexBufferPoints.get(), vertices);
               auto indexBufferEuler = mesh->addIndexBuffer(DrawType::Lines, ConnectivityType::Strip);
               auto indexBufferRK = mesh->addIndexBuffer(DrawType::Lines, ConnectivityType::Strip);
               if (propRk4.get()) {
                   stepsRk4 = Integrator::StreamLines(true, vectorField, startPoint,
                       propFwdDir, propStepSize, propNormalized, false, propNumSteps,
                       false, propArc, true, false, false, propVelocity,
                       propDisplayPoints.get() != 0, indexBufferPoints.get(), indexBufferRK.get(), vertices);
               }
               if (propEuler.get()) {
                   stepsEuler = Integrator::StreamLines(false, vectorField, startPoint,
                       propFwdDir, propStepSize, propNormalized, false, propNumSteps,
                       false, propArc, true, false, false, propVelocity,
                       propDisplayPoints.get() != 0, indexBufferPoints.get(), indexBufferEuler.get(), vertices);
               }
           }
        }
        else if (propSeedUniform.get() == true) { // UNIFORM GRID
            double x_size = (BBoxMax_[0] - BBoxMin_[0]) / (propNumXgrid+1);
            double y_size = (BBoxMax_[1] - BBoxMin_[1]) / (propNumYgrid+1);
            vec2 startPoint;
            for (int i = 0; i < propNumXgrid; i++) {
                for (int j = 0; j < propNumYgrid; j++) {
                    startPoint[0] = BBoxMin_[0] + ((i+1)* x_size);
                    startPoint[1] = BBoxMin_[1] + ((j+1)* y_size);
                    if (propDisplayPoints.get() != 0)
                        Integrator::drawPoint(startPoint, vec4(0, 0, 0, 1), indexBufferPoints.get(), vertices);
                    auto indexBufferEuler = mesh->addIndexBuffer(DrawType::Lines, ConnectivityType::Strip);
                    auto indexBufferRK = mesh->addIndexBuffer(DrawType::Lines, ConnectivityType::Strip);
                    if (propRk4.get()) {
                        stepsRk4 = Integrator::StreamLines(true, vectorField, startPoint,
                            propFwdDir, propStepSize, propNormalized, false, propNumSteps,
                            false, propArc, true, false, false, propVelocity,
                            propDisplayPoints.get() != 0, indexBufferPoints.get(), indexBufferRK.get(), vertices);
                    }
                    if (propEuler.get()) {
                        stepsEuler = Integrator::StreamLines(false, vectorField, startPoint,
                            propFwdDir, propStepSize, propNormalized, false, propNumSteps,
                            false, propArc, true, false, false, propVelocity,
                            propDisplayPoints.get() != 0, indexBufferPoints.get(), indexBufferEuler.get(), vertices);
                    }
                }         
                
            }

        }

    }

    mesh->addVertices(vertices);
    meshOut.setData(mesh);
}


float StreamlineIntegrator::randomValue(const float min, const float max) const {
    return min + uniformReal(randGenerator) * (max - min);
}



}  // namespace inviwo
