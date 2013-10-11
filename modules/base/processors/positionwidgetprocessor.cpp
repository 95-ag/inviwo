#include "positionwidgetprocessor.h"
#include <inviwo/core/interaction/pickingmanager.h>
#include <inviwo/core/datastructures/geometry/basemeshcreator.h>

namespace inviwo {

ProcessorClassName(PositionWidgetProcessor, "PositionWidget"); 
ProcessorCategory(PositionWidgetProcessor, "Geometry Rendering");
ProcessorCodeState(PositionWidgetProcessor, CODE_STATE_EXPERIMENTAL); 

PositionWidgetProcessor::PositionWidgetProcessor()
    : CompositeProcessorGL(),
      inport_("inport"),
      outport_("outport", COLOR_DEPTH_PICKING),
      widgetType_("widgetType", "Widget Type"),
      position_("position", "Position", vec3(0.0f), vec3(-1.f), vec3(1.f)),
      camera_("camera", "Camera", vec3(0.0f, 0.0f, -2.0f), vec3(0.0f, 0.0f, 0.0f), vec3(0.0f, 1.0f, 0.0f))
{
    addPort(inport_);
    addPort(outport_);

    widgetType_.addOption("cube", "Cube");
    widgetType_.addOption("sphere", "Sphere");
    widgetType_.set("sphere");
    addProperty(widgetType_);

    addProperty(position_);
    addProperty(camera_);

    widgetPickingObject_ = PickingManager::instance()->registerPickingCallback(this, &PositionWidgetProcessor::updateWidgetPositionFromPicking);

    vec3 posLLF = vec3(0.5f);
    vec3 posURB = vec3(1.5f);

    widget_ = new Geometry(BaseMeshCreator::rectangularPrism(posLLF, posURB, posLLF, posURB, vec4(1.0f), vec4(1.0f)));
}

PositionWidgetProcessor::~PositionWidgetProcessor() {}

void PositionWidgetProcessor::initialize() {
    CompositeProcessorGL::initialize();

    program_ = new Shader("picking.frag");
    program_->build();
}

void PositionWidgetProcessor::deinitialize() {
    CompositeProcessorGL::deinitialize();

    delete program_;
    program_ = 0;
}

void PositionWidgetProcessor::updateWidgetPositionFromPicking(){
    ivec2 move = widgetPickingObject_->getPickingMove();
    LogInfo("Picking Object with ID : " << widgetPickingObject_->getPickingId() << " moved with 2D vector (" << move.x << "," << move.y << ")");
}

void PositionWidgetProcessor::process() {    
    activateAndClearTarget(outport_);

    program_->activate();
    program_->setUniform("pickingColor_", widgetPickingObject_->getPickingColor());

    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
    glLoadMatrixf(glm::value_ptr(camera_.projectionMatrix()));
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadMatrixf(glm::value_ptr(camera_.viewMatrix()));

    widget_->getRepresentation<GeometryGL>()->render();

    glMatrixMode(GL_PROJECTION);
    glPopMatrix();
    glLoadIdentity();
    glMatrixMode(GL_MODELVIEW);
    glPopMatrix();
    glLoadIdentity();

    program_->deactivate();

    deactivateCurrentTarget();

    compositePortsToOutport(outport_, inport_);
}

} // namespace
