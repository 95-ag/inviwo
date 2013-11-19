#ifndef IVW_VOLUMESLICEGL_H
#define IVW_VOLUMESLICEGL_H

#include <modules/base/basemoduledefine.h>
#include <inviwo/core/common/inviwo.h>
#include <modules/opengl/inviwoopengl.h>
#include <modules/opengl/processorgl.h>
#include <inviwo/core/ports/volumeport.h>
#include <inviwo/core/ports/imageport.h>
#include <inviwo/core/properties/optionproperties.h>
#include <inviwo/core/properties/scalarproperties.h>
#include <inviwo/core/properties/transferfunctionproperty.h>
#include <inviwo/core/datastructures/geometry/geometrytype.h>
#include <modules/opengl/glwrap/shader.h>

namespace inviwo {

class IVW_MODULE_BASE_API VolumeSliceGL : public ProcessorGL {
public:
    VolumeSliceGL();
    ~VolumeSliceGL();
    
    InviwoProcessorInfo();

    void initialize();
    void deinitialize();

protected:
    virtual void process();

    void coordinatePlaneChanged();
    void volumeDimensionChanged();

private:
    VolumeInport inport_;
    ImageOutport outport_;

    TemplateOptionProperty<CoordinatePlane> coordinatePlane_;
    IntProperty sliceNumber_;

    TransferFunctionProperty transferFunction_;

    Shader* shader_;

    uvec3 volumeDimensions_;
};

}

#endif //IVW_VOLUMESLICEGL_H
