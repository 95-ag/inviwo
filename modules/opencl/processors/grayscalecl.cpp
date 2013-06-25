#include "grayscalecl.h"
#include <modules/opencl/inviwoopencl.h>
#include <modules/opencl/imagecl.h>
#include <modules/opencl/volumecl.h>
#include <modules/opencl/kernelmanager.h>

namespace inviwo {

ProcessorClassName(GrayscaleCL, "GrayscaleCL"); 
ProcessorCategory(GrayscaleCL, "Image Operation");
ProcessorCodeState(GrayscaleCL, CODE_STATE_EXPERIMENTAL);

GrayscaleCL::GrayscaleCL()
    : Processor(),
    inputPort_("color image"),
    outport_("outport"),
    kernel_(NULL)
{
    addPort(inputPort_, "ImagePortGroup1");
    addPort(outport_, "ImagePortGroup1");
}

GrayscaleCL::~GrayscaleCL() {}

void GrayscaleCL::initialize() {
    Processor::initialize();
    try {
        cl::Program* program = KernelManager::getRef().buildProgram(IVW_DIR+"modules/opencl/cl/grayscale.cl");
        kernel_ = KernelManager::getRef().getKernel(program, "grayscaleKernel");
    } catch (cl::Error&) {
        
    }
}

void GrayscaleCL::deinitialize() {
    Processor::deinitialize();
}

void GrayscaleCL::process() {
    if( kernel_ == NULL) {
        return;
    }
    Image* outImage = outport_.getData();
    const ImageCL* colorImageCL = inputPort_.getData()->getRepresentation<ImageCL>();
    outImage->resize(colorImageCL->getDimension());
    uvec2 outportDim = outImage->getDimension();
    ImageCL* outImageCL = outImage->getEditableRepresentation<ImageCL>();

    cl_uint arg = 0;
    kernel_->setArg(arg++, *colorImageCL);
    kernel_->setArg(arg++, *outImageCL);

    OpenCL::getInstance()->getQueue().enqueueNDRangeKernel(*kernel_, cl::NullRange, static_cast<glm::svec2>(outportDim));

}

} // namespace
