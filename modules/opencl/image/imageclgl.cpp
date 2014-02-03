 /*********************************************************************************
 *
 * Inviwo - Interactive Visualization Workshop
 * Version 0.6b
 *
 * Copyright (c) 2014 Inviwo Foundation
 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met: 
 * 
 * 1. Redistributions of source code must retain the above copyright notice, this
 * list of conditions and the following disclaimer. 
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimer in the documentation
 * and/or other materials provided with the distribution. 
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
 * ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * 
 * Main file author: Erik Sund�n
 *
 *********************************************************************************/

#include <modules/opencl/image/imageclgl.h>
#include <modules/opencl/image/layerclgl.h>

namespace inviwo {

ImageCLGL::ImageCLGL()
    : ImageRepresentation(), layerCLGL_(NULL)
{}

ImageCLGL::ImageCLGL(const ImageCLGL& rhs )
    : ImageRepresentation(rhs)
{}

ImageCLGL::~ImageCLGL() { 
}

ImageCLGL* ImageCLGL::clone() const {
    return new ImageCLGL(*this);
}

void ImageCLGL::initialize() {
}

void ImageCLGL::deinitialize() {
}

LayerCLGL* ImageCLGL::getLayerCLGL(){
    return layerCLGL_;
}

const LayerCLGL* ImageCLGL::getLayerCLGL() const {
    return layerCLGL_;
}

bool ImageCLGL::copyAndResizeRepresentation(DataRepresentation* targetRep) const {
    ImageCLGL* targetCLGL = dynamic_cast<ImageCLGL*>(targetRep);

    if (!targetCLGL) return false;
	
	return layerCLGL_->copyAndResizeLayer(targetCLGL->getLayerCLGL());
}

void ImageCLGL::update(bool editable) {
    //TODO: Convert more then just first color layer
    layerCLGL_ = NULL;
    if(editable){
        layerCLGL_ = owner_->getColorLayer()->getEditableRepresentation<LayerCLGL>();
    }
    else{
        layerCLGL_ = const_cast<LayerCLGL*>(owner_->getColorLayer()->getRepresentation<LayerCLGL>());
    }
}

} // namespace

namespace cl {

template <>
cl_int Kernel::setArg(cl_uint index, const inviwo::ImageCLGL& value) {
    return setArg(index, value.getLayerCLGL()->get());
}

} // namespace cl