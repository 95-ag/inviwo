/*********************************************************************************
 *
 * Inviwo - Interactive Visualization Workshop
 * Version 0.6b
 *
 * Copyright (c) 2012-2014 Inviwo Foundation
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
 * Main file authors: Erik Sund�n, Peter Steneteg
 *
 *********************************************************************************/

#include <inviwo/core/datastructures/image/imagedisk.h>
#include <inviwo/core/datastructures/image/layerdisk.h>

namespace inviwo {

ImageDisk::ImageDisk()
    : ImageRepresentation() {}

ImageDisk::ImageDisk(const ImageDisk& rhs)
    : ImageRepresentation(rhs) {
}
ImageDisk& ImageDisk::operator=(const ImageDisk& that) {
    if (this != &that)
        ImageRepresentation::operator=(that);

    return *this;
}
ImageDisk* ImageDisk::clone() const {
    return new ImageDisk(*this);
}

ImageDisk::~ImageDisk() {
}

void ImageDisk::initialize() {
}

void ImageDisk::deinitialize() {}


std::string ImageDisk::getClassName() const {
    return "ImageDisk";
}

bool ImageDisk::copyAndResizeRepresentation(DataRepresentation*) const {
    return false;
}

void ImageDisk::update(bool editable) {
    Image *owner = this->getOwner();
    if (editable) {
        for (size_t i=0; i<owner->getNumberOfColorLayers(); ++i)
            owner->getColorLayer(i)->getEditableRepresentation<LayerDisk>();

        Layer* depthLayer = owner->getDepthLayer();

        if (depthLayer)
            depthLayer->getEditableRepresentation<LayerDisk>();

        Layer* pickingLayer = owner->getPickingLayer();

        if (pickingLayer)
            pickingLayer->getEditableRepresentation<LayerDisk>();
    }
    else {
        for (size_t i=0; i<owner->getNumberOfColorLayers(); ++i)
            owner->getColorLayer(i)->getRepresentation<LayerDisk>();

        Layer* depthLayer = owner->getDepthLayer();

        if (depthLayer)
            depthLayer->getRepresentation<LayerDisk>();

        Layer* pickingLayer = owner->getPickingLayer();

        if (pickingLayer)
            pickingLayer->getRepresentation<LayerDisk>();
    }
}

} // namespace
