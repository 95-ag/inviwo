/*********************************************************************************
 *
 * Inviwo - Interactive Visualization Workshop
 * Version 0.6b
 *
 * Copyright (c) 2013-2014 Inviwo Foundation
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
 * Main file authors: Daniel J�nsson
 *
 *********************************************************************************/

#ifndef IVW_VOLUMECL_BASE_H
#define IVW_VOLUMECL_BASE_H

#include <inviwo/core/common/inviwo.h>
#include <inviwo/core/datastructures/image/layerrepresentation.h>
#include <inviwo/core/datastructures/volume/volume.h>
#include <modules/opencl/inviwoopencl.h>
#include <modules/opencl/openclmoduledefine.h>

namespace inviwo {
// Parent classes are responsible for creating the appropriate cl::Image type (Image2D, Image3D, Image2DGL/ImageGL and so forth)
// This class enables inviwo to use cl::Image(s) in a generic way (i.e. not caring if it is an Image2D or Image2DGL/ImageGL).
class IVW_MODULE_OPENCL_API VolumeCLBase {

public:
    VolumeCLBase();
    VolumeCLBase(const VolumeCLBase& other);
    virtual ~VolumeCLBase();

    virtual cl::Image& getEditable() { return *clImage_; }
    virtual const cl::Image& get() const { return *const_cast<const cl::Image*>(clImage_); }

    /** 
     * \brief Calculates scaling for 12-bit data dependent on internal OpenCL format.
     * Scaling will be applied using: dataValue * scaling
     * @return vec2 Offset in first component and scaling in second.
     */
    virtual vec2 getVolumeDataOffsetAndScaling(const Volume* volume) const;

protected:
    cl::Image* clImage_;
};

} // namespace

namespace cl {

// Kernel argument specializations for VolumeCLBase type
// (enables calling cl::Queue::setArg with VolumeCLBase)
template <>
IVW_MODULE_OPENCL_API cl_int Kernel::setArg(cl_uint index, const inviwo::VolumeCLBase& value);

} // namespace cl



#endif // IVW_VOLUMECL_BASE_H
