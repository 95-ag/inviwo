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
 * Main file authors: Viktor Axelsson, Timo Ropinski, Erik Sund�n
 *
 *********************************************************************************/

#ifndef IVW_TRANSFERFUNCTIONPROPERTY_H
#define IVW_TRANSFERFUNCTIONPROPERTY_H

#include <inviwo/core/properties/templateproperty.h>

#include <inviwo/core/datastructures/transferfunction.h>
#include <inviwo/core/datastructures/volume/volume.h>
#include <inviwo/core/datastructures/histogram.h>
#include <inviwo/core/ports/volumeport.h>

namespace inviwo {

class IVW_CORE_API TransferFunctionProperty 
    : public TemplateProperty<TransferFunction>
    , public TransferFunctionObserver {

public:
    TransferFunctionProperty(const std::string &identifier,
                             const std::string &displayName,
                             const TransferFunction &value = TransferFunction(),
                             VolumeInport* volumeInport = NULL,
                             PropertyOwner::InvalidationLevel invalidationLevel=PropertyOwner::INVALID_OUTPUT,
                             PropertySemantics semantics=PropertySemantics::Default);
        
    ~TransferFunctionProperty();
    InviwoPropertyInfo();

    const vec2 getMask() const;
    void setMask(float maskMin, float maskMax);

    const vec2& getZoomH() const;
    void setZoomH(float zoomHMin, float zoomHMax);

    const vec2& getZoomV() const;
    void setZoomV(float zoomVMin, float zoomVMax);

    void setShowHistogram(int type);
    int getShowHistogram();
    VolumeInport* getVolumeInport();

    virtual void setCurrentStateAsDefault();
    virtual void resetToDefaultState();

    virtual void serialize(IvwSerializer& s) const;
    virtual void deserialize(IvwDeserializer& d);

    // Override
    virtual void set(const TransferFunction& property);
    virtual void set(const Property *property);
    virtual void onControlPointAdded(TransferFunctionDataPoint* p);
    virtual void onControlPointRemoved(TransferFunctionDataPoint* p);
    virtual void onControlPointChanged(const TransferFunctionDataPoint* p);

private:
    vec2 zoomH_;
    vec2 defaultZoomH_;
    vec2 zoomV_;
    vec2 defaultZoomV_;
    int showHistogram_;
    int defaultShowHistogram_;

    VolumeInport* volumeInport_;
};

} // namespace inviwo

#endif // IVW_TRANSFERFUNCTIONPROPERTY_H