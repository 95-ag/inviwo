/*********************************************************************************
 *
 * Inviwo - Interactive Visualization Workshop
 *
 * Copyright (c) 2012-2015 Inviwo Foundation
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
 *********************************************************************************/

#ifndef IVW_IMAGEPORT_H
#define IVW_IMAGEPORT_H

#include <inviwo/core/common/inviwocoredefine.h>
#include <inviwo/core/ports/datainport.h>
#include <inviwo/core/ports/dataoutport.h>
#include <inviwo/core/datastructures/image/image.h>
#include <inviwo/core/interaction/events/eventhandler.h>
#include <inviwo/core/interaction/events/resizeevent.h>
#include <inviwo/core/util/imagecache.h>

namespace inviwo {

class ImageOutport;

class IVW_CORE_API ImageInport : public DataInport<Image> {
public:
    ImageInport(std::string identifier, bool outportDeterminesSize = false);
    virtual ~ImageInport();

    /**
     * Connects this inport to the outport. Propagates the inport size to the outport if the
     * processor is an end processor (Canvas) or any of the dependent outports of this inport are
     * connected.
     *
     * @note Does not check if the outport is an ImageOutport
     * @param Outport * outport ImageOutport to connect
     */
    virtual void connectTo(Outport* outport) override;
    const Image* getData() const override;
    virtual std::string getContentInfo() const override;

    // Actually returns the requested size... not size of the data.
    uvec2 getDimensions() const;

    /**
     * Handle resize event
     */
    void changeDataDimensions(ResizeEvent* resizeEvent);

    bool isOutportDeterminingSize() const;
    void setOutportDeterminesSize(bool outportDeterminesSize);
    
    void passOnDataToOutport(ImageOutport* outport) const;

private:
    uvec2 requestedDimensions_;
    bool outportDeterminesSize_;
};

class IVW_CORE_API ImageOutport : public DataOutport<Image>, public EventHandler {
    friend class ImageInport;

public:
    ImageOutport(std::string identifier, const DataFormatBase* format = DataVec4UINT8::get(),
                 bool handleResizeEvents = true);

    virtual ~ImageOutport();


    /**
     *	We will not handle resize event if we are not the data owner
     */
    virtual void setData(Image* data, bool ownsData = true) override;
    virtual void setConstData(const Image* data) override;
    const Image* getResizedImageData(uvec2 dimensions);
    
    /**
     * Handle resize event
     */
    void changeDataDimensions(ResizeEvent* resizeEvent);
    uvec2 getDimensions() const;
    /**
     * Set the dimensions of this port without propagating the size
     * through the network. Will resize the image contained within the port.
     */
    void setDimensions(const uvec2& newDimension);

    bool addResizeEventListener(EventListener*);
    bool removeResizeEventListener(EventListener*);

    /**
     * Determine if the image data should be resized during a resize event.
     * We will only resize if we own the data in the port.
     * @param handleResizeEvents True if data should be resized during a resize propagation,
     * otherwise false
     */
    void setHandleResizeEvents(bool handleResizeEvents);
    bool isHandlingResizeEvents() const;

protected:
    virtual void invalidate(InvalidationLevel invalidationLevel) override;

private:
    void updateImageFromInputSource();

    uvec2 dimensions_;
    bool handleResizeEvents_;  // True if data should be resized during a resize propagation,
                               // otherwise false

    ImageCache cache_;
};

}  // namespace

#endif  // IVW_IMAGEPORT_H
