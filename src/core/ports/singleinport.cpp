 /*********************************************************************************
 *
 * Inviwo - Interactive Visualization Workshop
 * Version 0.6b
 *
 * Copyright (c) 2013 Inviwo Foundation
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
 * Main file author: Daniel J�nsson
 *
 *********************************************************************************/

#include <inviwo/core/ports/singleinport.h>

#include <inviwo/core/processors/processor.h>

namespace inviwo {

SingleInport::SingleInport(std::string identifier,
                           PropertyOwner::InvalidationLevel invalidationLevel)
    : Inport(identifier)
    , connectedOutport_(NULL)
    , invalidationLevel_(invalidationLevel) {
}

SingleInport::~SingleInport() {}

//Inport should determine if we can connect to the outport
void SingleInport::connectTo(Outport* outport) {
    connectedOutport_ = outport;
    outport->connectTo(this);
    invalidate(invalidationLevel_);
}

void SingleInport::disconnectFrom(Outport* outport) {
    ivwAssert(connectedOutport_==outport, "Ports are not connected.");

    if (outport == connectedOutport_) {
        connectedOutport_ = NULL;
        outport->disconnectFrom(this);
        invalidate(invalidationLevel_);
    }
}

bool SingleInport::isConnected() const {
    return (connectedOutport_!=NULL);
}

bool SingleInport::isConnectedTo(Outport* outport) const {
    return connectedOutport_==outport;
}

void SingleInport::invalidate(PropertyOwner::InvalidationLevel invalidationLevel) {
    invalidationLevel_ = std::max(invalidationLevel_, invalidationLevel);
    //TODO: for port properties Port::invalidate() should be called here
    Port::invalidate(invalidationLevel);
}

} // namespace
