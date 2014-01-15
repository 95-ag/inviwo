/**********************************************************************
 * Copyright (C) 2013 Scientific Visualization Group - Link�ping University
 * All Rights Reserved.
 * 
 * Unauthorized copying of this file, via any medium is strictly prohibited
 * Proprietary and confidential
 * No part of this software may be reproduced or transmitted in any
 * form or by any means including photocopying or recording without
 * written permission of the copyright owner.
 *
 * Primary author : Peter Steneteg
 *
 **********************************************************************/

#include <inviwo/core/datastructures/buffer/bufferrepresentation.h>

namespace inviwo {

BufferRepresentation::BufferRepresentation(size_t size, const DataFormatBase* format, BufferType type, BufferUsage usage)
    : DataRepresentation(format), size_(size), type_(type), usage_(usage){
}

BufferRepresentation::BufferRepresentation(const BufferRepresentation& rhs) 
    : DataRepresentation(rhs), size_(rhs.size_), type_(rhs.type_), usage_(rhs.usage_){
}

BufferRepresentation& BufferRepresentation::operator=(const BufferRepresentation& that) {
    if(this != &that) {
        DataRepresentation::operator=(that);
        size_ = that.size_;
        type_ = that.type_;
        usage_ = that.usage_;
    }
    return *this;
}


} //namespace