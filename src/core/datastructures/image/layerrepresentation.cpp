/**********************************************************************
 * Copyright (C) 2012-2013 Scientific Visualization Group - Link�ping University
 * All Rights Reserved.
 * 
 * Unauthorized copying of this file, via any medium is strictly prohibited
 * Proprietary and confidential
 * No part of this software may be reproduced or transmitted in any
 * form or by any means including photocopying or recording without
 * written permission of the copyright owner.
 *
 * Primary author : Erik Sund�n
 *
 **********************************************************************/

#include <inviwo/core/datastructures/image/layerrepresentation.h>

namespace inviwo {

LayerRepresentation::LayerRepresentation(uvec2 dimensions, const DataFormatBase* format)
    : DataRepresentation(format), dimensions_(dimensions){
}

LayerRepresentation::~LayerRepresentation() {}

void LayerRepresentation::resize(uvec2 dimensions){
    dimensions_ = dimensions;
}    


} // namespace
