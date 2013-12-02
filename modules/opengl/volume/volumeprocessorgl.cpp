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
 * Primary author : Timo Ropinski
 *
 **********************************************************************/

#include "volumeprocessorgl.h"

namespace inviwo {

VolumeProcessorGL::VolumeProcessorGL()
    : ProcessorGL()
{}
VolumeProcessorGL::~VolumeProcessorGL() {}

void VolumeProcessorGL::initialize() {
    ProcessorGL::initialize();
}

void VolumeProcessorGL::deinitialize() {
    ProcessorGL::deinitialize();
}

} // namespace
