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
 * Primary author : Erik Sund�n
 *
 **********************************************************************/

#include <inviwo/core/util/settings.h>

namespace inviwo {

Settings::Settings(){}

Settings::~Settings() {
}

void Settings::initialize() {
}

void Settings::deinitialize() {
}

void Settings::invalidate() {
    PropertyOwner::invalidate(PropertyOwner::INVALID_OUTPUT);
}

void Settings::serialize(IvwSerializer& s) const {
    PropertyOwner::serialize(s);
}
void Settings::deserialize(IvwDeserializer& d) {
    PropertyOwner::deserialize(d);  
}

} // namespace
