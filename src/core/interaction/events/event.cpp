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
 * Primary author : Timo Ropinski
 *
 **********************************************************************/

#include <inviwo/core/interaction/events/event.h>

namespace inviwo {

Event::Event() {}
Event::~Event() {}

void Event::serialize(IvwSerializer& s) const {}
void Event::deserialize(IvwDeserializer& d) {}

} // namespace