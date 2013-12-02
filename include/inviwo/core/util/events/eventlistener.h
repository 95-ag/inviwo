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

#ifndef IVW_EVENTLISTENER_H
#define IVW_EVENTLISTENER_H

#include "inviwo/core/inviwocoredefine.h"
#include "inviwo/core/inviwo.h"

namespace inviwo {

class IVW_CORE_API EventListener {

public:
    EventListener();
    virtual ~EventListener();
};

} // namespace

#endif // IVW_EVENTLISTENER_H
