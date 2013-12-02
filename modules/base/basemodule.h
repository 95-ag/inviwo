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

#ifndef IVW_BASEMODULE_H
#define IVW_BASEMODULE_H

#include <modules/base/basemoduledefine.h>
#include <inviwo/core/common/inviwomodule.h>

namespace inviwo {

class IVW_MODULE_BASE_API BaseModule : public InviwoModule {

public:
    BaseModule();

};

} // namespace

#endif // IVW_BASEMODULE_H
