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
 * Primary author : Daniel J�nsson
 *
 **********************************************************************/

#ifndef IVW_BASECLMODULE_H
#define IVW_BASECLMODULE_H

#include <modules/basecl/baseclmoduledefine.h>
#include <inviwo/core/common/inviwomodule.h>

namespace inviwo {

class IVW_MODULE_BASECL_API BaseCLModule: public InviwoModule {

public:
    BaseCLModule();
	virtual ~BaseCLModule();

};

} // namespace

#endif // IVW_BASECLMODULE_H
