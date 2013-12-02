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
 * Primary author : Sathish Kottravel
 *
 **********************************************************************/

#ifndef IVW_VOLUMEEXPORT_H
#define IVW_VOLUMEEXPORT_H

#include <modules/base/basemoduledefine.h>
#include <inviwo/core/common/inviwo.h>
#include <inviwo/core/processors/processor.h>
#include <inviwo/core/ports/volumeport.h>

namespace inviwo {

class IVW_MODULE_BASE_API VolumeExport : public Processor {
public:
    VolumeExport();
    ~VolumeExport();

    InviwoProcessorInfo();

    virtual void initialize();
    virtual void deinitialize();
    void exportVolume();

protected:
    virtual void process();

private:
    VolumeInport volumePort_;
    FileProperty volumeFileName_;
    ButtonProperty exportVolumeButton_;
};

} // namespace

#endif // IVW_VOLUMEEXPORT_H
