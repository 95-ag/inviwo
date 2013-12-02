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

#include "inviwo/core/io/rawvolumewriter.h"


namespace inviwo {


RawVolumeWriter::RawVolumeWriter()
{}

void RawVolumeWriter::saveRawData(WriterSettings& writerSettings) {
    if (writerSettings.dataFormat_ == DataUINT8::str()) {
        saveData<DataUINT8::type>(writerSettings.rawFileAbsolutePath_, writerSettings.dimensions_, writerSettings.texels_);
    }
    else if (writerSettings.dataFormat_ == DataUINT16::str()) {
        saveData<DataUINT16::type>(writerSettings.rawFileAbsolutePath_, writerSettings.dimensions_, writerSettings.texels_ );
    }
}

} // namespace
