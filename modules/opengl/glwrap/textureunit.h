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

#ifndef IVW_TEXTUREUNIT_H
#define IVW_TEXTUREUNIT_H


#include <modules/opengl/openglmoduledefine.h>
#include <inviwo/core/common/inviwo.h>
#include <modules/opengl/inviwoopengl.h>
#include <inviwo/core/util/capabilities.h>

namespace inviwo {

class IVW_MODULE_OPENGL_API TextureUnit {
public:

    TextureUnit();
    virtual ~TextureUnit();

    GLint getEnum();
    GLint getUnitNumber();
    void requestUnit();
    void activate();

    static void setZeroUnit();

private:
    static bool initialized_;
    static unsigned int numRequestedUnits_;
    static std::vector<bool> requestedUnits_;

    static void initialize();

    bool requested_;
    GLint unitEnum_;
    GLint unitNumber_;
};

} // namespace

#endif // IVW_TEXTUREUNIT_H
