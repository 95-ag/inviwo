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

#ifndef IVW_PICKINGOBJECT_H
#define IVW_PICKINGOBJECT_H

#include <inviwo/core/common/inviwocoredefine.h>
#include <inviwo/core/common/inviwo.h>
#include <inviwo/core/util/callback.h>

namespace inviwo {

/** \class PickingObject 
*/
class IVW_CORE_API PickingObject {

public:
    PickingObject(size_t, DataVec3UINT8::type);

    virtual ~PickingObject();

    const size_t& getPickingId() const;
    const vec3& getPickingColor() const;
    const DataVec3UINT8::type& getPickingColorAsUINT8() const;
    const vec2& getPickingPosition() const;
    const vec2& getPickingMove() const;
    const float& getPickingDepth() const;

    void picked() const;

    void setPickingMove(vec2);
    void setPickingPosition(vec2);
    void setPickingDepth(float);

    SingleCallBack* getCallbackContainer();

private:
    size_t id_;
    DataVec3UINT8::type colorUINT8_;
    vec3 color_;
    vec2 pos_;
    float depth_;
    vec2 move_;
    SingleCallBack* onPickedCallback_;
};

} // namespace

#endif // IVW_PICKINGOBJECT_H