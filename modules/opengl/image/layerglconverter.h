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
 * Primary author : Daniel J�nsson
 *
 **********************************************************************/

#ifndef IVW_LAYERGLCONVERTER_H
#define IVW_LAYERGLCONVERTER_H

#include <modules/opengl/openglmoduledefine.h>
#include <modules/opengl/image/layergl.h>
#include <inviwo/core/datastructures/image/layerramconverter.h>

namespace inviwo {

class IVW_MODULE_OPENGL_API LayerRAM2GLConverter : public RepresentationConverterType<LayerGL> {

public:
    LayerRAM2GLConverter();
    virtual ~LayerRAM2GLConverter();

    inline bool canConvertFrom(const DataRepresentation* source) const {
        return dynamic_cast<const LayerRAM*>(source) != NULL;
    }

    DataRepresentation* createFrom(const DataRepresentation* source);
    void update(const DataRepresentation* source, DataRepresentation* destination);
};

class IVW_MODULE_OPENGL_API LayerGL2RAMConverter : public RepresentationConverterType<LayerRAM> {

public:
    LayerGL2RAMConverter();
    virtual ~LayerGL2RAMConverter();

    inline bool canConvertFrom(const DataRepresentation* source) const {
        return dynamic_cast<const LayerGL*>(source) != NULL;
    }

    DataRepresentation* createFrom(const DataRepresentation* source);
    void update(const DataRepresentation* source, DataRepresentation* destination);
};

class IVW_MODULE_OPENGL_API LayerDisk2GLConverter : public RepresentationConverterPackage<LayerGL> {

public:
    LayerDisk2GLConverter() : RepresentationConverterPackage<LayerGL>(){
        addConverter(new LayerDisk2RAMConverter());
        addConverter(new LayerRAM2GLConverter());
    };
    virtual ~LayerDisk2GLConverter() {};
};

} // namespace

#endif // IVW_LAYERGLCONVERTER_H
