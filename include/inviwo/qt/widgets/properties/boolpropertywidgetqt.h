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

#ifndef IVW_BOOLPROPERTYWIDGETQT_H
#define IVW_BOOLPROPERTYWIDGETQT_H

#include <inviwo/qt/widgets/inviwoqtwidgetsdefine.h>
#include <inviwo/qt/widgets/editablelabelqt.h>
#include <QCheckBox>

#include <inviwo/qt/widgets/properties/propertywidgetqt.h>

#include <inviwo/core/properties/boolproperty.h>

namespace inviwo {

class IVW_QTWIDGETS_API BoolPropertyWidgetQt : public PropertyWidgetQt {

    Q_OBJECT

public:
    BoolPropertyWidgetQt(BoolProperty* property);

    void updateFromProperty();

private:
    BoolProperty* property_;
    QCheckBox* checkBox_;
    EditableLabelQt* label_;
    void generateWidget();

public slots:
    void setPropertyValue();
    void setPropertyDisplayName();
};

} // namespace

#endif // IVW_BOOLPROPERTYWIDGETQT_H