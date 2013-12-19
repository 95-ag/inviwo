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

#ifndef IVW_GEOMETRYREPRESENTATION_H
#define IVW_GEOMETRYREPRESENTATION_H

#include <inviwo/core/datastructures/datagrouprepresentation.h>
#include <inviwo/core/datastructures/geometry/geometrytype.h>

namespace inviwo {

class IVW_CORE_API GeometryRepresentation : public DataGroupRepresentation {

public:
    GeometryRepresentation();
    virtual ~GeometryRepresentation();
    virtual void performOperation(DataOperation*) const;
    virtual DataRepresentation* clone() const = 0;
    virtual std::string getClassName() const;
};

} // namespace

#endif // IVW_GEOMETRYREPRESENTATION_H
