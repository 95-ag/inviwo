/*********************************************************************************
 *
 * Inviwo - Interactive Visualization Workshop
 * Version 0.6b
 *
 * Copyright (c) 2013-2014 Inviwo Foundation
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice, this
 * list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimer in the documentation
 * and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
 * ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * Main file author: Erik Sund�n
 *
 *********************************************************************************/

#ifndef IVW_DATAGROUPREPRESENTATION_H
#define IVW_DATAGROUPREPRESENTATION_H

#include <inviwo/core/common/inviwocoredefine.h>
#include <inviwo/core/common/inviwo.h>
#include <inviwo/core/datastructures/data.h>
#include <inviwo/core/datastructures/datarepresentation.h>

namespace inviwo {

/** \brief The base class for all DataGroupRepresentation objects.
 *
 *  It has reference to zero or many DataRepresentation objects, but never owns them,
 *  they are always owned by the Data object.
 *
 *  Differences between DataGroupRepresentation and DataRepresentation:
 *    - DataGroupRepresentation does not own DataRepresentation, does should never delete them.
 *    - DataGroupRepresentation becomes invalid when a child DataRepresentation is invalid.
 */

class DataGroup;

class IVW_CORE_API DataGroupRepresentation : public DataRepresentation {

    friend class DataGroup;

public:
    DataGroupRepresentation();
    DataGroupRepresentation(const DataGroupRepresentation& rhs);
    DataGroupRepresentation& operator=(const DataGroupRepresentation& that);
    virtual DataGroupRepresentation* clone() const = 0;
    virtual ~DataGroupRepresentation();

    virtual std::string getClassName() const;
    virtual void performOperation(DataOperation*) const = 0;

    void setAsInvalid();
    bool isValid();

    virtual void initialize() = 0;
    virtual void deinitialize() = 0;

protected:
    //Update representations_ with DataRepresentation from each Data and DataGroup object
    virtual void update(bool) = 0;

    void setAsValid();

private:
    bool valid_;
};

} // namespace

#endif // IVW_DATAGROUPREPRESENTATION_H
