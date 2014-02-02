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

#ifndef IVW_FILEDIRECTORY_H
#define IVW_FILEDIRECTORY_H

#include <inviwo/core/common/inviwocoredefine.h>
#include <inviwo/core/common/inviwo.h>

namespace inviwo {

class IVW_CORE_API URLParser {
public:
    URLParser() {}

    static std::string addBasePath(const std::string url);
    static bool fileExists(std::string fileName);
    static std::string getFileDirectory(const std::string url);
    static std::string getFileNameWithExtension(const std::string url);
    static std::string getFileNameWithoutExtension(const std::string url);
    static std::string getFileExtension(const std::string url);
    static std::string replaceFileExtension(const std::string url, const std::string newFileExtension);
    static std::string getRelativePath(const std::string& basePath, const std::string absolutePath);
};

} // namespace

#endif // IVW_FILEDIRECTORY_H
