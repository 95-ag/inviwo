/**********************************************************************
 * Copyright (C) 2012 Scientific Visualization Group - Link�ping University
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

#include "modules/glut/glutmodule.h"

#include "modules/glut/canvasglut.h"

namespace inviwo {

GLUTModule::GLUTModule() : InviwoModule() {
    setIdentifier("GLUT");
    setXMLFileName("glut/glutmodule.xml");
}

} // namespace
