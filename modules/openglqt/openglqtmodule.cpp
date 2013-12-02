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

#include <modules/openglqt/openglqtmodule.h>
#include <modules/openglqt/openglqtcapabilities.h>
#include <inviwo/qt/widgets/inviwoapplicationqt.h>
#include <modules/opengl/canvasprocessorgl.h>
#include <modules/openglqt/processors/canvasprocessorwidgetqt.h>

namespace inviwo {

OpenGLQtModule::OpenGLQtModule() : InviwoModule() {
    setIdentifier("OpenGLQt");
    setXMLFileName("openglqt/openglqtmodule.xml");   

    addProcessorWidgetAndAssociate<CanvasProcessorGL>(new CanvasProcessorWidgetQt());
    addCapabilities(new OpenGLQtCapabilities());
}

} // namespace
