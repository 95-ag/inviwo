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
 * Primary author : Peter Steneteg
 *
 **********************************************************************/

#include <modules/example/examplemodule.h>

#include <modules/example/exampleprocessor.h>

namespace inviwo {

ExampleModule::ExampleModule() : InviwoModule() {
    setIdentifier("Example");
    setXMLFileName("example/examplemodule.xml");

    registerProcessor(ExampleProcessor);
	
}

} // namespace
             