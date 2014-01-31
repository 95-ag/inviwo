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
 * Primary author : Rickard Englund
 *
 **********************************************************************/

#ifndef IVW_PYLISTMEHTODSINVIWO_H
#define IVW_PYLISTMEHTODSINVIWO_H



#include <modules/python/pythonmoduledefine.h>

#include "../pythoninterface/pymethod.h"

namespace inviwo {
PyObject* py_listProperties(PyObject* /*self*/, PyObject* /*args*/);
PyObject* py_listProcesoors(PyObject* /*self*/, PyObject* /*args*/);


    class IVW_MODULE_PYTHON_API PyListPropertiesMethod : public PyMethod{
    public:
        PyListPropertiesMethod();
        virtual ~PyListPropertiesMethod(){}
        virtual std::string getName()const{return "listProperties";}
        virtual std::string getDesc()const{return "List all properties for a processor";}
        virtual PyCFunction getFunc(){return py_listProperties;}
    private:
        PyParamString processor_;

    };
    class IVW_MODULE_PYTHON_API PyListProcessorsMethod : public PyMethod{
    public:
        virtual std::string getName()const{return "listProcessors";}
        virtual std::string getDesc()const{return "Lists all processors in the current network";}
        virtual PyCFunction getFunc(){return py_listProcesoors;}

    };

} //namespace


#endif // IVW_PYLISTMEHTODSINVIWO_H


