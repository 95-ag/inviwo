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
 * Primary author : Alexander Johansson
 *
 **********************************************************************/

#ifndef IVW_BUTTONPROPERTY_H
#define IVW_BUTTONPROPERTY_H

#include <inviwo/core/util/callback.h>
#include <inviwo/core/properties/property.h>

namespace inviwo {    
/** class ButtonProperty
* \brief  The Button property class provides buttons that you can bind functions to buttons.
* 
* The button property has a widget witch creates a button and register a function to it.
* You can only assign one function to the property.
* To bind a function to a button property use the registerClassMemberFunction.
* Example usage
* myButton_.registerClassMemberFunction(this, &MyButton::doSomethingFunction);
* A button property are normally used by a processor.
* The button property is created and added in the constructor of the processor.
*  
*  
* @see ButtonPropertyWidgetQt
*/
class ButtonProperty : public Property {

public:
	/** 
	 * \brief Constructs a button property
	 *
	 * <DESCRIBE THE METHOD>
	 * 
	 * @param std::string identifier <DESCRIBE ME>
	 * @param std::string displayName <DESCRIBE ME>
	 * @param PropertyOwner::InvalidationLevel invalidationLevel <DESCRIBE ME>
	 * @param PropertySemantics::Type semantics <DESCRIBE ME>
	 * @return  <DESCRIBE ME>
	 */ButtonProperty(std::string identifier, std::string displayName,
                   PropertyOwner::InvalidationLevel invalidationLevel=PropertyOwner::INVALID_OUTPUT, 
                   PropertySemantics::Type semantics = PropertySemantics::Default) : Property(identifier,displayName,invalidationLevel,semantics) {}

	virtual void serialize(IvwSerializer& s) const{
        Property::serialize(s);
    }

    virtual void deserialize(IvwDeserializer& d){
        Property::deserialize(d);
    }

    virtual std::string getClassName()  const { return "ButtonProperty"; }

};

} //namespace

#endif //IVW_BUTTONPROPERTY_H