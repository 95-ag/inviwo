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
 * Primary author : Viktor Axelsson
 *
 **********************************************************************/

#include <inviwo/core/properties/transferfunctionproperty.h>

namespace inviwo {

TransferFunctionProperty::TransferFunctionProperty(std::string identifier, std::string displayName, TransferFunction value, PropertyOwner::InvalidationLevel invalidationLevel, PropertySemantics::Type semantics )
: TemplateProperty<TransferFunction>(identifier, displayName, value, invalidationLevel, semantics), maskProperty_("mask_", "Mask", 0, 255, 0, 255), zoomProperty_("zoom_", "Zoom", 0, 255, 0, 255)
{}

void TransferFunctionProperty::serialize(IvwSerializer& s) const {
	Property::serialize(s);
	std::stringstream stream;
    s.serialize("size", (int)value_.getNumberOfDataPoints());
    s.serialize("zoom_", zoomProperty_);
    s.serialize("mask_", maskProperty_);
    //s.serialize("maskMin", (int)value_.getMaskMin());
    //s.serialize("maskMax", (int)value_.getMaskMax());

	for (int i = 0; i < static_cast<int>(value_.getNumberOfDataPoints()); i++){
		stream << "pos" << i;
		s.serialize(stream.str(), value_.getPoint(i)->getPos());
		stream.clear();
		stream.str(std::string());

		stream << "rgba" << i;
		s.serialize(stream.str(), value_.getPoint(i)->getRgba());
		stream.clear();
		stream.str(std::string());
	}
}

void TransferFunctionProperty::deserialize(IvwDeserializer& d) {
	Property::deserialize(d);
	int size;
    //float maskMin;
    //float maskMax;
	vec2 pos;
	vec4 rgba;
	std::stringstream stream;

	d.deserialize("size", size);
    d.deserialize("zoom_", zoomProperty_);
    //d.deserialize("maskMin", maskMin);
    //d.deserialize("maskMax", maskMax);

    //value_.setMaskMin(maskMin);
    //value_.setMaskMax(maskMax);

	for (int i = 0; i < size; i++){
		stream << "pos" << i;
		d.deserialize(stream.str(), pos);
		stream.clear();
		stream.str(std::string());

		stream << "rgba" << i;
		d.deserialize(stream.str(), rgba);
		stream.clear();
		stream.str(std::string());

		value_.addPoint(pos, rgba);
	}

    propertyModified();
}

} // namespace
