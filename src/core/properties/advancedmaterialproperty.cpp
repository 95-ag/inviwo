/*********************************************************************************
*
* Inviwo - Interactive Visualization Workshop
* Version 0.6b
*
* Copyright (c) 2012-2014 Inviwo Foundation
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
* Main file authors: Daniel J�nsson
*
*********************************************************************************/

#include <inviwo/core/properties/advancedmaterialproperty.h>

namespace inviwo {

PropertyClassIdentifier(AdvancedMaterialProperty, "org.inviwo.AdvancedMaterialProperty");

AdvancedMaterialProperty::AdvancedMaterialProperty(
    std::string identifier, std::string displayName,
    PropertyOwner::InvalidationLevel invalidationLevel, PropertySemantics semantics)
    : CompositeProperty(identifier, displayName, invalidationLevel,
                        semantics)
    , phaseFunctionProp("phaseFunction", "Phase function", PropertyOwner::VALID)
    , indexOfRefractionProp("IOR", "Index of refraction", 1.f, 1.f, 20.f)
    , roughnessProp("roughness", "Roughness", 0.1f, 0.01f, 1.f)
    , specularColorProp("specularColor", "Specular color", vec4(1.f))
    , anisotropyProp("anisotropy", "Anisotropy (g)", 0.f, -1.f, 1.f) {

    phaseFunctionProp.addOption("isotropic", "Isotropic");
    phaseFunctionProp.addOption("HenyeyGreenstein", "Henyey-Greenstein");
    phaseFunctionProp.addOption("Schlick", "Schlick");
    phaseFunctionProp.addOption("Ashikhmin", "Ashikhmin");
    phaseFunctionProp.addOption("BlinnPhong", "Blinn-Phong");
    phaseFunctionProp.addOption("Ward", "Ward");
    phaseFunctionProp.addOption("CookTorrance", "Cook-Torrance");
    phaseFunctionProp.addOption("ABCMicrofacet", "ABC microfacet");
    phaseFunctionProp.addOption("mix", "Mix");
    phaseFunctionProp.setCurrentStateAsDefault();
    phaseFunctionProp.onChange(this, &AdvancedMaterialProperty::phaseFunctionChanged);

    specularColorProp.setSemantics(PropertySemantics::Color);

    addProperty(phaseFunctionProp);
    addProperty(indexOfRefractionProp);
    addProperty(roughnessProp);
    addProperty(specularColorProp);
    addProperty(anisotropyProp);
}

AdvancedMaterialProperty::~AdvancedMaterialProperty() {}

vec4 AdvancedMaterialProperty::getCombinedMaterialParameters() const {
    ShadingFunctionEnum::Enum shadingType = getPhaseFunctionEnum();
    vec4 scaledMaterial(getIndexOfRefractionTerm(1.f, indexOfRefractionProp.get()),
                        roughnessProp.get(), anisotropyProp.get(), 0.f);

    if (shadingType == ShadingFunctionEnum::HENYEY_GREENSTEIN) {
        // scaledMaterial.z = anisotropyProp_.get();
    } else if (shadingType == ShadingFunctionEnum::SCHLICK) {
        scaledMaterial.z = HenyehGreensteinToSchlick(anisotropyProp.get());
    } else if (shadingType == ShadingFunctionEnum::BLINN_PHONG) {
        scaledMaterial.y = roughnessToShininess(roughnessProp.get());
    } else if (shadingType == ShadingFunctionEnum::WARD) {
        scaledMaterial.y = WardParameterMapping(roughnessProp.get());
    } else if (shadingType == ShadingFunctionEnum::COOK_TORRANCE) {
        // scaledMaterial.y = roughnessToShininess(material.y);
        scaledMaterial.y = WardParameterMapping(roughnessProp.get());
    } else if (shadingType == ShadingFunctionEnum::ABC_MICROFACET) {
        // TODO: Use better parameter naming for this
        vec2 bcParameter = BCParameterMapping(vec2(roughnessProp.get(), anisotropyProp.get()));
        scaledMaterial.y = bcParameter.x;
        scaledMaterial.z = bcParameter.y;
    } else if (shadingType == ShadingFunctionEnum::ASHIKHMIN) {
        scaledMaterial.y = roughnessToShininess(roughnessProp.get());
    } else if (shadingType == ShadingFunctionEnum::MIX) {
        scaledMaterial.y = roughnessToShininess(roughnessProp.get());
    }

    return scaledMaterial;
}

void AdvancedMaterialProperty::phaseFunctionChanged() {
    ShadingFunctionEnum::Enum shadingType = getPhaseFunctionEnum();
    indexOfRefractionProp.setVisible(false);
    roughnessProp.setVisible(false);
    specularColorProp.setVisible(false);
    anisotropyProp.setVisible(false);

    if (shadingType == ShadingFunctionEnum::HENYEY_GREENSTEIN) {
        anisotropyProp.setVisible(true);
    } else if (shadingType == ShadingFunctionEnum::SCHLICK) {
        anisotropyProp.setVisible(true);
    } else if (shadingType == ShadingFunctionEnum::BLINN_PHONG) {
        indexOfRefractionProp.setVisible(true);
        roughnessProp.setVisible(true);
        specularColorProp.setVisible(true);
    } else if (shadingType == ShadingFunctionEnum::WARD) {
        indexOfRefractionProp.setVisible(true);
        roughnessProp.setVisible(true);
        specularColorProp.setVisible(true);
    } else if (shadingType == ShadingFunctionEnum::COOK_TORRANCE) {
        indexOfRefractionProp.setVisible(true);
        roughnessProp.setVisible(true);
        specularColorProp.setVisible(true);
    } else if (shadingType == ShadingFunctionEnum::ABC_MICROFACET) {
        // TODO: Use better parameters
        anisotropyProp.setVisible(true);
        roughnessProp.setVisible(true);
        specularColorProp.setVisible(true);
    } else if (shadingType == ShadingFunctionEnum::ASHIKHMIN) {
        indexOfRefractionProp.setVisible(true);
        roughnessProp.setVisible(true);
        specularColorProp.setVisible(true);
    } else if (shadingType == ShadingFunctionEnum::MIX) {
        indexOfRefractionProp.setVisible(true);
        roughnessProp.setVisible(true);
        specularColorProp.setVisible(true);
    }

    propertyModified();
}

ShadingFunctionEnum::Enum AdvancedMaterialProperty::getPhaseFunctionEnum() const {
    const std::string& shadingFunction = phaseFunctionProp.get();

    if (shadingFunction == "HenyeyGreenstein") {
        return ShadingFunctionEnum::HENYEY_GREENSTEIN;
    } else if (shadingFunction == "Schlick") {
        return ShadingFunctionEnum::SCHLICK;
    } else if (shadingFunction == "Ashikhmin") {
        return ShadingFunctionEnum::ASHIKHMIN;
    } else if (shadingFunction == "mix") {
        return ShadingFunctionEnum::MIX;
    } else if (shadingFunction == "BlinnPhong") {
        return ShadingFunctionEnum::BLINN_PHONG;
    } else if (shadingFunction == "Ward") {
        return ShadingFunctionEnum::WARD;
    } else if (shadingFunction == "ABCMicrofacet") {
        return ShadingFunctionEnum::ABC_MICROFACET;
    } else if (shadingFunction == "CookTorrance") {
        return ShadingFunctionEnum::COOK_TORRANCE;
    } else {
        return ShadingFunctionEnum::ISOTROPIC;
    }
}

void AdvancedMaterialProperty::deserialize(IvwDeserializer& d) {
    CompositeProperty::deserialize(d);
    phaseFunctionChanged();
}

float AdvancedMaterialProperty::HenyehGreensteinToSchlick(float g) const {
    // Should work for intermediate values of k in Schlick
    return 1.55f * g - 0.55f * g * g * g;
    // return 1.5f*g-(1.f-1.5f)*g*g*g;
}

float AdvancedMaterialProperty::roughnessToShininess(float m) const {
    return 2.f / (std::max(1e-4f, m * m)) - 2.f;
}

float AdvancedMaterialProperty::WardParameterMapping(float x) const { 
    return std::max(x, 1e-4f); 
}

vec2 AdvancedMaterialProperty::BCParameterMapping(const vec2& bc) const {
    return bc * vec2(13000.f, 1.f);
}

float AdvancedMaterialProperty::getIndexOfRefractionTerm(float n1, float n2) const {
    // Compute index of refraction (IOR) f0 = ((n1-n2)/(n1+n2))^2
    return std::pow((n1 - n2) / (n1 + n2), 2.f);
}

}  // namespace
