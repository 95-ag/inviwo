/**********************************************************************
 * Copyright (C) 2012-2013 Scientific Visualization Group - Link�ping University
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

#include <inviwo/core/common/inviwocore.h>

//Data Structures
#include <inviwo/core/datastructures/volume/volumeramconverter.h>
#include <inviwo/core/datastructures/image/imageramconverter.h>

//Meta Data
#include <inviwo/core/metadata/metadata.h>
#include <inviwo/core/metadata/processormetadata.h>
#include <inviwo/core/metadata/processorwidgetmetadata.h>

//Utilizes
#include <inviwo/core/util/formatconversion.h>
#include <inviwo/core/util/systemcapabilities.h>
#include <inviwo/core/util/vectoroperations.h>

//Io
#include <inviwo/core/io/datvolumereader.h>
#include <inviwo/core/io/ivfvolumereader.h>
#include <inviwo/core/io/datvolumewriter.h>
#include <inviwo/core/io/ivfvolumewriter.h>

//Others
#include <inviwo/core/processors/canvasprocessor.h>
//Properties
#include <inviwo/core/properties/buttonproperty.h>


namespace inviwo {

InviwoCore::InviwoCore() : InviwoModule() {
    setIdentifier("Core");
    
    // Register Converters
    addRepresentationConverter(new VolumeDisk2RAMConverter());
    addRepresentationConverter(new ImageDisk2RAMConverter());

    // Register MetaData
    #define MetaDataMacroC(n, t, d) addMetaData(new n##MetaData());
	MetaDataMacroC(Bool, bool, false)
	MetaDataMacroC(Int, int, 0)
	MetaDataMacroC(Float, float, 0.0f)
	MetaDataMacroC(String, std::string, "")
	MetaDataMacroC(IVec2, ivec2, ivec2(0))
	MetaDataMacroC(IVec3, ivec3, ivec3(0))
	MetaDataMacroC(IVec4, ivec4, ivec4(0))
	MetaDataMacroC(Vec2, vec2, vec2(0))
	MetaDataMacroC(Vec3, vec3, vec3(0))
	MetaDataMacroC(Vec4, vec4, vec4(0))
	MetaDataMacroC(DVec2, dvec2, dvec2(0))
	MetaDataMacroC(DVec3, dvec3, dvec3(0))
	MetaDataMacroC(DVec4, dvec4, dvec4(0))
	MetaDataMacroC(UVec2, uvec2, uvec2(0))
	MetaDataMacroC(UVec3, uvec3, uvec3(0))
	MetaDataMacroC(UVec4, uvec4, uvec4(0))
	MetaDataMacroC(Mat2, mat2, mat2(0))
	MetaDataMacroC(Mat3, mat3, mat3(0))
	MetaDataMacroC(Mat4, mat4, mat4(0))
	
	//#include <inviwo/core/metadata/metadatadefinefunc.h>

    addMetaData(new VectorMetaData<2,float>());
    addMetaData(new VectorMetaData<3,float>());
    addMetaData(new VectorMetaData<4,float>());

    addMetaData(new VectorMetaData<2,int>());
    addMetaData(new VectorMetaData<3,int>());
    addMetaData(new VectorMetaData<4,int>());

    addMetaData(new VectorMetaData<2,unsigned int>());
    addMetaData(new VectorMetaData<3,unsigned int>());
    addMetaData(new VectorMetaData<4,unsigned int>());

    addMetaData(new MatrixMetaData<2,float>());
    addMetaData(new MatrixMetaData<3,float>());
    addMetaData(new MatrixMetaData<4,float>());

    addMetaData(new PositionMetaData());
    addMetaData(new ProcessorMetaData());
    addMetaData(new ProcessorWidgetMetaData());

    // Register Capabilities
    addCapabilities(new SystemCapabilities());

    // Register Data readers
    addDataReader(new DatVolumeReader());
    addDataReader(new IvfVolumeReader());

    // Register Data writers
    addDataWriter(new DatVolumeWriter());
    addDataWriter(new IvfVolumeWriter());

    allocTest_ = NULL;
}

void InviwoCore::setupModuleSettings(){
    if (getSettings()){

        OptionPropertyInt* viewMode_ = new OptionPropertyInt("viewMode","",0);
        viewMode_->addOption("developerMode","developerMode",0);
        viewMode_->addOption("applicationMode","applicationMode",1);
        getSettings()->addProperty(viewMode_);
        viewMode_->setVisibility(INVISIBLE);
        getSettings()->addProperty(new BoolProperty("txtEditor", "Use system text editor", true));

        getSettings()->addProperty(new BoolProperty("shaderReloading", "Automatically reload shaders", true));

        getSettings()->addProperty(new BoolProperty("enablePortInspectors", "Enable port inspectors", true));

        getSettings()->addProperty(new BoolProperty("enableSound", "Enable sound", true));

        getSettings()->addProperty(new BoolProperty("displayLinks", "Display links", true));

        getSettings()->addProperty(new IntProperty("useRAMPercent", "Max Use Mem %", 50, 1, 100));

        ButtonProperty* btnAllocTest = new ButtonProperty("allocTest", "Perform Allocation Test");
        btnAllocTest->onChange(this, &InviwoCore::allocationTest);
        getSettings()->addProperty(btnAllocTest);

        SystemCapabilities* sysInfo = getTypeFromVector<SystemCapabilities>(getCapabilities());
        if (sysInfo){
            ButtonProperty* btnSysInfo = new ButtonProperty("printSysInfo", "Print System Info");
            btnSysInfo->onChange(sysInfo, &SystemCapabilities::printInfo);
            getSettings()->addProperty(btnSysInfo);  
        }           
    }
}

void InviwoCore::allocationTest(){
    if (getSettings()){
        SystemCapabilities* sysInfo = getTypeFromVector<SystemCapabilities>(getCapabilities());
        if (sysInfo){
            if (allocTest_){
                delete allocTest_;
                LogInfo("Deleted previous test allocation");
            }
            IntProperty* useRAMPercent = dynamic_cast<IntProperty*>(getSettings()->getPropertyByIdentifier("useRAMPercent"));
            uint64_t memBytesAlloc = sysInfo->getAvailableMemory(); //In Bytes
            LogInfo("Maximum Available Memory is " << formatBytesToString(memBytesAlloc));
            memBytesAlloc /= 100; //1% of total available memory
            memBytesAlloc *= useRAMPercent->get(); //?% of total available memory
            try
            {
                allocTest_ = new uint32_t[static_cast<uint32_t>(memBytesAlloc/4)];
                LogInfo("Allocated " << formatBytesToString(memBytesAlloc) << ", which is " << useRAMPercent->get() << "% of available memory");
            }
            catch(std::bad_alloc&)
            {
                LogError("Failed allocation of " << formatBytesToString(memBytesAlloc) << ", which is " << useRAMPercent->get() << "% of available memory");
            }
        }
    }
}

} // namespace
