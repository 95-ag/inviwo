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

#include <inviwo/core/datastructures/diskrepresentation.h>

namespace inviwo {

DiskRepresentation::DiskRepresentation() : sourceFile_(""), reader_(NULL) {}

DiskRepresentation::DiskRepresentation(std::string srcFile) : sourceFile_(srcFile), reader_(NULL) {}

DiskRepresentation::DiskRepresentation(const DiskRepresentation& rhs) 
    : sourceFile_(rhs.sourceFile_)
    , reader_(rhs.reader_->clone()){
}

DiskRepresentation& DiskRepresentation::operator=(const DiskRepresentation& that) {
    if(this != &that) {
        sourceFile_ = that.sourceFile_;
        if(reader_) {
            delete reader_;
            reader_ = NULL;
        }
        reader_ = that.reader_->clone();
    }
    return *this;
}

DiskRepresentation* DiskRepresentation::clone() const {
    return new DiskRepresentation(*this);
}

DiskRepresentation::~DiskRepresentation(){
    if(reader_){
        delete reader_;
        reader_ = NULL;
    }
}

const std::string& DiskRepresentation::getSourceFile() const { 
    return sourceFile_; 
}

bool DiskRepresentation::hasSourceFile() const { 
    return !sourceFile_.empty(); 
}

void DiskRepresentation::setDataReader(DataReader* reader){
    if(reader_){
        delete reader_;
    }
    reader_ = reader;
}

void* DiskRepresentation::readData() const {
    if(reader_)
        return reader_->readData();

    return NULL;
}

void DiskRepresentation::readDataInto(void* dest) const {
    if(reader_)
        reader_->readDataInto(dest);
}

} // namespace
