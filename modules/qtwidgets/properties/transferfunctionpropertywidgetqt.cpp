/*********************************************************************************
 *
 * Inviwo - Interactive Visualization Workshop
 *
 * Copyright (c) 2013-2018 Inviwo Foundation
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
 *********************************************************************************/

#include <modules/qtwidgets/properties/transferfunctionpropertywidgetqt.h>
#include <modules/qtwidgets/properties/collapsiblegroupboxwidgetqt.h>
#include <modules/qtwidgets/editablelabelqt.h>
#include <modules/qtwidgets/properties/transferfunctionpropertydialog.h>
#include <modules/qtwidgets/inviwoqtutils.h>

#include <inviwo/core/datastructures/tfprimitive.h>

#include <warn/push>
#include <warn/ignore/all>
#include <QHBoxLayout>
#include <QWidget>
#include <warn/pop>

namespace inviwo {

TransferFunctionPropertyWidgetQt::TransferFunctionPropertyWidgetQt(
    TransferFunctionProperty* property)
    : PropertyWidgetQt(property)
    , label_{new EditableLabelQt(this, property_)}
    , btnOpenTF_{new TFPushButton(static_cast<TransferFunctionProperty*>(property_), this)} {

    QHBoxLayout* hLayout = new QHBoxLayout();
    hLayout->setContentsMargins(0, 0, 0, 0);
    hLayout->setSpacing(7);

    hLayout->addWidget(label_);

    connect(btnOpenTF_, &TFPushButton::clicked, [this]() {
        if (!transferFunctionDialog_) {
            transferFunctionDialog_ = util::make_unique<TransferFunctionPropertyDialog>(
                static_cast<TransferFunctionProperty*>(property_));
            transferFunctionDialog_->setVisible(true);
        } else {
            transferFunctionDialog_->setVisible(!transferFunctionDialog_->isVisible());
        }
    });

    {
        QWidget* widget = new QWidget(this);
        QSizePolicy sliderPol = widget->sizePolicy();
        sliderPol.setHorizontalStretch(3);
        widget->setSizePolicy(sliderPol);
        QGridLayout* vLayout = new QGridLayout();
        widget->setLayout(vLayout);
        vLayout->setContentsMargins(0, 0, 0, 0);
        vLayout->setSpacing(0);

        vLayout->addWidget(btnOpenTF_);
        hLayout->addWidget(widget);
    }

    setLayout(hLayout);
    updateFromProperty();

    QSizePolicy sp = sizePolicy();
    sp.setVerticalPolicy(QSizePolicy::Fixed);
    setSizePolicy(sp);
}

TransferFunctionPropertyWidgetQt::~TransferFunctionPropertyWidgetQt() {
    if (transferFunctionDialog_) transferFunctionDialog_->hide();
}

void TransferFunctionPropertyWidgetQt::updateFromProperty() { btnOpenTF_->updateFromProperty(); }

TransferFunctionPropertyDialog* TransferFunctionPropertyWidgetQt::getEditorWidget() const {
    return transferFunctionDialog_.get();
}

bool TransferFunctionPropertyWidgetQt::hasEditorWidget() const {
    return transferFunctionDialog_ != nullptr;
}

void TransferFunctionPropertyWidgetQt::setReadOnly(bool readonly) { label_->setDisabled(readonly); }

TFPushButton::TFPushButton(TransferFunctionProperty* property, QWidget* parent)
    : IvwPushButton(parent), tfProperty_(property) {}

void TFPushButton::updateFromProperty() {
    const QSize size = this->size() - QSize(2, 2);
    
    setIcon(utilqt::toQPixmap(*tfProperty_, size));
    setIconSize(size);
}

void TFPushButton::resizeEvent(QResizeEvent* event) {
    updateFromProperty();
    IvwPushButton::resizeEvent(event);
}

}  // namespace inviwo
