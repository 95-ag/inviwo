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
 * Primary author : Erik Sund�n
 *
 **********************************************************************/

#include <inviwo/qt/editor/settingswidget.h>
#include <inviwo/core/common/inviwoapplication.h>
#include <inviwo/qt/widgets/properties/propertywidgetfactoryqt.h>
#include <QLayout>
#include <QFrame>
#include <QSettings>

namespace inviwo {

SettingsWidget::SettingsWidget(QString title, QWidget* parent) : InviwoDockWidget(title, parent) {
    generateWidget();
}

SettingsWidget::SettingsWidget(QWidget* parent) : InviwoDockWidget(tr("Settings"), parent) {
    generateWidget();
}

void SettingsWidget::generateWidget() {
    setObjectName("SettingsWidget");
    setAllowedAreas(Qt::LeftDockWidgetArea | Qt::RightDockWidgetArea);

    //QFrame* frame = new QFrame();
    //vLayout_ = new QVBoxLayout(frame);
    vLayout_ = new QVBoxLayout();
	vLayout_->setSpacing(0);
    tabWidget_ = new QTabWidget();
    QWidget* tab1 = new QWidget();
    QWidget* tab2 = new QWidget();
    QWidget* tab3 = new QWidget();
    tab1->setLayout(vLayout_);
    tabWidget_->addTab(tab1,"Mainsettings");
    tabWidget_->addTab(tab2,"tab2");
    tabWidget_->addTab(tab3,"tab3");
    //frame->setLayout(vLayout_);
    setWidget(tabWidget_);
    //setWidget(frame);
}

SettingsWidget::~SettingsWidget() {}

//Load settings from QSettings
void SettingsWidget::loadSettings() {
    Settings* mainsettings = InviwoApplication::getRef().getSettings();
    std::vector<Property*> properties = mainsettings->getProperties();
    QSettings qmainsettings("Inviwo", "Inviwo");
    qmainsettings.beginGroup("mainsettings");
    QStringList keys = qmainsettings.allKeys();
    for (size_t i=0; i<properties.size(); i++) {
        Property* curProperty = properties[i];
        QString name = QString::fromStdString(curProperty->getIdentifier());
        if(keys.contains(name)){
            QVariant qval = qmainsettings.value(name);
            Variant val(std::string(qval.toString().toLocal8Bit().constData()));
            curProperty->setVariant(val);
        }

        PropertyWidgetQt* propertyWidget = PropertyWidgetFactoryQt::getRef().create(curProperty);
        vLayout_->addWidget(propertyWidget);
        curProperty->registerPropertyWidget(propertyWidget);
    }
    qmainsettings.endGroup();
    vLayout_->addStretch(0);
}

//Save application settings to QSettings
void SettingsWidget::saveSettings() {
    Settings* mainsettings = InviwoApplication::getRef().getSettings();
    std::vector<Property*> properties = mainsettings->getProperties();
    QSettings qmainsettings("Inviwo", "Inviwo");
    qmainsettings.beginGroup("mainsettings");
    for (size_t i=0; i<properties.size(); i++) {
        Property* curProperty = properties[i];
        qmainsettings.setValue(QString::fromStdString(curProperty->getIdentifier()), QVariant(QString::fromStdString(curProperty->getVariant().getString())));
    }
    qmainsettings.endGroup();
}

} // namespace