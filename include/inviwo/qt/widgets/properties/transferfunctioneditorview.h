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

#ifndef IVW_TRANSFERFUNCTIONEDITORVIEW_H
#define IVW_TRANSFERFUNCTIONEDITORVIEW_H

#include <inviwo/core/properties/transferfunctionproperty.h>
#include <inviwo/qt/widgets/inviwoqtwidgetsdefine.h>
#include <inviwo/qt/widgets/properties/transferfunctioneditor.h>
#include <inviwo/qt/widgets/inviwoqtwidgetsdefine.h>
#include <inviwo/core/datastructures/volume/volumeram.h>
#include <QtEvents>
#include <QGraphicsScene>
#include <QGraphicsView>
#include <QThread>

namespace inviwo {
class IVW_QTWIDGETS_API TransferFunctionEditorView : public QGraphicsView, public VoidObserver {

	Q_OBJECT

public:
    TransferFunctionEditorView(TransferFunctionProperty* tfProperty);

    void setMask(float maskMin, float maskMax) { if (maskMax<maskMin) maskMax=maskMin; mask_ = vec2(maskMin, maskMax); }
    virtual void notify();

signals:
    void resized();

public slots:
    void histogramThreadFinished();
    void zoomHorizontally(int zoomHMin, int zoomHMax);
    void zoomVertically(int zoomVMin, int zoomVMax);

protected:
    const NormalizedHistogram* getNormalizedHistogram();

    void updateZoom();
	void resizeEvent(QResizeEvent * event);

	void drawForeground(QPainter *painter, const QRectF &rect);
    void drawBackground(QPainter* painter, const QRectF& rect);

private:
    TransferFunctionProperty* tfProperty_;
    vec2 mask_;
    vec2 zoomH_;
    vec2 zoomV_;
    VolumeInport* volumeInport_;

    bool histogramTheadWorking_;
    QThread* workerThread_;
};

class IVW_QTWIDGETS_API HistogramWorkerQt : public QObject{
    Q_OBJECT
public:
    HistogramWorkerQt(const VolumeRAM* volumeRAM) : volumeRAM_(volumeRAM){}
    ~HistogramWorkerQt(){ volumeRAM_ = NULL; };

public slots:
    void process();

signals:
    void finished();

private:
    const VolumeRAM* volumeRAM_;
};

} // namespace

#endif // IVW_TRANSFERFUNCTIONEDITORVIEW_H