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

#ifndef IVW_TRANSFERFUNCTION_H
#define IVW_TRANSFERFUNCTION_H

#include <inviwo/core/common/inviwocoredefine.h>
#include <inviwo/core/datastructures/transferfunctiondatapoint.h>
#include <inviwo/core/util/observer.h>
#include <inviwo/core/util/fileextension.h>

namespace inviwo {

class Layer;

template <typename T>
class LayerRAMPrecision;

class IVW_CORE_API TransferFunctionObserver : public Observer {
public:
    virtual void onControlPointAdded(TransferFunctionDataPoint* p);
    virtual void onControlPointRemoved(TransferFunctionDataPoint* p);
    virtual void onControlPointChanged(const TransferFunctionDataPoint* p);
};
class IVW_CORE_API TransferFunctionObservable : public Observable<TransferFunctionObserver> {
protected:
    void notifyControlPointAdded(TransferFunctionDataPoint* p);
    void notifyControlPointRemoved(TransferFunctionDataPoint* p);
    void notifyControlPointChanged(const TransferFunctionDataPoint* p);
};

/**
 * \ingroup datastructures
 * \brief for holding 1D transfer function data.
 *  This class holds 1D transfer function data, currently one parameter in the variable data_.
 */
class IVW_CORE_API TransferFunction : public Serializable,
                                      public TransferFunctionObservable,
                                      public TransferFunctionPointObserver {

public:
    using Point = TransferFunctionDataPoint::Point;

    TransferFunction(size_t textureSize = 1024);
    TransferFunction(const std::vector<Point>& points, size_t textureSize = 1024);
    TransferFunction(const TransferFunction& rhs);
    TransferFunction& operator=(const TransferFunction& rhs);

    virtual ~TransferFunction();

    const Layer* getData() const;
    size_t getNumPoints() const;
    size_t getTextureSize();

    TransferFunctionDataPoint* getPoint(size_t i);
    const TransferFunctionDataPoint* getPoint(size_t i) const;

    /**
     * Add a transfer function point at pos with value color
     *
     * @param pos     position of TF point in range [0,1]
     * @param color   color and opacity, i.e. rgba, of the TF point
     * @throws RangeException if pos is outside [0,1]
     */
    void addPoint(const float& pos, const vec4& color);

    /**
     * Add a transfer function point
     *
     * @param point   TF point to be added
     * @throws RangeException if position of point is outside [0,1]
     */
    void addPoint(const Point& point);

    /**
     * Add a transfer function point at pos.x() where pos.y is used as alpha and the color is
     * interpolated from existing TF points before and after the given position
     *
     * @param pos     pos.x refers to the position of TF point in range [0,1], pos.y will be mapped
     *                to alpha
     * @throws RangeException if pos.x is outside [0,1]
     */
    void addPoint(const vec2& pos);

    /**
     * Add a transfer function points
     *
     * @throws RangeException if any of the given points is outside [0,1]
     */
    void addPoints(const std::vector<Point>& points);

    /**
     * Deprecated. Add a transfer function point at pos.x() with value color, pos.y is not used.
     */
    [[deprecated("was declared deprecated. Use `addPoint(const float& pos, const vec4& color)` instead")]]
    void addPoint(const vec2& pos, const vec4& color);

    void removePoint(TransferFunctionDataPoint* dataPoint);

    void clearPoints();

    float getMaskMin() const;
    void setMaskMin(float maskMin);
    float getMaskMax() const;
    void setMaskMax(float maskMax);

    /**
     * Notify that the layer data (texture) needs to be updated next time it is requested.
     */
    void invalidate();

    virtual void onTransferFunctionPointChange(const TransferFunctionDataPoint* p);

    virtual void serialize(Serializer& s) const;
    virtual void deserialize(Deserializer& d);

    /**
     * Sample the transfer function at position v and return the respective color and
     * opacity (rgba). The range of the transfer function is [0,1].
     *
     * @param v   sampling position, if v is outside the range [0,1] it is clamped to [0,1]
     * @return color and opacity at position v
     */
    vec4 sample(double v) const;
    /**
     * Sample the transfer function at position v and return the respective color and
     * opacity (rgba). The range of the transfer function is [0,1].
     *
     * @param v   sampling position, if v is outside the range [0,1] it is clamped to [0,1]
     * @return color and opacity at position v
     */
    vec4 sample(float v) const;

    friend bool operator==(const TransferFunction& lhs, const TransferFunction& rhs);

    void save(const std::string& filename, const FileExtension& ext = FileExtension()) const;
    void load(const std::string& filename, const FileExtension& ext = FileExtension());

protected:
    void addPoint(std::unique_ptr<TransferFunctionDataPoint> dataPoint);
    void removePoint(std::vector<std::unique_ptr<TransferFunctionDataPoint>>::iterator pos);
    void sort();
    void calcTransferValues() const;

private:
    float maskMin_;
    float maskMax_;
    std::vector<std::unique_ptr<TransferFunctionDataPoint>> points_;
    std::vector<TransferFunctionDataPoint*> sorted_;

    mutable bool invalidData_;
    std::shared_ptr<LayerRAMPrecision<vec4>> dataRepr_;
    std::unique_ptr<Layer> data_;
};

bool operator==(const TransferFunction& lhs, const TransferFunction& rhs);
bool operator!=(const TransferFunction& lhs, const TransferFunction& rhs);

}  // namespace inviwo
#endif  // IVW_TRANSFERFUNCTION_H
