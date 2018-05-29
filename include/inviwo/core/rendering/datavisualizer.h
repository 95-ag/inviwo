/*********************************************************************************
 *
 * Inviwo - Interactive Visualization Workshop
 *
 * Copyright (c) 2018 Inviwo Foundation
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

#ifndef IVW_DATAVISUALIZER_H
#define IVW_DATAVISUALIZER_H

#include <inviwo/core/common/inviwocoredefine.h>
#include <inviwo/core/common/inviwo.h>

#include <inviwo/core/ports/inport.h>
#include <inviwo/core/ports/outport.h>
#include <inviwo/core/processors/processor.h>
#include <inviwo/core/util/fileextension.h>

#include <inviwo/core/util/document.h>

namespace inviwo {

/**
 * A Data Visualizer can be used to easily add a visualization of a data type.
 * The DataVisualizerManager keeps a list of registered Data Visualizers.
 * An Inviwo module can create a new Data Visualizer by deriving from DataVisualizer and 
 * implementing all the virtual functions and then call registerDataVisualizer to register it with
*  the DataVisualizerManager
 * A Data Visualizer has two parts the "Source Processor" that can be used to load data from disk.
 * And a "Visualizer Network" that will get data from a Outport like the one in the 
 * "Source Processor" and generate some useful visualization. Both parts a optional, and can 
 * be used independently or together, like:
 *    1) Adding a source processor that can load the requested file
 *    2) Adding a set of processors that will take input from a existing port to create
 *       a visualization
 *    3) Both 1 and 2
 * @see DataVisualizerManager
 */
class IVW_CORE_API DataVisualizer {
public:
    DataVisualizer() = default;
    virtual ~DataVisualizer() = default;

    /**
     * A descriptive name of the visualizer
     */
    virtual std::string getName() const = 0;

    /**
     * A detailed description
     */ 
    virtual Document getDescription() const = 0;

    /**
     * A list of file formats that the visualizer can open
     */
    virtual std::vector<FileExtension> getSupportedFileExtensions() const = 0;

    /**
     * Test if the visualizer can be connected to the given outport
     */ 
    virtual bool isOutportSupported(const Outport* port) const = 0;

    /**
     * Does the visualizer provide a source processor
     */ 
    virtual bool hasSourceProcessor() const = 0;

    /**
     * Does the visualizer provide a visualization network
     */
    virtual bool hasVisualizerNetwork() const = 0;

    /**
     * Adds a source processor with the requested file to the network
     * hasSourceProcessor should return true for this function to be called.
     * @param filename The file to load in the source processor. The extension of the filename
     *        should be in the list of extensions returned by getSupportedFileExtensions.
     * @param network The network to add the processor to.
     * @return a pair or the added source processor and the outport in the source processor
     *         with data from the given file.
     */
    virtual std::pair<Processor*, Outport*> addSourceProcessor(const std::string& filename,
                                                               ProcessorNetwork* network) const = 0;
    /**
     * Adds a set of processors that generate a visualization of the data in the given outport.
     * Is isOutportSupported of outport should return true for this function to work. Other wise
     * call will likely not add any processors.
     * hasVisualizationNetwork should return true for this function to be called.
     * @param outport The port to get data from. It will be connected to the added network.
     * @param network The network to add the Visualizer Network into.
     * @return A list of added processors.
     */
    virtual std::vector<Processor*> addVisualizerNetwork(Outport* outport,
                                                         ProcessorNetwork* network) const = 0;

    /**
     * Adds a source processor with the requested file and a set of processors that generate
     * a visualization of the data.
     * hasSourceProcessor and hasVisualizationNetwork should return true for this function to be
     * called.
     * @param filename The file to load in the source processor. The extension of the filename
     *        should be in the list of extensions returned by getSupportedFileExtensions.
     * @param network The network to add the processor to.
     * @return A list of added processors.
     */
    virtual std::vector<Processor*> addSourceAndVisualizerNetwork(
        const std::string& filename, ProcessorNetwork* network) const = 0;
};

}  // namespace inviwo

#endif  // IVW_DATAVISUALIZER_H
