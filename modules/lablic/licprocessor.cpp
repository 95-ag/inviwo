/*********************************************************************
 *  Author  : Himangshu Saikia
 *  Init    : Monday, October 02, 2017 - 13:31:17
 *
 *  Project : KTH Inviwo Modules
 *
 *  License : Follows the Inviwo BSD license model
 *********************************************************************
 */

#include <inviwo/core/datastructures/volume/volumeram.h>
#include <lablic/licprocessor.h>
#include <labstreamlines/integrator.h>

namespace inviwo {

// The Class Identifier has to be globally unique. Use a reverse DNS naming scheme
const ProcessorInfo LICProcessor::processorInfo_{
    "org.inviwo.LICProcessor",  // Class identifier
    "LICProcessor",             // Display name
    "KTH Labs",                 // Category
    CodeState::Experimental,    // Code state
    Tags::None,                 // Tags
};

const ProcessorInfo LICProcessor::getProcessorInfo() const { return processorInfo_; }

LICProcessor::LICProcessor()
    : Processor()
    , volumeIn_("volIn")
    , noiseTexIn_("noiseTexIn")
    , licOut_("licOut")
    , debugOut_("debugOut")
// TODO: Register additional properties
    , propKernalSize("kernalSize", "Kernal Size", 3, 1, 100, 1)
    , propLICType("lICType", "LIC Type")
    , propContrast("contrast", "Constrast", false)
    , propMean("mean", "Mean", 0.5, 0.01, 1, 0.01)
    , propStdDev("sd", "SD", 0.1, 0.01, 1, 0.01)
    , propColor("color", "Color", false)
    , propIsoTransferFunc("isoTransferFunc", "Colors", &volumeIn_)
{
    // Register ports
    addPort(volumeIn_);
    addPort(noiseTexIn_);
    addPort(licOut_);
    addPort(debugOut_); 

    // Register properties
    // TODO: Register additional properties
    addProperty(propLICType); 
    propLICType.addOption("fastLIC", "Fast LIC", 0);
    propLICType.addOption("lIC", "LIC", 1);
    addProperty(propKernalSize);
    addProperty(propContrast);
    addProperty(propMean);
    addProperty(propStdDev);
    addProperty(propColor);
    addProperty(propIsoTransferFunc);

    // The default transfer function has just two blue points
    propIsoTransferFunc.get().clear();
    propIsoTransferFunc.get().add(0.0f, vec4(0.0f, 0.0f, 1.0f, 1.0f)); //blue
    propIsoTransferFunc.get().add(0.25f, vec4(0.0f, 1.0f, 1.0f, 1.0f)); //cyan
    propIsoTransferFunc.get().add(0.5f, vec4(0.0f, 1.0f, 0.0f, 1.0f)); //green
    propIsoTransferFunc.get().add(0.75f, vec4(1.0f, 1.0f, 0.0f, 1.0f)); //yellow
    propIsoTransferFunc.get().add(1.0f, vec4(1.0f, 0.0f, 0.0f, 1.0f)); //red
    propIsoTransferFunc.setCurrentStateAsDefault();

    util::hide(propStdDev, propMean, propKernalSize);
    propContrast.onChange([this]() {
        if (propContrast.get() == true) {
            util::show(propStdDev, propMean);
        }
        else {
            util::hide(propStdDev, propMean);
        }
    });
    propLICType.onChange([this]() {
        if (propLICType.get() == 0) {
            util::hide(propKernalSize);
        }
        else {
            util::show(propKernalSize);
        }
    });

    util::hide(propIsoTransferFunc);
    propColor.onChange([this]() {
        if (propColor.get() == 0) {
            util::hide(propIsoTransferFunc);
        }
        else {
            util::show(propIsoTransferFunc);
        }
        });
}

void LICProcessor::process() {
    // Get input
    if (!volumeIn_.hasData()) {
        return;
    }

    if (!noiseTexIn_.hasData()) {
        return;
    }

    auto vol = volumeIn_.getData();
    const VectorField2 vectorField = VectorField2::createFieldFromVolume(vol);
    vectorFieldDims_ = vol->getDimensions();
    const dvec2 BBoxMin_ = vectorField.getBBoxMin();
    const dvec2 BBoxMax_ = vectorField.getBBoxMax();
    const ivec2 nVertPerDim = vectorField.getNumVerticesPerDim();
    const dvec2 cellSize = vectorField.getCellSize();
    // Extract the minimum and maximum value from the input data
    //double minVal = glm::length(vectorField.getMinValue());
    //double maxVal = glm::length(vectorField.getMaxValue());
    

    auto tex = noiseTexIn_.getData();
    const RGBAImage texture = RGBAImage::createFromImage(tex);
    texDims_ = tex->getDimensions();

    double value = texture.readPixelGrayScale(size2_t(0, 0));

    LogProcessorInfo(value);

    // Prepare the output, it has the same dimensions as the texture and rgba values in [0,255]
    auto outImage = std::make_shared<Image>(texDims_, DataVec4UInt8::get());
    RGBAImage licImage(outImage);

    // Debug image 
    auto debugOutImage = std::make_shared<Image>(texDims_, DataVec4UInt8::get());
    RGBAImage debugImage(debugOutImage);


    std::vector<std::vector<double>> licTexture(texDims_.x, std::vector<double>(texDims_.y, 0.0));
    std::vector<std::vector<double>> licMag(texDims_.x, std::vector<double>(texDims_.y, 0.0));

    // Hint: Output an image showing which pixels you have visited for debugging
    std::vector<std::vector<int>> visited(texDims_.x, std::vector<int>(texDims_.y, 0));

    // TODO: Implement LIC and FastLIC
    // This code instead sets all pixels to the same gray value

    double pixelSize_x = (double) (BBoxMax_[0]-BBoxMin_[0]) / texDims_.x;
    double pixelSize_y = (double) (BBoxMax_[1]-BBoxMin_[1])/ texDims_.y; 
    //std::cout << "BBoxMax_x: "<< BBoxMax_[0] <<  "BBoxMax_y: "<< BBoxMax_[1] << "BBoxMin_x: "<< BBoxMin_[0] <<  "BBoxMin_y: "<< BBoxMin_[1] << "\n";
    //std::cout << "texDims_.x: "<< texDims_.x <<  "texDims_.y "<< texDims_.y << "\n"; 
    //std::cout << "pixelSize_x: "<< pixelSize_x <<  "pixelSize_y: "<< pixelSize_y << "\n"; 
    //std::cout << "vectorFieldDims_.x: "<< vectorFieldDims_.x <<  "vectorFieldDims_.y: "<< vectorFieldDims_.y << "\n"; 
    //std::cout << "nVertPerDim.x: "<< nVertPerDim[0] <<  "nVertPerDim.y: "<< nVertPerDim[1] << "\n";
    LogProcessorInfo("pixelSize_x: " << pixelSize_x << "pixelSize_y: " << pixelSize_y << "\n");

    float stepSize = std::min(pixelSize_x, pixelSize_y);
    float arc_length = propKernalSize * stepSize;
    vec2 startPoint;
    double minVal = INFINITY, maxVal = -INFINITY;
    // LIC
    if(propLICType == 1){
        for (int j = 0; j < (int)texDims_.y; j++) {
            for (int i = 0; i < (int)texDims_.x; i++) {
                startPoint[0] = BBoxMin_[0] + ((i)* pixelSize_x);
                startPoint[1] = BBoxMin_[1] + ((j)* pixelSize_y);

                int counter = 0;
                double sum = 0.0;
                if(vectorField.isInside(startPoint)){
                    //std::cout << " startPoint[0]: "<<  startPoint[0] <<  " startPoint[1]: "<<  startPoint[1] << "\n"; 
                    std::list<vec2> streamLines = Integrator::StreamLine(vectorField, startPoint, stepSize, true, arc_length);
                    //std::cout << "after streamLines" << "\n"; 
                    dvec2 point; 
                    while(!streamLines.empty()){
                        int new_i, new_j; 
                        point = streamLines.front();
                        //std::cout << "   point[0]: "<<  point[0] <<  " point[1]: "<<  point[1] << "\n"; 
                        new_i = int ((point[0] - BBoxMin_[0]) /pixelSize_x); 
                        new_j = int ((point[1] - BBoxMin_[1]) /pixelSize_y);
                        //std::cout << new_i << " " << new_j << "\n"; 
                        if(new_i >= 0 && new_i < (int)texDims_.x && new_j >= 0 && new_j < (int)texDims_.y){
                            sum += texture.readPixelGrayScale(size2_t(new_i, new_j));
                            counter++;
                            visited[new_i][new_j] = 1; 
                            //std::cout << new_i << " " << new_j << "\n"; 
                        }
                        streamLines.pop_front(); 
                    }
                    licTexture[i][j] = sum/counter;
                }
                licMag[i][j] = glm::length(vectorField.interpolate(startPoint));
                minVal = std::min(minVal, licMag[i][j]);
                maxVal = std::max(maxVal, licMag[i][j]);
            }
        }
    }
    // Fast LIC 
    else{
        for (int j = 0; j < (int)texDims_.y; j++) {
            for (int i = 0; i < (int)texDims_.x; i++) {
                startPoint[0] = BBoxMin_[0] + ((i)*pixelSize_x);
                startPoint[1] = BBoxMin_[1] + ((j)*pixelSize_y);
                if(visited[i][j] == 0){
                    int counter = 0;  
                    double sum = 0.0; 
                    if(vectorField.isInside(startPoint)){
                        //std::cout << " startPoint[0]: "<<  startPoint[0] <<  " startPoint[1]: "<<  startPoint[1] << "\n"; 
                        std::list<vec2> streamLines = Integrator::StreamLine(vectorField, startPoint, stepSize, false, arc_length);
                        //std::cout << "after streamLines" << "\n"; 
                        dvec2 point; 
                        std::list<vec2> streamLinesPixels; 
                        while(!streamLines.empty()){
                            int new_i, new_j; 
                            point = streamLines.front();
                            //std::cout << "   point[0]: "<<  point[0] <<  " point[1]: "<<  point[1] << "\n"; 
                            new_i = int ((point[0] - BBoxMin_[0]) /pixelSize_x); 
                            new_j = int ((point[1] - BBoxMin_[1]) /pixelSize_y);
                            //std::cout << new_i << " " << new_j << "\n"; 
                            if(new_i >= 0 && new_i < (int)texDims_.x && new_j >= 0 && new_j < (int)texDims_.y){
                                sum += texture.readPixelGrayScale(size2_t(new_i, new_j));
                                counter++;
                                streamLinesPixels.push_back(vec2(new_i, new_j)); 
                            }
                            streamLines.pop_front(); 
                        }
                        while(!streamLinesPixels.empty()){
                            vec2 ind = streamLinesPixels.front(); 
                            licTexture[ind[0]][ind[1]] = sum/counter;
                            visited[ind[0]][ind[1]] = 1; 
                            streamLinesPixels.pop_front(); 
                        }
                    }

                }
                licMag[i][j] = glm::length(vectorField.interpolate(startPoint));
                minVal = std::min(minVal, licMag[i][j]);
                maxVal = std::max(maxVal, licMag[i][j]);
            }
        }
    }

    // Adjust contrast
    double mean = 0.0, stdDev = 0.0, P = 0.0, c_sum = 0.0, colorMean;
    int nbPoints = 0;
    for (size_t j = 0; j < texDims_.y; j++) {
        for (size_t i = 0; i < texDims_.x; i++) {
            if (licTexture[i][j] > 0) { // get non-black pixel values
                c_sum += licTexture[i][j];
                P += (licTexture[i][j] * licTexture[i][j]);
                nbPoints++;
            }
        }
    }
    mean = c_sum / nbPoints;
    stdDev = sqrt((P - (nbPoints * mean * mean)) / (nbPoints - 1));
    colorMean = mean;
    if (propContrast.get()) {
        LogProcessorInfo("Contrast" << "pixelSize_x: " << pixelSize_x << "pixelSize_y: " << pixelSize_y << "\n")
        for (size_t j = 0; j < texDims_.y; j++) {
            for (size_t i = 0; i < texDims_.x; i++) {  
                if (licTexture[i][j] > 0) { // set new non-black pixel values
                    licTexture[i][j] = (propMean.get()*255) + (((propStdDev.get()*255) / stdDev) * (licTexture[i][j] - mean));
                }
            }
        }
        colorMean = (propMean.get() * 255);
    }

    // Color the images
    LogProcessorInfo("min " << minVal << "max " << maxVal << "\n");
    for (size_t j = 0; j < texDims_.y; j++) {
        for (size_t i = 0; i < texDims_.x; i++) {
            dvec4 color;
            int val = int(licTexture[i][j]);

            if (propColor.get()) {
                if (licTexture[i][j] > colorMean) {
                    double normalizedVal = (licMag[i][j] - minVal) / (maxVal - minVal);
                    auto scaledcolor = propIsoTransferFunc.get().sample(normalizedVal);
                    color[0] = int(scaledcolor[0] * 255);
                    color[1] = int(scaledcolor[1] * 255);
                    color[2] = int(scaledcolor[2] * 255);
                    color[3] = int(scaledcolor[3] * 255);
                }
                else {
                    color = dvec4(val, val, val, 255);
                }
            }
            else {
                color = dvec4(val, val, val, 255);
            }
                
            //licImage.setPixelGrayScale(size2_t(i, j), val);
            licImage.setPixel(size2_t(i, j), color);


            if(visited[i][j] == 1){
                debugImage.setPixel(size2_t(i, j), dvec4(255, 0, 0, 255)); 

            }
            else {
                debugImage.setPixel(size2_t(i, j), dvec4(0, 255, 0, 255)); 
            }
        }
    }

    debugOut_.setData(debugOutImage);
    licOut_.setData(outImage);
}

}  // namespace inviwo
