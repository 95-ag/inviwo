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
 * Primary author : Daniel J�nsson
 *
 **********************************************************************/

#ifndef LIGHTSOURCE_CL
#define LIGHTSOURCE_CL

typedef enum LightSourceType {
    LIGHT_AREA = 0,
    LIGHT_CONE,
    LIGHT_POINT,
    LIGHT_DIRECTIONAL
} LightSourceType;

// Note that largest variables should be placed first 
// in order to ensure struct size
typedef struct PhotonLightSource {
    float16 tm; // Transformation matrix from local to world coordinates
    float3 radiance; // Note that float3 occupies same space as float4
    float2 size; // width, height of light source, used by area light
    int type; // LightSourceType, use integer to handle size of struct easier
    float area; // area of light source
    float cosFOV;  // cos( (field of view)/2 ), used by cone light
    
    int padding[7];
} PhotonLightSource;

#endif // LIGHTSOURCE_CL