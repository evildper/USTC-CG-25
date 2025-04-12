// Blinn-Phong shader for deferred lighting
#version 430 core

// Define a uniform struct for lights
struct Light {
    // The matrices are used for shadow mapping. You need to fill it according to how we are filling it when building the normal maps (node_render_shadow_mapping.cpp). 
    // Now, they are filled with identity matrix. You need to modify C++ code innode_render_deferred_lighting.cpp.
    // Position and color are filled.
    mat4 light_projection;
    mat4 light_view;
    vec3 position;
    float radius;
    vec3 color; // Just use the same diffuse and specular color.
    int shadow_map_id;
};

layout(binding = 0) buffer lightsBuffer {
    Light lights[4];
};

uniform vec2 iResolution;

uniform sampler2D diffuseColorSampler;
uniform sampler2D normalMapSampler; // You should apply normal mapping in rasterize_impl.fs
uniform sampler2D metallicRoughnessSampler;
uniform sampler2DArray shadow_maps;//深度图采样
uniform sampler2D position;

// uniform float alpha;
uniform vec3 camPos;

uniform int light_count;

layout(location = 0) out vec4 Color;

void main() {
    vec2 uv = gl_FragCoord.xy / iResolution;//获取当前片元的uv坐标

    vec3 pos = texture2D(position, uv).xyz;//通过uv坐标获取当前片元的世界坐标
    vec3 normal = texture2D(normalMapSampler, uv).xyz;//通过uv坐标获取当前片元的法线贴图

    vec4 metalnessRoughness = texture2D(metallicRoughnessSampler, uv);//通过uv坐标获取当前片元的金属度和粗糙度贴图
    float metal = metalnessRoughness.x;//金属度
    float roughness = metalnessRoughness.y;//粗糙度

    vec3 textureColor = texture2D(diffuseColorSampler, uv).xyz;//通过uv坐标获取当前片元的漫反射贴图
    Color = vec4(0, 0, 0, 1);//初始化颜色为黑色

    for(int i = 0; i < light_count; i ++) {
        normal = normalize(normal);//法线归一化
        float shadow_map = texture(shadow_maps, vec3(uv, lights[i].shadow_map_id)).x;//通过uv坐标获取当前片元的阴影贴图
        vec3 L = normalize(lights[i].position - pos);
        vec3 V = normalize(camPos - pos);
        vec3 H = normalize(L + V);

        float NdotL = max(dot(normal, L), 0.0);
        float NdotH = max(dot(normal, H), 0.0);

        // 环境光
        vec3 ambient = 0.2 * lights[i].color*textureColor; // 适当降低环境光强度

        // 漫反射
        vec3 diffuse = (1.0 - 0.8 * metal) * NdotL * lights[i].color * textureColor;

        // 高光
        float shininess = mix(32.0, 128.0, 1.0 - roughness);
        vec3 specular = 0.8*metal * pow(NdotH, shininess) * lights[i].color;

        // Shadow 
        mat4 projection = lights[i].light_projection;
        mat4 view = lights[i].light_view;
        vec4 clipPos =  projection * view * (vec4(pos, 1.0));
        clipPos.w = max(clipPos.w, 0.0001);  //
        float depth_tmp = (clipPos.z / clipPos.w);
        float x_tmp = (clipPos.x / clipPos.w * 0.5) + 0.5;
        float y_tmp = (clipPos.y / clipPos.w * 0.5) + 0.5;
        
        // x_tmp = clamp(x_tmp, 0.0, 1.0);
        // y_tmp = clamp(y_tmp, 0.0, 1.0);
        float bias = max(0.005, 0.01 * length(lights[i].position - pos));
        float closestDepth = texture(shadow_maps, vec3(x_tmp, y_tmp, lights[i].shadow_map_id)).x;
        //float shadow = depth_tmp - bias > closestDepth ? 0.0 : 1.0;
        float shadow = 1.0;        
        vec2 texelSize =  1.0 / textureSize(shadow_maps, 0).xy;
        int filterX = 2;
        if(depth_tmp - bias > closestDepth) {
            float d_bl = 0.0;
            float count_b = 0.0;
            for(int x = -1; x <= 1; ++x) {
                for(int y = -1; y <= 1; ++y) {
                    float pcfDepth = texture(shadow_maps, vec3(x_tmp + x * texelSize.x, y_tmp + y * texelSize.y, lights[i].shadow_map_id)).x;
                    if(depth_tmp - bias > pcfDepth) {
                        count_b += 1;
                        d_bl += pcfDepth;
                    }
                }   
            }
            if(count_b > 0) {
                d_bl = d_bl / count_b;
                filterX = int(25.0 * d_bl / depth_tmp);
                filterX = min(filterX, 25);
                filterX = max(filterX, 0);
            }
        }
        for(int x = -filterX; x <= filterX; ++x) {
            for(int y = -filterX; y <= filterX; ++y) {
                float pcfDepth = texture(shadow_maps, vec3(x_tmp + x * texelSize.x, y_tmp + y * texelSize.y, lights[i].shadow_map_id)).x;
                shadow += depth_tmp - bias > pcfDepth ? 0.0 : 1.0;        
            }    
        }
        float count = 2 * filterX + 1;
        shadow = shadow / (count * count);
        if (x_tmp < 0.0 || x_tmp > 1.0 || y_tmp < 0.0 || y_tmp > 1.0) {
            shadow = 0.0;
        }      
        vec3 result = shadow* (specular + diffuse) + ambient;
        Color += vec4(result, 1.0);
    }
    Color.rgb = clamp(Color.rgb, 0.0, 1.0);
    Color.a = 1.0;
}
       