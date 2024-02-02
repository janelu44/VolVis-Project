#include "renderer.h"
#include <algorithm>
#include <algorithm> // std::fill
#include <cmath>
#include <functional>
#include <glm/common.hpp>
#include <glm/gtx/component_wise.hpp>
#include <iostream>
#include <tbb/blocked_range2d.h>
#include <tbb/parallel_for.h>
#include <tuple>

namespace render {

// The renderer is passed a pointer to the volume, gradinet volume, camera and an initial renderConfig.
// The camera being pointed to may change each frame (when the user interacts). When the renderConfig
// changes the setConfig function is called with the updated render config. This gives the Renderer an
// opportunity to resize the framebuffer.
Renderer::Renderer(
    const volume::Volume* pVolume,
    const volume::GradientVolume* pGradientVolume,
    const render::RayTraceCamera* pCamera,
    const RenderConfig& initialConfig)
    : m_pVolume(pVolume)
    , m_pGradientVolume(pGradientVolume)
    , m_pCamera(pCamera)
    , m_config(initialConfig)
{
    resizeImage(initialConfig.renderResolution);
}

// Set a new render config if the user changed the settings.
void Renderer::setConfig(const RenderConfig& config)
{
    if (config.renderResolution != m_config.renderResolution)
        resizeImage(config.renderResolution);

    m_config = config;
}

// Resize the framebuffer and fill it with black pixels.
void Renderer::resizeImage(const glm::ivec2& resolution)
{
    m_frameBuffer.resize(size_t(resolution.x) * size_t(resolution.y), glm::vec4(0.0f));
}

// Clear the framebuffer by setting all pixels to black.
void Renderer::resetImage()
{
    std::fill(std::begin(m_frameBuffer), std::end(m_frameBuffer), glm::vec4(0.0f));
}

// Return a VIEW into the framebuffer. This view is merely a reference to the m_frameBuffer member variable.
// This does NOT make a copy of the framebuffer.
gsl::span<const glm::vec4> Renderer::frameBuffer() const
{
    return m_frameBuffer;
}

// Main render function. It computes an image according to the current renderMode.
// Multithreading is enabled in Release/RelWithDebInfo modes. In Debug mode multithreading is disabled to make debugging easier.
void Renderer::render()
{
    resetImage();

    static constexpr float sampleStep = 1.0f;
    const glm::vec3 planeNormal = -glm::normalize(m_pCamera->forward());
    const glm::vec3 volumeCenter = glm::vec3(m_pVolume->dims()) / 2.0f;
    const Bounds bounds { glm::vec3(0.0f), glm::vec3(m_pVolume->dims() - glm::ivec3(1)) };

    // 0 = sequential (single-core), 1 = TBB (multi-core)
#ifdef NDEBUG
    // If NOT in debug mode then enable parallelism using the TBB library (Intel Threaded Building Blocks).
#define PARALLELISM 1
#else
    // Disable multi threading in debug mode.
#define PARALLELISM 0
#endif

#if PARALLELISM == 0
    // Regular (single threaded) for loops.
    for (int x = 0; x < m_config.renderResolution.x; x++) {
        for (int y = 0; y < m_config.renderResolution.y; y++) {
#else
    // Parallel for loop (in 2 dimensions) that subdivides the screen into tiles.
    const tbb::blocked_range2d<int> screenRange { 0, m_config.renderResolution.y, 0, m_config.renderResolution.x };
        tbb::parallel_for(screenRange, [&](tbb::blocked_range2d<int> localRange) {
        // Loop over the pixels in a tile. This function is called on multiple threads at the same time.
        for (int y = std::begin(localRange.rows()); y != std::end(localRange.rows()); y++) {
            for (int x = std::begin(localRange.cols()); x != std::end(localRange.cols()); x++) {
#endif
            // Compute a ray for the current pixel.
            const glm::vec2 pixelPos = glm::vec2(x, y) / glm::vec2(m_config.renderResolution);
            Ray ray = m_pCamera->generateRay(pixelPos * 2.0f - 1.0f);

            // Compute where the ray enters and exists the volume.
            // If the ray misses the volume then we continue to the next pixel.
            if (!instersectRayVolumeBounds(ray, bounds))
                continue;

            // Get a color for the current pixel according to the current render mode.
            glm::vec4 color {};
            switch (m_config.renderMode) {
            case RenderMode::RenderSlicer: {
                color = traceRaySlice(ray, volumeCenter, planeNormal);
                break;
            }
            case RenderMode::RenderMIP: {
                color = traceRayMIP(ray, sampleStep);
                break;
            }
            case RenderMode::RenderComposite: {
                color = traceRayComposite(ray, sampleStep);
                break;
            }
            case RenderMode::RenderCompositeEnhancedOp: {
                color = traceRayCompositeEnhancedOp(ray, sampleStep);
                break;
            }
            case RenderMode::RenderIso: {
                color = traceRayISO(ray, sampleStep);
                break;
            }
            case RenderMode::RenderTF2D: {
                color = traceRayTF2D(ray, sampleStep);
                break;
            }
            };
            // Write the resulting color to the screen.
            fillColor(x, y, color);

#if PARALLELISM == 1
        }
    }
});
#else
            }
        }
#endif
}

// ======= DO NOT MODIFY THIS FUNCTION ========
// This function generates a view alongside a plane perpendicular to the camera through the center of the volume
//  using the slicing technique.
glm::vec4 Renderer::traceRaySlice(const Ray& ray, const glm::vec3& volumeCenter, const glm::vec3& planeNormal) const
{
    const float t = glm::dot(volumeCenter - ray.origin, planeNormal) / glm::dot(ray.direction, planeNormal);
    const glm::vec3 samplePos = ray.origin + ray.direction * t;
    const float val = m_pVolume->getSampleInterpolate(samplePos);
    return glm::vec4(glm::vec3(std::max(val / m_pVolume->maximum(), 0.0f)), 1.f);
}

// ======= DO NOT MODIFY THIS FUNCTION ========
// Function that implements maximum-intensity-projection (MIP) raycasting.
// It returns the color assigned to a ray/pixel given it's origin, direction and the distances
// at which it enters/exits the volume (ray.tmin & ray.tmax respectively).
// The ray must be sampled with a distance defined by the sampleStep
glm::vec4 Renderer::traceRayMIP(const Ray& ray, float sampleStep) const
{
    float maxVal = 0.0f;

    // Incrementing samplePos directly instead of recomputing it each frame gives a measureable speed-up.
    glm::vec3 samplePos = ray.origin + ray.tmin * ray.direction;
    const glm::vec3 increment = sampleStep * ray.direction;
    for (float t = ray.tmin; t <= ray.tmax; t += sampleStep, samplePos += increment) {
        const float val = m_pVolume->getSampleInterpolate(samplePos);
        maxVal = std::max(val, maxVal);
    }

    // Normalize the result to a range of [0 to mpVolume->maximum()].
    return glm::vec4(glm::vec3(maxVal) / m_pVolume->maximum(), 1.0f);
}

// ======= TODO: IMPLEMENT ========
// This function should find the position where the ray intersects with the volume's isosurface.
// If volume shading is DISABLED then simply return the isoColor.
// If volume shading is ENABLED then return the phong-shaded color at that location using the local gradient (from m_pGradientVolume).
//   Use the camera position (m_pCamera->position()) as the light position.
// Use the bisectionAccuracy function (to be implemented) to get a more precise isosurface location between two steps.
glm::vec4 Renderer::traceRayISO(const Ray& ray, float sampleStep) const
{
    //static constexpr glm::vec3 isoColor { 0.7f, 0.3f, 0.5f };
    static constexpr glm::vec3 isoColor { 1.0f, 1.0f, 1.0f };

    glm::vec3 samplePos = ray.origin + ray.tmin * ray.direction;
    const glm::vec3 increment = sampleStep * ray.direction;
    for (float t = ray.tmin; t <= ray.tmax; t += sampleStep, samplePos += increment) {
        const float val = m_pVolume->getSampleInterpolate(samplePos);
        if (val > m_config.isoValue) {
            const auto bisectedPos = samplePos - (t - bisectionAccuracy(ray, t - sampleStep, t, m_config.isoValue)) * ray.direction;
            if (m_config.volumeShading) {
                return glm::vec4(
                    computePhongShading(
                        isoColor,
                        m_pGradientVolume->getGradientInterpolate(bisectedPos),
                        m_pCamera->position() - bisectedPos,
                        -ray.direction),
                    1.0f);
            } else if (m_config.toneShading) {
                return glm::vec4(computeToneShading(isoColor, m_pGradientVolume->getGradientInterpolate(bisectedPos), m_pCamera->position() - bisectedPos, -ray.direction), 1.0f);
            } else if (m_config.approxToneShading) {
                return glm::vec4(approxToneShading(isoColor, m_pGradientVolume->getGradientInterpolate(bisectedPos), m_pCamera->position() - bisectedPos, -ray.direction), 1.0f);
            } else if (m_config.approxToneShadingHighlights) {
                return glm::vec4(approxToneShadingHighlights(isoColor, m_pGradientVolume->getGradientInterpolate(bisectedPos), m_pCamera->position() - bisectedPos, -ray.direction), 1.0f);
            }
            
            return glm::vec4(isoColor, 1.0f);
        }
    }

    return glm::vec4(0.0f);
}

// ======= TODO: IMPLEMENT ========
// Given that the iso value lies somewhere between t0 and t1, find a t for which the value
// closely matches the iso value (less than 0.01 difference). Add a limit to the number of
// iterations such that it does not get stuck in degerate cases.
float Renderer::bisectionAccuracy(const Ray& ray, float t0, float t1, float isoValue) const
{
    const int maxIterations = m_config.bisection ? m_config.bisectionMaxIterations : 0;

    float lo = t0, t = t1, hi = t1;
    for (auto i = 0; i < maxIterations; i++) {
        t = (lo + hi) / 2.0f;
        const auto val = m_pVolume->getSampleInterpolate(ray.origin + t * ray.direction);
        if (glm::abs(val - isoValue) < m_config.bisectionErrorThreshold)
            return t;
        if (val > isoValue)
            hi = t;
        if (val < isoValue)
            lo = t;
    }

    return t;
}

// ======= TODO: IMPLEMENT ========
// Compute Phong Shading given the voxel color (material color), the gradient, the light vector and view vector.
// You can find out more about the Phong shading model at:
// https://en.wikipedia.org/wiki/Phong_reflection_model
//
// Use the given color for the ambient/specular/diffuse (you are allowed to scale these constants by a scalar value).
// You are free to choose any specular power that you'd like.
glm::vec3 Renderer::computePhongShading(const glm::vec3& color, const volume::GradientVoxel& gradient, const glm::vec3& L, const glm::vec3& V)
{   
    const auto k = glm::vec3(0.1f, 0.7f, 0.2f);
    const float alpha = 100.0f;

    auto cosTheta = glm::max(0.0f, glm::dot(glm::normalize(-L), glm::normalize(gradient.dir)));
    auto cosPhi = glm::max(0.0f, glm::dot(glm::normalize(glm::reflect(-L, gradient.dir)), glm::normalize(V)));

    return glm::dot(k, glm::vec3(1.0f, cosTheta, glm::pow(cosPhi, alpha))) * color;
}

// This function computes tone shading by approximating it using the Phong model
glm::vec3 Renderer::approxToneShading(const glm::vec3& color, const volume::GradientVoxel& gradient, const glm::vec3& L, const glm::vec3& V) {
    
    // Phong constants. Can be adjusted as needed
    float kd = 0.7f;

    // Blue and yellow colors. Adjust values to your case to determine the strength of temperature shift
    glm::vec3 blueColor(0.0f, 0.0f, 0.7f);
    glm::vec3 yellowColor(0.4f, 0.4f, 0.0f);

    // Proeminence of object color and luminance strength. These can be adjusted to your case
    float alpha = 0.25f;
    float beta = 0.55f;

    // Tone creation
    const glm::vec3 cool = blueColor + alpha * kd;
    const glm::vec3 warm = yellowColor + beta * kd;

  
    // Calculate light vector which is perpendicular to the gaze (using the gradient of the surface doesn't give good results, so I take a world vector). Paper suggests to have the light vector perpendicular to the gaze.
    // This vector direction can be adjusted as needed. 
    glm::vec3 worldUpVector = glm::vec3(0.0f, 1.0f, 0.0f);
    glm::vec3 lightVector = glm::cross(worldUpVector, glm::normalize(V));

    // Calculate light intensities and ambient term
    glm::vec3 light1 = 0.5f * (warm - cool);
    glm::vec3 light2 = 0.5f * (cool - warm);
    glm::vec3 ambient = 0.5f * (cool + warm);

    // Calculate the diffuse shading terms
    float cosTheta1 = glm::dot(glm::normalize(-lightVector), glm::normalize(gradient.dir));
    float diffuse1 = glm::max(0.0f, cosTheta1);

    float cosTheta2 = glm::dot(glm::normalize(lightVector), glm::normalize(gradient.dir));
    float diffuse2 = glm::max(0.0f, cosTheta2);

    // Combine the shading terms
    glm::vec3 shading = ambient + light1 * diffuse1 + light2 * diffuse2;

    // Paper assumes that for the approximation, the color of the object to be set to white
    // So multiplying with the color here does not make any major change
    // However you can see the difference if you uncomment first line in traceRayISO => less cool to warm shift => depth/shape not as easy to distinguish.
    return shading * color;

}


// This function computes tone shading by approximating it using the Phong model and adds highlights
// Tone shading is primarily based on the diffuse terms. The paper doesn't include highlights in their phong approximation.
// Paper mentiones that highlights could be added on systems with accumulation buffers
// So this is just an extra addition from my side just to check how it looks like (it is not supposed to give better results than the standard tone shading)
glm::vec3 Renderer::approxToneShadingHighlights(const glm::vec3& color, const volume::GradientVoxel& gradient, const glm::vec3& L, const glm::vec3& V)
{

    // Phong constants. Can be adjusted as needed
    float kd = 0.7f;
    float ks = 0.2f; 
    float shininess = 70.0f;

    // Blue and yellow colors. Adjust values to your case to determine the strength of temperature shift
    glm::vec3 blueColor(0.0f, 0.0f, 0.7f);
    glm::vec3 yellowColor(0.4f, 0.4f, 0.0f);

    // Proeminence of object color and luminance strength. These can be adjusted to your case
    float alpha = 0.25f;
    float beta = 0.55f;

    // Tone creation
    const glm::vec3 cool = blueColor + alpha * kd;
    const glm::vec3 warm = yellowColor + beta * kd;

    // Calculate light vector which is perpendicular to the gaze (using the gradient of the surface doesn't give good results, so I take a world vector)
    glm::vec3 worldUpVector = glm::vec3(0.0f, 1.0f, 0.0f);
    glm::vec3 lightVector = glm::cross(worldUpVector, glm::normalize(V));

    // Calculate light intensities and ambient term
    glm::vec3 light1 = 0.5f * (warm - cool);
    glm::vec3 light2 = 0.5f * (cool - warm);
    glm::vec3 ambient = 0.5f * (cool + warm);

    // Calculate the diffuse shading terms
    float cosTheta1 = glm::dot(glm::normalize(-lightVector), glm::normalize(gradient.dir));
    float diffuse1 = glm::max(0.0f, cosTheta1);

    float cosTheta2 = glm::dot(glm::normalize(lightVector), glm::normalize(gradient.dir));
    float diffuse2 = glm::max(0.0f, cosTheta2);

    // Calculate the specular terms
    float cosPhi1 = glm::max(0.0f, glm::dot(glm::normalize(glm::reflect(-lightVector, gradient.dir)), glm::normalize(V)));
    float specular1 = ks * glm::pow(cosPhi1, shininess);

    float cosPhi2 = glm::max(0.0f, glm::dot(glm::normalize(glm::reflect(lightVector, gradient.dir)), glm::normalize(V)));
    float specular2 = ks * glm::pow(cosPhi2, shininess);

    // Combine the shading terms
    glm::vec3 shading = ambient + light1 * diffuse1 + light2 * diffuse2 + specular1 + specular2;

    return shading * color;
}

// This function computes tone shading using new model
glm::vec3 Renderer::computeToneShading(const glm::vec3& color, const volume::GradientVoxel& gradient, const glm::vec3& L, const glm::vec3& V)
{

    // Blue and yellow colors. Adjust values to your case to determine the strength of temperature shift. These can be adjusted as needed.
    glm::vec3 blueColor(0.0f, 0.0f, 0.7f);
    glm::vec3 yellowColor(0.4f, 0.4f, 0.0f);

    // These can be adjusted to your case.
    float kd = 0.7f;
    float ka = 0.1f;

    // Proeminence of object color and luminance strength. These can be adjusted to your case.
    float alpha = 0.25f;
    float beta = 0.55f;

    // Tone creation
    const glm::vec3 cool = blueColor + alpha * kd;
    const glm::vec3 warm = yellowColor + beta * kd;

    // The light vector should be perpendicular to the gaze vector in order for the dot product between light and normal to fully vary between [-1,1], as suggested in the paper.
    // Vector direction can be adapted as you want.
    glm::vec3 worldUpVector = glm::vec3(0.0f, 1.0f, 0.0f);
    glm::vec3 lightVector = glm::cross(worldUpVector, glm::normalize(V));
    
    // Since we want a full variation between [-1,1], I will not cap the value of the dot product to 0. However, in comparison to the approxToneShading, when you switch to composite it won't look as good. This is because I do not cap negative values to 0. This is why the paper also suggests the approximation using Phong.
    float cosTheta = glm::dot(glm::normalize(lightVector), glm::normalize(gradient.dir));

    // In case you want to add some ambient light (it will make the image lighter) uncomment the next line and comment the penultimate one.
    //const glm::vec3 pixelColor = ka + ((((1 + cosTheta) / 2) * cool) + ((1 - ((1 + cosTheta) / 2)) * warm));

    // You can also take into account the original color of the object. However, now it's set to white (as suggested in the paper), so it doesn't make any difference.
    // If you want to see the effect, uncomment first line in traceRayISO which sets pink as the isocolor. Then uncomment the line below and comment the penultimate one.
    // The color of the object will combine with the cool-warm tone, which will still show a sense of depth, but not as good as the blue-yellow hue shift.
    // const glm::vec3 pixelColor = ((((1 + cosTheta) / 2) * cool) + ((1 - ((1 + cosTheta) / 2)) * warm)) * color;

    const glm::vec3 pixelColor = ((((1 + cosTheta) / 2) * cool) + ((1 - ((1 + cosTheta) / 2)) * warm));
    return pixelColor;
}


// ======= TODO: IMPLEMENT ========
// In this function, implement 1D transfer function raycasting.
// Use getTFValue to compute the color for a given volume value according to the 1D transfer function.
glm::vec4 Renderer::traceRayComposite(const Ray& ray, float sampleStep) const
{
    glm::vec3 samplePos = ray.origin + ray.tmin * ray.direction;
    const glm::vec3 increment = sampleStep * ray.direction;

    auto color = glm::vec4(0.0f);

    for (float t = ray.tmin; t <= ray.tmax; t += sampleStep, samplePos += increment) {
        const float val = m_pVolume->getSampleInterpolate(samplePos);

        auto tfValue = getTFValue(val);
 

        if (m_config.volumeShading) {
            tfValue = glm::vec4(
                computePhongShading(
                    tfValue,
                    m_pGradientVolume->getGradientInterpolate(samplePos),
                    m_pCamera->position() - samplePos,
                    -ray.direction),
                tfValue.w);
        } else if (m_config.toneShading) {
            tfValue = glm::vec4(
                computeToneShading(
                    tfValue,
                    m_pGradientVolume->getGradientInterpolate(samplePos),
                    m_pCamera->position() - samplePos,
                    -ray.direction),
                tfValue.w);
        } else if (m_config.approxToneShading) {
            tfValue = glm::vec4(
                approxToneShading(
                    tfValue,
                    m_pGradientVolume->getGradientInterpolate(samplePos),
                    m_pCamera->position() - samplePos,
                    -ray.direction),
                tfValue.w);
        } else if (m_config.approxToneShadingHighlights) {
            tfValue = glm::vec4(
                approxToneShadingHighlights(
                    tfValue,
                    m_pGradientVolume->getGradientInterpolate(samplePos),
                    m_pCamera->position() - samplePos,
                    -ray.direction),
                tfValue.w);
        }
        const auto tfColor = tfValue * glm::vec4(glm::vec3(tfValue.w), 1.0f);

        color += (1.0f - color.w) * tfColor;

        if (color.w > m_config.earlyRayTerminationThreshold)
            return color;
    }

    return color;
}

// 1D Transfer function raycasting with enhanced boundaries based on gradient magnitude
glm::vec4 Renderer::traceRayCompositeEnhancedOp(const Ray& ray, float sampleStep) const
{
    glm::vec3 samplePos = ray.origin + ray.tmin * ray.direction;
    const glm::vec3 increment = sampleStep * ray.direction;

    auto color = glm::vec4(0.0f);

    for (float t = ray.tmin; t <= ray.tmax; t += sampleStep, samplePos += increment) {
        const float val = m_pVolume->getSampleInterpolate(samplePos);

        auto tfValue = getTFValue(val);

        // Calculate gradient + final opacity value. Kgc, Kgs, Kge can be adjusted as needed. A higher Kgs will define boundaries even better, but will darken the inside of the fish.
        const volume::GradientVoxel& gradient = m_pGradientVolume->getGradientInterpolate(samplePos);
        double power = glm::pow(gradient.magnitude, 0.3);
        float opacity = tfValue.w * (0.2 + 1.4 * power);
        glm::vec4 tfColor;

        if (m_config.volumeShading) {
            tfValue = glm::vec4(
                computePhongShading(
                    tfValue,
                    m_pGradientVolume->getGradientInterpolate(samplePos),
                    m_pCamera->position() - samplePos,
                    -ray.direction),
                opacity);
           
        } else if (m_config.toneShading) {
            tfValue = glm::vec4(
                computeToneShading(
                    tfValue,
                    m_pGradientVolume->getGradientInterpolate(samplePos),
                    m_pCamera->position() - samplePos,
                    -ray.direction),
                opacity);
           
        } else if (m_config.approxToneShading) {
            tfValue = glm::vec4(
                approxToneShading(
                    tfValue,
                    m_pGradientVolume->getGradientInterpolate(samplePos),
                    m_pCamera->position() - samplePos,
                    -ray.direction),
                opacity);
            
        } else if (m_config.approxToneShadingHighlights) {
            tfValue = glm::vec4(
                approxToneShadingHighlights(
                    tfValue,
                    m_pGradientVolume->getGradientInterpolate(samplePos),
                    m_pCamera->position() - samplePos,
                    -ray.direction),
                opacity);          
        } 
            
            tfColor = tfValue * glm::vec4(glm::vec3(tfValue.w), 0.8f);
            // I made multiple alternatives. Uncomment each line to see the differences. This parameter can be adjusted based on your needs.
            //tfColor = tfValue * glm::vec4(glm::vec3(tfValue.w), 1.0f);
            //tfColor = tfValue * glm::vec4(glm::vec3(tfValue.w), 0.6f);
            //tfColor = tfValue * glm::vec4(glm::vec3(tfValue.w), 0.4f);
            //tfColor = tfValue * glm::vec4(glm::vec3(tfValue.w), opacity * 5);

        color += (1.0f - color.w) * tfColor;

        if (color.w > m_config.earlyRayTerminationThreshold)
            return color;
    }

    return color;
}

// ======= DO NOT MODIFY THIS FUNCTION ========
// Looks up the color+opacity corresponding to the given volume value from the 1D tranfer function LUT (m_config.tfColorMap).
// The value will initially range from (m_config.tfColorMapIndexStart) to (m_config.tfColorMapIndexStart + m_config.tfColorMapIndexRange) .
glm::vec4 Renderer::getTFValue(float val) const
{
    // Map value from [m_config.tfColorMapIndexStart, m_config.tfColorMapIndexStart + m_config.tfColorMapIndexRange) to [0, 1) .
    const float range01 = (val - m_config.tfColorMapIndexStart) / m_config.tfColorMapIndexRange;
    const size_t i = std::min(static_cast<size_t>(range01 * static_cast<float>(m_config.tfColorMap.size())), m_config.tfColorMap.size() - 1);
    return m_config.tfColorMap[i];
}

// ======= TODO: IMPLEMENT ========
// In this function, implement 2D transfer function raycasting.
// Use the getTF2DOpacity function that you implemented to compute the opacity according to the 2D transfer function.
glm::vec4 Renderer::traceRayTF2D(const Ray& ray, float sampleStep) const
{
    glm::vec3 samplePos = ray.origin + ray.tmin * ray.direction;
    const glm::vec3 increment = sampleStep * ray.direction;

    auto color = glm::vec4(0.0f);

    for (float t = ray.tmin; t <= ray.tmax; t += sampleStep, samplePos += increment) {
        const float val = m_pVolume->getSampleInterpolate(samplePos);
        const auto gradient = m_pGradientVolume->getGradientInterpolate(samplePos);

        auto tfValue = m_config.TF2DColor;
        tfValue.w *= getTF2DOpacity(val, gradient.magnitude);

        if (m_config.volumeShading) {
            tfValue = glm::vec4(
                computePhongShading(
                    tfValue,
                    gradient,
                    m_pCamera->position() - samplePos,
                    -ray.direction),
                tfValue.w);
        } else if (m_config.toneShading) {
            tfValue = glm::vec4(
                computeToneShading(
                    tfValue,
                    gradient,
                    m_pCamera->position() - samplePos,
                    -ray.direction),
                tfValue.w);
        } else if (m_config.approxToneShading) {
            tfValue = glm::vec4(
                approxToneShading(
                    tfValue,
                    gradient,
                    m_pCamera->position() - samplePos,
                    -ray.direction),
                tfValue.w);
        } else if (m_config.approxToneShadingHighlights) {
            tfValue = glm::vec4(
                approxToneShadingHighlights(
                    tfValue,
                    gradient,
                    m_pCamera->position() - samplePos,
                    -ray.direction),
                tfValue.w);
        }
        const auto tfColor = tfValue * glm::vec4(glm::vec3(tfValue.w), 1.0f);

        color += (1.0f - color.w) * tfColor;

        if (color.w > m_config.earlyRayTerminationThreshold)
            return color;
    }

    return color;
}

// ======= TODO: IMPLEMENT ========
// This function should return an opacity value for the given intensity and gradient according to the 2D transfer function.
// Calculate whether the values are within the radius/intensity triangle defined in the 2D transfer function widget.
// If so: return a tent weighting as described in the assignment
// Otherwise: return 0.0f
//
// The 2D transfer function settings can be accessed through m_config.TF2DIntensity and m_config.TF2DRadius.
float Renderer::getTF2DOpacity(float intensity, float gradientMagnitude) const
{
    const float normalizedMagnitude = gradientMagnitude / m_pGradientVolume->maxMagnitude();
    const float maxDistance = normalizedMagnitude * m_config.TF2DRadius;
    const float distance = glm::abs(intensity - m_config.TF2DIntensity);

    if (distance == 0.0f && maxDistance == 0.0f)
        return 1.0f;

    if (distance < maxDistance)
        return 1.0f - (distance / maxDistance);

    return 0.0f;
}

// This function computes if a ray intersects with the axis-aligned bounding box around the volume.
// If the ray intersects then tmin/tmax are set to the distance at which the ray hits/exists the
// volume and true is returned. If the ray misses the volume the the function returns false.
//
// If you are interested you can learn about it at.
// https://www.scratchapixel.com/lessons/3d-basic-rendering/minimal-ray-tracer-rendering-simple-shapes/ray-box-intersection
bool Renderer::instersectRayVolumeBounds(Ray& ray, const Bounds& bounds) const
{
    const glm::vec3 invDir = 1.0f / ray.direction;
    const glm::bvec3 sign = glm::lessThan(invDir, glm::vec3(0.0f));

    float tmin = (bounds.lowerUpper[sign[0]].x - ray.origin.x) * invDir.x;
    float tmax = (bounds.lowerUpper[!sign[0]].x - ray.origin.x) * invDir.x;
    const float tymin = (bounds.lowerUpper[sign[1]].y - ray.origin.y) * invDir.y;
    const float tymax = (bounds.lowerUpper[!sign[1]].y - ray.origin.y) * invDir.y;

    if ((tmin > tymax) || (tymin > tmax))
        return false;
    tmin = std::max(tmin, tymin);
    tmax = std::min(tmax, tymax);

    const float tzmin = (bounds.lowerUpper[sign[2]].z - ray.origin.z) * invDir.z;
    const float tzmax = (bounds.lowerUpper[!sign[2]].z - ray.origin.z) * invDir.z;

    if ((tmin > tzmax) || (tzmin > tmax))
        return false;

    ray.tmin = std::max(tmin, tzmin);
    ray.tmax = std::min(tmax, tzmax);
    return true;
}

// This function inserts a color into the framebuffer at position x,y
void Renderer::fillColor(int x, int y, const glm::vec4& color)
{
    const size_t index = static_cast<size_t>(m_config.renderResolution.x * y + x);
    m_frameBuffer[index] = color;
}
}