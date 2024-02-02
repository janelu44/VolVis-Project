#pragma once
#include <array>
#include <glm/vec2.hpp>
#include <glm/vec4.hpp>
#include <cstring> // memcmp  // macOS change TH

namespace render {

enum class RenderMode {
    RenderSlicer,
    RenderMIP,
    RenderIso,
    RenderComposite,
    RenderTF2D,
    RenderAmbientOcclusion
};

struct RenderConfig {
    RenderMode renderMode { RenderMode::RenderComposite };
    glm::ivec2 renderResolution;

    bool volumeShading { false };
    float isoValue { 95.0f };

    bool bisection { true };
    float bisectionErrorThreshold { 0.01f };
    int bisectionMaxIterations { 20 };

    float earlyRayTerminationThreshold { 0.99f };

    bool ambientOcclusionDebug { false };
    int ambientOcclusionDebugShellLevel { 1 };
    int ambientOcclusionMaxShells { 3 };
    float ambientOcclusionNormalFactor { 2.0f };

    // 1D transfer function.
    std::array<glm::vec4, 256> tfColorMap;
    // Used to convert from a value to an index in the color map.
    // index = (value - start) / range * tfColorMap.size();
    float tfColorMapIndexStart;
    float tfColorMapIndexRange;

    // 2D transfer function.
    float TF2DIntensity;
    float TF2DRadius;
    glm::vec4 TF2DColor;
};

// NOTE(Mathijs): should be replaced by C++20 three-way operator (aka spaceship operator) if we require C++ 20 support from Linux users (GCC10 / Clang10).
inline bool operator==(const RenderConfig& lhs, const RenderConfig& rhs)
{
    return std::memcmp(&lhs, &rhs, sizeof(RenderConfig)) == 0;
}
inline bool operator!=(const RenderConfig& lhs, const RenderConfig& rhs)
{
    return !(lhs == rhs);
}

}