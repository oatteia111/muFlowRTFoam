#ifndef MUFLOWRT_TIME_SCHEDULE_H
#define MUFLOWRT_TIME_SCHEDULE_H

#include <algorithm>
#include <cmath>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace muFlowTime
{

template<class Scalar>
struct BoundaryFrame
{
    Scalar nextTime;
    std::vector<float> data;
};

// The binary cell payload is unchanged; version 2 supplies precise frame times.
template<class Scalar>
BoundaryFrame<Scalar> readBoundaryFrame
(
    const std::string& path,
    const int frame,
    const int cellCount,
    const int version,
    const Scalar noNextTime
)
{
    const auto fail = [&path](const std::string& reason)
    {
        throw std::runtime_error(path + ": " + reason);
    };
    if (cellCount <= 0 || frame < 0)
    {
        fail("invalid boundary cell count or frame index");
    }
    if (version != 1 && version != 2)
    {
        fail("unsupported timeScheduleVersion");
    }
    const std::streamoff frameBytes = std::streamoff(cellCount)*4*sizeof(float);
    std::ifstream binary(path, std::ios::binary | std::ios::ate);
    if (!binary)
    {
        fail("cannot open boundary data");
    }
    const std::streamoff bytes = binary.tellg();
    if (bytes <= 0 || bytes % frameBytes != 0)
    {
        fail("empty or incomplete boundary frame");
    }
    const std::size_t count = static_cast<std::size_t>(bytes/frameBytes);
    if (static_cast<std::size_t>(frame) >= count)
    {
        fail("boundary frame is beyond end of file");
    }

    std::vector<Scalar> times;
    times.reserve(count);
    if (version == 2)
    {
        std::ifstream companion(path + ".times");
        if (!companion)
        {
            fail("timeScheduleVersion 2 requires a .times companion");
        }
        Scalar value;
        while (companion >> value)
        {
            times.push_back(value);
        }
        if (!companion.eof() || times.size() != count)
        {
            fail("invalid .times contents or frame count mismatch");
        }
    }
    else
    {
        for (std::size_t index = 0; index < count; ++index)
        {
            float value;
            binary.seekg(std::streamoff(index)*frameBytes);
            if (!binary.read(reinterpret_cast<char*>(&value), sizeof(value)))
            {
                fail("cannot read legacy boundary timestamp");
            }
            times.push_back(Scalar(value));
        }
    }
    if (times.front() != Scalar(0))
    {
        fail("first boundary frame must start at zero");
    }
    for (std::size_t index = 0; index < count; ++index)
    {
        if (!std::isfinite(times[index]) || times[index] >= noNextTime
            || (index && times[index] <= times[index-1]))
        {
            fail("boundary timestamps must be finite and strictly increasing");
        }
    }

    BoundaryFrame<Scalar> result;
    result.nextTime = static_cast<std::size_t>(frame)+1 < count
        ? times[frame+1] : noNextTime;
    result.data.resize(static_cast<std::size_t>(cellCount)*4);
    binary.seekg(std::streamoff(frame)*frameBytes);
    if (!binary.read(reinterpret_cast<char*>(result.data.data()), frameBytes))
    {
        fail("cannot read complete boundary frame");
    }
    for (int cell = 0; cell < cellCount; ++cell)
    {
        if (!std::isfinite(result.data[cell*4])
            || result.data[cell*4] != static_cast<float>(times[frame]))
        {
            fail("binary and companion frame timestamps disagree");
        }
    }
    return result;
}

template<class Scalar>
Scalar stepEnd
(
    const Scalar current,
    const Scalar numericalDelta,
    const Scalar writeTime,
    const Scalar boundaryTime,
    const Scalar reactionTime,
    const Scalar endTime
)
{
    if (!std::isfinite(numericalDelta) || numericalDelta <= 0)
    {
        throw std::runtime_error("numerical timestep must be finite and positive");
    }
    const Scalar target = std::min
    (
        std::min(current + numericalDelta, endTime),
        std::min(writeTime, std::min(boundaryTime, reactionTime))
    );
    if (!std::isfinite(target) || target <= current)
    {
        throw std::runtime_error("scheduled timestep cannot advance time");
    }
    return target;
}

}

#endif
