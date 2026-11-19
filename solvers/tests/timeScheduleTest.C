#include "../timeSchedule.H"
#include <filesystem>
#include <chrono>
#include <iomanip>
#include <iostream>

namespace fs = std::filesystem;
const double never = 1e30;

void require(bool condition, const char* message)
{
    if (!condition) {throw std::runtime_error(message);}
}

template<class Operation>
void rejects(Operation operation)
{
    bool rejected = false;
    try {operation();}
    catch (const std::runtime_error&) {rejected = true;}
    require(rejected, "invalid input was accepted");
}

void writeFrames(const fs::path& path, const std::vector<double>& times)
{
    std::ofstream data(path, std::ios::binary);
    std::ofstream companion(path.string()+".times");
    companion << std::setprecision(17);
    for (std::size_t index=0; index<times.size(); ++index)
    {
        const float values[4] = {float(times[index]), 0, float(10+index), 0};
        data.write(reinterpret_cast<const char*>(values), sizeof(values));
        companion << times[index] << '\n';
    }
}

void checkReactionCadence
(
    const int transportSteps,
    const std::vector<int>& expectedCalls
)
{
    const int reactionSteps = -5;
    const double endTime = transportSteps*0.4;
    int transportCount = 0;
    double lastReactionTime = 0;
    double reactedTime = 0;
    std::vector<int> calls;
    for (int step=1; step<=transportSteps; ++step)
    {
        const double time = step*0.4;
        ++transportCount;
        if (transportCount > -reactionSteps-1) {transportCount = 0;}
        if (muFlowTime::reactionDue
            (reactionSteps, transportCount, time, never, endTime, lastReactionTime))
        {
            calls.push_back(step);
            reactedTime += time-lastReactionTime;
            lastReactionTime = time;
            require(!muFlowTime::reactionDue
                (reactionSteps, transportCount, time, never, endTime, lastReactionTime),
                "chemistry must not repeat at an already reacted endpoint");
        }
    }
    require(calls == expectedCalls, "normal chemistry cadence changed");
    require(lastReactionTime == endTime, "final chemistry was omitted");
    require(std::abs(reactedTime-endTime) < 1e-12, "chemistry interval was skipped or repeated");
}

int main()
{
    const fs::path scratch = fs::temp_directory_path()/
        ("muflowrt-time-test-"+std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
    fs::create_directory(scratch);
    try
    {
        const fs::path path = scratch/"hfix";
        writeFrames(path, {0, 10, 100.125});
        auto read = [&](int frame, int version=2)
        {
            return muFlowTime::readBoundaryFrame<double>(path.string(), frame, 1, version, never);
        };
        require(read(0).nextTime == 10 && read(0).data[2] == 10, "initial boundary");
        require(read(1).nextTime == 100.125 && read(1).data[2] == 11, "second boundary");
        require(read(2).nextTime == never && read(2).data[2] == 12, "EOF must hold last frame");
        require(read(1, 1).nextTime == 100.125, "legacy binary reader");
        rejects([&] {read(3);});
        rejects([&] {read(-1);});
        rejects([&] {read(0, 3);});

        // These distinct event times collapse to one float32 timestamp.
        writeFrames(path, {0, 100.000001, 100.000002});
        require(read(0).nextTime == 100.000001, "precise first event");
        require(read(1).nextTime == 100.000002, "precise second event");
        rejects([&] {read(0, 1);});
        fs::remove(path.string()+".times");
        rejects([&] {read(0);});
        writeFrames(path, {0, 10});
        {std::ofstream companion(path.string()+".times"); companion << "0\n";}
        rejects([&] {read(0);});
        {std::ofstream companion(path.string()+".times"); companion << "0\n0\n";}
        rejects([&] {read(0);});
        {std::ofstream companion(path.string()+".times"); companion << "0\ninvalid\n";}
        rejects([&] {read(0);});
        writeFrames(path, {0, 10});
        fs::resize_file(path, fs::file_size(path)-1);
        rejects([&] {read(0);});

        auto step = [](double time, double dt, double output, double boundary, double reaction, double end)
        {
            return muFlowTime::stepEnd(time, dt, output, boundary, reaction, end);
        };
        require(step(0, 100, 1, 10, never, 800) == 1, "startup must not skip first output");
        require(step(7.94286, 3.82334, 50, 10, never, 800) == 10, "boundary before output");
        require(step(9, 5, 10, 10, 10, 800) == 10, "simultaneous events");
        require(step(100, 100, 150, never, 112.5, 800) == 112.5, "reaction subdivision");
        require(step(799.999999, 100, 800, never, never, 800) == 800, "final short step");
        require(step(100, 100, 150, never, never, 800) == 150, "restart pending output");
        require(step(0, 100, 100, never, never, 100) == 100, "single output");
        require(step(100.000001, 1, 100.000002, never, never, 800) == 100.000002, "close output times");
        rejects([&] {step(1, 0, 50, never, never, 800);});
        rejects([&] {step(1, 1, 1, never, never, 800);});
        rejects([&] {step(1e20, 1, never, never, never, 1e25);});

        checkReactionCadence(1, {1});
        checkReactionCadence(13, {4, 9, 13});
        checkReactionCadence(14, {4, 9, 14});
        require(!muFlowTime::reactionDue(-5, 0, 0.0, never, 0.0, 0.0),
            "zero-duration final chemistry");
        require(muFlowTime::reactionDue(-5, 1, 5.0, never, 5.0, 4.9116),
            "benchmark final chemistry remainder");

        // Positive subdivisions retain their events, including the final one.
        double chemistryTime = 0;
        double lastReactionTime = 0;
        int reactionCount = 1;
        std::vector<double> reactionTimes;
        while (chemistryTime < 1)
        {
            const double reactionTime = reactionCount*0.25;
            chemistryTime = step(chemistryTime, 0.1, 1.0, never, reactionTime, 1.0);
            if (muFlowTime::reactionDue
                (4, 0, chemistryTime, reactionTime, 1.0, lastReactionTime))
            {
                reactionTimes.push_back(chemistryTime);
                lastReactionTime = chemistryTime;
                ++reactionCount;
            }
        }
        require(reactionTimes == std::vector<double>({0.25, 0.5, 0.75, 1.0}),
            "positive chemistry events changed or final reaction was duplicated");

        // Exercise repeated adaptive advances through the flow3d write schedule.
        double time = 0;
        double dt = 1;
        for (double output=50; output<=800; output+=50)
        {
            while (time < output)
            {
                const double boundary = time < 10 ? 10 : never;
                time = step(time, dt, output, boundary, never, 800);
                dt = std::min(100.0, dt*1.4);
            }
            require(time == output, "accepted output drifted");
        }
        require(time == 800, "solver schedule overshot final output");
        fs::remove_all(scratch);
        std::cout << "time scheduling, chemistry cadence, and boundary reader checks passed\n";
    }
    catch (...)
    {
        fs::remove_all(scratch);
        throw;
    }
}
