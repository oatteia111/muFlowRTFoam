#ifndef MUFLOWRT_WINDOWS_CLI_H
#define MUFLOWRT_WINDOWS_CLI_H

#include "argList.H"
#include "foamVersion.H"
#include "windowsRuntime.H"
#include <iomanip>
#include <iostream>
#include <string>

namespace muflowrtCli
{
struct InformationOption
{
    const char* name;
    const char* aliases[3];
    const char* description;
    const char* link;
};

static const InformationOption informationOptions[] =
{
    {"help", {"-help", "--help", "-h"}, "Show commands and examples, then exit.", nullptr},
    {"version", {"-version", "--version", "-v"}, "Print the muFlowRT version, then exit.", nullptr},
    {"doc", {"-doc", nullptr, nullptr}, "Print the muFlowRT README link, then exit.",
        "https://github.com/oatteia111/muFlowRTFoam#readme"},
    {"srcDoc", {"-srcDoc", nullptr, nullptr}, "Print the muFlowRT source link, then exit.",
        "https://github.com/oatteia111/muFlowRTFoam/blob/main/solvers/muFlowRT.C"}
};

inline std::string canonicalOption(const std::string& argument)
{
    for (const auto& option : informationOptions)
    {
        for (const char* alias : option.aliases)
        {
            if (alias && argument == alias) return option.name;
        }
    }
    return argument.size() > 1 && argument[0] == '-' ? argument.substr(1) : "";
}

inline const InformationOption* informationOption(const std::string& name)
{
    for (const auto& option : informationOptions)
    {
        if (name == option.name) return &option;
    }
    return nullptr;
}

inline void registerOptions()
{
    using Foam::argList;
#ifdef MUFLOWRT_SERIAL_OPENFOAM
    argList::noParallel();
#endif
    for (const auto& option : informationOptions)
    {
        argList::addBoolOption(option.name, option.description);
    }
    argList::addOption("case", "directory",
        "Use this case for configuration, inputs and outputs. Default: current directory.");
    argList::addOption("fileHandler", "handler",
        "Select an OpenFOAM file I/O handler. Bundled default: uncollated; configuration can override it.");
    argList::addOption("libs", "list",
        "Preload compatible DLLs using OpenFOAM list syntax, for example '(\"example.dll\")'.");
    argList::addUsage("noFunctionObjects",
        "Disable configured OpenFOAM function objects (normally enabled).");
}

inline int error(const std::string& message)
{
    std::cerr << "muFlowRT: " << message << ".\nRun muFlowRT --help for usage.\n";
    return 2;
}

inline void printHelp()
{
    using Foam::argList;
    std::cout << "muFlowRT " << MUFLOWRT_VERSION
        << "\nOpenFOAM-" << Foam::FOAMversion << " (build " << Foam::FOAMbuild << ")\n";
#ifdef MUFLOWRT_SERIAL_OPENFOAM
    std::cout << "Serial OpenFOAM; MPI execution is unavailable.\n";
#else
    std::cout << "MPI-enabled OpenFOAM.\n";
#endif
#ifdef _OPENMP
    std::cout << "OpenMP chemistry; OMP_NUM_THREADS controls chemistry threads.\n";
#else
    std::cout << "Single-threaded chemistry.\n";
#endif
    std::cout << "\nUsage: muFlowRT [OPTIONS]\n\nInformational commands (no case required):\n";
    for (const auto& option : informationOptions)
    {
        std::string aliases;
        for (const char* alias : option.aliases)
        {
            if (!alias) continue;
            if (!aliases.empty()) aliases += ", ";
            aliases += alias;
        }
        std::cout << "  " << std::left << std::setw(27) << aliases
            << argList::optionUsage[option.name] << '\n';
    }
    std::cout << "\nSimulation options:\n";
    const Foam::wordList options = argList::validOptions.sortedToc();
    for (const auto& name : options)
    {
        if (informationOption(name)) continue;
        std::string syntax = "-" + std::string(name);
        const Foam::string& parameter = argList::validOptions[name];
        if (!parameter.empty()) syntax += " <" + std::string(parameter) + ">";
        std::cout << "  " << std::left << std::setw(27) << syntax
            << argList::optionUsage[name] << '\n';
    }
    std::cout
        << "\nRelative case paths resolve from the startup directory. Repeated -case uses the first value.\n"
        << "The selected case supplies system/controlDict and all case-relative chemistry files.\n"
        << "Global OpenFOAM configuration is separate; bundled defaults are discovered automatically.\n"
        << "Keep the complete EXE/DLL/configuration bundle together. No wrapper or environment setup is required.\n"
        << "Version and documentation commands must be used alone. Documentation commands print links.\n"
        << "\nPowerShell examples:\n"
        << "  & 'C:\\tools\\muFlowRT\\bin\\muFlowRT.exe' --help\n"
        << "  & 'C:\\tools\\muFlowRT\\bin\\muFlowRT.exe' -version\n"
        << "  & 'C:\\tools\\muFlowRT\\bin\\muFlowRT.exe' -doc\n"
        << "  & 'C:\\tools\\muFlowRT\\bin\\muFlowRT.exe' -case 'C:\\models\\Case B'\n";
}

// Returns -1 for simulation, otherwise the process exit code.
inline int handle()
{
    registerOptions();
    const auto& startup = Foam::windowsRuntime::startup();
    const auto& arguments = startup.arguments;
    const InformationOption* information = nullptr;
    for (std::size_t i = 1; i < arguments.size(); ++i)
    {
        const std::string name = canonicalOption(arguments[i]);
        if (name.empty() || !Foam::argList::validOptions.found(name))
        {
            return error("Unknown argument '" + arguments[i] + "'");
        }
        const auto* option = informationOption(name);
        if (option)
        {
            if (name != "help" && arguments.size() != 2)
                return error("Option '-" + name + "' must be used alone");
            information = option;
        }
        else if (!Foam::argList::validOptions[name].empty())
        {
            if (++i == arguments.size()
             || (!arguments[i].empty() && arguments[i][0] == '-'))
                return error("Option '-" + name + "' requires an argument");
        }
    }
    if (!startup.error.empty()) return error(startup.error);
    if (!information) return -1;
    if (std::string(information->name) == "help") printHelp();
    else if (information->link) std::cout << information->link << '\n';
    else std::cout << "muFlowRT " << MUFLOWRT_VERSION << '\n';
    return 0;
}
}
#endif
