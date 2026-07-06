#pragma once

class PreprocessApplication
{
public:
    // Top-level workflow coordinator. Keeps main.cpp small and leaves
    // mesh algorithms and output formats in their existing modules.
    int Run(int argc, char **argv);
};
