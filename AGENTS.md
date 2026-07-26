# Signalizer

Real-time audio visualization plugin targeting VST/VST3/AU/Standalone emphasizing gridless engineering precision and playful research & discovery of audio.

Built in C++ using JUCE 8, located in `External/juce`.
Signalizer is largely a test-bed / shell on top of a base cross-platform library located in `External/cpl`, where any reusable utility for audio plugin development goes.

Signalizer contains a vectorscope, an oscilloscope and a spectrum/spectrogram.

## Build

```
python prepare.py                # First-time setup: init/update submodules (+ Linux deps)
python build.py dev              # Debug Standalone — use this for iteration
python build.py release          # All formats, packaged into /Releases
```
By default output is quiet, full logs are printed and found in `Make/Logs/` or use `--verbose` for full compiler output (`--help` is also available).

VST2 formats need the SDK at `../SDKs/vstsdk2.4` (`release -s` skips VST2).

## Structure

- `Source/` — plugin source (see layout below)
- `Builds/` — IDE project files (per platform)
- `Make/` — build scripts (`common.py`, `build_win/osx/linux.py`)
- `External/` — submodules (`juce`, `cpl`)
- `JuceLibraryCode/` — generated JUCE glue (do not edit)
- `Signalizer.jucer` — Projucer project file; regenerates `Builds/` and `JuceLibraryCode/`

### `Source/` layout

- `Common/` — shared infra: `HostGraph`, `MixGraphListener`, `SentientViewState`, and the MVC base classes in `CommonSignalizer.h` (`ProcessorState`, `ProcessorStreamState`, `SystemView`, `StateEditor`)
- `Processor/` — `AudioProcessor` (the JUCE `AudioProcessor`, in `PluginProcessor.cpp`)
- `Editor/` — `MainEditor`, `GraphEditor`
- `Config/` — `SignalizerConfiguration.h`
- `Oscilloscope/`, `Spectrum/`, `Vectorscope/` — one folder per visualization
- `Unity/SignalizerSource.cpp` — unity-include build for embedding the visualizations elsewhere

Each visualization folder follows the same five-way split (mapping onto the MVC description below):

- `*Parameters.h` — the `…Content` persistent state class
- `*.cpp` / `*.h` — the visualization class
- `*Controller.cpp` — the editor UI
- `*Rendering.cpp` — OpenGL / 2D paint
- `*DSP.cpp` or `*DSP.inl` — signal processing

## Signal Flow and Architecture

1. `PluginProcessor` receives the signal from device or DAW
2. The signal is sent into one of two `cpl::AudioStream`s (`realtimeInput`)
3. The `HostGraph` connects to the `realtimeOutput` of the stream, publishing the data (timestamped) to later be available to any other Signalizer (including itself)
4. The `HostGraph` discovers and remembers other Signalizers (persistent IDs) and receives mapped input data (configured in `GraphEditor`) through `MixGraphListener` to the visualization audio stream (`presentationInput`), enabling sidechain/mixing of any Signalizer's output into another for comparison and overlay
5. Each visualization class is split into something similar to "model-view-controller", eg.
	- The persistent State (built from simple parameters exposed to the DAW): `OscilloscopeContent`
	- The Visualization: `Oscilloscope`
	- The Editor UI: `OscilloscopeController`
6. MainEditor orchestrates the tabbing system of editor UIs, global settings UI, and view instantiation by forwarding `presentationOutput` through `SentientViewState`, which defers view instantiation and manages serialization on creation/destruction.
7. Each Visualization has a `ProcessorShell` (an `AudioStream::Listener`) that processes audio asynchronously into `StreamState`, held as a separate `shared_ptr` that the `AudioStream` safely releases to avoid deadlocks during plugin instance destruction.
8. In an OpenGL thread, each Visualization's `onOpenGLRendering` ingests parameter updates via `handleFlagUpdates()`, renders audio data from the stream/`StreamState`, and paints 2D grids and text overlay.
9. Each Visualization typically does runtime CPU feature detection to dispatch to various template instantiations taking advantage of SIMD/FMA etc. for graphics and DSP.

## Notes
