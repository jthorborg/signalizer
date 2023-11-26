# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## 0.4.1 - 2023-11-26

### Fixed

- Issue where loading old Signalizer presets (or opening old projects) would not connect the input to the output in the graph editor
- Crash when sustaining Signalizer on a non-playing track in the Oscilloscope when loading an old preset or project

## 0.4.0 - 2023-11-24

### Added

- Graph editor window that can be opened through a new global button in the upper right. From here you can sidechain any Signalizer into another, freely routing channels to at most 16 inputs. This configuration is saved in the project and presets as well.
- All views can now visualize up 8 channel pairs (16 channels in total), colours beyond the first pair are automatically distinct but based on the primary pair
- Signalizer will now utilize more CPU cores when processing more than one stereo pair
- VST3 support and builds on Windows
- Ask to prune the exception log file if it is becoming large
- Option to enable Legends in the views displaying the source name and colour used to draw it
- A control for changing what channel the Oscilloscope triggers on

### Fixed

- Crash in Oscilloscope when running in mono
- Spurious crash in Oscilloscope on startup
- Spurious crash when channel layout changes
- Class of hangs or crashes complaining about "Acquired non-recursed mutex had non-zero ref count!"
- Class of hangs or crashes complaining about "Unsafe writer created, while reader exists!"
- Bug where stereo filters in Vectorscope would not update when "Peak" mode was selected for envelopes
- Bug where vector scope channels would go out of sync and looks noisy
- Bug where graphical glitches would happen when switching to spectrogram mode
- Decay rate in the spectrogram that was incorrectly affected by the frame rate
- Issue where discontinuities would be present in the spectrogram under heavy usage
- Issue with latency in the spectrogram when frame smoothing was enabled
- Issue where the spectrum would produce a slightly out of data view
- Flicker when reducing the window size in the spectrum
- Graphical issue when the Spectrum has a dB range of zero (there's now a small minimum enforced)
- Graphical stutter or lag in hosts on Windows 10 and 11
- Issue where wireframe in the vectorscope would be rendered with jagged lines

### Changed

- Polar mode in Vectorscope fills the whole screen instead of being incorrectly scaled (going forward)
- Optimized wireframe drawing code
- Further optimized Vectorscope DSP code 
- Constant angular resolution for arcs
- Spectrum is now computed in 32-bit instead of 64-bit
- All frequency domain math is now considerably faster
- Dramatically reduced drawcalls using batching instead, vastly improving performance and responsiveness especially on less capable hardware
- Ultimately hide tabs after 7 seconds, hiding a problem where they would not auto-hide anyway
- Reduced scale of scrolling when altering gain in multichannel viewing mode of the Oscilloscope
- Improved treshold offset for trigger detection in spectral mode when analysing harmonic but sinusoidal signals

### Removed

- Hi-res timer. Use V-Sync for precision instead.

## 0.3.8 - 2023-10-29

## 0.3.7 - 2023-09-17

## 0.3.5 - 2022-11-28

## 0.3.4 - 2021-11-15

## 0.3.3 - 2019-09-22

### Added

- Unity build option, you now add `Source/Unity/SignalizerSource.cpp` to build all symbols necessary for embedding Oscilloscope, Vectorscope or Spectrum into another project
- CPL is now a direct submodule instead of parallel checkout
- `prepare.py` for fetching submodules and setting up the dev environment

### Changed

- Visual studio updated to 2019.
- Xcode updated 10.2


## 0.3.2 (alpha) - 2017-07-12

### Added

- New+more presets, scaled down antialiasing defaults
- Smoothing for the values in the frequency tracker in the spectrum, adjustable on the last page
- Audio processing and 2D graphics run in parallel and vastly improves performance in certain cases

### Fixed

- bad recall of the vectorscope's window size parameter

### Changed

- Minimum VS dev environment raised to VS 2017
- Minimum XCode dev environment raised to 7.3, implies El Capitan


## 0.3.1 (alpha) - 2017-04-13

### Added

- Option for auto-hiding widgets (frequency trackers etc.) when the mouse leaves the window
- Double-clicking right mouse button switches freezing mode permanently, a single right click still only toggles freeze mode while holding down the button

### Fixed

- Colour knobs initialize properly to a preset's value (they would previously 'stick')
- Fixed speed and update timings in certain cases in the spectrum, if diagnostics isn't enabled
- RMS auto gain on vectorscope adjusted (it was errournously adjusting the gain too high)


## 0.3.0 (alpha) - 2017-01-08

### Added

- Oscilloscope
- Main settings now contain options for disabling/enabling async audio processing, so you can choose between saving CPU or having all views in sync in freezed mode

### Changed

- Auto gain modes now work disjoint from the input gain, so you can adjust a gain offset while having auto gain

### Removed

- Repo stripped of binaries, refs updated; reclone it


## 0.2.10 (alpha) - 2017-01-08

### Added

- Experimental Linux support

### Changed

- Minimum compilation environment on OS X raised to XCode 7.2

### Fixed

- Possible crash on OS X for machines with older processors
- Fixed possible missing labels on toggled buttons
- Compiles on Ubuntu 16.04 LTS using Code::Blocks and distributed GCC, check repo for more details


## 0.2.9 (alpha) - 2016-08-28

### Changed

- Abbreviated some controls' names that couldn't fit

### Fixed

- A bug when using the "Compare" function in Logic Pro X
- Buttons not having correct state adopted from parameters
- Control edit space titles on colour controls when controlling a parameter


## 0.2.8 (alpha) - 2016-08-24

### Added

- Control edit spaces' title is now the control's associated automation parameter name
- All parameters are now exported and visible, and can be automated
- Spectrum no longer confined to a minimum of 3 dBs
- Spectrum window size no longer confined to multiples of 8, and can be less than 16 now
- Optional auto-hide of tabs after a second (makes Signalizer completely borderless with no unneeded visual content)
- Updated presets and added one for a simple polar mode

### Changed

- Old presets *are* supported and will continue to be, but they don't take effect until you open the tabs for the view in question.

### Fixed

- Correct and threaded serialization (with the new parameter system as well) for get state/set state, that can now be called safely on background threads. This was broken on multithreaded hosts
- Ordering of names of RGBA channels in the edit space for colour controls are now in correct order
- Slope values not being saved & restored correctly
- Flat Top windows only working with periodic shapes
- Versioning on serialized states, so a "preset" can now correctly contain parts from different versions
- Bugs with window size being set to zero in almost all circumstances
- Correct conversions between time and samples in window sizes
- Assumed unit on window size interpreters (Vectorscope would assume samples, even though it's in time)
- The DB meter graph to be completely correct and floating now; also supports negative dBs (needed for free-floating automation)
- Roundings of colour conversions from 8-bit to 64-bit
- Removal of a content component when clicking on an icon in the vertical tab
- Spacing in the Vectorscope so it isn't off by 1/N
- Scanning of peaks in the spectrum to be inclusive of the last considered element (nyquist was never considered otherwise)
- Interpretation of peak scanning to contain less NaN's in case of complete zero response
- FFT transforms reporting nyquist bin to be +3 dB.
- Frequency graphs for complex channel configurations in linear and logarithmic modes
- Drawing of lines for the grid should now only disable if alpha channel in the colour is zero
- Most of black lines/segments in flood fills
- Flood fills is now less taxing on graphics cards
- Not storing versioning in some preset storages (fixed some parameters not being saved)
- Fix for crash when mouse reaches edge of window in some cases
- Possible crash when changing sample rates
- Not redrawing spectrum graph on sample rate changes
- Frequency tracker being wrongly offset relative to the cursor in its normal mode
- Frequency graph not being recalculated after changing between complex and linear channel configurations
- Update smoothing not having any effects
- Audio history size (in the global settings) not being saved nor recalled
- Changing audio history size now correctly truncates the window sizes in the views, and they're correctly restored on resets
- Scrolling on OS X in views to support mouses AND trackpads while holding shift
- A regression that caused stereo filters on the vectorscope to use the filtering coefficient from the balance filters
- Stereo filters not being updated when auto-gain mode is set to none.


## 0.2.7 (alpha) - 2016-05-25

### Added

- Ability to adjust reference tuning that controls conversions to/from musical notes
- Variable slopes to the graph and spectrum views (to implement pink-noise scaling and such)
- Ability to default reset controls and widgets by alt-clicking them. This will reset them to what the previous loaded preset was.
- Windows Vista & XP support

### Changed

- Halved cpu-usage on audio thread, depending on sample buffer sizes
- The analyser box can now be moved while freezing the spectrum view by rightclicking
- Analysis box now renders in monospace font, and text layout has been reworked so it's much easier to read
- Analysis box (the peaktracker) now displays seminotes and cents of both the cursor and the tracked frequency

### Fixed

- Crash on resizing while using the peak tracker
- OpenGL crash on switching views
- Complete auval pass now (fix for incompatible version headers)
- x64 binaries not being recognized and/or crashing on Windows
- DB meter graph in the spectrum view being rendered half-way offscreen (luckily noone has noticed this)
- Crashes in relation to switching viewing modes and channel configuration modes in the spectrum (this should also fix occasional drop-outs in the spectrogram)
- Colouring of textboxes (and their names) in the colour control editor space


## 0.2.6 (alpha) - 2016-05-21

### Changed

- Better debug output

### Fixed

- PDB files


## 0.2.5 (alpha) - 2016-05-21

### Added

- Initial release.