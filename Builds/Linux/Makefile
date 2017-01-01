# Automatically generated makefile, created by the Introjucer
# Don't edit this file! Your changes will be overwritten when you re-save the Introjucer project!

# (this disables dependency generation if multiple architectures are set)
DEPFLAGS := $(if $(word 2, $(TARGET_ARCH)), , -MMD)

ifndef CONFIG
  CONFIG=Debug
endif

ifeq ($(CONFIG),Debug)
  BINDIR := build
  LIBDIR := build
  OBJDIR := build/intermediate/Debug
  OUTDIR := build

  ifeq ($(TARGET_ARCH),)
    TARGET_ARCH := -march=native
  endif

  CPPFLAGS := $(DEPFLAGS) -D "LINUX=1" -D "DEBUG=1" -D "_DEBUG=1" -D "JUCER_LINUX_MAKE_7346DA2A=1" -D "JUCE_APP_VERSION=0.2.10" -D "JUCE_APP_VERSION_HEX=0x20a" -I /usr/include -I /usr/include/freetype2 -I ../../JuceLibraryCode
  CFLAGS += $(CPPFLAGS) $(TARGET_ARCH) -g -ggdb -fPIC -O0 -flax-vector-conversions -std=c++14 -I ../../../SDKs/
  CXXFLAGS += $(CFLAGS)
  LDFLAGS += $(TARGET_ARCH) -L$(BINDIR) -L$(LIBDIR) -shared -L/usr/X11R6/lib/ 
  LDDEPS :=
  RESFLAGS :=  -D "LINUX=1" -D "DEBUG=1" -D "_DEBUG=1" -D "JUCER_LINUX_MAKE_7346DA2A=1" -D "JUCE_APP_VERSION=0.2.10" -D "JUCE_APP_VERSION_HEX=0x20a" -I /usr/include -I /usr/include/freetype2 -I ../../JuceLibraryCode
  TARGET := Signalizer
  BLDCMD = $(CXX) -o $(OUTDIR)/$(TARGET) $(OBJECTS) $(LDFLAGS) $(RESOURCES) $(TARGET_ARCH)
  CLEANCMD = rm -rf $(OUTDIR)/$(TARGET) $(OBJDIR)
endif

ifeq ($(CONFIG),Release)
  BINDIR := build
  LIBDIR := build
  OBJDIR := build/intermediate/Release
  OUTDIR := build

  ifeq ($(TARGET_ARCH),)
    TARGET_ARCH := -march=native
  endif

  CPPFLAGS := $(DEPFLAGS) -D "LINUX=1" -D "NDEBUG=1" -D "JUCER_LINUX_MAKE_7346DA2A=1" -D "JUCE_APP_VERSION=0.2.10" -D "JUCE_APP_VERSION_HEX=0x20a" -I /usr/include -I /usr/include/freetype2 -I ../../JuceLibraryCode
  CFLAGS += $(CPPFLAGS) $(TARGET_ARCH) -fPIC -flax-vector-conversions -Os -std=c++14 -I ../../../SDKs/ -fmax-errors=2
  CXXFLAGS += $(CFLAGS)
  LDFLAGS += $(TARGET_ARCH) -L$(BINDIR) -L$(LIBDIR) -shared -fvisibility=hidden -L/usr/X11R6/lib/ 
  LDDEPS :=
  RESFLAGS :=  -D "LINUX=1" -D "NDEBUG=1" -D "JUCER_LINUX_MAKE_7346DA2A=1" -D "JUCE_APP_VERSION=0.2.10" -D "JUCE_APP_VERSION_HEX=0x20a" -I /usr/include -I /usr/include/freetype2 -I ../../JuceLibraryCode
  TARGET := Signalizer
  BLDCMD = $(CXX) -o $(OUTDIR)/$(TARGET) $(OBJECTS) $(LDFLAGS) $(RESOURCES) $(TARGET_ARCH)
  CLEANCMD = rm -rf $(OUTDIR)/$(TARGET) $(OBJDIR)
endif

OBJECTS := \
  $(OBJDIR)/COscilloscope_274f3e73.o \
  $(OBJDIR)/COscilloscopeRendering_59b4dd23.o \
  $(OBJDIR)/CSpectrum_8f1434b5.o \
  $(OBJDIR)/CSpectrumDSP_aefdb1d6.o \
  $(OBJDIR)/CSpectrumRendering_fc13d121.o \
  $(OBJDIR)/CVectorScope_1ac4a3f7.o \
  $(OBJDIR)/CVectorScopeRendering_695f6c1f.o \
  $(OBJDIR)/MainEditor_7dd2070f.o \
  $(OBJDIR)/PluginProcessor_a059e380.o \
  $(OBJDIR)/SignalizerDesign_7672022d.o \

.PHONY: clean

$(OUTDIR)/$(TARGET): $(OBJECTS) $(LDDEPS) $(RESOURCES)
	@echo Linking Signalizer
	-@mkdir -p $(BINDIR)
	-@mkdir -p $(LIBDIR)
	-@mkdir -p $(OUTDIR)
	@$(BLDCMD)

clean:
	@echo Cleaning Signalizer
	@$(CLEANCMD)

strip:
	@echo Stripping Signalizer
	-@strip --strip-unneeded $(OUTDIR)/$(TARGET)

$(OBJDIR)/COscilloscope_274f3e73.o: ../../Source/COscilloscope.cpp
	-@mkdir -p $(OBJDIR)
	@echo "Compiling COscilloscope.cpp"
	@$(CXX) $(CXXFLAGS) -o "$@" -c "$<"

$(OBJDIR)/COscilloscopeRendering_59b4dd23.o: ../../Source/COscilloscopeRendering.cpp
	-@mkdir -p $(OBJDIR)
	@echo "Compiling COscilloscopeRendering.cpp"
	@$(CXX) $(CXXFLAGS) -o "$@" -c "$<"

$(OBJDIR)/CSpectrum_8f1434b5.o: ../../Source/CSpectrum.cpp
	-@mkdir -p $(OBJDIR)
	@echo "Compiling CSpectrum.cpp"
	@$(CXX) $(CXXFLAGS) -o "$@" -c "$<"

$(OBJDIR)/CSpectrumDSP_aefdb1d6.o: ../../Source/CSpectrumDSP.cpp
	-@mkdir -p $(OBJDIR)
	@echo "Compiling CSpectrumDSP.cpp"
	@$(CXX) $(CXXFLAGS) -o "$@" -c "$<"

$(OBJDIR)/CSpectrumRendering_fc13d121.o: ../../Source/CSpectrumRendering.cpp
	-@mkdir -p $(OBJDIR)
	@echo "Compiling CSpectrumRendering.cpp"
	@$(CXX) $(CXXFLAGS) -o "$@" -c "$<"

$(OBJDIR)/CVectorScope_1ac4a3f7.o: ../../Source/CVectorScope.cpp
	-@mkdir -p $(OBJDIR)
	@echo "Compiling CVectorScope.cpp"
	@$(CXX) $(CXXFLAGS) -o "$@" -c "$<"

$(OBJDIR)/CVectorScopeRendering_695f6c1f.o: ../../Source/CVectorScopeRendering.cpp
	-@mkdir -p $(OBJDIR)
	@echo "Compiling CVectorScopeRendering.cpp"
	@$(CXX) $(CXXFLAGS) -o "$@" -c "$<"

$(OBJDIR)/MainEditor_7dd2070f.o: ../../Source/MainEditor.cpp
	-@mkdir -p $(OBJDIR)
	@echo "Compiling MainEditor.cpp"
	@$(CXX) $(CXXFLAGS) -o "$@" -c "$<"

$(OBJDIR)/PluginProcessor_a059e380.o: ../../Source/PluginProcessor.cpp
	-@mkdir -p $(OBJDIR)
	@echo "Compiling PluginProcessor.cpp"
	@$(CXX) $(CXXFLAGS) -o "$@" -c "$<"

$(OBJDIR)/SignalizerDesign_7672022d.o: ../../Source/SignalizerDesign.cpp
	-@mkdir -p $(OBJDIR)
	@echo "Compiling SignalizerDesign.cpp"
	@$(CXX) $(CXXFLAGS) -o "$@" -c "$<"

-include $(OBJECTS:%.o=%.d)
