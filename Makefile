# Signalizer Makefile
# Python-free build system for macOS

# Configuration (from config.ini)
VERSION_MAJOR := 0
VERSION_MINOR := 4
VERSION_BUILD := 3
VERSION_STRING := $(VERSION_MAJOR).$(VERSION_MINOR).$(VERSION_BUILD)
VERSION_INT := $(shell echo $$(($(VERSION_MAJOR) << 48 | $(VERSION_MINOR) << 32 | $(VERSION_BUILD))))

PRODUCT_NAME := Signalizer
COMPANY := Lightbridge
AUTHOR := Janus Lynggaard Thorborg
MANU4 := LbJt
SUB4 := Sign
DESCRIPTION := Real-time audio visualization plugin

# Build configuration
ARCH := arm64
BUILD_FOLDER := Signalizer_OSX
BUILD_DIR := $(BUILD_FOLDER)/$(ARCH)
RELEASE_DIR := Releases
ZIP_OUTPUT := $(RELEASE_DIR)/Signalizer\ OS\ X\ $(VERSION_STRING)
XCODE_PROJECT := Builds/MacOSX/Signalizer.xcodeproj
PLIST_FILE := Builds/MacOSX/Info.plist

# Architecture-specific library paths
ifeq ($(ARCH),x86_64)
    LIB_INCLUDE_PATH := /usr/local/include
    LIB_LIBRARY_PATH := /usr/local/lib
    LIB_PNG_PATH := /usr/local/lib/libpng.dylib
else ifeq ($(ARCH),arm64)
    LIB_INCLUDE_PATH := /opt/homebrew/include  
    LIB_LIBRARY_PATH := /opt/homebrew/lib
    LIB_PNG_PATH := /opt/homebrew/lib/libpng.dylib
endif

# Get build info
BUILD_TIME := $(shell date -u "+%Y-%m-%d %H:%M:%S")
BUILD_USER := $(shell whoami)
BUILD_SYSTEM := $(shell uname -s)
BUILD_RELEASE := $(shell uname -r)
GIT_BRANCH := $(shell git branch --show-current 2>/dev/null || echo "unknown")
GIT_COMMIT := $(shell git describe --always 2>/dev/null || echo "unknown")

.PHONY: all setup build build-debug clean increment-major increment-minor increment-patch install-macos-arm64 install-debug help build-projucer regenerate-project check-deps-x86 check-deps-arm64 build-x86 build-internal build-debug-x86 build-debug-arm64

all: setup regenerate-project build

help:
	@echo "Signalizer Build System"
	@echo ""
	@echo "Targets:"
	@echo "  all              - Setup and build (default)"
	@echo "  setup            - Initialize submodules"
	@echo "  build-projucer   - Build Projucer tool"
	@echo "  regenerate-project - Regenerate Xcode project from .jucer file"
	@echo "  build            - Build Signalizer for ARM64 (default)"
	@echo "  build-debug      - Build Signalizer with debug symbols for ARM64 (default)"
	@echo "  build-x86        - Build Signalizer for x86_64"
	@echo "  build-debug-x86  - Build Signalizer with debug symbols for x86_64"
	@echo "  check-deps-x86   - Check x86_64 dependencies"
	@echo "  check-deps-arm64 - Check ARM64 dependencies"
	@echo "  clean            - Clean build artifacts"
	@echo "  increment-major  - Increment major version and build"
	@echo "  increment-minor  - Increment minor version and build"
	@echo "  increment-patch  - Increment patch version and build"
	@echo "  install          - Install built plugins to user library"
	@echo "  install-debug    - Install debug VST3 with symbols for Instruments profiling"
	@echo ""
	@echo "Current version: $(VERSION_STRING)"
	@echo "Current architecture: $(ARCH)"
	@echo "Library paths: $(LIB_INCLUDE_PATH), $(LIB_LIBRARY_PATH)"

setup:
	@echo ">> Updating submodules..."
	git submodule update --init --recursive
	git submodule update
	@echo ">> Dev environment setup without errors."

build-projucer:
	@echo ">> Building Projucer..."
	@if [ ! -f "External/juce/extras/Projucer/Builds/MacOSX/build/Release/Projucer.app/Contents/MacOS/Projucer" ]; then \
		cd External/juce/extras/Projucer/Builds/MacOSX && \
		xcodebuild -project Projucer.xcodeproj -scheme "Projucer - App" -configuration Release; \
	fi

regenerate-project: build-projucer
	@echo ">> Regenerating Xcode project from .jucer file..."
	External/juce/extras/Projucer/Builds/MacOSX/build/Release/Projucer.app/Contents/MacOS/Projucer --resave Signalizer.jucer

check-deps-x86:
	@echo ">> Checking x86_64 dependencies..."
	@if [ ! -f "/usr/local/lib/libpng.dylib" ]; then \
		echo "❌ ERROR: x86_64 libpng not found at /usr/local/lib/libpng.dylib"; \
		echo ""; \
		echo "To install x86_64 Homebrew and libpng:"; \
		echo "1. Open Terminal with Rosetta:"; \
		echo "   arch -x86_64 /usr/bin/env bash"; \
		echo "2. Install x86_64 Homebrew:"; \
		echo "   /bin/bash -c \"\$$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)\""; \
		echo "3. Install libpng:"; \
		echo "   /usr/local/bin/brew install libpng"; \
		echo ""; \
		exit 1; \
	else \
		echo "✅ x86_64 libpng found at /usr/local/lib/libpng.dylib"; \
		lipo -info /usr/local/lib/libpng.dylib; \
	fi

check-deps-arm64:
	@echo ">> Checking ARM64 dependencies..."
	@if [ ! -f "/opt/homebrew/lib/libpng.dylib" ]; then \
		echo "❌ ERROR: ARM64 libpng not found at /opt/homebrew/lib/libpng.dylib"; \
		echo "Install with: brew install libpng"; \
		exit 1; \
	else \
		echo "✅ ARM64 libpng found at /opt/homebrew/lib/libpng.dylib"; \
		lipo -info /opt/homebrew/lib/libpng.dylib; \
	fi

build: setup check-deps-arm64
	@echo "------> Building Signalizer v. $(VERSION_STRING) release targets ($(VERSION_INT)) for arm64"
	@$(MAKE) build-internal ARCH=arm64 BUILD_DIR=$(BUILD_FOLDER)/arm64

build-x86: setup check-deps-x86
	@echo "------> Building Signalizer v. $(VERSION_STRING) release targets ($(VERSION_INT)) for x86_64"
	@$(MAKE) build-internal ARCH=x86_64 BUILD_DIR=$(BUILD_FOLDER)/x64

build-internal:
	@mkdir -p $(BUILD_DIR)
	@mkdir -p $(RELEASE_DIR)
	
	# Update version header
	@echo "Updating version header..."
	@echo '#define SIGNALIZER_MAJOR $(VERSION_MAJOR)' > Source/version.h
	@echo '#define SIGNALIZER_MINOR $(VERSION_MINOR)' >> Source/version.h
	@echo '#define SIGNALIZER_BUILD $(VERSION_BUILD)' >> Source/version.h
	@echo '#define SIGNALIZER_BUILD_INFO "$(BUILD_TIME): $(PRODUCT_NAME) $(VERSION_STRING) built on $(BUILD_SYSTEM) $(BUILD_RELEASE) by $(BUILD_USER)\\n$(GIT_BRANCH)\\n$(GIT_COMMIT)"' >> Source/version.h
	@echo '#define SIGNALIZER_VERSION_STRING "$(VERSION_STRING)"' >> Source/version.h
	@echo '#define SIGNALIZER_VST_VERSION_HEX 0x$(shell printf "%02x%02x%02x" $$(( $(VERSION_MAJOR) % 255 )) $$(( $(VERSION_MINOR) % 255 )) $$(( $(VERSION_BUILD) % 255 )))' >> Source/version.h
	
	# Update plist
	@echo "Updating plist..."
	/usr/libexec/PlistBuddy -c "Set :CFBundleIdentifier com.$(COMPANY).$(PRODUCT_NAME)" "$(PLIST_FILE)"
	/usr/libexec/PlistBuddy -c "Set :CFBundleShortVersionString $(VERSION_STRING)" "$(PLIST_FILE)"
	/usr/libexec/PlistBuddy -c "Set :CFBundleVersion $(VERSION_STRING)" "$(PLIST_FILE)"
	/usr/libexec/PlistBuddy -c "Set :NSHumanReadableCopyright Copyright (c) $(shell date +%Y) $(AUTHOR)" "$(PLIST_FILE)"
	/usr/libexec/PlistBuddy -c "Set :AudioComponents:0:description $(DESCRIPTION)" "$(PLIST_FILE)"
	/usr/libexec/PlistBuddy -c "Set :AudioComponents:0:manufacturer $(MANU4)" "$(PLIST_FILE)"
	/usr/libexec/PlistBuddy -c "Set :AudioComponents:0:name $(COMPANY): $(PRODUCT_NAME)" "$(PLIST_FILE)"
	/usr/libexec/PlistBuddy -c "Set :AudioComponents:0:subtype $(SUB4)" "$(PLIST_FILE)"
	/usr/libexec/PlistBuddy -c "Set :AudioComponents:0:type aufx" "$(PLIST_FILE)"
	/usr/libexec/PlistBuddy -c "Set :AudioComponents:0:version $(VERSION_INT)" "$(PLIST_FILE)"
	
	# Build with Xcode
	@echo "---------> Compiler invocation for $(ARCH):"
	xcodebuild \
		-project $(XCODE_PROJECT) \
		-scheme "Signalizer - All" \
		-configuration Release \
		CONFIGURATION_BUILD_DIR=$(shell pwd)/$(BUILD_DIR)/ \
		STRIP_INSTALLED_PRODUCT=YES \
		SEPARATE_STRIP=YES \
		COPY_PHASE_STRIP=YES \
		ARCHS=$(ARCH) \
		ONLY_ACTIVE_ARCH=NO \
		DYLIB_CURRENT_VERSION=$(VERSION_STRING) \
		GCC_FAST_MATH=NO \
		ENABLE_STRICT_OBJC_MSGSEND=NO \
		OTHER_CFLAGS="-DJUCE_INCLUDE_PNGLIB_CODE=0 -I$(LIB_INCLUDE_PATH) -Wno-writable-strings -DJUCE_SILENCE_XCODE_15_LINKER_WARNING=1 -DDONT_SET_USING_JUCE_NAMESPACE=1 -Wno-error -fno-lto" \
		OTHER_LDFLAGS="-L$(LIB_LIBRARY_PATH) -lpng -lz -lSignalizer -fno-lto -Wl,-v" \
		VERBOSE=1 \
		CLANG_ENABLE_OBJC_WEAK=NO \
		-destination "platform=macOS"
	
	# Create build log
	@mkdir -p $(BUILD_DIR)/Signalizer.component/Contents/Resources
	@echo "$(BUILD_TIME): Signalizer $(VERSION_STRING) built on $(BUILD_SYSTEM) $(BUILD_RELEASE) by $(BUILD_USER)" > $(BUILD_DIR)/Signalizer.component/Contents/Resources/Build.log
	@echo "$(GIT_BRANCH)" >> $(BUILD_DIR)/Signalizer.component/Contents/Resources/Build.log
	@echo "$(GIT_COMMIT)" >> $(BUILD_DIR)/Signalizer.component/Contents/Resources/Build.log
	@echo "" >> $(BUILD_DIR)/Signalizer.component/Contents/Resources/Build.log
	@git log -5 >> $(BUILD_DIR)/Signalizer.component/Contents/Resources/Build.log 2>/dev/null || true
	
	# Copy changelog and skeleton resources
	@cp CHANGELOG.md $(BUILD_DIR)/Signalizer.component/Contents/Resources/ || exit 1
	@cp -R Make/Skeleton/* $(BUILD_DIR)/Signalizer.component/Contents/Resources/ || exit 1
	# Re-sign after adding resources
	@codesign --force --sign - $(BUILD_DIR)/Signalizer.component
	
	@echo ""
	@echo "------> All builds finished, generating plugin permutations ..."
	
	# Create plugin variants
	@cp -R $(BUILD_DIR)/Signalizer.component $(BUILD_DIR)/Signalizer.vst
	# Re-sign VST copy
	@codesign --force --sign - $(BUILD_DIR)/Signalizer.vst
	# VST3 is built separately by Xcode, so we need to add resources to it
	@cp -R Make/Skeleton/* $(BUILD_DIR)/Signalizer.vst3/Contents/Resources/ || exit 1
	@cp CHANGELOG.md $(BUILD_DIR)/Signalizer.vst3/Contents/Resources/ || exit 1
	@cp $(BUILD_DIR)/Signalizer.component/Contents/Resources/Build.log $(BUILD_DIR)/Signalizer.vst3/Contents/Resources/ || exit 1
	# Re-sign after adding resources
	@codesign --force --sign - $(BUILD_DIR)/Signalizer.vst3
	
	@echo "------> Zipping output directories..."
	
	# Copy installation instructions and create zip
	@cp Make/macos_installation_advice.txt "$(BUILD_FOLDER)/HOW TO INSTALL.txt"
	@cd $(BUILD_FOLDER) && zip -r "../$(ZIP_OUTPUT).zip" .
	
	@echo "------> Built Signalizer successfully into:"
	@echo "------> $(ZIP_OUTPUT).zip"

# Debug build targets
build-debug: setup check-deps-arm64
	$(MAKE) build-debug-arm64

build-debug-x86:
	$(MAKE) ARCH=x86_64 BUILD_DIR=$(BUILD_FOLDER)/x64 build-debug-internal

build-debug-arm64: setup check-deps-arm64
	$(MAKE) ARCH=arm64 BUILD_DIR=$(BUILD_FOLDER)/arm64 build-debug-internal

build-debug-internal:
	@echo "========================================="
	@echo "Building Signalizer DEBUG for $(ARCH)..."
	@echo "========================================="
	@echo ""
	@echo "---------> Creating version.h with debugging info for $(ARCH):"
	@mkdir -p Source
	@echo "#define SIGNALIZER_MAJOR $(VERSION_MAJOR)" > Source/version.h
	@echo "#define SIGNALIZER_MINOR $(VERSION_MINOR)" >> Source/version.h
	@echo "#define SIGNALIZER_BUILD $(VERSION_BUILD)" >> Source/version.h
	@echo "#define SIGNALIZER_BUILD_INFO \"$(BUILD_TIME) $(BUILD_USER)@$(BUILD_SYSTEM) $(BUILD_RELEASE) [$(GIT_BRANCH):$(GIT_COMMIT)]\"" >> Source/version.h
	@cat Source/version.h
	@echo ""
	@echo "---------> Compiler invocation for $(ARCH) DEBUG:"
	xcodebuild \
		-project $(XCODE_PROJECT) \
		-scheme "Signalizer - All" \
		-configuration Debug \
		CONFIGURATION_BUILD_DIR=$(shell pwd)/$(BUILD_DIR)/ \
		STRIP_INSTALLED_PRODUCT=NO \
		SEPARATE_STRIP=NO \
		GCC_GENERATE_DEBUGGING_SYMBOLS=YES \
		DEBUG_INFORMATION_FORMAT=dwarf-with-dsym \
		GCC_OPTIMIZATION_LEVEL=0 \
		SWIFT_OPTIMIZATION_LEVEL=-Onone \
		COPY_PHASE_STRIP=NO \
		ARCHS=$(ARCH) \
		VALID_ARCHS=$(ARCH) \
		ENABLE_STRICT_OBJC_MSGSEND=NO \
		OTHER_CFLAGS="-g -O0 -DDEBUG=1 -DJUCE_INCLUDE_PNGLIB_CODE=0 -I$(LIB_INCLUDE_PATH) -Wno-writable-strings -DJUCE_SILENCE_XCODE_15_LINKER_WARNING=1 -DDONT_SET_USING_JUCE_NAMESPACE=1 -Wno-error" \
		OTHER_LDFLAGS="-L$(LIB_LIBRARY_PATH) -lpng -lz -lSignalizer $(LIB_PNG_PATH)" \
		VERBOSE=1 \
		CLANG_ENABLE_OBJC_WEAK=NO
	
	@echo ""
	
	# Copy resources and build metadata (same as release build but with debug symbols)
	@echo "---------> Preparing debug build directories ..."
	@mkdir -p $(BUILD_DIR)/Signalizer.component/Contents/Resources
	@echo "Build: $(BUILD_TIME)" > $(BUILD_DIR)/Signalizer.component/Contents/Resources/Build.log
	@echo "User: $(BUILD_USER)" >> $(BUILD_DIR)/Signalizer.component/Contents/Resources/Build.log
	@echo "System: $(BUILD_SYSTEM) $(BUILD_RELEASE)" >> $(BUILD_DIR)/Signalizer.component/Contents/Resources/Build.log
	@echo "Branch: $(GIT_BRANCH)" >> $(BUILD_DIR)/Signalizer.component/Contents/Resources/Build.log
	@echo "Commit: " >> $(BUILD_DIR)/Signalizer.component/Contents/Resources/Build.log
	@echo "$(GIT_COMMIT)" >> $(BUILD_DIR)/Signalizer.component/Contents/Resources/Build.log
	@echo "" >> $(BUILD_DIR)/Signalizer.component/Contents/Resources/Build.log
	@git log -5 >> $(BUILD_DIR)/Signalizer.component/Contents/Resources/Build.log 2>/dev/null || true
	
	# Copy changelog and skeleton resources
	@cp CHANGELOG.md $(BUILD_DIR)/Signalizer.component/Contents/Resources/ || exit 1
	@cp -R Make/Skeleton/* $(BUILD_DIR)/Signalizer.component/Contents/Resources/ || exit 1
	# Re-sign after adding resources
	@codesign --force --sign - $(BUILD_DIR)/Signalizer.component
	
	@echo ""
	@echo "------> All debug builds finished, generating plugin permutations ..."
	
	# Create plugin variants
	@cp -R $(BUILD_DIR)/Signalizer.component $(BUILD_DIR)/Signalizer.vst
	# Re-sign VST copy
	@codesign --force --sign - $(BUILD_DIR)/Signalizer.vst
	# VST3 is built separately by Xcode, so we need to add resources to it
	@cp -R Make/Skeleton/* $(BUILD_DIR)/Signalizer.vst3/Contents/Resources/ || exit 1
	@cp CHANGELOG.md $(BUILD_DIR)/Signalizer.vst3/Contents/Resources/ || exit 1
	@cp $(BUILD_DIR)/Signalizer.component/Contents/Resources/Build.log $(BUILD_DIR)/Signalizer.vst3/Contents/Resources/ || exit 1
	# Re-sign after adding resources
	@codesign --force --sign - $(BUILD_DIR)/Signalizer.vst3
	
	@echo "------> Built Signalizer DEBUG successfully into:"
	@echo "------> $(BUILD_FOLDER)/"

clean:
	@echo "Cleaning build artifacts..."
	@rm -rf $(BUILD_FOLDER)
	@rm -f Source/version.h

install-macos-arm64:
	@echo "Installing plugins to user library..."
# 	@cp -R $(BUILD_DIR)/Signalizer.component /Library/Audio/Plug-Ins/Components/
# 	@cp -R $(BUILD_DIR)/Signalizer.vst /Library/Audio/Plug-Ins/VST/
	rm -rf /Library/Audio/Plug-Ins/VST3/Signalizer.vst3
	cp -R Signalizer_OSX/arm64/Signalizer.vst3 /Library/Audio/Plug-Ins/VST3/
	xattr -rc /Library/Audio/Plug-Ins/VST3/Signalizer.vst3
	echo "Installation complete. You may need to restart your DAW."

install-debug:
	@echo "Installing debug VST3 plugin for Instruments profiling..."
	@if [ ! -d "$(BUILD_DIR)/Signalizer.vst3" ]; then \
		echo "Error: Debug build not found at $(BUILD_DIR)/Signalizer.vst3"; \
		echo "Run 'make build-debug' first."; \
		exit 1; \
	fi
	@rm -rf /Library/Audio/Plug-Ins/VST3/Signalizer.vst3
	@rm -rf /Library/Audio/Plug-Ins/VST3/Signalizer.vst3.dSYM
	@cp -R $(BUILD_DIR)/Signalizer.vst3 /Library/Audio/Plug-Ins/VST3/
	@if [ -d "$(BUILD_DIR)/Signalizer.vst3.dSYM" ]; then \
		cp -R $(BUILD_DIR)/Signalizer.vst3.dSYM /Library/Audio/Plug-Ins/VST3/; \
		echo "✅ Installed debug symbols: /Library/Audio/Plug-Ins/VST3/Signalizer.vst3.dSYM"; \
	fi
	@xattr -rc /Library/Audio/Plug-Ins/VST3/Signalizer.vst3
	@echo "✅ Installed debug plugin: /Library/Audio/Plug-Ins/VST3/Signalizer.vst3"
	@echo "✅ Architecture: $(ARCH)"
	@echo "✅ Build directory: $(BUILD_DIR)"
	@echo ""
	@echo "The debug plugin is now ready for Instruments profiling."
	@echo "Restart your DAW to load the new version."

increment-major:
	@echo "Incrementing major version..."
	@$(eval VERSION_MAJOR := $(shell echo $$(($(VERSION_MAJOR) + 1))))
	@sed -i '' 's/^major = .*/major = $(VERSION_MAJOR)/' Make/config.ini
	@$(MAKE) build

increment-minor:
	@echo "Incrementing minor version..."
	@$(eval VERSION_MINOR := $(shell echo $$(($(VERSION_MINOR) + 1))))
	@sed -i '' 's/^minor = .*/minor = $(VERSION_MINOR)/' Make/config.ini
	@$(MAKE) build

increment-patch:
	@echo "Incrementing patch version..."
	@$(eval VERSION_BUILD := $(shell echo $$(($(VERSION_BUILD) + 1))))
	@sed -i '' 's/^build = .*/build = $(VERSION_BUILD)/' Make/config.ini
	@$(MAKE) build

