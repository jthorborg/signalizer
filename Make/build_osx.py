import io
import configparser as cp
import os
import sys
import shutil as sh
import zipfile as zip
import subprocess
import common as cm

from datetime import date

def compiler_invoke(scheme, vstring, reloutdir):
	command = (
			   "xcodebuild "
			   "-project ../builds/macosx/signalizer.xcodeproj "
			   "clean"
			   )
	# new build system doesn't work with cleaning, hopefully you won't need it
	#if os.system(command) != 0:
	#	return -1
	
	command = (
			   "xcodebuild "
			   "-project ../builds/macosx/signalizer.xcodeproj "
			   "-scheme \"Signalizer - " + scheme + "\" "
			   "-configuration Release "
			   "CONFIGURATION_BUILD_DIR=\"" + cm.join(os.getcwd(), reloutdir) + "/\" "
			   # Following optional line removes nearly all symbol info, so it creates smaller packages but not really that great for debugging.
			   #"DEPLOYMENT_POSTPROCESSING=YES "
			   #"STRIP_INSTALLED_PRODUCT=YES "
			   #"SEPARATE_STRIP=YES "
			   #"COPY_PHASE_STRIP=YES "
			   "ONLY_ACTIVE_ARCH=NO "
			   "DYLIB_CURRENT_VERSION=" + vstring + " "
			   )
	print("---------> Compiler invocation: \n" + command)
	return os.system(command)

def set_plist_option(rel_plist_path, command):
	full_path = cm.join(os.getcwd(), rel_plist_path)
	os.system('/usr/libexec/PlistBuddy -c "' + command + '" "' + full_path + '"')

# parse config
config = cp.ConfigParser()
config.read("config.ini")

parameters = []
skip_vst2 = False
# handle cmd arguments
if len(sys.argv) > 1:
	for arg in sys.argv[1:]:
		inc = arg.find("-inc:")
		if inc != -1:
			parameters.append(arg[inc + 5:])
		else:
			skip_vst2 = skip_vst2 or arg.find("-skipvst2")

flush_parameters = False

# handle operations
for param in parameters:
	config.set("version", param, str(int(config.get("version", param)) + 1))
	print("------> Increasing " + param + " to " + config.get("version", param))

# write new configuration?
if len(parameters) > 0:
	flush_parameters = True

#configurations
major = config.get("version", "major")
minor = config.get("version", "minor")
build = config.get("version", "build")
desc = config.get("info", "description")
name = config.get("info", "productname")
company = config.get("info", "company")
author = config.get("info", "author")
manu4 = config.get("info", "manu4")
sub4 = config.get("info", "sub4")
version_string = major + "." + minor + "." + build
version_int = (int(major) << 48) | (int(minor) << 32) | int(build)
zipoutput = "../Releases/Signalizer macOS " + version_string

#diagnostic
print("------> Cleaning prior builds... ")

if os.system("xattr -w com.apple.xcode.CreatedByBuildSystem true ../Builds/MacOSX/build") != 0:
    print("------> Failed changing xattr for cleaning ...")
    exit(-2)
    
if os.system("xcodebuild -project ../builds/macosx/signalizer.xcodeproj -scheme \"Signalizer - All\" clean ONLY_ACTIVE_ARCH=NO") != 0:
    print("------> Failed cleaning...")
    exit(-3)

# xcodebuild -project ../builds/macosx/signalizer.xcodeproj -scheme "Signalizer - AU" -configuration Release clean archive ONLY_ACTIVE_ARCH=NO

print("------> Building Signalizer v. " + version_string + " release targets (" + str(version_int))

# rewrite program internal version

cm.rewrite_version_header("../Source/version.h", major, minor, build)

# rewrite build plist
root_plist = cm.join("../Builds/MacOSX/Info")
plist_variants = [root_plist + list + ".plist" for list in ["-AU", "-VST", "-VST3", "-VST3_Manifest_Helper", ""]]

for plist in plist_variants:
    print("------> rewriting plist " + plist)
    set_plist_option(plist, "Set :CFBundleIdentifier com." + company + "." + name)
    set_plist_option(plist, "Set :CFBundleShortVersionString " + version_string)
    set_plist_option(plist, "Set :CFBundleVersion " + version_string)
    set_plist_option(plist, "Set :NSHumanReadableCopyright Copyright (c) " + str(date.today().year) + " " + author)

# set the audio unit plugin description
aulist = plist_variants[0]
set_plist_option(aulist, "Set :AudioComponents:0:description " + desc)
set_plist_option(aulist, "Set :AudioComponents:0:manufacturer " + manu4)
set_plist_option(aulist, "Set :AudioComponents:0:name " + company + ": " + name)
set_plist_option(aulist, "Set :AudioComponents:0:subtype " + sub4)
set_plist_option(aulist, "Set :AudioComponents:0:type aufx")
set_plist_option(aulist, "Set :AudioComponents:0:version " + str(version_int))

build_variants = ["AU", "VST3"]

if not skip_vst2:
    build_variants = build_variants + ["VST"]
    
print(build_variants)

build_folder = "Signalizer macOS Universal"
if os.path.exists(build_folder):
	sh.rmtree(build_folder)

#run all targets
for build in build_variants:
	if compiler_invoke(build, version_string, build_folder) != 0:
		print("\n------> Error compiling for target " + build)
		sh.rmtree(build_folder)
		exit(-4)

#remove temporary junk
if os.path.exists(cm.join(build_folder, "libSignalizer.a")) != 0:
	os.remove(cm.join(build_folder, "libSignalizer.a"))

if os.path.exists(cm.join(build_folder, "juce_vst3_helper")) != 0:
	os.remove(cm.join(build_folder, "juce_vst3_helper"))

#append extra goodies
cm.create_build_file(cm.join(build_folder, "Build.log"), version_string)
sh.copyfile("../CHANGELOG.md", cm.join(build_folder, "CHANGELOG.md"))
sh.copyfile("macos_installation_advice.txt", cm.join(build_folder, "HOW TO INSTALL.txt"))

print("------> Zipping output directories...")

zx = sh.make_archive(zipoutput, "zip", build_folder)

print("------> Build Signalizer successfully to:")
print("------> " + zx)

# clean up dirs
sh.rmtree(build_folder)

# done, if we made it here, increase the conf build

if flush_parameters:
	with open("config.ini", "w") as f:
		config.write(f)
