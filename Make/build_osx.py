import io
import ConfigParser as cp
import os
import sys
import shutil as sh
import zipfile as zip
import subprocess
import common as cm

from datetime import date

def compiler_invoke(arch, vstring, reloutdir):
	command = (
			   "xcodebuild "
			   "-project ../builds/macosx/signalizer.xcodeproj "
			   "clean"
			   )
	
	if os.system(command) != 0:
		return -1
	
	command = (
			   "xcodebuild "
			   "-project ../builds/macosx/signalizer.xcodeproj "
			   "-scheme Signalizer "
			   "-configuration Release "
			   "CONFIGURATION_BUILD_DIR=" + cm.join(os.getcwd(), reloutdir) + "/ "
			   # Following optional line removes nearly all symbol info, so it creates smaller packages but not really that great for debugging.
			   #"DEPLOYMENT_POSTPROCESSING=YES "
			   "STRIP_INSTALLED_PRODUCT=YES "
			   "SEPARATE_STRIP=YES "
			   "COPY_PHASE_STRIP=YES "
			   "ARCHS=" + arch + " "
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

# handle cmd arguments
if len(sys.argv) > 1:
	for arg in sys.argv[1:]:
		inc = arg.find("-inc:")
		if inc != -1:
			parameters.append(arg[inc + 5:])

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
build_folder = "Signalizer_OSX"
zipoutput = "../Releases/Signalizer OS X " + version_string
#diagnostic
print("------> Building Signalizer v. " + version_string + " release targets (" + str(version_int))

# [0] = arg to clang, [1] = output folder
targets = [["i386", cm.join(build_folder, "x32")], ["x86_64", cm.join(build_folder, "x64")]]

# rewrite program internal version

cm.rewrite_version_header("../Source/version.h", major, minor, build)

# rewrite build plist
plist = cm.join("../Builds/MacOSX/Info.plist")
set_plist_option(plist, "Set :CFBundleIdentifier com." + company + "." + name)
set_plist_option(plist, "Set :CFBundleShortVersionString " + version_string)
set_plist_option(plist, "Set :CFBundleVersion " + version_string)
set_plist_option(plist, "Set :NSHumanReadableCopyright Copyright (c) " + str(date.today().year) + " " + author)

# set the audio unit plugin description
set_plist_option(plist, "Set :AudioComponents:0:description " + desc)
set_plist_option(plist, "Set :AudioComponents:0:manufacturer " + manu4)
set_plist_option(plist, "Set :AudioComponents:0:name " + company + ": " + name)
set_plist_option(plist, "Set :AudioComponents:0:subtype " + sub4)
set_plist_option(plist, "Set :AudioComponents:0:type aufx")
set_plist_option(plist, "Set :AudioComponents:0:version " + str(version_int))


#run all targets
for option in targets:
	if compiler_invoke(option[0], version_string, option[1]) != 0:
		print("\n------> Error compiling for target " + option[0])
		sh.rmtree(build_folder)
		exit(1)
	else:
		cm.create_build_file(cm.join(cm.join(cm.join(cm.join(option[1], "Signalizer.component"), "Contents"), "Resources"), "Build.log"), version_string)


print("\n------> All builds finished, generating plugin permutations ...")

# make build permutations and set up the plist file

for option in targets:
	original = cm.join(option[1], "Signalizer.component")
	# permute
	sh.copytree(original, cm.join(option[1], "Signalizer.vst"))
	sh.copytree(original, cm.join(option[1], "Signalizer.vst3"))


print("------> Zipping output directories...")

zx = sh.make_archive(zipoutput, "zip", build_folder)

print("------> Builded Signalizer successfully into:")
print("------> " + zx)

# clean up dirs
sh.rmtree(build_folder)

# done, if we made it here, increase the conf build

if flush_parameters:
	with open("config.ini", "w") as f:
		config.write(f)
