import io
import ConfigParser as cp
import os
import sys
import shutil as sh
import zipfile as zip
from time import gmtime, strftime
import getpass
import platform
import subprocess
from datetime import date

def compiler_invoke(arch, vstring, reloutdir):
	command = (
			   "xcodebuild "
			   "-project ../builds/macosx/signalizer.xcodeproj "
			   "-scheme Signalizer "
			   "-configuration Release "
			   "CONFIGURATION_BUILD_DIR=" + com_path(os.getcwd(), reloutdir) + "/ "
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

def rewrite_version_header(where, major, minor, build):
	contents = "#define SIGNALIZER_MAJOR " + major + "\n#define SIGNALIZER_MINOR " + minor + "\n#define SIGNALIZER_BUILD " + build + "\n"
	with open(where, "w") as out:
		out.writelines(contents)

def create_build_file(where, vstring):
	# add latest git commit to build log
	git = subprocess.Popen("git --git-dir ../.git log -2", shell=True, stdout=subprocess.PIPE)
	git_log = git.stdout.read();

	build_info = strftime("%Y-%m-%d %H:%M:%S", gmtime()) + ": Signalizer " + vstring + " built on " + platform.system() + " " + platform.release() + " by " + getpass.getuser() + "\n\n"

	with open(where, "w") as out:
		out.writelines(build_info)
		out.writelines(git_log)

def set_plist_option(rel_plist_path, command):
	full_path = com_path(os.getcwd(), rel_plist_path)
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


com_path = os.path.join

# [0] = arg to clang, [1] = output folder
targets = [["i386", com_path(build_folder, "x32")], ["x86_64", com_path(build_folder, "x64")]]

# rewrite program internal version

rewrite_version_header("../Source/version.h", major, minor, build)

#run all targets
for option in targets:
	if compiler_invoke(option[0], version_string, option[1]) != 0:
		print("\n------> Error compiling for target " + option[0])
		sh.rmtree(build_folder)
		exit(1)
	else:
		create_build_file(com_path(com_path(com_path(com_path(option[1], "Signalizer.component"), "Contents"), "Resources"), "Build.log"), version_string)


print("\n------> All builds finished, generating plugin permutations ...")

# make build permutations and set up the plist file

for option in targets:
	original = com_path(option[1], "Signalizer.component")
	plist = com_path(com_path(original, "Contents"), "Info.plist")
	# set bundle description
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
	
	# permute
	sh.copytree(original, com_path(option[1], "Signalizer.vst"))
	sh.copytree(original, com_path(option[1], "Signalizer.vst3"))


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
