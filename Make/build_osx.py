import io
import configparser as cp
import os
import sys
import shutil as sh
import zipfile as zip
import subprocess
import common as cm

from datetime import date

def compiler_invoke(scheme, vstring, reloutdir, config):
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
			   "-configuration " + config + " "
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

def build(program):

	version_int = (int(program.major) << 16) | (int(program.minor) << 8) | int(program.build)
	zipoutput = "../Releases/Signalizer_macOS_" + program.version_string

	#diagnostic
	print("------> Cleaning prior builds... ")

	if os.system("xattr -w com.apple.xcode.CreatedByBuildSystem true ../Builds/MacOSX/build") != 0:
		print("------> Failed changing xattr for cleaning ...")
		exit(-2)
		
	if os.system("xcodebuild -project ../builds/macosx/signalizer.xcodeproj -scheme \"Signalizer - All\" clean ONLY_ACTIVE_ARCH=NO") != 0:
		print("------> Failed cleaning...")
		exit(-3)

	# rewrite build plist
	root_plist = cm.join("../Builds/MacOSX/Info")
	plist_variants = [root_plist + list + ".plist" for list in ["-AU", "-VST", "-VST3", "-VST3_Manifest_Helper", ""]]

	for plist in plist_variants:
		print("------> rewriting plist " + plist)
		set_plist_option(plist, "Set :CFBundleIdentifier com." + program.company + "." + program.name)
		set_plist_option(plist, "Set :CFBundleShortVersionString " + program.version_string)
		set_plist_option(plist, "Set :CFBundleVersion " + program.version_string)
		set_plist_option(plist, "Set :NSHumanReadableCopyright Copyright (c) " + str(date.today().year) + " " + program.author)

	# set the audio unit plugin description
	aulist = plist_variants[0]
	set_plist_option(aulist, "Set :AudioComponents:0:description " + program.desc)
	set_plist_option(aulist, "Set :AudioComponents:0:manufacturer " + program.manu4)
	set_plist_option(aulist, "Set :AudioComponents:0:name " + program.company + ": " + program.name)
	set_plist_option(aulist, "Set :AudioComponents:0:subtype " + program.sub4)
	set_plist_option(aulist, "Set :AudioComponents:0:type aufx")
	set_plist_option(aulist, "Set :AudioComponents:0:version " + str(version_int))

	build_variants = ["AU", "VST3"]

	if not program.skipvst2:
		build_variants = build_variants + ["VST"]
		
	print(build_variants)
	
	build_folder = "Signalizer macOS Universal"
	program.make_release_folder_with_goodies(build_folder, "macos_installation_advice.txt")

	#run all targets
	for build in build_variants:
		if compiler_invoke(build, program.version_string, build_folder, program.configString) != 0:
			print("\n------> Error compiling for target " + build)
			sh.rmtree(build_folder)
			exit(-4)

	#remove temporary junk
	if os.path.exists(cm.join(build_folder, "libSignalizer.a")) != 0:
		os.remove(cm.join(build_folder, "libSignalizer.a"))

	if os.path.exists(cm.join(build_folder, "juce_vst3_helper")) != 0:
		os.remove(cm.join(build_folder, "juce_vst3_helper"))

	#append extra goodies
	print("------> Zipping output directories...")

	zx = sh.make_archive(zipoutput, "zip", build_folder)

	# clean up dirs
	sh.rmtree(build_folder)

	return zx
