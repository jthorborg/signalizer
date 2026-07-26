import io
import configparser as cp
import os
import sys
import shutil as sh
import zipfile as zip
import subprocess
import platform
import datetime
import re
import common as cm

from datetime import date

_ERROR_RE = re.compile(r': error:', re.IGNORECASE)

def build_dev(config):
	"""Standalone build for development/testing. Returns the path to the built executable."""
	xcode_arch = {"x64": "x86_64", "arm64": "arm64"}[config.arch]
	# ONLY_ACTIVE_ARCH restricts the build to the machine's own arch (fast, single-slice);
	# it's keyed off the actual host arch, not ARCHS, so requesting the *other* arch needs
	# ONLY_ACTIVE_ARCH=NO + an explicit ARCHS override instead.
	only_active = "YES" if xcode_arch == platform.machine() else "NO"

	command = [
		"xcodebuild",
		"-project", "../Builds/MacOSX/Signalizer.xcodeproj",
		"-scheme", "Signalizer - Standalone Plugin",
		"-configuration", config.configString,
		"ONLY_ACTIVE_ARCH=" + only_active,
		"ARCHS=" + xcode_arch,
	]

	if not config.verbose:
		command.append("-quiet")

	# Product > Archive (run manually in Xcode against this project) leaves a symlink at
	# build/<Config>/Signalizer.app pointing into DerivedData's ArchiveIntermediates, which
	# breaks mkdir -p for a normal build once DerivedData is cleaned. Clear it if present.
	app_bundle = cm.join("..", "Builds", "MacOSX", "build", config.configString, "Signalizer.app")
	if os.path.islink(app_bundle):
		os.unlink(app_bundle)

	os.makedirs(config.logs_dir, exist_ok=True)
	timestamp = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
	log_path = cm.join(config.logs_dir, f'build_{timestamp}.log')

	if config.verbose:
		print("---------> Compiler invocation: \n" + " ".join(command))
		proc = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
		lines = []
		for line in proc.stdout:
			print(line, end='')
			lines.append(line)
		proc.wait()
		returncode = proc.returncode
	else:
		proc = subprocess.run(command, capture_output=True, text=True)
		lines = proc.stdout.splitlines(keepends=True) + proc.stderr.splitlines(keepends=True)
		returncode = proc.returncode

	with open(log_path, 'w') as f:
		f.writelines(lines)

	if returncode != 0:
		print(f"------> Build failed (full log: {log_path})")
		if not config.verbose:
			errors = [l for l in lines if _ERROR_RE.search(l)]
			for line in errors[:config.max_errors]:
				print(line, end='')
		exit(1)

	build_dir = cm.join("..", "Builds", "MacOSX", "build", config.configString)
	binary = os.path.abspath(cm.join(build_dir, "Signalizer.app", "Contents", "MacOS", "Signalizer"))
	return binary, (None if config.verbose else log_path)

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
	zipoutput = "../Releases/Signalizer_MacOS_" + program.version_string

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
	plist_variants = [root_plist + list + ".plist" for list in ["-AU", "-VST", "-VST3", "-VST3_Manifest_Helper", "-Standalone_Plugin", ""]]

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

	build_variants = ["AU", "VST3", "Standalone Plugin"]

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
