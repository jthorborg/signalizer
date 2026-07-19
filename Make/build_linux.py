import io
import os
import sys
import shutil as sh
import zipfile as zip
import common as cm
import subprocess

make_dir = cm.join("..", "Builds", "LinuxMakefile")

def compiler_invoke(args):
	return os.system("make --directory=" + make_dir + " " + args)

def build_dev(config):
	"""Standalone build for development/testing. Returns the path to the built executable."""
	if compiler_invoke("CONFIG=" + config.configString + " Standalone") != 0:
		print("------> Error building...")
		exit(1)

	build_dir = os.path.abspath(cm.join(make_dir, "build"))

	# merge the resource skeleton (presets/resources/licenses) next to the binary
	# so the standalone is runnable straight out of its build location, mirroring
	# what the VS2022 exporter's postbuildCommand does on Windows.
	sh.copytree("Skeleton", build_dir, dirs_exist_ok=True)

	return cm.join(build_dir, "Signalizer"), None

def build(program):
	zipoutput = "../Releases/Signalizer Linux VST " + program.version_string

	#run targets
	if program.release:
		if compiler_invoke("clean") != 0:
			print("------> Error cleaning...")
			exit(-1)

	if compiler_invoke("CONFIG=" + program.configString) != 0:
		print("------> Error building...")

	print("\n------> All builds finished, generating skeletons...")

	# output dirs
	rootdir = "Signalizer Linux"

	# build skeleton
	program.make_release_folder_with_goodies(rootdir, "linux_installation_advice.txt")

	build_dir = cm.join(make_dir, "build")

	# VST2
	output_dir = cm.join(rootdir, "Signalizer.vst")

	sh.copytree("Skeleton", output_dir)
	sh.copyfile(cm.join(build_dir, "Signalizer.so"), cm.join(output_dir, "Signalizer.so"))

	# VST3 section
	output_dir = cm.join(rootdir, "Signalizer.vst3")

	sh.copytree(cm.join(build_dir, "Signalizer.vst3"), output_dir)
	sh.copytree("Skeleton", cm.join(output_dir, "Contents", "x86_64-linux"), dirs_exist_ok=True)

	# Standalone section
	output_dir = cm.join(rootdir, "Signalizer")

	sh.copytree("Skeleton", output_dir)
	sh.copy(cm.join(build_dir, "Signalizer"), cm.join(output_dir, "Signalizer"))

	print("------> Zipping output directories...")

	zx = sh.make_archive(zipoutput, "zip", rootdir)

	# clean up dirs
	sh.rmtree(rootdir)

	return zx