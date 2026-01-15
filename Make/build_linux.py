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

program = cm.ProgramConfig("config.ini")
zipoutput = "../Releases/Signalizer Linux VST " + program.version_string

#diagnostic
config = "Release" if program.release else "Debug"

print("------> Building Signalizer v. " + program.version_string + " " + config + " targets")

program.rewrite_version_header()

#run targets
if program.release:
	if compiler_invoke("clean") != 0:
		print("------> Error cleaning...")
		exit(-1)

if compiler_invoke("CONFIG=" + config) != 0:
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

print("------> Zipping output directories...")

zx = sh.make_archive(zipoutput, "zip", rootdir)

print("------> Built Signalizer successfully into:")
print("------> " + zx)

# clean up dirs
sh.rmtree(rootdir)

# done, if we made it here, increase the conf build
program.flush()
