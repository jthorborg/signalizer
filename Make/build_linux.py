import io
import ConfigParser
import os
import sys
import shutil as sh
import zipfile as zip
import common as cm
import subprocess

def compiler_invoke(args):
	return os.system("codeblocks ../Builds/CodeBlocks/Signalizer.cbp " + args + " --rebuild")

# parse config
config = ConfigParser.ConfigParser()
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


version_string = major + "." + minor + "." + build
zipoutput = "../Releases/Signalizer Linux VST " + version_string
#diagnostic
print("------> Building Signalizer v. " + version_string + " release targets")

cm.rewrite_version_header("../Source/version.h", major, minor, build)

#run targets

compiler_invoke("--target=Release")

# output dirs
release_dir = "Signalizer Linux Release " + version_string

cm.create_build_file("Build.log", version_string)

# build skeleton
sh.copytree("Skeleton", release_dir)
sh.copyfile("Build.log", cm.join(release_dir, "Build.log"))


print("\n------> All builds finished, generating skeletons...")

# copy in builds
sh.copy("../Builds/CodeBlocks/bin/Release/libSignalizer.so", cm.join(release_dir, "Signalizer.so"))

print("------> Zipping output directories...")

zx = sh.make_archive(zipoutput, "zip", release_dir)

print("------> Builded Signalizer successfully into:")
print("------> " + zx)

# clean up dirs
sh.rmtree(release_dir)
os.remove("Build.log")
# done, if we made it here, increase the conf build

if flush_parameters:
	with open("config.ini", "w") as f:
		config.write(f, True)
