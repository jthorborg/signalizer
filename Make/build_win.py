import io
import configparser
import os
import sys
import shutil as sh
import zipfile as zip
import common as cm
import subprocess

def setup_resource(outputfile, major, minor, build, name, description):
	version_comma = str(major) + "," + str(minor) + "," + str(build) + ",0"
	version_dot = "\"" + str(major) + "." + str(minor) + "." + str(build) + "\""

	contents = ("#ifdef JUCE_USER_DEFINED_RC_FILE\n"
				" #include JUCE_USER_DEFINED_RC_FILE\n"
				"#else\n"
				"#undef  WIN32_LEAN_AND_MEAN\n"
				"#define WIN32_LEAN_AND_MEAN\n"
				"#include <windows.h>\n"
				"VS_VERSION_INFO VERSIONINFO\n"
				"FILEVERSION " + version_comma + "\n"
				"PRODUCTVERSION " + version_comma + "\n"
				"BEGIN\n"
				"  BLOCK \"StringFileInfo\"\n" 
				"  BEGIN\n"
				"    BLOCK \"040904E4\"\n" 
				"    BEGIN\n"
				"      VALUE \"FileDescription\", \"" + description + "\"\n" 
				"      VALUE \"FileVersion\", " + version_dot + "\n"
				"      VALUE \"ProductName\", \"" + name + "\"\n"
				"      VALUE \"ProductVersion\", " + version_dot + "\n"
				"    END\n" 
				"  END\n"
				"  BLOCK \"VarFileInfo\"\n"
				"  BEGIN\n"
				"    VALUE \"Translation\", 0x409, 65001\n"
				"  END\n"
				"END\n"
				"#endif\n")

	with open(outputfile, "w") as out:
		out.writelines(contents)

def compiler_invoke(compiler_arch, target):
	return os.system("vscompile.bat " + compiler_arch + " " + target)

# parse config
config = configparser.ConfigParser()
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
vcxpath = "../builds/visualstudio2010"
zipoutput = "../Releases/Signalizer Windows VST " + version_string
#diagnostic
print("------> Building Signalizer v. " + version_string + " release targets")

#overwrite resource to embed version numbers
setup_resource(cm.join(vcxpath, "resources.rc"), major, minor, build, name, desc)

targets = [["x86", '"Release|win32"'], ["x64", '"Release|x64"']]



cm.rewrite_version_header("../Source/version.h", major, minor, build)

#run all targets

for option in targets:
	if compiler_invoke(option[0], option[1]) != 0:
		print("\n------> Error compiling for target " + option[0])
		exit(1)


# output dirs
release_dir = cm.join("Signalizer Windows", "Release " + version_string)
release_debug_dirs = [cm.join("Signalizer Windows", "Debug " + version_string + " " + targets[0][0]), cm.join("Signalizer Windows", "Debug " + version_string + " " + targets[1][0])]

cm.create_build_file("Build.log", version_string)

# build skeleton
for p in [release_debug_dirs[0], release_debug_dirs[1], release_dir]:
    sh.copytree("Skeleton", p)
    sh.copyfile("Build.log", cm.join(p, "Build.log"))


print("\n------> All builds finished, generating skeletons...")

# copy in builds
x32release = cm.join(cm.join(vcxpath, "Release"))
x64release = cm.join(cm.join(cm.join(vcxpath, "x64")), "Release")

# who needs for loops anyway
sh.copy(cm.join(x32release, "Signalizer.dll"), cm.join(release_dir, "Signalizer " + targets[0][0] + ".dll"))
sh.copy(cm.join(x64release, "Signalizer.dll"), cm.join(release_dir, "Signalizer " + targets[1][0] + ".dll"))

sh.copy(cm.join(x32release, "Signalizer.dll"), cm.join(release_debug_dirs[0], "Signalizer " + targets[0][0] + ".dll"))
sh.copy(cm.join(x32release, "Signalizer.pdb"), cm.join(release_debug_dirs[0], "Signalizer.pdb")) # important that its name is signalizer.pdb
sh.copy(cm.join(x64release, "Signalizer.dll"), cm.join(release_debug_dirs[1], "Signalizer " + targets[1][0] + ".dll"))
sh.copy(cm.join(x64release, "Signalizer.pdb"), cm.join(release_debug_dirs[1], "Signalizer.pdb")) # see above

print("------> Zipping output directories...")

zx = sh.make_archive(zipoutput, "zip", "Signalizer Windows")

print("------> Builded Signalizer successfully into:")
print("------> " + zx)

# clean up dirs
sh.rmtree("Signalizer Windows")
os.remove("Build.log")
# done, if we made it here, increase the conf build

if flush_parameters:
	with open("config.ini", "w") as f:
		config.write(f, True)