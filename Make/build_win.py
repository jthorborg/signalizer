import io
import configparser
import os
import sys
import shutil as sh
import zipfile as zip

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

# handle operations
for param in parameters:
	config.set("version", param, str(int(config.get("version", param)) + 1))
	print("------> Increasing " + param + " to " + config.get("version", param))

# write new configuration?
if len(parameters) > 0:
	with open("config.ini", "w") as f:
		config.write(f, True)

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
com_path = os.path.join

#overwrite resource to embed version numbers
setup_resource(com_path(vcxpath, "resource.rc"), major, minor, build, name, desc)

targets = [["x86", '"Release|win32"'], ["x64", '"Release|x64"']]

#run all targets

for option in targets:
	if compiler_invoke(option[0], option[1]) != 0:
		print("\n------> Error compiling for target " + option[0])
		exit(1)


# output dirs
release_dir = com_path("Signalizer Windows", "Release " + version_string)
release_debug_dir = com_path("Signalizer Windows", "Debug " + version_string)

# build skeleton
for p in [release_debug_dir, release_dir]:
	sh.copytree("Skeleton", p, option[0])

print("\n------> All builds finished, generating skeletons...")

# copy in builds
x32release = com_path(com_path(vcxpath, "Release"))
x64release = com_path(com_path(com_path(vcxpath, "x64")), "Release")

sh.copy(com_path(x32release, "Signalizer.dll"), com_path(release_dir, "Signalizer " + targets[0][0] + ".dll"))
sh.copy(com_path(x64release, "Signalizer.dll"), com_path(release_dir, "Signalizer " + targets[1][0] + ".dll"))

sh.copy(com_path(x32release, "Signalizer.dll"), com_path(release_debug_dir, "Signalizer " + targets[0][0] + ".dll"))
sh.copy(com_path(x64release, "Signalizer.dll"), com_path(release_debug_dir, "Signalizer " + targets[1][0] + ".dll"))
sh.copy(com_path(x64release, "Signalizer.pdb"), com_path(release_debug_dir, "Signalizer " + targets[0][0] + ".pdb"))
sh.copy(com_path(x32release, "Signalizer.pdb"), com_path(release_debug_dir, "Signalizer " + targets[1][0] + ".pdb"))

print("------> Zipping output directories...")

zx = sh.make_archive(zipoutput, "zip", "Signalizer Windows")

print("------> Builded Signalizer successfully into:")
print("------> " + zx)

# clean up dirs
sh.rmtree("Signalizer Windows")