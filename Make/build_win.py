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

		
print("\n------> All builds finished, generating skeletons...")

cm.create_build_file("Build.log", version_string)
sh.copyfile("windows_installation_advice.txt", cm.join("Signalizer Windows", "HOW TO INSTALL.txt"))

option_to_build = { "x86": cm.join(cm.join(vcxpath, "Release")), "x64": cm.join(cm.join(cm.join(vcxpath, "x64")), "Release") }

for option in targets:
	release_dir = cm.join("Signalizer Windows", "Release " + version_string + " " + option[0], "Signalizer")
	debug_dir = cm.join("Signalizer Windows", "Debug " + version_string + " " + option[0], "Signalizer")
	build = option_to_build[option[0]]

	for p in [release_dir, debug_dir]:
		sh.copytree("Skeleton", p)
		sh.copyfile("Build.log", cm.join(p, "Build.log"))
		
	# copy in builds
	sh.copy(cm.join(build, "Signalizer.dll"), cm.join(release_dir, "Signalizer.dll"))
	sh.copy(cm.join(build, "Signalizer.dll"), cm.join(debug_dir, "Signalizer.dll"))
	# important that its name is signalizer.pdb
	sh.copy(cm.join(build, "Signalizer.pdb"), cm.join(debug_dir, "Signalizer.pdb")) 


print("------> Zipping output directories...")

zx = sh.make_archive(zipoutput, "zip", "Signalizer Windows")

print("------> Built Signalizer successfully into:")
print("------> " + zx)

# clean up dirs
sh.rmtree("Signalizer Windows")
os.remove("Build.log")
# done, if we made it here, increase the conf build

if flush_parameters:
	with open("config.ini", "w") as f:
		config.write(f, True)