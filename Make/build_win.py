import io
import configparser
import os
import sys
import shutil as sh
import zipfile as zip
import common as cm
import subprocess

def kvsz(key, value):
	return "      VALUE " + "\"" + key + "\", \"" + value + "\\0\"\n"

def setup_resource(outputfile, major, minor, build, name, description, company):
	version_comma = str(major) + "," + str(minor) + "," + str(build) + ",0"
	version_dot = "\"" + str(major) + "." + str(minor) + "." + str(build) + "\""

	contents = ("#pragma code_page(65001)\n\n"
				"#ifdef JUCE_USER_DEFINED_RC_FILE\n"
				" #include JUCE_USER_DEFINED_RC_FILE\n"
				"#else\n"
				"#undef  WIN32_LEAN_AND_MEAN\n"
				"#define WIN32_LEAN_AND_MEAN\n"
				"#include <windows.h>\n"
				"VS_VERSION_INFO VERSIONINFO\n"
				"FILEVERSION " + version_comma + "\n"
				"BEGIN\n"
				"  BLOCK \"StringFileInfo\"\n" 
				"  BEGIN\n"
				"    BLOCK \"040904E4\"\n" 
				"    BEGIN\n" +
				kvsz("CompanyName", company) +
				kvsz("FileDescription", description) +
				kvsz("FileVersion", version_dot) +
				kvsz("ProductName", name) +
				kvsz("ProductVersion", version_dot) +
				"    END\n" 
				"  END\n"
				"  BLOCK \"VarFileInfo\"\n"
				"  BEGIN\n"
				"    VALUE \"Translation\", 0x409, 1252\n"
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
company = config.get("info", "company")
desc = config.get("info", "description")
name = config.get("info", "productname")


version_string = major + "." + minor + "." + build
vcxpath = "../Builds/VisualStudio2022"

zipoutput = "../Releases/Signalizer_Windows_VST_" + version_string
zippdboutput = "../Releases/Signalizer_Windows_Debug_PDBs_" + version_string

#diagnostic
print("------> Building Signalizer v. " + version_string + " release targets")

#overwrite resource to embed version numbers
setup_resource(cm.join(vcxpath, "resources.rc"), major, minor, build, name, desc, company)
cm.rewrite_version_header("../Source/version.h", major, minor, build)

#run all archs
archs = [["x64", '"Release|x64"']]
for option in archs:
	if compiler_invoke(option[0], option[1]) != 0:
		print("\n------> Error compiling for target " + option[0])
		exit(1)

print("\n------> All builds finished, generating skeletons...")

rootdir = "Signalizer Windows"

# VST2 section
build_dir = cm.join(vcxpath, "x64", "Release", "VST")
output_dir = cm.join(rootdir, "Signalizer.vst")

sh.copytree("Skeleton", output_dir)
sh.copyfile(cm.join(build_dir, "Signalizer.dll"), cm.join(output_dir, "Signalizer.dll"))
os.makedirs(cm.join("Symbols", "VST"))
sh.copy(cm.join(build_dir, "Signalizer.pdb"), cm.join("Symbols", "VST", "Signalizer.pdb")) 

# VST3 section
build_dir = cm.join(vcxpath, "x64", "Release", "VST3")
output_dir = cm.join(rootdir, "Signalizer.vst3")

sh.copytree(cm.join(build_dir, "Signalizer.vst3"), output_dir)
sh.copytree("Skeleton", cm.join(output_dir, "Contents", "x86_64-win"), dirs_exist_ok=True)
os.makedirs(cm.join("Symbols", "VST3"))
sh.copy(cm.join(build_dir, "Signalizer.pdb"), cm.join("Symbols", "VST3", "Signalizer.pdb")) 

# Shared
cm.create_build_file("Build.log", version_string)
sh.copyfile(cm.join("Skeleton", "READ ME.txt"), cm.join(rootdir, "READ ME.txt"))
sh.copyfile("Build.log", cm.join(rootdir, "Build.log"))
sh.copyfile("../CHANGELOG.md", cm.join(rootdir, "CHANGELOG.md"))

sh.copyfile("windows_installation_advice.txt", cm.join(rootdir, "HOW TO INSTALL.txt"))

print("------> Zipping output directories...")

zx = sh.make_archive(zipoutput, "zip", rootdir)
zxpdb = sh.make_archive(zippdboutput, "zip", "Symbols")

print("------> Built Signalizer successfully into:")
print("------> " + zx)
print("------> " + zxpdb)

# clean up dirs
if os.path.exists(rootdir):
	sh.rmtree(rootdir)
if os.path.exists("Symbols"):
	sh.rmtree("Symbols")

os.remove("Build.log")
# done, if we made it here, increase the conf build

if flush_parameters:
	with open("config.ini", "w") as f:
		config.write(f, True)