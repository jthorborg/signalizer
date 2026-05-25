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

def build(program):

	vcxpath = "../Builds/VisualStudio2022"

	zipoutput = "../Releases/Signalizer_Windows_VST_" + program.version_string
	zippdboutput = "../Releases/Signalizer_Windows_Debug_PDBs_" + program.version_string

	#overwrite resource to embed version numbers
	setup_resource(cm.join(vcxpath, "resources.rc"), program.major, program.minor, program.build, program.name, program.desc, program.company)

	#run all archs
	archs = [["x64", f'"{program.configString}|x64"']]
	for option in archs:
		if compiler_invoke(option[0], option[1]) != 0:
			print("\n------> Error compiling for target " + option[0])
			exit(1)

	print("\n------> All builds finished, generating skeletons...")

	rootdir = "Signalizer Windows"
	program.make_release_folder_with_goodies(rootdir, "windows_installation_advice.txt")

	# VST2 section
	build_dir = cm.join(vcxpath, "x64", program.configString, "VST")
	output_dir = cm.join(rootdir, "Signalizer.vst")

	sh.copytree("Skeleton", output_dir)
	sh.copyfile(cm.join(build_dir, "Signalizer.dll"), cm.join(output_dir, "Signalizer.dll"))
	os.makedirs(cm.join("Symbols", "VST"))
	sh.copy(cm.join(build_dir, "Signalizer.pdb"), cm.join("Symbols", "VST", "Signalizer.pdb")) 

	# VST3 section
	build_dir = cm.join(vcxpath, "x64", program.configString, "VST3")
	output_dir = cm.join(rootdir, "Signalizer.vst3")

	sh.copytree(cm.join(build_dir, "Signalizer.vst3"), output_dir)
	sh.copytree("Skeleton", cm.join(output_dir, "Contents", "x86_64-win"), dirs_exist_ok=True)
	os.makedirs(cm.join("Symbols", "VST3"))
	sh.copy(cm.join(build_dir, "Signalizer.pdb"), cm.join("Symbols", "VST3", "Signalizer.pdb")) 

	print("------> Zipping output directories...")

	zx = sh.make_archive(zipoutput, "zip", rootdir)
	zxpdb = sh.make_archive(zippdboutput, "zip", "Symbols")

	# clean up dirs
	if os.path.exists(rootdir):
		sh.rmtree(rootdir)
	if os.path.exists("Symbols"):
		sh.rmtree("Symbols")

	return zx + "\n" + zxpdb
