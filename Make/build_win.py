import io
import configparser
import os
import sys
import shutil as sh
import zipfile as zip
import common as cm
import subprocess
import re
import datetime

# Matches MSBuild/MSVC error lines: "path(line): error CODE: message"
_ERROR_RE = re.compile(r': error \w', re.IGNORECASE)


_VSWHERE = r'C:\Program Files (x86)\Microsoft Visual Studio\Installer\vswhere.exe'
_VS_EDITIONS = ('Community', 'Professional', 'Enterprise', 'BuildTools')


def find_vcvarsall():
	"""Locate vcvarsall.bat for VS2022, regardless of edition."""
	if os.path.exists(_VSWHERE):
		result = subprocess.run(
			[_VSWHERE, '-latest', '-products', '*',
			 '-requires', 'Microsoft.VisualStudio.Component.VC.Tools.x86.x64',
			 '-property', 'installationPath'],
			capture_output=True, encoding='oem'
		)
		install_path = result.stdout.strip()
		if install_path:
			candidate = os.path.join(install_path, 'VC', 'Auxiliary', 'Build', 'vcvarsall.bat')
			if os.path.exists(candidate):
				return candidate

	for edition in _VS_EDITIONS:
		candidate = rf'C:\Program Files\Microsoft Visual Studio\2022\{edition}\VC\Auxiliary\Build\vcvarsall.bat'
		if os.path.exists(candidate):
			return candidate

	raise RuntimeError(
		"Could not locate vcvarsall.bat for VS2022 (checked vswhere and common install paths)"
	)


def get_msvc_env(arch='x64'):
	"""Return os.environ copy with MSVC toolchain loaded for the given arch."""
	vcvarsall = find_vcvarsall()
	cmd = f'"{vcvarsall}" {arch} && set'
	result = subprocess.run(cmd, capture_output=True, encoding='oem', shell=True)
	env = {}
	for line in result.stdout.splitlines():
		if '=' in line:
			k, v = line.split('=', 1)
			env[k.upper()] = v
	return env


def run_msbuild(projects, config_str, arch, env, verbose, max_errors, logs_dir, solution_dir=None):
	"""
	Build a list of project/solution files in order using MSBuild.
	Returns (success: bool, error_lines: list[str], log_path: str).
	solution_dir is only needed when building .vcxproj files directly.
	In verbose mode, output is streamed and error_lines is always empty.
	"""
	os.makedirs(logs_dir, exist_ok=True)
	timestamp = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
	log_path = os.path.join(logs_dir, f'build_{timestamp}.log')

	msbuild = sh.which('msbuild', path=env.get('PATH', ''))
	if not msbuild:
		raise RuntimeError("msbuild not found in MSVC environment PATH")

	all_lines = []

	for project in projects:
		cmd = [
			msbuild, project,
			f'/p:Configuration={config_str}',
			f'/p:Platform={arch}',
			'/nologo',
			'/verbosity:quiet' if not verbose else '/verbosity:normal',
		]
		if solution_dir is not None:
			cmd.append(f'/p:SolutionDir={solution_dir}')

		if verbose:
			proc = subprocess.Popen(cmd, env=env, stdout=subprocess.PIPE,
			                        stderr=subprocess.STDOUT, text=True)
			log_lines = []
			for line in proc.stdout:
				print(line, end='')
				log_lines.append(line)
			proc.wait()
			all_lines.extend(log_lines)

			if proc.returncode != 0:
				with open(log_path, 'w') as f:
					f.writelines(all_lines)
				return False, [], log_path
		else:
			proc = subprocess.run(cmd, capture_output=True, text=True, env=env)
			lines = proc.stdout.splitlines() + proc.stderr.splitlines()
			all_lines.extend(lines)

			if proc.returncode != 0:
				with open(log_path, 'w') as f:
					f.write('\n'.join(all_lines))
				errors = [l for l in lines if _ERROR_RE.search(l)]
				return False, errors[:max_errors], log_path

	with open(log_path, 'w') as f:
		f.write('\n'.join(all_lines))

	return True, [], log_path


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

def build_dev(config):
	"""Standalone build for development. Returns the path to the built executable."""
	vcxpath = os.path.abspath(cm.join("..", "Builds", "VisualStudio2022"))
	solution_dir = vcxpath + "\\"

	projects = [
		cm.join(vcxpath, "Signalizer_SharedCode.vcxproj"),
		cm.join(vcxpath, "Signalizer_StandalonePlugin.vcxproj"),
	]

	env = get_msvc_env(config.arch)
	success, errors, log_path = run_msbuild(
		projects, config.configString, config.arch,
		env, config.verbose, config.max_errors, config.logs_dir,
		solution_dir=solution_dir
	)

	if not success:
		print(f"------> Build failed (full log: {log_path})")
		for line in errors:
			print(line)
		exit(1)

	return cm.join(vcxpath, config.arch, config.configString, "Standalone Plugin", "Signalizer.exe"), log_path


def build(program):

	vcxpath = "../Builds/VisualStudio2022"

	zipoutput = "../Releases/Signalizer_Windows_VST_" + program.version_string
	zippdboutput = "../Releases/Signalizer_Windows_Debug_PDBs_" + program.version_string

	#overwrite resource to embed version numbers
	setup_resource(cm.join(vcxpath, "resources.rc"), program.major, program.minor, program.build, program.name, program.desc, program.company)

	env = get_msvc_env('x64')
	success, errors, log_path = run_msbuild(
		[cm.join(vcxpath, "Signalizer.sln")], program.configString, 'x64',
		env, program.verbose, 5, 'Logs'
	)
	if not success:
		print(f"\n------> Build failed (full log: {log_path})")
		for line in errors:
			print(line)
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

	# Standalone section
	build_dir = cm.join(vcxpath, "x64", program.configString, "Standalone Plugin")
	output_dir = cm.join(rootdir, "Signalizer")

	sh.copytree("Skeleton", output_dir)
	sh.copyfile(cm.join(build_dir, "Signalizer.exe"), cm.join(output_dir, "Signalizer.exe"))
	os.makedirs(cm.join("Symbols", "Standalone"))
	sh.copy(cm.join(build_dir, "Signalizer.pdb"), cm.join("Symbols", "Standalone", "Signalizer.pdb"))

	print("------> Zipping output directories...")

	zx = sh.make_archive(zipoutput, "zip", rootdir)
	zxpdb = sh.make_archive(zippdboutput, "zip", "Symbols")

	# clean up dirs
	if os.path.exists(rootdir):
		sh.rmtree(rootdir)
	if os.path.exists("Symbols"):
		sh.rmtree("Symbols")

	return zx + "\n" + zxpdb, log_path
