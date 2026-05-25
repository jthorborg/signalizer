import os
import argparse

parser = argparse.ArgumentParser(prog='build.py')
sub = parser.add_subparsers(dest='mode', required=True)

dev_p = sub.add_parser('dev', help='Debug Standalone build for development/testing')
dev_p.add_argument('--arch', default='x64', choices=['x64', 'arm64'])
dev_p.add_argument('--verbose', action='store_true', help='Stream full build output')
dev_p.add_argument('--errors', type=int, default=3, metavar='N',
                   help='Max errors to print on failure (default: 3)')

rel_p = sub.add_parser('release', help='Full release build with packaging')
rel_p.add_argument('-d', '--debug', action='store_true', help='Debug code generation')
rel_p.add_argument('-s', '--skipvst2', action='store_true',
                   help='Build without VST2 SDK dependency')
rel_p.add_argument('--verbose', action='store_true', help='Stream full build output')
rel_p.add_argument('-j', '--increase-major', action='store_true')
rel_p.add_argument('-n', '--increase-minor', action='store_true')
rel_p.add_argument('-b', '--increase-build', action='store_true')

args = parser.parse_args()

os.chdir('Make')

try:
    if args.mode == 'dev':
        import Make.common as cm
        # TODO: Implement development builds on other platforms.
        import Make.build_win as bw
        config = cm.DevConfig(args)
        print(f"------> Building Signalizer Debug Standalone {config.arch} targets")
        binary, log_path = bw.build_dev(config)
        print("------> Built Signalizer successfully into:")
        print(binary)
        if not config.verbose:
            print(f"------> Build log: {log_path}")

    elif args.mode == 'release':
        import Make.common as cm
        program = cm.ProgramConfig(args, 'config.ini')
        print("------> Building Signalizer v. " + program.version_string + " " + program.configString + " targets")
        program.rewrite_version_header()

        log_path = None
        if cm.is_linux:
            import Make.build_linux as bl
            zx = bl.build(program)
        elif cm.is_windows:
            import Make.build_win as bw
            zx, log_path = bw.build(program)
        elif cm.is_mac:
            import Make.build_osx as bo
            zx = bo.build(program)

        print("------> Built Signalizer successfully into:")
        print(zx)
        if log_path and not program.verbose:
            print(f"------> Build log: {log_path}")
        program.flush()

finally:
    os.chdir('..')
