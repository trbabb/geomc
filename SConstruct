import itertools
import platform

import os

AddOption(
    '--wasm',
    action='store_true',
    help='Compile to WASM',
    default=False)
AddOption(
    '--sanitize',
    action='store_true',
    help='Enable address sanitization',
    default=False)

debug  = ARGUMENTS.get('debug',  False)
noopt  = ARGUMENTS.get('noopt',  False)

DOC_DIR = 'doc/gen'

def subst_wasm_flags(env, wasm_flags):
    wasm_opts = (('-s', f'{k}={v}') for k,v in wasm_flags.items())
    return list(itertools.chain(*wasm_opts))

wasm_cflags = {
    'NO_DISABLE_EXCEPTION_CATCHING' : 1,
}

wasm_linkflags = {
    'ALLOW_MEMORY_GROWTH' : 1,
}

pkg_config_path = os.environ.get('PKG_CONFIG_PATH') if 'PKG_CONFIG_PATH' in os.environ else []

if debug:
    wasm_cflags['ASSERTIONS'] = 1

env = Environment(
    CXX='clang++',
    CXXFLAGS=[
        '-O3' if not noopt else '-O0',
        '-std=c++20',
        '-fcolor-diagnostics',
        '-Wall',
        '-Werror',
        '-Wno-return-type-c-linkage',
        '-Wno-limited-postlink-optimizations',
        '-Wno-unknown-warning-option',
        # '-v'
    ],
    LIBS=[],
    CPPPATH=['#'],
    ENV={'PKG_CONFIG_PATH': pkg_config_path},
    COMPILATIONDB_USE_ABSPATH=True,
)
env.AddMethod(subst_wasm_flags)
env.Tool('compilation_db')

if env['PLATFORM'] == 'darwin':
    # homebrew / pkgconfig not in the path by default
    env.AppendENVPath('PATH', ['/opt/homebrew/bin','/usr/local/bin'])

if debug:
    env.Append(CXXFLAGS='-g')

if GetOption('sanitize'):
    env.Append(CXXFLAGS =['-fsanitize=address', '-fno-omit-frame-pointer'])
    env.Append(LINKFLAGS=['-fsanitize=address', '-fno-omit-frame-pointer'])

if GetOption('wasm'):
    arch           = 'wasm'
    env['CXX']     = 'em++'
    env['AR']      = 'emar'
    env['RANLIB']  = 'emranlib'
    # todo: make this more portable.
    env.Append(LIBPATH=['/usr/local/wasm/lib'])
    env.Append(CPPPATH=['/usr/local/wasm/include'])
    env.Append(CXXFLAGS=env.subst_wasm_flags(wasm_cflags))
    env.Append(LINKFLAGS=env.subst_wasm_flags(wasm_linkflags))
    env.Append(LIBS=[])
    env['LIBPATH'] = [f'#build/{arch}/geomc', '/usr/local/wasm/lib']
    if debug:
        env.Append(LINKFLAGS=['-gsource-map'])
else:
    arch = 'native'
    if env['PLATFORM'] == 'darwin':
        env.Append(LINKFLAGS=[
            '-Wl,-rpath,@loader_path',
            '-Wl,-rpath,/usr/local/lib',
        ])
    elif env['PLATFORM'] == 'linux':
        env.Append(LINKFLAGS=[
            '-Wl,-rpath,$ORIGIN',
            '-Wl,-rpath,/usr/local/lib'
        ])
    env.Append(CXXFLAGS='-march=native')
    env.Append(LIBPATH=['/usr/local/lib', f'#build/{arch}/geomc'])

env['ARCH'] = arch

docs = env.Command('docs', None, [Mkdir(DOC_DIR), 'doxygen'])
env.Alias('docs', docs)

Export("env")

# data = SConscript('objects/SConscript', variant_dir='build/objects')
lib  = SConscript('geomc/SConscript',      variant_dir=f'build/{arch}/geomc')
test = SConscript('regression/SConscript', variant_dir=f'build/{arch}/regression',
    exports={'lib': lib}
)
comp_db = env.CompilationDatabase(target='compile_commands.json')

env.Alias('lib', lib)
env.Alias('compile_commands', comp_db)

Default([lib, test, comp_db])
