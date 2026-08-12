#!/usr/bin/env python3
"""Build and launch the geomc frustum ray-intersection inspector."""

import argparse
import functools
import http.server
import pathlib
import shutil
import socketserver
import subprocess
import sys
import threading
import webbrowser


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--samples", type=int, default=64, help="samples to pre-generate (default: 64)")
    parser.add_argument("--rays", type=int, default=1000, help="lines per sample (default: 1000)")
    parser.add_argument("--seed", help="integer seed, decimal or 0x-prefixed (default: random)")
    parser.add_argument("--shapes", default="rectangle,circle,triangle",
                        help="comma-separated base shapes (default: rectangle,circle,triangle)")
    parser.add_argument("--port", type=int, default=8765, help="local HTTP port (default: 8765)")
    parser.add_argument("--no-browser", action="store_true", help="do not open a browser automatically")
    args = parser.parse_args()

    root = pathlib.Path(__file__).resolve().parent.parent
    build = root / "build" / "frustum_visualizer"
    build.mkdir(parents=True, exist_ok=True)
    executable = build / "frustum_visualizer_data"
    samples = build / "samples.json"

    compiler = shutil.which("clang++") or shutil.which("c++")
    if not compiler:
        parser.error("no C++ compiler found (expected clang++ or c++)")
    print("Building geomc…")
    subprocess.run(["scons", "lib"], cwd=root, check=True)
    compile_command = [compiler, "-std=c++23", "-O2", f"-I{root}",
                       str(root / "tools" / "frustum_visualizer.cpp"),
                       str(root / "build" / "native" / "geomc" / "libgeomc.a"),
                       "-o", str(executable)]
    print("Compiling sample generator…")
    subprocess.run(compile_command, check=True)

    generate_command = [str(executable), "--output", str(samples),
                        "--samples", str(args.samples), "--rays", str(args.rays),
                        "--shapes", args.shapes]
    if args.seed is not None:
        # Validate here while retaining the user's preferred spelling for the C++ tool.
        int(args.seed, 0)
        generate_command += ["--seed", args.seed]
    subprocess.run(generate_command, check=True)
    shutil.copy2(root / "tools" / "frustum_visualizer.html", build / "index.html")

    handler = functools.partial(http.server.SimpleHTTPRequestHandler, directory=build)
    class ReusableServer(socketserver.ThreadingTCPServer):
        allow_reuse_address = True
    try:
        with ReusableServer(("127.0.0.1", args.port), handler) as server:
            url = f"http://127.0.0.1:{args.port}/"
            print(f"Serving {url} — press Ctrl-C to stop")
            if not args.no_browser:
                threading.Timer(0.2, lambda: webbrowser.open(url)).start()
            server.serve_forever()
    except KeyboardInterrupt:
        print("\nStopped.")
    except OSError as error:
        parser.error(f"could not listen on port {args.port}: {error}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
