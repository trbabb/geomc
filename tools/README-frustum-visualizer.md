# Frustum ray-intersection inspector

This opt-in browser tool generates random `Frustum<Base>` objects and calls the
actual C++ `Frustum::intersect()` implementation for a dense set of parameterized
lines. It displays the returned intersection segments in orange, the remainder
of each displayed line in gray, and the frustum wireframe in cyan. A malformed
(non-finite) interval is shown as a magenta diagnostic line and counted in the
HUD rather than being written as invalid JSON.

No graphics packages are required. The launcher compiles a small C++ data
generator and serves a dependency-free Canvas viewer locally.

## Run

From the repository root:

```sh
python3 tools/frustum_visualizer.py
```

Useful options:

```sh
python3 tools/frustum_visualizer.py \
  --rays 3000 \
  --samples 100 \
  --shapes rectangle,circle,triangle \
  --seed 0x41e754ec0ffa0b67
```

Run `python3 tools/frustum_visualizer.py --help` for all options. The seed shown
in the viewer can be passed back with `--seed` to reproduce the generated set.

## Controls

- **Left drag:** orbit
- **Mouse wheel / trackpad scroll:** zoom
- **Space** or **Right arrow:** next sample
- **Left arrow:** previous sample
- **R:** reset the camera

The library's `Ray` intersection methods return an interval over the complete
parameterized line `origin + s * direction`, including negative `s`; they are
not clipped to a forward-only `s >= 0` ray. The viewer intentionally preserves
that behavior.

The circle wireframe is tessellated for display only. Colored segments always
come from the concrete C++ base shape and `Frustum::intersect()`, never from the
wireframe mesh.
