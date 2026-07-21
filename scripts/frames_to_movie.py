#!/usr/bin/env python3
"""
Assemble already-rendered PNG frames into a looping GIF and/or a truecolor MP4.

Pure Python, no PyMOL dependency -- this is a post-processing step you run
AFTER rendering frames (e.g. with rfd3_movie / style_fixed / apply_camera in
.pymolrc), not something PyMOL itself needs to know about.

Two gotchas this avoids (see docs/RFD3_FIGURES.md for the full writeups):

- A single global GIF palette sampled from only a handful of frames tends to
  wash out whichever frames it wasn't sampled from -- if your longest-held or
  most-important frame is near the end (e.g. a final "clean" panel), a sparse
  sample can miss it entirely. `frames_to_gif` quantizes each frame to its own
  local palette instead: slightly bigger files, much more faithful colors.
- ffmpeg's `-f concat` demuxer, when durations are set per source file, can
  write an empty edit-list entry for the first clip -- players honor that as
  a gap and the video opens on black, silently dropping the intended hold on
  frame 1. `frames_to_mp4` sidesteps this by writing out a plain numbered
  image sequence (held frames duplicated as symlinks -- instant, no pixel
  copy) and encoding that with the `image2` demuxer instead, which starts
  cleanly at t=0. Verify with `ffprobe -show_entries stream=start_time
  -show_entries format=duration out.mp4` -- start_time should read
  0.000000.

Usage as a library:

    from frames_to_movie import frames_to_gif, frames_to_mp4
    frames_to_gif(["a.png", "b.png", "c.png"], durations_ms=[100]*3, out_gif="out.gif")
    frames_to_mp4([("a.png", 100), ("b.png", 100), ("c.png", 400)], out_mp4="out.mp4")

Or from the command line, one frame per --frame flag (duration in ms):

    python3 frames_to_movie.py --gif out.gif --frame a.png 100 --frame b.png 100
    python3 frames_to_movie.py --mp4 out.mp4 --fps 25 --frame a.png 100 --frame b.png 400
"""
import argparse
import os
import shutil
import subprocess
import tempfile


def frames_to_gif(image_paths, durations_ms, out_gif, loop=0):
    """Assemble PNGs into a looping GIF. `durations_ms` is either one int
    (applied to every frame) or a list matching `image_paths`. Each frame gets
    its OWN locally-quantized 256-color palette (see module docstring) --
    slightly larger output than a single shared palette, but frames near the
    end don't wash out."""
    from PIL import Image

    if not image_paths:
        raise ValueError("frames_to_gif: no image_paths given")
    if isinstance(durations_ms, (int, float)):
        durations_ms = [int(durations_ms)] * len(image_paths)
    if len(durations_ms) != len(image_paths):
        raise ValueError("frames_to_gif: durations_ms length must match image_paths")

    quantized = [
        Image.open(p).convert("RGB").quantize(
            colors=256, method=Image.MEDIANCUT, dither=Image.Dither.FLOYDSTEINBERG
        )
        for p in image_paths
    ]
    quantized[0].save(
        out_gif, save_all=True, append_images=quantized[1:],
        duration=durations_ms, loop=loop, disposal=2,
    )
    return out_gif


def frames_to_mp4(specs, out_mp4, fps=25, crf=16, preset="medium"):
    """specs = [(image_path, duration_ms), ...]. Encodes a truecolor MP4 with
    correct hold timing by duplicating held frames as symlinks into a temp
    directory (instant, no pixel copy) and encoding that as a plain numbered
    image sequence -- avoids the concat-demuxer black-intro bug (see module
    docstring). Requires ffmpeg on PATH."""
    if shutil.which("ffmpeg") is None:
        raise RuntimeError("frames_to_mp4: ffmpeg not found on PATH")
    if not specs:
        raise ValueError("frames_to_mp4: no specs given")

    tmp = tempfile.mkdtemp(prefix="frames_to_mp4_")
    n = 0
    try:
        for path, duration_ms in specs:
            repeats = max(1, round(duration_ms / (1000.0 / fps)))
            for _ in range(repeats):
                os.symlink(os.path.abspath(path), os.path.join(tmp, "f_%05d.png" % n))
                n += 1
        subprocess.run(
            ["ffmpeg", "-y", "-v", "error", "-framerate", str(fps),
             "-i", os.path.join(tmp, "f_%05d.png"),
             "-c:v", "libx264", "-crf", str(crf), "-preset", preset,
             "-pix_fmt", "yuv420p", "-movflags", "+faststart", out_mp4],
            check=True,
        )
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    return out_mp4


def _main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--frame", nargs=2, action="append", metavar=("PATH", "DURATION_MS"),
                     required=True, help="one PNG + hold duration in ms; repeat per frame")
    ap.add_argument("--gif", metavar="OUT.gif", help="write a looping GIF here")
    ap.add_argument("--mp4", metavar="OUT.mp4", help="write a truecolor MP4 here")
    ap.add_argument("--fps", type=int, default=25, help="MP4 frame rate (default 25)")
    args = ap.parse_args()

    specs = [(path, int(dur)) for path, dur in args.frame]
    if not args.gif and not args.mp4:
        ap.error("specify --gif and/or --mp4")
    if args.gif:
        frames_to_gif([p for p, _ in specs], [d for _, d in specs], args.gif)
        print("wrote %s" % args.gif)
    if args.mp4:
        frames_to_mp4(specs, args.mp4, fps=args.fps)
        print("wrote %s" % args.mp4)


if __name__ == "__main__":
    _main()
