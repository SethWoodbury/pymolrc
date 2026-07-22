# Making RFdiffusion3 trajectory figures in PyMOL — principles & gotchas

Reference notes for building Figure-1-style stills and storyboard movies from
an RFdiffusion3 (RFd3) enzyme-design diffusion trajectory, on top of
`rfd3_movie()` (loading/gradient/defer-builds) and `style_fixed()` /
`apply_camera()` / `draw_connectors()` / `catalytic_sel_from_motif()` /
`add_custom_bond()` (all in `.pymolrc`, see the command reference in the main
README). Each item below cost real debugging time — read this before
re-discovering any of them the hard way.

## 1. `reverse=1` flips FRAMES, not STATES

`rfd3_movie(..., reverse=1)` remaps PyMOL's *movie frame* -> *state* order so
that frame 1 shows noise and frame N shows folded. It does **not** reorder the
underlying states: raw file state 1 is still the folded structure, and state
N is still noise, regardless of `reverse`.

Anything that takes a literal state index — above all
`cmd.create(new, obj, source_state, target_state)` — must use the real state
number, not the frame number. Copying `n_states` expecting "the folded frame"
silently grabs noise instead. Frame-based operations (`cmd.frame`, `cmd.png`,
`cmd.orient` at the current frame) follow the `reverse` mapping and are fine
as-is.

**This has only been confirmed on `_noisy_model_` trajectory files.**
`rfd3_movie()` does not special-case a `_denoised_model_` trajectory — it's
loaded and rendered exactly the same way, and the filename is only ever used
to auto-locate the design JSON — so don't assume the same state-1-is-folded
convention holds for a denoised trajectory without checking. `rfd3_movie` now
auto-runs `rfd3_check_frames(obj, reverse=...)` at the end of every call, which
measures the actual per-state geometry (item 3) and prints whether your
`reverse=` setting is correct for that specific file; call it by hand any time
you want to re-verify.

## 2. Disconnected fragments ignore transparency

RFd3 fixed-motif residues are frequently tiny, bond-less fragments (a lone
`CB` plus placeholder atoms, no backbone at all). PyMOL renders bond-less
atoms via the `nonbonded` (crosses) or `lines` representations, and **those
representations ignore `sphere_transparency` and `stick_transparency`
entirely.** No number of `set ..._transparency, 0` calls makes them look
solid if they're actually being drawn as `nonbonded`/`lines` under the hood.

Fix (this is what `style_fixed()` does): `hide nonbonded` + `hide lines` on
the fixed selection first, then `show sticks`/`show spheres` fresh, and only
set transparency to 0 *last*, after the representation is rebuilt as
sticks/spheres.

## 3. RFd3 noise is a COLLAPSED blob, not an expanded cloud

Measured on a real 100-state run: folded-state max radius (from centroid)
≈ 28 Å, noise-state max radius ≈ 10 Å. Noise is a tight Gaussian ball near
the origin; the chain *expands* as it folds. This is the opposite of the
intuitive "explosion" picture.

Consequence for cameras: the **folded frame is the largest state and the one
that can clip the render** — frame/zoom the camera to the folded frame, not
to the noise frame. `apply_camera(..., span_states=True, n_states=n_states)`
does this automatically (takes rotation + zoom from the folded frame, then
holds that fixed camera across the whole trajectory).

## 4. Don't toggle representations per-frame on the live multi-state object

`rfd3_movie` runs the trajectory object under `defer_builds_mode 3` (caches
per-state representations for fast scrubbing). Calling `show cartoon` (or any
other representation toggle) on the diffusing chain *inside a per-frame loop*
forces what looks like a full-trajectory rebuild on every call, and hard-hangs
PyMOL — unrecoverable with Ctrl-C; you need `pkill -9 -f pymol` (or
`killall -9 PyMOL`) from a separate terminal to recover.

Fix: never toggle representations per-frame on the deferred object. For any
one-off or per-state representation change (e.g. a cartoon on just the final
folded frame), pull that one state out first with
`cmd.create(tmp_name, obj, source_state, 1)` and work on the resulting
single-state copy instead — it's cheap and never touches the deferred
object's cache. Switching between pre-built alternatives via
`cmd.enable`/`cmd.disable` on separate objects is also cheap; it's the
*rebuild* that's expensive, not having multiple objects around.

## 5. Cartoon on an all-UNK chain traces through the sidechain cloud

Setting `cartoon_trace_atoms, 1` on a chain that's entirely `resn UNK` threads
the ribbon spline through *every* atom present per residue — including the V0..V8
placeholder "sidechain cloud" atoms, not just CA — producing a stringy mess
instead of a clean backbone ribbon.

Two fixes: (a) build the cartoon on a CA-only `cmd.create` copy (strip
everything but `name CA` first), or (b) — usually better — render any
"final, clean" panel from a real sequence-assigned output model (not the raw
`UNK` trajectory), where `show cartoon` builds correctly out of the box. Run
`cmd.dss()` first on that model for secondary-structure assignment.

## 6. The clean output model shares the trajectory's coordinate frame

Because the theozyme is *fixed* through diffusion, the motif residues sit at
identical (or extremely close) coordinates in both the raw noisy trajectory
and a real, sequence-assigned "clean" output model of the same design. That's
what makes `catalytic_sel_from_motif()` precise: it coordinate-matches each
trajectory motif residue to its nearest real residue in the clean model,
typically to well under 1 Å. Use it to identify catalytic residues in a clean
model instead of a raw distance cutoff off the ligand (e.g. "within 4.5 of
resn LIG"), which grabs the entire first shell and over-colors.

## 7. `.pymolrc` functions are called bare, not via `cmd.*`

`cmd.extend(name, function)` registers a function for PyMOL's command
language (so you can type it at the `PyMOL>` prompt), but it does **not**
attach it as a `cmd.<name>` Python attribute. A script run via PyMOL's `run`
command shares the same execution namespace as `.pymolrc` (both load through
the same mechanism), so call these as bare globals: `rfd3_movie(...)`,
`style_fixed(...)`, `apply_camera(...)`, not `cmd.rfd3_movie(...)`.

Also: `.pymolrc` is already auto-loaded by PyMOL at startup. Don't
`cmd.run(os.path.expanduser("~/.pymolrc"))` a second time inside a script —
PyMOL's `run` command tries to compile a file with no recognized extension as
plain Python, and chokes on the `python ... python end` block syntax that
only its own startup loader understands.

## 8. GIF assembly has two quantization/timing traps

- **A single global palette sampled sparsely (e.g. every 10th frame) can wash
  out whichever frames it wasn't built from** — if your longest-held or most
  important frame (a final "clean" panel, say) lands outside the sample, it
  gets mapped onto a palette built from other frames and looks off. Prefer a
  **per-frame local palette** (quantize each frame independently) — bigger
  files, much lower color error. `scripts/frames_to_movie.py::frames_to_gif`
  does this.
- **GIF is hard-capped at 256 colors per frame**, full stop. For anything
  where color fidelity matters, also emit a truecolor MP4 alongside the GIF —
  `scripts/frames_to_movie.py::frames_to_mp4`.

## 9. ffmpeg's concat demuxer can silently drop your first hold

Feeding ffmpeg's `-f concat` demuxer a file list with a per-entry `duration`
is fast to write, but can emit an empty edit-list entry for the first clip
(`media_time: -1` in the resulting container); players honor that as a gap
and open on black, dropping the intended hold on frame 1. Trying to patch
this with `setpts`/`settb` filters tends to collapse the holds instead of
fixing the gap.

More reliable: build a plain **numbered image sequence** — duplicate held
frames as **symlinks** (instant, no pixel copy) rather than real file copies
— and encode that with the `image2` demuxer + `-framerate` instead of
`concat`. This is what `frames_to_mp4` in `scripts/frames_to_movie.py` does.
Verify the result with `ffprobe -show_entries stream=start_time
-show_entries format=duration out.mp4`: `start_time` should read exactly
`0.000000`, and total frame count should equal `fps × total_duration`.

## 10. You can't see a headless render — verify offline first

Running `pymol -cq script.py` produces no visual feedback loop; you don't see
the image until it's written to disk. Before spending render time (especially
ray-traced, multi-panel batches), it's much cheaper to check things
numerically:

- Parse the structure file directly (`gemmi`, or plain mmCIF/PDB text) to
  confirm atom counts per selection, which residues sit near the
  ligand/metal, and per-state extents (radius from centroid) — this is how
  the noise-vs-folded size relationship in item 3 was actually confirmed, not
  guessed.
- Given a candidate `get_view` matrix, you can check whether every atom of
  interest projects inside the render frame at the actual render aspect
  ratio (which may differ from the interactive window's) before ray-tracing:
  project each atom through the view rotation, and check
  `|x|/(depth·tan(fov/2)·aspect) ≤ 1` and `|y|/(depth·tan(fov/2)) ≤ 1`, with
  `depth = camera_distance − z`, against the near/far clip planes.
- Sanity-check `catalytic_sel_from_motif` matches (item 6) by confirming the
  trajectory's fixed motif residues and the clean model's matched residues
  coincide to within a small tolerance before trusting the resulting
  selection for coloring/labeling.

Cheap math beats a slow ray-trace of the wrong thing.
