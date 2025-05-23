#! /bin/bash

# Make movie from initial frames of simulation file.
# Requires FFmpeg.

# DEFAULT PARAMETERS

__PYTHON_DEFAULT=$(which python)    # default python executable
__FFMPEG_DEFAULT=$(which ffmpeg)    # default FFmpeg executable

# HELP

usage() {
cat <<< "
Make movie from initial frames of simulation file using modules \`cells.read\`
and \`cells.plot\`.

Requires FFmpeg. (see https://ffmpeg.org/download.html)

SYNOPSIS

    [bash] $(basename "${BASH_SOURCE[0]}") [OPTIONS] SIMULATION_FILE_NAME
    [bash] $(basename "${BASH_SOURCE[0]}") [OPTIONS] -d FRAMES_DIRECTORY_NAME

OPTIONS

    -h  Display this help.

    PLOTTING MODE (NOTE: these do nothing if -d option is used)
    (without argument)
    -r  Make rainbow plot.
    -c  Clear the plot of all cell colouring.
    -v  Display velocities on vertices.
    -n  Highlight number of neighbours per cell.
    -H  Highlight hexatic bond orientational order parameter per cell.
    -T  Highlight translational order parameter per cell.
    -s  Highlight shape index per cell.
    -t  Automatically fit frame to figure using matplotlib.pyplot.tight_layout.
        NOTE:    Size of figures may be unconsistent from frame to frame.
    -W  Read input file with cells.read.ReadWYC rather than cells.read.Read.
        CAUTION: This object is expected to be prone to errors.
        NOTE:    Progress bars will most likely be erroneous.
    (with argument)
    -D  Set frame resolution in dots per inch.
        DEFAULT: (not specified)
    -S  Set font size.
        DEFAULT: (not specified)

    -f  Force check input file consistency (see cells.read.Read).
        NOTE:    This does nothing if -d option is used.

    -d  Set directory with frames from which to make movie.
        DEFAULT: (not specified)
        NOTE:    Frames must be *.png and ordered.

    -y  \"Yes\" to all FFmpeg prompts.

    -p  Set python executable.
        DEFAULT: $__PYTHON_DEFAULT
    -F  Set FFmpeg executable.
        DEFAULT: $__FFMPEG_DEFAULT

    -C  Make movie with H.265 compression.
        NOTE:    This compression is not compatible with all players.
                 FFmpeg must be configured with --enable-gpl --enable-libx265.
"
}

# OPTIONS

unset __FLAGS __RAINBOW __CLEAR __VELOCITIES __NEIGHBOURS __HEXATIC \
    __TRANSLATIONAL __TIGHT_LAYOUT __WYC __DPI __FONT __FORCE_CHECK __DIR \
    __YES __PYTHON __FFMPEG __H265
while getopts "hrcvnHTstWD:S:fd:yp:F:C" OPTION; do
    case $OPTION in
        h)  # help
            usage; exit 1;;
        r)  # rainbow plot
            __RAINBOW=true;
            __FLAGS=${__FLAGS}r;;
        c)  # clear plot
            __CLEAR=true;
            __FLAGS=${__FLAGS}c;;
        v)  # velocity plot
            __VELOCITIES=true;
            __FLAGS=${__FLAGS}v;;
        n)  # neighbours plot
            __NEIGHBOURS=true;
            __FLAGS=${__FLAGS}n;;
        H)  # hexatic plot
            __HEXATIC=true;
            __FLAGS=${__FLAGS}H;;
        T)  # translational plot
            __TRANSLATIONAL=true;
            __FLAGS=${__FLAGS}T;;
        s)  # shape index plot
            __SHAPE_INDEX=true;
            __FLAGS=${__FLAGS}s;;
        t)  # tight layout
            __TIGHT_LAYOUT=true;;
        W)  # use ReadWYC
            __WYC=true;;
        D)  # frame resolution
            __DPI="$OPTARG";;
        S)  # font size
            __FONT="$OPTARG";;
        f)  # force check
            __FORCE_CHECK=true;;
        d)  # frames directory
            __DIR="$OPTARG";;
        y)  # yes to FFmpeg
            __YES=true;;
        p)  # python executable
            if [ "$OPTARG" != "" ]; then __PYTHON=$OPTARG; fi;;
        F)  # FFmpeg executable
            if [ "$OPTARG" != "" ]; then __FFMPEG=$OPTARG; fi;;
        C)  # H.265 compression
            __H265=true;;
    esac
done
shift $(expr $OPTIND - 1)

if [[ -z "$1" && -z "$__DIR" ]]; then
    usage;
    exit 1;
fi

__PYTHON=${__PYTHON-$__PYTHON_DEFAULT}
__FFMPEG=${__FFMPEG-$__FFMPEG_DEFAULT}

# SCRIPT

# make frames
if [[ -z "$__DIR" ]]; then

    __DIR=$(mktemp -d);

    cat >> ${__DIR}/movie.py <<EOF
"""
Save frames from simulation file.

This script was generated by "${BASH_SOURCE[0]}".
"""

from cells.read import Read${__WYC:+WYC} as Read
from cells.read import _progressbar
${__VELOCITIES:+from cells.plot import plot_velocities}
${__NEIGHBOURS:+from cells.plot import plot_neighbours}
${__HEXATIC:+from cells.plot import plot_hexatic}
${__TRANSLATIONAL:+from cells.plot import plot_translational}
${__SHAPE_INDEX:+from cells.plot import plot_p0}

import os, sys, traceback
import numpy as np
${__TIGHT_LAYOUT:+from matplotlib.pyplot import tight_layout}
${__FONT:+from matplotlib import rcParams; rcParams[\"font.size\"] = $__FONT}

# simulation file
r = Read("$(realpath $1)"
    ${__FORCE_CHECK:+, check=True})
# plotted frames
assert(r.skip.size > 0)                                     # check there are saved frames
frames = r.frames[r.frames >= r.t0[0]]
if (np.diff(frames, n=2) != 0).any(): frames = r.t0         # logarithmically spaced saved frames
frames = list(frames[frames <= r.frames[r.skip.size - 1]])  # restrict to existing saved frames

print("Saving frames to temporary directory \"%s\"." % "$__DIR")
_progressbar(0)
for frame in frames:

    # plot
    r.plot(frame
        ${__RAINBOW:+, rainbow=r.t0[0]}
        ${__CLEAR:+, clear=True}
        ${__VELOCITIES:+, override=plot_velocities}
        ${__NEIGHBOURS:+, override=plot_neighbours}
        ${__HEXATIC:+, override=plot_hexatic}
        ${__TRANSLATIONAL:+, override=plot_translational}
        ${__SHAPE_INDEX:+, override=plot_p0})
    # save
    ${__TIGHT_LAYOUT:+tight_layout()}
    while True:
        try:
            r.fig.savefig(
                os.path.join("$__DIR", "%05d.png" % frames.index(frame))
                ${__DPI:+, dpi=$__DPI})
            break
        except SyntaxError:
            # dirty fix to "SyntaxError: not a PNG file" when using concurrent matplotlib instances
            print(traceback.format_exc(), file=sys.stderr)
            pass

    _progressbar((frames.index(frame) + 1)/len(frames))

EOF

    # save frames
    $__PYTHON ${__DIR}/movie.py || exit 1;

fi

# make movie
# https://stackoverflow.com/questions/20847674/ffmpeg-libx264-height-not-divisible-by-2
__MOVIE=${1:-movie}
__MOVIE=${__MOVIE%.*}${__FLAGS:+.$__FLAGS}.mp4
$__FFMPEG ${__YES:+-y} -r 5 -f image2 -s 1280x960 -pix_fmt yuv420p \
    -pattern_type glob -i "${__DIR}/*.png" ${__H265:+-vcodec libx265 -crf 28} \
    -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2:color=white" \
    $__MOVIE

